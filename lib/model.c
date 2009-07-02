/* model.c
 * 
 * Slightly higher-level data manipulation routines
 * These routines generally use both data.c and geometry.c.  They are
 * the machinery of the model itself, doing things like handling 
 * vertex neighbor location and trimming.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Craig DeForest, 2004-2007
 * 
 * You may modify and/or distribute this software under the temrs of
 * the Gnu Public License, version 2.  You should have received a copy
 * of the license with this software, in the file COPYING to be found
 * in the top level directory of the distribution.  You may obtain
 * additional copies of the licesnse via http://www.gnu.org or by
 * writing to the Free Software Foundation, 59 Temple Place - Suite
 * 330, Boston, MA 02111-1307 USA.
 *
 * The software comes with NO WARRANTY.
 * 
 * You may direct questions, kudos, gripes, and/or patches to the
 * author, Craig DeForest, at "deforest@boulder.swri.edu".
 *
 * This file is part of the FLUX 2.2 release (22-Nov-2008).
 * 
 */
#include "model.h"
#include <stdio.h>
#include <math.h>

#include <unistd.h> /* for isatty -- delete if needed */
extern FILE *gl_outf; /* for debugging -- delete */

#include "data.h"
#include "geometry.h"
#include "io.h"

char *code_info_model="%%%FILE%%%";

/**********************************************************************
 **********************************************************************
 ***** Neighbor routines
 ***** Vertices keep track of their neighbors.  Some reciprocal
 ***** actions are handy to have.  Here they are.
 *****
 ***** Neighbor distances are determined by closest approach of the
 ***** neighbor field line in the projected perpendicular plane.  
 */

/**********************************************************************
 * world_check
 * 
 * Performs some rudimentary consistency checks on the world, and enforces
 * some trivial conditions.    Returns 0 if nothing was done, a positive
 * number if the world was "fixed-up", and a negative number if a real 
 * (unfixable) problem was found.
 *
 * In particular:
 * 
 *  No two adjacent vertices on a fluxon can occupy the same location -- 
 *  if this is found, then one of them is deleted.
 *
 * This is intended to be a general cleanup pass for consistency checks
 * after ingestion or, if necessary, after a timestep. 
 * 
 * Other checks that would be useful include:
 *   - checks for invalid neighbors
 *   - checks for crossing a photospheric boundary
 * ...?
 *
 */
static int world_check_code;

static long w_c_springboard(FLUXON *fl, int lab, int link, int depth) {
  VERTEX *v;

  /* skip magic boundary fluxons */
  if( FL_ISDUMMY(fl) )
    return 0;


  /**********
   * Check that either there are no vertices, or that both start and end are defined.
   */
  if( (fl->start == 0) ^ (fl->end == 0) ) {
    fprintf(stderr, "world_check ERROR: fluxon %d has inconsistent start & end vertices (%d vs %d)!\n",fl->label, fl->start, fl->end);
    world_check_code = -1;
    return 0;
  }

  if( (fl->start != 0 && fl->start == fl->end) ) {
    fprintf(stderr, "world_check ERROR: fluxon %d has only one VERTEX (%d) - start & end must be different.\n",fl->label, fl->start->label);
    world_check_code = -1;
    return 0;
  }

  /**********
   * Check that no vertices are on top of one another.
   * This is a correctible error -- delete the second of the two unless it's the end -- then 
   * delete the first of the two.
   */
  for(v=fl->start; v && v->next && v != fl->end; v=v->next) {
    if( sqrt(cart2_3d(v->x, v->next->x) / (norm2_3d(v->x)+norm2_3d(v->next->x))) < EPSILON ) {

      if(!v->prev && !v->next->next) {
	fprintf(stderr,"world_check ERROR: fluxon %d has only two verts (%d & %d) and they are at the same location!\n\tCan't delete either without invalidating the fluxon.\n",fl->label, v->label, v->next->label);
	world_check_code = -1;
	return 0;
      }
      
      fprintf(stderr,"world_check WARNING: fluxon %4d: verts %4d (%s) & %4d (%s) are at the same location - deleted %4d...\n",fl->label, v->label, (v->prev ? "mid" : "beg"), v->next->label, (v->next->next ? "mid" : "end"), (v->next->next ? v->next->label : v->label));
      delete_vertex( v->next->next ? v->next : v );
      if(world_check_code == 0) 
	world_check_code = 1;
      return 0;
    }
  }
  
  return 0;
}

int world_check(WORLD *a) {
  world_check_code = 0;
  
  safe_tree_walker(a->lines, fl_lab_of, fl_all_ln_of, w_c_springboard, 0);
}


/**********************************************************************
 * world_update_ends
 * Updates the ends of each fluxon in turn...
 */
static long w_u_e_springboard(FLUXON *fl, int lab, int link, int depth) {
  /* Skip magic boundary fluxons */
  if( FL_ISDUMMY(fl) )
    return 0 ;

  fluxon_update_ends(fl);
  //if(fl->fc0->world->verbosity)
  //  printf("world_update_ends - fluxon %d\n",fl->label);
  return 0;
}

void world_update_ends(WORLD *a) {
  safe_tree_walker(a->lines,fl_lab_of, fl_all_ln_of, w_u_e_springboard,0);
  //tree_walker(a->lines,fl_lab_of, fl_all_ln_of, w_u_e_springboard,0);
}

/**********************************************************************
 * fluxon_update_ends
 * Checks and updates the magic boundary vertices at the end of a fluxon.
 * 
 * If the world has automatic open vertex handling ON, then 
 * fluxon_update_ends also checks to see if any new ends are
 * required, and creates 'em.
 *
 * It also handles U-loop exit (via the auto_open mechanism)
 */
void fluxon_update_ends(FLUXON *f) {
  WORLD *w = f->fc0->world;

  /*** End condition updates ***/
   if(f->fc0->bound) {	
	(*(f->fc0->bound))(f->start);
   } else if(f->fc0->world->default_bound) {
     (*(f->fc0->world->default_bound))(f->start);
   }

   if(f->fc1->bound) {
	(*(f->fc1->bound))(f->end);
   } else if(f->fc1->world->default_bound) {
     (*(f->fc1->world->default_bound))(f->end);
   }

   if(w->auto_open) {
     fluxon_auto_open(f);
   }

   if(f->plasmoid) {
     fluxon_plasmoid_cleanup(f);
   }
}

/**********************************************************************
 * fluxon_plasmoid_cleanup - given a FLUXON, check that it is a plasmoid,
 * and if it is, check for trivial plasmoid removal.  The problem to be
 * solved is that trivial loop plamoids in gas-free simulations tend to
 * shrink without limit.  If they shrink far enough, eliminate 'em.
 */
void fluxon_plasmoid_cleanup(FLUXON *f) {
  if(!f)
    return;
  if(!f->plasmoid)
    return;
  if(trivloop(f)) { // trivloop is in geometry.c
    delete_fluxon(f);
  }
  
}

/**********************************************************************
 * fluxon_auto_open - given a FLUXON, check whether any of its vertices
 * need to be opened (are outside the open-field boundary).
 *
 * If so, do so.
 * 
 * Does nothing unless the associated WORLD auto_open flag is set.
 */
void fluxon_auto_open(FLUXON *f) {

  WORLD *w = f->fc0->world;
  
  VERTEX *v = f->start;
  NUM x0[3];
  int v_ct = 0;
  NUM r2;
  NUM *xcen;
  
  if( ! (w->auto_open) )
    return;

  r2 = w->fc_ob->locale_radius;
  r2 *= r2;
  
  xcen = w->fc_ob->x;
  
  /**************************************************
   * Check for vanishing U-loops
   */

  if( (f->fc0 == w->fc_ob) &&
      (f->fc1 == w->fc_oe) ) {
    // It's a U-loop...

    // Trivial U-loops...
    if(f->v_ct < 4) {
      printf(" Fluxon %d is a trivial U-loop.  Deleting...\n",f->label);
      delete_fluxon(f);
      return;
    } 
  
    // Tiny U-loops...
    //    printf("Checking fluxon %d for triviality (start is at (%g,%g,%g))\n",
    //   f->label, f->start->x[0], f->start->x[1], f->start->x[2]);
    if(trivloop(f)) { // trivloop is in geometry.c
      //printf("Deleting %d...\n",f->label);
      delete_fluxon(f);
      return;
    }
    
  }

  /**************************************************
   * Check for newly opened lines
   */

  for(v=v->next; v && v->next; (v=v->next) && v_ct++) {
    diff_3d(x0,v->x,xcen);
    fflush(stdout);
    
    if( norm2_3d(x0) > r2 ) {
     printf("v%d (on %d) ",v->label,f->label);
      printf("Open!  ");
      fflush(stdout);
      
      /******************************
       * Cut the fluxon into two.  The new 
       * endpoints are placed on the open surface
       * at the points of intersection.
       */
      
      if( f->plasmoid ) {
	/* Plasmoid case - turn it into a U-loop.  This gets 
	 * a little complex becasue we have to rotate
	 * the endpoint of the plasmoid around to the new 
	 * true endpoints.
	 */

	VERTEX *vp = v->prev;  // vp gets the last one inside the boundary
	VERTEX *vn = v;        
	char stop = 0;

	/******************************
         * Walk around the plasmoid to find how many vertices are outside.
         */
	for( vn=v ; vn != vp  && !stop 
	       ; vn = ( (vn->next && vn->next->next) ? vn->next : f->start ) ) {
	  diff_3d(x0, vn->x, xcen);
	  stop = (norm2_3d(x0) <= r2);
	}
	
	/******************************
	 * All the vertices were outside -- delete the fluxon.
         */
	if(vn==vp) {
	  printf("Deleting fluxon %d\n",f->label);
	  delete_fluxon(f);
	  return;
	}

	/******************************
         * Normal case -- Splice the plasmoid into a single fully-open fluxon.
         */

	// Wrap around the end of the plasmoid if we have to...
	if(! vp->prev)
	  vp= f->end->prev;
	if(! vp->next)
	  vp= f->start->next;

	if(!vn->prev)
	  vn = f->end->prev;
	if(!vn->next)
	  vn = f->start->next;


	// Create a real loop, cutting out the dummy vertices
	f->end->prev->next   = f->start->next;
	f->start->next->prev = f->end->prev;

	// Unlink the dummy first and last vertices from the plasmoid
	f->start->next = NULL;
	f->start->prev = NULL;
	f->end->next=NULL;
	f->end->prev=NULL;

	  
	// Delete the dummy first and last vertices (also decrements v_ct)
	delete_vertex(f->start);
	delete_vertex(f->end);

	// Delete the intervening vertices between vp and vn
	// (also decrementing v_ct)
	{
	  VERTEX *vv;
	  for(vv=vp->next; vv != vn; vv=vp->next) 
	    delete_vertex(vv);
	}

	// Link in the proper (new) start and end
	f->start = vn;
	f->end = vp;

	// Remove the circular links
	f->start->prev = NULL;
	f->end->next = NULL;

	// clear plasmoid flag and relink flux concentrations...
	f->plasmoid = 0;
	{
	  WORLD *w;
	  w = f->fc0->world;
	  if(f->fc0 != w->fc_pb) 
	    printf("Hey! plasmoid beginning wasn't the world plasmoid start\n");
	  if(f->fc1 != w->fc_pe) 
	    printf("Hey! plasmoid end wasn't the world plasmoid end\n");

	  // Unlink from the plasmoid flux concentration placeholders
	  f->fc0->lines = tree_unlink(f, fl_lab_of, fl_start_ln_of);
	  f->fc1->lines = tree_unlink(f, fl_lab_of, fl_end_ln_of);

	  // Link to the open flux concentration placeholders
	  f->fc0 = w->fc_ob;
	  f->fc1 = w->fc_oe;
	  f->fc0->lines = tree_binsert(f->fc0->lines, f, fl_lab_of, fl_start_ln_of);
	  f->fc1->lines = tree_binsert(f->fc1->lines, f, fl_lab_of, fl_end_ln_of);
	}

	return;
	
      } else /* not a plasmoid */ {
	
	/****************
         * Not a plasmoid... cut into two separate fluxons
         */

	NUM x1[3];
	NUM newx[3];
	NUM rp_n, rn_n, r_n_diff, alpha;
	VERTEX *Vnext;
	VERTEX *Vprev;
	
	// Vnext is the next vertex located inside the sphere.  If there aren't any at all, 
	// then we simply open the fluxon at the current point.
	{
	  int inside = 0;
	  for(Vnext = v->next; Vnext && !inside;) {
	    NUM x2[3];
	    diff_3d(x2, Vnext->x, w->fc_oe->x);
	    rn_n = norm_3d(x2)/w->fc_oe->locale_radius;
	    inside = (rn_n < 1);
	    if (!inside) {
	      Vnext = Vnext->next;
	    }
	  }
	}
	
	if(!Vnext) {
	  
	  printf("Hmmmm, no Vnext --- purging to end of fluxon...\n");

	  /****************************************
	   * The whole rest of the fluxon is open.  Truncate it and connect it to the 
	   * open end vertex.
           *
           *
           * **** CASE: whole last half is outside the boundary **** 
	   */
	  while(v->next) {
	    delete_vertex(v->next);
	  }
	  f->end = v;

	  if(f->fc1 != w->fc_oe) {
	    f->fc1->lines = tree_unlink(f,fl_lab_of, fl_end_ln_of);
	    f->fc1 = w->fc_oe;
	    w->fc_oe->lines = tree_binsert(w->fc_oe->lines, f, fl_lab_of, fl_end_ln_of);

	    // Invoke the end boundary condition on the newly opened endpoint.
	    (*(f->fc1->bound))(f->end);
	  }
	  

	} else /* Vnext exists - normal case */ { 
	  
	  /****************************************
	   * There's at least one non-outside vertex following, so we need to 
	   * cut the current fluxon in two.
	   */
	  
	  NUM xp[3];
	  
	  // Calculate Vprev...
	  Vprev = v->prev;
	  diff_3d(xp, v->prev->x, w->fc_ob->x);
	  rp_n = norm_3d(xp)/w->fc_ob->locale_radius;
	  
	  if(rp_n > 1) {

	    printf("Hmmm .. no previous vertex.  Opening start of fluxon...\n");
	    
	    /* If there are no previous inside vertices, but there are following 
	     * inside vertices, then we need to open the beginning...
             *
             *
             *  **** CASE: whole first half is outside the boundary ****
             */
	    if(Vnext) {
	      while(Vnext->prev) {
		delete_vertex(Vnext->prev);
	      }
	      f->start = Vnext;
	    } else {
	      while(v->prev) {
		delete_vertex(v->prev);
	      }
	      f->start = v;
	    }

	    if(f->fc0 != w->fc_ob) {
	      f->fc0->lines = tree_unlink(f,fl_lab_of, fl_start_ln_of);
	      f->fc0 = w->fc_ob;
	      w->fc_ob->lines = tree_binsert(w->fc_ob->lines, f, fl_lab_of, fl_start_ln_of);
	    }
	    
	    
	  }  else /* Previous vertex is inside */ {

	    /*Vprev is inside, so we cut.
	     * Now v (and possible next links) is outside, Vprev is
	     * inside, and Vnext is inside.  We need exactly two
	     * vertices between Vnext and Vprev.  If there is only
	     * one, add one.  If there is more than one, excise
	     * 'em. */
	    
	    if(Vprev->next == Vnext->prev) {

	      // Add an extra vertex, colocated.  

	      VERTEX *Vnew = new_vertex(0, v->x[0],v->x[1],v->x[2], f);
	      add_vertex_after(f, v, Vnew);

	      {
		long n;
		for(n=0; n < v->neighbors.n; n++) {
		  VERTEX *vn = (VERTEX *)(v->neighbors.stuff[n]);
		  dumblist_add(&(Vnew->neighbors),(void *)vn);
		  dumblist_add(&(vn->nearby), (void *)Vnew);
		}
	      }
	    
	    } else if( Vprev->next->next != Vnext->prev ) {

	      // Trim vertices out of the middle

	      VERTEX *vv;
	      for( vv=Vprev->next->next->next; vv != Vnext; vv=vv->next ) {
		delete_vertex(vv->prev);
	      }
	    }
	    
	    // Now we should have Vprev->a->b->Vnext.
	    printf("Vprev is %d; next is %d; next is %d; next is %d; Vnext is %d\n",
	       Vprev->label, Vprev->next->label, Vprev->next->next->label, Vprev->next->next->next->label, Vnext->label);
	    
	    /* Now reposition the intermediate vertices approximately
	       on the sphere, by truncating their respective segments. */
	    {
	      NUM ra_n, rb_n;
	      NUM alpha_a, alpha_b;
	      NUM xx[3];
	      diff_3d(xx, Vprev->next->x, w->fc_ob->x);
	      ra_n = norm_3d(xx)/w->fc_ob->locale_radius;
	      
	      diff_3d(xx, Vnext->prev->x, w->fc_oe->x);
	      rb_n = norm_3d(xx)/w->fc_oe->locale_radius;
	      
	      alpha_a = (1 - rp_n)/(ra_n-rp_n) * 0.9;
	      alpha_b = (1 - rn_n)/(rb_n-rn_n) * 0.9;
	      //printf("alpha_a = %g; alpha_b=%g\n",alpha_a, alpha_b);
	      
	      diff_3d(xx, Vprev->next->x, Vprev->x);
	      scale_3d(xx, xx, alpha_a);
	      sum_3d(Vprev->next->x, xx, Vprev->x);
	      
	      diff_3d(xx, Vnext->prev->x, Vnext->x);
	      scale_3d(xx, xx, alpha_b);
	      sum_3d(Vnext->prev->x, xx, Vnext->x);
	      xx[0] = xx[0] * 1;
	    }
	    
	    /* Finally, create a new fluxon and cut the old one
	     * between the two intermediate vertices.  Convenience
	     * block...*/
	    {
	      int n;
	      VERTEX *vv;
	      FLUXON *Fnew = new_fluxon( f->flux,
					 w->fc_ob,
					 f->fc1,
					 0,
					 0
					 );
	      // Cut the linked lists.
	      Vnext->prev->prev->next = 0;
	      Vnext->prev->prev = 0;
	      Fnew->end = f->end;
	      f->end = Vprev->next;
	      Fnew->start = Vnext->prev;

	      
	      // Switch the original fluxon to its new endpoint
	      f->fc1->lines = tree_unlink(f, fl_lab_of, fl_end_ln_of);
	      f->fc1 = w->fc_oe;
	      w->fc_oe->lines = tree_binsert(w->fc_oe->lines, f, fl_lab_of, fl_end_ln_of);
	      
	      // Cut off the vertices in the middle
	      n=0;
	      for(vv=f->start; vv; vv=vv->next) n++;
	      f->v_ct = n;

	      // Set up the vertices...
	      n=0;
	      for(vv=Vnext->prev; vv; vv=vv->next) {
		n++;
		vv->line = Fnew;
	      }
	      Fnew->v_ct = n;

	      // Finally - clean up new fluxon (old will be cleaned below)
	      // and ensure we keep scanning down the new fluxon.
	      v=Vnext;
	      f=Fnew;
	    } /* End of convenience block */
	  } /* end of normal-case check (prev. vertex is inside) */
	} /* end of normal-case check (next vertex is inside) */
      } /* end of non-plasmoid opening code */
    } /* end of open-check loop */
  } /* end of vertex for-loop */
  //printf("\n");
} /* end of fluxon_auto_open */

  
/**********************************************************************
 * world_update_neighbors
 * 
 * Calls fluxon_update_neighbors to process the whole world!
 */
static WORLD *gl_a;  /* springboard world parameter */
static char gl_gl;   /* springboard global flag */
static void ((**gl_f_funcs)()); /* springboard functions list */

/* Helper routine for snarfing up the closest vertex in each fluxon to
   a particular vertex-- called via tree_walk, in data.c. It does not
   winnow, and it does not check the current fluxon. that should be
   taken care of in the next step when it does two regular neighbor
   searches.  */
static VERTEX *vertex_for_closest_v;
static DUMBLIST *snarfer_workspace;

static long closest_v_snarfer(FLUXON *f, int lab_of, int ln_of, long depth) {
  VERTEX *v2; //current vertex in the fluxon
  VERTEX *v_keep; //vertex we want to snarf
  int i; //dummy counter
  NUM r ; //distance b/w v2 and v
  NUM rmax=0; //distance accumulator
  
  v2 = f->start; //(don't do first and last)
  for (i=0; i < f->v_ct; i++){
    r = cart_3d(v2->x,vertex_for_closest_v->x);
      //sqrt(((v2->x[0] - vertex_for_closest_v->x[0])**2)+((v2->x[1] - vertex_for_closest_v->x[1])**2)+((v2->x[2] - vertex_for_closest_v->x[2])**2));
    if (r > rmax){
      v_keep = v2;
    }

    if (v2->next) { //for i=v_ct-1, there is no vnext
      v2 = v2->next;
    }

  }
  dumblist_quickadd (snarfer_workspace, v_keep);
  return 0;
}

static long fast_world_gather_neighbors(VERTEX *v, int lab_of, int ln_of, long depth){
  if( FL_ISDUMMY(v->line) ) //added. hopefully works.
    return 0 ;

  if (v == v->line->start || v == v->line->end)
    return 0;

  static DUMBLIST *workspace =0;
  void **foo;
  int i;
  int n;
  int verbosity = v->line->fc0->world->verbosity;
  long passno = ++(v->line->fc0->world->passno);
  DUMBLIST *vn;
  FLUXON *f;

  vertex_for_closest_v=v;
  
  f = v->line;
  while(f->all_links.up)
    f=f->all_links.up;
  
  if(!workspace) 
    workspace = new_dumblist();
  
  workspace->n = 0;
  
  snarfer_workspace = workspace;
  tree_walk(f, fl_lab_of, fl_all_ln_of,closest_v_snarfer);
  
  //is this the right way?
  vn = &(v->neighbors);
  vn->n=0;
  dumblist_snarf(vn,workspace);
  //printf("%d  ",v->label);
  
  return 0;
}


static long w_u_n_springboard(FLUXON *fl, int lab, int link, int depth) {
  /* Skip magic boundary fluxons */
  if( FL_ISDUMMY(fl) )
    return 0 ;

  if(fl->fc0->world->verbosity >= 5) dump_all_fluxon_tree(gl_a->lines);
  if(fl->fc0->world->verbosity) {
    printf(" %d",fl->label);
    fflush(stdout);
  }
  return fluxon_update_neighbors(fl, gl_gl);
}

void world_update_neighbors(WORLD *a, char global) {
  gl_a = a;
  gl_gl = global;
  if (gl_gl == -1) { 
    /*do fast neighbor search then do 2 iterations of regular neighbor search. */
    tree_walk(a->vertices,v_lab_of,v_ln_of,fast_world_gather_neighbors);
    printf("done with fast gather neighbors\n");
    gl_gl = 0;
    tree_walker(a->lines, fl_lab_of, fl_all_ln_of, w_u_n_springboard, 0);
    printf("done with first gather neighbors\n");
    tree_walker(a->lines, fl_lab_of, fl_all_ln_of, w_u_n_springboard, 0);
    printf("done with second gather neighbors\n");
    gl_gl = -1;
  } else {
    tree_walker(a->lines, fl_lab_of, fl_all_ln_of, w_u_n_springboard, 0);
  }
}




/**********************************************************************
 * fast_world_update_neighbors
 * 
 * Initializes the neighbors as the closest vertex to you from each
 * fluxon.
 */

static long f_w_u_n_springboard(FLUXON *fl, int lab, int link, int depth) {
  /* Skip magic boundary fluxons */
  if( FL_ISDUMMY(fl) )
    return 0 ;

  if(fl->fc0->world->verbosity >= 5) dump_all_fluxon_tree(gl_a->lines);
  if(fl->fc0->world->verbosity) {
    printf(" %d",fl->label);
    fflush(stdout);
  }
  return fast_fluxon_update_neighbors(fl, gl_gl);
}

void fast_world_update_neighbors(WORLD *a, char global) {
  gl_a = a;
  gl_gl = global;
  tree_walker(a->lines, fl_lab_of, fl_all_ln_of, f_w_u_n_springboard, 0);
}


/**********************************************************************
 * world_update_mag updates the magnetic forces for the whole world 
 */
static long w_u_m_springboard(FLUXON *fl, int lab, int link, int depth) {
  /* Skip magic boundary fluxons */
  if( FL_ISDUMMY(fl) )
    return 0 ;

  return fluxon_update_mag(fl,gl_gl, gl_f_funcs);
}

int world_update_mag(WORLD *a, char global) {
  gl_a = a;
  gl_gl = global;
  gl_f_funcs = a->f_funcs;

  init_minmax_accumulator(a);
  return tree_walker(a->lines, fl_lab_of, fl_all_ln_of, w_u_m_springboard, 0);
  finalize_minmax_accumulator(a);
}


/**********************************************************************
 * world_fluxon_length_check makes sure that there are at least
 * 4 vertices in any fluxon. We don't care about fluxons that only 
 * have 2 vertices and there aren't any that only have 1.
 */
static long check_fluxon_length(FLUXON *f, int lab, int link, int depth){
  //printf("in the function upper\n");
    //printf("fluxon %d, count is %d\n",f->label,f->v_ct);
    fflush(stdout);
  if (f->v_ct == 3) { //fix_curvature on the middle vertex
    fix_curvature(f->start->next,1e-9,1e-10);
    //printf("in the function lower\n");
    fflush(stdout);
  }
  return 0;
}

void world_fluxon_length_check(WORLD *w, char global){
  gl_a = w;
  gl_gl = global;  
  //printf("in world_fluxon_length_check\n");
  fflush(stdout);
  //tree_walker(tree, label offset, link offset, function, depth)
  tree_walker(w->lines, fl_lab_of, fl_all_ln_of, check_fluxon_length,0);
}


/**********************************************************************
 * world_collect_stats
 *
 * Collects neighbor and force stats about the world.  
 * Mainly useful as a relaxation diagnostic; see 'stats' in the PDL
 * XS routines.
 *
 */
static struct VERTEX_STATS gl_st;
static long w_c_s_springboard(FLUXON *fl, int lab, int link, int depth) {
  /* Skip magic boundary fluxons */
  if( FL_ISDUMMY(fl) )
    return 0 ;

  fluxon_collect_stats(fl,&gl_st);
  return 0;
}

VERTEX_STATS *world_collect_stats(WORLD *a) {
  gl_st.n = 0;
  gl_st.f_acc = 
    gl_st.f_max =
    gl_st.f_tot_acc =
    gl_st.f_tot_max =
    gl_st.n_acc =
    gl_st.n_max = 0;
  tree_walker(a->lines, fl_lab_of, fl_all_ln_of, w_c_s_springboard, 0);
  return &gl_st;
}

void fluxon_collect_stats(FLUXON *fl, VERTEX_STATS *st) {
  VERTEX *v = fl->start->next;
  while(v && v->next) {
    NUM f;
    st->n++;

    f = norm_3d(v->f_t);
    if( f > st->f_max ) 
      st->f_max = f;
    st->f_acc += f;

    f = v->f_s_tot + ( v->f_v_tot + v->prev->f_v_tot ) / 2;
    if( f > st->f_tot_max )
      st->f_tot_max = f;
    st->f_tot_acc += f;

    if(v->neighbors.n > st->n_max)
      st->n_max = v->neighbors.n;
    st->n_acc += v->neighbors.n;

    v = v->next;
  }
}
   

/**********************************************************************
 *
 * world_relax_step -- relax the world according to a given set of force
 * laws.  You feed in a world and a tau-timestep.
 */
static NUM gl_t;

static long w_ca_s_springboard(FLUXON *fl, int lab, int link, int depth) {
  /* Skip magic boundary fluxons */
  if( FL_ISDUMMY(fl) )
    return 0 ;

  fluxon_calc_step(fl, gl_t);
  return 0;
}

static long w_r_s_springboard(FLUXON *fl, int lab, int link, int depth) {
  /* Skip magic boundary fluxons */
  if( FL_ISDUMMY(fl) )
    return 0 ;

  fluxon_relax_step(fl, gl_t);
  return 0;
}

void world_relax_step(WORLD *a, NUM t) {
  gl_t=t;
  tree_walker(a->lines, fl_lab_of, fl_all_ln_of, w_ca_s_springboard, 0);
  gl_t=t;
  init_minmax_accumulator(a);
  tree_walker(a->lines, fl_lab_of, fl_all_ln_of, w_r_s_springboard, 0);
  finalize_minmax_accumulator(a);
}

/**********************************************************************
 * fluxon_update_neighbors
 * 
 * Calls vertex_update_neighbors to process a whole fluxon. 
 * This is wasteful, as it throws away the Voronoi vertex information
 * from each cell.  Use fluxon_update_mag instead -- that handles
 * physics too.
 */
int fluxon_update_neighbors(FLUXON *fl, char global) {
  int i=0;
  VERTEX *v =fl->start;
  int verbosity = fl->fc0->world->verbosity;

  if(verbosity>=2) printf("fluxon_update_neighbors... (gl=%d), fluxon %d\n",global,fl->label);

  if(FL_ISDUMMY(fl))
    return 0;

  /* First thing - make sure the end conditions are up to date... */
  if(fl->fc0->bound) {	
    (*(fl->fc0->bound))(fl->start);
  }
  if(fl->fc1->bound) {
    (*(fl->fc1->bound))(fl->end);
  }

  v=v->next;

  while(v->next) {   
    if(verbosity>=3)  printf("\tfluxon_update_neighbors... vertex %d\n",v->label);

    vertex_update_neighbors(v,global || (v->neighbors.n == 0));
    if(verbosity>=4)   fdump_fluxon(stdout,fl,0);
    if(verbosity>=3)   printf("=============\n\n");

    v=v->next;
    i++;
  }

  return 0;
}

/**********************************************************************
 * fluxon_update_mag
 * 
 * Calls vertex_update_neighbors and physics functions to process
 * a whole fluxon.  Differs from update_neighbors in that you feed
 * in a list of pointers to force functions that need calling at each
 * node.  At the end of the call, each fluxon contains not only 
 * updated neighbor information but also fresh force values.  
 * Returns the largest and smallest forces and 
 * force/neighbor-distance ratios (in that order) in an array of four
 * NUMs (useful for computing timesteps).
 */

int fluxon_update_mag(FLUXON *fl, char global, void ((**f_funcs)())) {
  int i=0;
  WORLD *w = fl->fc0->world;
  int verbosity = w->verbosity;
  VERTEX *v = fl->start;

  if(verbosity >= 2)  printf("fluxon_update_mag (fluxon %d): %c",fl->label,(verbosity==2?' ':'\n'));

  /*
   * Set forces on the end vertex to zero, since it has no 
   * bend and no following segment and hence does not get 
   * looped over.
   */
  {
    VERTEX *v = fl->end;
    v->f_v[0] = v->f_v[1] = v->f_v[2] = 0;
    v->f_v_ps[0] = v->f_v_ps[1] = v->f_v_ps[2] = 0;
    v->f_s[0] = v->f_s[1] = v->f_s[2] = 0;
    v->f_s_tot = v->f_v_tot = 0;
    v->r_s = v->r_v = -1;
  }

  
  for(i=0, v=fl->start; v->next; v=v->next, i++) {
    void (**f_func)();
    NUM r;
    int ii;
    HULL_VERTEX *vertices;
    if(verbosity >= 3) { printf("V%4d ",v->label); fflush(stdout); }

    /* Update neighbor map and establish neighbors' (r,a) variables, 
     *  which will be used by the physics funcs! */

    vertices = vertex_update_neighbors(v,global || (v->neighbors.n==0)); 

    /* Zero the vertex's force accumulators */
    v->f_v[0] = v->f_v[1] = v->f_v[2] = 0;
    v->f_s[0] = v->f_s[1] = v->f_s[2] = 0;
    v->f_s_tot = v->f_v_tot = 0;
    v->r_s = v->r_v = -1;

    /* Find smallest 'r' associated with this vertex */
    {
      int i;
      v->r_cl = -1;
      for(i=0;i<v->neighbors.n;i++)
	if( (v->r_cl < 0) 
	    || (v->r_cl > (((VERTEX **)(v->neighbors.stuff))[i])->r)
	    )
	  v->r_cl = (((VERTEX **)(v->neighbors.stuff))[i])->r;
    }
    v->r_ncl = v->r_cl;


    if(fl->fc0->world->verbosity >= 3) 
      printf("---forces:\n");

    {
      int ignore_segment_forces = 
	( ((!v->prev) && funky_fl_b(v->line->fc0))
	  ||
	  (  ((!v->next) || (!v->next->next))  && funky_fl_b(v->line->fc1) )
	  );

     /* Accumulate forces and relevant lengthscales */
    for(f_func = &f_funcs[0]; *f_func; f_func++) 
      (**f_func)(v,vertices,ignore_segment_forces);
    }

    if(fl->fc0->world->verbosity >= 3)
      printf("\n");

  }
  if(verbosity >= 3) printf("\n");


  /* Find closest approach radius, and calculate total forces */
  
  if(verbosity >= 2) printf("Radius - fluxon %d: %c",fl->label,(verbosity==2?' ':'\n'));

  for(v=fl->start;v->next;v=v->next) {
    NUM f[3], fns, fnv, fn;
    NUM a;
  
    diff_3d(f,v->next->x,v->x);
    a = norm_3d(f);
    if(a<v->r_cl && a>0)
      v->r_cl=a;
    
    if(v->prev) {
      diff_3d(f, v->x, v->prev->x);
      a = norm_3d(f);
      if(a<v->r_cl && a>0)
	v->r_cl = a;

      sum_3d(f,v->prev->f_s,v->f_s);
      scale_3d(f,f,0.5);
    } else {
      cp_3d(f,v->f_s);
    }
    
    fns = (v->r_s > 0) ? (norm_3d(f)) : 0;
    fnv = (v->r_v > 0) ? (norm_3d(v->f_v)) : 0;
    
    sum_3d(v->f_t,f,v->f_v);
    v->f_n_tot = norm_3d(v->f_t);

    /* Accumulate force and curvature data for the whole world.... */
    vertex_accumulate_f_minmax(v,w);

  }

  return 0;
}



  
/**********************************************************************
 * init_minmax_accumulator, vertex_accumulate_f_minmax
 *
 * Accumulates global stats by vertex.  You can call 
 * world_accumulate_minmax to get 'em all, or do it vertex-by-vertex
 * if you are already iterating through the arena.
 *
 * fluxon_update_mag calls the accumulator on a vertex-by-vertex
 * basis, which is currently a waste in the parallel case (because the 
 * daughter processes' accumulators will be discarded).
 *
 */
void init_minmax_accumulator(WORLD *w) {
  w->f_min = -1;
  w->f_max = -1;
  w->fr_min = -1;
  w->fr_max = -1;
  w->ca_min = -1;
  w->ca_max = -1;
  w->ca_acc = 0;
  w->ca_ct = 0;
}

void finalize_minmax_accumulator(WORLD *w) {
  w->max_angle  = RAD2DEG * acos(sqrt(fabs(w->ca_min)))            + 90 * (w->ca_min < 0);
  w->mean_angle = RAD2DEG * acos(sqrt(fabs(w->ca_acc / w->ca_ct))) + 90 * (w->ca_acc < 0);
}

void vertex_accumulate_f_minmax(VERTEX *v, WORLD *w) {
  NUM fr;
  NUM cangle;
  NUM b1hat[3];
  NUM b2hat[3];

  if(!w)
    w = v->line->fc0->world;

  /* Check max, min, etc. */
  if(w->f_min < 0 || v->f_n_tot < w->f_min)
    w->f_min = v->f_n_tot;
  if(w->f_max < 0 || v->f_n_tot > w->f_max)
    w->f_max = v->f_n_tot;
  
  fr = v->f_n_tot / v->r_cl;
  if(w->fr_min < 0 || fr < w->fr_min)
    w->fr_min = fr;
  if(w->fr_max < 0 || fr > w->fr_max)
      w->fr_max = fr;
  
  if (v->prev && v->next) {
    diff_3d(b1hat, v->x, v->prev->x);
    diff_3d(b2hat, v->next->x, v->x);
    cangle = inner_3d(b1hat, b2hat);
    cangle *= fabs(cangle);
    cangle /= (norm2_3d(b1hat) * norm2_3d(b2hat));
    
    // Calculate using only the cosine of the 
    // angle (faster than taking a gazillion acos's)
    if (w->ca_min < 0 || cangle < w->ca_min )
      w->ca_min = cangle;
    if (w->ca_max < 0 || cangle > w->ca_max )
      w->ca_max = cangle;
    
    w->ca_acc += cangle;
    w->ca_ct++;
  }
}

/**********************************************************************
 * fluxon_relax_step
 *ARD NOTES
 * fluxon relax step is now split into two parts. 
 * fluxon_calc_step and fluxon_relax_step.
 *
 *PREVIOUS NOTES
 * 
 * Relaxes the vertex positions of a fluxon.  You must have pre-loaded
 * the fluxon's forces with fluxon_update_mag.
 * 
 * You feed in an amount of time (really, a force-to-offset scaling 
 * factor) to advance, and everything advances toward it.  There is
 * no vertex_relax_step because it's such a simple operation -- though
 * one might want to write one later for extensibility.
 * 
 * The physics routines stuff forces into the f_s and f_v fields;
 * update_mag does the averaging for f_s and stuffs it into f_v for you,
 * so by the time you see this you should only have to look at f_v.
 *
 * Variables used for scaling:
 *     d  - harmonic mean of the various distances of neighbors and such
 *     b  - mean of the field strengths of the two segments on the vertex
 *     s  - stiffness coefficient from the force calculation
 *     ds - ratio of sum to difference of local vertex's forces and next/previous vertices'
 */


/*** fastpow is quite fast for small integer powers, slightly slower
 *** for fractional. Acceleration fails for large negative powers.
 *** Powers within 0.0001 of an integer get treated as integers.
 ***/
NUM fastpow( NUM num, NUM exponent ) {
  NUM out = 1.0;
  NUM recip;
  int exp = (exponent + 100.0001) - 100;

  if( exponent - exp < 0.0001 && exponent-exp > -0.0001) {
    if(exponent<0) {
      recip = 1.0/num;
      while(exp<0) {  exp++; out *= recip;  }
    } else 
      while(exp>=0.999) { exp--; out *= num; }
    return out;

  } else {

    return pow(num,exponent);

  }
}

NUM calc_stiffness(VERTEX *v) {
  NUM f_denom;
  if(!v->prev || !v->next)
    return 0 ;
  f_denom = v->f_v_tot + 0.5 * (v->f_s_tot + v->prev->f_s_tot);
  return (f_denom == 0) ? 1 : ( norm_3d(v->f_t)  / (1e-9 + f_denom));
}
    
void fluxon_calc_step(FLUXON *f, NUM dt) {
  VERTEX *v = f->start;
  WORLD *w = f->fc0->world;
  NUM a[3];
  NUM total[3];
  NUM stiffness;
  int verbosity = f->fc0->world->verbosity;

  a[2] = a[1] = a[0] = 0;


  for(v=v->next; v && v->next; v=v->next) {
    NUM r_cl= v->r_cl;
    NUM foo[3];
    NUM d,d1;



    // Calculate the harmonic mean of the relevant segment lengths
    // and the closest neighbor approach distance.
    
    diff_3d(foo,v->next->x,v->x);
    d = norm_3d(foo);
    
    diff_3d(foo,v->x,v->prev->x);
    d1 = norm_3d(foo);

    if(r_cl <= 0) {
      fprintf(stderr,"ASSERTION FAILED!  Negative distance %g on vertex %d!\n",r_cl,v->label);
      fprintf(stderr,"vertex has %d neighbors\n",v->neighbors.n);
      r_cl = 1; /* hope the problem goes away */
    }

    d = 4/(1/d + 1/d1 + 1/r_cl + 1/v->prev->r_cl);

    stiffness = calc_stiffness(v);
    if(stiffness == 0) stiffness = 1e-6;

    if(verbosity >= 3)  printf("fluxon %d, vertex %d: x=(%g,%g,%g).  v->r_cl=%g,  r_cl=%g,  stiffness = %g, f_t=(%g,%g,%g)[%g]\t",f->label,v->label, v->x[0],v->x[1],v->x[2], v->r_cl, r_cl, stiffness, norm_3d(v->f_t), v->f_t[0],v->f_t[1],v->f_t[2],norm_3d(v->f_t));
    
    if(stiffness > 1.00001) {
      fprintf(stderr,"fluxon %d, vertex %d: stiffness = %g, >1!  This is allegedly impossible! You've got trouble, gov\n",f->label,v->label,stiffness);
      fflush(stdout);
      fflush(stderr);
    }
    
    /* 
     * We now have the force per unit length on the vertex, with or without B scaling
     * depending on the force law.  Use the scaling information stored in the WORLD to 
     * generate a step.
     */

    {
      NUM fac = 1;
      NUM b_mag_mean;
      WORLD *w = f->fc0->world;

      if(w->step_scale.b_power) {
	b_mag_mean = v->b_mag;
	if(v->prev) {
	  b_mag_mean += v->prev->b_mag;
	  b_mag_mean *= 0.5;
	}
	
	fac *= fastpow(b_mag_mean,w->step_scale.b_power);
      }

      fac *= fastpow(d,             w->step_scale.d_power);
      fac *= fastpow(stiffness,  w->step_scale.s_power);

      /*
       * Handle acceleration of steps when everything's moving the same way
       */
      if(w->step_scale.ds_power && v->prev && v->prev->prev && v->next && v->next->next) {
	NUM a[3];
	NUM b[3];
	NUM pfrac,nfrac;

	sum_3d(a, v->f_t, v->prev->f_t);
	diff_3d(b, v->f_t, v->prev->f_t);
	pfrac = norm_3d(b)/norm_3d(a);
	
	sum_3d(a, v->f_t, v->next->f_t);
	diff_3d(b, v->f_t, v->next->f_t);
	nfrac = norm_3d(b)/norm_3d(a);
        fac *= fastpow( 2.0/(pfrac + nfrac), w->step_scale.ds_power );
      }

      if(finite(fac)) 
	scale_3d(a, v->f_t, dt * fac );	
      else 
	if(verbosity >= 3) 
	  printf("DS FAC NOT FINITE - SETTING to 1.0  ");

    }

    {
      /* Check step length -- absolute maximum is 
       * 0.49 of the distance to the nearest
       * boundary-condition neighbor, or all the distance to the nearest
       * non-boundary-condition neighbor.
       */
      if(v->next && v->prev) {
	NUM diff[3];
	NUM l1,l2,steplen;

	/*** Check distance to buddies on the same fluxon ***/
	diff_3d(diff,v->x,v->prev->x);
	l1 = norm_3d(diff);
	diff_3d(diff,v->next->x,v->x);
	l2 = norm_3d(diff);
	if(l2<l1)
	  l1=l2;

	/*** Check distance to image vertex ***/
	{
	  int i;
	  char flag = 0;
	  for (i=0; i<v->neighbors.n && !flag; i++) {
	    VERTEX *vn = ((VERTEX *)(v->neighbors.stuff[i]));
	    diff_3d(diff,v->x,vn->x);
	    l2 = norm_3d(diff);
	    //if( vn == v->line->fc0->world->image )
	      l2 *= 0.49;
	    if(l1<0 || l2 < l1)
	      l1 = l2;
	  }
	}

	steplen = norm_3d(a);
	if(l1>0 && steplen > l1)
	  scale_3d(a,a,l1/steplen);
      }
    }

    // OK we've finished calculating the step - stuff it into vertex->plan_step

    v->plan_step[0] = a[0];
    v->plan_step[1] = a[1];
    v->plan_step[2] = a[2];

    // Finally, accumulate force and curvature info
    vertex_accumulate_f_minmax(v, w);
  }
}

// ARD - Moved expand_lengthwise and expand_via_neighbors routines in file
// so that they are defined before they are called.

void expand_lengthwise(DUMBLIST *workspace, int start_idx, long passno) {
  int i;
  long n = workspace->n;
  for(i=start_idx;i<n;i++) {
    VERTEX *v = ((VERTEX *)(workspace->stuff[i]));
    if(v->prev && v->prev->passno != passno) {
      v->prev->passno = passno;
      dumblist_quickadd(workspace, v->prev);
    }
    if(v->next && v->next->passno != passno) {
      v->next->passno = passno;
      dumblist_quickadd(workspace,v->next);
    }
  }
}

void expand_via_neighbors(DUMBLIST *workspace, int start_idx, long passno) {
  int i;
  int j;
  long n = workspace->n;
  for(i=start_idx;i<n;i++) {
    VERTEX *v = ((VERTEX *)(workspace->stuff[i]));
    DUMBLIST *nlist;

    nlist = &(v->neighbors);
    for(j=0; j<nlist->n; j++) {
      VERTEX *vv;
      //printf("j=%d, ",j);
      vv = ((VERTEX *)(nlist->stuff[j]));
      //printf,("j=%d workspace.n=%d, neighbors.n=%d, neighbors.size=%d, workspace.size=%d \n",j,n,v->neighbors.n,v->neighbors.size, workspace->size);
      //printf("vv=%d, ",vv->label);
      //printf("vv->line=%d, ",vv->line->label);
      //printf("vv->passno=%d, ",vv->passno);
      //printf("passno=%d\n",passno);
      if(vv->line && vv->passno != passno) {
	//printf("here1 ");
	vv->passno = passno;
	//printf("here2 ");
	dumblist_quickadd(workspace, vv);
	//printf("here3 ");
      }
    }

    //printf("here4 ");
    nlist = &(v->nearby);
    for(j=0;j<nlist->n;j++) {
      VERTEX *vv;
      vv = ((VERTEX *)(nlist->stuff[j]));
      if(vv->line && vv->passno != passno) {
	vv->passno = passno;
	dumblist_quickadd(workspace, vv);
      }
    }
  }
}


/*******************************
 * vertex_enforce_photosphere 
 * Force a vertex to obey photospheric boundaries.
 * 
 * You feed in a VERTEX and a PHOTOSPHERE, and 
 * the VERTEX is forced into agreement with the PHOTOSPHERE.
 * 
 * There is no checking that the VERTEX agrees with all 
 * photospheres -- it is only moved in the most expedient 
 * way to agree with the current one.
 * 
 * Only planar and cylindrical photospheres are currently
 * supported.
 * 
 * This is necessary because, while normal photospheric 
 * neighborly action (image charges) usually works, occasional
 * quirks of relaxation or other boundaries can bring a 
 * point outside a defined PHOTOSPHERE.  This mops up those 
 * glitches. 
 *
 */
void vertex_enforce_photosphere(VERTEX *v, PHOTOSPHERE *p) {
  if(!p->type)
    return;

  switch(p->type) {
  case PHOT_PLANE: {
    /* Check for and eliminate photospheric plane crossings */
    NUM x1[3];
    NUM x1z;
    diff_3d(x1,v->x,p->plane->origin);
    
    x1z = inner_3d(x1,p->plane->normal);
    if(x1z < 0) {
      NUM delta[3];
      scale_3d(delta, p->plane->normal, - (x1z * 1.001) / norm2_3d(p->plane->normal));
      sum_3d(x1,x1,delta);
      sum_3d(v->x, x1, p->plane->origin);
    }
  }
    break;

  case PHOT_CYL: {
    /* Check for and eliminate crossings of the cylindrical boundary */
    NUM x1[3];
    NUM x1r2;
    NUM M[9];
    NUM xrot[3];
    NUM r2;
    r2 = norm2_3d(p->plane->normal);
    
    sum_3d(x1,p->plane->origin,p->plane->normal);

    /* We can't be outside the cylinder unless we're outside the inscribed */
    /* sphere -- which is much faster to check, so check that first */
    if(norm2_3d(x1) > r2) { 
      projmatrix(M,p->plane->origin,x1);
      
      /* Find the perpendicular-to-cylinder component of the vertex location */
      diff_3d(x1, v->x, p->plane->origin);
      mat_vmult_3d(xrot, M, x1);
      x1r2 = norm2_2d(xrot);
      
      /* If it's outside the cylinder, move it perpendicular to the axis  */
      /* (the Z axis in the rotated coordinates) until it's just inside the */
      /* cylinder. */
      if(x1r2 > r2) {
	scale_2d(xrot, xrot, (1.0 - 1e-5) * sqrt(r2/x1r2));
	vec_mmult_3d(x1, M, xrot);
	sum_3d( v->x, x1, p->plane->origin);
      }
      
    }
  }
    break;

    case PHOT_SPHERE: {
      NUM x1[3];
      NUM x1r;
      
      /* sphere center is at p->origin; radius is p->plane->normal[0]; outer/inner flag is sign of p->plane->normal[1]. */
      diff_3d(x1, v->x, p->plane->origin);
      x1r = norm_3d(x1);

      if( p->plane->normal[1] >= 0 ? (x1r < p->plane->normal[0]) : (x1r > p->plane->normal[0]) ) {
	if(p->plane->normal[1] >= 0) {
	  scale_3d(x1, x1, (1 + 1e-5) * p->plane->normal[0] / x1r);
	} else {
	  scale_3d(x1, x1, (1 - 1e-5) * p->plane->normal[0] / x1r);
	}
	sum_3d(v->x, x1, p->plane->origin);
      }
    }
    break;
	
    default: break;
      
    } /* end of photospheric checking switch */
}       
			      


// ARD - New fluxon_relax_step routine. Step calculation is now in 
// fluxon_calc_step

void fluxon_relax_step(FLUXON *f, NUM dt) {
  VERTEX *v = f->start;
  VERTEX *vtmp;
  NUM total[3];
  static DUMBLIST *workspace =0;
  long passno;
  int idx, j, i;
  WORLD *world = f->fc0->world;
  int verbosity = world->verbosity;
  int n_coeffs = world->n_coeffs;
  POINT3D step, tmp_step;
  VERTEX *vert_neigh;

  if(!workspace) 
    workspace = new_dumblist(); 

  if(!v) {
    fprintf(stderr,"WARNING: fluxon_relax_step called with a fluxon (%d) that has no start vertex!\n",f->label);
    return;
  }
    
  for(v=v->next; v && v->next; v=vtmp) {
    vtmp = v->next;
    /* vtmp is an assignment. This is here because sometimes the code
     * adds vertices because of the photopsheres. Those new vertices
     * don't have v->plan_step initialized so when those new vertices
     * get stepped they go to garbage-land. This way, the next vertex
     * you go to will always be one that was already initialized. */

    // ARD - Add in standard step - scaled in case we do something strange!

    scale_3d(step, v->plan_step, world->coeffs[0]);

    ///////////////////// ARD - Now add relative contributions of neighbors

    dumblist_clear(workspace);

    // Avoid duplication with the passno marking mechanism...
    passno = ++(world->passno);
    idx = 0;

    if (n_coeffs > 1) {

      // ARD - add current vertex as first element of dumblist
      dumblist_quickadd(workspace, v);

      for (i=1;i<n_coeffs;i++) {
	long newidx = workspace->n;
	long n_found = 0;

	expand_via_neighbors(workspace, idx, passno);
	expand_lengthwise(workspace, idx, passno);
	idx = newidx;

	// Accumulate temporary step for the current set of neighbors.
	// Start accumulating at the END of the previous batch, to keep
	// track of neighbor generational averages.
	// Ignore anyone whose stiffness is worse than 5%.
	tmp_step[0] = tmp_step[1] = tmp_step[2] = 0.0;

	for (j=idx;j<workspace->n;j++) {
	  vert_neigh = ((VERTEX *)(workspace->stuff[j]));
	  if(vert_neigh->next && vert_neigh->prev) {
	    NUM stiffness= calc_stiffness(vert_neigh);
	    if(stiffness < 0.05) {
	      sum_3d(tmp_step, tmp_step, vert_neigh->plan_step);
	      n_found++;
	    }
	  }
	}
	
	// Get the mean for this batch, scale it, and add.
	// The if clause avoids division-by-zero if none exist at this tier.
	if(n_found) {
	  scale_3d(tmp_step, tmp_step, (world->coeffs[i] / (n_found)));
    	
	//	printf("step:   %lf %lf %lf\ntmp_step: %lf %lf %lf\n\n", step[0], step[1], step[2], tmp_step[0], tmp_step[1], tmp_step[2]);
	
	  sum_3d(step, step, tmp_step);
	}
      }
    }
     
    //    printf("step:   %lf %lf %lf\n", step[0], step[1], step[2]);
    //fflush(stdout);
    
    if(finite(step[0]) && finite(step[1]) &&finite(step[2])) {
      sum_3d(v->x,v->x,step);
    } else {
	printf("NON_FINITE OFFSET! f_t=(%g,%g,%g), vertex=%d (ignoring)\n",v->f_t[0],v->f_t[1],v->f_t[2],v->label);
    }

    vertex_enforce_photosphere( v, &(world->photosphere) );
    vertex_enforce_photosphere( v, &(world->photosphere2) );

    if(v->prev && v->prev->prev) {
      vertex_accumulate_f_minmax(v->prev, world);
    }

    if(verbosity >= 3)    
      printf("after update: x=(%g,%g,%g)\n",v->x[0],v->x[1],v->x[2]);
  }
  
  /******************************
   * Done with mid-fluxon update -- now adjust the start and end positions 
   * depending on boundary condition.  (Currently only line-tied, open-sphere, and
   * open-plane boundary conditions are supported).
   */

  fluxon_update_ends(f);
}
    
  
/**********************************************************************
 * vertex_update_neighbors 
 * 
 * Updates a vertex's neighbor list using gather_neighbor_candidates
 * followed by winnow_neighbor_candidates and hull_neighbors.
 * 
 * Returns the hull vertex list that is used in the final winnowing
 * process, since it's needed for some of the force laws (see 
 * vertex_update_mag).
 * 
 * The 'global' flag can have several values, and the behavior changes
 * accordingly. fast, faster, and gonzo are defined in model.::
 * 
 *   -1 - quasi-global operation (requires cleanup afterward): the closest 
 *        vertex (using regular 3-D Cartesian distance) along each fluxon 
 *        is a candidate for the hull)
 *    0 - "normal" non-global operation (about 300 candidates for the hull)
 *    1 - global operation (every vertex is a candidate for every hull
 *    2 - reduced operation (neighbors of neighbors and neighbors' next/prev only: about 60-70 candidates) ->fast_neighbors
 *    3 - stripped-down operation (neighbors and their next/prev only: about 24 candidates) ->faster_neighbors
 *    4 - gonzo operation (neighbors only: about 8 candidates) -> gonzo_neighbors
 *
 * The stripped-down operations operate much more quickly but at the
 * expense of not finding neighbor candidates under evolution.
 */

static int ptr_cmp(void *a, void *b) { /* Helper routine for sorting by pointer */
  if(a<b) return -1;
  else if(a>b) return 1;
  return 0;
}

static int label_cmp(void *a, void *b) { /* Helper routine for sorting by pointer */

  if( ((VERTEX *)a)->label < ((VERTEX *)b)->label) return -1;
  else if( ((VERTEX *)a)->label > ((VERTEX *)b)->label) return 1;
  return 0;
}

HULL_VERTEX *vertex_update_neighbors(VERTEX *v, char global) {
  DUMBLIST *dl;
  HULL_VERTEX *hv;
  int i,j;
  char cpflag=0;
  DUMBLIST *vn;
  int verbosity = v->line->fc0->world->verbosity;

  if(!v) {
    fprintf(stderr,"Vertex_update_neighbors: null vertex ignored!\n");
    return 0;
  }

  if(!v->next) {
    return 0;
  }


  vn = &(v->neighbors);


#if DEBUG_DUPE_CHECK
  // Debugging check - remove for production
  // Check to see if the VERTEX has any duplicate neighbors!
  {
    int i,j;
    for (i=1;i<vn->n;i++)
      for(j=0;j<i;j++)
	if(vn->stuff[i]==vn->stuff[j]) {
	  printf("ZOINKS! Vertex %d came in with dup neighbor %d (in locations %d and %d)\n",v->label,((VERTEX **)(vn->stuff))[i]->label,i,j);
	}
  }
#endif


  dl = gather_neighbor_candidates(v,global);


  if(verbosity >= 3) {
    printf("vertex_update_neighbors: Vertex %d has %d candidates: ",v->label,dl->n);
    if(verbosity >= 4)
      for(i=0;i<dl->n;i++)
	printf("  %d",((VERTEX *)((dl->stuff)[i]))->label);
    printf("\n");
  }

  hv = hull_neighbors(v, dl); /* save hv for return */
  
  if(verbosity >= 3) {
    printf("Hull_neighbors returned %d neighbors: ",dl->n);
    if(verbosity >= 4)
      for(i=0;i<dl->n;i++) 
	printf("  %d",((VERTEX *)((dl->stuff)[i]))->label);
    printf("\n");
  }
  
  /******************************
   * Update 'nearby' lists through brute force
   * Not too expensive since these lists are typically under 10 elements long
   */
  for(i=0;i<vn->n;i++) {
    dumblist_delete( &(((VERTEX *)(vn->stuff[i]))->nearby), v);
  }
  for(i=0;i<dl->n;i++) {
    dumblist_add(    &(((VERTEX *)(dl->stuff[i]))->nearby), v);
  }
  cpflag = 1;
  
  /* Copy the new neighbor list into the vertex. */
  if(cpflag) {
    vn->n = 0;
    dumblist_snarf(vn,dl);
  }

#if DEBUG_DUPE_CHECK
  // Debugging check - remove for production
  // Check to see if the VERTEX has any duplicate neighbors!
  {
    int i,j;
    for (i=1;i<vn->n;i++)
      for(j=0;j<i;j++)
	if(vn->stuff[i]==vn->stuff[j]) {
	  printf("ZOINKS! Vertex %d left with dup neighbor %d (in locations %d and %d)\n",v->label,((VERTEX **)(vn->stuff))[i]->label,i,j);
	}
  }
#endif

  return hv;
}



/**********************************************************************
 * gather_neighbor_candidates
 * 
 * This finds the neighbor candidates of a vertex, by the 
 * neighbors-of-neighbors strategy.  (That is to say, the only 
 * candidates are the neighbor-candidate list of its neighbors
 * on its own field line, and the neighbors of its previous neighbors).  
 * If they have NO neighbors, then the entire
 * fluxon's neighbor set is generated recursively using the 
 * web of neighbors.  This can be slow, and it might miss a neighbor from
 * time to time, but it's not as bad as the prima facie alternative of 
 * always doing a global neighbor search. 
 * 
 * It returns a pointer to its own workspace -- so calling it kills
 * your last instance of its return value, unless you've gone and copied
 * it yourself.
 * 
 * gather_neighbor_candidates should normally be followed immediately by
 * winnow_neighbor_candidates, which sorts and winnows the list.
 * 
 */


/*
 * snarf_filtered
 * Helper routine snarfs up a dumblist but skips over already-grabbed vertices,
 * image vertices, and plasmoid final segments.
 */
static inline void snarf_filtered(DUMBLIST *ws, DUMBLIST *dl, long passno) {
  int i;
  for(i=0;i<dl->n;i++) {
    VERTEX *V = (VERTEX *)(dl->stuff[i]);
    if( V && V->next && V->passno != passno && V->line && ((!V->line->plasmoid) || V->next->next) ) {        
      V->passno = passno;
      dumblist_quickadd(ws, dl->stuff[i]);
    }
  }
}

/* snarf_list
 * Helper routine for snarfing up a set of vertices and its next and prev
 * elements.  You feed in a list of vertices (e.g. the "neighbors" list
 * from a vertex) and a workspace, and the vertices that are neighbors and 
 * nearby to each of the elements of the list are themselves added to the 
 * workspace.
 *
 */

static inline void snarf_list(DUMBLIST *workspace, VERTEX **foo, int n, long passno) {
  int i;

  for(i=0;i<n;i++) {
    VERTEX *V = ((VERTEX *)(foo[i]));
    if(V->passno != passno && V->line) {
      V->passno = passno;
      dumblist_quickadd(workspace,(void *)foo[i]);
    }
    snarf_filtered(workspace,&(foo[i]->neighbors), passno);
    snarf_filtered(workspace,&(foo[i]->nearby), passno);
  }
}
  

/* expand_list
 * Helper routine for expanding a workspace list out:  you feed in 
 * a workspace with a list of vertices, and it gets expanded to include
 * the next and previous item on each of those vertex lists.
 * It also includes the neighbors of all of those items.
 */
static inline void expand_list(DUMBLIST *workspace, long passno) {
  int i,n;
  n = workspace->n;
  for(i=0;i<n;i++) {
    VERTEX *v = ((VERTEX **)(workspace->stuff))[i];
    if(v) {
      if( v->next && v->next->passno != passno )
	v->next->passno = passno;
	dumblist_quickadd(workspace, v->next);
      if( v->prev && v->prev->passno != passno ) 
	v->prev->passno = passno;
	dumblist_quickadd(workspace, v->prev);
    }
  }
}

/* Helper routine for snarfing up whole fluxons -- called via
   tree_walk, in data.c */
static long line_snarfer(FLUXON *f, int lab_of, int ln_of, long depth) {
  VERTEX *v;
  for(v = f->start; v->next && ((!f->plasmoid) || v->next->next); v=v->next) {
    dumblist_quickadd(snarfer_workspace, v);
  }
  return 0;
}

DUMBLIST *gather_neighbor_candidates(VERTEX *v, char global){
  static DUMBLIST *workspace =0;
  void **foo;
  int i;
  int n;
  int verbosity = v->line->fc0->world->verbosity;
  long passno = ++(v->line->fc0->world->passno);
  //printf(". ");

  if(verbosity >= 2) {
    printf("passno=%d...  ",passno);
  }

  if(!v)        /* Paranoia */
    return;

  if(!workspace) 
    workspace = new_dumblist();

  workspace->n = 0;

  
  /* Gather the candidates together 
   * This step has had a longish history -- I started out 
   * grabbing lots and lots of stuff, but "just" neighbors-of-neighbors
   * (and next and prev links) seems to be sufficient.  In pathological
   * cases, it might take a couple of timesteps to walk to the appropriate
   * new neighbor, but that doesn't seem to cause problems in practice.
   * 
   * The various "fast" setting of global introduce pathologies of their own
   * by cutting down on the neighbor search.  Use with care but they can speed
   * up the process considerably!
   */

  /* Snarf neighbors */
  if(global != global_neighbors) {
    VERTEX *v1;
    int ilen,iwid;

    dumblist_quickadd(workspace, v);

    if(verbosity>=3) {
      printf("initial seed - %d; ",workspace->n);
    }

    if(global != gonzo_neighbors && global != faster_neighbors) {
      expand_lengthwise(workspace, 0, passno);
    
      if(verbosity >= 3) {
	printf("len - %d; ", workspace->n);
      }
    }
    ilen = workspace->n;
      
    expand_via_neighbors(workspace, 0, passno);
    iwid = workspace->n;
	
    if(verbosity >= 3) {
      printf("wid - %d; ", workspace->n);
    }

    if(global != gonzo_neighbors) {
      
      expand_lengthwise(workspace, ilen, passno);
      
      if(verbosity >= 3) {
	printf("len - %d; ", workspace->n);
      }
      
      
      if(global != faster_neighbors && global != fast_neighbors) {
	expand_via_neighbors(workspace, iwid, passno);
	
	if(verbosity >= 3) {
	  printf ("wid - %d; ", workspace->n);
	}
	
	expand_lengthwise(workspace, 0, passno);
	
	if(verbosity >= 3) {
	  printf("final - %d\n",workspace->n);
	}
      }
    }
  }
  
    
   /* Incremental neighbor searching didn't work -- do a global 
     grab (ouch!).
  */
  if(global==1 || (workspace->n == 0)) {
    FLUXON *f;

    if(verbosity >= 3)
      printf("Using global neighbors...");

    /* Snarf every vertex into the list! (sloooow) */
    f = v->line;

    while(f->all_links.up)
      f=f->all_links.up;

    snarfer_workspace = workspace;

    tree_walk(f, fl_lab_of, fl_all_ln_of, line_snarfer);


    if(verbosity >= 4)
          {
            int i;
            printf("gather_neighbor_candidates:  global neighbor add gives:\n\t");
            for(i=0;i<workspace->n;i++)
      	printf("%d ",((VERTEX *)(workspace->stuff[i]))->label);
            printf("\n");
          }


  }

  /* If there's a photosphere, add the dummy photospheric mirror
   * image. For most conditions, we don't let the first or last segments
   * interact with the photosphere, since they generally intersect the
   * photosphere.  Plasmoids do interact with the photosphere, since
   * they generally don't intersect it. */

  if(v->next && 
     ( ( v->prev && v->next->next )  || v->line->plasmoid )
     )
  {
    PHOTOSPHERE *phot;
    VERTEX *image;
    NUM a;

    /* Here, we find the image location via a helper routine below
     * called image_find. Use the correct photosphere and the vertices
     * image and image->next to get the two image positions. Then
     * stuff those image positions into the workspace which is the
     * list of possible neighbor candidates. Can use the same image
     * and image->next for both photspheres because in the neighbor
     * search, it doesn't care about the bit of memory called image,
     * it just cares about the neighbor list.*/

    if(v->line->fc0->world->photosphere.type) { /* assignment, photosphere1 */
      phot = &(v->line->fc0->world->photosphere);
      image = (v->line->fc0->world->image);
      if(verbosity >= 4){
	printf("using photosphere (type is %d)...",v->line->fc0->world->photosphere.type);
	fflush(stdout);
      }
      image_find(phot,image,v); /*helper routine found below*/
      dumblist_quickadd(workspace, image);
    }
    
    if(v->line->fc0->world->photosphere2.type) { /* assignment, photosphere2 */
      phot = &(v->line->fc0->world->photosphere2);
      image = (v->line->fc0->world->image3);
      if(verbosity >= 4){
	printf("using photosphere2 (type is %d)...",v->line->fc0->world->photosphere2.type);
	fflush(stdout);
      }
      image_find(phot,image,v);
      dumblist_quickadd(workspace, image);
    }
  }

  if(verbosity >= 3) printf("gather_neighbor_candidates: returning %d elements (may include dupes)\n",workspace->n);

  return workspace;
}

void image_find(PHOTOSPHERE *phot, VERTEX *image, VERTEX *v) {
  /* helper routine to generate the image point and stuff it into
   * the correct image point in the world.*/
  NUM a;
  PLANE *p = phot->plane;

  switch(phot->type) {
    PLANE pl;
    POINT3D pt;
    NUM radius;
    NUM M[9]; /*don't need * b/c it is a ptr to the first elem*/
    POINT3D vec2,xd,xrot,irot;
    
  case PHOT_CYL:
    /* arbitrary normal and origin. length of normal is the
     * radius. there are a number of steps involved in this. 1 - find
     * the projection matrix M that rotates the normal to the z-axis,
     * 2 - rotate v->x, 3 - find image charge, 4 - de-rotate image
     * position.  radius norm_3d(p->normal)
     */

    radius = norm_3d(p->normal);   
    sum_3d(vec2,p->origin,p->normal);
    projmatrix(M,p->origin,vec2);

    diff_3d(xd,v->x,p->origin);
    mat_vmult_3d(xrot,M,xd);
    a = norm_2d(xrot);

    pl.origin[0] = xrot[0] * radius / a;
    pl.origin[1] = xrot[1] * radius / a;
    pl.origin[2] = 0;
    pl.normal[0] = xrot[0];
    pl.normal[1] = xrot[1];
    pl.normal[2] = 0;
    scale_3d(pl.normal, pl.normal, 1.0/norm_3d(pl.normal));
    reflect(irot, xrot, &pl);
    vec_mmult_3d(image->x,M,irot);
    sum_3d(image->x,image->x,p->origin);
    
    diff_3d(xd,v->next->x,p->origin);
    mat_vmult_3d(xrot,M,xd);
    a = norm_2d(xrot);

    pl.origin[0] = xrot[0] * radius / a;
    pl.origin[1] = xrot[1] * radius / a;
    pl.origin[2] = 0;
    pl.normal[0] = xrot[0];
    pl.normal[1] = xrot[1];
    pl.normal[2] = 0;
    scale_3d(pl.normal, pl.normal, 1.0/norm_3d(pl.normal));
    reflect(irot, xrot, &pl);
    vec_mmult_3d(image->next->x,M,irot);
    sum_3d(image->next->x,image->next->x,p->origin);
    
    break;
  case PHOT_SPHERE:
    /***********
     * Spherical photosphere - sphere is located at the origin in the 
     * photospheric plane structure; radius is normal[0]. 
     * For external spheres, set normal[1] >= 0.  For internal spheres,
     * set normal[1] < 0.
     */
    
    /* Construct the sub-vertex point on the sphere */
    diff_3d(&(pt[0]),      v->x, &(p->origin[0]));  /* pt gets (x - origin) */
    scale_3d(&(pt[0]), &(pt[0]), p->normal[0]/norm_3d(&(pt[0])));  /* Scale to be on the sphere */
    sum_3d(&(pt[0]), &(pt[0]), p->origin);   /* Put in real space */
    scale_3d(&(pt[0]), &(pt[0]), 2.0);                             
    diff_3d(image->x, &(pt[0]), v->x);  /* Put reflection in image */
    
    /***** Do for the other point too *****/
    diff_3d(&(pt[0]), v->next->x, p->origin);
    scale_3d(&(pt[0]), &(pt[0]), p->normal[0]/norm_3d(&(pt[0])));
    sum_3d(&(pt[0]),&(pt[0]),p->origin);
    scale_3d(&(pt[0]), &(pt[0]), 2.0);
    diff_3d(image->next->x, &(pt[0]), v->x); 
    
    break;
    
  case PHOT_PLANE:
    reflect(image->x, v->x, p);
    reflect(image->next->x, v->next->x, p);
    break;
  default:
    fprintf(stderr,"Illegal photosphere type %d!\n",phot->type);
    exit(13);
  }
} /*end of image_find subroutine*/


/**********************************************************************
 * winnow_neighbor_candidates - Remove duplication for the massive
 * collection of neighbor candidates returned by gather_neighbor_candidates
 * (above).
 *
 * We'd like to sort by field line ID, to avoid duplicates along the 
 * same field line -- but unfortunately that breaks some of the more
 * interesting topologies.  In particular, "S"-shaped field lines may
 * interact with a second field line (or with themselves!) in more than
 * one place.  So we first remove the obvious duplications by sorting
 * on vertex pointer, then traverse each of the potential neighbors' field
 * lines to find the closest vertex to the current one, and sort
 * (thereby eliminating duplicates) again.  Note that some of the 
 * "neighbors" will turn out to be this vertex itself!  The current
 * vertex is, of course, omitted from the final list.
 * 
 * winnow_neighbor_candidates is likely the hot spot for the model
 * as a whole; this should be thought about carefully but is not
 * currently.
 *
 */

/* Helper function -- integer comparison for pointers */
int winnow_cmp_1(void *a, void *b) {
  if( (long)a < (long)b ) 
    return -1;
  if( (long)a > (long)b )
    return 1;
  return 0;
}

void winnow_neighbor_candidates(VERTEX *v, DUMBLIST *horde) {
  int i;
  VERTEX **vnp;
	
  /* Traverse each member of the horde to find the locally closest
   * vertex */
  for(i=0;i<horde->n;i++) {
    NUM d0,d_prev,d_next;
    VERTEX *vn;
    int j;

    vnp = ((VERTEX **)horde->stuff) + i;
    j=0;
    do {
      vn = *((VERTEX **)vnp);

      /* Find distance to the actual segment in the list */
      d0 = fl_segment_dist(v,vn);

      /* If the actual segment isn't valid, trash it. */
      if(d0<0) {
	dumblist_rm(horde,i);
	i--;
	break;
      }
	
      /* Find distances to the next and prev segments */
      d_next = fl_segment_dist(v,vn->next);
      d_prev = fl_segment_dist(v,vn->prev ? vn->prev : (vn->line->plasmoid ? vn->line->end->prev->prev : NULL) );
      /* Step forward or backward if necessary.  Fork if the field line
       * bends toward us.
       */
      if( (d_next > 0) && (d_next < d0) ) {
	if( (d_prev > 0) && (d_prev < d_next)) 
	  dumblist_quickadd(horde,vn->prev);
	(*vnp) = vn->next;
	if( (*vnp)->line->plasmoid && !((*vnp)->next->next) ) {
	  (*vnp) = (*vnp)->line->start;
	}
      } else if( (d_prev > 0) && (d_prev < d0)) {
	(*vnp) = (vn->prev ? vn->prev : (vn->line->plasmoid ? vn->line->end->prev->prev : NULL));
      } else
	break;     /* Only exit from loop! */
    } while(1);
  }

  /* Remove this vertex from the list, just in case it's there!*/
  dumblist_delete(horde,v);

  /* Sort again, to remove any duplicates again. */
  dumblist_sort(horde,winnow_cmp_1);
}


/**********************************************************************
 * hull_neighbors - Out of a list of neighbor candidates,
 * hull_neighbors generates a true neighbor list.  It also produces
 * planar hull information and a host of other useful 2-D projected
 * goodies.
 *
 * Each of the neighbor candidate segments is projected into the
 * current segment's plane.  Then the Voronoi cell is constructed
 * around the central point in the projection plane.  Fortunately,
 * routines in geometry.c take care of the geometric details.
 *
 * As long as a hull calculation is being done anyhow, hull_neigbors
 * keeps the valuable hull vertex information and returns it in a
 * buffer of (x,y) pairs.  The buffer is ephemeral -- it gets
 * overwritten each time hull_neighbors is called -- so use it while
 * you can.
 */

HULL_VERTEX *hull_neighbors(VERTEX *v, DUMBLIST *horde) {
  int i;
  int a,b;
  POINT3D x0;
  int verbosity = v->line->fc0->world->verbosity;

  static HULL_VERTEX *voronoi_buf = 0;
  static int voronoi_bufsiz = 0;

  if(v->line->fc0->world->verbosity >= 5) 
    printf("Entering hull_neighbors...\n");

  if(!v->next) {
    fprintf(stderr,"pick_neighbors: assertion failed -- vertex should have a next member!\n\tGiving up...\n");
    return 0;
  }

  /* Project the horde into the plane perpendicular to v's line segment.
   * The projected vectors go into the vertices' scratch space.    
   */
  if(v->line->fc0->world->verbosity>=5)
    printf("hull_neighbors calling project_n_fill...\n");

  project_n_fill(v, horde); /* in geometry.c */

  if(v->line->fc0->world->verbosity>=5)
    printf("hull_neighbors back from project_n_fill...\n");


  /* Grow the buffer if necessary. */
  if(voronoi_bufsiz <= horde->n*2) {
    voronoi_bufsiz = horde->n*4;

    //fprintf(stderr,"   expanding voronoi_buf: horde->n is %d; new bufsiz is %d   ",horde->n,voronoi_bufsiz);

    if(voronoi_buf)
      localfree(voronoi_buf);
    voronoi_buf = (HULL_VERTEX *)localmalloc((voronoi_bufsiz)*sizeof(HULL_VERTEX),MALLOC_VL);

    if(!voronoi_buf) {
      fprintf(stderr,"Couldn't get memory in hull_neighbors!\n");
      exit(99);
    } 

  }

  if(verbosity > 1) {
    printf(".");
    
    if(verbosity >= 5)
      printf("V%4d: (%7.3g, %7.3g, %7.3g) -- (%7.3g, %7.3g, %7.3g)\n",v->label,v->x[0],v->x[1],v->x[2], v->next->x[0],v->next->x[1],v->next->x[2]);
  }

  /* Find the 2-D hull.  Don't want rejects. */
  /*  hull_2d(voronoi_buf,horde,0); */
  hull_2d_us(voronoi_buf, horde, v);
  
  if(horde->n==0) {
    printf("VERTEX %d: hull_2d gave up! That's odd....\n",v->label);
  }

  if(verbosity >= 5){
    printf("hull_neighbors:  hull trimming gives:\n");
    for(i=0;i<horde->n;i++) 
      printf("%d ",((VERTEX *)(horde->stuff[i]))->label); fflush(stdout);
    printf("\n");
  }
  
  return voronoi_buf;
}

/**********************************************************************
 ** fix_proximity
 ** fluxon_fix_proximity
 ** global_fix_proximity
 ** 
 ** You feed in a VERTEX with updated neighbor list,
 ** and it checks all the neighbors for proximity.  If the closest 
 ** neighbor distance is closer than the threshold you supply times the 
 ** length of the line segment, then the line segments are doubled.  
 ** The doubling maintains the straightness of each line segment,
 ** just cuts it up into segments that are smaller than the smallest 
 ** distance to a neighbor.  
 ** 
 ** For all of the calls, you get back the number of violating vertices 
 ** that were found: ie if there were no proximity errors you get 0.
 ** 
 ** If you feed in 0 for the scale factor, you really get 1.0.
 **
 ** You must supply a VERTEX that is pre-loaded with neighbors and 
 ** their distances by having hull_neighbors called on it!
 ** After that, the problem is trivial but at least this encapsulates it.
 **
 ** If the line segment is too long for its proximity, it gets split in 
 ** half. 
 ** 
 ** scale_thresh is the ratio of the segment length to the 
 ** closest allowed approach:  0.1 requires more VERTEXes than
 ** does 10, consistent with the general trend that tighter
 ** thresholds mean more numerical work.  1.0 is "normal", I think.
 ** 
 */
int fix_proximity(VERTEX *V, NUM scale_thresh) {
  int i;
  NUM dist = -1;
  NUM l;
  NUM diff[3];

  if(!V || !V->next)
    return 0;

  if( V_ISDUMMY(V) || FL_ISDUMMY(V->line) )
    return 0;

  /* Find minimum neighbor-segment distance from the segment, and 
   * add vertices if necessary.
   */

  diff_3d(diff,V->next->x,V->x);
  l = norm_3d(diff);

  if((dist = V->r_cl) * scale_thresh  < l) { /* assignment */
    if(dist==0)
      return 0;

    add_vertex_after(V->line, V,
		     new_vertex(0
				,V->x[0] + diff[0]*0.5
				,V->x[1] + diff[1]*0.5
				,V->x[2] + diff[2]*0.5
				,V->line)
		     );
    V->next->r_cl = V->r_cl;
    return 1;
  }
  return 0;
}

int fluxon_fix_proximity(FLUXON *F, NUM scale_thresh) {
  int ret = 0;
  VERTEX *V = F->start;
  int verbosity = F->fc0->world->verbosity;

  if(verbosity >= 2) printf("fluxon_fix_proximity: fluxon %d\n",F->label);

  while(V && V != F->end) {
    VERTEX *Vnext = V->next;
     if(verbosity >= 3) printf("  vertex %d\n",V->label);
    ret += fix_proximity(V,scale_thresh);
    V=Vnext;
  }
  return ret;
}


/* tree walker hooks for global_fix_proximity... */
static NUM sc_thr;
static long sc_acc;
static long gfp_tramp(FLUXON *fl, int lab, int link, int depth) {
  sc_acc += fluxon_fix_proximity(fl,sc_thr);
  return 0;
}

int global_fix_proximity(WORLD *w, NUM scale_thresh) {
  sc_thr = scale_thresh;
  sc_acc = 0;
  tree_walker(w->lines,fl_lab_of,fl_all_ln_of,gfp_tramp,0);
  return sc_acc;
}
  
  
/**********************************************************************  
 ** fix_curvature
 ** fluxon_fix_curvature
 ** global_fix_curvature
 ** 
 ** You feed in a VEREX, and the curvature of the field line is
 ** compared with the upper and lower curvature threshold that you
 ** pass in (in radians). The vertex is split into 3 vertexes if the
 ** curvature is too high, and the vertexs is deleted if the curvature
 ** is too low. When fix_curvature returns, the VERTEX you passed in
 ** may have been split, so that the new internal angle is less than
 ** your specified threshold everywhere.  If you are stepping along a
 ** fluxon fixing curvatures, you should store V->next BEFORE calling
 ** fix_curvature, or you'll have to step through exactly half of any
 ** daughter VERTEXes that are created, or even crash if the vertex
 ** gets deleted.
 ** 
 ** Set 0 curvature to get 0.05 radian (~3 degrees) for splitting 
 ** and 1 degree for merging).  
 **
 ** If a vertex has less than the low threshold, and is sufficiently
 ** far from other fluxons in the area, then it is unlinked and deleted.
 **
 ** fix_curvature splits single angles into triple angles, but is
 ** recursive:  a vertex with a too-acute angle is split into three
 ** vertices with a 2/3/2 split in angle.  If any given sub-vertex 
 ** is still too acute, then fix_curvature recurses on down the 
 ** line.  The 2/3/2 split is there to prevent vertices coming only
 ** in powers of three -- this way, the series of possible splits is
 ** 1, 3, 7, 9, 17, ...  instead of only powers of 3.
 ** 
 ** To make the fundamental 1->3 split, we blow a particular fluxel
 ** out into a pair of equal sides.  The new external vertex angle 
 ** is 2/7 of the original vertex angle, ensuring the 2/3/2 split.
 **
 ** The proximity-curvature dual check may not be optimal for relaxation; 
 ** the two could stand to be tweaked a bit and merged somehow.
 **/

int fix_curvature(VERTEX *V, NUM curve_thresh_high, NUM curve_thresh_low) {
  int i;
  NUM d1[3];
  NUM d2[3];
  NUM cr[3];
  NUM sincurve,coscurve,curve;
  NUM scale;
  int ret=0;

  if(curve_thresh_high==0)
    curve_thresh_high = 0.05;

  if(!V || !V->next || !V->prev || !V->line) 
    return 0;


  diff_3d(d1,V->prev->x,V->x); /* d1 = segment A */
  diff_3d(d2,V->x,V->next->x); /* d2 = segment B */
  
  scale = 1.0/sqrt(norm2_3d(d1) * norm2_3d(d2));

  cross_3d(cr,d1,d2);
  sincurve = norm_3d(cr) * scale;
  coscurve = inner_3d(d1,d2) * scale;
  curve = ATAN2(sincurve,coscurve);


  if( curve > curve_thresh_high ) {
    VERTEX *Vnew;
    NUM P[3], P0[3],P1[3];

    /* Find location of new prior vertex: 3/4 along the previous vec */
    scale_3d(P0, V->x, 3);
    cp_3d(P1,P0);
    sum_3d(P0, P0, V->prev->x);
    scale_3d(P0, P0, 1.0/4);

    /* Find location of new next vertex: 1/4 along the current vec */
    sum_3d(P1, P1, V->next->x);
    scale_3d(P1, P1, 1.0/4);

    /* Find centroid of new triangle - new current vertex. */
    sum_3d(P,P1,P0);
    sum_3d(P, P, V->x);
    scale_3d(V->x,P,1.0/3);

    /* Make the new prior vertex, and link */
    add_vertex_after(V->line, V->prev,
		     new_vertex(0, P0[0],P0[1],P0[2], V->line));
    add_vertex_after(V->line, V,
		     new_vertex(0, P1[0],P1[1],P1[2], V->line));

    // Force global neighbor checks...
    vertex_clear_neighbors(V->prev);
    vertex_clear_neighbors(V);
    vertex_clear_neighbors(V->next);

    return 2;

  } else {
    /******************************
     * Check for deletion
     */

    if( curve < curve_thresh_low ) {

      /* Check distance to nearest neighbor... */
      NUM dnext, dprev;
      
      dnext = cart_3d(V->x, V->next->x);
      dprev = cart_3d(V->x, V->prev->x);
      if(dnext < V->r_ncl / 2 &&
	 dprev < V->prev->r_ncl / 2) {
      
	/* Check maximum displacement of the fluxel if straightened */
	NUM offset_dist;
	NUM pdist,ndist;
	
	offset_dist = p_ls_dist(V->x,V->prev->x,V->next->x);
	offset_dist *= 1.5;
	if(offset_dist >= 0 && 
	   offset_dist < V->prev->r_ncl && 
	   offset_dist < V->r_cl && 
	   offset_dist < V->next->r_ncl) {
	  if(V->line->fc0->world->verbosity > 3){
	    printf("fix_curvature: unlinking %d from line %d\n",V->label,V->line->label);
	  }
	  delete_vertex(V);
	  return -1;
	}
      }
    }
  }
  return 0;
}


int fluxon_fix_curvature(FLUXON *f, NUM curve_thresh_high, NUM curve_thresh_low) {
  int ret=0; //number of v's different from original
  VERTEX *V = f->start;
  int verbosity = f->fc0->world->verbosity;

  if(verbosity >= 2) printf("fluxon_fix_curvature: fluxon %d\n",f->label);
  
  while(V && V != f->end) {
    VERTEX *Vnext = V->next;
    if(verbosity >= 3) printf("  vertex %d\n",V->label);
    ret += fix_curvature(V,curve_thresh_high, curve_thresh_low);
    V=Vnext;
  }
  return ret;
}

/* tree walker hooks for global_fix_curvature... */
static NUM cu_thrh;
static NUM cu_thrl;
static int cu_acc;
static long cu_tramp(FLUXON *fl, int lab, int link, int depth) {
  cu_acc += fluxon_fix_curvature(fl, cu_thrh, cu_thrl);
  return 0;
}
int global_fix_curvature(WORLD *w, NUM curv_thresh_high, NUM curv_thresh_low) {
  cu_thrh = curv_thresh_high;
  cu_thrl = curv_thresh_low;
  cu_acc = 0;
  tree_walker(w->lines,fl_lab_of,fl_all_ln_of,cu_tramp,0);
  world_update_ends(w);
  return cu_acc;
}

/*********************************************************************
 * Reconnection code - check for a threshold rotation/distance and
 * reconnect if it is exceeded.
 */

/******************************
 * reconnect_vertices: given two vertices on different fluxons, 
 * swap their ->next connections.  No need to slosh vertices -- that is handled
 * adequately in gather_neighbor_candidates, on the next force iteration.
 *
 */
void reconnect_vertices( VERTEX *v1, VERTEX *v2, long passno ) {
  VERTEX *vv;
  FLUXON *f1, *f2;
  FLUX_CONCENTRATION *fc;
  int i,j;

  if(!v1 || !v2 || !v1->next || !v2->next ) {
    fprintf(stderr,"reconnect_vertices: error -- got a null or end vertex!\n");
    return;
  }

  if(v1->passno ==passno || 
     v2->passno == passno || 
     v1->next->passno == passno || 
     v2->next->passno == passno ||
     (v1->prev && v1->prev->passno==passno) ||
     (v2->prev && v2->prev->passno==passno) ||
     (v1->next->next && v1->next->next->passno==passno) ||
     (v2->next->next && v2->next->next->passno==passno)
     ) {
    fprintf(stderr,"reconnect_vertices: tried to reconnect an already-reconnected vertex -- nope.\n");
    return;
  }
  
  if(v1==v2) {
    fprintf(stderr,"reconnect_vertices: error -- tried to reconnect a vertex to itself!\n");
    return;
  }

  f1 = v1->line;
  f2 = v2->line;
  
  if(f1==f2) {
    /*********************************************
     * Self-reconnection: create a new plasmoid 
     */
    FLUXON *Fnew;
    VERTEX *firstv, *lastv, *v;
    WORLD *w;
    long ifirst, ilast,i;
    firstv = 0;
    lastv = 0;

    printf("(self case): Reconnecting (%d; l=%d%s) -- (%d; l=%d%s)\n",
	   v1->label, v1->line->label, v1->line->plasmoid?" P":"",
	   v2->label, v2->line->label, v2->line->plasmoid?" P":""
	   );

    w = f1->fc0->world;

    /* Sort the vertices into (first,last) order */
    for(i=0, v=f1->start; v && v->next && !lastv; i++, v=v->next) {
      if( v==v1 || v==v2 ) {
	if( !firstv ) {
	  ifirst = i;
	  firstv = v;
	} else {
	  ilast = i;
	  lastv = v;
	}
      }
    }

    if(!firstv || !lastv) {
      fprintf(stderr,"Error in self-reconnection code: one or more of the vertices weren't on the fluxon! Ignoring.\n");
      return;
    }

    if(firstv->next==lastv || firstv->next->next==lastv) {
      fprintf(stderr,"Error in plasmoid reconnection code: must be more than two vertices between the two reconnecting vertices! Ignoring (vertices = (%d, l=%d; %d, l=%d).\n",v1->label,v1->line->label,v2->label,v2->line->label);
      return;
    }

    /* Create the new plasmoid and link it into the world */
    Fnew = new_fluxon(f1->flux, 
			       w->fc_pb,
			       w->fc_pe,
			       0,
			       0
		      );
    Fnew->plasmoid = 1;
    

    /* Create the new plasmoid's end vertices and link them into the world */
    Fnew->start = new_vertex(0, 
			     0,0,0,
			     Fnew);

    Fnew->end = new_vertex(0,
			   0,0,0,
			   Fnew);

    
    /* Now cut the plasmoid out of the original fluxon... */
    v = firstv->next;
    firstv->next = lastv->next;
    firstv->next->prev = firstv;
    firstv->line->v_ct -= (ilast - ifirst);
    vertex_clear_neighbors(firstv);
    vertex_clear_neighbors(firstv->next);

    /* ... and link it into the new fluxon */
    Fnew->start->next = v;
    v->prev = Fnew->start;

    vertex_clear_neighbors(v);
    vertex_clear_neighbors(v->prev);
    vertex_clear_neighbors(v->next);

    Fnew->end->prev = lastv;
    lastv->next = Fnew->end;
    
    vertex_clear_neighbors(lastv);
    vertex_clear_neighbors(lastv->next);
    vertex_clear_neighbors(lastv->prev);
    
    Fnew->v_ct = 2 + (ilast - ifirst);
    
    for(i=0, v=Fnew->start; v; i++, v=v->next )
      v->line = Fnew;
    
    if(i != Fnew->v_ct) {
      fprintf(stderr, "Error in plasmoid reconnection case - v_ct is %d, i is %d!\n", Fnew->v_ct, i);
      Fnew->v_ct = i;
    }

    /* Now taint the vertices with the current passno */
    Fnew->start->next->passno = passno;
    firstv->passno = passno;

    return;

  } else if(f1->plasmoid || f2->plasmoid) {
    /*********************************************
     * Non-self reconnection involving a plasmoid: destroys a fluxon
     */
    VERTEX *v, *vpstart;

    if(f1->plasmoid && !f2->plasmoid) {
      // Swap the vertices to ensure that f2 is the plasmoid 
      v=v2; v2=v1; v1=v;
      f1=v1->line; 
      f2=v2->line;
    }

    printf("(plasmoid case): Reconnecting (%d; l=%d%s) -- (%d; l=%d%s)\n",
	   v1->label, v1->line->label, v1->line->plasmoid?" P":"",
	   v2->label, v2->line->label, v2->line->plasmoid?" P":""
	   );

    /******************************
     * Rectify the plasmoid linked list...
     */
    // If it's the very first VERTEX, switch to second-to-last (equivalent in a plasmoid)...
    if( ! v2->prev ) {
      printf("v2 was at the start of the plasmoid: v%d -> v%d\n",v2->label,v2->line->end->prev->label);
      v2 = v2->line->end->prev;
    } else if(!v2->next->next) {
      printf("v2 was at the end of the plasmoid: v%d -> v%d\n",v2->label,v2->line->end->prev->label);
      v2 = v2->line->start->next;
    }

    // splice a genuine loop out of the plasmoid
    vertex_clear_neighbors(f2->end->prev);
    vertex_clear_neighbors(f2->start->next);

    f2->end->prev->next   = f2->start->next;
    f2->start->next->prev = f2->end->prev;
    f2->start->next       = f2->end;
    f2->end->prev         = f2->start;

    // Cut the genuine loop and splice it into the other fluxon
    v=v1->next;
    v1->next = v2->next;
    v1->next->prev = v1;
    v2->next = v;
    v2->next->prev = v2;

    vertex_clear_neighbors(v1);
    vertex_clear_neighbors(v2);
    
    // Set the fluxon pointer of all vertices
    for(v=f1->start; v; v=v->next)
      v->line = f1;

    // Fix up the vertex counters of both fluxons
    f1->v_ct += f2->v_ct - 2;
    f2->start->next = f2->end;
    f2->end->prev = f2->start;
    f2->v_ct = 2;

    // Taint the reconnected vertices
    v1->passno = passno;
    v2->passno = passno;

    // Now delete the trivial plasmoid
    delete_fluxon(f2);
     
  } else {
    /*********************************************
     * Non-self-reconnection - normal case 
     */

    printf("(normal case): Reconnecting (%d; l=%d%s) -- (%d; l=%d%s)\n",
	   v1->label, v1->line->label, v1->line->plasmoid?" P":"",
	   v2->label, v2->line->label, v2->line->plasmoid?" P":""
	   );
    vv = v1->next;
    v1->next = v2->next;
    v1->next->prev = v1;
    
    v2->next = vv;
    v2->next->prev = v2;

    vertex_clear_neighbors(v1);
    vertex_clear_neighbors(v2);
    vertex_clear_neighbors(v1->next);
    vertex_clear_neighbors(v2->next);
    if(v1->prev)
      vertex_clear_neighbors(v1->prev);
    if(v2->prev)
      vertex_clear_neighbors(v2->prev);
    
    /* Now clean up the fluxons and the back-to-fluxon links in the individual vertices. */
    
    /* f1 */
    i=1;
    for( vv = f1->start; vv && vv->next; vv=vv->next ) {
      i++;
      vv->line = f1;
    }
    f1->end = vv;
    f1->end->line = f1;
    
    /* f2 */
    j=1;
    for( vv = f2->start; vv && vv->next; vv=vv->next ) {
      j++;
      vv->line = f2;
    }
    f2->end = vv;
    f2->end->line = f2;
    
    if(i+j != f1->v_ct + f2->v_ct) {
      fprintf(stderr,"Hmmm -- something's funny with the vertex count.  Reconnected %d(fl %d) to %d (fl %d), final total vertex count is %d, previously %d\n",v1->label,f1->label,v2->label,f2->label,i+j,f1->v_ct+f2->v_ct);
    }
    f1->v_ct = i;
    f2->v_ct = j;
    
    /* Switch the end flux concentrations... */
    f1->fc1->lines = tree_unlink(f1, fl_lab_of, fl_end_ln_of);
    f2->fc1->lines = tree_unlink(f2, fl_lab_of, fl_end_ln_of);
    
    fc = f1->fc1;
    f1->fc1 = f2->fc1;
    f2->fc1 = fc;
    
    f1->fc1->lines = tree_binsert(f1->fc1->lines, f1, fl_lab_of, fl_end_ln_of);
    f2->fc1->lines = tree_binsert(f2->fc1->lines, f2, fl_lab_of, fl_end_ln_of);

    /* Taint the reconnected vertices */
    v1->passno = passno;
    v2->passno = passno;
  }
}
  
/******************************
 * v_recon_check: checks the current reconnection condition between a 
 * VERTEX and each of its neighbors.  The *following* segment is 
 * (possibly) reconnected. Reconnection with an image charge
 * is not allowed!  The reconnection conditions are found via the rc_funcs table
 * in the world. 
 *
 * Cleanliness is achieved via the passno tainting mechanism, so you don't
 * get multiple reconnections of the same segment.
 */
int vertex_recon_check( VERTEX *v1, long passno ) {
  WORLD *w = v1->line->fc0->world;
  RC_FUNC **rcfuncp;
  VERTEX *rv = 0;
  int i;

  if(!v1)
    return 0;

  if(v1->passno == passno)
    return 0;
  
  for( i=0, rcfuncp = w->rc_funcs; 
       !rv && i<N_RECON_FUNCS && *rcfuncp ; 
       i++, rcfuncp++ ) 
    rv = (**rcfuncp)(v1,w->rc_params[i]);

  if(rv) {
    if( w->verbosity) {
      printf("Reconnecting vertices %d (on %d at [%.3g,%.3g,%.3g]) and %d (on %d at [%.3g,%.3g,%.3g]): satisfied condition %d\n",
	     v1->label,  v1->line->label,    v1->x[0], v1->x[1], v1->x[2],
	     rv->label,  rv->line->label,    rv->x[0], rv->x[1], rv->x[2],
	     i-1
	     );
    }
    if(rv->passno == passno)
      return 0;
    reconnect_vertices(v1,rv, passno);
    return 1;
  }
  return 0;
}

/******************************
 * fluxon_recon_check
 * Iterates over all the vertices in a fluxon.
 * 
 * DEPRECATED!  The fluxon can disappear in midstream.  The
 * only way to avoid that is to fall back to the vertex level.
 *
 */
long fluxon_recon_check( FLUXON *f, long passno ) {
  VERTEX *v;
  int retval = 0;
  if( !(f->start) || !(f->start->next) || !(f->start->next->next) ) 
    return retval;

  for(v=f->start->next; v->next && v->next->next; v=v->next) {
    retval += vertex_recon_check(v, passno);
    if(retval)
      return retval;
  }
  return retval;
}

static DUMBLIST *grc_vertexlist = 0;
static long grc_tramp(VERTEX *v, int lab, int link, int depth) {
  if( v->next && (v->prev || !(v->line->plasmoid)) )
    grc_vertexlist->stuff[ (grc_vertexlist->n)++ ] = (void *)v;
  return 0;
}

long global_recon_check(WORLD *w) {
  long recon_ct = 0;
  long passno = ++(w->passno);
  long i;

  if(!grc_vertexlist) {
    grc_vertexlist = new_dumblist();
  }
  if(grc_vertexlist->size < w->vertices->world_links.n) {
    dumblist_grow(grc_vertexlist, w->vertices->world_links.n);
  }
  grc_vertexlist->n = 0;

  // Accumulate a DUMBLIST of all existing vertices that aren't 
  // plasmoid starts or fluxon ends.  (No other vertex can be 
  // destroyed by the reconnection process...)

  tree_walker( w->vertices, v_lab_of, v_ln_of, grc_tramp, 0);

  printf("global_recon_check: found %d vertices to check...\n",grc_vertexlist->n);

  for(i=0;i<grc_vertexlist->n;i++) 
    recon_ct += vertex_recon_check((VERTEX *)(grc_vertexlist->stuff[i]) , passno);

  return recon_ct;
}


/******************************
 * fc_cancel - cancel two flux concentrations.  
 * The two flux concentrations must have opposite fluxes (one source, one sink).
 * Fluxons that connect one concentration to the other are deleted.  If, at the 
 * end of the deletion, there are any fluxons left in the weaker of the two, then
 * fluxons are reconnected (in tree [random] order) to connect the two, and then 
 * deleted.  The process stops when one or both concentrations have all their 
 * flux exhausted.
 *
 * Finally, if either of the concentrations has zero flux remaining, it is deleted.
 * 
 * The first flux concentration must be a source, the second a sink.
 *
 * Returns 0 on success, or nonzero on failure.
 */

/** FIXME: Will need to be updated to handle variable-flux fluxons... **/
int fc_cancel(FLUX_CONCENTRATION *fc0, FLUX_CONCENTRATION *fc1) {
  long n_min;
  WORLD *w = fc0->world;
  int v = w->verbosity;


  if(v) printf("fc_cancel: %d and %d\n",fc0->label,fc1->label);
  
  // Handle the trivial cases.
  if(!fc0->lines && !fc1->lines) {
    if(v) printf("fc_cancel: trivial case (both sides).\n");
    delete_flux_concentration( fc0 );
    delete_flux_concentration( fc1 );
    return 0;
  }

  if(!fc0->lines) {
    if(v) printf("fc_cancel: trivial case (fc0).\n");
    delete_flux_concentration(fc0);
    return 0;
  } 

  if(!fc1->lines) {
    if(v) printf("fc_cancel: trivial case (fc1).\n");
    delete_flux_concentration(fc1);
    return 0;
  }

  if(fc0->lines->fc0 != fc0) {
    fprintf(stderr,"cancel: fc0 is not the source of its root fluxon!");
    return 1;
  }

  if(fc1->lines->fc1 != fc1) {
    fprintf(stderr,"cancel: fc1 is not the sink of its root fluxon!");
    return 2;
  }


  n_min = fc0->lines->start_links.n;
  if(fc1->lines->end_links.n < n_min)
    n_min = fc1->lines->end_links.n;

  if(v) printf("fc_cancel: fc0.n is %d; fc1.n is %d; n_min is %d\n",
	       fc0->lines->start_links.n,
	       fc1->lines->end_links.n,
	       n_min
	       );

  while( n_min  ) {
    printf("%d: ",n_min);

    if(!fc0->lines || !fc1->lines) {
      fprintf(stderr, "cancel: accounting error! (n_min=%d, fc0 has %d, fc1 has %d)\n",n_min, (fc0->lines? fc0->lines->start_links.n : 0), (fc1->lines ? fc1->lines->end_links.n : 0) );
      return 3;
    }
    
    if( fc0->lines->fc1 == fc1 ) {
      if(v) printf("fc_cancel: fluxon %d (top of fc0) goes fc0->fc1\n",fc0->lines->label);
      delete_fluxon( fc0->lines );
      goto cancel_loop_end;
    } 

    if( fc1->lines->fc0 == fc0 ) {
      if(v) printf("fc_cancel: fluxon %d (top of fc1) goes fc0->fc1\n",fc1->lines->label);
      delete_fluxon( fc1->lines );
      goto cancel_loop_end;
    } 

    if(v) printf("No fc0->fc1 found. Reconnecting %d and %d...\n",fc0->lines->label,fc1->lines->label);
    reconnect_vertices( fc0->lines->start, fc1->lines->end->prev, ++(fc0->world->passno));
    if(fc0->lines->fc1 != fc1) {
      fprintf(stderr, "cancel: reconnection didn't work out like we hoped for fluxon %d...\n",fc0->lines->label);
      return 4;
    }
    delete_fluxon(fc0->lines);

  cancel_loop_end: n_min--;
  }

  if(!fc0->lines) {
    printf("Deleting fc0 (%d) ",fc0->label);
    delete_flux_concentration(fc0);
  }

  if(!fc1->lines) {
    printf("Deleting fc1 (%d) ",fc1->label);
    delete_flux_concentration(fc1);
  }

  return 0;
}
    
/********************************************************************
 * fluxon end-condition handlers - update end position of the fluxon
 * to be consistent with the boundary condition.
 * 
 * the utility routine funky_fl_b indicates whether it's OK for the final
 * segment of the fluxons on this flux concentration to end at the exact
 * same point in space.  If you add a new boundary condition that 
 * allows swivel-hook action, you need to update funky_ends.
 */
int funky_fl_b(FLUX_CONCENTRATION *fc) {
  void *bound;
  if(!fc)
    return 0;
  bound = fc->bound;
  if(!bound)
    bound = fc->world->default_bound;
  return (bound==fl_b_tied_inject || bound==fl_b_tied_force);
}

struct F_B_NAMES F_B_NAMES[] = {
  {fl_b_tied_inject, "fl_b_tied_inject", "Auto-inject new vertices to maintain one near the FC"},
  {fl_b_tied_force,  "fl_b_tied_force",  "Force the closest vertex to be near the FC"},
  {fl_b_open,        "fl_b_open",        "Enforce open (source surface) condition" },
  {fl_b_plasmoid,    "fl_b_plasmoid",    "Enforce plasmoid (ourobouros) condition" }, 
  {0, 0},
  {0, 0}
};
void *boundary_name_to_ptr(char *s) {
  int i;
  if(!s) 
    return 0;
  for(i=0; F_B_NAMES[i].func; i++) {
    if(!strcmp(s, F_B_NAMES[i].name)){
      return F_B_NAMES[i].func;
    }
  }
  return 0;
}

char *boundary_ptr_to_name(void *f) {
  int i;
  if(!f)
    return 0;
  for(i=0; F_B_NAMES[i].func; i++) {
    if(F_B_NAMES[i].func == f) {
      return F_B_NAMES[i].name;
    }
  }
}



void fl_b_tied_inject(VERTEX *v) {
  POINT3D a;
  NUM lr;
  NUM r;

  if(!v->prev) {
    if(!v->next) {
      fprintf(stderr,"HEY! fl_b_tied_inject got a loner vertex! Ignoring....\n");
      return;
    }

    // Handle start-of-fluxon case

    lr = v->line->fc0->locale_radius;
    if(lr>0) {
      diff_3d(a, v->next->x, v->x);
      r = norm_3d(a);
      
      if(r > lr) {
	VERTEX *vn;
	long i;
       
	// Inject a new vertex, colinear, halfway between the offending vertex and the start
	sum_3d(a, v->next->x, v->x);
	scale_3d(a, a, 0.5);
	vn = new_vertex(0, a[0], a[1], a[2], v->line);
	add_vertex_after(v->line, v, vn);

	vn->r_cl = v->r_cl;
	
	// Seed the new vertex with some neighbors...
	for(i=0;i<v->neighbors.n;i++) {
	  dumblist_add( &(vn->neighbors), v->neighbors.stuff[i] );
	  dumblist_add( &(((VERTEX *)(v->neighbors.stuff[i]))->nearby), vn );
	}

	// ... and stick it in the neighbors' candidate lists...
	for(i=0; i<v->nearby.n; i++) {
	  dumblist_add( &(((VERTEX *)(v->nearby.stuff[i]))->neighbors), vn); // Add vn as a neighbor to one of v's nearby vertices
	  dumblist_add( &(vn->nearby), v->nearby.stuff[i] );                 // Add v's nearby vertex as a neighbor candidate.

	}

	// Now do an extra neighbor calculation to make sure everyone's copascetic
	vertex_update_neighbors(vn,0);
	for(i=0;i<vn->neighbors.n;i++) {
	  vertex_update_neighbors( ((VERTEX **)(vn->neighbors.stuff))[i], 0);
	}

	// ... Finally, rig it up to do nothing this time step...
	cp_3d(vn->f_s, v->f_s);
	vn->f_s_tot = v->f_s_tot;
	vn->f_v[0] = vn->f_v[1] = vn->f_v[2] = 0;
	vn->f_v_ps[0] =vn->f_v_ps[1] = vn->f_v_ps[2];
	vn->f_v_tot = 1.0;
	cp_3d(vn->f_t, v->f_t);
	vn->r_v = v->r_v;
	vn->r_s = v->r_s;
	vn->b_mag = v->b_mag;
	cp_3d(vn->b_vec, v->b_vec);
	

      }
    }
    return;
  } 


  else if( !v->next ) {
    
    // Handle end-of-fluxon case

    lr = v->line->fc0->locale_radius; 
    if(lr>0) {
      diff_3d(a, v->prev->x, v->x);
      r = norm_3d(a);

      if(r > lr) {
	VERTEX *vn;
	long i;
	
	// Inject a new vertex, colinear, halfway between the offending vertex and the end
	sum_3d(a, v->prev->x, v->x);
	scale_3d(a, a, 0.5);
	vn = new_vertex(0, a[0], a[1], a[2], v->line);
	add_vertex_after(v->line, v->prev, vn);
	
	vn->r_cl = v->prev->r_cl;
	
	// Seed the new vertex with some neighbors...
	for(i=0;i<v->prev->neighbors.n;i++) {
	  dumblist_add( &(vn->neighbors), v->prev->neighbors.stuff[i] );
	  dumblist_add( &(((VERTEX *)(v->prev->neighbors.stuff[i]))->nearby), vn );
	}

	// ... and stick it in the neighbors' candidate lists...
	for(i=0; i<v->prev->nearby.n; i++) {
	  dumblist_add( &(vn->nearby), v->prev->nearby.stuff[i] );
	  dumblist_add( &(((VERTEX *)(v->prev->nearby.stuff[i]))->neighbors), vn);
	}

	// ... Finally, rig it up to do nothing this time step...
	cp_3d(vn->f_s, v->prev->f_s);
	vn->f_s_tot = v->prev->f_s_tot;
	vn->f_v[0] = vn->f_v[1] = vn->f_v[2] = 0;
	vn->f_v_tot = 1.0;
	cp_3d(vn->f_t, v->prev->f_t);
	vn->r_v = v->prev->r_v;
	vn->r_s = v->prev->r_s;
	vn->b_mag = v->prev->b_mag;
	cp_3d(vn->b_vec, v->prev->b_vec);
      }
    }
    return;
  }

  else {
    fprintf(stderr," fl_b_tied_inject: got a middle vertex; ignoring...\n");
    return;
  }

}
	
void fl_b_tied_force(VERTEX *v) {
  POINT3D a;
  NUM lr;
  NUM r;

  if(!v->prev) {
    if(!v->next) {
      fprintf(stderr,"HEY! fl_b_tied got a loner vertex! Ignoring...\n");
      return;
    }

    // Handle start-of-fluxon case

    lr = v->line->fc0->locale_radius;
    if( lr > 0 ) {
      diff_3d(a, v->next->x, v->x);
      r = norm_3d(a);
      
      if( r > lr ) {
	scale_3d(a, a, lr / r);
	sum_3d(v->next->x, a, v->x);
      }
    }
    return;

  }

  else if(!v->next) {
    
    // handle end-of-fluxon case 

    lr = v->line->fc0->locale_radius;
    if(lr > 0) {
      diff_3d(a, v->prev->x, v->x);
      r = norm_3d(a);
      
      if( r > lr ) {
	scale_3d( a, a, lr/r );
	sum_3d(v->prev->x, a, v->x);
      }
    }
    return;
  }


  else {
    fprintf(stderr,"HEY! fl_b_tied got a middle vertex! Ignoring...\n");
    return;
  }

}


void fl_b_open(VERTEX *v) {
  POINT3D a;
  if(!v->prev) {

    /** Check for loner vertex **/
    if(!v->next) {
      fprintf(stderr,"HEY! fl_b_tied got a loner vertex! Ignoring...\n");
      return;
    }
    
    /** Start vertex **/
    if(v->line->fc0->locale_radius > 0) {
      diff_3d( a,    v->next->x, v->line->fc0->x );
      scale_3d(a,    a,          v->line->fc0->locale_radius / norm_3d(a) );
      sum_3d(  v->x, a,          v->line->fc0->x );
    }
    return;

  }

  /** End vertex **/
  if(!v->next) {

    if(v->line->fc1->locale_radius > 0) {
      diff_3d( a,     v->prev->x,   v->line->fc1->x );
      scale_3d(a,     a,            v->line->fc1->locale_radius / norm_3d(a) );
      sum_3d(  v->x,  a,            v->line->fc1->x );
    }
    return;

  }

  /** Shouldn't get here (normally) **/
  fprintf(stderr,"HEY! fl_b_open got a middle vertex! Doing nothing...\n");
}    

void fl_b_plasmoid(VERTEX *v) {
  if(!v->prev) {

    /** Check for loner vertex **/
    if(!v->next) {
      fprintf(stderr,"HEY! fl_b_plasmoid got a loner vertex! Ignoring....\n");
      return;
    }
    
    /** Start vertex **/

    if(!v->next->next)
      return;
    
    cp_3d( v->x, v->line->end->prev->x );
    return;

  }

  if(!v->next) {
    
    /** End vertex **/

    if(!v->prev->prev)
      return;
    
    cp_3d( v->x, v->line->start->next->x );
    return;
  }

  fprintf(stderr,"HEY! fl_b_plasmoid got a middle vertex! Doing nothing...\n");
  
}

/**********************************************************************
 * First-stage parallelization stuff...
 *
 * First-stage parallelization of world_relax_step and
 * world_update_mag.  Each one Forks off the detailed calculation in
 * roughly vertex-equal chunks.  They work by keeping a DUMBLIST of
 * the fluxons to be processed; when enough are accumulated, each forks
 * off a processor and starts the DUMBLIST over again.  Outstanding
 * processes are handled via another DUMBLIST, and their results are
 * processed via the binary dump mechanism.
 *
 * 
 * The parallel_prep and parallel_finish set up initial operations to get the 
 * various globals set right.  
 * 
 */

/******************************
 * Globals for parallization bookkeeping
 */
static long v_ct;            // Number of vertices in current batch
static long v_thresh;        // Threshold vertex count above which to spawn a worker process
static FLUXON *last_fluxon;  // Gets the last fluxon in the tree (set in v_ct_springboard)

static DUMBLIST *fluxon_batch = 0;  // Dumblist gets the current fluxon batch to be processed.

typedef struct SUBPROC_DESC {   // Keeps track of what we spawned...
  long pid;
  long pipe;
  DUMBLIST batch;
} SUBPROC_DESC;
static SUBPROC_DESC *sbd, *sbdi;    // Used for bookkeeping subprocs.

int (*work_springboard)(FLUXON *f); // Contains the actual action 


/******************************
 * v_ct_springboard - counts up the total number of vertices in the sim.
 */
static long v_ct_springboard(FLUXON *fl, int lab, int link, int depth) {
  if( ! FL_ISDUMMY(fl) )
    v_ct += fl->v_ct;
  last_fluxon = fl;
  return 0;
}


/******************************
 * parallel_prep - call this to set up a parallelization run.
 * Counts vertices and sets the v_thresh, and zeroes out the 
 * subproc dumblist...
 */
void parallel_prep(WORLD *a) {
  int i;

  /* Make sure our dumblists exist */
  if(!fluxon_batch)
    fluxon_batch = new_dumblist();

  if(a->concurrency > 1000) {
    fprintf(stderr,"parallel_prep: WARNING - concurrency is ludicrous! (%d)\n",a->concurrency);
  }

  sbdi = sbd = (SUBPROC_DESC *)localmalloc(sizeof(SUBPROC_DESC) * a->concurrency, MALLOC_MISC);
  
  if(!sbd) {
    fprintf(stderr,"parallel_prep: sbd is 0 (malloc failed).  I'm about to crash....\n");
  }

  // Zero out the pointers and such in the dumblists, so as not to confuse the library.
  for(i=0;i<a->concurrency; i++) 
    dumblist_init(&(sbd[i].batch));
 
  /* Accumulate the total number of vertices, to divide 'em up about equally */
  v_ct = 0;
  tree_walker(a->lines, fl_lab_of, fl_all_ln_of, v_ct_springboard,0);


  v_thresh = v_ct / a->concurrency;
  //  printf("v_ct is %d, v_thresh is %d\n",v_ct,v_thresh);
  fflush(stdout);
  
  dumblist_clear(fluxon_batch);
  v_ct = 0; // Zero the vertex count for comparison against v_thresh...
}


/******************************
 * parallel_fluxon_spawn_springboard - call this from tree_walker 
 * to launch your task, parallelized by fluxon.  You determine the
 * task by setting the global work_springboard, which is a 
 * pointer to a function that accepts a FLUXON * and returns an
 * int.  
 */

static int parallel_daughter(int p); // handles processing in the daughters.

static long parallel_fluxon_spawn_springboard(FLUXON *fl, int lab, int link, int depth) {
  /* Skip magic boundary fluxons */
  if( FL_ISDUMMY(fl) ) 
    return 0;

  dumblist_add( fluxon_batch, fl );
  v_ct += fl->v_ct;

  if( v_ct > v_thresh+1 || fl==last_fluxon ) { // +1 ensures we run over, not under; last_fluxon check ensures we get the last one.
    int p[2], pid;

    if(pipe(p)) {
      fprintf(stderr,"Pipe creation failed! Giving up!\n");
      return 1;
    }
    pid = fork();

    if(pid) {

      /** parent process **/
      close(p[1]);

      if(fl->fc0->world->verbosity) {
	printf("spawned %d at %d (thresh %d)...",pid,v_ct,v_thresh);
	fflush(stdout);
      }
      
      /* Clear the batch list */
      v_ct = 0;
      dumblist_clear( fluxon_batch );
      
      /* Store the pid information including re-do batch list in case we have to, er, redo.
       */
      sbdi->pid = pid;
      sbdi->pipe = p[0];
      dumblist_snarf(&(sbdi->batch), fluxon_batch);

      sbdi++;
      
      return 0;

    } else {

      /*** daughter process ***/
      close(p[0]);
      exit(parallel_daughter(p[1]));

      /*** Never executes ***/
      fprintf(stderr,"You should never see this!\n");
      return 1;
    }
  }

  /* Normal exit case - just iterate again and spawn another time */
  return 0;

}

/******************************
 * parallel_daughter 
 * Handle springboard processing in the daughter processes.
 * If the pipe argument contains an actual fd, then 
 * it dumps the fluxon results to the pipe (which carries them
 * back to the parent process.
 */
int parallel_daughter(int p) {
  int i;
  int pid = getpid();
  
  for(i=0; i<fluxon_batch->n; i++) {
    FLUXON *f = ((FLUXON **)(fluxon_batch->stuff))[i];
    if(f->fc0->world->verbosity) {
      printf("pid %d: processing fluxon %d (%d of %d)\n",pid,f->label,i,fluxon_batch->n);
      fflush(stdout);
    }
    if((*work_springboard)( f )) {
      printf("pid %d: fluxon %d failed its springboard! Exiting early.\n",pid,f->label);
      fflush(stdout);
      return(1);
    }
  }
  
  if(p) {
    for(i=0; i<fluxon_batch->n; i++) {
      FLUXON *f = ((FLUXON **)(fluxon_batch->stuff))[i];
      binary_dump_fluxon_pipe( p, f );
    }
    binary_dump_end( p );
  }
  return(0);
}

/******************************
 * parallel_finish - call this to suck back in all the state
 * from the daughter subprocesses.  If there's a failure (e.g. an
 * interrupted system call) then re-try the daughter call.
 * 
 * Returns 0 on normal completion, 1 on error.
 */
int parallel_finish(WORLD *a) {
  SUBPROC_DESC *sbdii;
  int throw_err = 0;
  int i;
  init_minmax_accumulator(a);

  for(sbdii=sbd; sbdii < sbdi; sbdii++) {
    WORLD *brd_ret;
    if(a->verbosity) {
      printf("Reading dumpfile for pid %d...",sbdii->pid);
      fflush(stdout);
    }
    brd_ret = binary_read_dumpfile( sbdii->pipe, a );
    close(sbdii->pipe);

    // Uh oh -- we failed!  Re-try the daughter process task. 
    // Since this shouldn't happen too often, we just do it here 
    // in the parent rather than trying to respawn the daughter.

    if(!brd_ret) {
      dumblist_clear(fluxon_batch);
      dumblist_snarf(fluxon_batch, &(sbdii->batch));

      fprintf(stderr,"parallel_finish: Oh noes! Daughter failed for fluxons");
      for(i=0; i < fluxon_batch->n; i++) {
	fprintf(stderr," %d",( ((FLUXON **)(fluxon_batch->stuff))[i] )->label );
      }
      fprintf(stderr,"!\n     This is sometimes caused by a system interrupt damaging the pipe. Re-trying in the parent...\n");
      
      i = parallel_daughter(0); // sending 0 omits dumping
      if(i) {
	fprintf(stderr,"     task failed on re-try.  I will snarf all remaining data and return an error.\n");
	throw_err++;
      } else {
	fprintf(stderr,"     OK.  Continuing...\n");
      }
    }

    if(a->verbosity) {
      printf("finished pid %d dump..\n",sbdii->pid);
    }
  }

  finalize_minmax_accumulator(a);

  /* Clean up the mess of the daughters... */
  {
    int status;
    int pid;
    while(  (pid = wait(&status)) > 0) { // assignment
      if(status) {
	printf("Whoa! pid %d terminated with status %d (%s)\n",pid,status,strerror(status));
      }
    }
  }


  // Clear the batch dumblist but leave around for next time.
  dumblist_clear(fluxon_batch);

  // Clean up the cached batches and free them (a waste - we really
  // ought to keep 'em around for next time...)
  for(sbdii=sbd; sbdii<sbdi; sbdii++) 
    dumblist_clean(&(sbdii->batch));
  localfree(sbd);

  return throw_err;
}


/**********************************************************************
 *
 * Parallelized arena operations follow
 *
 */

/******************************
 * world_relax_step_parallel
 */

static int fluxon_calc_step_sb(FLUXON *f) {
  fluxon_calc_step(f, gl_t);
  return 0;
}

void world_relax_step_parallel(WORLD *a, NUM t) {
  gl_t = t;

  /* If no concurrency, don't bother doing it in parallel */
  if(a->concurrency < 1) {
    world_relax_step(a, t);
    return;
  }

  parallel_prep(a);                       // Get ready

  gl_t = t;                               // Set up global dtau variable for springboarding                                 
  work_springboard = fluxon_calc_step_sb; // This is the parallelized operation
  tree_walker(a->lines, fl_lab_of, fl_all_ln_of, parallel_fluxon_spawn_springboard, 0); 

  parallel_finish(a);     // Snarf up our state again

  /* Now advance the thing (single threaded just yet -- fix?).... */
  gl_t = t;
  tree_walker(a->lines, fl_lab_of, fl_all_ln_of, w_r_s_springboard, 0);

}


/******************************
 * world_update_mag_parallel
 */

static int fluxon_update_mag_sb(FLUXON *f) {
  return fluxon_update_mag(f, gl_gl, gl_f_funcs);
}

int world_update_mag_parallel(WORLD *a, char global) {


  if(a->concurrency < 1) {
    printf("Warning - using single threaded mag...\n");
    return world_update_mag(a,global);
  }

  gl_a = a;
  gl_gl = global;
  gl_f_funcs = a->f_funcs;
  a->max_angle = 1.0;
  a->mean_angle = 0.0;

  parallel_prep(a);

  work_springboard = fluxon_update_mag_sb;
  tree_walker(a->lines, fl_lab_of, fl_all_ln_of, parallel_fluxon_spawn_springboard, 0);
  
  parallel_finish(a);
  return 0;
}




















/**********************************************************************
 * Photosphere hull routines
 * 
 * Tese routines are to find the hull of the begin/end vertices on the
 * photosphere. They are only called through PDL (p
 * $w->vertex(#)->photohull), and none of the information is stored in
 * the vertex itself.
 */

/**********************************************************************
 * gather_photosphere_neighbor_candidates
 * 
 * This finds the neighbor candidates of a first or last vertex on a
 * planar photosphere. The list of neighnor candidates are all of the
 * other first/last vertices in the world. This can be slow. This
 * routine was developed for determining the exact hull (and hence
 * flux) on a photosphere for a better energy calculation.
 * 
 * It returns a pointer to its own workspace -- so calling it kills
 * your last instance of its return value, unless you've gone and copied
 * it yourself.
 * 
 * Use with other routines that have 'photosphere' in the name. 
 * 
 * Created by Laurel Rachmeler: Jan 2009.
 *
 */



/*********************************************************************
 * Helper routine for snarfing up first and last vertex of all fluxons
 * -- called via tree_walk, in data.c. Does not collect v because v is
 * the only vertex that has the updated passno, and only does real
 * fluxons, no images. It also checkes to make sure that only tied
 * vertices are added, not those that are on an open boundary.
 *
 * NOTE: it does not check which photosphere the vertex is on. We only
 * tie to one photopshere now even though we can do it to 2. you may
 * need to modify this somehow to check for that later.
*/
static DUMBLIST *snarfer_workspace;
static long begin_end_snarfer(FLUXON *f, int lab_of, int ln_of, long depth) {
  long passno = f->fc0->world->passno;
  if (f->label > 0) {
    if (funky_fl_b(f->fc0) && f->start->passno != passno)
      dumblist_quickadd(snarfer_workspace, f->start);
    if (funky_fl_b(f->fc1) && f->end->passno != passno)  
      dumblist_quickadd(snarfer_workspace, f->end);
  }
  return 0;
}


DUMBLIST *gather_photosphere_neighbor_candidates(VERTEX *v, char global){
  static DUMBLIST *workspace =0;
  void **foo;
  int i;
  int n;
  int verbosity = v->line->fc0->world->verbosity;
  long passno = ++(v->line->fc0->world->passno);
  v->passno = passno; //identifier for current vertex
  
  if(verbosity >= 2) {
    printf("passno=%d...  ",passno);
  }

  if(!v)        /* Paranoia */
    return;

  if(v != v->line->start && v != v->line->end){
    printf("Oops, gather_photosphere_neighbor_candidates got a non begin/end vertex (label %d) \n",v->label);
    return;
  }

  if(!workspace) 
    workspace = new_dumblist();

  workspace->n = 0;


   /* find all the end vertices world by doing a tree-walk trhough
      each line and grabbing the first/last vertex of each.  */

    FLUXON *f;

    if(verbosity >= 3)
      printf("searching for neighbors in gather_photosphere_neighbors.");

    f = v->line;

    while(f->all_links.up)
      f=f->all_links.up; //goes to the top of the tree

    snarfer_workspace = workspace;
    // What is the reason for the new pointer name? Just to make begin_end_snarfer as general as possible


    /* calls a tree (fluxons) and walks along all of the tree. uses
       begin_end_snarfer to add each begin and end vertex into
       snarfer_workspace.*/
    tree_walk(f, fl_lab_of, fl_all_ln_of, begin_end_snarfer);

    if(verbosity >= 4)
          {
            int i;
            printf("gather_photosphere_neighbor_candidates:  global neighbor add gives:\n\t");
            for(i=0;i<workspace->n;i++)
      	printf("%d ",((VERTEX *)(workspace->stuff[i]))->label);
            printf("\n");
          }

  if(verbosity >= 3) printf("gather_photosphere_neighbor_candidates: returning %d elements\n",workspace->n);

  return workspace;
}



/**********************************************************************
 * photosphere_hull_neighbors - Out of a list of neighbor candidates,
 * photosphere_hull_neighbors generates a true neighbor list.  It also
 * produces planar hull information and a host of other useful 2-D
 * projected goodies.
 *
 * Each of the neighbor candidate segments is projected onto the
 * photosphere's plane (photosphere MUST be a plane).  Then the
 * Voronoi cell is constructed around the central point in the
 * projection plane.  Fortunately, routines in geometry.c take care of
 * the geometric details.
 *
 * As long as a hull calculation is being done anyhow,
 * photosphere_hull_neigbors keeps the valuable hull vertex
 * information and returns it in a buffer of (x,y) pairs.  The buffer
 * is ephemeral -- it gets overwritten each time
 * photosphere_hull_neighbors or hull_neighbors is called -- so use it
 * while you can.
 */

HULL_VERTEX *photosphere_hull_neighbors(VERTEX *v, DUMBLIST *horde) {
  int i;
  int a,b;
  POINT3D x0;
  NUM pm[9];
  NUM or[3]={0,0,0};//how do i initialize this array? (syntax)
  NUM temp_p[3];
  int verbosity = v->line->fc0->world->verbosity;

  static HULL_VERTEX *voronoi_buf = 0;
  static int voronoi_bufsiz = 0;

  if(v->line->fc0->world->verbosity >= 5) 
    printf("Entering photosphere hull_neighbors...\n");

  if(v != v->line->start && v != v->line->end){
    printf("Oops, photosphere_hull_neighbors got a non begin/end vertex (label %d) \n",v->label);
    return 0;
  }

  /* Project the horde into the photospheric plane The projected
   * vectors go into the vertices' scratch space. In the case of a
   * horizonal plane, the effect should be to remove the z-coordinate.
   */
  if(v->line->fc0->world->verbosity>=5)
    printf("hull_neighbors calling project_n_fill_photosphere...\n");

  project_n_fill_photosphere(v, horde); /* in geometry.c, not
					   nonlinear, only used vertex
					   positions, not fluxel
					   paths */

  if(v->line->fc0->world->verbosity>=5)
    printf("hull_neighbors_photosphere back from project_n_fill_photosphere...\n");


  /* Grow the buffer if necessary. */
  if(voronoi_bufsiz <= horde->n*2) {
    voronoi_bufsiz = horde->n*4;

    if(voronoi_buf)
      localfree(voronoi_buf);
    voronoi_buf = (HULL_VERTEX *)localmalloc((voronoi_bufsiz)*sizeof(HULL_VERTEX),MALLOC_VL);

    if(!voronoi_buf) {
      fprintf(stderr,"Couldn't get memory in hull_neighbors!\n");
      exit(99);
    } 

  }

  if(verbosity > 1) {
    printf(".");
    
    if(verbosity >= 5)
      printf("V%4d: (%7.3g, %7.3g, %7.3g) -- (%7.3g, %7.3g, %7.3g)\n",v->label,v->x[0],v->x[1],v->x[2], v->next->x[0],v->next->x[1],v->next->x[2]);
  }

  /* Find the 2-D hull.  Don't want rejects. */
  /*  hull_2d(HULL_VERTEX out,VERTEX list for neighbors,central v),
      are they the same?; */
  hull_2d_us(voronoi_buf, horde, v);

  projmatrix(pm,or,v->line->fc0->world->photosphere.plane->normal); 
  //pm is the rotation matrix to the plane perpendicular to the photosphere

  for(i=0;i<horde->n;i++){ /* set the z-coord of the hull vertex */
    /* use temp_p to find the derotated and de-translated position of
       the hull point. copy this new position into the old position.*/
    voronoi_buf[i].p[2]=0;
    vec_mmult_3d(temp_p,pm,voronoi_buf[i].p);
    sum_3d(temp_p,temp_p,v->x);
    cp_3d(voronoi_buf[i].p,temp_p);
  }
  
  if(horde->n==0) {
    printf("VERTEX %d: hull_2d gave up! That's odd....\n",v->label);
  }

  if(verbosity >= 5){
    printf("hull_neighbors_photosphere:  hull trimming gives:\n");
    for(i=0;i<horde->n;i++) 
      printf("%d ",((VERTEX *)(horde->stuff[i]))->label); fflush(stdout);
    printf("\n");
  }
  
  return voronoi_buf;
}


/**********************************************************************
 * photosphere_vertex_update_neighbors 
 * 
 * Updates a begin/end vertex's neighbor list using
 * gather_photosphere_neighbor_candidates followed by
 * winnow_neighbor_candidates and hull_neighbors.
 * 
 * Returns the hull vertex list that is used in the final winnowing
 * process, since it's needed for some of the force laws (see
 * vertex_update_mag).  
 */

HULL_VERTEX *photosphere_vertex_update_neighbors(VERTEX *v, char global, int *n_ptr) {
  //n is the final number of hull/neighbor vertices
  DUMBLIST *dl;
  HULL_VERTEX *hv;
  int i,j;
  //DUMBLIST *vn; //vertex neighbors
  int verbosity = v->line->fc0->world->verbosity;

  if(!v) {
    fprintf(stderr,"Vertex_update_neighbors: null vertex ignored!\n");
    return 0;
  }

   if(v != v->line->start && v != v->line->end){
    printf("Oops, photosphere_vertex_update_neighbors got a non begin/end vertex (label %d) \n",v->label);
    return 0;
  }

 if(v->line->fc0->world->photosphere.type != 1) {
    printf("Oops, photosphere_vertex_update_neighbors got a vertex (label %d) that was not on a planar surface \n",v->label);
    return 0;
  }

 //vn = &(v->neighbors);

  dl = gather_photosphere_neighbor_candidates(v,global);


  if(verbosity >= 3) {
    printf("photosphere_vertex_update_neighbors: Vertex %d has %d candidates: ",v->label,dl->n);
    if(verbosity >= 4)
      for(i=0;i<dl->n;i++)
	printf("  %d",((VERTEX *)((dl->stuff)[i]))->label);
    printf("\n");
  }

  hv = photosphere_hull_neighbors(v, dl); /* save hv for return */
  
  if(verbosity >= 3) {
    printf("Hull_neighbors returned %d neighbors: ",dl->n);
    if(verbosity >= 4)
      for(i=0;i<dl->n;i++) 
	printf("  %d",((VERTEX *)((dl->stuff)[i]))->label);
    printf("\n");
  }

  *n_ptr = dl->n;

  return hv;
}
