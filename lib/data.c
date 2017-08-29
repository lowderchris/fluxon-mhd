/* data.c 
 *
 * Basic data manipulation routines and definitions for FLUX -- how to 
 * handle VERTEXs and FLUXONs etc.
 * 
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Southwest Research Institute, 2004-2008
 * 
 * You may modify and/or distribute this software under the terms of
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
 * Functions:
 *  barf - debugging function to find possible problem function(s)
 *  new_label - generates long-int label for a fluxon
 *  new_vertex_label - generates long-int label for a vertex
 *  new_fluxon - Generates a new, empty field line structure
 *  delete_fluxon - Removes a fluxon from the world, and frees it.
 *  new_vertex - Generates a vertex from X,Y,Z, and fluxon references.
 *  new_flux_concentration - generates a new flux concentration at a location
 *  new_world - creates a new world
 *  unlink_vertex - unlinks and vertex and links prev and next to eachother
 *  delete_vertex  - Unlink and delete a vertex. (unlinks neighbors too)
 *  add_vertex_pos - adds a vertex at a given position in a fluxon
 *  add_vertex_after - adds vertex after a given vertex in a fluxon
 *  vertex_renumber - changes the label of a vertex, fixing trees appropriately.
 *  clear_links - Initialize a tree-link data structure.
 *  tree_top - Returns a pointer to the top node of a tree.
 *  tree_walk - calls tree_walker
 *  tree_walker - performs a task to each node in a tree
 *  stw_helper - helper for safe_tree_walker
 *  safe_tree_walker - does what tree_walker does without fear that the 
 *                     tree will change in the middle of walking it.
 *  tree_find - Find a node in a tree (there are several kinds!).
 *  tree_binsert - (the call of choice) - Insert a node into a tree, 
 *                   and rebalance the tree if necessary.
 *  tree_insert - insert a node into a tree.
 *  tree_balance_check - checks if a tree is balanced
 *  tree_unlink - Remove a node from a tree and returns the tree root.
 *  tree_balance - Balances a tree and returns the new root.
 *  world_state_name - returns a string of the local world state name
 *  new_dumblist - initializes and empty, size zero dumblist
 *  free_dumblist - locally frees a dumblist
 *  dumblist_quickadd - adds an item to the end of a dumblist
 *  dumblist_add - Add an item to a dumblist (no duplicates allowed)
 *  dumblist_delete - Delete an item from a dumblist.
 *  dumblist_rm - removes the ith item and replaces it with the last item
 *  dumblist_grow - grows dumblist to (3/2) the requested size
 *  dumblist_qsort - qsort routine for a dumblist, used in dumblist_sort
 *  dumblist_shellsort - bubble sorts or shell sorts a dumblist, used in dumblist_sort
 *  dumblist_crunch - removes duplicates and nulls from a sorted dumblist 
 *  dumblist_sort - sorts a dumblist with qsort or shellsort
 *  ptr_cmp - compares two pointers
 *  dumblist_snarf - Copy a dumblist into another one.
 *  dumblist_clear - zeros out a dumblist
 *  dd_vertex_printer - prints out the label of a given vertex
 *  dumblist_dump - prints out the labels or the pointers of a dumblist
 *  fl_eq - tests whether two numbers are reasonably close
 * mallocs: these are switched in by compiler flags, and are used for debugging.
 *  flux_malloc - 
 *  flux_perl_malloc - 
 *  flux_padded_malloc - 
 *
 *
 *  find_vertex_by_location - just that.  Finds the closest vertex to you.
 *  interpolate_by_location - linearly interpolates a NUM valueto a location
 *  
 * This file is part of FLUX 2.2 (22-Nov-2008)
 *
 */
#include "data.h"
#include "physics.h" /* for declaration of force subroutines */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *code_info_data="%%%FILE%%%";

/******************************
 * Boundary condition names, by number
 */

char *BOUNDARY_NAMES[] = {
  "NONE",
  "PLANE",
  "SPHERE",
  "CYL"
};

/**********************************************************************
 * barf
 * Calls perror and exits, with some indication of where you started from.
 */
static char program_name[80];
void init_barf(char *whoami) {
  strncpy(program_name,whoami,79);
  program_name[79]=0;
}

void barf(int barf_on_op, char *where) {
  char buf[10240];
  switch(barf_on_op) {
  case BARF_MALLOC:
    fprintf(stderr,"Problems with Malloc!\n");
    break;
  default:
    break;
  }
  sprintf(buf,"%s(%s)",program_name,where);
  perror(buf);
  exit(barf_on_op);
}


/**********************************************************************
 * new_label
 * Generates and returns a new long-int label for a fluxon.  With 2^64 labels
 * available, we don't bother trying to conserve -- just plough through
 * them in sequence. Fieldlines and vertexs share a label numbering
 * system.  If you pass in 0, you get the next sequentially available
 * label. If you pass in a positive label number, you get that number back but
 * also the last_label cache is set to at least that number. There are 
 * potential collision problems here, but the solution is not to 
 * allocate any default-numbered nodes (e.g. from simulation) until after
 * all the set-numbered nodes (e.g. from a file) have been allocated.
 * 
 * Labels seem mainly useful for I/O (storing stuff to files), but
 * they come in handy for debugging too. It might be possible 
 * eventually to remove them from the code.
 */
static long max_label=0;

long new_label(long request){

  if(!request) 			/* request=0 */
    request = ++max_label;	/* take next available label */

  if(max_label < request)	/* else give 'em what they wanted */
    max_label=request;

  return request;
}

/**********************************************************************
 * new_vertex_label -- 
 * same as new_label, except for VERTEXes. Returns a long-int.
 */

static long max_vertex_label=0;
long new_vertex_label(long request) {
  long out;
  switch(request) {
  case 0:			/* request=0 */
    out = ++max_vertex_label;	/* take next available label */
    break;
  default:
    out = request;		/* else give 'em what they wanted */
    if(max_vertex_label < request)
      max_vertex_label = request;
    break;
  }
  return out;
}

static int randomized = 0;
long hash_vertex_label(long request, WORLD *w) {
  long new_label;
  long count = 0;

  if(randomized == 0) {
    randomized = 1;
    srandom( getpid() * time(NULL) );
  }

  if(w==0)
    return new_vertex_label(request);

  if(!request) {
    new_label = random() & 0x7fffffff;
  }  else {
    new_label = request;
  }

  // Try ten times to find a non-collision
  while(tree_find( w->vertices, new_label, v_lab_of, v_ln_of ) && count++ < 10) {
    new_label = random() & 0x7fffffff;
  }

  if(count>=10) {
    fprintf(stderr,"Hey! hash_vertex_label gave up after 10 tries.  I give up!\n");
    exit(123);
  }
  // probably not too smart -- md5 will fill up the space reasonably quickly - but what the heck.
  if(new_label > max_vertex_label)
    max_vertex_label = new_label;

  return new_label;
}
      
    

/***********************************************************************
 * new_fluxon 
 * Creates a new FLUXON data structure given the flux, begin and end
 * concentrations, and a label. Plasmoid isn't really used. The new fluxon
 * has no linked vertices. 
 */

FLUXON *new_fluxon(      NUM flux, 
			 FLUX_CONCENTRATION *c1, 
			 FLUX_CONCENTRATION *c2, 
			 long label,
			 char plasmoid)
{
  int f_no = 0;
  struct FLUXON *nf;			/* new fluxon, nf, allocate memory */

  if(c1 == NULL) {
    fprintf(stderr,"new_fluxon:  no starting flux conc.!  I refuse...\n");
    return NULL;
  } 
  if(c2 == NULL) {
    fprintf(stderr,"new_fluxon: no ending flux conc.!  I refuse...\n");
    return NULL;
  }
  if(c1==c2) {
    fprintf(stderr,"new_fluxon: can't start and end at the same flux concentration!\n");
    return NULL;
  }
  if(c1->lines && c1->lines->fc1==c1) {
    fprintf(stderr,"new_fluxon: can't start at a sink flux concentration!\n");
    return NULL;
  }
  if(c2->lines && c2->lines->fc0==c2) {
    fprintf(stderr,"new_fluxon: can't end at a source flux concentration!\n");
    return NULL;
  }
  
  nf = (struct FLUXON *)localmalloc(sizeof(struct FLUXON),MALLOC_FLUXON);

  if(!nf) barf(BARF_MALLOC,"new_fluxon");
  
  nf->flux            = flux;
  nf->label           = new_label(label);
  nf->start           = 0;		/* begin vertex */
  nf->end             = 0;		/* end vertex */
  nf->v_ct            = 0;		/* vertex count */
  
  clear_links(&(nf->all_links));	/* new fluxon not connected to others */
  nf->all_links.sum = nf->flux;
  

  clear_links(&(nf->start_links));
  nf->start_links.sum = nf->flux;
  nf->fc0 = c1;				/* start concentration */

  clear_links(&(nf->end_links));
  nf->end_links.sum = nf->flux;
  nf->fc1 = c2;				/* end concentration */

  nf->plasmoid = 0;

  nf->fc0->world->lines = tree_binsert(nf->fc0->world->lines, nf, fl_lab_of, fl_all_ln_of);
  nf->fc0->lines = tree_binsert(nf->fc0->lines, nf, fl_lab_of, fl_start_ln_of);
  nf->fc1->lines = tree_binsert(nf->fc1->lines, nf, fl_lab_of, fl_end_ln_of);

  return nf;
}

/***********************************************************************
 * new_vertex
 *
 * Returns a vertex given an X,Y,Z, and fluxon. 
 * The neighbor list is empty.
 * 
 * Mass, T, and momentum are set to 0, which is almost 
 * certainly Wrong.
 *
 */

VERTEX *new_vertex(long label, NUM x, NUM y, NUM z, FLUXON *fluxon) { 
  VERTEX *tp;				/* new vertex name */
  WORLD *w;
  int i;
  static int v_ct = 0;

  if(!fluxon || !fluxon->fc0 || !fluxon->fc0->world) {
    fprintf(stderr,"new_vertex: got a bad fluxon, flux concentration, or world! No new vertex for you!\n");
    return 0;
  }
  w = fluxon->fc0->world;

  tp = (VERTEX *)localmalloc(sizeof(VERTEX),MALLOC_VERTEX);

  tp->energy = 0;

  if(!tp) barf(BARF_MALLOC,"new_vertex");
  tp->line = fluxon;
  tp->prev = 0;				/* not linked to other vertices */
  tp->next = 0;

  tp->x[0] = x;				/* tp position */
  tp->x[1] = y;
  tp->x[2] = z;

  dumblist_init( &(tp->neighbors) );
  dumblist_init( &(tp->nearby) );

  // Pre-grow the neighbor lists to avoid having to do it later.
  dumblist_grow( &(tp->neighbors), N_NEIGHBORS_PREALLOC );
  dumblist_grow( &(tp->nearby),    N_NEIGHBORS_PREALLOC );

  tp->b_mag = 0;
  tp->b_vec[0] = tp->b_vec[1] = tp->b_vec[2] = 0;

  //  tp->label = label ? hash_vertex_label(label, fluxon->fc0->world) : hash_vertex_label(0, fluxon->fc0->world);
  /* DAL try hash_ --> new_ here (twice) to revert to the old behavior
     of making the new vertex labels sequential instead of randomized
     by md5.  Sequential is helpful for making runs reproducible when
     debugging.
   */
  tp->label = label ? new_vertex_label(label) : new_vertex_label(0);

  tp->rho = 0;
  tp->p[0] = tp->p[1] = tp->p[2] = 0;
  tp->T = 0;

  clear_links(&(tp->world_links));

  w->vertices = tree_insert( w->vertices, tp, v_lab_of, v_ln_of );
  if(v_ct++ % 100 == 0) {
    w->vertices = tree_balance(w->vertices, v_lab_of, v_ln_of);
  }

  return tp;
}

/**********************************************************************
 * new_flux_concentration
 * 
 * Returns a flux concentration at the given location within the world
 *
 */

FLUX_CONCENTRATION *new_flux_concentration(
					  WORLD *world,
					  NUM x, NUM y, NUM z,
					  NUM flux,
					  long label
					   ) {
  FLUX_CONCENTRATION *fc;
  
  fc = (FLUX_CONCENTRATION *)localmalloc(sizeof(FLUX_CONCENTRATION),MALLOC_VERTEX);
  if(!fc) barf(BARF_MALLOC,"new_flux_concentration");
  fc->world = world;
  fc->flux = flux;
  fc->label = new_label(label);

  fc->lines = (FLUXON *)0;	/* no fluxons in this fc */
  fc->x[0] = x;	 	 	/* position */
  fc->x[1] = y;
  fc->x[2] = z;
  fc->locale_radius = 0; 	/* notional size of the flux concentration... */
  
  fc->bound = 0; 
  fc->pbound = 0;

  clear_links(&(fc->links));
  fc->links.sum = fc->flux;

  world->concentrations = tree_binsert(world->concentrations, fc, fc_lab_of, fc_ln_of);

  return fc;
}

/**********************************************************************
 * check_special_concentration - checks if a FLUX_CONCENTRATION 
 * has one of the magic labels and links accordingly if it does.
 * Useful for IO.
 */
void check_special_concentration(FLUX_CONCENTRATION *f) {
  if(!f)
    return;
  if(!f->world) {
    fprintf(stderr,"check_special_concentration - only works for concentrations that are already associated with a world!\n");
    return;
  }
  switch(f->label) {
  case FC_OB_ID: 
    f->world->fc_ob = f;
    break;
  case FC_OE_ID:
    f->world->fc_oe = f;
    break;
  case FC_PB_ID:
    f->world->fc_pb = f;
    break;
  case FC_PE_ID:
    f->world->fc_pe = f;
    break;
  default:
    break;
  }
}


/**********************************************************************
 * new_world
 * Creates a new world with sensible beginning parameters, forces etc.
 */
WORLD *new_world() {
  int i;
  WORLD *a = (WORLD *)localmalloc(sizeof(WORLD),MALLOC_WORLD);
  if(!a) barf(BARF_MALLOC,"new_world");
  
  a->frame_number = 0;
  a->state = WORLD_STATE_NEW;
  a->refct = 0;

  a->concentrations = NULL;	/* nothing in the world yet */
  a->lines = NULL;
  a->vertices = NULL;

  a->photosphere.type = 0;
  a->photosphere.plane = NULL;
  a->photosphere2.type = 0;
  a->photosphere2.plane = NULL;

  /* Initialize the two dummy vertices for the mirroring for each
   * photosphere this requires two(4) dummy vertices connected onto a
   * dummy fluxon.  */

  a->fc_im10 = new_flux_concentration(a,0,0,0,0,FC_IM10_ID);
  a->fc_im11 = new_flux_concentration(a,0,0,0,0,FC_IM11_ID);
  a->fc_im20 = new_flux_concentration(a,0,0,0,0,FC_IM20_ID);
  a->fc_im21 = new_flux_concentration(a,0,0,0,0,FC_IM21_ID);
  a->fl_im = new_fluxon(0,a->fc_im10, a->fc_im11, FL_IM1_ID, 0);
  a->fl_im2 = new_fluxon(0,a->fc_im20, a->fc_im21, FL_IM2_ID, 0);
  a->image =  new_vertex(-1,0,0,0,a->fl_im);
  a->image2 = new_vertex(-2,0,0,0,a->fl_im);
  a->image3 = new_vertex(-3,0,0,0,a->fl_im2);
  a->image4 = new_vertex(-4,0,0,0,a->fl_im2);

  a->image->next = a->image2;
  a->image2->prev = a->image;

  a->image3->next = a->image4;
  a->image4->prev = a->image3;
  
  a->fl_im->start = a->image;
  a->fl_im->end = a->image2;
  a->fl_im->v_ct = 2;

  a->fl_im2->start = a->image3;
  a->fl_im2->end = a->image4;
  a->fl_im2->v_ct = 2;

  /*** By default, don't handle automatically open field lines ***/
  a->locale_radius = 0;
  a->auto_open = 0;

  /*** Init. the open-field and plasmoid pseudo flux concentrations ***/
  a->fc_ob = new_flux_concentration(a,0,0,0,1, FC_OB_ID);
  a->fc_ob->label = -1;
  a->fc_ob->bound = fl_b_open;
  a->fc_ob->locale_radius = 0;

  a->fc_oe = new_flux_concentration(a,0,0,0,-1, FC_OE_ID);
  a->fc_oe->label = -2;
  a->fc_oe->bound = fl_b_open;
  a->fc_oe->locale_radius = 0;
  
  a->fc_pb = new_flux_concentration(a,0,0,0,1, FC_PB_ID);
  a->fc_pb->label = -3;
  a->fc_pb->bound = fl_b_plasmoid;
  a->fc_pb->locale_radius = 0;
  
  a->fc_pe = new_flux_concentration(a,0,0,0,-1, FC_PE_ID);
  a->fc_pe->label = -4;
  a->fc_pe->bound = fl_b_plasmoid;
  a->fc_pe->locale_radius = 0;

  /*** By default, maintain tied boundaries by injecting vertices */
  a->default_bound = fl_b_tied_inject;

  a->verbosity = 0;  /* No verbose printouts by default */

  /* Put in a sensible default force list (physics.c) */
  a->f_funcs[0] = f_pressure_equi2b;
  a->f_funcs[1] = f_curvature;
  a->f_funcs[2] = f_vertex4;
  a->f_funcs[3] = 0;

  /* Put in NO mass list (physics.c) */
  a->m_funcs[0] = 0;

  /* Put in a blank reconnection list (physics.c) */
  a->rc_funcs[0] = 0;

  /* initialize scaling laws.  Default scaling is for b-normalized
   * forces, no acceleration
   */
  a->step_scale.b_power = 0;
  a->step_scale.d_power = 2;
  a->step_scale.s_power = 0;
  a->step_scale.ds_power = 0;

  a->passno = 0;      // initialize anticollision indicator
  a->handle_skew = 0; // by default don't skew...

  a->max_angle = 0;   // calculated value - initialize to 0
  a->mean_angle = 0;  // ditto

  a->dtau = 0.1;   // Reasonable default
  a->rel_step = 0; // No steps so far

  a->coeffs[0] = 1.0; // By default, just step naturally
  a->n_coeffs = 1;
  a->maxn_coeffs = MAXNUMCOEFFS;
  a->concurrency = 1;

  a->f_min = -1;
  a->f_max = -1;
  a->fr_min = -1;
  a->fr_max = -1;
  a->ca_min = -1;
  a->ca_max = -1;
  a->ca_acc = 0;
  a->ca_ct = 0;

  a->use_fluid = 0;  // Ignore mass by default
  a->k_b = 1.3807e-23 / (1.67e-27 * (1 * 0.87 + 4*0.13) / (2 * 0.87 + 3 * 0.13)); // particle-K.E / kelvin for coronal plasma, 13% He by number, fully ionized
  a->gravity_type = 0; // no gravity, by default.  1 = sphere, 2 = plane
  a->gravity_origin[0] = 0;
  a->gravity_origin[1] = 0;
  a->gravity_origin[2] = 0;
  a->g = 274.2 / 6.96e8 / 6.96e8 / 6.96e8 * 3600 * 3600; // default: time in hours, distance in solar radii
  

  return a;
}

/**********************************************************************
 * delete_fluxon - deletes a fluxon by unlinking it and its vertices from the 
 * surrounding structures, and then frees it.
 */
void delete_fluxon ( FLUXON *f ) {
  VERTEX *v;

  if(f->fc0->world->verbosity>=2) 
    printf("deleting fluxon %ld (%ld under it)...\n",f->label,f->all_links.n);

  if(f->start) {
    for(v=f->start->next; v && v != f->end; ) {
      delete_vertex(v);
      v=f->start->next;
    }
  }
  delete_vertex(f->start); f->start = 0;
  delete_vertex(f->end);   f->end = 0;
  f->v_ct = 0;
  f->fc0->world->lines = tree_unlink(f, fl_lab_of, fl_all_ln_of);
  f->fc0->lines = tree_unlink(f, fl_lab_of, fl_start_ln_of);
  f->fc1->lines = tree_unlink(f, fl_lab_of, fl_end_ln_of);
  {
    int v = f->fc0->world->verbosity;
    if(v>=3)
      printf("freeing the fluxon (%ld)...\n",f->label);
    if(f->label >= 0 || f->label <= -10)
      localfree(f);
    if(v>=3)
      printf("ok\n");
  }
}

/**********************************************************************
 * delete_flux_concentration
 */
void delete_flux_concentration ( FLUX_CONCENTRATION *fc ) {
  WORLD *w = fc->world;

  while(fc->lines) {
    delete_fluxon(fc->lines);
  }

  if(w->verbosity>=2)
    printf("delete_flux_concentration: deleting myself (%ld)...\n",fc->label);
  
  fc->world->concentrations = tree_unlink( fc, fc_lab_of, fc_ln_of );

  if(w->verbosity>=2)
    printf("...\n");

  if(fc->label >= 0 || fc->label <= -10) 
    localfree(fc);

  if(w->verbosity)
    printf("delete_flux_concentration: done...\n");
}
    
  

/**********************************************************************
 * free_world
 * Deallocates a world and everything in it
 */
void free_world( WORLD *w ) {
  VERTEX *v;
  FLUX_CONCENTRATION *fc;


  /////// Delete everything in the main fc tree...
  // this trickles down to all vertices etc. 
  while(w->concentrations) {
    delete_flux_concentration(w->concentrations);
  }

  ///////// Free the world...
  localfree(w);
}


/**********************************************************************
 **********************************************************************
 ***** 
 ***** List routines: handle the doubly linked list structures in the
 ***** fluxons themselves
 *****
 */

/**********************************************************************
 * unlink_vertex
 * Remove the given vertex from its list and its neighborhood.  
 * Neighborhood information gets sloshed around.
 */

#ifdef DEBUGGING_DEREF
    /**********
     * helper for debugging -- check a fluxon over for references to a particular vertex.
     * This is a static helper routine for tree_walker, used below. 
     */
     static VERTEX *srch_vtx;
     static long srch_vtx_lab;
     static char nflag;
     static long dbg_vsearch( FLUXON *f, int lb, int ln, int dp) {
       VERTEX *v;
       int i;
       for(v=f->start; v; v=v->next) {
	 if(nflag==0) {
	   for(i=0;i<v->neighbors.n;i++) {
	     char lbmatch = (((VERTEX *)(v->neighbors.stuff[i]))->label == srch_vtx_lab );
	     char vtmatch = (((VERTEX *)(v->neighbors.stuff[i])) == srch_vtx );
	     if(lbmatch && vtmatch) {
	       printf("Vertex %d <-- neighbor %d\n",srch_vtx_lab,v->label);
	     } else if(lbmatch) {
	       printf("Vertex %d (but not the one we want) <-- neighbor %d\n",srch_vtx_lab,v->label);
	     } else if(vtmatch) {
	       printf("Vertex %d <-- neighbor %d but some sort of corruption has damaged the link\n",srch_vtx_lab,v->label);
	     }
	   }
	 } else {
	   for(i=0;i<v->nearby.n;i++) {
	     char lbmatch = (((VERTEX *)(v->nearby.stuff[i]))->label == srch_vtx_lab );
	     char vtmatch = (((VERTEX *)(v->nearby.stuff[i])) == srch_vtx );
	     if(lbmatch && vtmatch) {
	       printf("Vertex %d <-- nearby %d\n",srch_vtx_lab,v->label);
	     } else if(lbmatch) {
	       printf("Vertex %d (but not the one we want) <-- neighbor %d\n",srch_vtx_lab,v->label);
	     } else if(vtmatch) {
	       printf("Vertex %d <-- nearby %d but some sort of corruption has damaged the link\n",srch_vtx_lab,v->label);
	     }
	   }
	 }
       }
       return 0;
     }
#endif


static int ptr_cmp(void *a, void *b) { /* Helper routine for sorting by pointer */
  if(a<b) return -1;
  else if(a>b) return 1;
  return 0;
}

void unlink_vertex(VERTEX *v) {
  int i;


  //printf("unlink_vertex: %d has %d neighbors and %d nearby\n",v->label,v->neighbors.n, v->nearby.n);
#ifdef NEVER
  {
    int ii;
    printf("unlink_vertex: unlinking %ld.\n\tNeighbors are: ",v->label);
    for(ii=0;ii<v->neighbors.n;ii++) {
      printf("%ld ",v->neighbors.stuff[ii]?(((VERTEX *)(v->neighbors.stuff[ii]))->label):0);
    }
    printf("\n\tNearby is: ");
    for(ii=0;ii<v->nearby.n;ii++) {
      printf("%ld ",v->nearby.stuff[ii]?(((VERTEX *)(v->nearby.stuff[ii]))->label):0);
    }
    printf("\n");
  }
#endif
  
  
  if(!v){			/* do nothing if not given a v */
    fprintf(stderr,"unlink_vertex: warning - got a null vertex!\n");
    return;
  }

  if((!v->prev || !v->next) && v->line->fc0->world->verbosity) {	/* if one of the links is missing */
    fprintf(stderr,"unlink_vertex: warning - vertex %ld is a fluxon endpoint. Proceeding anyway...\n",v->label);
  }

  if(!v->line) {		/* if it doesn't belong to a fluxon */
    fprintf(stderr,"unlink_vertex: warning - vertex %ld has no fluxon!  Proceeding anyway...!\n",v->label);
    return;
  }

  if(v->prev) 
    v->prev->next = v->next;
  else 
    v->line->start = v->next;

  if(v->next)
    v->next->prev = v->prev;
  else
    v->line->end = v->prev;
      
  v->line->v_ct--;		/* decrease vertex count in the fluxon */


  /* Shuffle around the neighbor and nearby vertex lists to remove v
   * and put in v's next and previous vertices. 
   */

  dumblist_crunch(&(v->neighbors));


  for(i=0;i<v->neighbors.n;i++) {
    int j;
    VERTEX *a = ((VERTEX **)(v->neighbors.stuff))[i];

    if( a && a != v ) { //could we get rid of both of these conditions? We just ran dumblist_crunch on v->neighbors(and also in neary loop below)  Keeping in for now, with the NOTICE below.
      if(v->line->fc0->world->verbosity>=3)
	printf(" neighbor %d: purging %ld's nearby link to %ld\n",i,a->label,v->label);

      long n = a->nearby.n;
      dumblist_delete( &(a->nearby), v);
      if(n == a->nearby.n) {
  	printf("WHOA THERE! doomed vertex %ld's neighbor %ld had no nearby link back!\n",v->label,a->label);

	//DAL print extra diagnostics
	printf("%ld's nearby list is: ",a->label);
	for(j=0; j<a->nearby.n; j++) {
	  printf("%ld ",( ((VERTEX **)(a->nearby.stuff))[j] )->label);
	}
	printf("\n");
	printf("doomed vertex %ld has neighbor list: ",v->label);
	for(j=0; j<v->neighbors.n; j++) {
	  printf("%ld ",( ((VERTEX **)(v->neighbors.stuff))[j] )->label);
	}
	printf("\n");
	printf("(we're on the %dth element of the neighbor list)\n",i);
      }

      if(v->next && a != v->next && a != v->next->next ) {
      	dumblist_add(    &(v->next->neighbors), a);
      	dumblist_add(    &(a->nearby), v->next);
      }

      if(v->prev && a != v->prev && a != v->prev->prev ) {
	dumblist_add(    &(v->prev->neighbors), a);
	dumblist_add(    &(a->nearby), v->prev);
      }
    } else {
      printf("NOTICE: In unlink_vertex, we reached an else condition that probably shouldn't have been possible: we just ran dumblist_crunch on v->neighbors.  Now v's neighbor #%d (vertex 'a') is at memory location %p and v is at memory location %p\n",i,a,v);
    }
    //v->neighbors.stuff[i]=0;//????????DAL added for symmetry with below, but probably non-functional
  }
  v->neighbors.n=0;

  dumblist_crunch(&(v->nearby));

  for(i=0;i<v->nearby.n;i++) {
    int j;
    VERTEX *a = ((VERTEX **)(v->nearby.stuff))[i];
    if( a && a != v ) { //kept for now with NOTICE below.
     if(v->line->fc0->world->verbosity>=3)
       printf(" nearby %d: purging %ld's neighbor link to %ld\n",i,a->label,v->label);

      long n = a->neighbors.n;
      dumblist_delete(&(a->neighbors), v);
      if(n == a->neighbors.n) {
	printf("WHOA THERE! doomed vertex %ld's nearby %ld had no neighbor link back!\n",v->label,a->label);

	//DAL print extra diagnostics
	printf("%ld's neighbor list is: ",a->label);
	for(j=0; j<a->neighbors.n; j++) {
	  printf("%ld ",( ((VERTEX **)(a->neighbors.stuff))[j] )->label);
	}
	printf("\n");
	printf("doomed vertex %ld has nearby list: ",v->label);
	for(j=0; j<v->nearby.n; j++) {
	  printf("%ld ",( ((VERTEX **)(v->nearby.stuff))[j] )->label);
	}
	printf("\n");
	printf("(we're on the %dth element of the nearby list)\n",i);

      }
      
      if(v->next && a != v->next && a != v->next->next ) {
	dumblist_add(   &(a->neighbors), v->next);
	dumblist_add(   &(v->next->nearby), a);
      }
      
      if(v->prev && a != v->prev && a != v->prev->prev ) {
	dumblist_add(   &(a->neighbors), v->prev);
	dumblist_add(   &(v->prev->nearby), a);
      }
    }  else {
      printf("NOTICE: In unlink_vertex, we reached an else condition that probably shouldn't have been possible: we just ran dumblist_crunch on v->nearby.  Now v's nearby #%d (vertex 'a') is at memory location %p and v is at memory location %p\n",i,a,v);
    }
    v->nearby.stuff[i]=0;
  }
  v->nearby.n=0;
  v->prev = 0;
  v->next = 0;

#ifdef DEBUGGING_DEREF
  {
    FLUXON *f = v->line->fc0->world->lines;
    
    printf("==== Unlinked vertex %d. post-unlink locations it's found:\n",v->label);
    srch_vtx_lab = v->label;
    srch_vtx = v;
    nflag = 0;
    tree_walker(f,fl_lab_of,fl_all_ln_of,dbg_vsearch,0);
    nflag = 1;
    tree_walker(f,fl_lab_of,fl_all_ln_of,dbg_vsearch,0);
    printf("\n\n");
    fflush(stdout);
  }
#endif

  if(v->line->fc0->world->verbosity>=3) {
    printf("Unlinking vertex %ld (up=%ld)...\n",v->label,v->world_links.up?((VERTEX *)(v->world_links.up))->label:0);
  }

  {
    void *root;
    WORLD *w = v->line->fc0->world;

    /* Normalize world root... */
    for(root=w->vertices;
	root && ((LINKS *)(root+v_ln_of))->up; 
	root=((LINKS *)(root+v_ln_of))->up)
      printf("WHOA!  Normalizing world root in unlink_vertex -- looks crazy from here!  Walking upwards...\n");
      ;
    v->line->fc0->world->vertices = root;
    
    /* Compare with actual root */
    for(root=v; root && ((LINKS *)(root+v_ln_of))->up; root=((LINKS *)(root+v_ln_of))->up);
    if(root==w->vertices) {
      if(w->verbosity>=4)
	printf("calling tree_unlink...\n");
      w->vertices = tree_unlink(v, v_lab_of, v_ln_of);
      if(w->verbosity>=4)
	printf("\n");  
    } else {
      if(w->verbosity >= 3)
	printf("no unlinking here... (derived root is %ld; actual root is %ld)\n",root?((VERTEX *)root)->label:0,v->line->fc0->world->vertices?v->line->fc0->world->vertices->label:0);
    }

    if(w->verbosity>=4){
      printf("Finished unlinking vertex %ld\n",v->label);
    }
  }
}

/**********************************************************************
 * vertex_clear_neighbors
 * Clears out the neighbor list of a vertex, removing it from the appropriate
 * nearby lists as necessary.
 */
void vertex_clear_neighbors(VERTEX *v) {
  int i;
  VERTEX *vn;

  // Clear out the nearby links
  for(i=0;i<v->neighbors.n;i++) {
    vn = ((VERTEX **)(v->neighbors.stuff))[i];
    dumblist_delete( &(vn->nearby), v );
  }

  // Clear out the neighbors
  v->neighbors.n = 0;
}

/**********************************************************************
 * vertex_add_neighbor
 * Adds a neighbor and links this vertex to the neighbor's nearby list
 */
void vertex_add_neighbor(VERTEX *v, VERTEX *neighbor) {
  dumblist_add( &(v->neighbors), neighbor );
  dumblist_add( &(neighbor->nearby), v );
}

/**********************************************************************
 * delete_vertex
 * Unlinks and then zero's out a given vertex. 
 */

void delete_vertex(VERTEX *v) {
  if(!v) {
    fprintf(stderr,"delete_vertex - got a null vertex!\n");
    return;
  }
  //  printf("delete_vertex: deleting %d\n",v->label);

  unlink_vertex(v);
  if(v->neighbors.stuff) {		/* zero out the neighbor links */
    localfree(v->neighbors.stuff);
    v->neighbors.stuff=0;
    v->neighbors.size=0;
    v->neighbors.n=0;
  }
  if(v->nearby.stuff) {			/* zero out the nearby links */
    localfree(v->nearby.stuff);
    v->nearby.stuff=0;
    v->nearby.size=0;
    v->nearby.n=0;
  }

  if(v->world_links.up || v->world_links.left || v->world_links.right) {
    fprintf(stderr,"WARNING: NOT freeing vertex %ld (still linked into the world tree!)\n",v->label);
    return;
  }

  localfree(v);
}

/**********************************************************************
 * add_vertex_pos 
 * Stick vertex <v> into fluxon <f> at position <pos> within the fluxon.
 * Return int status 0=success, 1=error.  
 *
 * Pos Input:  
 * 	0	trivial initial node 
 * 	-1 	trivial last node
 * 	-2 	adds to end of the middle nodes (just before the last node).  
 *	1 	adds to the beginning of the list of middle nodes.  
 * Numbers are truncated such that it's impossible to scrozzle the final
 * node by feeding in a positive number -- you have to use -1 to get that.
 *
 * Calls add_vertex_after in most cases
 *
 */

int add_vertex_pos(FLUXON *f, long pos, VERTEX *v) {
  
  if(!f || !v) 			/* error if either of them is not given */
    return 1;

  if(v->next != 0 || v->prev != 0) {
    fprintf(stderr,"Hey! add_vertex_pos got a *list*, not a vertex.  Not supported yet.\n");
    return 1;
  }


  /* Handle the empty-fluxon case: start and end vertices of f become v */
  if(f->start == 0 || f->end == 0){
    if(f->start == 0 && f->end == 0) {
      f->start = v;
      f->end = v;
      f->v_ct = 1;
      v->next = 0;
      v->prev = 0;
      return 0;
    } else {
      fprintf(stderr, "Encountered a fluxon with at START but no END, or vice versa!\n");
      return 1;
    }
  }

 
  if(pos >= 0) {		/* pos is positive */
    int i;
    VERTEX *v0;

    /* count backwards if pos is near the end */
    if(pos > f->v_ct/2) {
      for(i=f->v_ct-1, v0 = f->end; i>pos; i--)
	v0=v0->prev;
    } else {
      for(i = 0, v0 = f->start; v0->next && i<pos; i++)
	v0 = v0->next;
    }

    /* else count forwards */    
    if( f->v_ct <= pos ) {
      if(!f->end->prev) {
	fprintf(stderr,"add_vertex_pos got a 1-node fl with a large pos. position.  Not kosher!\n");
	return 1;
      }
      return add_vertex_after(f, f->end->prev, v);
    }
    return add_vertex_after(f, v0->prev, v);
  }

  /* if pos is negative */
  {
    int i;
    VERTEX *v0;

    if( f->v_ct < -pos ) {
      return add_vertex_after(f, NULL, v);
    }
    
    /* Add optional counting-forwards here for minor speed improvement */
    for(i = -1, v0 = f->end; v0->prev && i>pos; i--) 
      v0 = v0->prev;
    
    return add_vertex_after(f, v0, v);
  }
}

/**********************************************************************
 * add_vertex_after: 
 * Stick vertex <v> into fluxon <f>, just AFTER <lucy>.  If <lucy> 
 * is 0, then <v> goes at the beginning of the list.  Return an error flag
 * (0=success, 1=error). The newly inserted vertex gets the prior and next
 * vertices' neighbor lists, too.
 */

int add_vertex_after(FLUXON *f, VERTEX *lucy, VERTEX *v) {

  /* Insert at beginning of list if <lucy> is not given 
   * (<v> inserted after <lucy>) */
  if(!lucy) {
    v->prev = 0;
    v->next = f->start;
    if(f->start) {	   /* There was at least one other vertex */
      f->start->prev = v;  
      f->v_ct++;
    }
    else {		   /* This is the first vertex */
      f->end = v;
      f->v_ct=1;           /* Enforce the v_ct to be correct */
    }
    f->start = v;
  }

  else {		   /* if there is a <lucy> */
    if(lucy && (lucy->line != f)) {  /* if <lucy> isn't on same fluxon */
      fprintf(stderr,"add_vertex_after: adding %ld to f%ld after %ld: %ld is on f%ld (0x%x), not f%ld (0x%x)!\n",v->label,f->label,lucy?lucy->label:0,lucy?lucy->label:0,lucy?lucy->line->label:0,(unsigned int)(lucy?lucy->line:0),f->label,(unsigned int)f);
      return 1;
    }
    
    /* Insert at end of list */
    if(!lucy->next) {	   /* if <lucy> is the last vertice in the fluxon */
      if(f->end != lucy) { /* Print an error and abort if inconsistent */
	fprintf(stderr,"Strange vertex list found by add_vertex_after\n");
	exit(1021);
      }
      v->prev = lucy;      	/* link <lucy> in */
      v->next = 0;
      lucy->next = v;
      f->end = v;
      f->v_ct++;
    }
    
    else {			/* Insert in middle of list */
      v->next = lucy->next;	/* link <lucy> in */
      v->prev = lucy;
      lucy->next->prev = v;
      lucy->next = v;
      f->v_ct++;
    }  
  }

  /* Fix up neighbor and nearby lists lists */
  if(v->prev) 
    dumblist_snarf(&(v->neighbors),&(v->prev->neighbors));
  if(v->next)
    dumblist_snarf(&(v->neighbors),&(v->next->neighbors));
  {
    int i;
    for(i=0;i<v->neighbors.n;i++) {
      dumblist_add( &(((VERTEX *)(v->neighbors.stuff[i]))->nearby), (void *)v);
      dumblist_add( &(((VERTEX *)(v->neighbors.stuff[i]))->neighbors), (void *)v);
      dumblist_add( &(v->nearby), v->neighbors.stuff[i]);
    }
  }

  return 0;
}

/**********************************************************************
 * vertex_renumber - change the label of a vertex. 
 * Returns 0 on success, 1 on failure (which happens if a duplicate 
 * vertex is found).  If the newlab is 0, assigns a new label automatically.
 */
int vertex_renumber(VERTEX *v, long newlab) {
  VERTEX *v2;
  WORLD *w;

  if(!v) {
    fprintf(stderr,"vertex_renumber: needs a vertex!\n");
    return 1;
  }
  if(!v->line) {
    fprintf(stderr,"vertex_renumber: vertex %ld has no fluxon!\n",v->label);
    return 2;
  }
  if(!v->line->fc0) {
    fprintf(stderr,"vertex_renumber: vertex %ld's fluxon (%ld) needs a flux concentration!\n",v->label,v->line->label);
    return 3;
  } 
  if(!v->line->fc0->world) {
    fprintf(stderr,"vertex_renumber: vertex %ld (fluxon %ld, fc %ld) seems not to have a world!",v->label, v->line->label, v->line->fc0->label);
    return 4;
  }
  w = v->line->fc0->world;

  if(!newlab)
    newlab = new_label(0);


  v2 = (VERTEX *)(tree_find(w->vertices, newlab, v_lab_of, v_ln_of));
  if(v2) {
    fprintf(stderr,"vertex_renumber: %ld -> %ld no can do: vertex %ld already exists!\n",v->label,newlab,newlab);
    return -1;
  }
  
  w->vertices = tree_unlink(v, v_lab_of, v_ln_of);
  v->label = newlab;
  w->vertices = tree_binsert(w->vertices, v, v_lab_of, v_ln_of);

  return 0;
}

int fluxon_renumber(FLUXON *f, long newlab) {
  FLUXON *f2;
  WORLD *w;
  FLUX_CONCENTRATION *fc0, *fc1;
  if(!f) {
    fprintf(stderr,"fluxon_renumber: needs a fluxon!\n");
    return 1;
  } 
  if(!f->fc0) {
    fprintf(stderr,"fluxon_renumber: fluxon %ld needs a flux concentration!\n",f->label);
    return 2;
  }
  if(!f->fc0->world) {
    fprintf(stderr,"fluxon_renumber: fluxon %ld (fc %ld) seems not to have a world!",f->label, f->fc0->label);
    return 3;
  }
  w = f->fc0->world;
  fc0 = f->fc0;
  fc1 = f->fc1;

  if(!newlab) {
    newlab = new_label(0);
  }

  f2 = (FLUXON *)(tree_find(w->lines, newlab, fl_lab_of, fl_all_ln_of));
  if(f2) {
    fprintf(stderr,"fluxon_renumber: %ld->%ld no can do: fluxon %ld already exists!\n",f->label,newlab, newlab);
    return -1;
  }
  
  w->lines = tree_unlink(f, fl_lab_of, fl_all_ln_of);
  fc0->lines = tree_unlink(f, fl_lab_of, fl_start_ln_of);
  fc1->lines = tree_unlink(f, fl_lab_of, fl_end_ln_of);
  f->label = newlab;
  w->lines = tree_binsert(w->lines, f, fl_lab_of, fl_all_ln_of);
  fc0->lines = tree_binsert(fc0->lines, f, fl_lab_of, fl_start_ln_of);
  fc1->lines = tree_binsert(fc1->lines, f, fl_lab_of, fl_end_ln_of);
  
  return 0;
}

int concentration_renumber(FLUX_CONCENTRATION *fc, long newlab) {
  FLUX_CONCENTRATION *fc2;
  WORLD *w;
  
  if(!fc) {
    fprintf(stderr,"concentration_renumber: needs a flux concentration!\n");
    return 1;
  }
  if(!fc->world) {
    fprintf(stderr,"concentration_renumber: fc %ld seems to have no world!\n",fc->label);
    return 2;
  }
  w = fc->world;
  
  fc2 = (FLUX_CONCENTRATION *)(tree_find(w->concentrations, newlab, fc_lab_of, fc_ln_of));
  if(fc2) {
    fprintf(stderr,"concentration_renumber: %ld->%ld no can do: concentration %ld already exists!\n",fc->label, newlab, newlab);
    return -1;
  }
  
  w->concentrations = tree_unlink( fc, fc_lab_of, fc_ln_of );
  fc->label = newlab;
  w->concentrations = tree_binsert( w->concentrations, fc, fc_lab_of, fc_ln_of );
  return 0;
}

  
  

  

/**********************************************************************
 **********************************************************************
 *****
 *****  The tree routines:  handle multiply linked tree structures
 *****  with long-int labels and [left,right,up] link structures.
 *****
 *****  If you use a different struct type with these routines,
 *****  make sure that the first element is a NUM that can be
 *****  running-summed over.
 *****
 *****
 */

/**********************************************************************
 * clear_links
 * Clears out the links for a new data structure.
 */
void clear_links(LINKS *links) {
  links->left  = NULL;
  links->right = NULL;
  links->up    = NULL;
  links->n     = 1;
  links->sum   = 0;
  links->balanced =1;
}

/**********************************************************************
 * tree_top
 * Returns the top node of a tree
 */
void *tree_top(void *tree, int link_offset) {
  long i;
  for(i=0; ((LINKS *)(tree+link_offset))->up && i<MAX_TREE_DEPTH; i++)
    tree= ((LINKS *)(tree+link_offset))->up;
  if(i==MAX_TREE_DEPTH) {
    fprintf(stderr,"Exceeded maximum tree depth!  Something's wrong here...\n");
    exit(-1);
  }
  return tree;
}

/**********************************************************************
 * tree_walker 
 * Recursive function to walk along a tree, calling a function once
 * per node in label order. It gets the node pointer and a depth
 * meter.  Returns int 0 on normal completion, int 1 on error.
 *
 * tree_walk is a convenience function for tree_walker, which does the
 * work and also carts along a depth meter.
 */
long tree_walker(void *tree
		, int label_offset
		, int link_offset
		, long((*func)())
		, int depth){
  long err=0;

  LINKS *l = (LINKS *)(tree+link_offset);

  if(!tree)			/* not given a tree */
    return 1;

  if(l->left) 			/* go left and call tree_walker */
    err = tree_walker(l->left,label_offset, link_offset, func, depth+1);
  
  if(!err)			/* perform func on node (if err=0) */
    err = (*func)(tree, label_offset, link_offset,depth);

  if(!err && l->right)		/* go right and call tree_walker (if err=0) */
    err = tree_walker(l->right, label_offset, link_offset, func, depth+1);

  return err;
}
  
long tree_walk(void *tree
	       , int label_offset, int link_offset
	       , long ((*func)())){
  return tree_walker(tree,label_offset,link_offset,func,0);
}

/**********************************************************************
 * safe_tree_walker

 * Recursive function to walk along a tree, calling a function once
 * per node in label order. It gets the node pointer and a depth
 * meter.  Returns int 0 on normal completion, int 1 on error.
 *
 * Difference from regular tree_walker: If the function func changes
 * the tree such that it has to be re-balanced, this routine does not
 * care. It makes a copy of the tree at the first step and works with
 * that tree as opposed to the current balanced tree.
 *
 * tree_walk is a convenience function for tree_walker, which does the
 * work and also carts along a depth meter.
 */
static void **glist;
//sets each node to the right pointer in glist and hence list
long stw_helper(void *node, int lab, int link, int depth){
  *(glist++)=node;
  return 0;
}

long safe_tree_walker(void *tree
		, int label_offset
		, int link_offset
		, long((*func)())
		, int depth){
 
  void **list;
  int n = ((LINKS *)(tree+link_offset))->n;
  int i;

  if(!tree)			/* not given a tree */
    return 1;

  //make an n-size array of pointers, 
  list = (void **)localmalloc(sizeof(void *)*n, MALLOC_MISC);

  glist = list; //point global list to same place as list

  //make each pointer in glist point to the right thing.
  tree_walker(tree,label_offset,link_offset,stw_helper,depth);

  //do the function calling to each element of the list
  for( i=0; i<n; i++){
    int err;
    err = (*func)(list[i], label_offset, link_offset,1);
    if (err){
      localfree(list);
      return err;
    }
  }
  localfree(list);
  return 0;
}


/**********************************************************************
 * tree_find
 * Given a data type that includes tree links and an offset to the label
 * field and the tree link array, find the lebelled item in the tree,
 * or NULL. Recursive function.
 *
 * Returns a pointer to the desired node with that label.
 */
void *tree_find(void *tree, long label, int label_offset, int link_offset){
  long node_label;
  LINKS *links;   

  if(tree==NULL) return NULL;	/* if no tree */
  node_label = *((long *)(tree+label_offset));
  if(node_label == label) {
    return tree;
  }

  links = (LINKS *)((void *)tree + link_offset);

  /* go left if desired label is bigger than current one, call tree_find */
  if( node_label > label) {
    if( links->left == NULL ) 
      return NULL;

    return tree_find( links->left,
		      label,
		      label_offset, 
		      link_offset 
		      );
  } 

  if( links->right == NULL ) 
    return NULL; /* if it isn't on the tree */


  /* go right if desired label is smaller than current one, call tree_find */
  return tree_find( links->right,
		    label,
		    label_offset,
		    link_offset
		    );
}

/**********************************************************************
 * tree_insert
 * Given a tree and a datum, insert the datum into the tree.  
 * Balancing is given proper thought, but not done.  The root of the tree
 * is returned (because it might change!)
 *
 * Because this is a single-node insertion only, the tree 
 * links from the item should both be NULL before insertion. It will not
 * insert an item that is already in the tree, but it will intert and item
 * that has the same label as an item already in the tree (with a warning). 
 *
 * As an efficiency measure (in part because tree_balance uses tree_insert
 * to shuffle things around), tree_insert does not call tree_balance on 
 * branches that may become imbalanced.  It's best to call tree_binsert
 * if that's what you want, which calls tree_intert followed by tree_balace.
 */

void *tree_binsert(void *root, void *item, int label_offset, int link_offset) {
  return tree_balance(tree_insert(root, item, label_offset, link_offset),
		      label_offset, link_offset);
}
  
void *tree_insert(void *root, void *item, int label_offset, int link_offset) {
  long node_label;
  LINKS *item_links = (LINKS *)(item+link_offset);
  void *foo;
  long foo_label;

  if(!item) {
    fprintf(stderr,"HEY! tree_insert can't insert null items! I'm giving this tree back unchanged.\n");
    return root;
  }

  node_label = *((long *)(item+label_offset));

  /* Check assertion that this is a simple item and not a tree */
  if( (item_links->left != NULL) ||
      (item_links->right != NULL) ||
      (item_links->up != NULL) ||
      (item_links->n > 1)) {
    fprintf(stderr,"HEY!  tree_insert got item #%ld, which is a tree! I refuse to insert this.\n",(long)item_links);
    return root;
  }
  
  /* Handle the no-root-at-all case */
  if(root == NULL) {
    (((LINKS *)(item+link_offset))->n)=1;
    return item;
  }

  if( ((LINKS *)(root+link_offset))->up != NULL ) {
    NUM a;
    NUM b;
    printf("tree_insert -- warning: up from here is nonzero!  Throwing an arithmetic exception...\n");
    a=1.0;
    b=0.0;
    a/= b;
  }


  /* Find the proper parent location for the new item */
  for( foo=root; foo; ) {
    long foo_label = *( (long *)(foo+label_offset) );
    LINKS *foo_links = (LINKS *)(foo+link_offset);

    if( foo_label < node_label ) {	/* move to right */
      if(foo_links->right == NULL) {
	foo_links->right = item;
	break;
      } else {
	if( ((LINKS *)(foo_links->right+link_offset))->up != foo) {
	  NUM a,b;
	  printf("Right link's up doesn't link back!  Killing myself with an arithmetic exception...\n");
	  a=1;  b=0;  a/=b;
	}
	foo = foo_links->right;
      }
    }

    else if( foo_label > node_label) {	/* move to left */
      if(foo_links->left == NULL) {
	foo_links->left = item;
	break;
      } else {
	if( ((LINKS *)(foo_links->left+link_offset))->up != foo) {
	  NUM a,b;
	  printf("Left link's up doesn't link back!  Killing myself with an arithmetic exception...\n");
	  a=1;  b=0;  a/=b;
	}
	foo = foo_links->left;
      } 
    }

    else if(foo == item) 
      break; /* Item is already in tree... */
    else {
      /* Item is not already in tree but label is not unique. */
      fprintf(stderr,"tree_insert: Hey! Label %ld(%ld) is not unique! Proceeding anyway, but this is a real problem for you.\n",foo_label, (long)foo);

      if(foo_links->left == NULL) {
	foo_links->left = item;
	break;
      } else {
	if( ((LINKS *)(foo_links->left+link_offset))->up != foo) {
	  NUM *a = 0;
	  NUM b;
	  printf("Left link's up doesn't link back!  Killing myself...\n");
	  b = *a;
	  printf("%g",b);
	}
	foo = foo_links->left;
      }
    }
  }

  if(foo == NULL) {		/* full exit */
    fprintf(stderr,"tree_insert should never get here.  You lose.\n");
    exit(2);
  }

  /* If it's already in the tree, return. */
  if(item_links->up == foo) {
    printf("Hey! tree_insert found that this item is already in the tree!\n");
    return root;
  }

  /* Link it up! (and notate the item as a leaf) */
  item_links->up = foo;
  item_links->n = 1;
  item_links->sum = *((NUM *)item);
  item_links->balanced = 1;

  /* Increment tree-size figures and set imbalance flags for the branch */
  for(; foo; foo = ((LINKS *)(foo+link_offset))->up) {
    LINKS *fool = (LINKS *)(foo+link_offset);
    if( fool->up && 
       ((LINKS *)( fool->up + link_offset))->left != foo && 
       ((LINKS *)( fool->up + link_offset))->right != foo) {
      NUM *a = 0;
      NUM b;
      printf("Hey: tree_insert - big trouble in insertion!  Killing myself.\n");
      b = *a;
      printf("%g",b);
    }
    (fool->n)++;
    (fool->sum) += *((NUM *)item);
    tree_balance_check(foo,link_offset);
  }

  return root;
}

/**********************************************************************
 * tree_balance_check
 * Given a node in a tree, check whether it is balanced (non-recursively)
 * using its node counts and its descendents' balanced flags.  Return
 * the balanced flag (0=unbalanced, 1=balanced), but also put it into the node.
 * Allows 30% slop in the size of the two branches before it complains.
 */

char tree_balance_check(void *foo, int link_offset) {
  LINKS *l;
  long ln, rn;
  char lb, rb;

  if(!foo) return 1;		/* no tree is a balanced tree */

  l = (LINKS *)(foo+link_offset);

  ln = l->left  ? ((LINKS *)(l->left  + link_offset))->n : 0;
  rn = l->right ? ((LINKS *)(l->right + link_offset))->n : 0;

  lb = l->left  ? ((LINKS *)(l->left  + link_offset))->balanced : 1;
  rb = l->right ? ((LINKS *)(l->right + link_offset))->balanced : 1;

  return 
    l->balanced = 
    (lb && rb)
    ? ( ( ((ln - rn) <= 1) && ((rn - ln) <= 1) )
	|| ( fabs((float)ln - (float)rn)/((float)ln+(float)rn) < 0.3) ) 
    : 0;	/* are both nonzero? yes, do branches have similar n,
		 * then tree is balanced. no, tree not balanced */

}
  
/**********************************************************************
 * tree_unlink
 * Given an item in a tree, unlink the item from that tree.  The item
 * is not freed -- only unlinked from the tree and cleared!
 *
 * As always, you have to supply the offsets to the label and 
 * link fields of the tree data structure.
 * 
 * Returns the new root, which might change.
 */

void *tree_unlink(void *data, int label_offset, int link_offset) {
  LINKS *data_links = (LINKS *)(data+link_offset);
  LINKS *upline_links;
  void *foo, *root;
  void *downline = (data_links->left ? 
		    data_links->left : 
		    data_links->right); /* downline is left if possible, else right */
  
  /* Initial guess at new root is one of the branches of this node */
  root = data_links->left ? data_links->left : data_links->right;

  /* Walk up to the top, decrementing link counts as we go and finding the new root */
  for(foo=data_links->up; 
      foo; 
      foo = ((LINKS *)(foo+link_offset))->up) {
    
    ((LINKS *)(foo+link_offset))->sum -= *((NUM *)data);
    if( (--(((LINKS *)(foo+link_offset))->n))  < 0)
      fprintf(stderr,"You hoser -- negative node count found!\n");
    
    root = foo; /* Gets last nonzero link */
  }

  upline_links = ( data_links->up ? 
		    (LINKS *)(((void *)data_links->up)+link_offset) :
		    0 ); 	/* upline_links, up or null */


  /* If neither branch is null, then the right-hand branch gets 
   * grafted onto the rightmost leaf of the left-hand branch. */
  if(data_links->left && data_links->right) {
      void *leaf;
      LINKS *leaf_links;
      long n;
      NUM sum;

      /* Grab the n and sum from the right-hand tree; prepare to
       * accumulate them into the left-hand side. */
      n =   ((LINKS *)(data_links->right+link_offset))->n;
      sum = ((LINKS *)(data_links->right+link_offset))->sum;

      /* Walk down to the leaf, adding sum and n at each level */
      for(leaf = downline;
	  (leaf_links = (LINKS *)((leaf+link_offset)))->right;
	  leaf = leaf_links->right
	  ) {
	leaf_links->n   += n;
	leaf_links->sum += sum;
      }

      /* Fix up the leaf: sum, n, links */
      leaf_links->n   += n;
      leaf_links->sum += sum;
      leaf_links->right = data_links->right;
      ((LINKS *)((data_links->right)+link_offset))->up = leaf;

      /* Mark the right-hand branch as fully grafted to the left */
      data_links->right = 0;
  }


  /* Now there are either 0 or 1 branches below this one, and 
   * downline is set appropriately (if there were initially two
   * then downline is pointing to the left-hand one).
   * Link them in to the parent. */

   if(upline_links) {
     if(upline_links->left == data) 
       upline_links->left = downline;
     else if(upline_links->right == data)
       upline_links->right = downline;
     else
       fprintf(stderr,"tree_delete: HEY!  The tree is unglued...\n");
   } else 
     root = downline;
   
   if(downline) 
     ((LINKS *)(downline+link_offset))->up = data_links->up;
   
   clear_links(data_links);
   return root;
}

/**********************************************************************
 * tree_balance
 *
 * If there's an imbalance between left and right sides, 
 * unlink the rightmost <n> leaves from the left side or the leftmost
 * <n> leaves from the right side, and stick 'em on the opposite side.
 * Then recursively rebalance both subtrees.  Return the new tree root.
 *
 * This could be made more efficient by optimizing the tree_unlink
 * call a bit (and keeping it local):  rather than finding the rightmost
 * leaf every time, it's possible to "hold one's place".  
 * 
 * Note to self:  copy a real tree algorithm someday.  (not that this
 * bit of code will be the slowest in the bunch -- the "hot spot" is
 * in the geometry section.)  2-oct-2000
 *
 */

void *tree_balance(void *tree, int label_offset, int link_offset) {
  int i;
  void *l_tree, *r_tree;
  void *upnode;
  long l_n, r_n;
  LINKS *l_link, *r_link, *t_link;
  
  if(tree == NULL) 			/* a null tree is balanced */
    return tree;

  t_link = (LINKS *)(tree+link_offset);
  l_tree = t_link->left;
  r_tree = t_link->right;

  l_link = (LINKS *)(l_tree + link_offset);
  r_link = (LINKS *)(r_tree + link_offset);

  l_n = l_tree ? l_link->n : 0;
  r_n = r_tree ? r_link->n : 0;

  /* Some shuffling has to be done:  sever the head of the tree
   * from the branches if it isn't balanced. */
  if( ! t_link->balanced ) {
    char l_tweak=0;
    char r_tweak=0;

    if(l_tree) l_link->up = NULL;
    if(r_tree) r_link->up = NULL;
    
    t_link->left = NULL;
    t_link->right = NULL;
    t_link->up = NULL;
    t_link->n = 1;
    t_link->sum = *((NUM *)tree);
    upnode = t_link->up;

    
    /* Shuffle from left to right as necessary */
    for(i=0; i < (l_n - r_n)/2 ; i++) {
      void *foo, *f;
      l_tweak=1;
      
      r_tree = tree_insert(r_tree, tree, label_offset, link_offset);


      for(foo = l_tree; (f = ((LINKS *)(foo + link_offset))->right); (foo = f)) ;
      l_tree = tree_unlink(foo, label_offset, link_offset);

      tree = foo;
    }

    /* Shuffle from right to left as necessary */
    for(i=0; i< (r_n - l_n)/2 ; i++) {
      void *foo, *f;
      r_tweak = 1;
      
      l_tree = tree_insert(l_tree, tree, label_offset, link_offset);
      
      for(foo = r_tree; (f = ((LINKS *)(foo + link_offset))->left); (foo = f) );
      r_tree = tree_unlink(foo, label_offset, link_offset);

      tree = foo;
    }


    if(l_tweak && r_tweak) {
      fflush(stdout);
      fprintf(stderr,"balance_tree:  Got an impossible dual-tweak!\n");
      fflush(stderr);
    }
  }
   
  l_link = (LINKS *)(l_tree+link_offset);
  r_link = (LINKS *)(r_tree + link_offset);

  if(l_tree && (!l_link->balanced)) 	/* balance left subtree */
    l_tree = tree_balance(l_tree, label_offset, link_offset);

  if(r_tree && (!r_link->balanced))	/* balance right subtree */
    r_tree = tree_balance(r_tree, label_offset, link_offset);
  
  /* Mend the links and reconnect left and right trees */
  t_link = (LINKS *)(tree+link_offset);
  t_link->up = NULL;
  t_link->left = l_tree;
  t_link->right = r_tree;
  t_link->n = 1;
  t_link->sum = *((NUM *)tree);
  t_link->balanced = 1;
    
  if(l_tree) {
    l_link = (LINKS *)(l_tree + link_offset);
    l_link->up = tree;
    t_link->n += l_link->n;
    t_link->sum += l_link->sum;
  }
  
  if(r_tree) {
    r_link = (LINKS *)(r_tree + link_offset);
    r_link->up = tree;
    t_link->n += r_link->n;
    t_link->sum += r_link->sum;
  }
  
  return tree;
  
}


/**********************************************************************
 * world_state_name
 * dumb -- just produces a string with the name of the state the world
 * (the given one, not the big one) is in.  
 */

const char *world_state_name(WORLD *a){
  if(a==NULL) {
    return "NULL WORLD";
  }
  switch(a->state) {
  case WORLD_STATE_NEW:     return "NEW";
  case WORLD_STATE_LOADING: return "LOADING";
  case WORLD_STATE_LOADED:  return "LOADED";
  case WORLD_STATE_WORKING: return "WORKING";
  case WORLD_STATE_RELAXED: return "RELAXED";
  case WORLD_STATE_READY:   return "READY";
  default: return "UNKNOWN";
  }
}

/**********************************************************************
 **********************************************************************
 ***** Dumblist routines
 ***** dumb lists keep track of vertices' neighbors.  They're pretty
 ***** dumb, but here's a generic implementation.  Each is an 
 ***** array of unsorted pointers to objects, and knows how many things
 ***** are in the list and how big the allocated space is.  
 *****
 ***** They grow to the power of 1.5 closest to the high-water mark and 
 ***** stay there till you get rid of 'em.  That's useful for short lists 
 ***** of things (say under 50) where you don't want lots of malloc 
 ***** overhead.  
 *****
 ***** For sanity reasons, these things crash if you try to make a dumblist
 ***** with more than 2^15 things in it.
 *****/


/**********************************************************************
 * new_dumblist
 * initializes a new, emptry dumblist with size zero 
 */

DUMBLIST *new_dumblist() {
  DUMBLIST *foo = (DUMBLIST *)localmalloc(sizeof(DUMBLIST),MALLOC_DL);
  foo->n = 0;
  foo->size = 0;
  foo->stuff = 0;
  return foo;
}

/**********************************************************************
 * free_dumblist
 * locally frees the pointers in a dumblist and then the dumblist itself
 */ 

void free_dumblist(DUMBLIST *foo) {
  if(foo->stuff) 
    localfree (foo->stuff);
  localfree (foo);
}

/**********************************************************************
 * dumblist_quickadd
 * adds an element after the last used element in a dumblist. sometimes
 * increases the size of the dumblist if needed
 */ 

void dumblist_quickadd(DUMBLIST *dl, void *a) {
  int i;
  void **foo;
  if(!dl || !a) return;
  if(!(dl->stuff)) {	/* if its a new dumblist, give it 16 elements */
    dl->stuff = (void **)localmalloc( 16*sizeof(void *),MALLOC_DL_L);
    dl->size = 16;
    dl->n = 1;
    dl->stuff[0] = a;
    return;
  }

  if(dl->n >= dl->size)		                /* note that if n>size, you have bigger problems */
    dumblist_grow(dl, dl->size + 10);    	/* grows to more or less (3/2)*size */
  
  dl->stuff[dl->n] = a;
  dl->n++;
}

/**********************************************************************
 * dumblist_add
 * adds an element after the last used element in a dumblist. sometimes
 * increases the size of the dumblist if needed. checks to make sure that
 * the element isn't already in the dumblist
 */ 
    
void dumblist_add(DUMBLIST *dl, void *a) {
  int i;
  void **foo;
  if(!dl || !a) return;

  /* If it's a null list, allocate some space for a 16 element list (this could 
     actually fit in the extend-list conditions, but what the heck) */
  if(!(dl->stuff)) {
    dl->stuff = (void **)localmalloc( 16 * sizeof(void *),MALLOC_DL_L);
    dl->size = 16;
    dl->n = 1;
    dl->stuff[0] = a;
    return;
  }

  /* Traverse the list to see if a is already in it! (this is what makes it dumb...) */
  for(foo=dl->stuff,i=0;i<dl->n && *foo != a;foo++,i++)
    ;
  if(i<dl->n) 	/* don't add it if it is already there */
    return;
  
  /* OK, a isn't in the list -- add it on the end, increase size by a few. */
  if(dl->n >= dl->size)
    dumblist_grow(dl,dl->size  + 10);

  dl->stuff[dl->n] = a;
  dl->n++;

  return;
}

/**********************************************************************
 * dumblist_delete
 * attempt to remove an item from a dumb list. Finds each occurence 
 * of the item and blasts 'em.  If no such item is found, then
 * nothing happens.  
 */

 void dumblist_delete(DUMBLIST *dl, void *a) {
  int i;
  void **foo;

  if(!dl || !a) return;

  for(foo=dl->stuff,i=0; i<dl->n; foo++, i++) 
    if(*foo == a) {
      dumblist_rm(dl,i);
      i--;
      foo--;
    }
}
    
/**********************************************************************
 * dumblist_rm
 * Removes the ith item from a dumblist.  Move the last
 * item in the dumblist to the vanished position.
 */

 void dumblist_rm(DUMBLIST *dl, int i) {
  if(i < dl->n) {
    dl->n--;
    if( i < dl->n )	/* with the new n */
      dl->stuff[i] = dl->stuff[dl->n];
  } else {
    fprintf(stderr,"Uh oh!  dumblist_rm tried to remove off the end of the list...\n");
    {
      float a = 5;
      a /= 0.; // throw arithmetic exception
    }
  }
}


/**********************************************************************
 * dumblist_grow
 * Grows a dumblist to (3/2) the requested size or 16, whichever is bigger.
 */
void dumblist_grow(DUMBLIST *dl, int size) {
  void **foo;
  void **a,**b;
  int i;
  long newsize;

  newsize = (3 * size)/2;

  if(newsize < 16)         /* Minimum size */
    newsize = 16;

  if(newsize <= dl->size)  /* No extra work please */
    return;
  
  if(newsize > 32768) {
    if(newsize > 1073741820) {
      fprintf(stderr,"dumblist_grow:  size is too big!\n");
      exit(-30);
    }
    else {
      int a,b=0;
      fprintf(stderr,"dumblist_grow: Warning -- size is %ld!\n",newsize);
      /* j      a /= b;  throw an arithmetic exception */
    }
  }
  
  /* This seems to corrupt the arena on my linux box -- replaced 
   * with more pedestrian version below.  --CED 07-Oct-2004
   *
   * dl->stuff = (void **)realloc(dl->stuff,newsize * sizeof(void *));
   */

  {
    void **oldstuff, **ost, **newstuff;
    int i;
    oldstuff = ost = dl->stuff;
    dl->stuff = newstuff = (void **)localmalloc(newsize * sizeof(void *),MALLOC_DL_L);
    if(ost) {
      for(i=0;i<dl->size; i++)
      *(newstuff++) = *(oldstuff++);
      localfree(ost);
    }
  }
    
    dl->size = newsize;
}
    
/**********************************************************************
 * static variables */

static void **dls_wk = 0;
static long dls_size = 0;         /* Size of dumblist we can currently qsort */
static long dls_sz = 0;           /* Size of allocated space for qsorting */
static int ((*dls_cmp)(void *a, void *b));

/**********************************************************************
 * dumblist_qsort
 * helper function for dumblist_sort 
 * Sort the given dumblist -- just the n elements starting at l --
 * using the workspace starting at element wk_start.  Return the 
 * start of the sorted list in the workspace.  Eat up workspace linearly --
 * there should be plenty (never call this routine from outside; call the 
 * interface routine dumblist_sort instead!)
 * 
 * This is just your basic qsort routine -- the main reasons to 
 * reimplement it here rather than using the libraries are:
 *    (A) it's kind of fun to mess with qsort after all these years;
 *    (B) I wanted one that would munch up a known workspace with
 *        known-max-size lists, to avoid lots of mallocation.
 * 
 * This routine is doomed to be part of the 'hot spot' for the model --
 * it should probably be retweaked for maximum speed, or rewritten 
 * to omit the internal function call for each compare (pre-process the
 * list to use floating point size indices?)
 * 
 */

void **dumblist_qsort(void **l, int n, void **wk_start) {
#if SORT_WARN
  printf("Dumblist_qsort...\n");
#endif
  if(n<=1) {
    if(!n)        /* Paranoia */
      return l;

    *wk_start = *l;  /* n=1: return the input list (it's sorted!) */
    return wk_start;
  }

  {  /* Block avoids variable initialization for the trivial case */
    int nl=n/2;
    int nr=n-nl;
    int n_copied; 
    void **out,**out1;
    void **left, **left0;	/* left0 is start of left branch, left is current */
    void **right, **right0;	/* right0 is start of right branch, right is current */
    
    /* Divide the list in half (nl and nr), and sort each half. */
    if(n>2) {
      left0  = left  = dumblist_qsort(l    , nl, wk_start      );
      right0 = right = dumblist_qsort(l+nl , nr, left + nl + 1 );
    } else {  /* n==2 */
      /* No need to iterate to the hairy end for 2 trivial lists. */
      left0 = left = wk_start;	/* start out with left,left0 same place */
      right0 = right = left+1;	/* start out with right,right0 same place */
      *left = *l;
      *right = l[1];
    }
 
    /* out1 is beginning of out sorted list, the 1 is a buffer space */    
    out1 = out = right + nr + 1;
    
    /* Merge the lists, skipping duplicates and copying back into place. */
    for(n_copied = 0; 
	n_copied < n;
	n_copied++) {
      /* compare the positions of left, left0 and right, right0. they have 
       * been sorted if they are not at the same place  */
      int l_done = (left  - left0 >= nl)   || (*left == NULL)  ;
      int r_done = (right - right0 >= nr)  || (*right == NULL) ;
      int foo;

      if(!l_done) {
	if(!r_done) {
	  switch (foo=((*dls_cmp)(*left,*right))) {
	    
	  case -1:			/* *left < *right */
	    *(out++) = *(left++);
	    break;
	    
	  case 1:			/* *left > *right */
	    *(out++) = *(right++);
	    break;
	    
	  case 0:			/* *left == *right */
	    *(out++) = *(right++);
	    left++;
	    break;
	    
	  default:
	    fprintf(stderr,"Unknown response %d from sort function in dumblist_qsort.  Giving up.\n",foo);
	    exit(-10);
	    break;
	    
	  }

	} else  /* r_done && !l_done */
	  *(out++) = *(left++);
      }
      else if(!r_done) /* !r_done && l_done */
	*(out++) = *(right++);
      else /* l_done and r_done are true */ {
	break; /* Don't waste any more iterations! */
      }
    }

    *out = 0; /* Lay down a fence, buffer space */

    return out1;
  } /* end of convenience block */
} /* end of dumblist_qsort */


/******************************
 * dumblist_shellsort
 * is a basic comb/shell sorter.
 * It's useful for lists too small to need quicksort (up to, say, 50 elements).
 */

void dumblist_shellsort( DUMBLIST *dl, int ((*cmp)(void *a, void *b)) ) {
  int i,j,increment;
  void *temp;
  int n = dl->n;
  int iter=0;

#if SORT_WARN
  printf("dumblist_shellsort\n");
#endif

#ifdef testing_with_bubblesort
  /*  Bubblesort to make sure shellsort isn't scrozzling memory, if option is set */
  char done;
  int ret;

  /* come at it from both ends at once, going all the way through to the opp end */
  do {
    done = 1;
    for(i=0;i<n-1;i++) {
      if(  (ret = ((*cmp)(dl->stuff[i],dl->stuff[i+1]))) > 0) {
	void *foo = dl->stuff[i+1]; /* switch places if i > i+1 */
	dl->stuff[i+1] = dl->stuff[i];
	dl->stuff[i] = foo;
	done = 0;
      }
      if(((*cmp)(dl->stuff[n-2-i],dl->stuff[n-1-i]))>0) {
	void *foo = dl->stuff[n-2-i]; /* switch places if n-2-i > n-1-i */
	dl->stuff[n-2-i] = dl->stuff[n-1-i];
	dl->stuff[n-1-i] = foo;
	done = 0;
      }
    }
  } while(!done);
  return;

#endif

  increment = (dl->n / 3) || 1;

  while(increment>0)
    { 
      for(i=0; i < dl->n; i++) {

	j = i;
	temp = dl->stuff[i];
	while((j >= increment) && ( (*cmp)(dl->stuff[ j - increment ],temp) > 0) ) {
	  dl->stuff[j] = dl->stuff[j-increment];
	  j = j - increment;
	}
	dl->stuff[j] = temp;
      }
      increment >>= 1; /* binary /2, drop remainder, basically shorten the increment*/
    }
}

/******************************
 * dumblist_crunch 
 * crunches a dumblist, removing duplicate elements and nulls. duplicate
 * elements must be pointing to same thing to be considered duplicate.
 * 
 * The dumblist need not be sorted.
 *
 * The algorithm is O(n^2), so don't use it to crunch huge dumblists -- it's great, though,
 * for stuff the size of a typical neighbor list.
 *
 */

void dumblist_crunch(DUMBLIST *dl) {
  void **a, **b;
  int i,j;

  if(dl->n <= 1)
    return;

  // Loop over the whole list...
  for(i=0; i<dl->n; i++) {
    
    
    if(dl->stuff[i] == NULL) {
      // If the element is null, remove it.

      dumblist_rm( dl, i );
      i--;

    } else {

      // If the element is not null, sweep the rest of the dumblist for dups
      for(j=i+1; j<dl->n; j++) {

	if( dl->stuff[i] == dl->stuff[j] ) {
	  dumblist_rm(dl, j);
	  j--;
	}
      }
    }
  }
}


/**********************************************************************
 * dumblist_sort
 * calls dumblist_qsort, dumblist_shellsort and dumblist_crunch
 * Sort a dumblist by some external criterion.  You feed in 
 * a function (*cmp) that returns -1, 0, or 1 depending on whether the
 * first argument is less than, equal to, or greater than the 
 * second, and the list is sorted for you. 
 * 		cmp (a,b)
 *		-1	a<b
 *		0	a=b
 *		1	a>b
 *
 * The sorting takes place in a workspace that is static, so that
 * it isn't always being remalloc'ed. 
 *
 * See notes in dumblist_qsort -- this is doomed to be a hot-spot
 * routine, but isn't particularly well optimized.
 *
 * If the number of elements is small, then shellsort is used as
 * it's faster than qsort for the small stuff.
 *
 */

void dumblist_sort(DUMBLIST *dl, int ((*cmp)(void *a, void *b))) {
  void ** sorted;
  void **dl_a;
  int i;
  /* **dls_wk, dls_size, dls_sz and *dls_cmp are static variables */
#if SORT_WARN
  printf("dumblist_sort...\n");
#endif

  if(dl->n <= 1) 
    return;

  if( dl->n <= 150 ) {             /* use shellsort for the small cases */
    //    printf("s");fflush(stdout);
    dumblist_shellsort(dl,cmp);
    dumblist_crunch(dl);
  } else {
    int odls_size;
    //    printf("q");fflush(stdout);
    /* Make sure that there's enough scratch space */
    if(dls_size < dl->size) {
      long levels = 2;
      unsigned long a;
      /* This is kind of wasteful but should be tiny compared to the
	 main data structure -- so don't be makin' global gather_neighbor
	 calls if you can possibly avoid it.  

	 Reserve enough space for two copies at each level, and enough
	 levels to cover the whole depth.

      */
      odls_size = dls_size;
      dls_size = a = dl->size * 2;
      while(a >>= 1)   /* shift and test for nonzero */
	levels++;
      dls_sz = (2*dls_size)*levels;
      
      /* realloc seems to bust the arena -- so manually reallocate...
       * dls_wk = (void **)realloc(dls_wk, dls_sz * sizeof(void *));
       */
      if(dls_wk)
	localfree(dls_wk);
      dls_wk = (void **)localmalloc(dls_sz * sizeof(void *),MALLOC_DLS_ARENA);
    }
    
    /* Launch the sort */
    dls_cmp = cmp;
    sorted = dumblist_qsort(dl->stuff,dl->n,dls_wk);
    
    /* Copy the sorted list back into the vertex. */
    dl_a = dl->stuff;
    for(i=0;i<dl->n && *sorted;i++) {
      *(dl_a++) = *(sorted++);
    }
    dl->n = i;
  }
}
  

/**********************************************************************
 * dumblist_snarf 
 * Copy every item from the source dumblist into the destination.
 * Dumblists are not sorted.
 */

void dumblist_snarf(DUMBLIST *dest, DUMBLIST *source) {
  long i;
  void **b;
  void **b1;
  void **c;
  
  if( !dest || !source || dest==source) {
    fprintf(stderr,"dumblist_snarf: Got null src or dest, or dest==src! (d=%ld, s=%ld)\n",(long)dest,(long)source);
    fflush(stderr);
    return;
  }
  if(dest->size <= dest->n + source->n) { /* make dest bigger if need be */
    dumblist_grow(dest,dest->n + source->n);
  }
  
  for(b=source->stuff, c=&(dest->stuff[dest->n]), i=0; i<source->n; i++) {
    *(c++) = *(b++);			/* do the copying */
  }
  
  dest->n += source->n;			/* increase dest n */

}

/**********************************************************************
 * dumblist_clear
 * Zeroes out a dumblist.  Especially useful for the dumblists in VERTEXes,
 * because they're included in the base struct instead of being pointers out.
 */

 void dumblist_clear(DUMBLIST *foo) {
  foo->n = 0;
}

/**********************************************************************
 * dumblist_init
 * Prep a freshly-allocated dumblist for life
 */
 void dumblist_init(DUMBLIST *foo) {
  foo->n = foo->size = 0;
  foo->stuff = 0;
}

/**********************************************************************
 * dumblist_clean
 * De-allocate the stuff[] array and clean up the dumblist for freeing.
 */
void dumblist_clean(DUMBLIST *foo) {
  foo->n = foo->size = 0;
  if(foo->stuff) {
    localfree(foo->stuff);
  }
  foo->stuff = 0;
}

/**********************************************************************
 * dd_vertex_printer
 * prints out the vertex label of the given vertex
 */

void dd_vertex_printer(void *a) {
  printf("vertex #%ld",((VERTEX *)a)->label);
}

/**********************************************************************
 * dumblist_dump
 * prints out the labels or the pointers of a dumblist
 */

 void dumblist_dump(DUMBLIST *foo, void ((*printer)(void *a))) {
  int i;

  if(!foo) {
    printf("dumblist_dump got a null list!\n");
  } else {
    printf("list has %d elements\n",foo->n);
      for (i=0;i<foo->n;i++) {
	printf("\t%d: ",i);
	if(printer)
	  (*printer)(foo->stuff[i]);
        else 
	  printf("element ptr=%ld",(long)(foo->stuff[i]));
	printf("\n");
      }
  }
}

/**********************************************************************
 * fl_eq
 * Test whether two NUMs are reasonably close.  
 * Currently, 100ppm is the slop.
 */

 char fl_eq(NUM a, NUM b) {
  return ( (a==0 && b==0) ||
	   ( !(a==0 || b==0) &&
	     (a/b > 0) &&
	     ( (a/b < 1.0001) && (b/a < 1.0001) )
	     )
	   );
}

/**********************************************************************
 *
 * flux_malloc
 *
 * A wrapper for malloc that includes optional fencing.
 * Only gets compiled if you asked for debugging malloc with
 * -DUSE_DEBUGGING_MALLOC at compile time.
 *
 * flux_malloc allocates a fence of size fencecount*sizeof(long) 
 * on either side of each newly allocated bit of memory, then fills
 * the fences with 0xDEADBEEF, and monitors them on every memory operation.
 * If the fences get touched, a warning gets printed.
 * 
 * Of course, there is a rather severe hit on performance! 
 * 
 * With small fences (4 longs), memory contention issues still
 * appear in the perl side -- notably, the malloc arena gets corrupted
 * and the code crashes.  But none of the fences ever sustains a hit.
 * With large fences (100 longs), the crashes go away -- but none
 * of the fences ever seem to get touched. 
 * 
 */
#ifdef USE_DEBUGGING_MALLOC

static char *flux_malloc_types[MALLOC_MAXTYPENO] = 
  {"World","Fluxon","Vertex","Flux Concentration",
   "Dumblist","Dumblist stuff field","Vertex list","DLS arena",
   "Photospheric plane", "Misc"
  };

#define MALLOC_MAXNO (10000)
static struct flux_malloc_alloc {
    void *where;
    long size;
    int type;
  } flux_malloc_alloc[MALLOC_MAXNO];
  static int flux_malloc_next = 0;

#define fencecount (100)
#define fencesize  (4*fencecount)

char *flux_malloc(long size, int what_for) {
  char *p;
  long *p2;

  flux_memcheck();
  
  if(flux_malloc_next >= MALLOC_MAXNO){
    fprintf(stderr,"Filled malloc table; ending test\n");
    exit(1);
  }

  if(!valid_malloc_type(what_for)) {
    fprintf(stderr,"Invalid type received by flux_malloc\n");
    /*  Throw an exception */
    {
    int i = 0;
    int j = 2;
    j /= i;
    }
    return (void *)0; 
  }
  
  /* call malloc and pad for the fence.*/
  p = (char *)malloc(size + 2*fencesize);
  p2 = (long *)p;
  { int i;
  for(i=0;i<fencecount;i++) *(p2++) = 0xDEADBEEF;
  p = (char *)p2;
  p2 = (long *)(p+size);
  for(i=0;i<fencecount;i++) *(p2++) = 0xDEADBEEF;
  }

  flux_malloc_alloc[flux_malloc_next].where = p;
  flux_malloc_alloc[flux_malloc_next].size = size;
  flux_malloc_alloc[flux_malloc_next].type = what_for;
  flux_malloc_next++;
  return (void *)p;
}

void flux_dump_memblocks() {
  int i;
  printf("\n\n---snapshot of allocated blocks:\n");
  for(i=0;i<flux_malloc_next;i++) {
    if(flux_malloc_alloc[i].type) {
      printf("\t%.4d: type=%d, size=%d, loc=%d\n",i,flux_malloc_alloc[i].type,flux_malloc_alloc[i].size,flux_malloc_alloc[i].where);
    }
  }
  printf("\n");
}

void flux_free(void *p) {
  int i;
  int ok = -1;
  flux_memcheck();


  //  flux_dump_memblocks();


  for(i=0;ok<0 && i<flux_malloc_next;i++) 
    if (p==flux_malloc_alloc[i].where) 
      ok = i;
  
  if(ok>=0) {
    flux_malloc_alloc[ok].where = 0;
    flux_malloc_alloc[ok].size = 0;
    flux_malloc_alloc[ok].type = 0;
    free((char *)p-fencesize);
    return;
  } else {
    fprintf(stderr,"%%%%%%flux_free: got an unknown block to free\n");
    fprintf(stderr,"block was %d\n",p);
    fprintf(stderr,"Dump of blocks:\n");
    for(i=0;i<flux_malloc_next;i++) {
      fprintf(stderr,"\t%d (type=%d, size=%d)\n",flux_malloc_alloc[i].where,flux_malloc_alloc[i].type,flux_malloc_alloc[i].size);
    }
    /* Throw an exception */
    {
      int i=0,j=2;
      fprintf(stderr,"flux_free: Throwing a floating point exception...\n");
      j/=i;
    }
    free((char *)p-fencesize);
    return;
  }    
}

void flux_memcheck() {
  long badct = 0;
  long freed = 0;
  int i;
  char flag =0;
  static FILE *foo=0;

  if(!foo) {
    foo = fopen("foo.txt","w");
  }
  fprintf(foo,"%c======memcheck table:\n",(char)12);
  for(i=0;i<flux_malloc_next;i++) {
    char type_bad=0;
    char fence_bad=0;
    char *p;
    long *p2;
    fprintf(foo,"%4.4d: 0x%8.8x size=%6.6d (%4.4x) type=%d (%s)\n",
	   i,flux_malloc_alloc[i].where,
	   flux_malloc_alloc[i].size,flux_malloc_alloc[i].size,
	    flux_malloc_alloc[i].type,
	    malloc_types[flux_malloc_alloc[i].type]
	   );
    if(!flux_malloc_alloc[i].where) {
      freed++;
    } else {
      int ii;
      p = flux_malloc_alloc[i].where;
      p2 = (long *)(p-fencesize);
      
      for(ii=0; ii<fencecount;ii++) 
	fence_bad += (*(p2++) != 0xDEADBEEF);
      p2 = (long *)(p+flux_malloc_alloc[i].size);
      if(!valid_malloc_type(flux_malloc_alloc[i].type)) 
	type_bad = 1;
      for(ii=0; ii<fencecount; ii++)
	fence_bad += (*(p2++) != 0xDEADBEEF);
    }
    if(fence_bad) {
      int ii;
      char buf[fencecount*2+10];
      char *s=buf;
      /* Announce the problem */
      fprintf(stderr,"%%%%%%flux_memcheck: slot %d, 0x%x (%s) has a bad fence! (%d bad longs)\n",i,p,type_bad?"UNKNOWN_TYPE":flux_malloc_types[flux_malloc_alloc[i].type],fence_bad);
      for(ii=0;ii<fencecount;ii++) {
	*(s++) = ( ((long *)(flux_malloc_alloc[i].where - fencesize))[ii]==0xDEADBEEF ) ? 'X' : '.';
      }
      *(s++) = ' '; *(s++) = '|'; *(s++) = ' ';
      for(ii=0;ii<fencecount;ii++) {
	*(s++) = ( ((long *)(flux_malloc_alloc[i].where + flux_malloc_alloc[i].size))[ii]==0xDEADBEEF ) ? 'X' : '.';
      }
      fprintf(stderr,"  MAP: %s\n",buf);
      /* Should mend the fence here... */
    }
    
    if(fence_bad || type_bad)
      badct++;
  }
  fprintf(foo, "%sflux_memcheck (both-sides fencing): %d table entries (%d freed) OK, %d good, %d bad\n",(badct==0?"   ":"%%%"),i,freed,i-freed-badct,badct);
  fflush(foo);
}

char *malloc_types[MALLOC_MAXTYPENO+1] = {
  "Nothing",
  "World","Fluxon","Vertex","Flux Conc.",
  "Dumblist","Dumblist_stuff","Vertex List","DLS Arena",
  "Plane"};

 /* end of debugging-malloc code */

#else
#ifdef USE_PERL_MALLOC

/******************************
 * perl malloc uses the perl allocator; this was an idea 
 * to avoid memory contention issues that cropped up with the 
 * perl interface.  It doesn't help.
 * --CED 27-Oct-2004
 */
char *flux_perl_malloc(long size) {
  char *out;
  New(0,out,size,char);
  if(!out) {
    int i,j;
    fprintf(stderr,"Out of memory!\n");
    i=2; j=0;
    i /= j; /* throw arithmetic exception */
    exit(3498);
  }
  printf("Alloc'ed %d bytes (%x)\n",size,out);
  fflush(stdout);
  return out;
}

void flux_perl_free(void *foo) {
  /*  Safefree(foo);    do nothing for now */
}

#else
#ifdef USE_PADDED_MALLOC

/********** 
 * Padded malloc allocates more space but doesn't keep track of it.  That
 * seems to make the perl allocator happy, I'm not sure why. (allocation pattern
 * is similar to debugging malloc, but we don't bother with keeping track of the
 * fencing).  
 *  --CED 27-Oct-2004
 */
#define padding (500)
char *flux_padded_malloc(long size) {
  return malloc(size+2*padding)+padding;
}

void flux_padded_free( void *foo ) {
  free(foo-padding);
}
#endif  
#endif
#endif

