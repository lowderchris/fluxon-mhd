/* model.c
 * 
 * Slightly higher-level data manipulation routines
 * These routines generally use both data.c and geometry.c.  They are
 * the machinery of the model itself, doing things like handling 
 * vertex neighbor location and trimming.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Southwest Research Institute, 2004
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
 */
#include "model.h"
#include <stdio.h>
#include <math.h>

#include <unistd.h> /* for isatty -- delete if needed */
extern FILE *gl_outf; /* for debugging -- delete */

#include "data.h"
#include "geometry.h"
#include "io.h"


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
 * world_update_neighbors
 * 
 * Calls fluxon_update_neighbors to process the whole world!
 */
static WORLD *gl_a;  /* springboard world parameter */
static char gl_gl;   /* springboard global flag */
static void ((**gl_f_funcs)()); /* springboard functions list */
static NUM *gl_minmax; /* Keeps track of minimum and maximum forces */

static long w_u_n_springboard(FLUXON *fl, int lab, int link, int depth) {
  if(fl->fc0->world->verbosity >= 5) dump_all_fluxon_tree(gl_a->lines);
  fluxon_update_neighbors(fl, gl_gl);
  return 0;
}

void world_update_neighbors(WORLD *a, char global) {
  gl_a = a;
  gl_gl = global;
  tree_walker(a->lines, fl_lab_of, fl_all_ln_of, w_u_n_springboard, 0);
}

/**********************************************************************
 * world_update_mag updates the magnetic forces for the whole world 
 */
static long w_u_m_springboard(FLUXON *fl, int lab, int link, int depth) {
  gl_minmax = fluxon_update_mag(fl,gl_gl, gl_f_funcs, gl_minmax);
  return 0;
}

NUM *world_update_mag(WORLD *a, char global) {
  gl_a = a;
  gl_gl = global;
  gl_f_funcs = a->f_funcs;
  gl_minmax = 0;
  tree_walker(a->lines, fl_lab_of, fl_all_ln_of, w_u_m_springboard, 0);
  return gl_minmax;
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
static long w_r_s_springboard(FLUXON *fl, int lab, int link, int depth) {
  fluxon_relax_step(fl, gl_t);
  return 0;
}

void world_relax_step(WORLD *a, NUM t) {
  gl_t=t;
  tree_walker(a->lines, fl_lab_of, fl_all_ln_of, w_r_s_springboard, 0);
}

/**********************************************************************
 * fluxon_update_neighbors
 * 
 * Calls vertex_update_neighbors to process a whole fluxon. 
 * This is wasteful, as it throws away the Voronoi vertex information
 * from each cell.  Use fluxon_update instead -- that handles
 * physics too.
 */
void fluxon_update_neighbors(FLUXON *fl, char global) {
  int i=0;
  VERTEX *v = fl->start;
  int verbosity = fl->fc0->world->verbosity;

  if(verbosity>=2) printf("fluxon_update_neighbors... (gl=%d), fluxon %d\n",global,fl->label);

  while(v->next) {
    if(verbosity>=3)  printf("\tfluxon_update_neighbors... vertex %d\n",v->label);
    vertex_update_neighbors(v,global);
    if(verbosity>=4)   fdump_fluxon(stdout,fl,0);
    if(verbosity>=3)   printf("=============\n\n");
    v=v->next;
    i++;
  }
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

NUM *fluxon_update_mag(FLUXON *fl, char global, void ((**f_funcs)()), NUM *minmax) {
  static NUM output[4];
  NUM *f_min, *f_max, *fr_min, *fr_max;
  int i=0;
  int verbosity = fl->fc0->world->verbosity;
  VERTEX *v = fl->start;

  if(!minmax) {
    minmax = output;
    minmax[0] = minmax[1] = minmax[2] = minmax[3] = -1;
  }

  f_min = &minmax[0];
  f_max = &minmax[1];
  fr_min = &minmax[2];
  fr_max = &minmax[3];
  
  if(verbosity >= 2)  printf("fluxon_update_mag (fluxon %d): %c",fl->label,(verbosity==2?' ':'\n'));

  for(i=0, v=fl->start; v->next; v=v->next, i++) {
    void (**f_func)();
    NUM r;
    int ii;
    HULL_VERTEX *vertices;
    if(verbosity >= 3) { printf("V%4d ",v->label); fflush(stdout); }

    /* Update neighbor map and establish neighbors' (r,a) variables, 
     *  which will be used by the physics funcs!  */
    vertices = vertex_update_neighbors(v,global);

    /* Zero force accumulators */
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


    if(fl->fc0->world->verbosity >= 3) 
      printf("---forces:\n");

     /* Accumulate forces and relevant lengthscales */
    for(f_func = &f_funcs[0]; *f_func; f_func++) 
      (**f_func)(v,vertices);


    if(fl->fc0->world->verbosity >= 3)
      printf("\n");

  }
  if(verbosity >= 3) printf("\n");


  /* Find closest approach radius, and calculate total forces */
  
  if(verbosity >= 2) printf("Radius - fluxon %d: %c",fl->label,(verbosity==2?' ':'\n'));

  for(v=fl->start->next;v->next;v=v->next) {
    NUM f[3], fns, fnv, fn;
    NUM r = v->r_cl;
    NUM a;
    
    if(v->prev) {
      diff_3d(f, v->x, v->prev->x);
      a = norm_3d(f);
      if(a<r && a>0)
	r=a;
    }

    if(v->next) {
      diff_3d(f,v->next->x,v->x);
      a = norm_3d(f);
      if(a<r && a>0)
	r=a;
    }

    sum_3d(f,v->prev->f_s,v->f_s);
    scale_3d(f,f,0.5);

    fns = (v->r_s > 0) ? (norm_3d(f)) : 0;
    fnv = (v->r_v > 0) ? (norm_3d(v->f_v)) : 0;
    
    sum_3d(v->f_t,f,v->f_v);

    fn = norm_3d(v->f_t);

    /* Check max, min, etc. */
    if(*f_min < 0 || fn < *f_min)
      *f_min = fn;
    if(*f_max < 0 || fn > *f_max)
      *f_max = fn;
    if(*fr_min < 0 || fn/r < *fr_min)
      *fr_min = fn/r;
    if(*fr_max < 0 || fn/r > *fr_max)
      *fr_max = fn/r;

    if(verbosity >= 3){ printf("V%4d: r=%g, fmax=%g, frmax=%g ",v->label, r, *f_max, *fr_max); fflush(stdout);}
  }

  if(verbosity >= 3)  printf("\n");
  return minmax;
}

/**********************************************************************
 * fluxon_relax_step
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
 */
void fluxon_relax_step(FLUXON *f, NUM t0) {
  NUM t = t0;
  VERTEX *v = f->start;
  NUM a[3];
  NUM total[3];
  NUM force_factor, f_denom;
  int verbosity = f->fc0->world->verbosity;

  a[2] = a[1] = a[0] = 0;


  for(v=v->next;v && v->next; v=v->next) {
    NUM r_cl= v->r_cl;
    NUM foo[3];
    NUM d,d1;
    /* Take a step proportional to the force and the harmonic mean of the fluxon length and local approach distance. */
    
    diff_3d(foo,v->next->x,v->x);
    d = norm_3d(foo);
    
    diff_3d(foo,v->x,v->prev->x);
    d1 = norm_3d(foo);

    
    if(r_cl <= 0) {
      fprintf(stderr,"ASSERTION FAILED!  Negative distance %g on vertex %d!\n",r_cl,v->label);
      fprintf(stderr,"vertex has %d neighbors\n",v->neighbors.n);
      r_cl = 1; /* hope the problem goes away */
    }
    
    
    d = 4/(1/d + 1/d1 + 1/r_cl + (v->prev ? 1/v->prev->r_cl : 1/r_cl));

    f_denom = v->f_v_tot + 0.5 * (v->f_s_tot + v->prev->f_s_tot);
    force_factor = (f_denom == 0) ? 1 : ( norm_3d(v->f_t)  / f_denom);
    if(force_factor == 0) force_factor = 1e-3;

    if(verbosity >= 3)  printf("fluxon %d, vertex %d: x=(%g,%g,%g).  v->r_cl=%g,  r_cl=%g,  force_factor = %g (%g / %g), f_t=(%g,%g,%g)[%g]\t",f->label,v->label, v->x[0],v->x[1],v->x[2], v->r_cl, r_cl, force_factor, norm_3d(v->f_t), f_denom, v->f_t[0],v->f_t[1],v->f_t[2],norm_3d(v->f_t));
    
    if(force_factor > 1.00001) {
      fprintf(stderr,"fluxon %d, vertex %d: force_factor = %g, >1!  This is allegedly impossible! You've got trouble, gov\n",f->label,v->label,force_factor);
      fflush(stdout);
      fflush(stderr);
    }
    

    //    {
    //    NUM f2 = 2/(1+1/(force_factor));
    //    scale_3d(a,v->f_t, t * r_cl * f2 * f2 / force_factor);
    //    }

    //    scale_3d(a,v->f_t, t * r_cl * 1/(1/sqrt(force_factor) + 1/force_factor));
    
    /* `forces' are forces per unit length, so multiply times d to bring them down to normal-force level. */
    /* Throw in another factor of d to account for the fact that you have to take bigger steps, not just */
    /* equal-sized steps, where d is large. */
    scale_3d(a,v->f_t, t * d * d * force_factor );
    //scale_3d(a,v->f_t, t * r_cl * r_cl * force_factor);

    sum_3d(v->x,v->x,a);	     

    if(verbosity >= 3)    printf("after update: x=(%g,%g,%g)\n",v->x[0],v->x[1],v->x[2]);
  }

  /* Update of start and end positions goes here! */

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
 */

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
  dl = gather_neighbor_candidates(v,global);


  if(verbosity >= 3) {
    printf("vertex_update_neighbors: Vertex %d has %d candidates: ",v->label,dl->n);
    if(verbosity >= 4)
      for(i=0;i<dl->n;i++)
	printf("  %d",((VERTEX *)((dl->stuff)[i]))->label);
    printf("\n");
  }


  /* This winnow_neighbor_candidates call isn't strictly necessary --
     hull_neighbors should do the job.  But winnow_neighbor_candidates 
     or something like it should be called eventually to do walkalongs.
     (comment left in as a placeholder!) */
  /*winnow_neighbor_candidates(v,dl);*/
  
  hv = hull_neighbors(v, dl); /* save hv for return */

    if(verbosity >= 3) {
      printf("Hull_neighbors returned %d neighbors: ",dl->n);
      if(verbosity >= 4)
	for(i=0;i<dl->n;i++) 
	  printf("  %d",((VERTEX *)((dl->stuff)[i]))->label);
      printf("\n");
    }

  /* Walk through both dumblists, updating the neighbors' 
     "nearby" information as needed. */
  for(i=j=0; i<dl->n || j<vn->n;) {
    if( (j >= vn->n) || (i < dl->n && (dl->stuff[i] < vn->stuff[j]) ) ) {
      if(verbosity >= 6) {
	printf ("[i=%d,j=%d,nbor:%x(%d)]  ",i,j,(dl->stuff[i]),((VERTEX *)(dl->stuff[i]))->label);
	fflush(stdout);
      }
      dumblist_add( &( ((VERTEX *)(dl->stuff[i]))->nearby ), v );
      i++;
      cpflag=1;
    } else if( (i >= dl->n) || (j < vn->n && dl->stuff[i] > vn->stuff[j] )) {
      dumblist_delete( &( ((VERTEX *)(vn->stuff[j]))->nearby), v);
      cpflag=1;
      j++;
    } else /* if(dl->stuff[i] == vn->stuff[j]) */ {
      i++;
      j++;
    }
  }
  
  /* Copy the neighbor lists, if necessary */
  if(cpflag) {
    vn->n = 0;
    dumblist_snarf(vn,dl);
  }
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
 * The "global" flag indicates whether global search (SLOOOW) is required --
 * good for initial conditions or if a pathology is detected.
 */

/* snarf_list
 * Helper routine for snarfing up a set of vertices and its next and prev
 * elements.  You feed in a list of vertices (e.g. the "neighbors" list
 * from a vertex) and a workspace, and the vertices that are neighbors and 
 * nearby to each of the elements of the list are themselves added to the 
 * workspace.
 *
 */

static inline void snarf_list(DUMBLIST *workspace, VERTEX **foo, int n) {
  int i,j;
  for(i=0;i<n;i++) {
    dumblist_add(workspace,(void *)foo[i]);
    dumblist_snarf(workspace,&(foo[i]->neighbors));
    dumblist_snarf(workspace,&(foo[i]->nearby));
  }

  /* Purge image pseudo-vertices from the snarfed list */
  for(i=j=0;i<workspace->n;i++) {
    if( (((VERTEX **)(workspace->stuff))[i])->line ) {
      if(j!=i)
	workspace->stuff[j] = workspace->stuff[i];
      j++;
    }
  }
  workspace->n = j;

}


/* expand_list
 * Helper routine for expanding a workspace list out:  you feed in 
 * a workspace with a list of vertices, and it gets expanded to include
 * the next and previous item on each of those vertex lists.
 */
static inline void expand_list(DUMBLIST *workspace) {
  int i,n;
  n = workspace->n;
  for(i=0;i<n;i++) {
    VERTEX *v = ((VERTEX **)(workspace->stuff))[i];
    if(v) {
      if( v->next )
	dumblist_quickadd(workspace, v->next);
      if( v->prev ) 
	dumblist_quickadd(workspace, v->prev);
    }
  }
}

/* Helper routine for snarfing up whole fluxons -- called via
   tree_walk, in data.c */
static DUMBLIST *snarfer_workspace;
static long line_snarfer(FLUXON *f, int lab_of, int ln_of, long depth) {
  VERTEX *v;
  for(v = f->start; v->next; v=v->next) {
    dumblist_quickadd(snarfer_workspace, v);
  }
  return 0;
}

DUMBLIST *gather_neighbor_candidates(VERTEX *v, char global){
  static DUMBLIST *workspace =0,*workspace2 = 0;
  void **foo;
  int i;
  int n;
  int verbosity = v->line->fc0->world->verbosity;

  if(!v)        /* Paranoia */
    return;

  if(!workspace) 
    workspace = new_dumblist();

  if(!workspace2)
    workspace2 = new_dumblist();

  workspace->n = workspace2->n = 0;

  /**********************************************************************/
  /* Gather the candidates together */

  /* Snarf neighbors */

  if(!global) {
    VERTEX *v1;

    /* Grab neighbors & nearby from vertex & its siblings */
    snarf_list(workspace2,&v,1); 

    //    if(v->next)     snarf_list(workspace2,&(v->next),1);
    //    if(v->prev)     snarf_list(workspace2,&(v->prev),1);

    /* Grab siblings of all the neighbors & nearby */
    expand_list(workspace2);

    /* Grab neighbors & nearby of all *those* guys */
    snarf_list(workspace,(VERTEX **)workspace2->stuff,workspace2->n);

    /* Remove duplicates to save time in the next step */
    dumblist_sort(workspace,winnow_cmp_1);

    //    /* Grab siblings of all those guys */
    //    expand_list(workspace);

  }


  /* Incremental neighbor searching didn't work -- do a global 
     grab (ouch!).
  */
  if(global || (workspace->n == 0)) {
    FLUXON *f;

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

  /* If there's a photosphere, add the dummy photospheric mirror image */
  if(v->next) {
    PLANE *phot;
    if(phot = v->line->fc0->world->photosphere) { /* assignment */
      /* Generate image and stuff it into the image point in the world
	 space */
      reflect(v->line->fc0->world->image->x, v->x, phot);
      reflect(v->line->fc0->world->image->next->x, v->next->x, phot);
      dumblist_add(workspace, v->line->fc0->world->image);
    }
  }

  /* Sort by vertex number, to avoid simple duplication */
  dumblist_sort(workspace, winnow_cmp_1);

  if(verbosity >= 3) printf("gather_neighbor_candidates: sorting yielded %d elements\n",workspace->n);
  return workspace;
}


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
      d_prev = fl_segment_dist(v,vn->prev);
      /* Step forward or backward if necessary.  Fork if the field line
       * bends toward us.
       */
      if( (d_next > 0) && (d_next < d0) ) {
	if( (d_prev > 0) && (d_prev < d_next)) 
	  dumblist_quickadd(horde,vn->prev);
	(*vnp) = vn->next;
      } else if( (d_prev > 0) && (d_prev < d0)) {
	(*vnp) = vn->prev;
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
 * planar hull information and a host of other useful 2-D projected goodies.
 *
 * Each of the neighbor candidate segments is projected into the current segment's 
 * plane.  Then the Voronoi cell is constructed around the central point in the
 * projection plane.  Fortunately, routines in geometry.c take care of the
 * geometric details.
 *
 * As long as a hull calculation is being done anyhow, hull_neigbors
 * keeps the valuable hull vertex information and returns it in 
 * a buffer of (x,y) pairs.  The buffer is ephemeral -- it gets overwritten
 * each time hull_neighbors is called -- so use it while you can.
 */

HULL_VERTEX *hull_neighbors(VERTEX *v, DUMBLIST *horde) {
  int i;
  int a,b;
  POINT3D x0;
  int verbosity = v->line->fc0->world->verbosity;

  static HULL_VERTEX *voronoi_buf = 0;
  static int voronoi_bufsiz = 0;
  static DUMBLIST *rejects = 0;

  if(!v->next) {
    fprintf(stderr,"pick_neighbors: assertion failed -- vertex should have a next member!\n\tGiving up...\n");
    return 0;
  }

  /* Project the horde into the plane perpendicular to v's line segment.
   * The projected vectors go into the vertices' scratch space.    
   */

  project_n_fill(v, horde); /* in geometry.c */

  /* Grow the buffer if necessary. */
  if(voronoi_bufsiz < horde->n) {
    voronoi_bufsiz = horde->n*2;

    voronoi_buf = (HULL_VERTEX *)realloc(voronoi_buf, voronoi_bufsiz*sizeof(HULL_VERTEX));

    if(!voronoi_buf) {
      fprintf(stderr,"Couldn't get memory in hull_neighbors!\n");
      exit(99);
    } 

  }


  
  if(verbosity >= 5) {
    printf("V%4d: (%7.3g, %7.3g, %7.3g) -- (%7.3g, %7.3g, %7.3g)\n",v->label,v->x[0],v->x[1],v->x[2], v->next->x[0],v->next->x[1],v->next->x[2]);
    for(i=0;i<horde->n;i++) {
      NUM p0[3], p1[3];
      fl_segment_deluxe_dist(p0, p1, v, ((VERTEX *)(horde->stuff[i])));

      printf("\tV%4d:  3d=%g\tproj dist=%7.3g (%7.3g, %7.3g, %7.3g) -- (%7.3g, %7.3g, %7.3g).  Closest approach: (%7.3g, %7.3g, %7.3g) -- (%7.3g, %7.3g, %7.3g)\n"
	     ,((VERTEX *)(horde->stuff[i]))->label
	     ,((VERTEX *)(horde->stuff[i]))->r_cl
	     ,((VERTEX *)(horde->stuff[i]))->r
	     ,((VERTEX *)(horde->stuff[i]))->x[0]
	     ,((VERTEX *)(horde->stuff[i]))->x[1]
	     ,((VERTEX *)(horde->stuff[i]))->x[2]
	     ,((VERTEX *)(horde->stuff[i]))->next->x[0]
	     ,((VERTEX *)(horde->stuff[i]))->next->x[1]
	     ,((VERTEX *)(horde->stuff[i]))->next->x[2]
	     ,p0[0],p0[1],p0[2]
	     ,p1[0],p1[1],p1[2]
	     );
    }
  }


  /* Find the 2-D hull.  For now, ask for rejects in order to plot 'em. */
  hull_2d(voronoi_buf,horde,0);

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
  tree_walker(w->lines,fl_lab_of,fl_all_ln_of,gfp_tramp);
  return sc_acc;
}
  
  
/**********************************************************************  
 ** fix_curvature
 ** fluxon_fix_curvature
 ** global_fix_curvature
 ** 
 ** You feed in a VERTEX, and the curvature of the field line
 ** is compared with the curvature threshold that you
 ** pass in (in radians).  When fix_curvature returns, the 
 ** VERTEX you passed in may have been split, so that the 
 ** new internal angle is less than your specified threshold
 ** everywhere.  If you are stepping along a fluxon fixing 
 ** curvatures, you should store V->next BEFORE calling fix_curvature, 
 ** or you'll have to step through exactly half of any daughter VERTEXes 
 ** that are created.
 ** 
 ** Set 0 curvature to get 0.05 radian (~3 degrees) for splitting 
 ** and 1 degree for merging).  
 **
 ** For now, no merging is done!
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
int fix_curvature(VERTEX *V, NUM curve_thresh) {
  int i;
  NUM d1[3];
  NUM d2[3];
  NUM cr[3];
  NUM sincurve,coscurve,curve;
  NUM scale;
  int ret=0;

  if(curve_thresh==0)
    curve_thresh = 0.05;

  if(!V || !V->next || !V->prev || !V->line) 
    return 0;


  diff_3d(d1,V->prev->x,V->x); /* d1 = segment A */
  diff_3d(d2,V->x,V->next->x); /* d2 = segment B */
  
  cross(cr,d1,d2);

  scale = 1.0/sqrt(norm2_3d(d1) * norm2_3d(d2));

  sincurve = norm_3d(cr) * scale;
  coscurve = inner_3d(d1,d2) * scale;
	
  /* Do a bilinear, Q&D check to see if splitting is necessary -- this */
  /* avoids having to do any trig for most non-splitting cases! */
  /* This slows down a bit for curvature thresholds above 30 degrees, but */
  /* adds no overhead for small curvature thresholds.  theta/sin(theta) is */
  /* 1.047 for 30 degrees, 1.209 for 60 degrees, and pi/2 for 90 degrees -- */
  /* hence the magic numbers.  */
  if( ( (sincurve * 1.048 > curve_thresh) ||
	(curve_thresh > (M_PI * 1.0/3)  && 
	 ( (sincurve * 1.210 > curve_thresh) ||
	   (curve_thresh > (M_PI * 2.0/3) && 
	    ( (sincurve * (M_PI / 2)  > curve_thresh) ||
	      (coscurve < 0)
	      )
	    )
	   )
	 )
	)
	&& 
      /* This is the actual check but the trig operation is masked out */
      /* by the above mess for most cases.  That's an assignment.      */
      (curve = atan2(sincurve,coscurve)) > curve_thresh
      ) {
    VERTEX *Vnew;
    NUM P[3], P0[3],P1[3];

    /* Find location of new prior vertex: 1/3 along the previous vec */
    sum_3d(P0, V->prev->x, V->prev->x);
    sum_3d(P0, P0, V->x);
    scale_3d(P0, P0, 1.0/3);

    /* Find location of new next vertex: 2/3 along the current vec */
    sum_3d(P1, V->next->x, V->next->x);
    sum_3d(P1, P1, V->x);
    scale_3d(P1, P1, 1.0/3);

    /* Find centroid of triangle */
    sum_3d(P,V->prev->x, V->next->x);
    sum_3d(P, P, V->x);
    scale_3d(P,P,1.0/3);

    /* Find location of new current vertex: halfway between P and V->x */
    sum_3d(P, P, V->x);
    scale_3d(V->x, P, 0.5);

    /* Make the new prior vertex, and link */
    add_vertex_after(V->line, V->prev,
		     new_vertex(0, P0[0],P0[1],P0[2], V->line));
    add_vertex_after(V->line, V,
		     new_vertex(0, P1[0],P1[1],P1[2], V->line));
    return 1;

  }
  
  return 0;
}

int fluxon_fix_curvature(FLUXON *f, NUM curve_thresh) {
  int ret=0;
  VERTEX *V = f->start;
  int verbosity = f->fc0->world->verbosity;

  if(verbosity >= 2) printf("fluxon_fix_curvature: fluxon %d\n",f->label);
  
  while(V && V != f->end) {
    VERTEX *Vnext = V->next;
    if(verbosity >= 3) printf("  vertex %d\n",V->label);
    ret += fix_curvature(V,curve_thresh);
    V=Vnext;
  }
  return ret;
}

/* tree walker hooks for global_fix_curvature... */
static NUM cu_thr;
static int cu_acc;
static long cu_tramp(FLUXON *fl, int lab, int link, int depth) {
  cu_acc += fluxon_fix_curvature(fl, cu_thr);
  return 0;
}
int global_fix_curvature(WORLD *w, NUM curv_thresh) {
  cu_thr = curv_thresh;
  cu_acc = 0;
  tree_walker(w->lines,fl_lab_of,fl_all_ln_of,cu_tramp);
  return cu_acc;
}
