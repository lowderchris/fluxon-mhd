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
 * This is model.c version 1.1 - part of the FLUX 1.1 release.
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
 * world_update_ends
 * Updates the ends of each fluxon in turn...
 */
static long w_u_e_springboard(FLUXON *fl, int lab, int link, int depth) {
  fluxon_update_ends(fl);
  return 0;
}

void world_update_ends(WORLD *a) {
  tree_walker(a->lines,fl_lab_of, fl_all_ln_of, w_u_e_springboard,0);
}

/**********************************************************************
 * fluxon_update_ends
 * Checks and updates the magic boundary vertices at the end of a fluxon
 */
void fluxon_update_ends(FLUXON *f) {
   if(f->fc0->bound) {	
	(*(f->fc0->bound))(f->start);
   }
   if(f->fc1->bound) {
	(*(f->fc1->bound))(f->end);
   }
}  


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
  if(fl->fc0->world->verbosity) {
    printf(" %d",fl->label);
    fflush(stdout);
  }
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
 * from each cell.  Use fluxon_update_mag instead -- that handles
 * physics too.
 */
void fluxon_update_neighbors(FLUXON *fl, char global) {
  int i=0;
  VERTEX *v = fl->start;
  int verbosity = fl->fc0->world->verbosity;

  if(verbosity>=2) printf("fluxon_update_neighbors... (gl=%d), fluxon %d\n",global,fl->label);


  /* First thing - make sure the end conditions are up to date... */
  if(fl->fc0->bound) {	
    (*(fl->fc0->bound))(fl->start);
  }
  if(fl->fc1->bound) {
    (*(fl->fc1->bound))(fl->end);
  }

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

  /*
   * Set forces on the end vertex to zero, since it has no 
   * bend and no following segment and hence does not get 
   * looped over.
   */
  {
    VERTEX *v = fl->end;
    v->f_v[0] = v->f_v[1] = v->f_v[2] = 0;
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

     vertices = vertex_update_neighbors(v,global); 

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

  for(v=fl->start;v->next;v=v->next) {
    NUM f[3], fns, fnv, fn;
    NUM r = v->r_cl;
    NUM a;
    
    diff_3d(f,v->next->x,v->x);
    a = norm_3d(f);
    if(a<r && a>0)
      r=a;
    
    if(v->prev) {
      diff_3d(f, v->x, v->prev->x);
      a = norm_3d(f);
      if(a<r && a>0)
	r=a;

      sum_3d(f,v->prev->f_s,v->f_s);
      scale_3d(f,f,0.5);
    } else {
      f[0] = v->f_s[0];
      f[1] = v->f_s[1];
      f[2] = v->f_s[2];
    }
    
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
 *
 * Variables used for scaling:
 *     d  - harmonic mean of the various distances of neighbors and such
 *     b  - mean of the field strengths of the two segments on the vertex
 *     s  - stiffness coefficient from the force calculation
 */


/*** fastpow is quite fast for small integer powers, slightly slower for fractional ***/
NUM fastpow( NUM num, NUM exponent ) {
  NUM out = 1.0;
  NUM recip;
  
  if(exponent<0) {
    recip = 1.0/num;
    while(exponent<=-0.9999) {  exponent++; out *= recip;  }
  } else 
    while(exponent>=0.9999) { exponent--; out *= num; }

  if(exponent>0.0001 || exponent<-0.0001) { out *= pow(num,exponent);  }
  return out;
}
    
void fluxon_relax_step(FLUXON *f, NUM dt) {
  VERTEX *v = f->start;
  NUM a[3];
  NUM total[3];
  NUM force_factor, f_denom;
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

    f_denom = v->f_v_tot + 0.5 * (v->f_s_tot + v->prev->f_s_tot);

    force_factor = (f_denom == 0) ? 1 : ( norm_3d(v->f_t)  / (1e-9 + f_denom));
    if(force_factor == 0) force_factor = 1e-3;

    if(verbosity >= 3)  printf("fluxon %d, vertex %d: x=(%g,%g,%g).  v->r_cl=%g,  r_cl=%g,  force_factor = %g (%g / %g), f_t=(%g,%g,%g)[%g]\t",f->label,v->label, v->x[0],v->x[1],v->x[2], v->r_cl, r_cl, force_factor, norm_3d(v->f_t), f_denom, v->f_t[0],v->f_t[1],v->f_t[2],norm_3d(v->f_t));
    
    if(force_factor > 1.00001) {
      fprintf(stderr,"fluxon %d, vertex %d: force_factor = %g, >1!  This is allegedly impossible! You've got trouble, gov\n",f->label,v->label,force_factor);
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
      fac *= fastpow(force_factor,  w->step_scale.s_power);

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

    if(finite(a[0]) && finite(a[1]) &&finite(a[2])) 
      sum_3d(v->x,v->x,a);	     
    else
      if(verbosity >= 3) 
	printf("NON_FINITE OFFSET! f_s=(%g,%g,%g), f_v=(%g,%g,%g), f_t=(%g,%g,%g)",v->f_s[0],v->f_s[1],v->f_s[2],v->f_v[0],v->f_v[1],v->f_v[2],v->f_t[0],v->f_t[1],v->f_t[2]);

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


/*
 * snarf_filtered
 * Helper routine snarfs up a dumblist but skips over already-grabbed vertices and
 * image vertices.
 */
static inline void snarf_filtered(DUMBLIST *ws, DUMBLIST *dl, long passno) {
  int i;
  for(i=0;i<dl->n;i++) {
    VERTEX *V = (VERTEX *)(dl->stuff[i]);
    if( V->passno != passno && V->line ) {        
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
static DUMBLIST *snarfer_workspace;
static long line_snarfer(FLUXON *f, int lab_of, int ln_of, long depth) {
  VERTEX *v;
  for(v = f->start; v->next; v=v->next) {
    dumblist_quickadd(snarfer_workspace, v);
  }
  return 0;
}


/******/
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
      vv = ((VERTEX *)(nlist->stuff[j]));
      if(vv->line && vv->passno != passno) {
	vv->passno = passno;
	dumblist_quickadd(workspace, vv);
      }
    }

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
      

DUMBLIST *gather_neighbor_candidates(VERTEX *v, char global){
  static DUMBLIST *workspace =0;
  void **foo;
  int i;
  int n;
  int verbosity = v->line->fc0->world->verbosity;
  long passno = ++(v->line->fc0->world->passno);
  
  if(verbosity >= 2) {
    printf("passno=%d...  ",passno);
  }

  if(!v)        /* Paranoia */
    return;

  if(!workspace) 
    workspace = new_dumblist();

  workspace->n = 0;

  

  /**********************************************************************/
  /* Gather the candidates together 
   * This step has had a longish history -- I started out 
   * grabbing lots and lots of stuff, but "just" neighbors-of-neighbors
   * (and next and prev links) seems to be sufficient.  In pathological
   * cases, it might take a couple of timesteps to walk to the appropriate
   * new neighbor, but that doesn't seem to cause problems in practice.

  /* Snarf neighbors */
  if(!global) {
    VERTEX *v1;
    int ilen,iwid;

    dumblist_quickadd(workspace, v);

    if(verbosity>=3) {
      printf("initial seed - %d; ",workspace->n);
    }

    expand_lengthwise(workspace, 0, passno);
    ilen = workspace->n;

    if(verbosity >= 3) {
      printf("len - %d; ", workspace->n);
    }

    expand_via_neighbors(workspace, 0, passno);
    iwid = workspace->n;

    if(verbosity >= 3) {
      printf("wid - %d; ", workspace->n);
    }
    
    expand_lengthwise(workspace, ilen, passno);

    if(verbosity >= 3) {
      printf("len - %d; ", workspace->n);
    }

    expand_via_neighbors(workspace, iwid, passno);
    
    if(verbosity >= 3) {
      printf ("wid - %d; ", workspace->n);
    }

    expand_lengthwise(workspace, 0, passno);
    
    if(verbosity >= 3) {
      printf("final - %d\n",workspace->n);
    }
    
  }
    
   /* Incremental neighbor searching didn't work -- do a global 
     grab (ouch!).
  */
  if(global || (workspace->n == 0)) {
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

  /* If there's a photosphere, add the dummy photospheric mirror image */
  if(v->next) {
    PHOTOSPHERE *phot;
    PLANE *p;
    NUM a;
    if(v->line->fc0->world->photosphere.type) { /* assignment */
      /* Generate image and stuff it into the image point in the world
	 space */
      
      if(verbosity >= 4){
	printf("using photosphere (type is %d)...",v->line->fc0->world->photosphere.type);
	fflush(stdout);
      }
      phot = &(v->line->fc0->world->photosphere);
      p = phot->plane;
      switch(phot->type) {
	PLANE pl;
	POINT3D pt;
      case PHOT_CYL:
	/* Special case: mirror segment is reflected through a cylinder, radius p[0],
	 * aligned along the z axis.
	 */
	a = norm_2d(v->x);
	
	pl.origin[0] = v->x[0] * p->origin[0] / a;
	pl.origin[1] = v->x[1] * p->origin[0] / a;
	pl.origin[2] = 0;
	pl.normal[0] = v->x[0];
	pl.normal[1] = v->x[1];
	pl.normal[2] = 0;
	scale_3d(pl.normal, pl.normal, 1.0/norm_3d(pl.normal));
	reflect(v->line->fc0->world->image->x, v->x, &pl);
	
	a = norm_2d(v->next->x);
	pl.origin[0] = v->next->x[0] * p->origin[0] / a;
	pl.origin[1] = v->next->x[1] * p->origin[0] / a;
	pl.origin[2] = 0;
	pl.normal[0] = v->next->x[0];
	pl.normal[1] = v->next->x[1];
	pl.normal[2] = 0;
	scale_3d(pl.normal, pl.normal, 1.0/norm_3d(pl.normal));
	reflect(v->line->fc0->world->image->next->x, v->next->x, &pl);
	
	dumblist_quickadd(workspace, v->line->fc0->world->image);

	break;
      case PHOT_SPHERE:
	/***********
         * Spherical photosphere - sphere is located at the origin in the 
         * photospheric plane structuure; radius is normal[0]. */
	
	/* Construct the sub-vertex point on the sphere */
	diff_3d(&(pt[0]),      v->x, &(p->origin[0]));                 /* pt gets (x - origin) */
	scale_3d(&(pt[0]), &(pt[0]), p->normal[0]/norm_3d(&(pt[0])));  /* Scale to be on the sphere */
	sum_3d(&(pt[0]), &(pt[0]), p->origin);                    /* Put in real space */
	scale_3d(&(pt[0]), &(pt[0]), 2.0);                             
	diff_3d(v->line->fc0->world->image->x, &(pt[0]), v->x); /* Put reflection in image */

	/***** Do for the other point too *****/
	diff_3d(&(pt[0]), v->next->x, p->origin);
	scale_3d(&(pt[0]), &(pt[0]), p->normal[0]/norm_3d(&(pt[0])));
	sum_3d(&(pt[0]),&(pt[0]),p->origin);
	scale_3d(&(pt[0]), &(pt[0]), 2.0);
	diff_3d(v->line->fc0->world->image->next->x, &(pt[0]), v->x); 
	 
	dumblist_quickadd(workspace, v->line->fc0->world->image);

	break;
	  
      case PHOT_PLANE:
	reflect(v->line->fc0->world->image->x, v->x, p);
	reflect(v->line->fc0->world->image->next->x, v->next->x, p);
	dumblist_quickadd(workspace, v->line->fc0->world->image);
	break;
      default:
	fprintf(stderr,"Illegal photosphere type %d!\n",phot->type);
	exit(13);
      }
    }
  }

  /* Sort by vertex number, to avoid simple duplication */
  //  dumblist_sort(workspace, winnow_cmp_1);
  //if(verbosity >= 3) printf("gather_neighbor_candidates: sorting yielded %d elements\n",workspace->n);
  if(verbosity >= 3) printf("gather_neighbor_candidates: returning %d elements (may include dupes)\n",workspace->n);

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

  if(verbosity >= 1) {
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
 ** You feed in a VERTEX, and the curvature of the field line is
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
  if( ( (sincurve * 1.048 > curve_thresh_high) ||
	(curve_thresh_high > (M_PI * 1.0/3)  && 
	 ( (sincurve * 1.210 > curve_thresh_high) ||
	   (curve_thresh_high > (M_PI * 2.0/3) && 
	    ( (sincurve * (M_PI / 2)  > curve_thresh_high) ||
	      (coscurve < 0)
	      )
	    )
	   )
	 )
	)
	&& 
      /* This is the actual check but the trig operation is masked out */
      /* by the above mess for most cases.  That's an assignment.      */
      (curve = ATAN2(sincurve,coscurve)) > curve_thresh_high
      ) {

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
    return 2;

  } else {
    /******************************
     * Check for deletion
     */
    
    if( curve_thresh_low != 0 && 
	! ( (sincurve * 1.048 > curve_thresh_low) ||
	    (curve_thresh_low > (M_PI * 1.0/3)  && 
	     ( (sincurve * 1.210 > curve_thresh_low) ||
	       (curve_thresh_low > (M_PI * 2.0/3) && 
		( (sincurve * (M_PI / 2)  > curve_thresh_low) ||
		  (coscurve < 0)
		  )
		)
	       )
	     )
	    )
	&& 
	/* This is the actual check but the trig operation is masked out */
	/* by the above mess for most cases.  That's an assignment.      */
	(curve = ATAN2(sincurve,coscurve)) < curve_thresh_low
	) {
      NUM offset_dist;
      NUM pdist,ndist;
      /* Check maximum displacement of the fluxel if straightened */
      offset_dist = p_ls_dist(V->x,V->prev->x,V->next->x);
      offset_dist *= 1.333;
      if(offset_dist >= 0 && 
	 offset_dist < V->prev->r_cl && 
	 offset_dist < V->next->r_cl) {
	if(V->line->fc0->world->verbosity > 3){
	  printf("fix_curvature: unlinking %d from line %d\n",V->label,V->line->label);
	}
	delete_vertex(V);
	return -1;
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
 */
void reconnect_vertices( VERTEX *v1, VERTEX *v2 ) {
  VERTEX *vv;
  FLUXON *f1, *f2;
  FLUX_CONCENTRATION *fc;
  int i,j;

  if(!v1 || !v2 || !v1->next || !v2->next) {
    fprintf(stderr,"reconnect_vertices: error -- got a null or end vertex!\n");
    return;
  }

  f1 = v1->line;
  f2 = v2->line;
  
  if(f1==f2) {
    fprintf(stderr,"reconnect_vertices:  reconnection from %d to %d would create a plasmoid; not yet supported.  ignoring!\n",v1->label,v2->label);
    return;
  }
  
  /* Do the reconnection */
  vv = v1->next;
  v1->next = v2->next;
  v1->next->prev = v1;

  v2->next = vv;
  v2->next->prev = v2;

  /* Now clean up the fluxons and the back-to-fluxon links in the individual vertices. */
  
  /* f1 */
  i=1;
  for( vv = f1->start; vv->next; vv=vv->next ) {
    i++;
    vv->line = f1;
  }
  f1->end = vv;
  
  /* f2 */
  j=1;
  for( vv = f2->start; vv->next; vv=vv->next ) {
    i++;
    vv->line = f2;
  }
  f2->end = vv;
  
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
  
}

/******************************
 * v_recon_check: checks the current reconnection condition between a 
 * vertex and each of its neighbors.  Reconnection with an image charge
 * is not allowed!  The reconnection conditions are found via the rc_funcs table
 * in the world. 
 */
int vertex_recon_check( VERTEX *v1 ) {
  WORLD *w = v1->line->fc0->world;
  RC_FUNC **rcfuncp;
  VERTEX *rv = 0;
  int i;

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
    reconnect_vertices(v1,rv);
    return 1;
  }
  return 0;
}

/******************************
 * fluxon_recon_check
 * Iterates over all the vertices in a fluxon
 */
long fluxon_recon_check( FLUXON *f ) {
  VERTEX *v;
  int retval = 0;
  if( !(f->start) || !(f->start->next) || !(f->start->next->next) ) 
    return retval;

  for(v=f->start->next; v->next && v->next->next; v=v->next) {
    retval += vertex_recon_check(v);
  }
  return retval;
}

static long grc_acc;
static long grc_tramp(FLUXON *fl, int lab, int link, int depth) {
  grc_acc += fluxon_recon_check(fl);
}

long global_recon_check(WORLD *w) {
  grc_acc = 0;
  tree_walker( w->lines, fl_lab_of, fl_all_ln_of, grc_tramp, 0);
  return grc_acc;
}
    
/********************************************************************
 * fluxon end-condition handlers - update end position of the fluxon
 * to be consistent with the boundary condition
 */

void fl_b_start_open(VERTEX *v) {
  POINT3D a;
  if(!v->next) {
    fprintf(stderr,"HEY! fl_b_start_open got an end vertex!  Doing nothing...\n");
    return;
  }
  if(v->prev) {
    fprintf(stderr,"HEY! fl_b_start_open got a middle vertex!  Doing nothing...\n");
    return;
  } 
  
  diff_3d( a,    v->next->x, v->line->fc0->x );
  scale_3d(a,    a,          v->line->fc0->locale_radius / norm_3d(a) );
  sum_3d(  v->x, a,          v->line->fc0->x );
}

void fl_b_end_open(VERTEX *v) {
  POINT3D a;
  if(!v->prev) {
    fprintf(stderr,"HEY! fl_b_end_open got a beginning vertex! Doing nothing...\n");
    return;
  }
  if(v->next) {
    fprintf(stderr,"HEY! fl_b_end_open got a middle vertex!  Doing nothing...\n");
    return;
  }
  
  diff_3d( a,     v->prev->x,   v->line->fc1->x );
  scale_3d(a,     a,            v->line->fc1->locale_radius / norm_3d(a) );
  sum_3d(  v->x,  a,            v->line->fc1->x );
} 

void fl_b_start_plasmoid(VERTEX *v) {
  if(!v->next) {
    fprintf(stderr,"HEY! fl_b_start_plasmoid got an end vertex! Doing nothing...\n");
    return;
  }
  if(v->prev) {
    fprintf(stderr,"HEY! fl_b_start_plasmoid got a middle vertex! Doing nothing...\n");
    return;
  }
  if(!v->next->next) {
    fprintf(stderr,"HEY! plasmoid conditions require more than one middle vertex! Doing nothing...\n");
    return;
  }
  
  cp_3d( v->x, v->line->end->prev->x );
}

void fl_b_end_plasmoid(VERTEX *v) {
  if(!v->prev) {
    fprintf(stderr,"HEY! fl_b_end_plasmoid got a start vertex! Doing nothing...\n");
    return;
  }
  if(v->next) {
    fprintf(stderr,"HEY! fl_b_end_plasmoid got a middle vertex! Doing nothing...\n");
    return;
  }
  if(!v->prev->prev) {
    fprintf(stderr,"HEY! plasmoid conditions require real fluxons, you git!...\n");
    return;
  }
  
  cp_3d( v->x, v->line->start->next->x );
}

