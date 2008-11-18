/*
 * physics.c
 * 
 * Contains force laws and other directly physics-related 
 * calculations.  The 'force laws' manipulate vertex-local 
 * quantities; you may include additional local simulation steps
 * by adding additional routines that carry them out. 
 *
 * The force laws actually calculate force per unit length along
 * the fluxon; the relaxation step should multiply the f/l times
 * the length of each segment.
 * 
 * Quasi-local simulation steps such as reconnection that require
 * physical quantities in the neighbors will require a second
 * pass through the data (to update all local calculations such as
 * magnetic field), but can also be implemented using the same 
 * force-law structure.
 *
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Craig DeForest, 2004-2007.
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
 * This file is part of the FLUX 2.0 release (31-Oct-2007).
 * 
 */

#include <math.h>
#include "physics.h"
#include "geometry.h"

char *code_info_physics="%%%FILE%%%";

/* global force name table */
struct FLUX_FORCES FLUX_FORCES[] = {
  {"b_simple","inverse-area B field (breaks for open cells)",b_simple},
  {"b_eqa","2004 angular equipartion B calculation (B required at start of list)",b_eqa},
  {"e_simple", "B energy per vertex using simple inverse-area (breaks for open cells)",e_simple},
  {"e_simple2", "B energy per vertex using simple inverse-area ( doesnt break for open cells)",e_simple2},
  {"e_open", "B energy per vertex with angular equipartition, includes open cells",e_open},
  {"e_eqa","B energy per vertex using angular equipartition",e_eqa},
  {"f_curvature","(OLD) simple curvature force (B-normalized)",f_curvature},
  {"f_curvature2","(OLD) simple curvature force (B-normalized) with stiffener",f_curvature2},
  {"f_curvature3","simple curvature force (B-normalized) with stiffener and distance scaling",f_curvature3},
  {"f_pressure_equi","(OLD) 2001 pressure law (B-normalized)",f_pressure_equi},
  {"f_pressure_equi2","(OLD) 2001 pressure law (B-normalized)",f_pressure_equi2},
  {"f_pressure_equi2a","(OLD) 2001 pressure law (B-normalized; patched for open field)",f_pressure_equi2a},
  {"f_pressure_equi2b","(OLD) 2001 pressure law (B-normalized; patched for open field)",f_pressure_equi2b},
  {"f_curv_hm","Curvature force law (harmonic mean curvature)",f_curv_hm},
  {"f_curv_m","Curvature force law (mean curvature)",f_curv_m},
  {"f_p_eqa_radial","Angular equipartition pressure, radial forces",f_p_eqa_radial},
  {"f_vertex","Vertex distribution pseudo-force (DEPRECATED)",f_vertex},
  {"f_vertex2","Vertex distribution pseudo-force",f_vertex2},
  {"f_vertex3","Vertex distribution pseudo-force",f_vertex3},
  {"f_vertex4","Vertex distribution pseudo-force (r^2 repulsion)",f_vertex4},
  {"f_vertex5","Vertex distribution pseudo-force",f_vertex5},
  {"f_vert","Vertex distribution pseudo-force",f_vert},
  //  {"m_hydrostatic","Quasi-hydrostatic mass loading with N_0 and T_0 from powers of B_base; requires a b_",m_hydrostatic},
  {0,0,0}
};

char *force_ptr_to_str(void *f) {
  int i;
  for(i=0;*(FLUX_FORCES[i].name); i++) 
    if(FLUX_FORCES[i].func == f)
      return FLUX_FORCES[i].name;
  return 0;
}

void *force_str_to_ptr(char *s) {
  int i;
  if(!s)
    return 0;
  for(i=0; *(FLUX_FORCES[i].name); i++) 
    if( !strcmp(FLUX_FORCES[i].name, s) )
      return FLUX_FORCES[i].func;
  return 0;
}

/* global reconnection-condition calculator table */
struct FLUX_RECON FLUX_RECON[] = {
  {"rc_a_ad2","Threshold angle per d^2 (J proxy avoids explicit B calc)", rc_a_ad2, "(min. angle), (min A/D^2 [or 0 for none])"},
  {"rc_a_ad2_h","Similar to rc_a_ad2, but J threshold adjusts with a scale height", rc_a_ad2_h,"(min. angle),(min A/ln(z/h)/d^2), (h)"},
  {"","",0,""}
};


char *recon_ptr_to_str(void *f) {
  int i;
  for(i=0;*(FLUX_RECON[i].name); i++) 
    if(FLUX_RECON[i].func == f)
      return FLUX_RECON[i].name;
  return 0;
}

void *recon_str_to_ptr(char *s) {
  int i;
  if(!s)
    return 0;
  for(i=0; *(FLUX_RECON[i].name); i++) 
    if( !strcmp(FLUX_RECON[i].name, s) )
      return FLUX_RECON[i].func;
  return 0;
}


/* fluxon physics routines get called like this:
 *    force = routine( VERTEX *v );
 * where, on return, force is 0 (error) or a pointer to a static area 
 * containing the force 3-vector.
 * 
 * Later, there will be other routines for more complex models
 * handled along the fluxon. 
 */

/**********************************************************************
 **********************************************************************
 * The older forces used for field line manipulation had |B| divided out of
 * them, to simplify their calculation.   They are kept here, warts and
 * all, for reference.  But they are deprecated.  So instead of 
 * f_curvature and f_pressure_equi, use f_curv and f_p_eqa
 *
 *

/**********************************************************************
 * f_curvature
 * Curvature `force':   (B^ . \del) B^ 
 *      = b^_x ( db^_x / dx ) + b^_y ( db^_y / dy ) + b^_z ( db^_z / dz )
 * 
 * but b^ is just l^.
 *
 * The curvature force is a vertex-centered force.  There's actually 
 * a little too much information:  not only do we get an angle and a 
 * length, we get an extra length too (there are *two* legs out of the
 * segment, not just one, and they have differing lengths).  I use 
 * the harmonic mean because that's dominated by the smaller of the
 * two lengths rather than by the larger.
 * 
 */
void f_curvature(VERTEX *V, HULL_VERTEX *verts) {
  NUM recip_l1, recip_l2, recip_len, len;
  NUM b1hat[3];
  NUM b2hat[3];
  NUM curve[3];
  NUM force[3];


  if(!V->next || !V->prev) 
    return;

  diff_3d(b2hat,V->next->x,V->x);
  scale_3d(b2hat,b2hat, (recip_l2 = 1.0 / norm_3d(b2hat)));

  diff_3d(b1hat,V->x,V->prev->x);
  scale_3d(b1hat,b1hat, (recip_l1 = 1.0 /norm_3d(b1hat)));

  recip_len =  ( recip_l1 + recip_l2 ) * 0.5;
  
  diff_3d(curve,b2hat,b1hat);
  scale_3d(curve,curve,recip_len); /* convert to force per unit length */

  sum_3d(V->f_v,V->f_v,curve);

  V->f_v_tot += norm_3d(curve);


}

/**********************************************************************
 * f_curvature
 * Curvature `force':   (B^ . \del) B^ 
 *      = b^_x ( db^_x / dx ) + b^_y ( db^_y / dy ) + b^_z ( db^_z / dz )
 * 
 * but b^ is just l^.
 *
 * The curvature force is a vertex-centered force.  There's actually 
 * a little too much information:  not only do we get an angle and a 
 * length, we get an extra length too (there are *two* legs out of the
 * segment, not just one, and they have differing lengths).  I use 
 * the harmonic mean because that's dominated by the smaller of the
 * two lengths rather than by the larger.
 *
 * This version 2 of the force law goes to infinity as the angle goes to 
 * 180 degrees -- hopefully the extra stiffness far from the small angle
 * condition will prevent certain kinds of runaway.  The stiffening factor
 * is secant(theta/2), scaled down by a factor of 2.  
 * 
 * Stiffening fudge factor at selected angles:
 *   15 deg - 1.008
 *   30 deg - 1.036
 *   45 deg - 1.086
 *   60 deg - 1.167
 *   75 deg - 1.294
 *   90 deg - 1.500
 *  105 deg - 1.849
 *  120 deg - 2.500
 *  135 deg - 3.914
 *  150 deg - 7.965
 *  165 deg -29.85 
 *  
 */
void f_curvature2(VERTEX *V, HULL_VERTEX *verts) {
  NUM recip_l1, recip_l2, recip_len, len;
  NUM b1hat[3];
  NUM b2hat[3];
  NUM curve[3];
  NUM force[3];
  // NUM cross[3];
  // NUM sin2;
  NUM sec;
  NUM fudge;

  if(!V->next || !V->prev) 
    return;

  diff_3d(b2hat,V->next->x,V->x);
  scale_3d(b2hat,b2hat, (recip_l2 = 1.0 / norm_3d(b2hat)));

  diff_3d(b1hat,V->x,V->prev->x);
  scale_3d(b1hat,b1hat, (recip_l1 = 1.0 /norm_3d(b1hat)));

  recip_len =  ( recip_l1 + recip_l2 ) * 0.5;
  

  // Calculate fudge factor -- near unity where we like the angles to be, 
  // but extra powerful near 180 degrees.  
  sec = 2.0/(1 + inner_3d(b1hat, b2hat)); // secant(theta/2)
  fudge = (sec+1)/2;
  
  diff_3d(curve,b2hat,b1hat);
  scale_3d(curve,curve,recip_len * fudge); // force/length, with fudge
  
  sum_3d(V->f_v,V->f_v,curve);

  V->f_v_tot += norm_3d(curve);


}

void f_curvature3(VERTEX *V, HULL_VERTEX *verts) {
  NUM recip_l1, recip_l2, recip_len, len;
  NUM b1hat[3];
  NUM b2hat[3];
  NUM curve[3];
  NUM force[3];
  // NUM cross[3];
  // NUM sin2;
  NUM sec;
  NUM fudge;
  NUM dist_recip;

  if(!V->next || !V->prev) 
    return;

  diff_3d(b2hat,V->next->x,V->x);
  scale_3d(b2hat,b2hat, (recip_l2 = 1.0 / norm_3d(b2hat)));

  diff_3d(b1hat,V->x,V->prev->x);
  scale_3d(b1hat,b1hat, (recip_l1 = 1.0 /norm_3d(b1hat)));

  recip_len =  ( recip_l1 + recip_l2 ) * 0.5;
  

  // Calculate fudge factor -- near unity where we like the angles to be, 
  // but extra powerful near 180 degrees.  
  sec = 2.0/(1 + inner_3d(b1hat, b2hat)); // secant(theta/2)
  fudge = (sec+1)/2;
  
  diff_3d(curve,b2hat,b1hat);
  
  // Proximity
  dist_recip = 0.5 * (1/V->r_cl + 1/V->prev->r_cl);

  scale_3d(curve,curve,recip_len * fudge * dist_recip ); // force/length, with fudge
  sum_3d(V->f_v,V->f_v,curve);
  V->f_v_tot += norm_3d(curve);


}
  
/**********************************************************************
 * f_pressure_equi
 * Magnetic pressure `force' : - \del_perp(|B|) / |B|
 *
 * Magnetic pressure force, using the equipartition method 
 * which distributes a given fluxon's flux equally by delta-phi.
 * This method distributes current throughout each fluxon's cross-section.
 * 
 * This uses the Kankelborg/DeForest discretization:
 * del_perp(|B|) / |B| ~ \sum_i ( n^_i delta\phi_i / (pi r_i) )
 * where n^_i is the vector from the current point to the ith neighbor,
 * delta\phi is the angular extent of the ith line segment in the 
 * Voronoi hull, and r_i is the radius to the closest approach of the
 * ith line segment (see DeForest notebook V, pp. 69-70)
 */
void f_pressure_equi(VERTEX *V, HULL_VERTEX *verts) {
  NUM force[3]; /* force accumulator */
  static HULL_VERTEX *point_cache = 0;  /* Workspace grows and sticks */
  static int pc_size = 0;               /* Workspace size */
  int i;
  NUM pmatrix[9];
  NUM len;

  force[0] = force[1] = force[2] = 0;

  if(!V->next) 
    return;

  /* Get the perpendicular-plane projection matrix */
  /* (This is wasteful -- it's already generated in geometry.c.  Is
     there an easy way to send it back?) */
  projmatrix(pmatrix, V->x, V->next->x);

  /*printf("f_pressure_equi:  V%4d,  n=%d\n",V->label, V->neighbors.n);*/

  for(i=0;i<V->neighbors.n;i++) {
    NUM f_i[3];
    NUM phi_l, phi_r, deltaphi, f;
    HULL_VERTEX *left, *right;
    VERTEX *N = (VERTEX *)V->neighbors.stuff[i];

    left = &verts[i];
    right = (i==0) ? &verts[V->neighbors.n - 1] : &verts[i-1];
    
    /* Calculate delta-phi by calculating the maximal phis and subtracting. 
     * In the open case, the endpoint phi is the angle at which the 
     * perpendicular bisector to r_i goes away.  In the closed case, 
     * it's the angle to the vertex.
     */

    deltaphi = (left->a_l - right->a_r);
    TRIM_ANGLE(deltaphi);
    // printf("vertex %d (%d): deltaphi=%g\t",V->label,N->label,deltaphi*180/3.14159);

    if(deltaphi < -EPSILON) {
      fprintf(stderr,"Assertion failed!  deltaphi <0 in f_pressure_equi, vertex %d, neighbor %d (%18.12g deg); correcting to %18.12g\n",V->label, i, deltaphi*180/M_PI,deltaphi*180/M_PI + 360);
      deltaphi += M_PI+M_PI;

      printf("(Vertex %d --  neighbors are:\n",V->label);
      {
	int ii;
	for(ii=0;ii<V->neighbors.n;ii++) {
	  VERTEX *vv= (VERTEX *)(V->neighbors.stuff[ii]);
	  printf("\t%d -\tx=(%g,%g,%g),\tprojx = (%g,%g), a=%g,\tr=%g, lefthull a_l=%g, righthull a_r=%g, open: l=%d, r=%d\n",
		 vv->label,
		 vv->x[0],vv->x[1],vv->x[2],
		 vv->scr[0],vv->scr[1],
		 vv->a * 180/3.1415926, vv->r,
		 verts[ii].a_l * 180/3.1415926, verts[ ii==0 ? V->neighbors.n-1 : ii-1 ].a_r * 180/3.1415926,
		 verts[ii].open, verts[ ii==0 ? V->neighbors.n-1 : ii-1 ].open
		 );
	  printf("\t\thull[%d]: p=(%g,%g)\n",ii,verts[ii].p[0],verts[ii].p[1]);
	}
	printf(")\n");
      }
    }
    /*    fprintf(stderr," VERTEX #%d, i=%d: deltaphi=%g\tr=%g\n",N->label,i,deltaphi,N->r);*/


    f = deltaphi / M_PI / N->r;

    N->scr[2] = 0; /* Make sure we're operating in the perp. plane */
    
    vec_mmult_3d(f_i,pmatrix,N->scr);       /* Put back into 3-D */
    /*    printf("  in plane: (%g,%g,%g); projects back to (%g,%g,%g)\n",N->scr[0],N->scr[1],N->scr[2],f_i[0],f_i[1],f_i[2]);*/

    {
      /* Scale to the proper force */
      NUM a = norm_3d(f_i);
      if(a==0) {
	fprintf(stderr,"Assertion failed in physics.c!  f_i has 0 norm!\n");
	exit(99);
      }
      scale_3d(f_i,f_i, - f / a);  /* Scale to the calc. force */
    }

    sum_3d(force,force,f_i);
    V->f_s_tot += fabs(f);
  }


  sum_3d(V->f_s,V->f_s,force);

  /*printf("f_pressure_equi -- force is (%8g,%8g,%8g), total is (%8g,%8g,%8g)\n\n",force[0],force[1],force[2],V->f_s[0],V->f_s[1],V->f_s[2]);*/

  {
    NUM r;
    int ii;
    /* Accumulate closest neighbor */
    /* This code will need updating if fluxons get variable flux! */
    r = -1;
    for(ii=0;ii<V->neighbors.n;ii++) {
      NUM a;
      a = ((VERTEX *)(V->neighbors.stuff[ii]))->r;
      if(r<0 || a<r)
	r = a;
    }
    V->r =  r / 2;  /* If all fluxons have the same amount of flux, then the min.
		       radius is always half of the closest neighbor's distance. */
  }
  
  return;
}
/**********************************************************************
 * f_pressure_equi2
 * Magnetic pressure `force' : - \del_perp(|B|) / |B|
 * Corrects for the factor-of-two error in f_pressure_equi
 *
 * Magnetic pressure force, using the equipartition method 
 * which distributes a given fluxon's flux equally by delta-phi.
 * This method distributes current throughout each fluxon's cross-section.
 * 
 * This uses the Kankelborg/DeForest discretization:
 * del_perp(|B|) / |B| ~ \sum_i ( n^_i delta\phi_i / (pi r_i) )
 * where n^_i is the vector from the current point to the ith neighbor,
 * delta\phi is the angular extent of the ith line segment in the 
 * Voronoi hull, and r_i is the radius to the closest approach of the
 * ith line segment (see DeForest notebook V, pp. 69-70)
 */
void f_pressure_equi2(VERTEX *V, HULL_VERTEX *verts) {
  NUM force[3]; /* force accumulator */
  static HULL_VERTEX *point_cache = 0;  /* Workspace grows and sticks */
  static int pc_size = 0;               /* Workspace size */
  int i;
  NUM pmatrix[9];
  NUM len;

  force[0] = force[1] = force[2] = 0;

  if(!V->next) 
    return;

  /* Get the perpendicular-plane projection matrix */
  /* (This is wasteful -- it's already generated in geometry.c.  Is
     there an easy way to send it back?) */
  projmatrix(pmatrix, V->x, V->next->x);

  /*printf("f_pressure_equi:  V%4d,  n=%d\n",V->label, V->neighbors.n);*/

  for(i=0;i<V->neighbors.n;i++) {
    NUM f_i[3];
    NUM phi_l, phi_r, deltaphi, f;
    HULL_VERTEX *left, *right;
    VERTEX *N = (VERTEX *)V->neighbors.stuff[i];

    left = &verts[i];
    right = (i==0) ? &verts[V->neighbors.n - 1] : &verts[i-1];
    
    /* Calculate delta-phi by calculating the maximal phis and subtracting. 
     * In the open case, the endpoint phi is the angle at which the 
     * perpendicular bisector to r_i goes away.  In the closed case, 
     * it's the angle to the vertex.
     */

    deltaphi = (left->a_l - right->a_r);
    TRIM_ANGLE(deltaphi);
    //    printf("vertex %d (%d): deltaphi=%g\t",V->label,N->label,deltaphi*180/3.14159);

    if(deltaphi < -EPSILON) {
      //      fprintf(stderr,"Assertion failed!  deltaphi <0 in f_pressure_equi (%18.12g deg); correcting to %18.12g (vertex %d, neighbor %d, hull on left: %g, hull on right: %g)\n",deltaphi*180/M_PI,deltaphi*180/M_PI + 360,V->label,i, left->a_l*180/3.14159, right->a_r*180/3.14159);
      deltaphi += M_PI+M_PI;
#ifdef NEVER
	{
	  int j,k;
	  VERTEX *vn;
	  for(j=0;j<V->neighbors.n;j++) {
	    vn = ((VERTEX *)(V->neighbors.stuff[j]));
	    fprintf(stderr,"Neighbor #%d: label is %d, a is %g, r is %g, hull angles: %g, %g\n",j,vn->label,vn->a*180/3.1415926,vn->r, verts[j].a_l*180/3.14159,verts[(j+1) % V->neighbors.n].a_r*180/3.14159);
	    //    for(k=0;k<vn->nearby.n;k++) {
	    //  fprintf(stderr,"  nearby #%d: %d\n",k,((VERTEX *)(vn->nearby.stuff[k]))->label);
	    //}
	  }
	}
#endif
    }
    
    /*    fprintf(stderr," VERTEX #%d, i=%d: deltaphi=%g\tr=%g\n",N->label,i,deltaphi,N->r);*/


    f = 2 * deltaphi / ( M_PI * N->r );

    N->scr[2] = 0; /* Make sure we're operating in the perp. plane */
    
    vec_mmult_3d(f_i,pmatrix,N->scr);       /* Put back into 3-D */
    /*    printf("  in plane: (%g,%g,%g); projects back to (%g,%g,%g)\n",N->scr[0],N->scr[1],N->scr[2],f_i[0],f_i[1],f_i[2]);*/

    {
      /* Scale to the proper f  orce */
      NUM a = norm_3d(f_i);
      if(a==0) {
	fprintf(stderr,"Assertion failed in physics.c!  f_i has 0 norm (vertex %d on %d (next is %d, prev is %d), neighbor %d on %d (next is %d, prev is %d), i=%d, V->neighbors.n=%d)\n",V->label,V->line->label, V->next?V->next->label:0, V->prev?V->prev->label:0, N->label, N->line->label , N->next?N->next->label:0, N->prev?N->prev->label:0,V->neighbors.n);
	{
	  int j,k;
	  VERTEX *vn;
	  for(j=0;j<V->neighbors.n;j++) {
	    vn = ((VERTEX *)(V->neighbors.stuff[j]));
	    fprintf(stderr,"Neighbor #%d: label is %d\n",j,vn->label);
	    //	    for(k=0;k<vn->nearby.n;k++) {
	    //  fprintf(stderr,"  nearby #%d: %d\n",k,((VERTEX *)(vn->nearby.stuff[k]))->label);
	    //}
	  }
	}
	exit(99);
      }
      scale_3d(f_i,f_i, - f / a);  /* Scale to the calc. force */
    }

    sum_3d(force,force,f_i);
    V->f_s_tot += fabs(f);
  }


  sum_3d(V->f_s,V->f_s,force);

  /*printf("f_pressure_equi -- force is (%8g,%8g,%8g), total is (%8g,%8g,%8g)\n\n",force[0],force[1],force[2],V->f_s[0],V->f_s[1],V->f_s[2]);*/

  {
    NUM r;
    int ii;
    /* Accumulate closest neighbor */
    /* This code will need updating if fluxons get variable flux! */
    r = -1;
    for(ii=0;ii<V->neighbors.n;ii++) {
      NUM a;
      a = ((VERTEX *)(V->neighbors.stuff[ii]))->r;
      if(r<0 || a<r)
	r = a;
    }
    V->r =  r / 2;  /* If all fluxons have the same amount of flux, then the min.
		       radius is always half of the closest neighbor's distance. */
  }
  
  return;
}
/**********************************************************************
 * f_pressure_equi2a
 * Magnetic pressure `force' : - \del_perp(|B|) / |B|
 * Corrects for the factor-of-two error in f_pressure_equi
 *
 * Magnetic pressure force, using the equipartition method 
 * which distributes a given fluxon's flux equally by delta-phi.
 * This method distributes current throughout each fluxon's cross-section.
 * 
 * This uses the Kankelborg/DeForest discretization:
 * del_perp(|B|) / |B| ~ \sum_i ( n^_i delta\phi_i / (pi r_i) )
 * where n^_i is the vector from the current point to the ith neighbor,
 * delta\phi is the angular extent of the ith line segment in the 
 * Voronoi hull, and r_i is the radius to the closest approach of the
 * ith line segment (see DeForest notebook V, pp. 69-70)
 *
 * This version applies the force from the median angle of the cell wall,
 * rather than from the normal.
 */
void f_pressure_equi2a(VERTEX *V, HULL_VERTEX *verts) {
  NUM force[3]; /* force accumulator */
  static HULL_VERTEX *point_cache = 0;  /* Workspace grows and sticks */
  static int pc_size = 0;               /* Workspace size */
  int i;
  NUM pmatrix[9];
  NUM len;

  force[0] = force[1] = force[2] = 0;

  if(!V->next) 
    return;

  /* Get the perpendicular-plane projection matrix */
  /* (This is wasteful -- it's already generated in geometry.c.  Is
     there an easy way to send it back?) */
  projmatrix(pmatrix, V->x, V->next->x);

  /*printf("f_pressure_equi:  V%4d,  n=%d\n",V->label, V->neighbors.n);*/

  for(i=0;i<V->neighbors.n;i++) {
    NUM f_i[3];
    NUM f_j[3];
    NUM phi_l, phi_r, deltaphi, fn, fp, phi;
    NUM nhat_x, nhat_y;
    HULL_VERTEX *left, *right;
    VERTEX *N = (VERTEX *)V->neighbors.stuff[i];

    left = &verts[i];
    right = (i==0) ? &verts[V->neighbors.n - 1] : &verts[i-1];
    
    deltaphi = left->a_l - right->a_r;
    TRIM_ANGLE(deltaphi);
    //    printf("vertex %d (%d): deltaphi=%g\t",V->label,N->label,deltaphi*180/3.14159);
    
    if(deltaphi < -EPSILON) {
      fprintf(stderr,"Assertion failed!  deltaphi <0 in f_pressure_equi (%18.12g deg); correcting to %18.12g\n",deltaphi*180/M_PI,deltaphi*180/M_PI + 360);
      deltaphi += M_PI+M_PI;
    }
    
    /*    fprintf(stderr," VERTEX #%d, i=%d: deltaphi=%g\tr=%g\n",N->label,i,deltaphi,N->r);*/


    fn = 2 * deltaphi / (M_PI * N->r);
    fp = 2 * (cos(right->a_r - N->a) - cos(left->a_l - N->a)) / (M_PI * N->r);
    {
      NUM r = norm_2d(N->scr);
      nhat_x = N->scr[0]/r;
      nhat_y = N->scr[1]/r;
    }
    f_i[0] = - nhat_x * fn   +   nhat_y * fp ;
    f_i[1] = - nhat_y * fn   -   nhat_x * fp ;

    vec_mmult_3d( f_j, pmatrix, f_i );      /* Put back into 3-D */

    sum_3d(force,force,f_j);
    V->f_s_tot += sqrt( fn*fn + fp*fp );
  }

  sum_3d(V->f_s,V->f_s,force);

  /*printf("f_pressure_equi -- force is (%8g,%8g,%8g), total is (%8g,%8g,%8g)\n\n",force[0],force[1],force[2],V->f_s[0],V->f_s[1],V->f_s[2]);*/

  {
    NUM r;
    int ii;
    /* Accumulate closest neighbor */
    /* This code will need updating if fluxons get variable flux! */
    r = -1;
    for(ii=0;ii<V->neighbors.n;ii++) {
      NUM a;
      a = ((VERTEX *)(V->neighbors.stuff[ii]))->r;
      if(r<0 || a<r)
	r = a;
    }
    V->r =  r / 2;  /* If all fluxons have the same amount of flux, then the min.
		       radius is always half of the closest neighbor's distance. */
  }
  
  return;
}

/**********************************************************************
 * f_pressure_equi2b
 * Magnetic pressure `force' : - \del_perp(|B|) / |B|
 * Corrects for the factor-of-two error in f_pressure_equi
 *
 * Magnetic pressure force, using the equipartition method which
 * distributes a given fluxon's flux equally by delta-phi.  This
 * method distributes current throughout each fluxon's cross-section.
 * 
 * This uses the Kankelborg/DeForest discretization: del_perp(|B|) /
 * |B| ~ \sum_i ( n^_i delta\phi_i / (pi r_i) ) where n^_i is the
 * vector from the current point to the ith neighbor, delta\phi is the
 * angular extent of the ith line segment in the Voronoi hull, and r_i
 * is the radius to the closest approach of the ith line segment (see
 * DeForest notebook V, pp. 69-70)
 *
 * This version applies the force from the median angle of the cell
 * wall, rather than from the normal.
 * 
 * In this version only the normal wall force, and not the parallel
 * wall force, is doubled.
 */

void f_pressure_equi2b(VERTEX *V, HULL_VERTEX *verts) {
  NUM force[3]; /* force accumulator */
  static HULL_VERTEX *point_cache = 0;  /* Workspace grows and sticks */
  static int pc_size = 0;               /* Workspace size */
  int i;
  NUM pmatrix[9];
  NUM len;

  force[0] = force[1] = force[2] = 0;

  if(!V->next) 
    return;

  /* Get the perpendicular-plane projection matrix (This is wasteful
     -- it's already generated in geometry.c.  Is there an easy way to
     send it back?) */

  projmatrix(pmatrix, V->x, V->next->x);

  /*printf("f_pressure_equi: V%4d, n=%d\n",V->label,
    V->neighbors.n);*/

  for(i=0;i<V->neighbors.n;i++) {
    NUM f_i[3];
    NUM f_j[3];
    NUM phi_l, phi_r, deltaphi, fn, fp, phi;
    NUM nhat_x, nhat_y;
    HULL_VERTEX *left, *right;
    VERTEX *N = (VERTEX *)V->neighbors.stuff[i];

    left = &verts[i];
    right = (i==0) ? &verts[V->neighbors.n - 1] : &verts[i-1];
    
    deltaphi = left->a_l - right->a_r;

    TRIM_ANGLE(deltaphi); //normalize to between -PI and PI
    
    if(deltaphi < -EPSILON) {
      fprintf(stderr,"Assertion failed!  deltaphi <0 in f_pressure_equi (%18.12g deg; vertex %d; neighbor %d); correcting to %18.12g\n",deltaphi*180/M_PI,V->label,((VERTEX *)(V->neighbors.stuff[i]))->label,deltaphi*180/M_PI + 360);
      deltaphi += M_PI+M_PI;
    }
    
    fn = 2 * deltaphi / (M_PI * N->r); //N->r is radius to neighbor
    fp = (cos(right->a_r - N->a) - cos(left->a_l - N->a)) / (M_PI * N->r);
    {
      NUM r = norm_2d(N->scr);
      nhat_x = N->scr[0]/r;
      nhat_y = N->scr[1]/r;
    }
    f_i[0] = - nhat_x * fn   +   nhat_y * fp ;
    f_i[1] = - nhat_y * fn   -   nhat_x * fp ;

    vec_mmult_3d( f_j, pmatrix, f_i );      /* Put back into 3-D */

    sum_3d(force,force,f_j);
    V->f_s_tot += sqrt( fn*fn + fp*fp );
  }

  sum_3d(V->f_s,V->f_s,force);

  /*printf("f_pressure_equi -- force is (%8g,%8g,%8g), total is
    (%8g,%8g,%8g)\n\n",force[0],force[1],force[2],V->f_s[0],V->f_s[1],V->f_s[2]);*/

  {
    NUM r;
    int ii;
    /* Accumulate closest neighbor */
    /* This code will need updating if fluxons get variable flux! */
    r = -1;
    for(ii=0;ii<V->neighbors.n;ii++) {
      NUM a;
      a = ((VERTEX *)(V->neighbors.stuff[ii]))->r;
      if(r<0 || a<r)
	r = a;
    }
    V->r =  r / 2;  /* If all fluxons have the same amount of flux,
		       then the min.  radius is always half of the
		       closest neighbor's distance. */
  }
  
  return;
}

/**********************************************************************
 * f_vert
 * Pseudoforce that moves vertices along field lines.  Used to keep them 
 * nicely evened out along the lines and to attract them to 
 * curvature.
 * 
 * Strictly speaking, this is two forces:  a vertex repulsive force
 * and another force that attracts vertices toward curvature.  But
 * the two should always be used together, so they're included together.
 *
 * f_vert is optimized for use with the newer forces that are straight
 * force-per-unit-length: it is artificially strengthened by the 
 * local B field magnitude.
 * 
 */
void f_vert(VERTEX *V, HULL_VERTEX *verts) {
  NUM force[3];
  NUM d1[3], d2[3];
  NUM d1nr,d2nr,fn,d1n,d2n;
  NUM r_clp, r_cln;
  NUM Bmag;

  /* Exclude endpoints */
  if(!V->next || !V->prev)
    return;

  /* Scale forces to the maximum of the B values in our vicinity */
  Bmag = (V->b_mag > V->prev->b_mag) ? V->b_mag : V->prev->b_mag ;

  /* Repulsive force from nearest neighbors.  The force drops as 
   *  1/r.
   */
  diff_3d(d1,V->x,V->prev->x);
  scale_3d(d1, d1, (d1nr = 1.0 / norm_3d(d1))); /* assignment */

  diff_3d(d2,V->next->x, V->x);
  scale_3d(d2, d2, (d2nr = 1.0 / norm_3d(d2)));  /* assignment */
  
  fn = (d1nr - d2nr);
  V->f_v_tot += fabs(fn*Bmag);

  /* Proximity-attractive force.  This attracts vertices toward places
   * where field lines are interacting. 
   */
  r_clp = 1.0 / V->prev->r_cl;
  r_cln = 1.0 / V->r_cl;
  
  V->f_v_tot += fabs(0.5 * (r_cln - r_clp) * Bmag);
  fn += 0.5 * (r_cln - r_clp);
  
  
  /* Generate a vector along the field line and scale it to the 
   * calculated force.  
   *
   * Finally, stick the force where it belongs in the VERTEX's force vector.
   */
  
  sum_3d(force,d1, d2);
  scale_3d(force, force, 0.5 * fn * Bmag / norm_3d(force) );
  
  sum_3d(V->f_v, V->f_v, force);

  return;
}

/**********************************************************************
 * f_vertex
 * Pseudoforce that moves vertices along field lines.  Used to keep them 
 * nicely evened out along the lines and to attract them to 
 * curvature.
 * 
 * Strictly speaking, this is two forces:  a vertex repulsive force
 * and another force that attracts vertices toward curvature.  But
 * the two should always be used together, so they're included together.
 *
 * f_vertex is optimized for the older forces (f_pressure_equi and f_curvature)
 * that are b-normalized.  Use f_vert with the newer ones.
 * 
 */
void f_vertex(VERTEX *V, HULL_VERTEX *verts) {
  NUM force[3];
  NUM d1[3], d2[3];
  NUM d1nr,d2nr,fn,d1n,d2n;
  NUM l1, l2;
  /* Exclude endpoints */
  if(!V->next || !V->prev)
    return;


  /* Repulsive force from nearest neighbors.  The force drops as 
   *  1/r.
   */
  diff_3d(d1,V->x,V->prev->x);
  scale_3d(d1, d1, (d1nr = 1.0 / (l1=norm_3d(d1)))); /* assignment */

  diff_3d(d2,V->next->x, V->x);
  scale_3d(d2, d2, (d2nr = 1.0 / (l2=norm_3d(d2))));  /* assignment */
  
  fn = (d1nr - d2nr);
  V->f_v_tot += fabs(fn);

  /* Proximity-attractive force.  This attracts vertices toward places
   * where field lines are interacting. 
   */
    if(V->prev && V->next) {
      NUM r_clp, r_cln;
  
      r_clp = 2.0 / V->prev->r_cl;
      r_cln = 2.0 / V->r_cl;
  
      V->f_v_tot += fabs(0.5 * (r_cln - r_clp));
      fn += 0.5 * (r_cln - r_clp);
    }
    

  /* Generate a vector along the field line and scale it to the 
   * calculated force.  
   *
   * Finally, stick the force where it belongs in the VERTEX's force vector.
   */

  sum_3d(force,d1, d2);
  scale_3d(force, force, 0.5 * fn / norm_3d(force) );

  /*printf("f_vertex: V%4d, force=(%g,%g,%g)\n",V->label,force[0],force[1],force[2]);*/

  sum_3d(V->f_v, V->f_v, force);

  return;
}
 


 /**********************************************************************
 * f_vertex2
 * Pseudoforce that moves vertices along field lines.  Used to keep them 
 * nicely evened out along the lines and to attract them to 
 * curvature.
 * 
 * Strictly speaking, this is threeo forces: a vertex repulsive force,
 * a force that attracts vertices toward curvature, and one that
 * attracts vertices towards a high density of fluxons.  But the three
 * should always be used together, so they're included together.
 *
 * f_vertex2 is optimized for the older forces (f_pressure_equi and
 * f_curvature) that are b-normalized.  Use f_vert with the newer
 * ones. This force is unitless.
 * 
 */

void f_vertex2(VERTEX *V, HULL_VERTEX *verts) {
  NUM force[3];
  NUM d1[3], d1n [3], d2[3], d2n[3];
  NUM d1nr,d2nr,fn1,fn2,fn3; /* fn's are scalers */
  NUM l1, l2, l_fac;
  NUM dpp[3], dnn[3];
  NUM alpha_p, alpha_n, lnn, lpp, i1; 
  /* Exclude endpoints */
  if(!V->next || !V->prev)
    return;


  /* Repulsive force from nearest neighbors.  The force drops as 
   *  1/r.
   */

  diff_3d(d1,V->x,V->prev->x);
  scale_3d(d1n, d1, (d1nr = 1.0 / (l1=norm_3d(d1)))); /* assignment */

  diff_3d(d2,V->next->x, V->x);
  scale_3d(d2n, d2, (d2nr = 1.0 / (l2=norm_3d(d2)))); /* assignment */
  
  l_fac = 3. / ( d1nr + d2nr );

  fn1 = (d1nr - d2nr) * l_fac;
  V->f_v_tot += fabs(fn1); 

  /* Proximity-attractive force.  This attracts vertices toward places
   * where field lines are interacting.
   */

  if(V->prev && V->next) {
    NUM r_clp, r_cln;
  
    r_clp = 2.0 / V->prev->r_cl;
    r_cln = 2.0 / V->r_cl;
  
    fn2 = 0.2 * (r_cln - r_clp) * l_fac; 
    V->f_v_tot += fabs(fn2);
      
  } else {

    fn2 = 0.0;

  }
    
  /* Curvature-attractive force.  This attracts vertices toward
   * curvature so that sharp angles attract vertices to smooth
   * themselves out.
   */

  if(V->prev->prev && V->next->next) {

    diff_3d(dnn,V->next->next->x,V->next->x);
    lnn = norm_3d(dnn);

    diff_3d(dpp,V->prev->x,V->prev->prev->x);
    lpp = norm_3d(dpp);

    alpha_p = M_PI - acos( (inner_3d(dpp, d1)) / (l1 * lpp) );
    alpha_n = M_PI - acos( (inner_3d(dnn, d2)) / (l2 * lnn) );
    
    fn3 = 5. * (alpha_p - alpha_n);
    V->f_v_tot += fabs(fn3);
     
  } else {

    fn3 = 0.0;
  
  }


  /* Generate a unit vector along the field line and scale it to the
   * calculated force.
   */

  sum_3d(force,d1n, d2n);
  scale_3d(force, force, 0.5 * (fn1 + fn2 + fn3) / norm_3d(force) );

  /* Stick the force where it belongs in the VERTEX's force vector.*/

  sum_3d(V->f_v, V->f_v, force);

  //printf("f_vertex2: V%4d, f=(%g,%g,%g), |f|=%g \n",V->label,force[0],force[1],force[2], norm_3d(force) );

  return;
}
 

 /**********************************************************************
 * f_vertex3
 * Pseudoforce that moves vertices along field lines.  Used to keep them 
 * nicely evened out along the lines and to attract them to 
 * curvature.
 * 
 * Strictly speaking, this is threeo forces: a vertex repulsive force,
 * a force that attracts vertices toward curvature, and one that
 * attracts vertices towards a high density of fluxons.  But the three
 * should always be used together, so they're included together.
 *
 * f_vertex2 is optimized for the older forces (f_pressure_equi and
 * f_curvature) that are b-normalized.  Use f_vert with the newer
 * ones. This force is unitless.
 * 
 */

void f_vertex3(VERTEX *V, HULL_VERTEX *verts) {
  NUM force[3];
  NUM d1[3], d1n [3], d2[3], d2n[3];
  NUM d1nr,d2nr,fn1,fn2,fn3; /* fn's are scalers */
  NUM l1, l2, l_fac;
  NUM dpp[3], dnn[3];
  NUM alpha_p, alpha_n, lnn, lpp, i1; 
  /* Exclude endpoints */
  if(!V->next || !V->prev)
    return;


  /* Repulsive force from nearest neighbors.  The force drops as 
   *  1/r.
   */

  diff_3d(d1,V->x,V->prev->x);
  scale_3d(d1n, d1, (d1nr = 1.0 / (l1=norm_3d(d1)))); /* assignment */

  diff_3d(d2,V->next->x, V->x);
  scale_3d(d2n, d2, (d2nr = 1.0 / (l2=norm_3d(d2)))); /* assignment */
  
  l_fac = 3. / ( d1nr + d2nr );

  fn1 = (d1nr - d2nr) * l_fac;
  V->f_v_tot += fabs(fn1); 

  fn2 = 0.0;
    
  /* Curvature-attractive force.  This attracts vertices toward
   * curvature so that sharp angles attract vertices to smooth
   * themselves out.
   */

  if(V->prev->prev && V->next->next) {

    diff_3d(dnn,V->next->next->x,V->next->x);
    lnn = norm_3d(dnn);

    diff_3d(dpp,V->prev->x,V->prev->prev->x);
    lpp = norm_3d(dpp);

    alpha_p = M_PI - acos( (inner_3d(dpp, d1)) / (l1 * lpp) );
    alpha_n = M_PI - acos( (inner_3d(dnn, d2)) / (l2 * lnn) );
    
    fn3 = 5. * (alpha_p - alpha_n);
    V->f_v_tot += fabs(fn3);
     
  } else {

    fn3 = 0.0;
  
  }


  /* Generate a unit vector along the field line and scale it to the
   * calculated force.
   */

  sum_3d(force,d1n, d2n);
  scale_3d(force, force, 0.25 * (fn1 + fn2 + fn3) / norm_3d(force)  );

  /* Stick the force where it belongs in the VERTEX's force vector.*/

  sum_3d(V->f_v, V->f_v, force);

  //printf("f_vertex2: V%4d, f=(%g,%g,%g), |f|=%g \n",V->label,force[0],force[1],force[2], norm_3d(force) );

  return;
}
 

/**********************************************************************
 * f_vertex4 
 *
 * Use f_vertex4 with the newer forces ones. This force has units of
 * 1/length. TThis one doesn't have a attraction to fluxon proximity
 * of other fluxons force, but it does have repulsion from vertices
 * (could stand to be a little stronger?) and attraction to curvature.
 * 
 */

void f_vertex4(VERTEX *V, HULL_VERTEX *verts) {
  NUM force[3];
  NUM d1[3], d1n [3], d2[3], d2n[3];
  NUM d1nr,d2nr,fn1,fn2,fn3; /* fn's are scalers */
  NUM l1, l2, l_fac;
  NUM dpp[3], dnn[3];
  NUM alpha_p, alpha_n, lnn, lpp, i1; 
  /* Exclude endpoints and image charges */
  if(!V->next || !V->prev || V->label < 0)
    return;


  // Find the normal vectors pointing along the previous and next segments,
  // and collect the lengths of the original vectors...

  diff_3d(d1,V->x,V->prev->x);
  scale_3d(d1n, d1, (d1nr = 1.0 / (l1=norm_3d(d1)))); /* assignment */

  diff_3d(d2,V->next->x, V->x);
  scale_3d(d2n, d2, (d2nr = 1.0 / (l2=norm_3d(d2)))); /* assignment */


  /* Repulsive force from nearest neighbors.  The force drops as 
   *  1/r^2 but is normalized l to yield something like 1.
   */
  fn1 = (d1nr*d1nr-d2nr*d2nr) * ( (l1+l2)*0.5 ) * ( (l1+l2)*0.5 );
  //fn1 *= 0.5;

  /* Curvature-attractive force.  This attracts vertices toward
   * curvature so that sharp angles attract vertices to smooth
   * themselves out.
   */

  if(V->prev->prev && V->next->next) {

    diff_3d(dnn,V->next->next->x,V->next->x);
    lnn = norm_3d(dnn);

    diff_3d(dpp,V->prev->x,V->prev->prev->x);
    lpp = norm_3d(dpp);

    alpha_p = M_PI - acos( (1-1e-6) * (inner_3d(dpp, d1)) / (l1 * lpp) );
    alpha_n = M_PI - acos( (1-1e-6) * (inner_3d(dnn, d2)) / (l2 * lnn) );
    
    fn3 = 1.0 * (alpha_p - alpha_n);
    V->f_v_tot += fabs(fn3);
     
    //fn3 *= 0.5;
  } else {

    fn3 = 0.0;
  
  }

  /* Generate a unit vector along the field line and scale it to the
   * calculated force.
   */

  sum_3d(force,d1n, d2n);
  scale_3d(force, force, (fn1 + fn3) / norm_3d(force) / (l1+l2) / 2);

  /* Stick the force where it belongs in the VERTEX's force vector.*/

  if (!finite(norm_3d(force))) {
    // printf("somethign is wrong with the forces on vertex %d, fn1=, fn3=%g, force=(%g,%g,%g) \n", V->label, fn1,fn3,force[0],force[1],force[2]);
    //printf("somethign is wrong with the forces on vertex %d, l1=%g,l2=%g,d1n=%g,d2n=%g,fn1=%g,fn3=%g,force=(%g,%g,%g) \n",V->label, l1,l2,d1nr,d2nr, fn1,fn3,force[0],force[1],force[2]);
    //printf("somethign is wrong with the forces on vertex %d, lnn=%g,lpp=%g,acos_p=%g, acos_n=%g,alpha_p=%g,alpha_n=%g,fn1=%g,fn3=%g,force=(%g,%g,%g) \n",V->label, lnn,lpp,acos( (inner_3d(dpp, d1)) / (l1 * lpp) ),acos( (inner_3d(dnn, d2)) / (l2 * lnn) ),alpha_p,alpha_n, fn1,fn3,force[0],force[1],force[2]);
    //printf("somethign is wrong with the forces on vertex %d, l1=%g,l2=%g,lnn=%g,lpp=%g,inner_n=%g,inner_p=%ginside_cos_p=%g, inside_cos_n=%g,acos_p=%g, acos_n=%g,alpha_p=%g,alpha_n=%g,fn3=%g \n",V->label, l1,l2,lnn,lpp,inner_3d(dnn, d2),inner_3d(dpp, d1),(inner_3d(dpp, d1)/(l1 * lpp)),(inner_3d(dnn, d2)/(l2 * lnn)), acos( (inner_3d(dpp, d1)) / (l1 * lpp) ),acos( (inner_3d(dnn, d2)) / (l2 * lnn) ),alpha_p,alpha_n,fn3);

    printf("somethign is wrong with the forces on vertex %d on fl%d, ln is %d on %d, nn is %d on %d, l1=%g,l2=%g,lnn=%g,lpp=%g,x=(%g,%g,%g),xn=(%g,%g,%g),xnn=(%g,%g,%g)\n",V->label,V->line->label,V->next->label,V->next->line->label, V->next->next->label,V->next->next->line->label,l1,l2,lnn,lpp,V->x[0],V->x[1],V->x[2],V->next->x[0],V->next->x[1],V->next->x[2],V->next->next->x[0],V->next->next->x[1],V->next->next->x[2]);
      }


  sum_3d(V->f_v, V->f_v, force);
  V->f_v_tot+= norm_3d(force);

  return;
}

/**********************************************************************
 * f_vertex4 
 *
 * Similar to f_vertex3, it has vertex repulstion, curvature
 * attraction, and attraction to close fluxons. It is modified
 * slightly for more stability. (f_vertex3 has vertex flow)
 *
 */
void f_vertex5(VERTEX *V, HULL_VERTEX *verts) {
  NUM force[3];
  NUM d1[3], d1n[3], d2[3], d2n[3];
  NUM l1, l2;
  NUM d1nr,d2nr,fn1,fn2,fn3; /* fn's are scalers */
  NUM dpp[3], dnn[3];
  NUM lnn, lpp, i1; 
  NUM alpha_p, alpha_n;
  NUM fn1_scale=.01,fn2_scale=2., fn3_scale=1.;

  /* Exclude endpoints */
  if(!V->next || !V->prev || V->label < 0)
    return;

  // Find the normal vectors pointing along the previous and next segments,
  // and collect the lengths of the original vectors...

  diff_3d(d1,V->x,V->prev->x);
  scale_3d(d1n, d1, (d1nr = 1.0 / (l1=norm_3d(d1)))); /* assignment */

  diff_3d(d2,V->next->x, V->x);
  scale_3d(d2n, d2, (d2nr = 1.0 / (l2=norm_3d(d2)))); /* assignment */

  /* Repulsive force from nearest neighbors.  
   */
  fn1 = (d1nr*d1nr-d2nr*d2nr) * ( (l1+l2)*0.5 ) * ( (l1+l2)*0.5 );
  fn1 *= fn1_scale;

  /* Proximity-attractive force.  This attracts vertices toward places
   * where field lines are interacting.  
   */

  if(V->prev && V->next) {
    NUM r_clp, r_cln;
  
    r_clp = l1 / V->prev->r_cl;
    r_clp *= r_clp;

    r_cln = l2 / V->r_cl;
    r_cln *= r_cln;
  
    fn2 = (r_cln - r_clp);
    fn2 *= fn2_scale;
      
  } else {

    fn2 = 0.0;

  }

  /* Curvature-attractive force.  This attracts vertices toward
   * curvature so that sharp angles attract vertices to smooth
   * themselves out.
   */

  if(V->prev->prev && V->next->next) {

    diff_3d(dnn,V->next->next->x,V->next->x);
    lnn = norm_3d(dnn);

    diff_3d(dpp,V->prev->x,V->prev->prev->x);
    lpp = norm_3d(dpp);

    alpha_p = M_PI - acos( (1-1e-6) * (inner_3d(dpp, d1)) / (l1 * lpp) );
    alpha_n = M_PI - acos( (1-1e-6) * (inner_3d(dnn, d2)) / (l2 * lnn) );
    
    fn3 = 1. * (alpha_p - alpha_n);
     
    fn3 *= fn3_scale;
  
  } else {
    alpha_p = 0; alpha_n=0;
    fn3 = 0.0;
  
  }

  /* Generate a unit vector along the field line and scale it to the
   * calculated force.  Remember, we want force per unit length (hence the
   * final division by the harmonic mean of l1 and l2).
   */

  sum_3d(force,d1n, d2n);
  scale_3d(force, force, (fn1 + fn2 + fn3) / norm_3d(force) / (l1+l2) / 2);

  if( !finite(norm_3d(force)) ) {
    //printf("vertex %d: vertex force was (%g,%g,%g); f_vertex5 is (%g,%g,%g) (total %g); i1 is %g\n",V->label,V->f_v[0],V->f_v[1],V->f_v[2],force[0],force[1],force[2],fabs((fn2+fn3)/i1),i1);
    // printf("\tfn2 is %g; fn3 is %g; V->prev->r_cl is %g; V->next->r_cl is %g; l1 is %g; l2 is %g; lpp is %g; lnn is %g; alpha_p is %g; alpha_n is %g; alpha_p arg is %g; alpha_n arg is %g\n",fn2, fn3, V->prev->r_cl, V->next->r_cl, l1, l2, lpp, lnn, alpha_p, alpha_n, (inner_3d(dpp,d1))/(l1*lpp),(inner_3d(dnn,d2))/(l2*lnn));
    printf("vertex %d: f_vertex5 is (%g,%g,%g) (total %g);fn1 is %g, fn2 is %g fn3 is %g \n",V->label,force[0],force[1],force[2],norm_3d(force), fn1,fn2,fn3);
	   }
  
    
  sum_3d(V->f_v, V->f_v, force);
  V->f_v_tot += norm_3d(force);


  return;
}
 

/**********************************************************************
 **********************************************************************
 *****
 ***** newer forces are NOT divided by B.
 * 
 * Note that you need a B calculation _in_advance_ of any of the other
 * force or geometrical calculations.  The B routines must calculate
 * not only the magnetic field but also the projected distances and 
 * angles to the relevant other points.
 *
 */

/**********************************************************************
 * b_simple
 * 
 * Simplest possible B calculation: just 1/area.  Fails when the Voronoi 
 * cell is open, as the area is then infinity and the B field is zero.
 * Useful for test cases and other constrained solutions.
 * 
 *
 * (august 2005)
 */
void b_simple (VERTEX *V, HULL_VERTEX *verts) {
  NUM Bmag = 0;
  NUM vec1[3];
  int n = V->neighbors.n;
  int i;
  int bad = 0;

  for(i=0;i<n && !bad;i++) {
    NUM A;
    HULL_VERTEX *left = &(verts[i]);
    HULL_VERTEX *right = (i==0) ? &(verts[n-1]) : &(verts[i-1]);

    if(left->open || right->open){
      bad=1;
    } else {
      Bmag += cross_2d( left->p, right->p );
    }
  }
  
  if(bad) {
    V->b_mag = 0;
    V->b_vec[0] = 0;
    V->b_vec[1] = 0;
    V->b_vec[2] = 0;
    return;
  }

  Bmag = 1.0/Bmag;
  Bmag *= V->line->flux * (1/PI);
  V->b_mag = Bmag;
  diff_3d(vec1, V->next->x, V->x);
  scale_3d(V->b_vec, vec1, Bmag/norm_3d(vec1));
  return;
}
    
/**********************************************************************
 * e_simple
 * 
 * Very simple magnetic-energy calculation per vertex.  Follows the
 * b_simple method -- if the cell is open then the energy associated
 * with that vertex is zero, as the area is then infinity and the B
 * field is zero.
 * 
 */
void e_simple (VERTEX *V, HULL_VERTEX *verts) {
  int flux = 1;
  NUM ds = 0;
  NUM Area = 0;
  NUM Energy = 0;
  int n = V->neighbors.n;
  int i;
  int bad = 0;

  if (!V->next) { // if ds=0, then energy=0
    V->energy = 0;
    return;
  } else {
    ds = cart_3d(V->x, V->next->x); //length to next vertex
  }

  for(i=0;i<n && !bad;i++) {
    NUM A;
    HULL_VERTEX *left = &(verts[i]);
    HULL_VERTEX *right = (i==(n-1)) ? &(verts[0]) : &(verts[i+1]); 
    /* right is the next hull vertex or the first one if 
       the ith is the last */

    if(left->open || right->open){
      bad=1;
    } else {
      A = 0.5*cross_2d(left->p,right->p);
      Area += A;
    }
  }

  if(Area < 0) {
    fprintf(stderr,"Hey!  Area is less than zero in e_simple (vertex %d!)\n",V->label);
    fflush(stderr);
    Area = fabs (Area);
  }

  if(bad) {
    V->energy = 0;
    return;
  }

  Energy = ds * flux * flux / (8 * PI * Area);
  V->energy = Energy;

  return;
}

/**********************************************************************
 * e_simple2
 * 
 * Simple magnetic-energy calculation per vertex. Only calculates
 * magnetic energy. Accounts for open hulls by calculating fluxes
 * of hull segments relative to angle and sets the open segment 
 * energy to zero.
 * 
 */

void e_simple2 (VERTEX *V, HULL_VERTEX *verts) {
  NUM flux = 1;
  NUM iflux = 0;
  NUM ds = 0;
  NUM angle = 0;
  NUM pseudo_energy =0;
  NUM Energy = 0;
  int n = V->neighbors.n;
  int i;

  if (!V->next) { // if ds=0, then energy=0
    V->energy = 0;
    return;
  } else {
    ds = cart_3d(V->x, V->next->x); //length to next vertex
  }

  for(i=0; i<n; i++) { //i from 0 to n-1
    NUM A;
    HULL_VERTEX *left = &(verts[i]);
    HULL_VERTEX *right = (i==(n-1)) ? &(verts[0]) : &(verts[i+1]); 
    /* right is the next hull vertex or the first one if 
       the ith is the last */

    if(!left->open && !right->open) {
      A = 0.5 * cross_2d(left->p,right->p);              // |left->p| * |right->p| * sin(theta) (cross-sectional area)
      angle = -atan2(left->p[1],left->p[0]) + atan2(right->p[1],right->p[0]); // theta
      /*fprintf(stderr,"angle=%g deg, \n",angle*180./PI);*/

      if (angle < -PI){
	/*fprintf(stderr,"angle=%g deg, vertex=%d \n",angle*180/PI, V->label);*/
	angle += 2.*PI;
	}


      iflux = flux * ( angle / (2.*PI) );  // Percentage of flux in this segment
      pseudo_energy += (iflux * iflux / A);  // sum segments

    }
  }

  if(pseudo_energy < 0) {
    fprintf(stderr,"Hey!  energy is less than zero in e_simple2 (vertex %d!)\n",V->label);
    fflush(stderr);
    pseudo_energy = fabs (pseudo_energy);
  }


  Energy= pseudo_energy * (ds / (8. * PI) );
  V->energy = Energy;

  return;
}

/**********************************************************************
 * e_open
 * 
 * Magnetic-energy calculation per vertex. Uses angular equipartition
 * of flux. Accounts for open hulls by calculating fluxes of hull
 * segments relative to angle and integrates the area of the open
 * cells.
 * 
 */

void e_open (VERTEX *V, HULL_VERTEX *verts) {
  NUM flux = 1;
  NUM ds = 0;
  NUM ienergy =0;
  NUM pseudo_energy =0;
  NUM Energy = 0;
  NUM factor=8.*PI;
  int n = V->neighbors.n;
  int i;

  if (!V->next) { // if ds=0, then energy=0
    V->energy = 0;
    return;
  } else {
    ds = cart_3d(V->x, V->next->x); //length to next vertex
  }

  for(i=0; i<n; i++) { //i from 0 to n-1

    HULL_VERTEX *current_h = &(verts[i]);
    HULL_VERTEX *next_h = (i==(n-1)) ? &(verts[0]) : &(verts[i+1]); 
    /* next is the next hull vertex or the first one if 
       the ith is the last */
    
    // n=h+1
    VERTEX *current_n = (i==(n-1)) ? 
                        ((VERTEX *)(V->neighbors.stuff[0])) : 
                        ((VERTEX *)(V->neighbors.stuff[i+1]));
         //to get to current neighbor, do i+1 from current hull 
    
    NUM ro = .5 * current_n->r;  //half distance to corresponding neighbor
    NUM nhull = (next_h->a_l < 0) ? next_h->a_l + (2.*PI) : next_h->a_l;
    NUM chull = (current_h->a_r < 0) ? current_h->a_r + (2.*PI) : current_h->a_r;
    NUM cneigh = (current_n->a < 0) ? current_n->a + (2.*PI) : current_n->a;
    NUM ang_dif = 0;

    ang_dif = nhull - chull;

    if (ang_dif < 0) {
      ang_dif += 2. * PI;
    }

    ienergy =  ang_dif + 
	       .5*(sin(2. * (nhull - cneigh)) - 
                   sin(2. * (chull - cneigh)) );
    pseudo_energy += ienergy/(ro*ro);

    if (ienergy/(ro*ro) < 0) {
      fprintf(stderr,"Hey!  energy is less than zero in e_open at ith hull (vertex %d!)\n",V->label);
      fprintf(stderr,"cur_h=%3.3f,%d;next_h=%3.3f,%d;cur_n=%3.3f;ro=%3.3f,ei=%3.3f,i=%d,vertex=%d \n",
	      chull,   current_h->open,
	      nhull,         next_h->open,
	      cneigh,
	      ro, ienergy/ro,i,V->label);
      fprintf(stderr,"sin part = %3.3f \n",.5*(sin(2. * (nhull - cneigh)) - sin(2. * (chull - cneigh))));
      fprintf(stderr,"first part = %3.3f \n", nhull - chull);
      fflush(stderr);
    }

    /*        if (current_h->open){
	  fprintf(stderr,"this one open,i=%d, p_energy=%f,vertex %d\n",i,pseudo_energy,V->label);
      fflush(stderr);
      } */
  }

  Energy= pseudo_energy * ((ds * flux * flux) / (4. * PI * PI * factor) );
  V->energy = Energy;
  //fprintf(stderr,"n=%d, energy=%f,vertex %d\n",n,Energy,V->label);
  //  fflush(stderr);
  
  return;
}


/**********************************************************************
 * e_eqa
 *
 * Calculates the magnetic energy contained in a fluxel, using angular 
 * equipartition of flux.  See DeForest Notebox X, pp. 122-123
 */
void e_eqa(VERTEX *V, HULL_VERTEX *verts) {
  NUM E_accum = 0;
  static NUM R;
  int n = V->neighbors.n;
  int i;

  for(i=0; i<n; i++) {
    HULL_VERTEX *left = &(verts[i]);
    HULL_VERTEX *right = (i==0) ? &(verts[n-1]) : &(verts[i-1]);
    NUM righta, lefta,deltaphi;

    righta = right->a_r;
    lefta = left->a_l;
    deltaphi = lefta-righta;
    TRIM_ANGLE(deltaphi);

    if(deltaphi<0) {
      printf("ASSERTION FAILED in e_eqa - lefta should be >= righta, is not. vertex %d, hullpoint %d, lefta=%g deg, righta=%g deg\n",V->label,i,lefta*180/3.14159,righta*180/3.14159);
    }
    
    R = ((VERTEX *)(V->neighbors.stuff[i]))->r * 0.5;
    E_accum += (deltaphi) + sin(2 * (deltaphi))/ ( 2 * R * R );
  }

  /* FIXME - flux per fluxon! */
  V->energy = E_accum  /* * V->line->flux */  * (4/PI/PI) * cart_3d(V->next->x, V->x);
}
  

/**********************************************************************
 * b_eqa
 *  
 * Not truly a force, but required for the others.  Calculates 
 * B and |B| and inserts them into the VERTEX.  b_eqa uses the 
 * angular equiparition estimate, and delivers a phi-averaged
 * magnitude.  See DeForest notebooks v. VIII, p. 70.
 *
 * As a side effect, also accumulates the smallest neighbor distance.
 *
 * The flux is taken from the fluxon itself, in preparation for 
 * future non-uniform-flux-per-fluxon simulations although now 
 * (august 2004) it's not here.
 * 
 */
void b_eqa(VERTEX *V, HULL_VERTEX *verts) {
  NUM Bmag = 0;
  static NUM vec1[3],vec2[3];
  VERTEX **nv = (VERTEX **)(V->neighbors.stuff);
  int n = V->neighbors.n;
  int i;

  for(i=0;i<n;i++) {
    HULL_VERTEX *left = &(verts[i]);
    HULL_VERTEX *right = (i==0) ? &(verts[n-1]) : &(verts[i-1]);
    VERTEX *v = nv[i];
    NUM righta, lefta;
    /* open left vertices are OK; open right vertices must be calculated 
     * from the current vertex.
     */
    righta = right->a_r - v->a;
    TRIM_ANGLE(righta);

    lefta = left->a_l - v->a;
    TRIM_ANGLE(lefta);

    if(righta > lefta) {
      printf("ASSERTION FAILED in b_eqa: left rel. angle %7.3g(%c) > right rel. angle %7.3g(%c), vertex %d (neighbor %d, %d of %d)\n",lefta*180/PI,left->open?'o':'c', righta*180/PI, right->open?'o':'c',V->label, v->label, i, n);

      {
	int j;
	for(j=0;j<n;j++) 
	  printf("\t%s corner %2d (%c): a_l=%7.3g,a_r=%7.3g (rel. %7.3g)\n",(j==i)?"==>":"   ",j,(verts[j].open?'o':'c'),verts[j].a_l*180/PI,verts[j].a_r*180/PI,verts[j].a_l*180/PI-v->a*180/PI);
      }
    }
    
    if(V->line->fc0->world->verbosity >= 4) {printf("b_eqa force: VERTEX %5d, neighbor %5d, a=%7.3g, r=%7.3g, righta = %7.3g(%c), lefta=%7.3g(%c)\n",V->label, v->label,v->a*180/PI, v->r, righta*180/PI,right->open?'o':'c',lefta*180/PI,left->open?'o':'c');
    }

    /* formula wants \frac{1}{2r_p^2}; but VERTEX r is 2r, so
     * we instead use \frac{2}{r^2}...
     */
    Bmag +=  
      ( 0.5 * ( sin( 2 * lefta ) - sin( 2 * righta ) )
	  + lefta - righta
       ) 
      / 
      (nv[i]->r * nv[i]->r) ;

  }
  
  Bmag *= V->line->flux * PI * PI;

  V->b_mag = Bmag;

  /* The B field along the segment points in the same direction as the segment itself. */
  diff_3d(vec1,V->next->x,V->x);
  scale_3d(V->b_vec,vec1,Bmag/norm_3d(vec1));

  if(V->line->fc0->world->verbosity >= 3) { fprintf(stderr,"b_eqa: VERTEX %4d (fluxon %4d); magnetic field mag. is %7.3g (vec is %7.3g, %7.3g, %7.3g)\n", V->label, V->line->label, V->b_mag, V->b_vec[0],V->b_vec[1],V->b_vec[2]); }


}


/**********************************************************************
 * f_curv_hm
 * Curvature `force':   (B^ . \del) B^ 
 *      = b^_x ( db^_x / dx ) + b^_y ( db^_y / dy ) + b^_z ( db^_z / dz )
 * 
 * but b^ is just l^.
 *
 * The curvature force is a vertex-centered force.  There's actually 
 * a little too much information:  not only do we get an angle and a 
 * length, we get an extra length too (there are *two* legs out of the
 * segment, not just one, and they have differing lengths).  
 *
 * f_curv_hm uses the harmonic mean to find the curvature at the vertex,
 * which lets the smaller radius dominate the curvature.
 *
 * f_curv is multiplied by the magnetic field magnitude, unlike
 * its deprecated predecessor, f_curvature.
 * 
 */
void f_curv_hm(VERTEX *V, HULL_VERTEX *verts) {
  NUM recip_l1, recip_l2, recip_len, len;
  NUM b1hat[3];
  NUM b2hat[3];
  NUM curve[3];
  NUM force[3];


  if(!V->next || !V->prev) 
    return;

  diff_3d(b2hat,V->next->x,V->x);
  scale_3d(b2hat,b2hat, (recip_l2 = 1.0 / norm_3d(b2hat)));

  diff_3d(b1hat,V->x,V->prev->x);
  scale_3d(b1hat,b1hat, (recip_l1 = 1.0 /norm_3d(b1hat)));

  recip_len =  ( recip_l1 + recip_l2 ) / 2;
  
  diff_3d(curve,b2hat,b1hat);
  scale_3d(curve,curve, recip_len * (2 / (1/V->prev->b_mag + 1/V->b_mag) ) );

  sum_3d(V->f_v,V->f_v,curve);

  V->f_v_tot += norm_3d(curve);

  if(V->a < 0 || V->a > (len = ( (recip_l2 > recip_l1) ? (1.0/recip_l2) : (1.0/recip_l1) ) ) );
    V->a = len;
}

/**********************************************************************
 * f_curv_m
 * Curvature `force':   (B^ . \del) B^ 
 *      = b^_x ( db^_x / dx ) + b^_y ( db^_y / dy ) + b^_z ( db^_z / dz )
 * 
 * but b^ is just l^.
 *
 * The curvature force is a vertex-centered force.  There's actually 
 * a little too much information:  not only do we get an angle and a 
 * length, we get an extra length too (there are *two* legs out of the
 * segment, not just one, and they have differing lengths).  
 *
 * f_curv_m uses the standard mean to find the curvature at the vertex;
 * this treats the curvature as spread throughout the fluxon.
 *
 * f_curv is multiplied by the magnetic field magnitude, unlike
 * its deprecated predecessor, f_curvature.
 * 
 */
void f_curv_m(VERTEX *V, HULL_VERTEX *verts) {
  NUM l1, l2, recip_len, len;
  NUM b1hat[3];
  NUM b2hat[3];
  NUM curve[3];
  NUM force[3];


  if(!V->next || !V->prev) 
    return;

  diff_3d(b2hat,V->next->x,V->x);
  scale_3d(b2hat,b2hat, 1.0 / (l2 = norm_3d(b2hat)));

  diff_3d(b1hat,V->x,V->prev->x);
  scale_3d(b1hat,b1hat, 1.0 / (l1 = norm_3d(b1hat)));

  recip_len =  2 / ( l1 + l2 );
  
  diff_3d(curve,b2hat,b1hat);
  scale_3d(curve, curve, recip_len * V->line->flux);

  //  scale_3d(curve,curve, recip_len * 0.5 * (V->prev->b_mag + V->b_mag));

  sum_3d(V->f_v,V->f_v,curve);

  V->f_v_tot += norm_3d(curve);

  if(V->a < 0 || V->a > ( len = ( (l2  < l1) ? l2 : l1) ) ) /* assignment */
    V->a = len;
}

/**********************************************************************
 *
 * f_p_eqa_radial
 *
 * Magnetic pressure force, acting radially all the way around the 
 * Voronoi cell boundary (acts like springs on a corrugated boundary)
 *
 * See DeForest notebooks, Vol. VIII, pp 75-76.
 *
 */

void f_p_eqa_radial(VERTEX *V, HULL_VERTEX *verts) {
  static NUM scr2d[3],scr3d[3];
  int i,n;
  NUM sina, cosa;
  NUM pmatrix[9];
  VERTEX **nv = (VERTEX **)(V->neighbors.stuff);
  NUM fac = V->line->flux * V->line->flux / (PI*PI*PI);
  
  if(!(V->next))
    return;

  /* Get the perpendicular-plane projection matrix */
  /* (This is wasteful - it is generated in geometry.c as well. */
  /* Is there an easy way to pass it around? */
  projmatrix(pmatrix,V->x,V->next->x);
  
  scr2d[0] = scr2d[1] = scr2d[2] = 0;

  n = V->neighbors.n;
  for(i=0; i< n; i++) {
    NUM fperp;
    NUM fpar;
    NUM factor;
    HULL_VERTEX *left = &(verts[i]);
    HULL_VERTEX *right = (i!=0) ? &(verts[i-1]) : &(verts[n-1]);
    VERTEX *v = nv[i];
    NUM righta, lefta;

    righta = right->a_r - v->a;   
    TRIM_ANGLE(righta);

    lefta = left->a_l - v->a;
    TRIM_ANGLE(lefta);

    
    if(lefta<righta) 
      lefta += 2*PI;

    if(  (lefta - righta) > PI ) 
      printf("ASSERTION FAILED: lefta-righta >0 (lefta = %5.3g, righta=%5.3g)\n",lefta,righta);
    
    /* Assemble perpendicular and parallel force components */
    factor = -fac / (v->r * v->r * v->r);
    
    fperp= factor * (0.25 * (sin(2*lefta) - sin(2*righta)) 
		     + 0.5 * (lefta-righta)
		     );
    fpar = factor * ( cos(2*righta) - cos(2*lefta) );

    {
      NUM s, c;
      s = sin(v->a);
      c = cos(v->a);
      scr2d[0] += c * fperp - s * fpar;
      scr2d[1] += s * fperp +  c * fpar;
    }

    if(V->line->fc0->world->verbosity >= 4) {printf("f_p_eqa_radial: VERTEX %4d, neighbor %5d, a=%7.3g, r=%7.3g, righta=%7.3g(%c), lefta=%7.3g(%c), fperp=%7.3g, fpar=%7.3g\n",V->label,v->label,v->a*180/PI,v->r,lefta*180/PI,left->open?'o':'c',righta*180/PI,right->open?'o':'c',fperp,fpar);
    }


  } /* end of neighbor loop */
  
   vec_mmult_3d(scr3d,pmatrix,scr2d); /* Convert to 3-D */
   //  mat_vmult_3d(scr3d,pmatrix,scr2d); /* COnvert to 3-D */

  if(V->line->fc0->world->verbosity >= 3) { printf("f_p_eqa_radial: VERTEX %4d (fluxon %4d); total force is %7.3g, %7.3g, %7.3g\n", V->label, V->line->label, scr3d[0],scr3d[1],scr3d[2]); }


  sum_3d(V->f_s,V->f_s,scr3d);

  V->f_s_tot += norm_3d(scr3d);
  
  return;
}
	
/**********************************************************************
 * Reconnection criteria: each accepts a vertex and a list of parameters whose
 * meaning is defined here and in the table at top.  It tests each of the vertex's
 * neighbors against the criterion, and if any of them succeed it immediately returns
 * a pointer to the guilty party, for subsequent reconnection elsewhere.
 */


/******************************
 * rc_a_ad2: threshold angle and angle/d-squared
 * param 0: angle (in radians)
 * param 1: a/d-squared (
 */
VERTEX *rc_a_ad2(VERTEX *v, NUM *params) {
  NUM ath  = params[0];
  NUM ad2th = params[1];
  int i;
  VERTEX *v2;
  NUM d;
  NUM a;
  NUM p1[3], p2[3];
  NUM pm[9];
  NUM l1l2;
  int verbosity = v->line->fc0->world->verbosity;
  
  if(verbosity) {
    printf("aad2: ");
  }

  /* Don't even try unless we've got an element to look at */
  if( !v || !(v->next) )
    return 0;

  if(verbosity > 1) {
    printf("Vertex %d: neighbors are: %d",v->label, ((VERTEX *)(v->neighbors.stuff[0]))->label);
    for(i=1;i<v->neighbors.n; i++) 
      printf(", %d",((VERTEX *)(v->neighbors.stuff[i]))->label);
    printf("\n");
  }

  for(i=0;i<v->neighbors.n;i++) {

    v2 = (VERTEX *)(v->neighbors.stuff[i]);

    if(!V_ISDUMMY(v2) && v2->next) {

      /* p1 and p2 get the closest approach points */
      d = fl_segment_deluxe_dist(p1, p2, v, v2);
	       
      /* Find matrix projecting them into the z axis */
      projmatrix(pm, p1, p2);
      
      /* Project the segments into that plane */
      diff_3d(p1, v->next->x, v->x);
      mat_vmult_3d(v->scr, pm, p1);

      diff_3d(p2, v2->next->x, v2->x);
      mat_vmult_3d( v2->scr, pm, p2);


      /* Now v->scr and v2->scr contain the projected 2-vecs of the corresponding
       * fluxel directions.  Now calculate the angle.  (Could avoid square roots
       * entirely using the half-angle formula, but no sense fixing that since this
       * isn't the hot spot for the whole relaxation).
       */
      l1l2 = sqrt( norm2_2d( v->scr ) * norm2_2d( v2->scr )  );
      a = acos( inner_2d( v->scr , v2->scr )  /  l1l2 );

      if(verbosity>1) {
	printf( "vertices %d-%d: norm2(v->scr)=%.2g, norm2(v2->scr)=%.2g, l1l2=%.2g, inner=%.2g, a=%.2g, d=%.2g, ad2=%.2g, ath=%.2g, ad2th=%.2g\n",
		v->label,
		v2->label,
		norm2_2d(v->scr),
		norm2_2d(v2->scr),
		l1l2,
		inner_2d(v->scr,v2->scr),
		a,
		d,
		a/d/d,
		ath,
		ad2th
		);
      }
      if(a     >       ath      &&
	 a / d / d >   ad2th)
	return v2;
    }
  }
  return 0;
}

VERTEX *rc_a_ad2_h(VERTEX *v, NUM *params) {
  NUM ath  = params[0];
  NUM ad2th = params[1];
  NUM h = params[2];
  int i;
  VERTEX *v2;
  NUM d;
  NUM a;
  NUM p1[3], p2[3];
  NUM pm[9];
  NUM l1l2;
  int verbosity = v->line->fc0->world->verbosity;
  
  if(verbosity) {
    printf("aad2: ");
  }

  /* Don't even try unless we've got an element to look at */
  if( !v || !(v->next) )
    return 0;

  if(verbosity > 1) {
    printf("Vertex %d: neighbors are: %d",v->label, ((VERTEX *)(v->neighbors.stuff[0]))->label);
    for(i=1;i<v->neighbors.n; i++) 
      printf(", %d",((VERTEX *)(v->neighbors.stuff[i]))->label);
    printf("\n");
  }

  for(i=0;i<v->neighbors.n;i++) {

    v2 = (VERTEX *)(v->neighbors.stuff[i]);

    if(!V_ISDUMMY(v2) && v2->next) {

      /* p1 and p2 get the closest approach points */
      d = fl_segment_deluxe_dist(p1, p2, v, v2);
	       
      /* Find matrix projecting them into the z axis */
      projmatrix(pm, p1, p2);
      
      /* Project the segments into that plane */
      diff_3d(p1, v->next->x, v->x);
      mat_vmult_3d(v->scr, pm, p1);

      diff_3d(p2, v2->next->x, v2->x);
      mat_vmult_3d( v2->scr, pm, p2);


      /* Now v->scr and v2->scr contain the projected 2-vecs of the corresponding
       * fluxel directions.  Now calculate the angle.  (Could avoid square roots
       * entirely using the half-angle formula, but no sense fixing that since this
       * isn't the hot spot for the whole relaxation).
       */
      l1l2 = sqrt( norm2_2d( v->scr ) * norm2_2d( v2->scr )  );
      a = acos( inner_2d( v->scr , v2->scr )  /  l1l2 );

      if(verbosity>1) {
	printf( "vertices %d-%d: norm2(v->scr)=%.2g, norm2(v2->scr)=%.2g, l1l2=%.2g, inner=%.2g, a=%.2g, d=%.2g, ad2=%.2g, ath=%.2g, ad2th=%.2g\n",
		v->label,
		v2->label,
		norm2_2d(v->scr),
		norm2_2d(v2->scr),
		l1l2,
		inner_2d(v->scr,v2->scr),
		a,
		d,
		a/d/d,
		ath,
		ad2th
		);
      }
      if(a     >       ath      &&
	 exp(v->x[2]/h) * a / d / d >   ad2th)
	return v2;
    }
  }
  return 0;
}
