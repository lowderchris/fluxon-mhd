/*
 * physics.c
 * 
 * Contains force laws and other directly physics-related 
 * calculations.
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

#include <math.h>
#include "physics.h"
#include "geometry.h"

/* global force name table */
struct FLUX_FORCES FLUX_FORCES[] = {
  {"f_curvature","simple curvature force over B",f_curvature},
  {"f_pressure_equi","Kankelborg/DeForest 2001 pressure law over B",f_pressure_equi},
  {"f_vertex","Vertex distribution pseudo-force",f_vertex},
  {0,0,0}
};


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
 * The forces used for field line manipulation have |B| divided out of
 * them, to simplify their calculation.  
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

  recip_len =  ( recip_l1 + recip_l2 ) / 2;
  
  diff_3d(curve,b2hat,b1hat);
  scale_3d(curve,curve, recip_len);
  sum_3d(V->f_v,V->f_v,curve);
  V->f_v_tot += norm_3d(curve);

  if(V->a < 0 || V->a > (len = ( (recip_l2 > recip_l1) ? (1.0/recip_l2) : (1.0/recip_l1) ) ) );
    V->a = len;

    /*  printf("f_curvature: V%4d, force=(%g,%g,%g)\n",V->label,force[0],force[1],force[2]);*/

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
    HULL_VERTEX *right, *left;
    VERTEX *N = (VERTEX *)V->neighbors.stuff[i];

    right = &verts[i];
    left = (i==0) ? &verts[V->neighbors.n - 1] : &verts[i-1];
    
    /* Calculate delta-phi by calculating the maximal phis and subtracting. 
     * In the open case, the endpoint phi is the angle at which the 
     * perpendicular bisector to r_i goes away.  In the closed case, 
     * it's the angle to the vertex.
     */
    if(right->open){
      phi_r = atan2( N->scr[1], N->scr[0] ) + M_PI_2 ;
    } else 
      phi_r = atan2( right->p[1], right->p[0] );
    if(phi_r < 0)
      phi_r += 2 * M_PI;
    
    if(left->open) {
      phi_l = atan2( N->scr[1], N->scr[0] ) - M_PI_2 ;
    } else
      phi_l = atan2( left->p[1], left->p[0] );
    if(phi_l < 0)
      phi_l += 2 * M_PI;
    
    deltaphi = (phi_r - phi_l);
    if(deltaphi < -1e-4 ) /* Should be 0; allow for slop */
      deltaphi += 2 * M_PI;

    /*    printf("from %d - %d: deltaphi is %g deg; left->open=%d; right->open=%d\n",((VERTEX *)(V->neighbors.stuff[i]))->label, ((VERTEX *)(V->neighbors.stuff[(i+1)%V->neighbors.n]))->label,deltaphi*180/M_PI,left->open,right->open); */

    if(deltaphi > M_PI) 
      fprintf(stderr,"Assertion failed!  deltaphi > M_PI in f_pressure_equi (%18.12g deg)\n",deltaphi*180/M_PI);
    
    /*    fprintf(stderr," VERTEX #%d, i=%d: deltaphi=%g\tr=%g\n",N->label,i,deltaphi,N->r);*/


    f = deltaphi / M_PI / N->r;

    N->scr[2] = 0; /* Make sure we're operating in the perp. plane */
    
    vec_mmult_3d(f_i,pmatrix,N->scr);       /* Put back into 3-D */
    /*    printf("  in plane: (%g,%g,%g); projects back to (%g,%g,%g)\n",N->scr[0],N->scr[1],N->scr[2],f_i[0],f_i[1],f_i[2]); */

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
 * fluxon_vertex_force
 * Pseudoforce that moves vertices along field lines.  Used to keep them 
 * nicely evened out along the lines and to attract them to 
 * curvature.
 * 
 * Strictly speaking, this is two forces:  a vertex repulsive force
 * and another force that attracts vertices toward curvature.  But
 * the two should always be used together, so they're included together.
 * 
 */
void f_vertex(VERTEX *V, HULL_VERTEX *verts) {
  NUM force[3];
  NUM d1[3], d2[3];
  NUM d1nr,d2nr,fn,d1n,d2n;

  /* Exclude endpoints */
  if(!V->next || !V->prev)
    return;


  /* Repulsive force from nearest neighbors.  The force drops as 
   *  1/r.
   */
  diff_3d(d1,V->x,V->prev->x);
  scale_3d(d1, d1, (d1nr = 1.0/(d1n = norm_3d(d1)))); /* assignment */

  diff_3d(d2,V->next->x, V->x);
  scale_3d(d2, d2, (d2nr = 1.0/(d2n = norm_3d(d2)))); /* assignment */
  
  fn = (d1nr - d2nr);
  V->f_v_tot += fabs(fn);

  /* Curvature-attractive force.  This keeps vertices from wandering
   * away from zones of more curvature.  It's meant to balance the
   * repulsive force, allowing more concentration of vertices where
   * they're needed.  The force is the secant of the vertex angle on each side.
   */
  
  /* Exclude endpoint-neighbors  (endpoints already excluded) */
  if(0 && V->next->next  && V->prev->prev) {
    NUM d0[3], d3[3];
    NUM d0n, d3n, sec1, sec2;
    
    diff_3d(d0,V->prev->x,V->prev->prev->x);
    d0n = norm_3d(d0);
    
    diff_3d(d3,V->next->next->x,V->next->x);
    d3n = norm_3d(d3);
    
    sec1 = d0n / inner_3d(d0,d1);
    sec2 = d3n / inner_3d(d2,d3);

    V->f_v_tot += fabs(sec2 - sec1);
    fn +=  sec2 - sec1;
  }


  /* Proximity-attractive force.  This attracts vertices toward places
   * where field lines are interacting. 
   */
  if(V->prev && V->next) {
    NUM r_clp, r_cln;

    r_clp = 1.0 / V->prev->r_cl;
    r_cln = 1.0 / V->r_cl;

    V->f_v_tot += fabs(0.5 * (r_cln - r_clp));
    fn += 0.5 * (r_cln - r_clp);
  }
    

  /* Generate a vector along the field line and scale it to the 
   * calculated force.  
   *
   * Finally, stick the force where it belongs in the VERTEX's force vector.
   */

  sum_3d(force,d1, d2);
  scale_3d(force, force, 0.1 * fn / norm_3d(force));

  /*printf("f_vertex: V%4d, force=(%g,%g,%g)\n",V->label,force[0],force[1],force[2]);*/

  sum_3d(V->f_v, V->f_v, force);

  return;
}
 
 
