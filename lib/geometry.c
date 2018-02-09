/**********************************************************************
 * Geometry.c -- routines that embody geometrical operations
 * used for the fieldine model.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Southwest Research Institute, 2004-2008
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
 * This file is part of the FLUX 2.2 release (22-Nov-2008)
 *
 */
#include <math.h>
#include <stdio.h>
#include "data.h"
#include "geometry.h"

// model.h is needed only for find_simplex, which pulls expand_lengthwise and expand_via_neighbors.
// Should find_simplex be moved to model.c?
#include "model.h"

#ifndef NAN
#define NAN (nan("NAN"))
#endif

char *code_info_geometry="%%%FILE%%% (new unsorted hull routine)";

/**********************************************************************
 **********************************************************************
 *****  Vector basics
 *****  norm, inner product, cross product, scalar product, sum,
 *****  difference, copy.
 */

/**********************************************************************
 * norm - find the length of a vector.  2-d or 3-d.
 */
NUM norm2_2d(NUM *x) {
  NUM out;
  out = *x * *x; x++;;
  out += *x * *x;
  return out;
}

NUM norm_2d(NUM *x) {
  NUM out;
  out = *x * *x; x++;
  out += *x * *x;
  return sqrt(out);
}

NUM norm2_3d(NUM *x) {
  NUM out;
  out = *x * *x; x++;
  out += *x * *x; x++;
  out += *x * *x; x++;
  return out;
}

NUM norm_3d(NUM *x) {
  NUM out;
  out = *x * *x; x++;
  out += *x * *x; x++;
  out += *x * *x;
  return sqrt(out);
}


/**********************************************************************
 * inner - Find the inner product of two vectors.  2-d or 3-d
 */
NUM inner_2d(NUM *p0, NUM *p1) {
  NUM out = *(p0++) * *(p1++);
  return out + ( *p0 * *p1 );
}

NUM inner_3d(NUM *p0, NUM *p1) {
  NUM out;
  out = *(p0++) * *(p1++);
  out += *(p0++) * *(p1++);
  out += *(p0) * *(p1);
  return out;
}

/**********************************************************************
 * cross - Cross product of two vectors.  In 2-D, this is a
 * scalar, in 3-D it's a vector and hence doesn't return the value
 * directly.
 */

NUM cross_2d(NUM *p0, NUM *p1) {
  return (p1[1]*p0[0] - p1[0]*p0[1]);
}

/* cross_3d is defined to cross in geometry.h */
 void cross(NUM *out, NUM *p0, NUM *p1) {
  *(out++) = p0[1]*p1[2]-p0[2]*p1[1];
  *(out++) = p0[2]*p1[0]-p0[0]*p1[2];
  *out     = p0[0]*p1[1]-p0[1]*p1[0];
}

/**********************************************************************
 * scale - Multiply a vector by a scalar in-place, and stick the
 * result in the given destination. OK to have the destination be
 * the source vector.
 */
void scale_3d(NUM *out, NUM *a, NUM b) {
  *(out++) = *(a++) * b;
  *(out++) = *(a++) * b;
  *(out) = *(a) * b;
}

void scale_2d(NUM *out, NUM *a, NUM b) {
  *(out++) = *(a++) * b;
  *(out) = *(a++) * b;
}


/**********************************************************************
 * sum - Add two 3-vectors and stick 'em into the destination array.
 * OK to have the destination be one of the sources.
 */
void diff_2d(NUM *out, NUM *a, NUM *b){
  *(out++) = *(a++) + *(b++);
  *(out) = *(a) + *(b);
}
void sum_3d(NUM *out, NUM *a, NUM *b) {
  *(out++) = *(a++) + *(b++);
  *(out++) = *(a++) + *(b++);
  *(out) = *(a) + *(b);
}

/**********************************************************************
 * diff -- Subtract two 3 vectors (a - b) and put the result in the
 * destination array. OK to have the destination be one of the sources.
 */
void diff_3d(NUM *out, NUM *a, NUM *b) {
  *(out++) = *(a++) - *(b++);
  *(out++) = *(a++) - *(b++);
  *(out) = *(a) - *(b);
}
#define diff_3d_macro(out,a,b) ( ( (out)[0] = (a)[0]-(b)[0] ), ( (out)[1] = (a)[1]-(b)[1]), ((out)[2]=(a)[2]-(b)[2]))


/**********************************************************************
 * centroid - Find the centroid of the plane formed from 3 3-vectors
 */
void centroid(NUM *out, NUM *a, NUM *b, NUM *c) {
  *(out++) = *(a++) + *(b++) + *(c++);
  *(out++) = *(a++) + *(b++) + *(c++);
  *(out) = (*(a) + *(b) + *(c))/3;
}


/**********************************************************************
 * mpdist - Find the distance from the point x to the midpoint of the
 * vector from a-b
 */
NUM mpdist(NUM *a, NUM*b, NUM *x) {
  NUM out;
  out = (0.5 * (*a - *b)) - *x; a++; b++; x++;
  out += (0.5 * (*a - *b)) - *x; a++; b++; x++;
  out += (0.5 * (*a - *b)) - *x;
  return out;
}

/**********************************************************************
 * cp_3d - Copy a 3-vector
 */
void cp_3d(NUM *a, NUM *b) {
  *(a++) = *(b++);
  *(a++) = *(b++);
  *(a) = *(b);
}



/**********************************************************************
 **********************************************************************
 ***** Matrix handling routines
 */

/**********************************************************************
 * rotmat_2d - Return a 2x2 matrix containing the specified rotation.
 * The matrix comes as an array of 4 elements:  top row, then bottom row.
 * It gets put in the specified place.
 */
void rotmat_2d(NUM *out, NUM alpha) {
  out[0] = out[3]  = cos(alpha);
  out[2] = - (out[1] = sin(alpha));
}

void rotmat_2d_fr_slope(NUM *out,NUM dy, NUM dx) {
  NUM r;
  r = sqrt(dy*dy+dx*dx);
  out[0] = out[3] = dx/r;
  out[2] = - (out[1] = dy/r);
}

/**********************************************************************
 * mat_mult_2d - Multiply two 2x2 matrices, returning a 2x2 matrix
 * in the destination.  Buffered, so the destination may be one
 * of the sources.  Each element has its own variable to facilitate
 * registerification.
 */
void mat_mult_2d(NUM *out, NUM *a, NUM *b) {
  NUM lt,rt,lb,rb;
  lt = a[0]*b[0] + a[1]*b[2];
  rt = a[0]*b[1] + a[1]*b[3];
  lb = a[2]*b[0] + a[3]*b[2];
  rb = a[2]*b[1] + a[3]*b[3];
  *out++=lt;
  *out++=rt;
  *out++=lb;
  *out=rb;
}

/**********************************************************************
 * Transpose a 2-D matrix.  This inverts rotations.  It is also
 * trivial. Yes, I know the ^= trick, but it's too much hassle with
 * the switchable NUMs.
 */
void transpose_2x2(NUM *mat) {
  NUM a = mat[1];
  mat[1] = mat[2];
  mat[2] = a;
}

/**********************************************************************
 * transpose_3x3
 * Transpose a 3x3 matrix in situ
 */
void transpose_3x3(NUM *mat) {
  NUM a;
  a=mat[1]; mat[1] = mat[3]; mat[3] = a;
  a=mat[2]; mat[2] = mat[6]; mat[6] = a;
  a=mat[5]; mat[5] = mat[7]; mat[7] = a;
}

/**********************************************************************
 * mat_mult_3x3
 * Multiply two 3x3 matrices, returning a 3x3 matrix in the destination.
 * Not buffered.
 */
void mat_mult_3d(NUM *out, NUM *a, NUM *b) {
  *(out++) = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
  *(out++) = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
  *(out++) = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];

  *(out++) = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
  *(out++) = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
  *(out++) = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];

  *(out++) = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
  *(out++) = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
  *(out)   = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}

/**********************************************************************
 * mat_vmult_2d - Multply a 2x2 matrix by a 2-vector on the right.
 * Buffered, so the destination may be the vector source.
 *
 * The arguments are treated as pointers not POINTs so that they can
 * be written to!
 */
void mat_vmult_2d(NUM *out, NUM *mat, NUM *v) {
  NUM t,b;
  t = mat[0]*v[0]+mat[1]*v[1];
  b = mat[2]*v[0]+mat[3]*v[1];
  *(out++)=t;
  *(out)=b;
}

/**********************************************************************
 * mat_vmult_3d - Multiply a 3x3 matrix by a 3-vector on the right.
 */
void mat_vmult_3d(NUM *out, NUM *mat, NUM *v) {
  *(out++) = v[0]*mat[0] + v[1]*mat[1] + v[2]*mat[2];
  *(out++) = v[0]*mat[3] + v[1]*mat[4] + v[2]*mat[5];
  *(out  ) = v[0]*mat[6] + v[1]*mat[7] + v[2]*mat[8];
}

#define mat_vmult_3d_macro(out,mat,v) (((out)[0] = (mat)[0]*(v)[0] + (mat)[1]*(v)[1] + (mat)[2]*(v)[2]), ((out)[1] = (mat)[3]*(v)[0] + (mat)[4]*(v)[1] + (mat)[5]*(v)[2]), ((out)[2] = (mat)[6]*(v)[0] + (mat)[7]*(v)[1] + (mat)[8]*(v)[2]))


/***********************************************************************
 * vec_mmult_3d - Multiply a 3x3 matrix by a 3-vector on the left.
 * This is the same as transposing the matrix, which is the same as
 * inverting it if it's a rotation matrix!
 */
void vec_mmult_3d(NUM *out, NUM *mat, NUM *v) {
  *(out++) = v[0] * mat[0] + v[1] * mat[3] + v[2] * mat[6];
  *(out++) = v[0] * mat[1] + v[1] * mat[4] + v[2] * mat[7];
  *(out)   = v[0] * mat[2] + v[1] * mat[5] + v[2] * mat[8];
}

/**********************************************************************
 * det_2d - Deteminant of a 2x2 matrix
 */
NUM det_2d(NUM *mat) {
  return mat[0]*mat[3] - mat[1]*mat[2];
}

/**********************************************************************
 * det_3d - Determinant of a 3x3 matrix
 */
NUM det_3d(NUM *mat) {
  return(  mat[0]*mat[4]*mat[8] + mat[1]*mat[5]*mat[6] + mat[2]*mat[3]*mat[7]
          -mat[0]*mat[5]*mat[7] - mat[1]*mat[3]*mat[8] - mat[2]*mat[4]*mat[6]);
}

/**********************************************************************
 **********************************************************************
 *** Plane & topology support
 */

/**********************************************************************
 * points2plane - given three noncolinear points, return the plane that
 * contains them as a PLANE structure.
 */
void points2plane(PLANE *plane, NUM *p0, NUM *p1, NUM *p2) {
  NUM x1[3];
  NUM x2[3];
  NUM cr[3];

  diff_3d(x1, p1, p0);
  diff_3d(x2, p2, p0);
  cross_3d(cr, x1, x2);
  scale_3d(plane->normal, cr, 1.0/norm_3d(cr));
  cp_3d(plane->origin, p0);
}

/**********************************************************************
 * lp2plane - given a line segment P0-P1 and a third point P2,
 * return the plane that passes through P0 and P1, with normal
 * through P2.
 */
void lp2plane(PLANE *plane, NUM *p0, NUM *p1, NUM *p2) {
  NUM x1[3];
  NUM x2[3];
  NUM alpha;

  diff_3d(x1,p1,p0);
  diff_3d(x2,p2,p0);

  // offset plane origin to the point under the normal
  alpha = inner_3d(x1,x2) / norm2_3d(x1);
  scale_3d(plane->origin, x1, alpha );
  sum_3d(plane->origin, plane->origin, p0);

  // copy the normal and, er, normalize it.
  diff_3d(plane->normal, p2, plane->origin);
  scale_3d(plane->normal, plane->normal, 1.0/norm_3d(plane->normal));
}


/**********************************************************************
 * p_l_intersection - given a plane and a line (specified by two points)
 * return the intersection between them.
 * Returns a flag indicating whether the intersection is between
 * the two points (1=between, 0=not between)
 */
int p_l_intersection(NUM *out, PLANE *plane, NUM *p0, NUM *p1) {
  NUM x0[3];
  NUM x1[3];
  NUM zeta0, zeta1;
  diff_3d(x0, p0, plane->origin);
  diff_3d(x1, p1, plane->origin);

  zeta0 = inner_3d(x0, plane->normal);
  zeta1 = inner_3d(x1, plane->normal);

  if(zeta1==0) {
    cp_3d(out,p1);
    return 1;
  } else if(zeta0==0){
    cp_3d(out,p0);
    return 1;
  } else {
    scale_3d(out, x1, -(zeta0/zeta1));
    sum_3d(out, out, x0);
    sum_3d(out, out, plane->origin);
    return (zeta0 * zeta1 <= 0);
  }
  // Never get here...
}

/**********************************************************************
 * xy_l_intersection - find the intersection between a line and the xy plane
 */
int xy_l_intersection(NUM *out, NUM *p0, NUM *p1) {
  if(p1[2]==0) {
    cp_3d(out, p1);
    return 1;
  } else if(p0[2]==0) {
    cp_3d(out, p0);
    return 1;
  } else {
    scale_3d(out, p1, -p0[2]/p1[2]);
    sum_3d(out, out, p0);
    return (p0[2]*p1[2] <= 0);
  }
  // Never get here...
}

/**********************************************************************
 * p_inside_tri - given a triangle and a point in the plane, determine
 * whether the point is strictly inside the triangle (1) or outside it (0).
 */
int p_inside_tri(NUM *tri0, NUM *tri1, NUM *tri2, NUM *p) {
  NUM c1, c2, c3;
  NUM x0[2],x1[2],x2[2];
  diff_2d(x0, tri0, p);
  diff_2d(x1, tri1, p);
  diff_2d(x2, tri2, p);

  c1 = cross_2d(x0,x1);
  c2 = cross_2d(x1,x2);
  c3 = cross_2d(x2,x0);

  return ( (c1 > 0  &&  c2 > 0  && c3 > 0) ||
           (c1 < 0  &&  c2 < 0  && c3 < 0)
           );
}


/**********************************************************************
 * trivloop - given a fluxon, return whether it is a trivial loop, that
 * is to say if it can be contained in a neighborhood that is small compared
 * to the distance to the nearest neighbor that is not part of the fluxon.
 *
 * Trivloop tells you nothing about the topology of the loop -- e.g. whether
 * it is the unknot.  Hence it is not so useful for determining whether a
 * plasmoid should vanish in ideal MHD -- but it is pretty good for deciding
 * whether a small U-loop should vanish in open boundary conditions.
 *
 * A fluxon is a trivial loop if all of its internal distances vanish
 * compared to the corresponding distances to each neighbor.
 *

 * Requires that the neighbor function have already been executed.
 */
#define trivloop_factor 3

int trivloop(FLUXON *f) {
  NUM max_internal_dist = 0;;
  NUM min_neighbor_dist = 1e100;
  VERTEX *v;
  VERTEX *v1;

  for(v=f->start; v && v->next; v=v->next) {
    int i;

    for(i=0; i<v->neighbors.n; i++) {
      NUM r;
      VERTEX *vn = (VERTEX *)(v->neighbors.stuff[i]);
      if(vn->next && vn->line != f) {
        r=p_ls_dist(v->x, vn->x, vn->next->x);
        if(r<min_neighbor_dist) {
          min_neighbor_dist = r;
          if(min_neighbor_dist < max_internal_dist * trivloop_factor) {
            //printf("trivloop: not trivial (%g < %g); v=%d,vn=%d\n",min_neighbor_dist,max_internal_dist * trivloop_factor,v->label,vn->label);
            return 0;
          }
        }
      }
    }

    for(i=0;i<v->nearby.n; i++) {
      NUM r;
      VERTEX *vn = (VERTEX *)(v->nearby.stuff[i]);
      if(vn->next && vn->line != f) {
        r=p_ls_dist(v->x, vn->x, vn->next->x);
        if(r<min_neighbor_dist) {
          min_neighbor_dist = r;
          if(min_neighbor_dist < max_internal_dist * trivloop_factor) {
            //printf("trivloop: not trivial; case 2 (%g < %g); v=%d,vn=%d\n",min_neighbor_dist,max_internal_dist * trivloop_factor,v->label,vn->label);
            return 0;
          }
        }
      }
    }


    for(v1=v; v1&&v1->next; v1=v1->next) {
      NUM r = cart_3d(v->x, v1->x);
      if(r > max_internal_dist) {
        max_internal_dist = r;
        if( min_neighbor_dist < max_internal_dist * trivloop_factor) {
          //printf("trivloop: not trivial; case 3 (%g < %g); v=%d\n",min_neighbor_dist,max_internal_dist * trivloop_factor,v->label);
          return 0;
        }
      }
    }
  }

  return 1;
}

/**********************************************************************
 **********************************************************************
 *** Distances.   Between points, lines, and linesegments...
 ***
 */

/***********************************************************************
 * cart2 - return the squared distance between two points, 3d
 * (avoids the square root if you only need a monotonically increasing
 * function!)
 */
NUM cart2_2d(NUM *x1, NUM *x2) {
  NUM out;
  NUM a;
  a = *(x2++) - *(x1++);
  out = a*a;
  a = *(x2++) - *(x1++);
  return (out + a*a);
}

NUM cart2_3d(NUM *x1, NUM *x2) {
  NUM out;
  NUM a;
  a = *(x2++) - *(x1++);
  out = a*a;
  a = *(x2++) - *(x1++);
  out += a*a;
  a = *(x2++) - *(x1++);
  return (out + a*a);
}


/***********************************************************************
 * cart - Return the distance between two points, 'N' dimensions
 */
NUM cart_2d(NUM *x1, NUM *x2) {
  NUM out;
  NUM a;
  a = *(x2++) - *(x1++);
  out = a*a;
  a = *(x2++) - *(x1++);
  return sqrt(out + a*a);
}

NUM cart_3d(NUM *x1, NUM *x2) {
  NUM out;
  NUM a;
  a = *(x2++) - *(x1++);
  out = a*a;
  a = *(x2++) - *(x1++);
  out += a*a;
  a = *(x2++) - *(x1++);
  return sqrt(out + a*a);
}


/**********************************************************************
 * p_l_dist
 * Find the distance between a point (in 3-space) and a line
 * (also in 3-space).   Relatively simple geometric construction.
 */
NUM p_l_dist(NUM *p0, NUM *x0, NUM *x1) {
  NUM p[3],x[3],u[3];

  diff_3d(x,x1,x0);      /* x = x1 - x0 */
  diff_3d(p,p0,x0);      /* p = p0 - x0 */

  cross(u,x,p);           /* u = p x x  */
  return norm_3d(u)/norm_3d(x);
}

/**********************************************************************
 * p_ls_dist
 * Find the distance between a point (in 3-space) and a line segment
 * (also in 3-space).  This is a little harder because of the conditionals
 * around the endpoints.
 */
NUM p_ls_dist(NUM *p0, NUM *x0, NUM *x1) {
  NUM p[3],x[3],u[3];
  NUM l,a;

  diff_3d(x,x1,x0);

  l = norm_3d(x);
  if(!l)                    /* Degenerate case -- x0==x1 */
    return cart_3d(p0,x0);

  diff_3d(p,p0,x0);

  a = inner_3d(p,x)/l;         /* a = ( (p-x0) . (x1-x0) ) / |x1 - x0| */

  if(a<0)
    return norm_3d(p);     /* Already did the subtraction... */
  else if(a>l)
    return cart_3d(p0,x1);

  cross(u,x,p);
  return norm_3d(u)/l;
}

/**********************************************************************
 * l_l_dist
 * Find the closest approach between two skew lines (in 3-space).
 * You supply two points (A and B) on one line segment and
 * two points (C and D) on the other.   You get back the closest
 * distance between the two LINES, not the segments.
 *
 * This could probably be cleaned up via some vector identities,
 * but at first blush they didn't look to actually improve the computational
 * load much.
 *
 * Formula:  dist = | C' x (D'-C') | / |(D'-C')|
 *
 * (where C' = (C - A) x (B - A)^
 *        D' = (C - A) x (B - A)^
 * )
 */
NUM l_l_dist(NUM a0[3], NUM b0[3], NUM c0[3], NUM d0[3]) {
  NUM b[3],foo[3],c[3],d[3];
  NUM cdpl;

  diff_3d(b,b0,a0);             /* b gets (B - A)-hat */
  scale_3d(b,b,1.0/norm_3d(b));

  diff_3d(foo,c0,a0);
  cross(c,foo,b);      /* c gets (C - A) x (B - A) */

  diff_3d(foo,d0,a0);
  cross(d,foo,b);      /* d gets (D - A) x (B - A) */

  diff_3d(d,d,c);      /* d gets (D' - C') */

  cross(foo,c,d);      /* foo gets   C' x (D' - C') */

  cdpl = norm_3d(d);
  if(cdpl)
    return norm_3d(foo)/cdpl;

  /* Degenerate case -- cd is parallel to ab, or c==d!  Makes no  *
   * difference either way... use the point/line formula.         */
  return p_l_dist(c0, a0, b0);
}

/**********************************************************************
 * p_l_closest_approach and p_ls_closest_approach --
 * Given a line segment and point C, return the location of the
 * extended line's closest approach to C, or the line segment's
 * closest approach to C.  Scant degeneracy checking is done, because
 * the algorithm works more-or-less gracefully in nearly-degenerate cases.
 *
 */

void p_ls_closest_approach(NUM p0[3], NUM a0[3], NUM b0[3], NUM c0[3]) {
  NUM c[3];
  NUM r;

  diff_3d(p0, b0, a0);
  if( p0[0]==0 && p0[1]==0 && p0[2]==0 ) {
    cp_3d(p0,a0);
  } else {

    diff_3d(c, c0, a0);
    r = inner_3d(p0, c) / inner_3d(p0,p0); /* r gets (b.c/b.b) */
    if(r<0)                                /* less than 0: before a */
      cp_3d(p0,a0);
    else if(r > 1)                         /* greater than 1: after b */
      cp_3d(p0,b0);
    else {                                 /* between 0-1: on ab line segment */
      scale_3d(p0, p0, r);
      sum_3d(p0, p0, a0);                  /* put p0 back into original coordinates. */
    }
  }
  return;
}


/**********************************************************************
 * ls_closest_approach
 *
 * You feed in the endpoints of two line segments, and two POINTs.
 * The two points get the closest-approach points of the two line
 * segments.
 *
 * There are nine cases to consider: four endpoint-to-endpoint cases,
 * four endpoint-to-line cases, and one line-to-line case.
 *
 * Method:
 *   - handle degenerate cases.
 *
 *   - Consider the line-to-line case.  If it's valid, use it and return.
 *
 *   - Consider the endpoint-to-line cases. If one is valid, use it and return.
 *
 *   - Consider the endpoint-to-endpoint cases.
 *
 */


void ls_closest_approach(NUM p0[3], NUM p1[3], NUM a0[3], NUM b0[3], NUM c0[3], NUM d0[3]) {
  NUM q1[3],q2[3], b[3], c[3], d[3], scr[3];
  NUM mat[9];
  NUM alpha,beta,d2;

  /******************************
   * Check for degenerate cases: two points, or point-line
   */
  {
    char ab0,cd0;
    ab0=(a0[0]==b0[0] && a0[1]==b0[1] && a0[2]==b0[2]);
    cd0=(c0[0]==d0[0] && c0[1]==d0[1] && c0[2]==d0[2]);

    if(ab0) {
      cp_3d(p0,a0);
      if(cd0) {
        cp_3d(p1,c0);
        /* printf("ls_closest: a==b, c==d (%g,%g,%g)\n",a0[0],a0[1],a0[2]); */
        return;
      } else {
        p_ls_closest_approach(p1,c0,d0,a0);
        /* printf("ls_closest: a==b (%g,%g,%g); closest on cd is (%g,%g,%g)\n",a0[0],a0[1],a0[2],p1[0],p1[1],p1[2]); */
        return;
      }
    }
    if(cd0){
      cp_3d(p1,c0);
      p_ls_closest_approach(p0,a0,b0,c0);
/*       printf("c==d (%g,%g,%g); closest on ab is (%g,%g,%g)\n",c0[0],c0[1],c0[2],p0[0],p0[1],p0[2]); */

    }
  }

  /******************************
   * Find the closest approach of the *lines* AB and CD.
   */

  /***********
   * Place a at the origin and rotate so that b points in the positive z direction
   */

  projmatrix(mat,b0,a0);

  diff_3d(scr,b0,a0);
  mat_vmult_3d(b, mat, scr);

  diff_3d(scr,c0,a0);
  mat_vmult_3d(c, mat, scr);

  diff_3d(scr,d0,a0);
  mat_vmult_3d(d, mat, scr);

  //  printf("ls_closest_approach (in ab frame):\n\tb is %g,%g,%g\n\tc is %g,%g,%g\n\td is %g,%g,%g\n",b[0],b[1],b[2],c[0],c[1],c[2],d[0],d[1],d[2]);

  /**********
   * Find the closest approach of the line cd to the origin.  Works by
   * displacing c to the origin and projecting a perpendicular to the
   * cd line from a until it lands on cd.  The vector is (cd) (
   * (cd).(ca)/(cd.cd) ) But cd.ca is just -(cd.ac), so we can avoid
   * having to actually displace...
   */
  diff_3d(scr,d,c);
  d2= norm2_2d(scr);
  if(d2==0)
    goto parallel;

  /********************
   * now calculate the alpha along cd of the closest approach point (in the plane).
   */
  alpha = - inner_2d(scr,c) / d2;

  /**********
   * Using alpha, find the altitude in the current space
   * of the closest approach, and use that to calculate beta.
   */
  beta = ( (1.0-alpha) * c[2] + alpha * d[2] ) / b[2];



  if(alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1) {
    /* The line-line case is valid: copy and return. */

    //    printf("ls_closest_approach: line-line case. alpha=%g; beta=%g\n",alpha,beta);

    scale_3d(scr, c0, (1.0-alpha) );
    scale_3d(p1,  d0, alpha );
    sum_3d(p1, p1, scr);

    scale_3d(scr, a0, (1.0-beta) );
    scale_3d(p0,  b0, beta );
    sum_3d(p0, p0, scr);
    return;
  }

  /**********
   * line-line is invalid.  Try the four endpoint-line cases
   * (which includes thhe four endpoint-endpoint cases too).
   * Boneheaded but certain.
   */
  {
    char found_one = 0;
    NUM d2,c2;

    /* Check AB - C */
    p_ls_closest_approach(scr, a0, b0, c0);
    d2 = cart2_3d(c0,scr);
    cp_3d(p0,scr);
    cp_3d(p1,c0);
    //        printf("AB-C: dist is %g\n",d2);

    /* Check AB - D */
    p_ls_closest_approach(scr, a0, b0, d0);
    c2 = cart2_3d(d0,scr);
    //        printf("AB-D: dist is %g\n",c2);
    if( c2<d2 ) {
      d2 = c2;
      cp_3d(p0,scr);
      cp_3d(p1,d0);
    }

    /* Check A - CD */
    p_ls_closest_approach(scr, c0, d0, a0);
    c2 = cart2_3d(a0,scr);
    //        printf("A-CD: dist is %g\n",c2);
    if( c2<d2 ) {
      d2 = c2;
      cp_3d(p1,scr);
      cp_3d(p0,a0);
    }

    /* Check B - CD */
    p_ls_closest_approach(scr, c0, d0, b0);
    c2 = cart2_3d(b0,scr);
    //        printf("B-CD: dist is %g\n",c2);
    if( c2<d2 ) {
      d2=c2;
      cp_3d(p1,scr);
      cp_3d(p0,b0);
    }
  }

  return;


/******************************
 * parallel -- parallel-lines case
 */
 parallel:
  //  printf("ls_dist: parallel case...");

  /* Both lines are parallel along the z axis.  b>0, a=0. */

  if( c[2] < 0 && d[2] < 0 ) {
    cp_3d(p0, a0);
    cp_3d(p1, (c[2] >= d[2]) ? c0 : d0);
    return;
  }

  if( c[2] > b[2]  && d[2] > b[2] ) {
    cp_3d(p0, b0);
    cp_3d(p1, (c[2] <= d[2]) ? c0 : d0);
    return;
  }

  if( c[2] >=0 && c[2] <= b[2] ) {
    cp_3d( p1, c0 );
    alpha = c[2]/b[2];
    scale_3d( scr, a0, (1.0-alpha) );
    scale_3d(  p0, b0, alpha );
    sum_3d( p0, p0, scr );
    return;
  }

  if( d[2] >= 0 && d[2] <= b[2] ) {
    cp_3d( p1, d0 );
    alpha = d[2]/b[2];
    scale_3d( scr, a0, (1.0-alpha) );
    scale_3d( p0, b0, alpha );
    sum_3d( p0, p0, scr );
    return;
  }

  /* got here: cd must straddle ab, so a is a right answer. */
  cp_3d(p0, a0);
  alpha = (0 - c[2]) / (d[2]-c[2]);
  scale_3d( scr, c0, (1.0-alpha) );
  scale_3d( p1,  d0, alpha );
  sum_3d(p1, p1, scr);
  return;

}



/**********************************************************************
 * ls_ls_dist
 * Find the distance of closest approach between two line segments in
 * 3-space.  You supply end endpoints A and B of the first
 * segment and C and D of the second segment.  You get back the
 * closest distance between the segments (which might be a
 * l-l distance, a p-p distance, or a l-p distance).
 *
 * This is especially complex because of the multiple cases.
 * This one finds the location of the closest point on
 * each line segment to the other line, and then gets the distance
 * between those two points explicitly.
 *
 */

NUM ls_ls_dist(NUM a0[3], NUM b0[3], NUM c0[3], NUM d0[3]) {
  NUM foo[3],bar[3];
  NUM dist;
  ls_closest_approach(foo,bar,a0,b0,c0,d0);
  dist = cart_3d(foo,bar);
  if(!isfinite(dist))
    printf( "dist is nan\n");
  if(!dist)
    printf("dist is 0\n");
  return dist;
}

/**********************************************************************
 * fl_segment_masked_dist
 * Given two pointers to VERTEXes, return the closest approach of their
 * two line segments, masking space to include only the demesnes of the
 * first line segment.  That is to say, return the distance to the
 * closest point of the second line segment that is within the volume of space
 * defined by the two face planes of the segment.  If there is only one face
 * plane (either this VERTEX or its next neighbor is an endpoint),
 * then half of space is allowed.
 *
 * If there is no point within the space, or if the VERTEX is at the end
 * of a fluxon, then return -1.
 *
 * It uses a helper routine, mask_halfspace, that masks and truncates two
 * points to within one of the two bounding halfspaces of a line segment.
 *
 * The "sign" argument to mask_halfspace should be 1 or -1.
 *
 * mask_halfspace returns -1, 0, or 1 depending on whether it rejects both
 * points, does nothing, or truncates a point.  Only the "-1" value matters
 * at the moment.
 */

static int  mask_halfspace(VERTEX *v, NUM X0[3], NUM X1[3], int sign){
  NUM A[3], B[3];
  NUM Y0[3], Y1[3];
  NUM d0, d1;

  if(!v->prev || !v->next)
    return 0;

  /* Find the normal vector to the plane */
  diff_3d(A,v->next->x,v->prev->x);
  scale_3d(A,A,1.0/norm_3d(A));

  /* Translate the origin to the vertex */
  diff_3d(Y0,X0,v->x);
  diff_3d(Y1,X1,v->x);

  /* Get the dot products */
  d0 = sign<0 ? -inner_3d(Y0,A) : inner_3d(Y0,A);
  d1 = sign<0 ? -inner_3d(Y1,A) : inner_3d(Y1,A);

  if(d0 <= 0 && d1 <= 0)
    return -1;

  if(d0 > 0 && d1 > 0)
    return 0;

  diff_3d(B, Y1,Y0);
  scale_3d(B,B,d0/(d0-d1));

  if(d0 <= 0)
    sum_3d(X0,X0,B);
  else /* if(d1 <= 0) */
    sum_3d(X1,X0,B);

  return 1;
}

NUM fl_segment_masked_dist(VERTEX *v0, VERTEX *v1) {
  NUM P0[3],P1[3];
  return fl_segment_masked_deluxe_dist(P0, P1, v0, v1);
}


/**********************************************************************
 * fl_segment_masked_deluxe_dist
 * You feed in an X0 and X1, and they get stuffed with the points of closest
 * approach of the segments.  You also get the distance back.
 *
 * This function is deprecated -- the smoothly masked fl_segment_deluxe_dist is
 * preferred.
 */
  NUM X0[3], X1[3];
NUM fl_segment_masked_deluxe_dist(NUM P0[3], NUM P1[3], VERTEX *v0, VERTEX *v1){
  int i;

  /* Make sure that we've got a valid line segment */
  if(!v0 || !v1 || !v0->next || !v1->next)
    return -1;

  if(v0==v1 || v0==v1->next || v0->next == v1)
    return 0;

  /* Make local copies of the points before truncating! */
  cp_3d(X0,v1->x);
  cp_3d(X1,v1->next->x);

  if( (mask_halfspace(v0,       X0,X1, 1) < 0) ||
  (mask_halfspace(v0->next, X0,X1,-1) < 0) )
  return -1;

  ls_closest_approach(P0, P1, v0->x, v0->next->x, X0, X1);

  return cart_3d(P0, P1);
}

/**********************************************************************
 * fl_segment_dist
 * fl_segment_deluxe_dist
 * Given two pointers to VERTEXes, return the closest approach of their
 * two line segments.  If either is invalid, return 1e50.
 * The deluxe_dist accepts two points, which get filled with the
 * points of closest approach of the two line segments.
 *
 * Does a projection in the inverse sine to get the output
 * distance. P0 and P1 are not projected.
 *
 */
NUM fl_segment_dist(VERTEX *v0, VERTEX *v1) {
  NUM P0[3],P1[3];
  return fl_segment_deluxe_dist(P0,P1,v0,v1);
}

NUM fl_segment_deluxe_dist(NUM P0[3],NUM P1[3], VERTEX *v0, VERTEX *v1) {
  NUM a,b,Plen;
  NUM P[3], seg[3], cr[3];
  WORLD *w = v0->line->fc0->world;

  /* Exclude trivial cases, adjacent-segment case, and next-nearest-segment case. */
  /* (farther segments aren't as likely to cause trouble) */
  if(!v0 || !v1 || !v0->next || !v1->next)
    return 1e50;

  if(v0==v1 ||   (v0->next == v1) ||  (v1->next == v0) )
    return 1e50;

  ls_closest_approach(P0,P1,v0->x,v0->next->x,v1->x,v1->next->x);

  diff_3d(P,P1,P0);
  Plen = norm_3d(P);
  if(Plen==0)
    return 1e50;

  // inverse square sine...
  /* Scale by the inverse square sine of the projection angle:  */
  /* AB x (closest).  Note that this makes the metric asymmetric! */

  /* See also reflect_deluxe, which undoes the sin^2 for reflected (image) */
  /* segment distances -- if you change this you must also change that...  */

  scale_3d(P,P,1.0/Plen);

  diff_3d(seg,v0->next->x,v0->x);
  scale_3d(seg,seg,1.0/norm_3d(seg));

  cross(cr,P,seg);
  a = norm2_3d(cr); // this is sin^2 theta...
  //  a *= a;          // sin^2 --> sin^4

  // Now scale by the limited skew angle between the two segments if our World is set for that.
  if(w->handle_skew) {
    NUM alpha;
    NUM skewlimit_recip = 4;
    NUM s1[3];
    NUM mat[9];


    projmatrix(mat, P0, P1);

    diff_3d(P, v1->next->x, v1->x);
    mat_vmult_3d(s1, mat, P);

    mat_vmult_3d(P, mat, seg);

    alpha = 2.0  / ( fabs(   ( norm_2d(s1) * norm_2d(P))   / inner_2d(s1, P) ) + skewlimit_recip ) ;
    if(a == 0)
      return 1e50;

    Plen *= alpha/a;
  } else {
    if(a == 0)
      return 1e50;
    Plen /= a;
  }


  return Plen;
}

/**********************************************************************
 **********************************************************************
 ***** Projection
 */

/**********************************************************************
 * projmatrix
 * Given two points in 3-space, return a rotation matrix that
 * rotates the vector connecting them into the Z axis (for projecting
 * stuff down to the plane perpendicular to the given line segment).
 * The returned matrix has a determinant of unity, by construction, so
 * distances are preserved.
 *
 * You supply the output location, pre-allocated with 9 NUMs.
 *
 * The product of two rotation matrices is hand-assembled here:
 * first the vector is rotated about the Z axis into the ZX plane;
 * then it's rotated about the Y axis into the Z axis.  See
 * DeForest notebooks, vol. IV-A, p. 177.
 */
void projmatrix(NUM *out, NUM *x0_3, NUM *x1_3) {
  NUM a[3];
  NUM r2d,r3d;
  NUM m1[9];
  NUM m2[9];
  NUM *m;

  if(x0_3) {
    if(x1_3) {
      diff_3d(a,x1_3,x0_3);
    } else {
      a[0] = -x0_3[0]; a[1]= -x0_3[1]; a[2]= -x0_3[2];
    }
  } else {
    if(x1_3) {
      a[0] = x1_3[0]; a[1]= x1_3[1]; a[2]=x1_3[2];
    } else {
      int i;
      for(i=0;i<9;i++) {out[i]=((i % 4) == 0);}
      return;
    }
  }

  r2d = norm_2d(a);
  r3d = norm_3d(a);
  if(r2d==0) {
    m = out;
    if(a[2]>=0) {
      *(m++) = 1; *(m++) = 0; *(m++) = 0;
      *(m++) = 0; *(m++) = 1; *(m++) = 0;
      *(m++) = 0; *(m++) = 0; *(m++) = 1;
    } else {
      *(m++) = -1; *(m++) = 0; *(m++) =  0;
      *(m++) =  0; *(m++) = 1; *(m++) =  0;
      *(m++) =  0; *(m++) = 0; *(m++) = -1;
    }
    return;
  }


  /* Generate the rotation matrices independently and multiply them.
     Slightly slower than direct assembly, but robust. */

  /* rotate to remove the y component */
  m = m1;
  *(m++) = a[0] / r2d;  *(m++) = a[1] / r2d;   *(m++) = 0;
  *(m++) = -a[1] / r2d;  *(m++) = a[0] / r2d;    *(m++) = 0;
  *(m++) = 0;           *(m++) = 0;             *(m)   = 1;

  /* rotate the vector (now in the xz plane) to point along the z axis */
  m=m2;
  *(m++) = a[2] / r3d;  *(m++) = 0;            *(m++) = -r2d / r3d;
  *(m++) = 0;            *(m++) = 1;            *(m++) =  0;
  *(m++) = r2d / r3d;   *(m++) = 0;            *(m)   = a[2]/r3d;
  mat_mult_3d(out,m2,m1);

}

/**********************************************************************
 * reflect
 * Given a PLANE and a point, reflect the point through the PLANE.
 * The reflected point is stuffed into the target array.
 */
void reflect(NUM out[3], NUM point[3], PLANE *plane) {
  NUM foo[3];
  NUM l;
  diff_3d(foo,point,plane->origin);
  scale_3d(foo, plane->normal, -2 * inner_3d(foo, plane->normal));
  sum_3d(out, foo, point);
}


/**********************************************************************
 **********************************************************************
 ***** Voronoi-cell math
 */

/**********************************************************************
 * perp_bisector_2d
 *
 * Find the bisector of the line between (P) and (Q) [this is trivial]
 * and return the bisector and slope in the given temporary 3-NUM array.
 * Return  0 on normal completion; 1 if the slope is infinite.
 *
 * If you feed in NULL for P, then P is treated as the origin.
 */
int perp_bisector_2d(NUM *out, NUM *P, NUM *Q) {
  NUM *delta;
  NUM scratch[2];

  if( !P ) {
    /* P-is-the-origin case */
    out[0] = Q[0]/2;
    out[1] = Q[1]/2;

    out[2] = - Q[0]/Q[1]; /* Perp. slope = (- delta-X / delta-Y). */
    return !isfinite(out[2]);
  }

  /* P-not-the-origin case */
  out[0] = ( Q[0] + P[0] ) / 2;
  out[1] = ( Q[1] + P[1] ) / 2;

  /* Slope is negative reciprocal */
  out[2] = - ( (Q[0]-P[0]) / (Q[1]-P[1]) );
  return !isfinite(out[2]);
}


/**********************************************************************
 * intersection_2d
 * Given L1 and L2 as 2-D lines (3-vectors including a 2-vector
 * of a point on the line, and a slope (the infinite case is handled)),
 * stuff the coordinates of the intersection into the return location and
 * return 0, or, if the lines are parallel, stuff NaN's into the coordinates
 * and return 1.
 *
 * There is probably a faster way to do this -- this one is the first
 * one that came to mind.  I calculate the difference in slopes, then
 * calculate the offset from the base point for L1.  The infinite-slope
 * cases are handled separately.
 */
int intersection_2d(NUM *out, NUM *L1, NUM *L2) {
  if( L1[2]<1e100 && L1[2]>-1e100 && L2[2]<1e100 && L2[2]>-1e100 ) {
    NUM delta_y;
    NUM delta_x;
    NUM delta_s;

    delta_s = L1[2] - L2[2];  /* Delta slope */

    if(!delta_s) {            /* Delta is zero -- parallel lines */
      out[1] = out[0] = NAN;
      return 1;
    }

    delta_y = L2[1] + L2[2] * ( L1[0] - L2[0] ) - L1[1]; /* Delta-Y at L1's X0 */
    delta_x = delta_y / delta_s; /* Offset from L1's X0 to intersection */

    out[0] = L1[0] + delta_x;
    out[1] = L1[1] + delta_x * L1[2];
    return 0;
  }
  else if( isfinite(L1[2]) ) {  /* L2 is vertical; L1 is not */
    out[0] = L2[0];
    out[1] = L1[1] + (L2[0] - L1[0]) * L1[2];
    return 0;
  } else if( isfinite(L2[2]) ) { /* L1 is vertical; L2 is not */
    out[0] = L1[0];
    out[1] = L2[1] + (L1[0] - L2[0]) * L2[2];
    return 0;
  }

  /* Parallel vertical lines */
  out[1] = out[0] = NAN;
  return 1;
}


/**********************************************************************
 * project_n_fill
 *
 * Given a VERTEX and a DUMBLIST of VERTICES, project the DUMBLIST
 * down to 2-D using fl_segment_deluxe_dist (closest approach of
 * fluxels associated with vertices, not the vertices themselves) to
 * find the projected radius, and fill the angle and radius fields of
 * the members of the DUMBLIST.
 *
 * Also accumulate the closest-real-approach distance between neighbors
 * and the current VERTEX.  By examining only neighbors that have been
 * winnowed with the special metric, we risk missing the actual
 * minimum distance -- but not very often.  Hopefully it won't lead to
 * pathological evolution problems!
 *
 */
void project_n_fill(VERTEX *v, DUMBLIST *horde) {
  int i;
  NUM pm[9];
  NUM X0[3];
  NUM p0[3], p1[3];
  NUM len, r;
  VERTEX *v1;
  char crunch = 0;

  if(v->line->fc0->world->verbosity >= 5)
    printf("Entered project_n_fill... got a horde of %d candidates...\n",horde->n);

  if(!v->next) {
    fprintf(stderr,"Hey!  Don't feed end vertices to project_n_fill....\n");
    fflush(stderr);
    return;
  }

  if(v->line->fc0->world->verbosity >= 5)
    printf("calling projmatrix with (%g,%g,%g) and (%g,%g,%g)\n",v->x[0],v->x[1],v->x[2],v->next->x[0],v->next->x[1],v->next->x[2]);

  projmatrix(pm,v->x,v->next->x);

  if(v->line->fc0->world->verbosity >= 5)
    printf("back\n");

  v->r_cl = -1; /* Initialize closest-approach accumulator */

  for(i=0;i<horde->n;i++) {

    if(v->line->fc0->world->verbosity >= 5) {
      printf("i=%d ",i);
      fflush(stdout);
    }


    v1 = (VERTEX *)((horde->stuff)[i]);

    r = fl_segment_deluxe_dist(p0, p1, v, v1); /*This is where the
                                                 non-linear projection
                                                 is done. It is for r
                                                 only, not for p0 and
                                                 p1*/

    if((v->line->fc0->world->verbosity - (r>=0))>=6)
      printf("\nfl_segment_deluxe_dist returned %g for vertex %ld\n",r,v1->label);

    if(r<=0 || !isfinite(r)) {
      horde->stuff[i] = 0;
      crunch=1;
    }
    else {

      if(r<v->r_cl || v->r_cl < 0 ) /* Accumulate closest approach distance */
        v->r_cl = r;

      diff_3d(X0, p1, p0);

      mat_vmult_3d(v1->scr,pm,X0); /* Project into the perpendicular plane */

      len = norm_2d(v1->scr);      /* * 2-D * length of vector */

      if(v->line->fc0->world->verbosity >= 5)
        printf("len=%g for vertex %ld\n",len,v1->label);

      if(len <= 0) {

        if(v->line->fc0->world->verbosity == 4)
          printf("len=%g for vertex %ld\n",len,v1->label);

        horde->stuff[i] = 0;
        crunch=1;
      }
      else {
        scale_3d(v1->scr,v1->scr,r/len); /*so scr is the vector of the
                                           nearest approach with a
                                           length of the inverese
                                           sin^2 projected distance
                                           that has been projected
                                           into a plane perpendicular
                                           to the v-fluxel.*/

        v1->a = ATAN2(v1->scr[1],v1->scr[0]); //store angles
        v1->r = r; //store non-linear projected distances.
      }
    }

  }

  /* If some of the items were invalidated, crunch down the horde... */
  if(crunch) {
    int j;
    for(i=j=0;i<horde->n;i++) {
      for(;i<horde->n && !horde->stuff[i]; i++)
        ;
      if(i<horde->n)
        horde->stuff[j++] = horde->stuff[i];
    }
    horde->n = j;
  }
}


/**********************************************************************
 * sort_by_angle_2d
 *
 * Given a DUMBLIST of VERTICES containing 2-D points in their
 * scratch spaces, sort them by angle from the origin.
 *
 * No actual checking of the vector in ->scr is done -- you must
 * pre-set the angle in ->a!  (This is normally done by project_n_fill.)
 */

static int angle_cmp(void *a, void *b) {
  if(a==b) return 0;

  TRIM_ANGLE( ((VERTEX *)a)->a );
  TRIM_ANGLE( ((VERTEX *)b)->a );

  if(((VERTEX *)a)->a == ((VERTEX *)b)->a ) {


    if( ((VERTEX *)a)->r==((VERTEX *)b)->r ) {
      return 0;
    }
    else {
      return ( ((VERTEX *)a)->r < ((VERTEX *)b)->r ) ? -1 : 1;
    }


  } else {
    return ( ((VERTEX *)a)->a < ((VERTEX *)b)->a ) ? -1 : 1;
  }
}

void sort_by_angle_2d(DUMBLIST *horde) {
  int i;

  /* Do finiteness checking... */
  for(i=0;i<horde->n;i++) {
    VERTEX *v1 = (VERTEX *)(horde->stuff[i]);
    if(!v1 || !isfinite(v1->a) || !isfinite(v1->r)) {
      //      printf("sort_by_angle_2d: v1 is bad (label %d)\n",v1?v1->label:0);
      /*** Break the symmetry ***/
      v1->r=1e50; v1->a=0;
      dumblist_rm(horde,i);
      //      printf("X: v=%d\tx=(%g,%g,%g)\tr=%g\ta=%g\n",v1->label,v1->x[0],v1->x[1],v1->x[2],v1->r,v1->a); fflush(stdout);
      i--;
    }
  }

  /* Do the actual sort... */
  dumblist_sort(horde,angle_cmp);
  //dumblist_crunch(horde); (crunching is handled by the sorter)
}


/**********************************************************************
 * hv_cp - copy a HULL_VERTEX
 */
 void hv_cp(HULL_VERTEX *to, HULL_VERTEX *from) {
  to->p[0] = from->p[0];
  to->p[1] = from->p[1];
  to->a_l = from->a_l;
  to->a_r = from->a_r;
  to->bisector[0] = from->bisector[0];
  to->bisector[1] = from->bisector[1];
  to->bisector[2] = from->bisector[2];
  to->open = from->open;
}


/**********************************************************************
 * hull_2d
 *
 * Given an angle-sorted DUMBLIST of VERTICES containing 2-D points in their
 * scratch spaces, find the minimal convex hull of the origin.
 * (note that this is different than the conventional "maximal convex hull"
 * problem). The original DUMBLIST winds up containing the points that were
 * used for the hull; you can also pass in a DUMBLIST
 * that gets the rejects.  On exit, the output array (which you
 * have to supply, pre-allocated) contains N HULL_VERTEXes
 * (N is the number of vertices in the final hull).
 * Each voronoi vertex is the one to the left of the corresponding
 * neighbor point.
 *
 * Typical usage:
 *    <collect dumblist>
 *    vertices = malloc(sizeof(NUM)*2*horde->n);
 *    hull_2d(vertices,horde,rejects=dumblist_new());
 *
 * If you don't want a list of rejects, pass in NULL for the
 * third parameter.
 *
 * Note:  (To avoid the hassle of frequently removing fields from the
 * horde DUMBLIST, the code nulls individual rejected elements, then
 * crunches 'em all at the end.)
 *
 * This is probably the hot spot for the code as a whole, so some
 * perhaps misguided trades have been made against readability and
 * for (hopefully) speed.
 *
 * The sort-by-angle step requires that the projected radius and
 * angle fields in the VERTEX be filled; this is performed by
 * project_n_fill, which should be called before hull_2d.
 * *
 */

#define MOD_INC(a,n) ( ( (++(a)) < (n) ) ? (a) : ((a)=0) )
#define MOD_DEC(a,n) ( ( (--(a)) >= 0  ) ? (a) : ((a)=(n)-1) )
#define MOD_NEXT(a,n) ( ((a)<((n)-1)) ? ((a)+1) : 0 )
#define MOD_PREV(a,n) ( ((a)>0) ? ((a)-1) : ((n)-1) )

void hull_2d(HULL_VERTEX *out, DUMBLIST *horde, DUMBLIST *rejects) {
  int n = horde->n;
  int i, i_r, i_l; /* Current, right, and left indices */
  int terminus; /* When we get here, we are done (moving target) */
  char been_there;
  VERTEX *iv, *rv, *lv;
  NUM a;
  char abort = 0;

  /* Do nothing if we get a trivial list */
  if(horde->n < 1)
    return;

  if(horde->n == 1) {
    perp_bisector_2d( &(out[0].bisector[0]),0,((VERTEX *)(horde->stuff[0]))->scr);
    out[0].open = 1;
    out[0].p[0] = out[0].p[1] = 0;
    out[0].a_l = ((VERTEX *)(horde->stuff[0]))->a + PI/2;
    out[0].a_r = ((VERTEX *)(horde->stuff[0]))->a - PI/2;
    return;
  }

  sort_by_angle_2d(horde);

  /* Get ready for main loop -- on entry, the first two candidates
   * need their perpendicular bisectors calculated...
   */
  n=horde->n;
  terminus = i = been_there = 0;
  i_r = n-1;
  iv = (VERTEX *)(horde->stuff[ i ]);
  rv = (VERTEX *)(horde->stuff[i_r]);

  perp_bisector_2d( &(out[ i ].bisector[0]), 0, iv->scr );
  perp_bisector_2d( &(out[i_r].bisector[0]), 0, rv->scr );
  intersection_2d( out[i_r].p, out[i_r].bisector, out[ i ].bisector );
  a = iv->a - rv->a; TRIM_ANGLE(a);
  if(a < 0 ||
     !isfinite(out[i_r].p[0]) || !isfinite(out[i_r].p[1])) {
    out[i_r].p[0] = out[i_r].p[1] = 0;
    out[i_r].open = 1;
  } else
    out[i_r].open = 0;

  /* Main loop */
  do {
    NUM a,b;
    if(horde->stuff[i]) {    /* Skip over already-zeroed elements */

      char open=0;

      /*** Find next guy and prep his bisector if necessary ***/
      for(i_l=MOD_NEXT(i,n);  i_l != i && !horde->stuff[i_l]; MOD_INC(i_l,n))
        ;
      lv = (VERTEX *)(horde->stuff[i_l]);

      // Remove dupes
      if( iv==rv )
        goto reject;

      if(!been_there)
        perp_bisector_2d( out[i_l].bisector, 0, lv->scr );

      /*** Check for pathologies ***/
      if(i_r == i || i_l == i) {
        abort = 1;
        //        fflush(stdout);
        //        fprintf(stderr,"hull_2d: eliminated all vertices! I give up.\n");
        //        horde->n=0;
        //        return;
      }

      if(!isfinite(iv->r) || !isfinite(iv->a))
        goto reject;

      if(iv->r == 0)
        goto reject;

      /*** Check colinearity and opentude on right & left ***/
      /* (right could be handled by caching, but it's pretty cheap) */
      b = iv->a - rv->a;
      TRIM_ANGLE(b);
      if( b < EPSILON && b > -EPSILON && rv->r < iv->r ) {
        goto reject; /* Colinear with (and farther than) right-side vertex */
      }

      a = lv->a - iv->a;

      TRIM_ANGLE(a);
      if( a < EPSILON && a > -EPSILON && lv->r <= iv->r ) {
        goto reject;  /* Colinear with (and farther than) left-side vertex */
      }

      if( a > EPSILON && a < PI-EPSILON) {
        /* Line isn't open on the left, so it intersects on the left.   */
        intersection_2d( out[ i ].p, out[ i ].bisector, out[i_l].bisector );
        out[i].open = 0;


        if( b > EPSILON && b < PI-EPSILON) {

          /* Line isn't open on the right either, so compare intersections */
          a = cross_2d( out[i_r].p, out[ i ].p );

          if( !isfinite(a) ) {
            printf("ASSERTION FAILED in hull_2d: intersection should exist (but doesn't)\n\t\tleft_p: %5.2g, %5.2g ;   right_p: %5.2g, %5.2g",out[i].p[0],out[i].p[1],out[i_r].p[0],out[i_r].p[1]);
            printf("\n\t\tleft_b: %5.2g, %5.2g, %5.2g;  this_b: %5.2g, %5.2g, %5.2g; right_b: %5.2g, %5.2g, %5.2g\n",out[i_l].bisector[0],out[i_l].bisector[1],out[i_l].bisector[2],out[i].bisector[0],out[i].bisector[1],out[i].bisector[2],out[i_r].bisector[0],out[i_r].bisector[1],out[i_r].bisector[2]);
            printf("\t\tlv->a=%5.2g, iv->a=%5.2g, rv->a=%5.2g,   lv-iv: %5.2g,   iv-rv: %5.2g\n",lv->a*180/PI,iv->a*180/PI,rv->a*180/PI,(lv->a-iv->a)*180/PI,(iv->a-rv->a)*180/PI);
            goto reject;
          }

          /*** If intersections were in the wrong order, reject the point ***/
          if(a < 0) {
            goto reject;
          }
        }

      } else {
        /* Line is open on the left -- set the left point to zero */
        out[i].p[0] = out[i].p[1] = 0;
        out[i].open = 1;
      }


      /*** Acceptance code ***/
      been_there = (been_there || (i_l < i));
      i_r = i;
      i = i_l;
      rv = iv;
      iv = lv;

      /*** The only way to get to the rejection code is to goto it! ***/
      if(0) {
      reject:
        /*** rejection code ***/

        /* "Normally", we back up when we hit a rejection -- but the
         * first item in the list is an exception -- in that case we
         * keep the same right side but walk forward one.  That requires
         * updating the intersection information for the right side.
         */
        if(rejects)
          dumblist_add(rejects,horde->stuff[i]);

        horde->stuff[i] = 0;

        if(been_there || (i > i_r)) {
          i = i_r;

          iv = rv;

          for(i_r=MOD_PREV(i,n); i_r!=i && !horde->stuff[i_r]; MOD_DEC(i_r,n))
            ;
          terminus =i_r;
          rv = (VERTEX *)(horde->stuff[i_r]);


        } else {

          i = i_l;
          terminus = i_r;
          iv = lv;

          been_there = 0;

          a = iv->a - rv->a; TRIM_ANGLE(a);
          if(a<0) {
            out[i_r].p[0] = out[i_r].p[1] = 0;
            out[i_r].open = 1;
          }  else {
            intersection_2d( out[i_r].p, out[i_r].bisector, out[ i ].bisector );
            out[i_r].open = 0;
          }
        }
      }
    } /* end of non-zeroed-out check */

  } while(!abort && (!been_there || i != terminus));  /* End of main loop */

  /* Finished -- now crunch the dumblist and, simultaneously, the output. */
  {
    int j;

    for(i=j=0; i<horde->n; i++) {
      VERTEX *hsi = (VERTEX *)(horde->stuff[i]);
      if(hsi) {
        if(i != j) {
          horde->stuff[j] = hsi;
          out[j] = out[i];
        }
        j++;
      }
    }
    horde->n = j;
  }

  /* Finally - insert appropriate atan2 angular fields into the output */
  {
    HULL_VERTEX *hv = out;
    for(i=0;i<horde->n;i++,hv++) {
      if(!hv->open)
        hv->a_r = hv->a_l = ATAN2(hv->p[1], hv->p[0]);
      else {
        hv->a_l =(((VERTEX *)(horde->stuff[i]))->a ) + PI/2;
        hv->a_r =(((VERTEX *)(horde->stuff[MOD_NEXT(i,horde->n)]))->a) - PI/2;
      }
    }
  }
}

/**********************************************************************
 * atan2_oct - sleazy-but-fast atan2
 * (not quite a drop-in replacement for atan2 -- be careful when using it.
 * Should only be used via the ATAN2 macro!)
 * Faster than atan2.  Most of what we use atan2 for is sorting by angle, but
 * for that we only need a monotonically increasing function.  atan2_slz uses
 * the octagonal approximation -- the unit circle is approximated as an octagon
 * with vertices on the axes and on the 45-degree lines.
 *
 * The range is -pi/2 to pi/2.
 *
 */

NUM atan2_oct( NUM y, NUM x) {
  NUM ay = fabs(y);
  NUM ax = fabs(x);
  NUM out;

  if(ax==0)
    return(y>0?PI/2: (y<0)?-PI/2 : nan("1"));

  /* fast atan2 in one quadrant using octagonal approx */
  out = 0.25 * PI * ((ax >= ay) ? (ay/ax) : (2 - ax/ay));

  if(x>0) {
    if(y<0)
      out *= -1;
  } else {
    if(y>=0) {
      out = PI - out;
    } else {
      out -= PI;
    }
  }
  return out;
}


/******************************
 * check_hullpoint - helper routine for hull_2d_us
 * Given a VERTEX and HULL_VERTEX pointer for next and prev, and
 * a VERTEX to test, returns true if it should be kept and false if not.
 * The return value is negative if the vertex is closed, 1 if it is open
 * on the next side, 2 if it is open on the prev side.
 * If the hullpoint is to be kept, then the previous hull_vertex has its
 * open and/or p fields updated as appropriate.
 */

int check_hullpoint(VERTEX *v,
                    VERTEX *pv,
                    HULL_VERTEX *phv,
                    VERTEX *nv,
                    HULL_VERTEX *nhv,
                    HULL_VERTEX *scr
                    ) {
  NUM point[2];
  NUM a;
  char open = 0;

  /* Check for identity -- very easy fix if we're identical */
  // (Why is this here?  The passno mechanism, below, should prevent ever triggering it)
  if(v==pv || v==nv) {
    if (v->line->fc0->world->verbosity >= 1) {
    fprintf(stderr,"Hmmm, check_hullpoint duplication check activated, pv:%ld, v:%ld, nv:%ld. Could be nearly cospatial.\n",pv->label,v->label,nv->label);
    /* This gets caught in for the third horde vertex. The problem
       initially occured in the second horde bvertex where the first
       and second horde vertices were nearly cospatial. The colinear
       check didn't work correcly in that case because both sides
       needed to be open.*/
    }
    return 0;
  }

  /* Check for colinearity -- easy fix if we're colinear */
  /* The colinear_keep flag gets return value that will be given if
   * the vertex turns out to be open - this is sort of the opposite
   * sense of the return value.
   */

  {
    int colinear_keep = 0;
    NUM eps = EPSILON * v->r;
    /* if the angle between them in small (sin(theta) near zero,
       cos(theta) is positive) */
    if( fabs( cross_2d(v->scr, pv->scr) ) < eps && (inner_2d(v->scr, pv->scr) > 0) ) {
      if(  pv->r >= pv->r  ) {
        return 0;
      }
      else
        colinear_keep = 1;
    }

    if( fabs( cross_2d( v->scr, nv->scr) ) < eps && (inner_2d(v->scr, nv->scr) > 0) )  {
      if(  v->r >= nv->r  ) {
        return 0;
      }
      else
        colinear_keep = 2;
    }

    if(colinear_keep) {
      perp_bisector_2d(scr->bisector, 0, v->scr);
      intersection_2d(point, scr->bisector, phv->bisector);
      intersection_2d(scr->p, scr->bisector, nhv->bisector);
      phv->p[0] = point[0];
      phv->p[1] = point[1];
      /* Don't touch phv->open: it should already have the correct value */
      return ( phv->open ? colinear_keep : -1 );
    }
  } /* end of colinear check convenience block */

  /* Not colinear -- find perpendicular bisector and start crunching... */

  perp_bisector_2d(scr->bisector, 0, v->scr);

  /* scr->p gets the hull intersection point with the next guy; */
  /* point gets the hull intersection point with the previous guy */

  intersection_2d( scr->p, scr->bisector, nhv->bisector );
  intersection_2d( point, scr->bisector, phv->bisector );

  /**** Open on next side ? ****/
  if(          cross_2d(  v->scr, nv->scr) <= 0) {
    if(        cross_2d(  pv->scr, v->scr) <= 0) {
      // Antiparallel case (open both sides) -- keep
      open = 3;
    } else {
      // Noncolinear and open on next side -- keep
      open = 1;
    }
  }
  /**** Open on previous side ? ****/
  else if(     cross_2d( pv->scr,  v->scr) <= 0) {
    // Noncolinear and open on prev side -- keep
    open = 2;
  }

  /**** Closed and OK? ****/
  else if(  cross_2d( point, scr->p ) >= 0  ) {
    // Noncolinear and closed, but OK -- keep
    open = -1;
  }

  /**** Closed and not OK? ****/
  else {
    /* closed both sides, but intersections in wrong order -- reject */
      return 0;
  }

  /* are the points the same? If so, reject it. */
  /*if (scr->p[0] == point[0] && scr->p[1] == point[1]) {
    return 0;
    }*/

  /*** If we get here we're keeping, so we have to copy the ***/
  /*** intersection point into the previous hull_vertex.    ***/
  phv->p[0] = point[0];
  phv->p[1] = point[1];
  phv->open = (open==2)?1:0;

  return open;

}



/**********************************************************************
 * hull_2d_us - hull routine for unsorted points
 * (see if we can avoid the sorting step)
 *
 * Keep a current-hull dumblist and add points to it by brute-force
 * copying and binary searching.  Brute force should be OK since the
 * typical hull size is only 6.
 *
 * We keep our own workspace, and build the hull in it, then copy the
 * elements by brute force into the horde (thereby truncating it considerably!)
 *
 * See notes on hull_2d, above -- the points should be pre-treated with
 * project_n_fill, which sets the projected 2-D coordinates in the scr field of
 * each candidate VERTEX, and also the a and r fields.  (a may contain the
 * sleazy arctan rather than the true arctan -- it only has to be monotonic).
 *
 */
static DUMBLIST *ws = 0;


// Local helper plasmoid_conjugate identifies nodes that are "adjacent"
// to the start and end nodes of a plasmoid fluxon
static inline char plasmoid_conjugate(VERTEX *a, VERTEX *b) {
  if((a->line != b->line) || (!a->line->plasmoid))
    return 0;
  if( a == a->line->start )
    return (b == b->line->end->prev || b==b->line->end->prev->prev);
  if( a == a->line->start->next )
    return (b == b->line->end->prev);
  if( b == a->line->start )
    return (a == a->line->end->prev || a==a->line->end->prev->prev);
  if( b == a->line->start->next )
    return (a == a->line->end->prev);
  return 0;
}

void hull_2d_us(HULL_VERTEX *hull, DUMBLIST *horde, VERTEX *central_v) {
  int n = horde->n;
  int horde_i;
  long passno;
   WORLD *w = 0;
  int i,j;

  if(ws==0) {
    ws = new_dumblist();
    dumblist_grow(ws,16);
  }

  if(! horde->n)
    return;

  /*go through each potential neighbor and make sure it is in the
    world still. Old problem that is now fixed, can do it an easier
    way now.
  for(i=0;!w && i<horde->n; i++)
    if( (((VERTEX **)(horde->stuff))[i])->line )
    w = (((VERTEX **)(horde->stuff))[i])->line->fc0->world;*/

  w = (((VERTEX **)(horde->stuff))[0])->line->fc0->world;
  if(!w) {
    fprintf(stderr,"You're in trouble, guv! exiting (horde->n was %d)...\n",horde->n);
    exit(99);
  }

  // use passno to prevent multiple looks at the same candidate.
  passno = ++w->passno;

  /* Set our output list to zero length */
  ws->n = 0;

  /* Nothing in the horde?  We're done! */
  if(horde->n == 0)
    return;


  /* always add the first eligible VERTEX in the list... */
  for(i=0;i<horde->n &&
        ( horde->stuff[i] == central_v ||              // skip the main vertex if present
          horde->stuff[i] == central_v->prev ||        // skip the following segment if present
          horde->stuff[i] == central_v->next ||        // skip the previous segment if present
          !( ((VERTEX *)(horde->stuff[i]))->next ) ||  // skip endpoint vertices if present
          plasmoid_conjugate( (VERTEX *)(horde->stuff[i]), central_v )
          );
      i++) //find first candidate
    ;

  /*first look at this candidate and add it, it will be removed later
    if it is not really a neighbor. */
  if(i<horde->n) {
    ((VERTEX **)(horde->stuff))[i]->passno = passno;
    ws->stuff[0] = horde->stuff[i];
    ws->n = 1;
    hull[0].open = 1; //because there are no others yet
    perp_bisector_2d(hull[0].bisector, 0, ((VERTEX *)(horde->stuff[i]))->scr);
  } else {
    horde->n = 0;
    return;
  }


  /** Loop over the horde and check vertices as we go... **/
  for(horde_i=1; horde_i<horde->n; horde_i++) {
    VERTEX *v =  ((VERTEX **)(horde->stuff))[horde_i];
    VERTEX *nv, *pv;
    static HULL_VERTEX scrhv;
    HULL_VERTEX *nh, *ph, *scr;
    int next_idx;
    char keep;
    int  flag;

    /* Check for skip conditions. This has the same conditions as
       above, but it rejects ws->stiff[0] because that passno has
       already been updated. */
    if(v->passno == passno ||                   // Already been here
       v == central_v ||                        // Skip ourselves if we encounter us
       !v->next ||                              // Must have a following segment to be valid
       ( v->line==central_v->line &&            // Several end conditions for same fluxon:
         ( v->next==central_v ||                // Skip next segment
           v->prev==central_v ||                // Skip previous segment
           plasmoid_conjugate( v, central_v)
           )
         )
       ) {
      // Skip (do nothing; just move on to the next candidate)
    } else {
      v->passno = passno; // Mark it processed (avoid duplicated effort)

      keep = 0;

      /* Linear searches are slow, but not so bad -- there are
       generally under 10 elements in ws, so binary search probably
       isn't worth the effort.  */
      /* Sort through the current candidates to find the first one
         that has a larger angle than v, this candidate will be nv. nv
         is assigned here. This loop determines next_idx which is also
         the index for nv.*/
      for(next_idx = 0;
          next_idx < ws->n &&
            (nv = (((VERTEX **)(ws->stuff))[next_idx]) )->a < v->a;  // assign to nv
          next_idx++
          )
        ;


      /******************************
       * Now set nv to the next vertex after this one and pv to the previous one.
       * next_idx is the true insertion point if this element is kept.
       */

      if( nv->a < v->a ) {
        /* We went off the end, so put ti at the end */

        next_idx = ws->n;
        nv = ((VERTEX **)(ws->stuff))[0];
        nh = hull;
        pv = ((VERTEX **)(ws->stuff))[ws->n-1];
        ph = hull + ws->n-1;
      } else {
        /* We didn't, it is at the beginning */
        if(next_idx==0) {
          nh = hull;
          pv = ((VERTEX **)(ws->stuff))[ws->n-1];
          ph = hull + ws->n-1;
        } else {
          /* It is somewhere in the middle */
          nh = hull + next_idx;
          pv = ((VERTEX **)(ws->stuff))[next_idx - 1];
          ph = hull + next_idx-1;
        }
      }

      flag = check_hullpoint(v, pv,ph, nv,nh, &scrhv);

      if(flag) {
        /* Make sure we have room in the workspace */
        if(ws->size <= ws->n)
          dumblist_grow(ws, ws->size * 1.5 + 10);

        /* Make room in the workspace dumb list and the hull
           list. backward loop, search backward and shift everything
           +1 until we get to next_idx. */
        for( j = ws->n; j>next_idx; j--){
          ws->stuff[j] = ws->stuff[j-1];
          hv_cp(hull+j, hull+j-1);
        }

        ws->n++;

        /* Set the insertion-point values */
        ws->stuff[next_idx] = v;
        hv_cp(hull+next_idx, &(scrhv));

        hull[next_idx].open = ( flag>0 && (flag & 1) ); /* logical && and bitwise & */
        hull[   ( next_idx ? next_idx : ws->n ) - 1  ].open =
          ( flag>0 && (flag & 2) ); /* logical && and bitwise & */

        /**************************
         * Now that we've put the new guy in, make sure he doesn't invalidate
         * his neighbors on the prev side...
         */
        {
          int n, p, pp, nn, nnn, check_flag;

          n = next_idx;              //n
          p = n; MOD_DEC(p, ws->n);  //n-1
          pp= p; MOD_DEC(pp, ws->n); //n-2
          nn = n; MOD_INC(nn, ws->n);//n+1
          nnn=nn; MOD_INC(nnn,ws->n);//n+2

          /* Walk backward until we find a still-keep-worthy vertex */
          while( p != n &&
                 ! (check_flag = check_hullpoint( (VERTEX *)(ws->stuff[p]),
                                                  (VERTEX *)(ws->stuff[pp]), hull + pp,
                                                  (VERTEX *)(ws->stuff[n]),  hull + n,
                                                  hull + p))
                 ) {

            p=pp;
            MOD_DEC(pp, ws->n); //decrement p and pp
          }

          hull[p].open = (check_flag > 0) && (check_flag & 1);
          hull[pp].open = (check_flag > 0) && (check_flag & 2);

          /* walk forward until we find a still-keep-worthy vertex */
          while( n != p &&
                 ! (check_flag = check_hullpoint( (VERTEX *)(ws->stuff[nn]),
                                                  (VERTEX *)(ws->stuff[n]),  hull + n,
                                                  (VERTEX *)(ws->stuff[nnn]), hull + nnn,
                                                  hull+nn
                                                  ) )
                 ) {

            nn=nnn;
            MOD_INC(nnn, ws->n); //increment nn and nnn
          }

          hull[nn].open = (check_flag > 0) && (check_flag & 1);
          hull[n].open = (check_flag > 0) && (check_flag & 2);

          /* Now p is the first keepworthy vertex, walking backward
             from n, and nn is the first keepworthy vertex, walking
             forward from n. Crunch down the list. */
          /* Case nomeclature: .=no longer neighbor, *=true
             neighbor. the line represents each vertex currently in
             ws.*/
          {

            int ii,jj, diff;

            if( p <= n-1 ) {
              if( nn <= n-1 ) {

                // ...n*****p...|.....  case
                for(ii=0, jj=nn; jj<=p; jj++,ii++) {
                  hv_cp( hull+ii, hull+jj );
                  ws->stuff[ii] = ws->stuff[jj];
                }
                hv_cp( hull+ii, hull +n );
                ws->stuff[ii] = ws->stuff[n];
                ii++;

              } else {

                // *****p...|...n**** case

                ii=p+1;
                if(ii == n) {
                  ii++;
                } else {
                  hv_cp( hull + ii, hull + n );
                  ws->stuff[ii] = ws->stuff[n];
                  ii++;
                }
                if(ii<nn) {
                  for( jj=nn; jj < ws->n; ii++,jj++) {
                    hv_cp(hull + ii, hull + jj);
                    ws->stuff[ii] = ws->stuff[jj];
                  }
                } else if(ii==nn) {
                  ii = ws->n;
                }
              }
            } else {

              // ....|...n****p.. case
              if(n>0) {
                hv_cp(hull, hull+n);
                ws->stuff[0] = ws->stuff[n];
              }
              ii=1;

              if( nn>1 ) {
                for( jj=nn; jj<=p; ii++, jj++) {
                  hv_cp( hull + ii, hull + jj );
                  ws->stuff[ii] = ws->stuff[jj];
                }
              } else {
                ii= p + 1;
              }

            }/* end of case testing */
            ws->n = ii;

          } /* end of crunching convenience block */
        } /* end of neighbor-testing convenience block */
      } else {
        /* no keeping -- just ignore this point and start the next loop */
      }
    } /* end of non-duplication test block */
  } /* end of horde loop */

  /* Calculate the (more expensive) true-atan angles for the hull points and for the
   * actual neighbors, since there are so few of them.
   */

  for(  i=0; i<ws->n; i++ ) {
    VERTEX *vv = ((VERTEX **)(ws->stuff))[i];

    if(hull[i].open) {
      VERTEX *vn = ((VERTEX **)(ws->stuff))[MOD_NEXT(i,ws->n)];
      hull[i].a_r = atan2(vn->scr[1],vn->scr[0]) - PI/2;
      hull[i].a_l = atan2(vv->scr[1],vv->scr[0]) + PI/2;
    } else {
      hull[i].a_r = hull[i].a_l = atan2(hull[i].p[1],hull[i].p[0]);
    }

    vv->a = atan2(vv->scr[1],vv->scr[0]);
    horde->stuff[i] = vv;
  }
  horde->n = i;
}





/**********************************************************************
 * find_vertex_by_location
 *
 * Given a location in 3-space, a world, and an (optional) VERTEX,
 * find the VERTEX closest to the point you want.
 *
 * If you use the "global" flag then an explicit search is performed
 * (slow); otherwise a simple minimization-by-step is performed.
 *
 */


/* A helper routine for the tree walk in the explicit case */
static VERTEX *fvbl_stash_v;
static POINT3D fvbl_stash_x;
static NUM fvbl_current_min;
static long fvbl_helper( void *v, int depth) {
  NUM dist;
  dist = cart2_3d(fvbl_stash_x, ((VERTEX *)v)->x);
  if(!V_ISDUMMY((VERTEX *)v) &&
     (fvbl_current_min < 0 || dist < fvbl_current_min)) {
    fvbl_stash_v = v;
    fvbl_current_min = dist;
  }
  return 0;
}



VERTEX *find_vertex_by_location(POINT3D x, WORLD *w, VERTEX *v, int global) {

  if(!w){
    fprintf(stderr,"FLUX find_vertex_by_location needs a WORLD!\n");
    return 0;
  }
  if(!v){
    v = w->vertices;
    if(!v) {
      fprintf(stderr,"This WORLD apparently has no vertices! Help!\n");
      return 0;
    }
  }

  /*** Global search -- explicitly walk through the whole sim ***/
  if(global) {
    cp_3d(fvbl_stash_x,x);
    fvbl_current_min = -1;
    fvbl_stash_v = 0;
    tree_walker(w->vertices,v_lab_of, v_ln_of, fvbl_helper, 0);
    return fvbl_stash_v;
  }

  /*** Local search -- minimize from the vertex ***/
  else {
    //printf("Starting at vertex %ld",v->label);
    NUM dist2 = cart2_3d(x,v->x);
    int done = 0;

    do {
      NUM d2;


      if( v->next &&
          ((d2 = cart2_3d(x, v->next->x)) < dist2)  /* assignment */
          ) {
        /* walk forward */
        v = v->next;
        dist2 = d2;
      } else if( v->prev &&
                 ((d2 = cart2_3d(x, v->prev->x)) < dist2) /* assignment */
                 ) {
        /* walk backward */
        v = v->prev;
        dist2 = d2;
      } else {
        /* search neighbors & nearby */
        int i;
        VERTEX *v_candidate = 0;
        NUM vcd2 = dist2;
        for(i=0; i < v->neighbors.n; i++) {
          if( (d2 = cart2_3d(x, ((VERTEX **)(v->neighbors.stuff))[i]->x)) < vcd2 ) {
            v_candidate = ((VERTEX **)(v->neighbors.stuff))[i];
            vcd2 = d2;
          }
        }
        for(i=0; i<v->nearby.n; i++) {
          if( (d2 = cart2_3d(x, ((VERTEX **)(v->nearby.stuff))[i]->x)) < vcd2 ) {
            v_candidate = ((VERTEX **)(v->nearby.stuff))[i];
            vcd2 = d2;
          }
        }

        if(v_candidate) {
          v = v_candidate;
          dist2 = vcd2;
        } else
          done = 1;
      }
    } while(!done);
  }
  return v;
}

/**********************************************************************
 * above_plane
 *
 * Given three noncolinear points (which define an oriented plane)
 * determine whether a 4th point is above or below the plane.  Returns
 * 1 if the 4th point is above or on, 0 if it is  below.
 */
int above_plane(POINT3D A, POINT3D B, POINT3D C, POINT3D X) {
  POINT3D AB, AC, AX;
  POINT3D ABxAC;
  diff_3d(AB, B, A);
  diff_3d(AC, C, A);
  diff_3d(AX, X, A);
  cross_3d(ABxAC, AC, AB);
  return (  inner_3d(ABxAC, AX) >= 0  );
}

/**********************************************************************
 * above_plane_ratio
 *
 * Given three noncolinear points (which define an oriented plane)
 * determine whether a 4th point is above or below the plane.
 * Returns the projection of this vector onto AX, as a fraction of the
 * lengths of AX and ABxAC. Positive / negative values are above / below the plane.
 */
int above_plane_ratio(POINT3D A, POINT3D B, POINT3D C, POINT3D X) {
  POINT3D AB, AC, AX;
  POINT3D ABxAC;
  NUM AXL, ABxACL;
  diff_3d(AB, B, A);
  diff_3d(AC, C, A);
  diff_3d(AX, X, A);
  AXL = norm_3d(AX);
  ABxACL = norm_3d(ABxAC);
  cross_3d(ABxAC, AC, AB);
  printf("\n %g", inner_3d(ABxAC, AX) / (AXL * ABxACL));
  return (  inner_3d(ABxAC, AX) / (AXL * ABxACL)  );
}

/**********************************************************************
 * opposite_plane
 * Given three points forming a plane,
 * determine if a fourth and fifth points are on opposite sides
 * of this plane. Return 1 for opposite sides, 0 for same sides (or *on* the plane).
 */
int opposite_plane(POINT3D A, POINT3D B, POINT3D C, POINT3D X, POINT3D Y) {
  POINT3D AB, AC, AX, AY;
  POINT3D ABxAC;
  diff_3d(AB, B, A);
  diff_3d(AC, C, A);
  diff_3d(AX, X, A);
  diff_3d(AY, Y, A);
  cross_3d(ABxAC, AC, AB);
  return ( (inner_3d(ABxAC, AX) > 0) ^ (inner_3d(ABxAC, AY) > 0) );
}

/**********************************************************************
 * in_simplex
 * Given four noncoplanar points (which define a 3-D simplex),
 * determine whether a 5th point is inside or outside the simplex.
 */
int in_simplex( POINT3D P0, POINT3D P1, POINT3D P2, POINT3D P3, POINT3D X) {
  POINT3D P01, P02, P03;

  return ( ! (  (above_plane(P0,P1,P2,P3) ^ above_plane(P0,P1,P2,X) ) ||
                (above_plane(P1,P2,P3,P0) ^ above_plane(P1,P2,P3,X) ) ||
                (above_plane(P2,P3,P0,P1) ^ above_plane(P2,P3,P0,X) ) ||
                (above_plane(P3,P0,P1,P2) ^ above_plane(P3,P0,P1,X) )
              )
           );
}

/**********************************************************************
 * in_simplex_ratio
 * Given four noncoplanar points (which define a 3-D simplex),
 * determine whether a 5th point is inside or outside the simplex,
 * utilizing the means of above_plane_ratio
 *
 * CL - Perhaps only apply the ratio version to the X points
 */
int in_simplex_ratio( POINT3D P0, POINT3D P1, POINT3D P2, POINT3D P3, POINT3D X) {
  POINT3D P01, P02, P03;
  NUM tol;

  tol = 0.2;

  return ( ! (  (above_plane(P0,P1,P2,P3) ^ (above_plane_ratio(P0,P1,P2,X)>=tol) ) ||
                (above_plane(P1,P2,P3,P0) ^ (above_plane_ratio(P1,P2,P3,X)>=tol) ) ||
                (above_plane(P2,P3,P0,P1) ^ (above_plane_ratio(P2,P3,P0,X)>=tol) ) ||
                (above_plane(P3,P0,P1,P2) ^ (above_plane_ratio(P3,P0,P1,X)>=tol) )
              )
           );
}



/**********************************************************************
 * find_simplex_by_location
 *
 * Given a location in 3-space, build a minimal simplex (tetrahedron) of
 * vertices around it. Attempts to provide a 'nice' simplex to comfortably
 * envelop the provided location.
 *
 * Takes the same parameters as find_vertex_by_location; returns a DUMBLIST
 * containing the VERTEXes that were found.  The DUMBLIST is statically
 * allocated and hence is only temporary storage.
 */

static DUMBLIST *fsbl_cache = 0;

// f_s_calc_stuff: shorthand for some geometrical calculations:
//  - pc gets the point of the argument VERTEX
//  - ac gets the offset vector to that VERTEX from the current origin (x).
//  - acl gets the norm of the offset vector.
#define f_s_calc_stuff(pc,ac,acl,v) do { (pc)=(v->x); diff_3d((ac),(pc),x); (acl)=norm_3d(ac); } while(0)

// f_s_copy_stuff: shorthand for copying a bunch of things.
//   - p gets pc;
//   - a gets ac;
//   - al gets acl;
//   - v gets vc.
#define f_s_copy_stuff(p,a,al,v,pc,ac,acl,vc) do { (p)=(pc); cp_3d((a),(ac)); (al)=(acl); (v)=(vc);} while(0)

DUMBLIST *find_simplex_by_location(POINT3D x, WORLD *w, VERTEX *v, int global) {
    VERTEX *vv;
    static VERTEX *(simplex[6]);                // Gets the VERTEXes of the simplex
    POINT3D a0, a1, a2a, a2b, a3a, a3b;         // These get the vectors relative to the desired point
    static POINT3D origin = {0,0,0};
    NUM *p0, *p1, *p2a, *p2b, *p3a, *p3b;       // These get the x vectors of the VERTEXes we're dealing with
    NUM a0l, a1l, a2al, a2bl, a3al, a3bl;

    // Initialize some things

    simplex[0] = simplex[1] = simplex[2] = simplex[3] = simplex[4] = simplex[5] = 0;

    if(!fsbl_cache) {
    fsbl_cache = new_dumblist();
    }
    dumblist_clear(fsbl_cache);

    // *************************************
    // Find the closest vertex, P0

    vv = find_vertex_by_location(x,w,v,global);
    if(!vv) {
    fprintf(stderr,"find_simplex_by_location: no vertices at all found for location (%g,%g,%g)!\n",x[0],x[1],x[2]);
    return 0;
    }

    // Put the vertex in the simplex, and calculate its ancillary vectors
    simplex[0] = vv;
    f_s_calc_stuff(p0, a0, a0l, vv);
    dumblist_add(fsbl_cache,simplex[0]);

    // End early if the point is exactly atop a vertex
    if(simplex[0]->x[0]==x[0] && simplex[0]->x[1]==x[1] && simplex[0]->x[2]==x[2]){
    //printf("The interpolated point is exactly at one of the vertices: returning with only that one point\n");
    return fsbl_cache;
    }

    // *************************************
    // Find a vertex, P1, close to perpendicular to the x-P0 line, penalizing for distance (?)
    {
    static DUMBLIST *cache = 0;
    NUM costheta, cosphi, l0l1sintheta;
    NUM *pc;
    POINT3D ac,scr;
    NUM acl;
    int i;
    long passno;
    NUM costheta_min, cosphi_min, cosphi_max;
    p1 = 0;

    if( !cache ) {
      cache = new_dumblist();
    }
    dumblist_clear(cache);

    passno = ++(vv->line->fc0->world->passno);
    dumblist_add(cache, vv);
    vv->passno = passno;
    expand_lengthwise(cache, 0, passno);
    expand_via_neighbors(cache, 0, passno);

    costheta_min = -1;
    for(i=1;i<cache->n;i++){ // skip simplex[0]
      vv = ((VERTEX **)(cache->stuff))[i];
      if (V_ISDUMMY(vv))
        continue;
      f_s_calc_stuff(  pc, ac, acl, vv );
      // Find the cosine of the angle between the first VERTEX found and the current one;
      // retain the lowest-cosine (highest angle) VERTEX.
      // An additional factor of vector lengths is provided to penalize distant points
        costheta = fabs(inner_3d(ac,a0) * (acl * a0l));
      if((costheta_min < 0) || (costheta < costheta_min)) { // assignment
          f_s_copy_stuff( p1, a1, a1l, simplex[1],        pc, ac, acl, vv );
          costheta_min = costheta;
      }
    }

    if(!p1) {
      fprintf(stderr, "find_simplex_by_location: FAILED to find a second neighbor! (Should never happen...?\n");
      return fsbl_cache;
    }
    dumblist_add(fsbl_cache,simplex[1]);

    // *************************************
    // Find two vertices, P2a and P2b, close to perpendicular to the x-P0-P1 plane, penalizing for distance (?)
    p2a = 0;
    p2b = 0;

    // Assemble candidates
    dumblist_clear(cache);
    passno = ++(vv->line->fc0->world->passno);
    dumblist_add(cache, simplex[0]);
    dumblist_add(cache, simplex[1]);
    simplex[0]->passno = simplex[1]->passno = passno;
    expand_lengthwise(cache, 0, passno);
    expand_via_neighbors(cache, 0, passno);

    // Find the normal vector to the x-p0-p1 plane
    cross_3d(scr,a0,a1);
    l0l1sintheta = norm_3d(scr);
    cosphi_max = 0;
    cosphi_min = 0;

    // Search through points to find two candidate points for p2
    for(i=2; i<cache->n; i++) { // skip simplex[0] and simplex[1]
      NUM triple;
      vv = ((VERTEX **)(cache->stuff))[i];
      if (V_ISDUMMY(vv))
        continue;

      f_s_calc_stuff(pc, ac, acl, vv);
      triple = inner_3d( ac, scr );                // triple is the triple product (ac . a0 x a1)
      cosphi = triple / l0l1sintheta / acl;   // Divide out to get cos(phi); get absolute value

      if( (cosphi > 0) && (cosphi > cosphi_max) ) {
          f_s_copy_stuff( p2a, a2a, a2al, simplex[2],      pc, ac, acl, vv );
          cosphi_max = cosphi;
      }
      if( (cosphi < 0) && (cosphi < cosphi_min) ) {
          f_s_copy_stuff( p2b, a2b, a2bl, simplex[3],      pc, ac, acl, vv );
          cosphi_min = cosphi;
      }
        // CL - Need to set up some logic if either of these points is not defined
    }
    if(!p2a && !p2b) {
      fprintf(stderr,"find_simplex_by_location: FAILED to find a third neighbor!\n");
      return fsbl_cache;
    }
    dumblist_add(fsbl_cache,simplex[2]);
    dumblist_add(fsbl_cache,simplex[3]);

    // *************************************
    // Find two vertices, P3a and P3b, close to parallel with the vector connecting
    //   x to the centroid of the plane(s) formed by P0-P1-P2a/b. Penalize for distance (?)
    POINT3D c012a, c012b, ax012a, ax012b;
    NUM axl012a, axl012b;
    NUM acosgamma, bcosgamma, acosgamma_max, bcosgamma_max;

    p3a = 0;
    p3b = 0;

    // Assemble candidates
    dumblist_clear(cache);
    passno = ++(vv->line->fc0->world->passno);
    dumblist_add(cache, simplex[0]);
    dumblist_add(cache, simplex[1]);
    dumblist_add(cache, simplex[2]);
    dumblist_add(cache, simplex[3]);
    simplex[0]->passno = simplex[1]->passno = passno;
    if (simplex[2] != 0) simplex[2]->passno = passno;
    if (simplex[3] != 0) simplex[3]->passno = passno;
    expand_lengthwise(cache, 0, passno);
    expand_via_neighbors(cache, 0, passno);
    expand_lengthwise(cache,0,passno);
    simplex[0]->passno = simplex[1]->passno = passno;
    if (simplex[2] != 0) simplex[2]->passno = passno;
    if (simplex[3] != 0) simplex[3]->passno = passno;
    {
      long n = cache->n;
      expand_lengthwise(cache,n,passno);
      expand_via_neighbors(cache,n,passno);
    }

    // Grab the centroids of the plane defined by P0, P1, and P2a/b
    centroid(c012a, a0, a1, a2a);
    centroid(c012b, a0, a1, a2b);

    // Define the vector between these centroids and the point x
    diff_3d(ax012a, x, c012a);
    diff_3d(ax012b, x, c012b);
    axl012a = norm_3d(ax012a);
    axl012b = norm_3d(ax012b);

    // Search through candidate points
    simplex[4] = 0;
    simplex[5] = 0;
    acosgamma_max = 0;
    bcosgamma_max = 0;
    for(i=4; i<cache->n; i++){ // Skipping the first four simplex candidates

        // Define any variables
        int oka, okb;

        // Check the vertex
        vv = ((VERTEX **)(cache->stuff))[i];
        if (V_ISDUMMY(vv))
            continue;

        // Calculate the angle between the candidate point and the {P0,P1,P2} centroid -> x vector
        // Normalize with lengths to penalize distant points
        f_s_calc_stuff( pc, ac, acl, vv);
        acosgamma = fabs(inner_3d(ac, ax012a)) / (acl * axl012a);
        bcosgamma = fabs(inner_3d(ac, ax012b)) / (acl * axl012b);

        // Check that this encloses the point x
        oka = in_simplex( a0, a1, a2a, ac, origin );
        okb = in_simplex( a0, a1, a2b, ac, origin );

        // If the simplex contains x,
        //   and if either the simplex hasn't been filled or the weighted angle exceeds
        //   the current maximum, copy things over.
        if ( oka && (simplex[2] != 0) && ((!simplex[4]) || (acosgamma > acosgamma_max ))) {
            f_s_copy_stuff(p3a, a3a, a3al, simplex[4], pc, ac, acl, vv);
            acosgamma_max = acosgamma;
        }
        if ( okb && (simplex[3] != 0) && ((!simplex[5]) || (bcosgamma > bcosgamma_max ))) {
            f_s_copy_stuff(p3b, a3b, a3bl, simplex[5], pc, ac, acl, vv);
            bcosgamma_max = bcosgamma;
        }
    }
    dumblist_add(fsbl_cache,simplex[4]);
    dumblist_add(fsbl_cache,simplex[5]);
  }

    // *************************************
    // Keeping the simplex points P0 and P1, find the combination of remaining candidates that:
    //  - enclose the point x
    //  ~ minimize simplex volume
    //  - maintain a close distribution of distances from x to simplex faces
    // If no adequate fourth point is found, return the best combination of three
    // Define some lengths
    NUM ed01, ed02a, ed02b, ed12a, ed12b;
    POINT3D c012a, c012b, c013a, c013b, c123a, c123b, c023a, c023b;
    POINT3D ax012a, ax012b, ax013a, ax013b, ax123a, ax123b, ax023a, ax023b;
    NUM axl012a, axl012b, axl013a, axl013b, axl123a, axl123b, axl023a, axl023b;
    NUM avga, avgb, vara, varb;

    // Output all of the simplex points for diagnostics
    // dumblist_clear(fsbl_cache);
    // dumblist_add(fsbl_cache, simplex[0]);
    // dumblist_add(fsbl_cache, simplex[1]);
    // dumblist_add(fsbl_cache, simplex[2]);
    // dumblist_add(fsbl_cache, simplex[3]);
    // dumblist_add(fsbl_cache, simplex[4]);
    // dumblist_add(fsbl_cache, simplex[5]);

    // If both final simplex points are defined...
    if (simplex[4] && simplex[5]){
        // Calculate some distances to face centroids

        // Define some centroids
        centroid(c012a, a0,  a1, a2a);
        centroid(c012b, a0,  a1, a2b);
        centroid(c013a, a0,  a1, a3a);
        centroid(c013b, a0,  a1, a3b);
        centroid(c123a, a1, a2a, a3a);
        centroid(c123b, a1, a2b, a3b);
        centroid(c023a, a0, a2a, a3a);
        centroid(c023b, a0, a2b, a3b);

        // Define the vector between these centroids and the point x
        diff_3d(ax012a, x, c012a);
        diff_3d(ax012b, x, c012b);
        diff_3d(ax013a, x, c013a);
        diff_3d(ax013b, x, c013b);
        diff_3d(ax123a, x, c123a);
        diff_3d(ax123b, x, c123b);
        diff_3d(ax023a, x, c023a);
        diff_3d(ax023b, x, c023b);
        axl012a = norm_3d(ax012a);
        axl012b = norm_3d(ax012b);
        axl013a = norm_3d(ax013a);
        axl013b = norm_3d(ax013b);
        axl123a = norm_3d(ax123a);
        axl123b = norm_3d(ax123b);
        axl023a = norm_3d(ax023a);
        axl023b = norm_3d(ax023b);

        // Calculate variances
        avga = (axl012a + axl013a + axl123a + axl023a) / 4;
        avgb = (axl012b + axl013b + axl123b + axl023b) / 4;
        vara = (pow(axl012a - avga, 2) + pow(axl013a - avga, 2) + pow(axl123a - avga, 2) + pow(axl023a - avga, 2)) / 4;
        varb = (pow(axl012b - avgb, 2) + pow(axl013b - avgb, 2) + pow(axl123b - avgb, 2) + pow(axl023b - avgb, 2)) / 4;

        // Write these best simplex points out
        dumblist_clear(fsbl_cache);
        dumblist_add(fsbl_cache, simplex[0]);
        dumblist_add(fsbl_cache, simplex[1]);
        if (vara < varb) {
            dumblist_add(fsbl_cache, simplex[2]);
            dumblist_add(fsbl_cache, simplex[4]);
        } else{
            dumblist_add(fsbl_cache, simplex[3]);
            dumblist_add(fsbl_cache, simplex[5]);
        }
    // If only the final a-simplex point is found...
    } else if (!simplex[5] && (simplex[4] != 0)){
        dumblist_clear(fsbl_cache);
        dumblist_add(fsbl_cache, simplex[0]);
        dumblist_add(fsbl_cache, simplex[1]);
        dumblist_add(fsbl_cache, simplex[2]);
        dumblist_add(fsbl_cache, simplex[4]);
    // If only the final b-simplex point is found...
    } else if (!simplex[4] && (simplex[5] != 0)){
        dumblist_clear(fsbl_cache);
        dumblist_add(fsbl_cache, simplex[0]);
        dumblist_add(fsbl_cache, simplex[1]);
        dumblist_add(fsbl_cache, simplex[3]);
        dumblist_add(fsbl_cache, simplex[5]);
    // If neither the final a/b-simplex points are found...
    } else if (!simplex[4] && !simplex[5]){
        // Calculate some distances to edges for this case
        ed01  = mpdist(a0,  a1, x);
        ed02a = mpdist(a0, a2a, x);
        ed02b = mpdist(a0, a2b, x);
        ed12a = mpdist(a1, a2a, x);
        ed12b = mpdist(a1, a2b, x);

        // Calculate variances
        avga = (ed01 + ed02a + ed12a) / 3;
        avgb = (ed01 + ed02b + ed12b) / 3;
        vara = (pow(ed01 - avga, 2) + pow(ed02a - avga, 2) + pow(ed12a - avga, 2)) / 3;
        varb = (pow(ed01 - avgb, 2) + pow(ed02b - avgb, 2) + pow(ed12b - avgb, 2)) / 3;

        dumblist_clear(fsbl_cache);
        dumblist_add(fsbl_cache, simplex[0]);
        dumblist_add(fsbl_cache, simplex[1]);
        if (vara < varb) {
            dumblist_add(fsbl_cache, simplex[2]);
        } else{
            dumblist_add(fsbl_cache, simplex[3]);
        }
    }

    return fsbl_cache;
}

/**********************************************************************
 * find_nsimplex_by_location
 *
 * Given a location in 3-space, build a minimal simplex (tetrahedron) of
 * vertices around it.  Works by finding the nearest VERTEX, then attempting
 * to find other nearby vertices.   If the point is outside the cloud of
 * vertices, deliver only a plane.
 *
 * This offshoot of the simplex searcher fetches the neighbor expander
 * to broaden the neighbor tree for an ideal fourth simplex point.
 *
 * Takes the same parameters as find_vertex_by_location; returns a DUMBLIST
 * containing the VERTEXes that were found.  The DUMBLIST is statically
 * allocated and hence is only temporary storage.
 */

// static DUMBLIST *fsbl_cache = 0;

// f_s_calc_stuff: shorthand for some geometrical calculations:
//  - pc gets the point of the argument VERTEX
//  - ac gets the offset vector to that VERTEX from the current origin (x).
//  - acl gets the norm of the offset vector.
#define f_s_calc_stuff(pc,ac,acl,v) do { (pc)=(v->x); diff_3d((ac),(pc),x); (acl)=norm_3d(ac); } while(0)

// f_s_copy_stuff: shorthand for copying a bunch of things.
//   - p gets pc;
//   - a gets ac;
//   - al gets acl;
//   - v gets vc.
#define f_s_copy_stuff(p,a,al,v,pc,ac,acl,vc) do { (p)=(pc); cp_3d((a),(ac)); (al)=(acl); (v)=(vc);} while(0)

DUMBLIST *find_nsimplex_by_location(POINT3D x, WORLD *w, VERTEX *v, int global) {
  VERTEX *vv;
  static VERTEX *(simplex[4]);  // Gets the VERTEXes of the simplex
  POINT3D a0, a1, a2, a3;       // These get the vectors relative to the desired point
  static POINT3D origin = {0,0,0};
  NUM *p0, *p1, *p2, *p3;       // These get the x vectors of the VERTEXes we're dealing with
  NUM a0l, a1l, a2l, a3l;

  simplex[0] = simplex[1] = simplex[2] = simplex[3] = 0;

  if(!fsbl_cache) {
    fsbl_cache = new_dumblist();
  }
  dumblist_clear(fsbl_cache);

  vv = find_vertex_by_location(x,w,v,global);
  if(!vv) {
    fprintf(stderr,"find_simplex_by_location: no vertices at all found for location (%g,%g,%g)!\n",x[0],x[1],x[2]);
    return 0;
  }

  // Put the vertex in the simplex, and calculate its ancillary vectors
  simplex[0] = vv;
  f_s_calc_stuff(p0, a0, a0l, vv);
  dumblist_add(fsbl_cache,simplex[0]);

  //  printf("x is (%g,%g,%g)\n",x[0],x[1],x[2]);
  //  printf("p0 is vertex %ld (%g,%g,%g)\n",simplex[0]->label, simplex[0]->x[0], simplex[0]->x[1], simplex[0]->x[2]);

  if(simplex[0]->x[0]==x[0] && simplex[0]->x[1]==x[1] && simplex[0]->x[2]==x[2]){
    //printf("The interpolated point is exactly at one of the vertices: returning with only that one point\n");
    return fsbl_cache;
  }
  /************
   * Now search for the neighbor of P0 that has the highest angle relative to the P0 - x line
   */
  {
    static DUMBLIST *cache = 0;
    NUM costheta, cosphi, l0l1sintheta;
    NUM *pc;
    POINT3D ac,scr;
    NUM acl;
    int i;
    long passno;
    NUM costheta_min, cosphi_max;
    p1 = 0;

    if( !cache ) {
      cache = new_dumblist();
    }
    dumblist_clear(cache);

    passno = ++(vv->line->fc0->world->passno);
    dumblist_add(cache, vv);
    vv->passno = passno;
    expand_lengthwise(cache, 0, passno);
    expand_via_neighbors(cache, 0, passno);

    costheta_min = -1;
    for(i=1;i<cache->n;i++){ // skip simplex[0]
        vv = ((VERTEX **)(cache->stuff))[i];
      if (V_ISDUMMY(vv))
        continue;
      f_s_calc_stuff(  pc, ac, acl, vv );
      // Find the cosine of the angle between the first VERTEX found and the current one;
      // retain the lowest-cosine (highest angle) VERTEX.
      // An additional factor of vector lengths is provided to penalize distant points
        costheta = fabs(inner_3d(ac,a0) * (acl * a0l * acl * a0l));
      if((costheta_min < 0) || (costheta < costheta_min)) { // assignment
          f_s_copy_stuff( p1, a1, a1l, simplex[1],        pc, ac, acl, vv );
          costheta_min = costheta;
      }
//        printf("%d-%d:%g ",vv->label,simplex[0]->label,costheta );
//            if(costheta<1e-100)
//                printf("\n\t%d:%g,%g,%g;  %fd:%g,%g,%g\n", vv->label,ac[0],ac[1],ac[2],simplex[0]->label,a0[0],a0[1],a0[2]);
//        fflush(stdout);
    }

    // printf("p1 is vertex %ld (%g,%g,%g)\n",simplex[1]->label, simplex[1]->x[0], simplex[1]->x[1], simplex[1]->x[2]);

    if(!p1) {
      fprintf(stderr, "find_simplex_by_location: FAILED to find a second neighbor! (Should never happen...?\n");
      return fsbl_cache;
    }
    dumblist_add(fsbl_cache,simplex[1]);

    /* Now we have two points (p0 and p1) and their ancillary info.  Search the neighbors of
     * p0 and p1 to find the maximum angle out of the (x,p0,p1) plane.  We do that by trying
     * to find a vector that is (parallel | antiparallel) to the (x,p0,p1) normal vector.
     *
     * That means maximizing abs(cos(phi)), the angle between the vector to the trial point and
     * the normal fo (x,p0,p1).
     */

    p2 = 0;

    // Assemble candidates
    dumblist_clear(cache);
    passno = ++(vv->line->fc0->world->passno);
    dumblist_add(cache, simplex[0]);
    dumblist_add(cache, simplex[1]);
    simplex[0]->passno = simplex[1]->passno = passno;
    expand_lengthwise(cache, 0, passno);
    expand_via_neighbors(cache, 0, passno);


    // scr gets the normal vector to the x,p0,p1 plane. l0l1sintheta gets its length.
    cross_3d(scr,a0,a1);
    l0l1sintheta = norm_3d(scr);
    cosphi_max = -0.1; /* impossibly tiny since we take absolute value later */

    for(i=2; i<cache->n; i++) { // skip simplex[0] and simplex[1]
      NUM triple;
      vv = ((VERTEX **)(cache->stuff))[i];
      if (V_ISDUMMY(vv))
        continue;

      f_s_calc_stuff(pc, ac, acl, vv);
      triple = inner_3d( ac, scr );                // triple is the triple product (ac . a0 x a1)
      cosphi = fabs(triple / l0l1sintheta / acl / l0l1sintheta / acl);   // Divide out to get cos(phi); get absolute value

      if( cosphi > cosphi_max ) {
          f_s_copy_stuff( p2, a2, a2l, simplex[2],      pc, ac, acl, vv );
          cosphi_max = cosphi;
      }
    }
    if(!p2) {
      fprintf(stderr,"find_simplex_by_location: FAILED to find a third neighbor!\n");
      return fsbl_cache;
    }
    dumblist_add(fsbl_cache,simplex[2]);

    // printf("p2 is vertex %ld (%g,%g,%g)\n",simplex[2]->label, simplex[2]->x[0], simplex[2]->x[1], simplex[2]->x[2]);

    /* We've got three points that are hopefully noncoplanar with the original point.
     * Now try to find a fourth VERTEX that makes an enclosing simplex.
     */


    // A New Hope... find the vertex close to the vector defined by the normal
    //   of the plane P0, P1, P2 and the point x.

    POINT3D c012, ax012;
    NUM axl012;
    NUM cosgamma, cosgamma_max;
    p3 = 0;
    long n;
    n = 0;

    // Assemble candidates
    dumblist_clear(cache);
    passno = ++(vv->line->fc0->world->passno);
    dumblist_add(cache, simplex[0]);
    dumblist_add(cache, simplex[1]);
    dumblist_add(cache, simplex[2]);
    simplex[0]->passno = simplex[1]->passno = simplex[2]->passno = passno;
    expand_lengthwise(cache, 0, passno);
    expand_via_neighbors(cache, 0, passno);
    expand_lengthwise(cache,0,passno); // Additional neighbor search(?) from original version
    simplex[0]->passno = simplex[1]->passno = simplex[2]->passno = passno;
    //printf("\n Neighbor search - n:%ld passno:%ld", n, passno);
    {
      long n = cache->n;
      expand_lengthwise(cache,n,passno);
      expand_via_neighbors(cache,n,passno);
    }

    // Grab the centroid of the plane defined by P0, P1, and P2
    centroid(c012, a0, a1, a2);

    // Define the vector between this centroid and the point x
    diff_3d(ax012, x, c012);
    axl012 = norm_3d(ax012);

    // Search through candidate points
    simplex[3] = 0;
    cosgamma_max = 0;
    for(i=3; i<cache->n; i++){ // Skipping the first three simplexes

        // Define any variables
        int ok;
        int opp;
        //int plnchk;

        // Check the vertex
        vv = ((VERTEX **)(cache->stuff))[i];
        if (V_ISDUMMY(vv))
            continue;

        // Calculate the angle between the candidate point and the {P0,P1,P2} centroid -> x vector
        // Normalize with lengths to penalize distant points
        f_s_calc_stuff( pc, ac, acl, vv);
        cosgamma = fabs(inner_3d(ac, ax012)) / (acl * acl * axl012 * axl012);

        // To force some out of plane movement, could we utilize above_plane here?
        // Perhaps checking that P3 lies above the plane formed by a combination of
        //   x and two of P0, P1, P2
        // Although the point might be *barely* above one of these planes...
        //plnchk = above_plane(a0, a1, x, ac);

        // Check that the distance of x from the plane centroid is within a tolerance ratio
        //   of the distance from the centroid to the remaining simplex point. Rinse and repeat.

        // Check that this encloses the point x
        ok = in_simplex( a0, a1, a2, ac, origin );

        // Check that this fourth point is opposite of P2
        opp = opposite_plane(a0, a1, origin, ac, a2);

        // CL - Print out fluxon ID
        //printf("\n okay: %d , acl: %g , axl012: %g, cosgamma: %g, flxn: %ld", ok, acl, axl012, cosgamma, vv->line->label);

        //if (vv->line->label == 300){
        //    printf("\n \t $a0 = pdl(%g, %g, %g)", a0[0]+x[0], a0[1]+x[1], a0[2]+x[2]);
        //    printf("\n \t $a1 = pdl(%g, %g, %g)", a1[0]+x[0], a1[1]+x[1], a1[2]+x[2]);
        //    printf("\n \t $a2 = pdl(%g, %g, %g)", a2[0]+x[0], a2[1]+x[1], a2[2]+x[2]);
        //    printf("\n \t $ac = pdl(%g, %g, %g)", ac[0]+x[0], ac[1]+x[1], ac[2]+x[2]);
        //    printf("\n \t $as = pdl($ac, $a0, $ac, $a1, $ac, $a2)");
        //    printf("\n \t $win->replot({trid=>1}, {with=>'lines'},$as->using(0,1,2))");
        //    printf("\n");
        //}

        // If the simplex contains x,
        //   and if either the simplex hasn't been filled or the weighted angle exceeds
        //   the current maximum, copy things over.
        if ( ok && opp && ((!simplex[3]) || (cosgamma > cosgamma_max ))) {
        //if ( ok && (cosgamma > 0) && ((!simplex[3]) || (cosgamma > cosgamma_max ))) {
            f_s_copy_stuff(p3, a3, a3l, simplex[3], pc, ac, acl, vv);
            cosgamma_max = cosgamma;
            //printf("\n Good point : okay: %d , acl: %g , axl012: %g, cosgamma: %g, flxn: %ld", ok, acl, axl012, cosgamma, vv->line->label);
        }
    }

    // Expand the neighbor search if no fourth point is immediately found
    while ((!simplex[3]) && (n < 10) && (0)){
        // Declare any variables
        int ok;
        int opp;
        int i;

        //CL - Assemble a few more candidates.
        passno = ++(vv->line->fc0->world->passno);
        n = ++n;
        i = n;
        simplex[0]->passno = simplex[1]->passno = simplex[2]->passno = passno;
        expand_lengthwise(cache, n, passno);
        expand_via_neighbors(cache, n, passno);
        expand_lengthwise(cache,n,passno); // Additional neighbor search(?) from original version
        simplex[0]->passno = simplex[1]->passno = simplex[2]->passno = passno;
        //printf("\n Neighbor search - n:%ld passno:%ld", n, passno);

        // Check the vertex
        vv = ((VERTEX **)(cache->stuff))[i];
        if (V_ISDUMMY(vv))
            continue;

        // Calculate the angle between the candidate point and the {P0,P1,P2} centroid -> x vector
        // Normalize with the lengths to penalize distant points
        f_s_calc_stuff( pc, ac, acl, vv);
        cosgamma = fabs(inner_3d(ac, ax012)) / (acl * acl * axl012 * axl012);

        // Check that this encloses the point x
        ok = in_simplex( a0, a1, a2, ac, origin );

        // Check that this fourth point is opposite of P2
        opp = opposite_plane(a0, a1, origin, ac, a2);

        // If the simplex contains x,
        //   and if either the simplex hasn't been filled or the weighted angle exceeds
        //   the current maximum, copy things over.
        if ( ok && opp && ((!simplex[3]) || (cosgamma > cosgamma_max ))) {
            f_s_copy_stuff(p3, a3, a3l, simplex[3], pc, ac, acl, vv);
            cosgamma_max = cosgamma;
            //printf("\n Good point : okay: %d , acl: %g , axl012: %g, cosgamma: %g, flxn: %ld", ok, acl, axl012, cosgamma, vv->line->label);
        }
    }

    //// Original simplex volume minimization code
    //passno = ++(vv->line->fc0->world->passno);
    //dumblist_clear(cache);
    //dumblist_add(cache,simplex[0]);
    //dumblist_add(cache,simplex[1]);
    //dumblist_add(cache,simplex[2]);

    //// printf("passno is %ld...",passno);
    //simplex[0]->passno = simplex[1]->passno = simplex[2]->passno = passno;
    //expand_lengthwise(cache,0,passno);
    //expand_via_neighbors(cache,0,passno);
    //expand_lengthwise(cache,0,passno);
    //{
    //  long n = cache->n;
    //  expand_lengthwise(cache,n,passno);
    //  expand_via_neighbors(cache,n,passno);
    //}
    ////    printf("Testing %d candidates for the final point of the simplex...\n",cache->n - 3);

    //simplex[3] = 0;
    //for(i=3; i<cache->n; i++) {
    //  int ok;
    //  vv =  ((VERTEX **)(cache->stuff)) [i] ;
    //  if (V_ISDUMMY(vv))
    //    continue;

    //  // acl gets the volume of the simplex formed by the three prior points and the candidate.
    //  // We want the smallest simplex that encloses the sample point.
    //  f_s_calc_stuff(pc, ac, acl, vv);

    //  // printf(" v%ld ",vv->label);
    //  ok =  in_simplex( a0, a1, a2, ac, origin );
    //
    //  // Some diagnostics for funky simplex finding
    //  //if( !ok && x[0] > 0 && x[0] < 0.5 && x[1] > 0 && x[1] < 0.5) {
    //    //printf("\n \t %ld", vv->line->label);
    //    //printf("\n \t $a0 = pdl(%g, %g, %g)", a0[0], a0[1], a0[2]);
    //    //printf("\n \t $a1 = pdl(%g, %g, %g)", a1[0], a1[1], a1[2]);
    //    //printf("\n \t $a2 = pdl(%g, %g, %g)", a2[0], a2[1], a2[2]);
    //    //printf("\n \t $ac = pdl(%g, %g, %g)", ac[0], ac[1], ac[2]);
    //    //printf("\n \t $as = pdl($a0, $a1, $a2, $ac)");
    //    //printf("\n \t $win->plot({trid=>1}, {with=>'points'},$as->using(0,1,2))");
    //    //printf("\n");
    //    //fflush(stdout);
    //  //}
    //
    //  if( ok ) {
    //    //printf("%d ",ok);
    //    //printf("\n \t %ld", vv->line->label);
    //    //printf("\n \t $a0 = pdl(%g, %g, %g)", a0[0], a0[1], a0[2]);
    //    //printf("\n \t $a1 = pdl(%g, %g, %g)", a1[0], a1[1], a1[2]);
    //    //printf("\n \t $a2 = pdl(%g, %g, %g)", a2[0], a2[1], a2[2]);
    //    //printf("\n \t $ac = pdl(%g, %g, %g)", ac[0], ac[1], ac[2]);
    //    //printf("\n \t $as = pdl($a0, $a1, $a2, $ac)");
    //    //printf("\n \t $win->plot({trid=>1}, {with=>'points'},$as->using(0,1,2))");
    //    //printf("\n");
    //    //fflush(stdout);
    //    // If the simplex hasn't been filled yet (always true on the first OK) or if
    //    // we're better than the last simplex-filler, copy the fourth point to the simplex.
    //    if(  (!simplex[3]) || (acl < a3l )  )
    //      f_s_copy_stuff(p3, a3, a3l, simplex[3],          pc, ac, acl, vv);
    //  }
    //}
    //  printf("\n");

    if(!simplex[3]) {
      // Not finding a 4th neighbor is not unusual if you are outside the sim, so we don't throw an error.
      // Kept as a comment because it is useful for debugging.
      //      fprintf(stderr,"find_simplex_by_location: FAILED to find a fourth neighbor!\n");
      //      fprintf(stderr,"\tChecking planarity condition...\n");
      //      fprintf(stderr,"\tx is (%g,%g,%g)\n\tp0 is (%g,%g,%g) (v%ld)\n\tp1 is (%g,%g,%g) (v%ld)\n\tp2 is (%g,%g,%g) (v%ld)\n",x[0],x[1],x[2],
      //              simplex[0]->x[0],simplex[0]->x[1],simplex[0]->x[2],simplex[0]->label,
      //              simplex[1]->x[0],simplex[1]->x[1],simplex[1]->x[2],simplex[1]->label,
      //              simplex[2]->x[0],simplex[2]->x[1],simplex[2]->x[2],simplex[2]->label
      //              );

      return fsbl_cache;
    }
    //printf("p3 is vertex %ld (%g,%g,%g)\n",simplex[3]->label, simplex[3]->x[0], simplex[3]->x[1], simplex[3]->x[2]);

    dumblist_add(fsbl_cache,simplex[3]);

    /* Should we try to shrink the simplex here? */

  }
  return fsbl_cache;
}


/**********************************************************************
 * interpolate_lin_3d - given a set of 1 <= n <= 4 points (p) in
 * 3-space at which values (val) are known, and a location (x) to
 * which to interpolate, do so.
 *
 * Scrozzles the contents of the original arrays.
 */
NUM interpolate_lin_3d(POINT3D x, NUM p[4*3], NUM val[4], int n) {
  POINT3D xx;
  POINT3D a, aa;
  PLANE pl;
  NUM alpha;
  NUM wgt_acc = 1.0;
  NUM acc = 0;
  int i;

  /*  printf("interpolate_lin_3d: number of simplex points is %d, x is %f, %f, %f\n",n,x[0],x[1],x[2]);
  printf("interpolate_lin_3d: val is %f, %f, %f, %f\n",val[0],val[1],val[2],val[3]);
  printf("p (simplex points):\n");
  for (i=0;i<3*n;i++){
    if(i%3==0)
      printf("(");
    printf("%f, ",p[i]);
    if(i%3==2)
      printf(") val = %f\n",val[i/3]);
  }
  */
  switch(n) {
  case 4:
    // Full simplex: Find intersection of the line P3 -> X with the P0,P1,P2 plane,
    // and reduce to the triangular (planar) case.
    points2plane(&pl, &p[3*0], &p[3*1], &p[3*2]);
    //printf("PLANE: origin (%f, %f, %f), normal (%f, %f, %f)\n",pl.origin[0],pl.origin[1],pl.origin[2],pl.normal[0],pl.normal[1],pl.normal[2]);
    p_l_intersection(a, &pl, x, &p[3*3]);
    //printf("OUT: a= (%f, %f, %f)\n",a[0],a[1],a[2]);

    // Find relative length
    diff_3d(xx, x, &p[3*3]); //xx gets the vector difference between x and the 4th simplex point.
    diff_3d(aa, a, &p[3*3]); //aa gets the vector difference between a and the 4th simplex point.
    alpha = inner_3d(xx,aa) / norm2_3d(aa);
    if(isfinite(alpha)){
      acc += wgt_acc * (1-alpha) * val[3];
      wgt_acc *= alpha;
      /*printf("a= %f, %f, %f\n",a[0],a[1],a[2]);
      printf("xx= %f, %f, %f\n",xx[0],xx[1],xx[2]);
      printf("aa= %f, %f, %f\n",aa[0],aa[1],aa[2]);
      printf("n=4; alpha=%f,wgt_acc=%f,acc=%f\n",alpha,wgt_acc,acc);*/
      cp_3d(x,a); // Move X into the plane and fall through to planar case
    }
  case 3:
    // Triangle (plane): Find the intersection of the line from P2->X with the
    // P0,P1 "line" (project everything down to the P0,P1,P2 plane), and reduce
    // to the linear case.
    lp2plane(&pl, &p[3*0], &p[3*1], &p[3*2]);
    p_l_intersection(a, &pl, x, &p[2*3]);

    // find relative length
    diff_3d(xx, x, &p[2*3]);
    diff_3d(aa, a, &p[2*3]);
    alpha = inner_3d(xx,aa) / norm2_3d(aa);
    if(isfinite(alpha)){
      acc += wgt_acc * (1-alpha) * val[2];
      wgt_acc *= alpha;
      //printf("n=3; alpha=%f,wgt_acc=%f,acc=%f\n",alpha,wgt_acc,acc);
      cp_3d(x,a); // Move X onto the P0-P1 line and fall through to the linear case
    }
  case 2:
    // Line: direct interpolation
    diff_3d(xx, x, &p[3*1]);
    diff_3d(a,  &p[3*0], &p[3*1]);

    alpha = inner_3d( xx, a ) / norm2_3d(a);
    if(isfinite(alpha)){
      acc += wgt_acc * (1 - alpha) * val[1];
      wgt_acc *= alpha;
      //printf("n=2; alpha=%f,wgt_acc=%f,acc=%f\n",alpha,wgt_acc,acc);
    }
  case 1:

    acc += wgt_acc * val[0];
    //printf("n=1; alpha=%f,wgt_acc=%f,acc=%f\n",alpha,wgt_acc,acc);
    break;
  }
  return acc;
}


/**********************************************************************
 * interpolate_value_simplex
 *
 * Given a location and a simplex of VERTEXes that surround it,
 * interpolate the value of a NUM between the vertices.  Any NUM
 * stored in a VERTEX can be interpolated.
 *
 */

NUM interpolate_value_simplex( POINT3D x, DUMBLIST *dl, int val_offset) {
  POINT3D scr;
  NUM acc;
  NUM wgt_acc;
  NUM wgt;
  int i;
  NUM p[4*3];
  NUM val[4];
  acc = 0;
  wgt_acc = 0;


  for(i=0;i<dl->n;i++) {
    cp_3d( &p[i*3], ((VERTEX **)(dl->stuff))[i]->x );
    val[i] = * ( (NUM *)( (void *)(dl->stuff[i]) + val_offset ) );
  }
  cp_3d(scr, x);
  return interpolate_lin_3d(scr, p, val, dl->n);
}

/**********************************************************************
 * interpolate_value
 * same as interpolate_value_simplex, but you don't supply the simplex --
 * find_simplex_by_location gets called for you.
 *
 * The calling parameters are the same as for find_simplex_by_location.
 */

NUM interpolate_value( POINT3D x, WORLD *w, VERTEX *v, int global, int val_offset ) {
  DUMBLIST *dl;
  int i;
  dl = find_simplex_by_location(x, w, v, global);
  printf("found %d vertices: ",dl->n);
  for(i=0;i<dl->n;i++) {
    printf(" %ld ",((VERTEX **)dl->stuff)[i]->label);
  }
  printf("\n");

  return interpolate_value_simplex(x, dl, val_offset);
}



















/**********************************************************************
 * project_n_fill_photosphere
 *
 * Similar to project_n_fill, but instead of using the fluxel to
 * derive the plane, we use photosphere1. Also, there are no
 * non-linear projections, and all distances are from vertices instead
 * of fluxels.
 *
 * Given a VERTEX and a DUMBLIST of VERTICES, project the DUMBLIST
 * down to 2-D using fl_segment_deluxe_dist to find the projected radius,
 * and fill the angle and radius fields of the members of the DUMBLIST.
 *
 * Also accumulate the closest-real-approach distance between neighbors
 * and the current VERTEX.  By examining only neighbors that have been
 * winnowed with the special metric, we risk missing the actual
 * minimum distance -- but not very often.  Hopefully it won't lead to
 * pathological evolution problems!
 *
 */
void project_n_fill_photosphere(VERTEX *v, DUMBLIST *horde) {
  int i;
  NUM pm[9];
  NUM X0[3];
  //NUM p0[3], p1[3]; //points of closest approach (essentially v and v1)
  NUM or[3]={0,0,0};//how do i initialize this array? (syntax)
  NUM len, r;
  VERTEX *v1;
  char crunch = 0;

  if(v->line->fc0->world->verbosity >= 5)
    printf("Entered project_n_fill_phototsphere... got a horde of %d candidates...\n",horde->n);

  if(v->line->fc0->world->verbosity >= 5)
    printf("calling projmatrix with (%g,%g,%g) and (%g,%g,%g)\n",or[0],or[1],or[2],v->line->fc0->world->photosphere.plane->normal[0],v->line->fc0->world->photosphere.plane->normal[1],v->line->fc0->world->photosphere.plane->normal[2]);

  projmatrix(pm,or,v->line->fc0->world->photosphere.plane->normal); //pm is the projection matrix

  if(v->line->fc0->world->verbosity >= 5)
    printf("back\n");

  v->r_cl = -1; /* Initialize closest-approach accumulator */

  for(i=0;i<horde->n;i++) {

    if(v->line->fc0->world->verbosity >= 5) {
      printf("i=%d ",i);
      fflush(stdout);
    }


    v1 = (VERTEX *)((horde->stuff)[i]);

    //p0 = v->x;

    //p1 = v1->x;

    diff_3d(X0, v->x, v1->x); /*X0 is the line of closest approach
                                between v and neighbor v1*/

    r = norm_3d(X0);

    if(r<=0 || !isfinite(r)) { /*not sure what these circumstance would be*/
      horde->stuff[i] = 0;
      crunch=1;
    }
    else {

      if(r<v->r_cl || v->r_cl < 0 ) /* Accumulate closest approach distance */
        v->r_cl = r;

      mat_vmult_3d(v1->scr,pm,X0);
      //v1->scr is X0 rotated to perpendicular to the photosphere.

      len = norm_2d(v1->scr);      /* 2-D length of vector, should be
                                    the same as r because length is
                                    preserved for rotation
                                    matrices*/

      if (len != r){
        printf("Oops! len!=r in project_n_fill_photsphere, v %ld, neighbor %ld,len=%f,r=%f \n",v->label, v1->label,len,r);
        printf("pm is ((%f,%f,%f),(%f,%f,%f),(%f,%f,%f)) \n",pm[0],pm[1],pm[2],pm[3],pm[4],pm[5],pm[6],pm[7],pm[8]);

      } /*else {
        printf("len=r for neighbor %d \n",v1->label);
        }*/

      if(v->line->fc0->world->verbosity >= 5)
        printf("len=%g for vertex %ld\n",len,v1->label);

      if(len <= 0) {

        if(v->line->fc0->world->verbosity == 4)
            printf("len=%g for vertex %ld\n",len,v1->label);

        horde->stuff[i] = 0;
        crunch=1;
      }
      else {
        //scale_3d(v1->scr,v1->scr,r/len); /* not this b/c should be the same.*/

        v1->a = ATAN2(v1->scr[1],v1->scr[0]); //accumulate angles
        v1->r = r; //accumulate non-linear projected distance.
      }
    }

  }

  /* If some of the items were invalidated, crunch down the horde... */
  if(crunch) {
    int j;
    for(i=j=0;i<horde->n;i++) {
      for(;i<horde->n && !horde->stuff[i]; i++)
        ;
      if(i<horde->n)
        horde->stuff[j++] = horde->stuff[i];
    }
    horde->n = j;
  }
}


/**********************************************************************
 * The only difference for the photosphere version is that it doesn't
 * skip end vertices.
 */
void hull_2d_us_photosphere(HULL_VERTEX *hull, DUMBLIST *horde, VERTEX *central_v) {
  int n = horde->n;
  int horde_i;
  long passno;
  WORLD *w = 0;
  int i,j;

  if(ws==0) {
    ws = new_dumblist();
    dumblist_grow(ws,16);
  }

  if(! horde->n)
    return;

  for(i=0;!w && i<horde->n; i++) //go through each potential neighbor
    if( (((VERTEX **)(horde->stuff))[i])->line )
      w = (((VERTEX **)(horde->stuff))[i])->line->fc0->world;
  if(!w) {
    fprintf(stderr,"You're in trouble, guv! exiting (horde->n was %d)...\n",horde->n);
    exit(99);
  } //somehow your potential neighbor was not connected to a fluxon

  // use passno to prevent multiple looks at the same candidate.
  passno = ++w->passno;

  /* Set our output list to zero length */
  ws->n = 0;

  /* Nothing in the horde?  We're done! */
  if(horde->n == 0)
    return;


  /* always add the first eligible VERTEX in the list... */
  for(i=0;i<horde->n &&
        ( horde->stuff[i] == central_v ||              // skip the main vertex if present
            plasmoid_conjugate( (VERTEX *)(horde->stuff[i]), central_v )
        );
      i++) //find first candidate
    ;
  if(i<horde->n) { //first look at this candidate and add it.
    ((VERTEX **)(horde->stuff))[i]->passno = passno;
    ws->stuff[0] = horde->stuff[i];
    ws->n = 1;
    hull[0].open = 1;
    perp_bisector_2d(hull[0].bisector, 0, ((VERTEX *)(horde->stuff[i]))->scr);
  } else {
    horde->n = 0;
    return;
  }//will we ever remove it?


  /** Loop over the horde and check vertices as we go... **/
  for(horde_i=1; horde_i<horde->n; horde_i++) {
    VERTEX *v =  ((VERTEX **)(horde->stuff))[horde_i];
    VERTEX *nv, *pv;
    static HULL_VERTEX scrhv;
    HULL_VERTEX *nh, *ph, *scr;
    int next_idx;
    char keep;
    int  flag;

    // Check for skip conditions...
    if(v->passno == passno ||                   // Already been here
       v == central_v ||                        // Skip ourselves if we encounter us
       ( v->line==central_v->line &&            // Several end conditions for same fluxon:
        ( v->next==central_v ||                // Skip next segment
            v->prev==central_v ||                // Skip previous segment
            plasmoid_conjugate( v, central_v)
        )
       )
       ) {
      // Skip (do nothing; just move on to the next candidate)
    } else {
      v->passno = passno; // Mark it processed (avoid duplicated effort)

      keep = 0;

      /* Linear searches are slow, but not so bad -- there are
       generally under 10 elements in ws, so binary search probably
       isn't worth the effort.  */
      for(next_idx = 0;
            next_idx < ws->n &&
            (nv = (((VERTEX **)(ws->stuff))[next_idx]) )->a < v->a;  // assign to nv
          next_idx++
          )
        ;//not sure about this loop


      /******************************
       * Now set nv to the next vertex after this one and pv to the previous one.
       * next_idx is the true insertion point if this element is kept.
       */

      if( nv->a < v->a ) {
        /* We went off the end... */

        next_idx = ws->n;
        nv = ((VERTEX **)(ws->stuff))[0];
        nh = hull;
        pv = ((VERTEX **)(ws->stuff))[ws->n-1];
        ph = hull + ws->n-1;
      } else {
        /* We didn't... */
        if(next_idx==0) {
            nh = hull;
            pv = ((VERTEX **)(ws->stuff))[ws->n-1];
            ph = hull + ws->n-1;
        } else {
            nh = hull + next_idx;
            pv = ((VERTEX **)(ws->stuff))[next_idx - 1];
            ph = hull + next_idx-1;
        }
      }

      flag = check_hullpoint(v, pv,ph, nv,nh, &scrhv);

      if(flag) {
        /* Make sure we have room in the workspace */
        if(ws->size <= ws->n)
         dumblist_grow(ws, ws->size * 1.5 + 10);

        /* Make room in the workspace dumb list and the hull list */
        for( j = ws->n; j>next_idx; j--){
            ws->stuff[j] = ws->stuff[j-1];
            hv_cp(hull+j, hull+j-1);
        }

        ws->n++;

        /* Set the insertion-point values */
        ws->stuff[next_idx] = v;
        hv_cp(hull+next_idx, &(scrhv));

        hull[next_idx].open = ( flag>0 && (flag & 1) ); /* logical && and bitwise & */
        hull[   ( next_idx ? next_idx : ws->n ) - 1  ].open =
            ( flag>0 && (flag & 2) ); /* logical && and bitwise & */

        /**************************
        * Now that we've put the new guy in, make sure he doesn't invalidate
        * his neighbors on the prev side...
        */
        {
          int n, p, pp, nn, nnn, check_flag;

          n = next_idx;
          p = n; MOD_DEC(p, ws->n);
          pp= p; MOD_DEC(pp, ws->n);
          nn = n; MOD_INC(nn, ws->n);
          nnn=nn; MOD_INC(nnn,ws->n);

          /* Walk backward until we find a still-keep-worthy vertex */
          while( p != n &&
            ! (check_flag = check_hullpoint( (VERTEX *)(ws->stuff[p]),
                                                            (VERTEX *)(ws->stuff[pp]), hull + pp,
                                                            (VERTEX *)(ws->stuff[n]),  hull + n,
                                                            hull + p))
          ) {

            p=pp;
            MOD_DEC(pp, ws->n);
          }

          hull[p].open = (check_flag > 0) && (check_flag & 1);
          hull[pp].open = (check_flag > 0) && (check_flag & 2);

          /* walk forward until we find a still-keep-worthy vertex */
          while( n != p &&
                 ! (check_flag = check_hullpoint( (VERTEX *)(ws->stuff[nn]),
                                                  (VERTEX *)(ws->stuff[n]),  hull + n,
                                                  (VERTEX *)(ws->stuff[nnn]), hull + nnn,
                                                  hull+nn
                                                  ) )
                 ) {

            nn=nnn;
            MOD_INC(nnn, ws->n);
          }

          hull[nn].open = (check_flag > 0) && (check_flag & 1);
          hull[n].open = (check_flag > 0) && (check_flag & 2);

          /* Now p is the first keepworthy vertex, walking backward from n, and */
          /* nn is the first keepworthy vertex, walking forward from n.         */
          /* Crunch down the list... */
          {

            int ii,jj, diff;

            if( p <= n-1 ) {
              if( nn <= n-1 ) {

                // ...*****...|.....  case
                for(ii=0, jj=nn; jj<=p; jj++,ii++) {
                  hv_cp( hull+ii, hull+jj );
                  ws->stuff[ii] = ws->stuff[jj];
                }
                hv_cp( hull+ii, hull +n );
                ws->stuff[ii] = ws->stuff[n];
                ii++;

              } else {

                // *****...|...**** case

                ii=p+1;
                if(ii == n) {
                  ii++;
                } else {
                  hv_cp( hull + ii, hull + n );
                  ws->stuff[ii] = ws->stuff[n];
                  ii++;
                }
                if(ii<nn) {
                  for( jj=nn; jj < ws->n; ii++,jj++) {
                    hv_cp(hull + ii, hull + jj);
                    ws->stuff[ii] = ws->stuff[jj];
                  }
                } else if(ii==nn) {
                  ii = ws->n;
                }
              }
            } else {

              // ....|...****.. case
              if(n>0) {
                hv_cp(hull, hull+n);
                ws->stuff[0] = ws->stuff[n];
              }
              ii=1;

              if( nn>1 ) {
                for( jj=nn; jj<=p; ii++, jj++) {
                  hv_cp( hull + ii, hull + jj );
                  ws->stuff[ii] = ws->stuff[jj];
                }
              } else {
                ii= p + 1;
              }

            }/* end of case testing */
            ws->n = ii;

          } /* end of crunching convenience block */
        } /* end of neighbor-testing convenience block */
      } else {
        /* no keeping -- just ignore this point and start the next loop */
      }
    } /* end of non-duplication test block */
  } /* end of horde loop */

  /* Calculate the (more expensive) true-atan angles for the hull points and for the
   * actual neighbors, since there are so few of them.
   */

  for(  i=0; i<ws->n; i++ ) {
    VERTEX *vv = ((VERTEX **)(ws->stuff))[i];

    if(hull[i].open) {
      VERTEX *vn = ((VERTEX **)(ws->stuff))[MOD_NEXT(i,ws->n)];
      hull[i].a_r = atan2(vn->scr[1],vn->scr[0]) - PI/2;
      hull[i].a_l = atan2(vv->scr[1],vv->scr[0]) + PI/2;
    } else {
      hull[i].a_r = hull[i].a_l = atan2(hull[i].p[1],hull[i].p[0]);
    }

    vv->a = atan2(vv->scr[1],vv->scr[0]);
    horde->stuff[i] = vv;
  }
  horde->n = i;
}
