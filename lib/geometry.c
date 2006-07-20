#define DEBUG_HULL 0
/**********************************************************************
 * Geometry.c -- routines that embody geometrical operations
 * used for the fieldine model.
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
 * This is geometry.c 1.1 - part of the FLUX 1.1 release. 
 *
 */
#include <math.h>
#include <stdio.h>
#include "data.h"
#include "geometry.h"


#ifndef NAN
#define NAN (nan("NAN"))
#endif


char *code_info_geometry="%%%FILE%%% (new unsorted hull routine)";

/**********************************************************************
 **********************************************************************
 *****  Vector basics
 *****  norm, inner product, cross product, scalar product, sum, 
 *****  difference, copy.
 
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
void *cross(NUM *out, NUM *p0, NUM *p1) {
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
  
/**********************************************************************
 * sum - Add two 3-vectors and stick 'em into the destination array.
 * OK to have the destination be one of the sources.
 */
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
 * returns the alpha coefficient, which can often be safely ignored.
 */

NUM p_ls_closest_approach(NUM p0[3], NUM a0[3], NUM b0[3], NUM c0[3]) {
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
  if(!finite(dist))
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
 */
NUM fl_segment_dist(VERTEX *v0, VERTEX *v1) {
  NUM P0[3],P1[3];
  return fl_segment_deluxe_dist(P0,P1,v0,v1);
}

NUM fl_segment_deluxe_dist(NUM P0[3],NUM P1[3], VERTEX *v0, VERTEX *v1) {
  NUM a,b,Plen;
  NUM P[3], seg[3], cr[3];

  /* Exclude trivial cases, adjacent-segment case, and next-nearest-segment case. */
  /* (farther segments aren't as likely to cause trouble) */
  if(!v0 || !v1 || !v0->next || !v1->next)
    return 1e50;

  if(v0==v1 ||   (v0->next == v1) ||  (v1->next == v0) )
    return 1e50;

  if(v0->line->fc0->world->verbosity >= 5) 
    printf("fl_segment_dist: calling ls_closest_approach...\n");

  ls_closest_approach(P0,P1,v0->x,v0->next->x,v1->x,v1->next->x);

  if(v0->line->fc0->world->verbosity >= 5) 
    printf("fl_segment_dist: got back... P1=%d; P0=%d\n",P1,P0);

  diff_3d(P,P1,P0);
  Plen = norm_3d(P);
  if(Plen==0)
    return 1e50;

  if (v0->line->fc0->world->verbosity >= 6) {
    printf("\nfl_segment_deluxe_dist: v0=%d,v1=%d; Plen is %g\n",v0->label,v1->label,Plen);
    printf("P0=(%g,%g,%g); P1=(%g,%g,%g)\n",P0[0],P0[1],P0[2],P1[0],P1[1],P1[2]);
  }

  // inverse FOURTH POWER sine!
  /* Scale by the inverse square sine of the projection angle:  */
  /* AB x (closest).  Note that this makes the metric asymmetric! */
   
  scale_3d(P,P,1.0/Plen);

  diff_3d(seg,v0->next->x,v0->x);
  scale_3d(seg,seg,1.0/norm_3d(seg));

  cross(cr,P,seg);
  a = norm2_3d(cr);
  a *= a;               // sin^2 --> sin^4

  if(a == 0) 
    return 1e50;

  Plen /= a;

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
    return !finite(out[2]);
  }    

  /* P-not-the-origin case */
  out[0] = ( Q[0] + P[0] ) / 2;
  out[1] = ( Q[1] + P[1] ) / 2;
  
  /* Slope is negative reciprocal */
  out[2] = - ( (Q[0]-P[0]) / (Q[1]-P[1]) );
  return !finite(out[2]);
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
  else if( finite(L1[2]) ) {  /* L2 is vertical; L1 is not */
    out[0] = L2[0];
    out[1] = L1[1] + (L2[0] - L1[0]) * L1[2];
    return 0;
  } else if( finite(L2[2]) ) { /* L1 is vertical; L2 is not */
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
    
    r = fl_segment_deluxe_dist(p0, p1, v, v1);

    if((v->line->fc0->world->verbosity - (r>=0))>=6)
      printf("\nfl_segment_deluxe_dist returned %g for vertex %d\n",r,v1->label);

    if(r<=0 || !finite(r)) {
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
	printf("len=%g for vertex %d\n",len,v1->label);

      if(len <= 0) {

	if(v->line->fc0->world->verbosity == 4) 
	  printf("len=%g for vertex %d\n",len,v1->label);

	horde->stuff[i] = 0;
	crunch=1;
      }
      else {
	scale_3d(v1->scr,v1->scr,r/len);
	
	v1->a = ATAN2(v1->scr[1],v1->scr[0]);
	v1->r = r;
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

  if(((VERTEX *)a)->a == ((VERTEX *)b)->a )
    if( ((VERTEX *)a)->r==((VERTEX *)b)->r )
      return 0;
    else 
      return ( ((VERTEX *)a)->r < ((VERTEX *)b)->r ) ? -1 : 1;
	    
  return ( ((VERTEX *)a)->a < ((VERTEX *)b)->a ) ? -1 : 1;
}

void sort_by_angle_2d(DUMBLIST *horde) {
  int i;

  /* Do finiteness checking... */
  for(i=0;i<horde->n;i++) {
    VERTEX *v1 = (VERTEX *)(horde->stuff[i]);
    if(!v1 || !finite(v1->a) || !finite(v1->r)) {
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
  dumblist_crunch(horde,angle_cmp);
}  


/**********************************************************************
 * hv_cp - copy a HULL_VERTEX 
 */
inline void hv_cp(HULL_VERTEX *to, HULL_VERTEX *from) {
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
  int verbosity;
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

  /* Sort out debugging level. Kind of a kludge, but what the hey... */
  { 
    WORLD *w = 0;
    for(i=0;!w && i<horde->n; i++) 
      if( (((VERTEX **)(horde->stuff))[i])->line )
	w = (((VERTEX **)(horde->stuff))[i])->line->fc0->world;
    verbosity = w?w->verbosity:0;
  }

  sort_by_angle_2d(horde);

  if(verbosity >= 4) 
    printf("hull_2d: got %d candidates...\n",horde->n);
  if(verbosity >= 5) {
    printf("candidate order:\n");
    for (i=0;i<horde->n;i++) {
      VERTEX *v = ((VERTEX **)(horde->stuff))[i];
      printf("pos %3.3d: (x,y)=(%5.2g,%5.2g); a=%5.2g; r=%5.2g; label=%d\n",
	     i,
	     v->scr[0],
	     v->scr[1],
	     v->a*180/3.14159,
	     v->r,
	     v->label
	     );
    }
  }



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
     !finite(out[i_r].p[0]) || !finite(out[i_r].p[1])) {
    out[i_r].p[0] = out[i_r].p[1] = 0;
    out[i_r].open = 1;
  } else
    out[i_r].open = 0;

  if(verbosity >= 5) {
    printf("Entering main loop: n=%d; horde->n=%d; i_r=%d; i_l=%d\n",n,horde->n,i_r,i_l);
  }
  /* Main loop */
  do {
    NUM a,b;
    if(horde->stuff[i]) {    /* Skip over already-zeroed elements */

      char open=0;
      
      /*** Find next guy and prep his bisector if necessary ***/
      for(i_l=MOD_NEXT(i,n);  i_l != i && !horde->stuff[i_l]; MOD_INC(i_l,n))
	;
      lv = (VERTEX *)(horde->stuff[i_l]);
      if(!been_there)
	perp_bisector_2d( out[i_l].bisector, 0, lv->scr );
      
      if(verbosity >= 5) 
	printf("\tpos %3d: i_r=%3d, i_l=%3d, l_scr:%5.3g,%5.3g (%5.3g deg), scr:%5.3g,%5.3g (%5.3g deg), r_scr:%5.3g,%5.3g (%5.3g deg)\n\t\t",i,i_r,i_l,lv->scr[0],lv->scr[1],lv->a*180/3.14159,iv->scr[0],iv->scr[1],iv->a*180/3.14159,rv->scr[0],rv->scr[1],rv->a*180/3.14159);


      /*** Check for pathologies ***/
      if(i_r == i || i_l == i) {
	abort = 1;
	//	fflush(stdout);
	//	fprintf(stderr,"hull_2d: eliminated all vertices! I give up.\n");
	//	horde->n=0;
	//	return;
      }

      if(!finite(iv->r) || !finite(iv->a))
	goto reject;

      if(iv->r == 0)
	goto reject;

      /*** Check colinearity and opentude on right & left ***/
      /* (right could be handled by caching, but it's pretty cheap) */
      b = iv->a - rv->a;
      TRIM_ANGLE(b);
      if( b < EPSILON && b > -EPSILON && rv->r < iv->r ) {
	if(verbosity>=5)
	  printf("REJECT: colinear and more distant than right side");
	goto reject; /* Colinear with (and farther than) right-side vertex */
      }

      a = lv->a - iv->a;

      TRIM_ANGLE(a);
      if( a < EPSILON && a > -EPSILON && lv->r <= iv->r ) {
	if(verbosity>=5)
	  printf("REJECT: colinear and more distant than left side");
	goto reject;  /* Colinear with (and farther than) left-side vertex */
      }

      if( a > EPSILON && a < PI-EPSILON) {
	/* Line isn't open on the left, so it intersects on the left.   */
	intersection_2d( out[ i ].p, out[ i ].bisector, out[i_l].bisector );
	out[i].open = 0;
	

	if( b > EPSILON && b < PI-EPSILON) {

	  /* Line isn't open on the right either, so compare intersections */
	  a = cross_2d( out[i_r].p, out[ i ].p );

	  if( !finite(a) ) {
	    printf("ASSERTION FAILED in hull_2d: intersection should exist (but doesn't)\n\t\tleft_p: %5.2g, %5.2g ;   right_p: %5.2g, %5.2g",out[i].p[0],out[i].p[1],out[i_r].p[0],out[i_r].p[1]);
	    printf("\n\t\tleft_b: %5.2g, %5.2g, %5.2g;  this_b: %5.2g, %5.2g, %5.2g; right_b: %5.2g, %5.2g, %5.2g\n",out[i_l].bisector[0],out[i_l].bisector[1],out[i_l].bisector[2],out[i].bisector[0],out[i].bisector[1],out[i].bisector[2],out[i_r].bisector[0],out[i_r].bisector[1],out[i_r].bisector[2]);
	    printf("\t\tlv->a=%5.2g, iv->a=%5.2g, rv->a=%5.2g,   lv-iv: %5.2g,   iv-rv: %5.2g\n",lv->a*180/PI,iv->a*180/PI,rv->a*180/PI,(lv->a-iv->a)*180/PI,(iv->a-rv->a)*180/PI);
	    goto reject;
	  }

	  /*** If intersections were in the wrong order, reject the point ***/
	  if(a < 0) {
	    if(verbosity>=5)
	      printf("REJECT: wrong order");
	    goto reject;
	  }
	}

      } else {
	/* Line is open on the left -- set the left point to zero */
	if(verbosity >= 5) printf(" looks open... ");
	out[i].p[0] = out[i].p[1] = 0;
	out[i].open = 1;
      }


      /*** Acceptance code ***/
      been_there = (been_there || (i_l < i));
      i_r = i;
      i = i_l;
      rv = iv;
      iv = lv;

      if(verbosity>=5) 
	printf("Accept.");
      
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

      if(verbosity >= 5) 
	printf(" terminus=%d; been_there=%d\n",terminus,been_there);

    } /* end of non-zeroed-out check */
    
  } while(!abort && (!been_there || i != terminus));  /* End of main loop */

  /* Finished -- now crunch the dumblist and, simultaneously, the output. */
  { 
    int j;

    for(i=j=0; i<horde->n; i++) {
      VERTEX *hsi = (VERTEX *)(horde->stuff[i]);
      if(hsi) {
	if(i != j) {
	  if(verbosity >5) printf("  %d:%d->%d ",hsi,i,j);
	  horde->stuff[j] = hsi;
	  out[j] = out[i];
	}
	j++;
      }
    }
    horde->n = j;
    if(verbosity>5) printf("\n");
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
      if(verbosity>5) 
	printf(" %d(%c):a_l=%.3g,a_r=%.3g,p1=%.3g,p2=%.3g  ",i,(hv->open)?'o':'c',180/PI*hv->a_l,180/PI*hv->a_r,hv->p[0],hv->p[1]);
    }
    if(verbosity>5)
      printf("\n");
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
#if DEBUG_HULL
static int verbosity = 0;
#endif

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

  /* Check for colinearity -- easy fix if we're colinear */
  /* The colinear_keep flag gets return value that will be given if
   * the vertex turns out to be open - this is sort of the opposite sense 
   * of the return value.
   */


  { 
    int colinear_keep = 0;
    NUM eps = EPSILON * v->r;
    if( fabs( cross_2d(v->scr, pv->scr) ) < eps && (inner_2d(v->scr, pv->scr) > 0) ) {
      if(  v->r >= pv->r  ) {
#if DEBUG_HULL
	if(verbosity >= 5) 
	  printf("\tColinear & more distant than prior -- reject\n");
#endif
	return 0;
      }
      else
	colinear_keep = 1;
    }

    if( fabs( cross_2d( v->scr, nv->scr) ) < eps && (inner_2d(v->scr, nv->scr) > 0) )  {
      if(  v->r >= nv->r  ) {
#if DEBUG_HULL
	if(verbosity >= 5) 
	  printf("\tColinear & more distant than next -- reject\n");
#endif
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
#if DEBUG_HULL
      if(verbosity >= 5) 
	printf("\tColinear & closer than next or prior - keep\n");
#endif
      return ( phv->open ? colinear_keep : -1 );
    }
  } /* end of colinear check convenience block */

  /* Not colinear -- find perpendicular bisector and start crunching... */
  
  perp_bisector_2d(scr->bisector, 0, v->scr);

  /* scr->p gets the hull intersection point with the next guy; */
  /* point gets the hull intersection point with the previous guy */
  
  intersection_2d( scr->p, scr->bisector, nhv->bisector );
  intersection_2d( point, scr->bisector, phv->bisector );


#if DEBUG_HULL
  if(verbosity >= 5) {
    printf("   perp_bisector_2d: v->scr is (%g,%g), scr->bisector is (%g,%g,%g)\n",v->scr[0],v->scr[1],scr->bisector[0],scr->bisector[1],scr->bisector[2]);
    printf("   intersections: prev bisector is (%g,%g,%g), intersection is (%g,%g); next bisector is (%g,%g,%g); intersection is (%g,%g)\n",
	   phv->bisector[0],phv->bisector[1],phv->bisector[2],
	   point[0],point[1],
	   nhv->bisector[0],nhv->bisector[1],nhv->bisector[2],
	   scr->p[0],scr->p[1]
	   );
  }
#endif
  

  /**** Open on next side ? ****/
  if(          cross_2d(  v->scr, nv->scr) <= 0) {
    if(        cross_2d(  pv->scr, v->scr) <= 0) {
#if DEBUG_HULL
      if(verbosity >= 5) 
	printf("\tAntiparallel case (open both sides) -- keep\n");
#endif
      open = 3;
    } else {
#if DEBUG_HULL
      if(verbosity >= 5) 
	printf("\tNoncolinear and open on next side -- keep\n");
#endif
      open = 1;
    }
  } 

  /**** Open on previous side ? ****/
  else if(     cross_2d( pv->scr,  v->scr) <= 0) {
#if DEBUG_HULL
    if(verbosity >= 5) 
      printf("\tNoncolinear and open on prev side -- keep\n");
#endif
    open = 2;
  }
  
  /**** Closed and OK? ****/
  else if(  cross_2d( point, scr->p ) >= 0  ) {
#if DEBUG_HULL
    if(verbosity >= 5) 
      printf("\tNoncolinear and closed -- keep (ok)\n");
#endif
    open = -1;
  } 

  /**** Closed and not OK? ****/
  else { 
    /* closed both sides, but intersections in wrong order -- reject */
#if DEBUG_HULL
    if(verbosity >= 5) {
      printf("\tNoncolinear and wrong vertex order -- reject\n");
      printf("\t ( prev intersection = (%g,%g), next intersection = (%g,%g) ) \n",point[0],point[1],scr->p[0],scr->p[1]);
      printf("scr->p is (%g,%g); nhv->p is (%g,%g)\n",scr->p[0], scr->p[1],nhv->p[0],nhv->p[1]);
      printf("scr->bisector is (%g,%g,%g); nhv->bisector is (%g,%g,%g)",scr->bisector[0],scr->bisector[1],scr->bisector[2], nhv->bisector[0],nhv->bisector[1],nhv->bisector[2]);
    }
#endif
      return 0;
  }
  
  
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
 */
static DUMBLIST *ws = 0;

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


  for(i=0;!w && i<horde->n; i++) 
    if( (((VERTEX **)(horde->stuff))[i])->line )
      w = (((VERTEX **)(horde->stuff))[i])->line->fc0->world;
  if(!w) {
    fprintf(stderr,"You're in trouble, guv! exiting...\n");
    exit(99);
  }

#if DEBUG_HULL
  verbosity = w->verbosity;
  if(!verbosity && central_v->label==32)
    verbosity = 5;
#endif

  passno = ++w->passno;

#if DEBUG_HULL
  if(verbosity >= 5) {
    printf("\n\n********************************* hull_2d_us: vertex %d (0x%x)\n",central_v->label, central_v);
  }
#endif

  /* Set our output list to zero length */
  ws->n = 0;

  /* Nothing in the horde?  We're done! */
  if(horde->n == 0)
    return;

  /* always add the first eligible VERTEX in the list... */
  for(i=0;i<horde->n && 
	( horde->stuff[i] == central_v ||             // skip the main vertex if present
	  !( ((VERTEX *)(horde->stuff[i]))->next ) || // skip endpoint vertices if present
	  !( ((VERTEX *)(horde->stuff[i]))->line)     // skip image vertices
	  );
      i++) 
    ;
  if(i<horde->n) {
    ((VERTEX **)(horde->stuff))[i]->passno = passno;
    ws->stuff[0] = horde->stuff[i];
    ws->n = 1;
    hull[0].open = 1;
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

    if(v->passno == passno || v == central_v || !v->next ) {
 #if DEBUG_HULL
      if(verbosity >= 5)
	printf("vertex %d has already been inspected (or is the same as %d), or has no next field (0x%x) - skipping...\n",v->label, central_v->label, v->next);
#endif
    } else {
      v->passno = passno; // Mark it processed (avoid duplicated effort)
      
      keep = 0;

#if DEBUG_HULL
      if(verbosity >= 5) {
	int i;
	printf("Considering vertex %d; a is %g degrees, r is %g...",v->label,v->a * 180/3.1415926, v->r);
	printf("\t(current list: ");
	for(i=0;i<ws->n;i++) {
	  VERTEX *v = (VERTEX *)(ws->stuff[i]);
	  printf(" %d:L%d(A%g,R%g) ",i,v->label,180/PI*v->a,v->r);
	}
	printf(")\n");
      }
#endif

      /* FIXME: Binary search later -- linear search now */
      for(next_idx = 0; 
	  next_idx < ws->n && 
	    (nv = (((VERTEX **)(ws->stuff))[next_idx]) )->a < v->a;  // assign to nv
	  next_idx++ 
	  )
	;

#if DEBUG_HULL      
      if(verbosity >= 5) 
	printf("position %d in current hull (out of %d)\n",next_idx,ws->n);
#endif

      /******************************
       * Now set nv to the next vertex after this one and pv to the previous one.
       * next_idx is the true insertion point if this element is kept.
       */

      if( nv->a < v->a ) {
	/* We went off the end... */
#if DEBUG_HULL	
	if(verbosity >= 5) 
	  printf("   off the end...\n");
#endif
	next_idx = ws->n;
	nv = ((VERTEX **)(ws->stuff))[0];
	nh = hull;
	pv = ((VERTEX **)(ws->stuff))[ws->n-1];
	ph = hull + ws->n-1;
      } else {
	/* We didn't... */
	if(next_idx==0) {
#if DEBUG_HULL
	  if(verbosity >= 5)
	    printf("   idx=0...\n");
#endif
	  pv = ((VERTEX **)(ws->stuff))[ws->n-1];
	  ph = hull + ws->n-1;
	  nh = hull;
	} else {
#if DEBUG_HULL
	  if(verbosity >= 5)
	    printf("   idx>0...\n");
#endif
	  pv = ((VERTEX **)(ws->stuff))[next_idx - 1];
	  ph = hull + next_idx-1;
	  nh = hull + next_idx;
	}
      }
      
      flag = check_hullpoint(v, pv,ph, nv,nh, &scrhv);

#if DEBUG_HULL
      if(verbosity>=5) 
	printf("   %s\n",(flag ? "KEEP" : "REJECT"));
#endif
      
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
#if DEBUG_HULL
	if(verbosity >= 5) 
	  printf("Inserting into position %d\n",next_idx);
#endif

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

#if DEBUG_HULL	  
	  if(verbosity >= 5) 
	    printf("Walking backward; n=%d... ",n);
#endif
	  /* Walk backward until we find a still-keep-worthy vertex */
	  while( p != n && 
		 ! (check_flag = check_hullpoint( (VERTEX *)(ws->stuff[p]), 
						  (VERTEX *)(ws->stuff[pp]), hull + pp, 
						  (VERTEX *)(ws->stuff[n]),  hull + n,
						  hull + p))
		 ) {
#if DEBUG_HULL
	    if(verbosity >= 5) 
	      printf("p=%d ");
#endif

	    p=pp;
	    MOD_DEC(pp, ws->n);
	  }

	  hull[p].open = (check_flag > 0) && (check_flag & 1);
	  hull[pp].open = (check_flag > 0) && (check_flag & 2);

#if DEBUG_HULL
	  if(verbosity >= 5) 
	    printf("final p =%d  prev_intersection=(%g,%g) \nWalking forward; n=%d... ",p,hull[p].p[0],hull[p].p[1],n);
#endif
	  
	  /* walk forward until we find a still-keep-worthy vertex */
	  while( n != p &&
		 ! (check_flag = check_hullpoint( (VERTEX *)(ws->stuff[nn]),
						  (VERTEX *)(ws->stuff[n]),  hull + n,
						  (VERTEX *)(ws->stuff[nnn]), hull + nnn,
						  hull+nn
						  ) )
		 ) {

#if DEBUG_HULL
 	    if(verbosity >= 5)
	      printf("nn=%d ",nn);
#endif

	    nn=nnn;
	    MOD_INC(nnn, ws->n);
	  }

	  hull[nn].open = (check_flag > 0) && (check_flag & 1);
	  hull[n].open = (check_flag > 0) && (check_flag & 2);

#if DEBUG_HULL
	  if(verbosity >= 5) {
	    printf("final nn=%d;  next_intersection=(%g,%g)   \n",nn,hull[n].p[0],hull[n].p[1]);
	    {
	      NUM pp[2];
	      NUM bis1[3];
	      NUM bis2[3];
	      perp_bisector_2d(bis1,0, ((VERTEX *)(ws->stuff[n]))->scr);
	      perp_bisector_2d(bis2,0, ((VERTEX *)(ws->stuff[nn]))->scr);
	      intersection_2d(pp, bis1, bis2);
	      printf("\tdirect intersection: bisector[n]=(%g,%g,%g); bisector[nn]=(%g,%g,%g); point is: (%g,%g)\n",
		     bis1[0],bis1[1],bis1[2],
		     bis2[0],bis2[1],bis2[2],
		     pp[0],pp[1]);
	    }
	  }
#endif

	  /* Now p is the first keepworthy vertex, walking backward from n, and */
	  /* nn is the first keepworthy vertex, walking forward from n.         */
	  /* Crunch down the list... */

#if DEBUG_HULL
	  if(verbosity >= 5) {
	    int ii;
	    printf("Kept & trimmed.  List is: ");
	    for(ii=0;ii<ws->n;ii++) {
	      printf(" %d ",((VERTEX *)(ws->stuff[ii]))->label);
	    }
	    printf("\n\tp=%d, n=%d, nn=%d\n",p,n,nn);
	  }
#endif

	  {

	    int ii,jj, diff;
	    
	    if( p <= n-1 ) {
	      if( nn <=
 n-1 ) {

#if DEBUG_HULL
		if(verbosity >= 5) printf("...*****...|.....    case\n");
#endif

		for(ii=0, jj=nn; jj<=p; jj++,ii++) {
		  hv_cp( hull+ii, hull+jj );
		  ws->stuff[ii] = ws->stuff[jj];
		}
		hv_cp( hull+ii, hull +n );
		ws->stuff[ii] = ws->stuff[n];
		ii++;

	      } else {

#if DEBUG_HULL
		if(verbosity >= 5) printf("*****...|...**** case\n");
#endif
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
#if DEBUG_HULL
		if(verbosity >= 5) printf("\tii=%d\n",ii);
#endif
	      }
	    } else {
	      
#if DEBUG_HULL
	      if(verbosity >= 5) printf("....|...****.. case\n");
#endif

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

#if DEBUG_HULL
	    if(verbosity >= 5) {
	      int ii;
	      printf("After post-keep trimming, we have: ");
	      for(ii=0;ii<ws->n;ii++) {
		printf(" %d %",((VERTEX *)(ws->stuff[ii]))->label);
	      }
	      printf("\n");
	    }
#endif
	  } /* end of crunching convenience block */
	} /* end of neighbor-testing convenience block */
      } else {
	/* no keeping -- just ignore this point and start the next loop */
      }
    } /* end of non-duplication test block */

#if DEBUG_HULL
    if(w->verbosity >= 3) {
      printf("horde-loop: finished %d; %d neighbors in list...\n",v->label,ws->n);
    }
#endif
    
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

#if DEBUG_HULL
  if(w->verbosity >= 5)  {
    int i;
    
    for(i=0;i<ws->n;i++) {
      VERTEX *vx = ((VERTEX *)(ws->stuff[i]));
      printf("\t%d: a=%g, r=%g, point=(%g,%g), hullpt=(%g,%g), open=%d, a_l = %g, a_r = %g\n",
	     vx->label,
	     vx->a,
	     vx->r, 
	     vx->scr[0],
	     vx->scr[1],
	     hull[i].p[0],
	     hull[i].p[1],
	     hull[i].open,
	     hull[i].a_r,
	     hull[i].a_l
	     );
    }
  }
#endif
}


