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
#define NAN (nan("NAN"))

char *code_info_geometry="%%%FILE%%%";

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
  if( finite(L1[2]) && finite(L2[2]) ) {
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
	
	v1->a = atan2(v1->scr[1],v1->scr[0]);
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
      dumblist_rm(horde,i);
      printf("X"); fflush(stdout);
      i--;
    }
  }

  /* Do the actual sort... */
  dumblist_sort(horde,angle_cmp);
  dumblist_crunch(horde,angle_cmp);
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
  NUM epsilon=1e-19; // double-precision roundoff
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
	hv->a_r = hv->a_l = atan2(hv->p[1], hv->p[0]);
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

#ifdef never 
/**********************************************************************
 * hull_2d_older
 * 
 * Given an angle-sorted DUMBLIST of VERTICES containing 2-D points in their 
 * scratch spaces, find the convex hull of the origin.  
 * The original DUMBLIST winds up containing the points that were used 
 * for the hull; you can also pass in a DUMBLIST
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
 *
 * hull_2d_older is an older version of hull_2d that is kept around for 
 * reference.
 * *
 */


void hull_2d_older(HULL_VERTEX *out, DUMBLIST *horde, DUMBLIST *rejects) {
  NUM epsilon=1e-19; // double-precision roundoff
  NUM left[2], right[2], b1[3], *b0;
  NUM a,b;
  int outdex = 0;
  char cache_ok = 0;
  char reject_collinear;
  char redo = 0;
  int j = -1;         /* j gets -1 or the next rightward VERTEX */
  int k = -1;         /* k gets -1 or the next leftward VERTEX */
  int i;

  int verbosity;
  
  if(horde->n < 1)
    return;

  sort_by_angle_2d(horde);

  verbosity = (((VERTEX **)(horde->stuff))[0])->line ? 
    (((VERTEX **)(horde->stuff))[0])->line->fc0->world->verbosity :
    (horde->n>1) ? 
       (((VERTEX **)(horde->stuff))[1])->line->fc0->world->verbosity : 
    0;
  verbosity = 5;

  if( verbosity >= 4 ) {
    int i;
    printf("\nhull_2d: after sort, the horde has %d elements:\n\t\t",horde->n);
    for(i=0;i<horde->n;i++) {
      printf("%5.2g(V%4d)   ",180/PI*( ((VERTEX **)(horde->stuff))[i] )->a, ( ((VERTEX **)(horde->stuff))[i])->label);
    }
    printf("\n");
  }

  /* Assemble the intersection points with neighbors, rejecting as necessary
   * on the way.  Algorithm:
   *   - Calculate left vertex.
   *   - Check if you have a cached right vertex and if not, calculate one.
   *   - Compare cached right vertex and newly calculated left vertex
   *   - Reject if necessary, voiding the cache.
   *   - Repeat until you get to an old, accepted vertex.
   */
  
  i=0;
  
  if(verbosity >= 5) 
    printf("hull_2d: processing...\n");
  
  do {
    NUM *ivec,*jvec,*kvec;
    int missing;
    
    for(missing=i=0;i<horde->n;i++) {
      if(!horde->stuff[i]){   /* On redos, skip missing elements */
	if(missing >= horde->n){
	  fprintf(stderr,"hull bug:  eliminated all vertices!\n");
	  exit(53);
	}
      } else if(!finite( ((VERTEX *)horde->stuff[i])->a)) {
	/* Skip degenerate points too -- this should't be necessary
	   if they've been pre-winnowed; but it's here as a safety check.
	*/ 
	printf("degen or nan trigger in hull_2d\n"); 
	horde->stuff[i] = 0;
      } else {
	
                     if(verbosity>=5)  printf("\t i=%3d, cache %s ",i,(cache_ok?"ok ":"bad"));
	
	/* Neither infinite nor zero -- handle normally */
	
	ivec = ((VERTEX *)(horde->stuff[i]))->scr;
	
	/* Grab perpendicular bisector for this guy */
	b0 = &(out[outdex].bisector[0]);
	perp_bisector_2d(b0, 0, ((VERTEX *)((horde->stuff[i])))->scr);
	
	/* Check cache and if necessary calculate the right side */
	if(cache_ok) { 
	  /* Last item was kept -- copy its vertex info into the LHS. */
	  right[0] = left[0];
	  right[1] = left[1];
	  jvec = ((VERTEX *)(horde->stuff[j]))->scr;
	} else {
	  /* Last item wasn't kept -- find the previous vertex and 
	     calculate the intersection for the LHS. */
	  if(j < 0) {
	    for(j = ( (i>0) ? (i-1) : (horde->n - 1) );
		(!horde->stuff[j]) && j!=i;
		j = ( (j>0) ? (j-1) : (horde->n - 1) )
		)
	      ;
	  }
	  jvec = ((VERTEX *)(horde->stuff[j]))->scr;
	  if(j==i) 
	    right[1] = right[0] = NAN;
	  else {
	    perp_bisector_2d(b1, 0, jvec);
	    intersection_2d(right, b0, b1);
	  }
	}
	
	if(verbosity>=5) printf("j=%3d ",j);
      
      /* Find next VERTEX */
      for( k= (i+1)%(horde->n);
	   (!horde->stuff[k]) && k != j;
	     k = ( (k+1) % horde->n )
	   )
	;
      
      kvec = ((VERTEX *)(horde->stuff[k]))->scr;

      if(verbosity >=5) printf("k=%3d ",k);
      
      if(i==j)
	reject_collinear=0;
      else {
	NUM inorm = norm_2d(ivec); 
	NUM knorm = norm_2d(kvec);
	NUM jnorm = norm_2d(jvec);
	perp_bisector_2d(b1, 0, kvec);
	intersection_2d(left,b0,b1);
	
	/* Check for collinearity with neighbors, either the 
	 *  left or right side. 
	 * It would shave a few usec to stash knorm somewhere and
	 * reuse for inorm -- but I can't be bothered. */
	
	reject_collinear = 
	  (inorm < 1e-6) 
	  ||
	  ((jnorm <= inorm) && jnorm != 0 && 
	   fl_eq(inner_2d( ivec,jvec )
		 , jnorm * inorm
		 )
	   )
	  ||
	  ((knorm <= inorm) && knorm != 0 &&
	   fl_eq(inner_2d( ivec,kvec )
		 , knorm * inorm
		 ) 
	   );
	if(verbosity >= 5) printf("%s collinear ",reject_collinear?"***":"not");
      }
      
      /* Only reject if the neighbors are less than 180 deg from each 
       * other and the right and left neighboring vertices are in the wrong
       * order (the neighbor bisectors crossed before hitting this one). 
       * 
       * Alternatively, reject if either neighbor is collinear and shorter
       * than the current one.
       */
      
      a = cross_2d(jvec,kvec);
      b = cross_2d(right,left);
      //	if( ( (!(j==k)) && (a * b < 0) )
      if( ( (!(j==k)) &&
	    finite(a) && finite(b) && (a > 0 && b < 0))
	  ||
	  reject_collinear
	  ) {

	if(verbosity >= 5) printf(" REJECT ");

	/* reject */
	if(rejects) 
	  dumblist_add(rejects,horde->stuff[i]);
	
	horde->stuff[i] = 0;
	cache_ok = 0;
	
	if(i == (horde->n - 1)) {
	  horde->n--;
	  redo = 1;
	}
	
	/* Back up to the previous accepted vertex, if there is one. */
	if(j<i) {
	  i=j-1; /* the one gets added back in by the for loop */
	  outdex--;
	}
	j=-1;
      } else {  /* End of rejection code; start of acceptance code */
	
	if(verbosity >= 5) printf(" ACCEPT ");

	/* accept */
	out[outdex].p[0] = left[0];
	out[outdex].p[1] = left[1];
	
	/* assignment below */
	if ( out[outdex].open = !(finite(left[0]) && finite(left[1]) 
				  && fabs(cross_2d(ivec,kvec)>epsilon)) ) {
	  out[outdex].a_l = ((VERTEX *)(horde->stuff[i]))->a + PI/2;
	  {
	    int k;
	    for(k=MOD_NEXT(i,horde->n);k!=i && !horde->stuff[k]; MOD_INC(k,horde->n))
	      ;
	    out[outdex].a_r = ((VERTEX *)(horde->stuff[k]))->a - PI/2;
	  }
	  if(out[outdex].a_l > PI) 
	    out[outdex].a_l -= 2*PI;
	  if(out[outdex].a_r < -PI)
	    out[outdex].a_r += 2*PI;

	} else 
	  out[outdex].a_l = out[outdex].a_r = atan2(out[outdex].p[1],out[outdex].p[0]);
	
	outdex++;
	cache_ok = 1;
	j=i;
	
	if(redo) {
	  /* To get here, we have to be redoing the first part because
	   * the last item in the list was rejected.  We only have to 
	   * redo stuff until the first acceptance, whereupon we are done.
	   */
	  redo = 0;
	  break;
	}
      }          /* End of acceptance case*/
      
      if(verbosity >= 5) printf("\n");
      }            /* End of evaluation case (non-skip) */
    }              /* End of for loop */
    
    
    
    /* Fix up the horde -- crunch it down to size. */
    for(j=i=0;i<horde->n;i++) {
      if(horde->stuff[i]) {
	if(j<i)
	  horde->stuff[j] = horde->stuff[i];
	j++;
      }
    }
    horde->n = j;
  } while(redo && (horde->n > 0));
  
  if(horde->n ==0) {
    fprintf(stderr,"Eliminated all vertices!\n");
    exit(57);
  }
  
  /* Assertion test -- this should never happen. */
  if(horde->n != outdex) {
    fprintf(stderr,"Assertion failed!  outdex != horde->n (%d != %d)!\n",outdex,horde->n);
  }
  
  /* Clean up doubly-open 0-element case */
  if(horde->n == 1 || ( out[0].open && out[horde->n-1].open ) ) {
    out[0].p[0] = ((VERTEX *)(horde->stuff[0]))->scr[0] / 2;
    out[0].p[1] = ((VERTEX *)(horde->stuff[0]))->scr[1] / 2;
  }
  
  
  return;
}
#endif

