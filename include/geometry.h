/* Geometry header files -- for handling point locations &c.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Craig DeForet, 2004-2007
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
 * This file is part of FLUX 2.2 (22-Nov-2008).
 */

#ifndef FLEM_GEOMETRY
#define FLEM_GEOMETRY 1

#include "data.h"


NUM norm_2d(NUM *x);
NUM norm_3d(NUM *x);

NUM norm2_2d(NUM *x);
NUM norm2_3d(NUM *x);

NUM inner_2d(NUM *p0, NUM *p1);
NUM inner_3d(NUM *p0, NUM *p1);

NUM cross_2d(NUM *p0, NUM *p1);
#define cross_3d cross
void cross(NUM *out, NUM *p0, NUM *p1);

void scale_2d(NUM *out, NUM *a, NUM alpha);
void scale_3d(NUM *out, NUM *a, NUM alpha);

void diff_2d(NUM *out, NUM *a, NUM *b);

void sum_3d(NUM *out, NUM *a, NUM *b);
void diff_3d(NUM *out, NUM *a, NUM *b);

void cp_3d(NUM *out, NUM *a);

/**********************************************************************
 * Matrix manipulation & rotation matrix generation
 */
 void rotmat_2d(NUM *out, NUM alpha);
 void rotmat_2d_fr_slope(NUM *out, NUM dy, NUM dx);
 void mat_mult_2d(NUM *out, NUM *a, NUM *b);

 void transpose_2x2(NUM *mat);
 void transpose_3x3(NUM *mat);

 void mat_mult_3d(NUM *out, NUM *a, NUM *b);
 void mat_vmult_2d(NUM *out, NUM *mat, NUM *v);
 void mat_vmult_3d(NUM *out, NUM *mat, NUM *v);
 void vec_mmult_3d(NUM *out, NUM *mat, NUM *v);

 NUM det_2d(NUM *mat);
 NUM det_3d(NUM *mat);

/**********************************************************************
 * Planar & topology functions
 */
 void points2plane(PLANE *plane, NUM *p0, NUM *p1, NUM *p2);
 int p_l_intersection(NUM *out, PLANE *plane, NUM *p0, NUM *p1);
 int xy_l_intersection(NUM *out, NUM *p0, NUM *p1);
 int p_inside_tri(NUM *tri0, NUM *tri1, NUM *tri2, NUM *p);
 int trivloop(FLUXON *f);

/**********************************************************************
 * Projection functions
 */
 void projmatrix(NUM *out, NUM *x0_3, NUM *x1_3);
 void reflect(NUM out[3], NUM point[3], PLANE *plane);
/**********************************************************************
 * Distance functions
 */
 NUM cart2_2d(NUM *x1, NUM *x2);
 NUM cart2_3d(NUM *x1, NUM *x2);

 NUM cart_2d(NUM *x1, NUM *x2);
 NUM cart_3d(NUM *x1, NUM *x2);

 NUM p_l_dist(NUM *p0, NUM *x0, NUM *x1);  /* point-line dist */
 NUM p_ls_dist(NUM *p0, NUM *x0, NUM *x1); /* point-line segment dist */
 NUM l_l_dist( NUM a0[3], NUM b0[3], NUM c0[3], NUM d0[3]); /* line-line dist */
 void p_ls_closest_approach(NUM p0[3], NUM a0[3], NUM b0[3], NUM c0[3]);
 void ls_closest_approach(NUM p0[3], NUM p1[3], NUM a0[3], NUM b0[3], NUM c0[3], NUM d0[3]);
 NUM ls_ls_dist(NUM a0[3], NUM b0[3], NUM c0[3], NUM d0[3]); /* seg-seg dist */

 NUM fl_segment_masked_dist(VERTEX *v0, VERTEX *v1); 
 NUM fl_segment_masked_deluxe_dist(NUM P0[3], NUM P1[3], VERTEX *v0, VERTEX *v1);

 NUM fl_segment_dist(VERTEX *v1, VERTEX *v2);        
 NUM fl_segment_deluxe_dist(NUM P0[3], NUM P1[3], VERTEX *v0, VERTEX *v1);





/********************************************************************** 
 * Voronoi-cell functions
 */


 int perp_bisector_2d(NUM *out, NUM *P, NUM *Q);           /* Find the perpendicular bisector
							     of a line segment */
 int intersection_2d(NUM *out, NUM *L1, NUM *L2);          /* Find the intersection of two
							     line segments */
void project_n_fill(VERTEX *v, DUMBLIST *horde);


NUM neighbor_triangle(NUM *dAdp, NUM *L, NUM *M, NUM *N);  /* Characterize the 
							 partial-neighborhood triangle
							 described by L, M, and N (with
							 central point at the origin).
							  */

int find_neighbors_from_list(VERTEX *tp, 
			     int n, 
			     int depth, 
			     VERTEX **where);


/* The a_left and a_right differ when the point is open... */

/* Generally you only need 2 points for the hull position because all
   are projected onto a perpendicular plane, but I changed it to 3
   because of the hull vertices on the photosphere specifically. For
   the regular relaxing code, it should only use the first two
   coordinates. */
typedef struct HULL_VERTEX {
  NUM p[3];  /* point coordinates */
  NUM a_l;   /* absolute point angle when seen on left */
  NUM a_r;   /* absolute point angle when seen on right */
  NUM bisector[3]; /* bisector line description -- cached here for use later */
  char open;
} HULL_VERTEX;

void hull_2d(HULL_VERTEX *out, DUMBLIST *horde, DUMBLIST *rejects);
void hull_2d_us(HULL_VERTEX *hull, DUMBLIST *horde, VERTEX *central_v);


VERTEX   *find_vertex_by_location(POINT3D x, WORLD *w, VERTEX *v, int global);

int above_plane(POINT3D A, POINT3D B, POINT3D C, POINT3D X);
int in_simplex(POINT3D P0, POINT3D P1, POINT3D P2, POINT3D P3, POINT3D X);
DUMBLIST *find_simplex_by_location(POINT3D x, WORLD *w, VERTEX *v, int global);

NUM interpolate_lin_3d( POINT3D x, NUM p[12], NUM val[4], int n);
NUM interpolate_value_simplex( POINT3D x, DUMBLIST *dl, int val_offset);
NUM interpolate_value( POINT3D x, WORLD *w, VERTEX *v, int global, int val_offset);
/********************************************************************** 
 *Photosphere only functions
 */
void project_n_fill_photosphere(VERTEX *v, DUMBLIST *horde);
void hull_2d_us_photosphere(HULL_VERTEX *hull, DUMBLIST *horde, VERTEX *central_v);

/**********************************************************************/

NUM atan2_oct(NUM x, NUM y);
#define ATAN2 atan2_oct


#define PI 3.141592653589793238462643383279502
#define DEG2RAD (PI/180.)
#define RAD2DEG (180./PI)
#define EPSILON 1e-6


/* regularize an angle that is between -3PI and 3PI, to the range
 * -PI to PI.
 */
#define TRIM_ANGLE(a) do { if((a)<-PI) { (a)+=2*PI; } else if((a)>=PI) { (a) -= 2*PI; } } while(0)

#endif /* overall file include */
