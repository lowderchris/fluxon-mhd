/* Geometry header files -- for handling point locations &c.
 *
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
 */

#ifndef FLEM_GEOMETRY
#define FLEM_GEOMETRY 1

#include "data.h"

inline NUM norm_2d(NUM *x);
inline NUM norm_3d(NUM *x);

inline NUM norm2_3d(NUM *x);

inline NUM inner_2d(NUM *p0, NUM *p1);
inline NUM inner_3d(NUM *p0, NUM *p1);

inline NUM cross_2d(NUM *p0, NUM *p1);
#define cross_3d cross
inline void *cross(NUM *out, NUM *p0, NUM *p1);

inline void scale_3d(NUM *out, NUM *a, NUM alpha);

inline void sum_3d(NUM *out, NUM *a, NUM *b);
inline void diff_3d(NUM *out, NUM *a, NUM *b);

inline void cp_3d(NUM *out, NUM *a);


/**********************************************************************
 * Matrix manipulation & rotation matrix generation
 */
inline void rotmat_2d(NUM *out, NUM alpha);
inline void rotmat_2d_fr_slope(NUM *out, NUM dy, NUM dx);
inline void mat_mult_2d(NUM *out, NUM *a, NUM *b);

inline void transpose_2x2(NUM *mat);
inline void transpose_3x3(NUM *mat);

inline void mat_mult_3d(NUM *out, NUM *a, NUM *b);
inline void mat_vmult_2d(NUM *out, NUM *mat, NUM *v);
inline void mat_vmult_3d(NUM *out, NUM *mat, NUM *v);
inline void vec_mmult_3d(NUM *out, NUM *mat, NUM *v);

inline NUM det_2d(NUM *mat);
inline NUM det_3d(NUM *mat);

/**********************************************************************
 * Projection functions
 */
inline void projmatrix(NUM *out, NUM *x0_3, NUM *x1_3);

/**********************************************************************
 * Distance functions
 */
inline NUM cart2_2d(NUM *x1, NUM *x2);
inline NUM cart2_3d(NUM *x1, NUM *x2);

inline NUM cart_2d(NUM *x1, NUM *x2);
inline NUM cart_3d(NUM *x1, NUM *x2);

inline NUM p_l_dist(NUM *p0, NUM *x0, NUM *x1);  /* point-line dist */
inline NUM p_ls_dist(NUM *p0, NUM *x0, NUM *x1); /* point-line segment dist */
inline NUM l_l_dist( NUM a0[3], NUM b0[3], NUM c0[3], NUM d0[3]); /* line-line dist */
inline void p_ls_closest_approach(NUM p0[3], NUM a0[3], NUM b0[3], NUM c0[3]);
inline void ls_closest_approach(NUM p0[3], NUM p1[3], NUM a0[3], NUM b0[3], NUM c0[3], NUM d0[3]);
inline NUM ls_ls_dist(NUM a0[3], NUM b0[3], NUM c0[3], NUM d0[3]); /* seg-seg dist */

inline NUM fl_segment_masked_dist(VERTEX *v0, VERTEX *v1); 
inline NUM fl_segment_masked_deluxe_dist(NUM P0[3], NUM P1[3], VERTEX *v0, VERTEX *v1);

inline NUM fl_segment_dist(VERTEX *v1, VERTEX *v2);        
inline NUM fl_segment_deluxe_dist(NUM P0[3], NUM P1[3], VERTEX *v0, VERTEX *v1);





/********************************************************************** 
 * Voronoi-cell functions
 */


inline int perp_bisector_2d(NUM *out, NUM *P, NUM *Q);           /* Find the perpendicular bisector
							     of a line segment */
inline int intersection_2d(NUM *out, NUM *L1, NUM *L2);          /* Find the intersection of two
							     line segments */

NUM neighbor_triangle(NUM *dAdp, NUM *L, NUM *M, NUM *N);  /* Characterize the 
							 partial-neighborhood triangle
							 described by L, M, and N (with
							 central point at the origin).
							  */

int find_neighbors_from_list(VERTEX *tp, 
			     int n, 
			     int depth, 
			     VERTEX **where);


typedef struct HULL_VERTEX {
  NUM p[2];  /* point coordinates */
  NUM a;     /* point angle -- cached here to avoid multiple atan2 calls */
  NUM bisector[3]; /* bisector point -- cached here for use later */
  char open; /* indicates point is at infinity */
} HULL_VERTEX;

void hull_2d(HULL_VERTEX *out, DUMBLIST *horde, DUMBLIST *rejects);
#endif /* overall file include */



#define PI 3.141592653589793238462643383279502
#define DEG2RAD (PI/180.)
#define RAD2DEG (180./PI)
