/* Nitty-gritty modeling routines for FLUX -- routines that find
 * neighbors, calculate neighborhoods, and such.  These routines 
 * use the libraries in data.c, geometry.c, and io.c.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Craig DeForest, 2004-2007
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
 * This file is part of FLUX 2.0 (31-Oct-2007)
 */

#ifndef FLEM_MODEL
#define FLEM_MODEL 1

#include "data.h"
#include "geometry.h"
#include "io.h"

/**********************************************************************
 * "Magic" end-of-line boundary conditions
 */
void world_update_ends(WORLD *a);
void fluxon_update_ends(FLUXON *f);
void fluxon_auto_open(FLUXON *f);

/**********************************************************************
 * Neighborhood handling routines
 */

enum neighbor_global {
  normal_neighbors=0,
  global_neighbors=1,
  fast_neighbors=2,
  faster_neighbors  =3,
  gonzo_neighbors =4
};

int world_check(WORLD *a);  
void world_update_neighbors(WORLD *a, char global);
NUM *world_update_mag(WORLD *a, char global);
void world_relax_step(WORLD *a, NUM t);
void world_fluxon_length_check(WORLD *a, char global);

void fluxon_update_neighbors(FLUXON *fl, char global);
NUM *fluxon_update_mag(FLUXON *fl, char global, void ((**f_funcs)()), NUM *minmax);

void fluxon_calc_step(FLUXON *fl, NUM t);
void fluxon_relax_step(FLUXON *fl, NUM t);

HULL_VERTEX *vertex_update_neighbors(VERTEX *v, char global);

DUMBLIST *gather_neighbor_candidates(VERTEX *v,char global);
void image_find(PHOTOSPHERE *phot, VERTEX *image, VERTEX *v);

int winnow_cmp_1(void *a, void *b);
void winnow_neighbor_candidates(VERTEX *v, DUMBLIST *horde);

HULL_VERTEX *hull_neighbors(VERTEX *v, DUMBLIST *horde);

void voronoi_vertex(NUM P[2], NUM B[2], NUM A[2]);
int fix_proximity(VERTEX *V, NUM scale_thresh);
int fluxon_fix_proximity(FLUXON *F, NUM scale_thresh);
int global_fix_proximity(WORLD *w, NUM scale_thresh);

int fix_curvature(VERTEX *V, NUM curve_thresh_high, NUM curve_thresh_low);
int fluxon_fix_curvature(FLUXON *F, NUM curve_thresh_high, NUM curve_thresh_low);
int global_fix_curvature(WORLD *w, NUM curve_thresh_high, NUM curve_thresh_low);

void reconnect_vertices( VERTEX *v1, VERTEX *v2, long passno );
int vertex_recon_check( VERTEX *v1, long passno );
long fluxon_recon_check( FLUXON *f, long passno );
long global_recon_check( WORLD *w );

int fc_cancel( FLUX_CONCENTRATION *fc0, FLUX_CONCENTRATION *fc1 );

typedef struct VERTEX_STATS {
  long n;         /* Number of vertices included in stats */
  NUM f_acc;      /* average resultant magnitude of force */
  NUM f_max;      /* maximum resultant magnitude */
  NUM f_tot_acc;  /* average sum-of-magnitudes of forces  */
  NUM f_tot_max;  /* maximum resultant force */
  NUM n_acc;    /* average number of neighbors */
  NUM n_max;    /* maximum number of neighbors */
} VERTEX_STATS;

VERTEX_STATS *world_collect_stats(WORLD *a);
void fluxon_collect_stats(FLUXON *fl, VERTEX_STATS *st);

/**********************************************************************
 * Fluxon end-condition handlers and names
 */
extern struct F_B_NAMES {
  void (*func)(VERTEX *v);
  char name[80];
  char summary[80];
} F_B_NAMES[];

void fl_b_tied_force(VERTEX *v);
void fl_b_tied_inject(VERTEX *v);
void fl_b_open(VERTEX *v);
void fl_b_plasmoid(VERTEX *v);

#endif /* overall file include */

