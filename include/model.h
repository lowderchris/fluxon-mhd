/* Nitty-gritty modeling routines for flem -- routines that find
 * neighbors, calculate neighborhoods, and such.  These routines 
 * use the libraries in data.c, geometry.c, and io.c.
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

#ifndef FLEM_MODEL
#define FLEM_MODEL 1

#include "data.h"
#include "geometry.h"
#include "io.h"

/**********************************************************************
 * Neighborhood handling routines
 */

void world_update_neighbors(WORLD *a, char global);
NUM *world_update_mag(WORLD *a, char global, void ((**f_funcs)()));
void world_relax_step(WORLD *a, NUM t);

void fluxon_update_neighbors(FLUXON *fl, char global);
NUM *fluxon_update_mag(FLUXON *fl, char global, void ((**f_funcs)()), NUM *minmax);

void fluxon_relax_step(FLUXON *fl, NUM t);

HULL_VERTEX *vertex_update_neighbors(VERTEX *v, char global);

DUMBLIST *gather_neighbor_candidates(VERTEX *v,char global);

int winnow_cmp_1(void *a, void *b);
void winnow_neighbor_candidates(VERTEX *v, DUMBLIST *horde);

HULL_VERTEX *hull_neighbors(VERTEX *v, DUMBLIST *horde);

void voronoi_vertex(NUM P[2], NUM B[2], NUM A[2]);
int fix_proximity(VERTEX *V, NUM scale_thresh);
int fluxon_fix_proximity(FLUXON *F, NUM scale_thresh);
int global_fix_proximity(WORLD *w, NUM scale_thresh);

int fix_curvature(VERTEX *V, NUM curve_thresh);
int fluxon_fix_curvature(FLUXON *F, NUM curve_thresh);
int global_fix_curvature(WORLD *w, NUM curve_thresh);

#endif /* overall file include */

