/**********************************************************************
 * geomview.h -- I/O routine headers to talk with geomview 
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
#ifndef FLEM_GEOMVIEW
#define FLEM_GEOMVIEW 1

#include "data.h"
#include <stdio.h>

NUM *expand_vertex_rim(VERTEX *v, NUM r, int n_faces);
void OOGL_dump_fluxon(FILE *f, FLUXON *fl, NUM radius, int n_faces, char label, NUM *col, NUM *vcol);

void OOGL_vect(FILE *f, int count, NUM **x, NUM *color);
void OOGL_segment(FILE *f, NUM *x0, NUM *x1, NUM *color);

void OOGL_vertex(FILE *f, NUM *x, NUM radius, NUM *color,char redefine);
void OOGL_neighbor_lines(FILE *f, VERTEX *v, DUMBLIST *neighbors, NUM *c1, NUM *c2);
void OOGL_neighborhood_cell(FILE *f, VERTEX *v, NUM *color,NUM *c2);
void OOGL_draw_force(FILE *f, VERTEX *v, NUM *color);


void OOGL_label(FILE *f, char *string, NUM x[3], NUM dir[3], NUM normal[3], NUM height, char *align, NUM *color);

#endif
