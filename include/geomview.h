/**********************************************************************
 * geomview.h -- I/O routine headers to talk with geomview 
 * (part of FLEM)
 *
 * Craig DeForest, this datestamp last modified 26-Jan-2001
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
