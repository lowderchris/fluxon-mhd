/**********************************************************************
 * physics.h -- physics routine headers for FLEM
 * 
 * Craig DeForest, 14-Mar-2001
 */
 
#ifndef FLUX_PHYSICS
#define FLUX_PHYSICS 1

#include "data.h"
#include "geometry.h"
#include "model.h"

#include <stdio.h>

void f_curvature(VERTEX *V, HULL_VERTEX *verts);
void f_pressure_equi(VERTEX *V, HULL_VERTEX *verts);
void f_vertex(VERTEX *V, HULL_VERTEX *verts);

/* F_CONV_TABLE is an array that associates function names with 
 * jumptable entries.  The actual array is defined in physics.c.
 */
extern struct FLUX_FORCES {
  char name[80];
  char summary[80];
  void ((*(func))());
} FLUX_FORCES[];

#endif
