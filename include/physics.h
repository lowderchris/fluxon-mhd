/**********************************************************************
 * physics.h -- physics routine headers for FLEM
 * 
 * Craig DeForest, 14-Mar-2001
 */
 
#ifndef FLEM_PHYSICS
#define FLEM_PHYSICS 1

#include "data.h"
#include "geometry.h"
#include "model.h"

#include <stdio.h>

void f_curvature(VERTEX *V, HULL_VERTEX *verts);
void f_pressure_equi(VERTEX *V, HULL_VERTEX *verts);
void f_vertex(VERTEX *V, HULL_VERTEX *verts);

#endif
