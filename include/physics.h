/**********************************************************************
 * physics.h -- physics routine headers for FLEM
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
 */
 
#ifndef FLUX_PHYSICS
#define FLUX_PHYSICS 1

#include "data.h"
#include "geometry.h"
#include "model.h"

#include <stdio.h>

void b_vec(VERTEX *V, HULL_VERTEX *verts); /* calculate B field at vertex */

void f_curvature(VERTEX *V, HULL_VERTEX *verts);     /* deprecated ( f/|B| ) */
void f_pressure_equi(VERTEX *V, HULL_VERTEX *verts); /* deprecated ( f/|B| ) */

void f_p_eqa(VERTEX *V, HULL_VERTEX *verts); /* ang. equipartition pressure */
void f_curv(VERTEX *V, HULL_VERTEX *verts);  /* curvature force */ 

void f_vertex(VERTEX *V, HULL_VERTEX *verts); /* vertex pseudoforce */


                                                   


/* F_CONV_TABLE is an array that associates function names with 
 * jumptable entries.  The actual array is defined in physics.c.
 */
extern struct FLUX_FORCES {
  char name[80];
  char summary[80];
  void ((*(func))());
} FLUX_FORCES[];

#endif
