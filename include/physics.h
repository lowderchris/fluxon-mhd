/**********************************************************************
 * physics.h -- physics routine headers for FLUX
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
 * This file is part of FLUX 2.0 (31-Oct-2007).
 */
 
#ifndef FLUX_PHYSICS
#define FLUX_PHYSICS 1

#include "data.h"
#include "geometry.h"
#include "model.h"

#include <stdio.h>

void b_eqa(VERTEX *V, HULL_VERTEX *verts); /* calculate B field at vertex */
void b_simple(VERTEX *V, HULL_VERTEX *verts); /* B field at vertex */

void e_simple(VERTEX *V, HULL_VERTEX *verts); /* energy, breaks if hull is open */
void e_simple2(VERTEX *V, HULL_VERTEX *verts); /* energy assoc. with vertex */
void e_open(VERTEX *V, HULL_VERTEX *verts); /* energy assoc. with vertex */
void e_eqa(VERTEX *V, HULL_VERTEX *verts);

void f_curvature(VERTEX *V, HULL_VERTEX *verts);     /* deprecated ( f/|B| ) */
void f_curvature2(VERTEX *V, HULL_VERTEX *verts);     /* deprecated ( f/|B| ) */
void f_curvature3(VERTEX *V, HULL_VERTEX *verts);   
void f_pressure_equi(VERTEX *V, HULL_VERTEX *verts); /* deprecated ( f/|B| ) */
void f_pressure_equi2(VERTEX *V, HULL_VERTEX *verts); /* deprecated ( f/|B| ) */
void f_pressure_equi2a(VERTEX *V, HULL_VERTEX *verts); /* deprecated ( f/|B| ) */
void f_pressure_equi2b(VERTEX *V, HULL_VERTEX *verts); /* deprecated ( f/|B| ) */

void f_p_eqa_radial(VERTEX *V, HULL_VERTEX *verts); /* ang. equipart. pressure */
void f_curv_hm(VERTEX *V, HULL_VERTEX *verts); /* harmonic-mean curvature  */ 
void f_curv_m(VERTEX *V, HULL_VERTEX *verts);  /* normal mean curvature */ 

void f_vertex(VERTEX *V, HULL_VERTEX *verts); /* vertex pseudoforce */
void f_vertex2(VERTEX *V, HULL_VERTEX *verts); /* vertex pseudoforce */
void f_vertex3(VERTEX *V, HULL_VERTEX *verts); /* vertex pseudoforce */
void f_vertex4(VERTEX *V, HULL_VERTEX *verts); /* vertex pseudoforce */
void f_vertex5(VERTEX *V, HULL_VERTEX *verts); /* vertex pseudoforce */
void f_vert(VERTEX *V, HULL_VERTEX *verts);   /* vertex pseudoforce */


/* F_CONV_TABLE is an array that associates function names with 
 * jumptable entries.  The actual array is defined in physics.c.
 */
extern struct FLUX_FORCES {
  char name[80];
  char summary[80];
  void ((*(func))());
} FLUX_FORCES[];

void *force_str_to_ptr(char *s);
char *force_ptr_to_str(void *f);

/* R_COND_TABLE is an array that associates function names with jumptable 8008ies

 * entries for reconnection conditions.  The actual array is defined in 
 * physics.c.
 */

extern struct FLUX_RECON {
  char name[80];
  char summary[80];
  RC_FUNC *func;
  char pnames[80];
} FLUX_RECON[];

void *recon_str_to_ptr(char *s);
char *recon_ptr_to_str(void *f);

RC_FUNC rc_a_ad2;
RC_FUNC rc_a_ad2_h;

#endif
