/**********************************************************************
 * physics.h -- physics routine headers for FLUX
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Craig DeForest, 2004-2008
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
 * This file is part of FLUX 2.2 (21-Nov-2008).
 */
 
#ifndef FLUX_PHYSICS
#define FLUX_PHYSICS 1

#include "data.h"
#include "geometry.h"
#include "model.h"

#include <stdio.h>

typedef void PHYSICS_FUNC(VERTEX *V, HULL_VERTEX *verts, int segflag);

PHYSICS_FUNC b_eqa;       // calculate B field at vertex
PHYSICS_FUNC b_simple;    // B field at vertex

PHYSICS_FUNC e_simple;    // energy, breaks if hull is open
PHYSICS_FUNC e_simple2;   // energy, better at open hulls
PHYSICS_FUNC e_open;      // energy, best with open hulls (still not great)
PHYSICS_FUNC e_eqa;       // equal-flux-per-angle formula for energy

PHYSICS_FUNC f_curvature;  
PHYSICS_FUNC f_curvature2; 
PHYSICS_FUNC f_curvature3; 

PHYSICS_FUNC f_pressure_equi;
PHYSICS_FUNC f_pressure_equi2;
PHYSICS_FUNC f_pressure_equi2a;
PHYSICS_FUNC f_pressure_equi2b;

PHYSICS_FUNC f_p_eqa_radial;
PHYSICS_FUNC f_curv_hm;
PHYSICS_FUNC f_curv_m;

PHYSICS_FUNC f_vertex;
PHYSICS_FUNC f_vertex2;
PHYSICS_FUNC f_vertex3;
PHYSICS_FUNC f_vertex4;
PHYSICS_FUNC f_vertex5;
PHYSICS_FUNC f_vert;

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
RC_FUNC rc_a_ad2_h_ad2hmax;

#endif
