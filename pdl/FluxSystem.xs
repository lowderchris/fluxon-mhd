/* 
 * FluxSystem.xs - glue code for the FluxSystem object
 * in perl.
 *
 *
 * Codes coverd here:
 *   PERL INTERFACE  FLUX SUBROUTINE / FUNCTION
 *
 *   read_world        read_world (io.c)       
 *   write_world       write_world (io.c)
 *
 *   fix_proximity     global_fix_proximity (model.c)
 *   fix_curvature     global_fix_curvature (model.c)
 *  
 *   update_neighbors  world_update_neighbors (model.c)
 *
 */
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"


/* Includes for libflux.a - don't bother with geomview stuff */
#include <flux/data.h>
#include <flux/geometry.h>
#include <flux/io.h>
#include <flux/model.h>
#include <flux/physics.h>

#include <stdio.h>

#define SvWorld(wsv,name) (   ( SvROK(wsv) && sv_derived_from((wsv),"FluxSystem") )  ? (WORLD *)(SvIVx(SvRV(wsv))) : ( croak("%s requires a FluxSystem.\n",(name)), (WORLD *)0 ) )


MODULE = FluxSystem    PACKAGE = FluxSystem 

SV *
read_world(s)
 char *s
PREINIT:
 FILE *f;
 WORLD *w;
 SV *sv;
/**********************************************************************/
/* read_world glue */
/**/
CODE:
	if(  !(f = fopen(s,"r"))  ) 
	    croak("Couldn't open file '%s' to read a FluxSystem",s);

	w = read_world(f,(WORLD *)0);

	sv = newSViv((IV)(w));
	RETVAL = newRV_noinc(sv);
	(void)sv_bless(RETVAL,gv_stashpv("FluxSystem",TRUE));
	fclose(f);
OUTPUT:
	RETVAL


void
write_world(wsv,s)
  SV *wsv
  char *s
PREINIT:
  WORLD *w;
  FILE *f;
  SV *sv;
/**********************************************************************/
/* write_world glue */
/**/
CODE:
  w = SvWorld(wsv,"write_world");
  if( !(f=fopen(s,"w")) ) 
    croak("Couldn't open file '%s' to write a FluxSystem",s);
  fprint_world(f,w,"");
  fclose(f);


int
fix_proximity(wsv,alpha)
 SV *wsv
 NV alpha
PREINIT:
 WORLD *w;
/**********************************************************************/
/* fix_proximity glue */
/**/
CODE:
  w = SvWorld(wsv,"fix_proximity");
  RETVAL = global_fix_proximity(w,alpha);
OUTPUT:
  RETVAL

int
fix_curvature(wsv,alpha)
 SV *wsv
 NV alpha
PREINIT:
 WORLD *w;
/**********************************************************************/
/* fix_curvature glue */
/**/
CODE:
  w = SvWorld(wsv,"fix_curvature");
  RETVAL = global_fix_curvature(w,alpha);
OUTPUT:
  RETVAL

void
update_neighbors(wsv, globalflag)
 SV *wsv
 IV globalflag
PREINIT:
 WORLD *w;
/**********************************************************************
 * update_neighbors glue
 */
CODE:
  w = SvWorld(wsv,"update_neighbors");
  world_update_neighbors(w,globalflag);


IV
verbosity(wsv,verbosity=0)
 SV *wsv
 IV verbosity
PREINIT:
 WORLD *w;
 IV v;
/**********************************************************************
 * verbosity glue - set the verbosity of activity with this World
 */
CODE:
  w = SvWorld(wsv,"verbosity");
  if(items>1)
    w->verbosity = verbosity;
  RETVAL = w->verbosity;
OUTPUT:
  RETVAL

