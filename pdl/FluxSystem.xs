/* 
 * FluxSystem.xs - glue code for the FluxSystem object
 * in perl.
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


MODULE = FluxSystem    PACKAGE = FluxSystem 



SV *
read_world(s)
 char *s
PREINIT:
 FILE *f;
 WORLD *w;
 SV *sv;
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
CODE:
  if(!SvROK(wsv) || !sv_derived_from(wsv,"FluxSystem"))
    croak("FluxSystem::write_world needs a FluxSystem to write.");
  if( !(f=fopen(s,"w")) ) 
    croak("Couldn't open file '%s' to write a FluxSystem",s);
  w = (WORLD *)(SvIVx(SvRV(wsv)));
  fprint_world(f,w,"");
  fclose(f);


int
fix_proximity(wsv,alpha)
 SV *wsv
 NV alpha
PREINIT:
 WORLD *w;
CODE:
  if(!SvROK(wsv)||!sv_derived_from(wsv,"FluxSystem"))
    croak("FluxSystem::global_proximity needs a FluxSystem");
  w = (WORLD *)(SvIVx(SvRV(wsv)));
  RETVAL = global_fix_proximity(w,alpha);
OUTPUT:
  RETVAL

int
fix_curvature(wsv,alpha)
 SV *wsv
 NV alpha
PREINIT:
 WORLD *w;
CODE:
  if(!SvROK(wsv)||!sv_derived_from(wsv,"FluxSystem"))
    croak("FluxSystem::global_curvature needs a FluxSystem");
  w = (WORLD *)(SvIVx(SvRV(wsv)));
  RETVAL = global_fix_curvature(w,alpha);
OUTPUT:
  RETVAL

void
update_neighbors(wsv, globalflag)
 SV *wsv
 IV globalflag
PREINIT:
 WORLD *w;
CODE:
  if(!SvROK(wsv)||!sv_derived_from(wsv,"FluxSystem"))
    croak("FluxSystem::update_neighbors needs a FluxSystem");
  w = (WORLD *)(SvIVx(SvRV(wsv)));
  world_update_neighbors(w,globalflag);
