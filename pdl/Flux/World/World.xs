/* 
 * World.xs - glue code for the Flux::World object
 * in perl.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) 2004 Craig DeForest.  You may distribute this
 * file under the terms of the Gnu Public License (GPL), version 2.
 * You should have received a copy of the GPL with this file.
 * If not, you may retrieve it from "http://www.gnu.org".
 *
 * Codes coverd here:
 *   PERL INTERFACE  FLUX SUBROUTINE / FUNCTION
 *
 *   _read_world        read_world (io.c) (leading '_' for ref in World.pm)
 *   write_world       write_world (io.c)
 *
 *   fix_proximity     global_fix_proximity (model.c)
 *   fix_curvature     global_fix_curvature (model.c)
 *  
 *   update_neighbors  world_update_neighbors (model.c)
 *
 *   fluxons           ----  Return a perl list of fluxon IDs in the world
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
#include <flux/perl.h>

#include <stdio.h>
#include <string.h>

/******************************
 * Some helper routines that are used internally
 */

/** fluxons_helper: callback to stuff a fluxon label into a list **/
static AV *fluxons_av;
static long fluxons_helper(void *tree, int lab_of, int ln_of, int depth) {
 SV *sv = newSViv(((FLUXON *)tree)->label);
 if(  !(  av_store(fluxons_av, av_len(fluxons_av)+1, sv)  )   )
	svREFCNT_dec(sv);
 return 0;
}

/**********************************************************************
 **********************************************************************
 ***
 *** XS definitions follow...
 ***/

MODULE = Flux::World    PACKAGE = Flux::World

SV *
_read_world(s)
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
	    croak("Couldn't open file '%s' to read a Flux::World",s);

	w = read_world(f,(WORLD *)0);

	sv = newSViv((IV)(w));
	RETVAL = newRV_noinc(sv);
	(void)sv_bless(RETVAL,gv_stashpv("Flux::World",TRUE));
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
    croak("Couldn't open file '%s' to write a Flux::World",s);
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



AV *
_fluxon_ids(wsv)
 SV *wsv
PREINIT:
 WORLD *w;
/**********************************************************************
 * _fluxon_ids - generate a perl list of fluxon IDs associated with this
 * world 
 * 
 * Hands back an array ref instead of a list, so I don't have to hassle
 * with the perl stack...
 */
CODE:
 w = SvWorld(wsv,"Flux::World::fluxon_ids");
 RETVAL = fluxons_av = newAV();  /* fluxons_av is static, at top */
 av_clear(RETVAL);
 if(w->lines) {
	av_extend(RETVAL, (w->lines->all_links).n);
 	tree_walker(w->lines,fl_lab_of,fl_all_ln_of,fluxons_helper);
 }
OUTPUT:
 RETVAL

SV *
fluxon(wsv,id)
 SV *wsv
 IV id
PREINIT:
 WORLD *w;
 FLUXON *f;
 SV *sv;
/**********************************************************************
 * fluxon - generate and return a Flux::Fluxon object associated 
 * with the given id
 */
CODE:
  w = SvWorld(wsv,"Flux::World::fluxon");
  f = (FLUXON *)tree_find(w->lines,id,fl_lab_of,fl_all_ln_of);
  if(f) {
     sv = newSViv((IV)(f));
     RETVAL = newRV_noinc(sv);
     (void)sv_bless(RETVAL,gv_stashpv("Flux::Fluxon",TRUE));
  } else {
     RETVAL = &PL_sv_undef;
  }
OUTPUT:
  RETVAL


AV *
_forces(wsv)
 SV *wsv
PREINIT:
 WORLD *w;
 SV *sv;
 int i;
/**********************************************************************
 * forces
 * Retrieve the force list in the World, using strings and the 
 * conversion array defined at the top of physics.c
 */
CODE:
 w=SvWorld(wsv,"Flux::World::forces");

 av_clear(RETVAL = newAV());

 for(i=0;i<N_FORCE_FUNCS && w->f_funcs[i]; i++) {
   int j;
   for(j=0; 
       FLUX_FORCES[j].func && FLUX_FORCES[j].func != w->f_funcs[i];
       j++
       );
   if(FLUX_FORCES[j].func) {
     sv = newSVpv(FLUX_FORCES[j].name,strlen(FLUX_FORCES[j].name));
   } else {
     char s[80];
     sprintf(s,"0x%x",(unsigned long)(w->f_funcs[i]));
     sv = newSVpv(s,strlen(s));
   }
   av_store(RETVAL, av_len(RETVAL)+1, sv) || svREFCNT_dec(sv);
   }
OUTPUT:
 RETVAL

void
_set_force(wsv,where,what)
SV *wsv
IV where
char * what
PREINIT:
 WORLD *w;
 SV *sv;
 int j;
/**********************************************************************
 * _set_force
 * Try to interpret a string as a force name and set the 
 * force pointer to it.
 */
CODE:
 w = SvWorld(wsv,"Flux::World::_set_force");
 if(where >= N_FORCE_FUNCS-1)
   croak("force index %d not allowed 0<i<%d\n",where,N_FORCE_FUNCS-1);
 for(j=0;
     FLUX_FORCES[j].func && strcmp(FLUX_FORCES[j].name,what);
     j++)
       ;
     if(!FLUX_FORCES[j].func) {
       if(what[0]=='0' && what[1]=='x') {
	 unsigned long ul;
	 sscanf(what+2,"%x",&ul);
	 w->f_funcs[where] = ul;
       } else {
         croak("Unknown force function '%s'\n",what);
       }
     } else {
        w->f_funcs[where] = FLUX_FORCES[j].func;
     }
