/* Fluxon.xs - glue code for the Flux::Fluxon object
 * in perl.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) 2004 Craig DeForest.  You may distribute this
 * file under the terms of the Gnu Public License (GPL), version 2.
 * You should have received a copy of the GPL with this file.
 * If not, you may retrieve it from "http://www.gnu.org".
 *
 * Codes covered here:
 *  PERL INTERFACE     FLUX SUBROUTINE / FUNCTION
 *
 *  _stringify         <NA>
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

MODULE = Flux::Fluxon      PACKAGE = Flux::Fluxon

char *
_stringify(flx)
 SV *flx
PREINIT:
 FLUXON *f;
 char str[BUFSIZ];
/**********************************************************************
 * _stringify - generate a summary string about a fluxon
 */
CODE: 
  f = SvFluxon(flx,"Flux::Fluxon::_stringify");
  sprintf(str,"Fluxon %5.5d: start-fc %d, end-fc %d   %d vertices\n",f->label,f->fc0->label,f->fc1->label,f->v_ct);
  RETVAL = str;
OUTPUT:
  RETVAL

SV *
vertex(flx,vno)
 SV *flx
 IV vno
PREINIT:
 FLUXON *f;
 VERTEX *v;
 SV *sv;
 long i;
/**********************************************************************
 * vertex - given a vertex index location, return a Flux::Vertex object
 * pointing to it.
 */
CODE:
  f = SvFluxon(flx,"Flux::Fluxon::vertex");
  v = (VERTEX *)0;
  if(vno >= 0) 
    for(i=0, v=f->start; i<vno && v; i++, v=v->next)
      ;
  if(v) {
      sv = newSViv((IV)(v));
      RETVAL = newRV_noinc(sv);
      (void)sv_bless(RETVAL,gv_stashpv("Flux::Vertex",TRUE));
  } else {
      RETVAL = &PL_sv_undef;
  }
OUTPUT:
  RETVAL
     
AV *
_vertices(flx)
 SV *flx
PREINIT:
 FLUXON *f;
 VERTEX *v;
 SV *sv,*rv;
/**********************************************************************
 * _vertices - return a perl list of all the vertices in a particular 
 * fluxon. Returns a list ref so I don't have to hassle with the perl 
 * stack.
 */
CODE:
 f = SvFluxon(flx,"Flux::Fluxon::vertices");
 RETVAL = newAV();                 /* Initialize array */
 av_clear(RETVAL);
 av_extend(RETVAL,f->v_ct);      /* Pre-extend for efficiency */
 for(v=f->start;   v;   v=v->next) {  /* Loop over vertices */
   sv = newSViv((IV)(v));             
   rv = newRV_noinc(sv);
   (void)sv_bless(rv,gv_stashpv("Flux::Vertex",TRUE));
   if( !(  av_store(RETVAL, av_len(RETVAL)+1, rv) ) ) {
     svREFCNT_dec(rv); 
     v=0;
   }
 }
OUTPUT:
  RETVAL


    






