/* Vertex.xs - glue code for the Flux::Vertex object
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

MODULE = Flux::Vertex      PACKAGE = Flux::Vertex

char *
_stringify(vrt)
 SV *vrt
PREINIT:
 VERTEX *v;
 char str[BUFSIZ];
/**********************************************************************
 * _stringify - generate a summary string about a fluxon
 */
CODE: 
  v = SvVertex(vrt,"Flux::Vertex::_stringify");
  if(v->line) {
	sprintf(str,"vertex %5d (fl %5d): xyz=%5.3g,%5.3g,%5.3g, n=%5d,p=%5d, out/in:%2d/%2d\n",
	v->label,
	v->line->label,
	v->x ? v->x[0]:-1e64,
	v->x ? v->x[1]:-1e64,
	v->x ? v->x[2]:-1e64,
	v->next ? v->next->label : 0,
	v->prev ? v->prev->label : 0, 
	v->neighbors.n,
        v->nearby.n
        );
  } else {
	sprintf(str,"vertex %5d (IMAGE; xyz may be invalid): xyz=%5.3g,%5.3g,%5.3g\n",v->label, v->x?v->x[0]:-1e64, v->x?v->x[1]:-1e64, v->x?v->x[2]:-1e64);
  }

  RETVAL = str;
OUTPUT:
  RETVAL

SV *
next(vrt)
 SV *vrt
PREINIT:
 VERTEX *v;
 SV *sv;
/**********************************************************************
 * next - hop to the next vertex on the fluxon
 */
CODE:
 v = SvVertex(vrt,"Flux::Vertex::next");
 if(v->next) {	
   sv = newSViv((IV)(v->next));
   RETVAL = newRV_noinc(sv);
   (void)sv_bless(RETVAL,gv_stashpv("Flux::Vertex",TRUE));
 } else {
  RETVAL = &PL_sv_undef;
 }
OUTPUT:
 RETVAL

SV *
prev(vrt)
 SV *vrt
PREINIT:
 VERTEX *v;
 SV *sv;
/**********************************************************************
 * prev - hop to the previous vertex on the fluxon
 */
CODE:
 v = SvVertex(vrt,"Flux::Vertex::prev");
 if(v->prev) {
   sv = newSViv((IV)(v->prev));
   RETVAL = newRV_noinc(sv);
   (void)sv_bless(RETVAL,gv_stashpv("Flux::Vertex",TRUE));
 } else {
  RETVAL = &PL_sv_undef;
 }
OUTPUT:
 RETVAL

AV *
_adjacent(svrt,nearby)
SV *svrt
IV nearby
PREINIT:
 VERTEX *v;
 DUMBLIST *dl;
 int i;
 SV *sv, *rv;
/**********************************************************************
 * _adjacent - return a perl list of the neighbors or nearby list for a
 * vertex. 
 */
CODE:
 v = SvVertex(svrt,"Flux::Vertex::_adjacent");
 RETVAL = newAV(); /* initialize array */
 
 if(nearby)
  dl = &(v->nearby);
 else 
  dl = &(v->neighbors);
 
 av_clear(RETVAL);
 av_extend(RETVAL,dl->n);
 for(i=0; i<dl->n; i++) {
   VERTEX *v = (VERTEX *)(dl->stuff[i]);
   sv = newSViv((IV)v);
   rv = newRV_noinc(sv);
   (void)sv_bless(rv,gv_stashpv("Flux::Vertex",TRUE));
   if( !( av_store(RETVAL,i,rv) ) ) {
	svREFCNT_dec(rv);
	fprintf(stderr,"Warning: problems with array in _adjacent...\n");
   }
 }
OUTPUT:
 RETVAL

 
  