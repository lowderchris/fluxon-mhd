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







