/* 
 * Concentration.xs - glue code for flux concentrations
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) 2004-2007 Craig DeForest.  You may distribute this
 * file under the terms of the Gnu Public License (GPL), version 2.
 * You should have received a copy of the GPL with this file.
 * If not, you may retrieve it from "http://www.gnu.org".
 *
 * This file is part of FLUX 2.2 (22-Nov-2008).
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
#include <flux/fluxperl.h>

#include <stdio.h>
#include <string.h>


#include "pdl.h"
#include "pdlcore.h"

static FluxCore* FLUX; /* FLUX core functions (run-time linking) */
static SV *FluxCoreSV;

static Core* PDL; /* PDL core functions (run-time linking)     */
static SV* CoreSV;/* gets perl var holding the core structures */

/******************************
 * Some helper routines used internally to this file...
/** fluxons_helper: callback to stuff a fluxon label into a list **/
static AV *fluxons_av;
static long fluxons_helper(void *tree, int lab_of, int ln_of, int depth) {
 SV *sv = newSViv(((FLUXON *)tree)->label);
 av_store(fluxons_av, av_len(fluxons_av)+1, sv) || svREFCNT_dec(sv);
 return 0;
}

/**********************************************************************
 **********************************************************************
 ***
 *** XS definitions follow...
 ***/
MODULE = Flux::Concentration    PACKAGE = Flux::Concentration

AV *
_fluxon_ids(wsv)
 SV *wsv
PREINIT:
 FLUX_CONCENTRATION *fc;
/**********************************************************************
 * _fluxon_ids - generate a perl list of fluxon IDs associated with this
 * world 
 * 
 * Hands back an array ref instead of a list, so I don't have to hassle
 * with the perl stack...
 */
CODE:
 fc = SvConc(wsv,"Flux::Concentration::_fluxon_ids",1);
 RETVAL = fluxons_av = newAV();  /* fluxons_av is static, at top */
 av_clear(RETVAL);
 if(fc->lines && fc->lines->fc0 == fc) {
	av_extend(RETVAL, (fc->lines->start_links).n);
 	FLUX->tree_walker(fc->lines,fl_lab_of,fl_start_ln_of,fluxons_helper,0);
 } else if(fc->lines && fc->lines->fc1 == fc) {
	av_extend(RETVAL, (fc->lines->end_links).n);
	FLUX->tree_walker(fc->lines,fl_lab_of,fl_end_ln_of,fluxons_helper,0);
 } else if(fc->lines) {
	croak("Flux::Concentration::_fluxon_ids: concentration seems screwed up!\n");
 }
OUTPUT:
 RETVAL

void
delete(fcsv)
 SV *fcsv
PREINIT:
 FLUX_CONCENTRATION *fc;
CODE:
 fc = SvConc(fcsv,"Flux::Concentration::delete",1);
 if(fc->label < 0) {
   croak("Flux::Concentration::delete: I refuse to delete a negative numbered concentration, that is too dangerous.");
 }
 FLUX->delete_flux_concentration(fc);

IV
cancel(fc0sv, fc1sv)
 SV *fc0sv
 SV *fc1sv
PREINIT:
 FLUX_CONCENTRATION *fc0, *fc1;
/**********************************************************************
 * cancel - Perl hook for fc_cancel in model.c
 */
CODE:
 fc0 = SvConc(fc0sv, "Flux::Concentration::cancel (0)",1);
 fc1 = SvConc(fc1sv, "Flux::Concentration::cancel (1)",1);
 RETVAL = FLUX->fc_cancel(fc0, fc1);
 if(RETVAL)
   croak("Flux::Concentration::cancel - fc_cancel threw an error!\n");
OUTPUT:
 RETVAL


BOOT:
/**********************************************************************
 **********************************************************************
 **** bootstrap code -- load-time dynamic linking to pre-loaded PDL
 **** modules and core functions.   **/
 require_pv("PDL/Core.pm");
 CoreSV = perl_get_sv("PDL::SHARE",FALSE);
 if(CoreSV==NULL)     Perl_croak(aTHX_ "Can't load PDL::Core module (required by Flux::Concentration)");

 PDL = INT2PTR(Core*, SvIV( CoreSV ));  /* Core* value */
 if (PDL->Version != PDL_CORE_VERSION)
    Perl_croak(aTHX_ "Flux::Concentration needs to be recompiled against the newly installed PDL");

 
 require_pv("Flux/Core.pm");
 FluxCoreSV = perl_get_sv("Flux::Core::FLUX",FALSE);
 if(FluxCoreSV == NULL)      Perl_croak(aTHX_ "Can't load Flux::Core module (required by Flux::Concentration)");
 FLUX = INT2PTR(FluxCore*, SvIV(FluxCoreSV));
 if(FLUX->CoreVersion != FLUX_CORE_VERSION) {
	printf("FLUX->CoreVersion is %ld; FLUX_%s is %d\n",FLUX->CoreVersion,"CORE_VERSION",FLUX_CORE_VERSION);
	Perl_croak(aTHX_ "Flux needs to be recompiled against the newly installed FLUX libraries");
}
