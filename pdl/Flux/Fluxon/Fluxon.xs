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

#include "pdl.h"
#include "pdlcore.h"


static Core* PDL;  /* PDL core functions (run-time linking) */
static SV* CoreSV; /* gets perl var holding the core structures */

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

SV *
polyline(flx)
 SV *flx
PREINIT:
 FLUXON *f;
 VERTEX *v;
 PDL_Double *d;
 pdl *p;
 SV *psv;
 int i;
 PDL_Long dims[2];
/**********************************************************************
 * polyline - return a 3xn PDL containing the coordinates of each 
 * VERTEX in the fluxon.  Useful for rendering.
 */
CODE:
 f = SvFluxon(flx,"Flux::Fluxon::polyline");

 /* Create the PDL and allocate its data */
 dims[0] = 3;
 dims[1] = f->v_ct;
 p = PDL->create(PDL_PERM);
 PDL->setdims(p,dims,2);
 p->datatype = PDL_D;
 PDL->allocdata(p);
 PDL->make_physical(p);
 d = p->data;

 /* Loop along the vertices, adding coordinates as we go */
 for(i=0, v=f->start; 
     i<f->v_ct && v;
     v=v->next, i++) {
  *(d++) = v->x[0];
  *(d++) = v->x[1];
  *(d++) = v->x[2];
 }
 RETVAL = NEWSV(546,0); /* 546 is arbitrary tag */
 PDL->SetSV_PDL(RETVAL, p);
OUTPUT:
 RETVAL

SV *
bfield(flx)
 SV *flx
PREINIT:
 FLUXON *f;
 VERTEX *v;
 PDL_Double *d;
 pdl *p;
 SV *psv;
 int i;
 PDL_Long dims[2];
/**********************************************************************
 * bfield - return a 6xN PDL containing the coordinates of each VERTEX
 * in the fluxon, together with the B-field components at that
 * location.
 */
CODE:
 f = SvFluxon(flx,"Flux::Fluxon::bfield");
 /* Create the PDL and allocate its data */
 dims[0] = 6;
 dims[1] = f->v_ct;
 p = PDL->create(PDL_PERM);
 PDL->setdims(p,dims,2);
 p->datatype = PDL_D;
 PDL->allocdata(p);
 PDL->make_physical(p);
 d = p->data;

 /* Loop along the vertices, adding coordinates as we go */
 for(i=0, v=f->start; 
     i<f->v_ct && v;
     v=v->next, i++) {
  *(d++) = v->x[0];
  *(d++) = v->x[1];
  *(d++) = v->x[2];
  *(d++) = v->b_vec[0];
  *(d++) = v->b_vec[1];
  *(d++) = v->b_vec[2];
 }
 RETVAL = NEWSV(547,0); /* 546 is arbitrary tag */
 PDL->SetSV_PDL(RETVAL, p);
OUTPUT:
 RETVAL


SV *
dump_vecs(flx)
 SV *flx
PREINIT:
 FLUXON *f;
 VERTEX *v;
 PDL_Double *d;
 pdl *p;
 SV *psv;
 int i;
 PDL_Long dims[2];
/**********************************************************************
 * dump_vecs - return a 17xN PDL containing the coordinates of each VERTEX
 * in the fluxon, together with the B-field components at that
 * location and forces.
 */
CODE:
 f = SvFluxon(flx,"Flux::Fluxon::dump_vecs");
 /* Create the PDL and allocate its data */
 dims[0] = 17;
 dims[1] = f->v_ct;
 p = PDL->create(PDL_PERM);
 PDL->setdims(p,dims,2);
 p->datatype = PDL_D;
 PDL->allocdata(p);
 PDL->make_physical(p);
 d = p->data;

 /* Loop along the vertices, adding coordinates as we go */
 for(i=0, v=f->start; 
     i<f->v_ct && v;
     v=v->next, i++) {
  *(d++) = v->x[0];
  *(d++) = v->x[1];
  *(d++) = v->x[2];
  *(d++) = v->b_vec[0];
  *(d++) = v->b_vec[1];
  *(d++) = v->b_vec[2];
  *(d++) = v->f_s[0];
  *(d++) = v->f_s[1];
  *(d++) = v->f_s[2];
  *(d++) = v->f_v[0];
  *(d++) = v->f_v[1];
  *(d++) = v->f_v[2];
  *(d++) = v->f_s_tot;
  *(d++) = v->f_v_tot;
  *(d++) = v->r_s;
  *(d++) = v->r_v;
  *(d++) = v->r_cl;
 }
 RETVAL = NEWSV(547,0); /* 546 is arbitrary tag */
 PDL->SetSV_PDL(RETVAL, p);
OUTPUT:
 RETVAL


BOOT:
/**********************************************************************
 **********************************************************************
 **** bootstrap code -- load-time dynamic linking to pre-loaded PDL
 **** modules and core functions.   **/
 perl_require_pv("PDL::Core");
 CoreSV = perl_get_sv("PDL::SHARE",FALSE);
 if(CoreSV==NULL)     Perl_croak(aTHX_ "Can't load PDL::Core module (required by Flux::Fluxon)");

 PDL = INT2PTR(Core*, SvIV( CoreSV ));  /* Core* value */
 if (PDL->Version != PDL_CORE_VERSION)
    Perl_croak(aTHX_ "Flux::Fluxon needs to be recompiled against the newly installed PDL");

 