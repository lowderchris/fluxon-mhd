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
 *  id                 Returns the label of the vertex - obviated by the tied hash interface
 *  fluxon             Returns the fluxon to which this vertex belongs - obviated by tied-hash
 *  next               Returns the next vertex on the fluxon - obviated by tied-hash
 *  prev               Returns the prevoius vertex on the fluxon - obviated by tied-hash
 *  _adjacent	       Returns a perl list containing either the neighbors or nearby dumblist.
 *  hull               Returnns the projected hull of the vertex as a 7xN PDL
 *  projmatrix         Returns the projection matrix used for hull, as a 3x3 PDL
 *  proj_neighbors     Calls vertex_update_neighbors.
 *  x                  Returns the coordinates of the vertex - obviated by tied-hash
 * 
 * This is Vertex.xs version 1.1 - part of the FLUX 1.1 release.
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

#include "pdl.h"
#include "pdlcore.h"

static Core* PDL; /* PDL core functions (run-time linking)     */
static SV* CoreSV;/* gets perl var holding the core structures */


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
	sprintf(str,"vertex %5d (fl %5d): xyz=%7.3g,%7.3g,%7.3g, |B|=%7.3g, n=%5d,p=%5d, out/in:%2d/%2d\n",
	v->label,
	v->line->label,
	v->x ? v->x[0]:-1e64,
	v->x ? v->x[1]:-1e64,
	v->x ? v->x[2]:-1e64,
	v->b_mag,
	v->next ? v->next->label : 0,
	v->prev ? v->prev->label : 0, 
	v->neighbors.n,
        v->nearby.n
        );
  } else {
	sprintf(str,"vertex %5d (IMAGE; xyz may be invalid): xyz=%7.3g,%7.3g,%7.3g\n",v->label, v->x?v->x[0]:-1e64, v->x?v->x[1]:-1e64, v->x?v->x[2]:-1e64);
  }

  RETVAL = str;
OUTPUT:
  RETVAL

long 
id(vrt)
 SV *vrt
PREINIT:
 VERTEX *v;
 SV *dv;
/**********************************************************************
 * id - return the longint ID as a perl scalar 
 */
CODE:
 v = SvVertex(vrt,"Flux::Vertex::id");
 RETVAL = v->label;
OUTPUT:
 RETVAL

SV *
fluxon(vrt)
 SV *vrt
PREINIT:
 FLUXON *f;
 VERTEX *v;
 SV *sv;
/**********************************************************************
 * fluxon - return the containing fluxon, as a Flux::Fluxon object
 */
CODE:
  v = SvVertex(vrt,"Flux::Vertex::fluxon");
  f = v->line;
  sv = newSViv((IV)(f));
  RETVAL = newRV_noinc(sv);
  (void)sv_bless(RETVAL,gv_stashpv("Flux::Fluxon",TRUE));
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

 if(nearby)
  dl = &(v->nearby);
 else 
  dl = &(v->neighbors);
 

 RETVAL = newAV(); /* initialize array */
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

 
SV *
hull(svrt,global=0)
SV *svrt
char global
PREINIT:
 VERTEX *v;
 pdl *p;
 HULL_VERTEX *hull_verts;
 PDL_Long dims[2];
 int i;
 PDL_Double *d;
CODE:
 v = SvVertex(svrt,"Flux::Vertex::hull");
 hull_verts = vertex_update_neighbors(v, global);
 
 /* Now stuff the hull vertices into the first PDL. */
 p = PDL->create(PDL_PERM);
 dims[0] = 7;
 dims[1] = v->neighbors.n;
 PDL->setdims(p,dims,2);	
 p->datatype = PDL_D;
 PDL->allocdata(p);
 PDL->make_physical(p);
 d =  p->data;

 for(i=0;i<v->neighbors.n; i++) {

   *(d++) = ((VERTEX *)(v->neighbors.stuff[i]))->scr[0];
   *(d++) = ((VERTEX *)(v->neighbors.stuff[i]))->scr[1];
   *(d++) = hull_verts[i].p[0];
   *(d++) = hull_verts[i].p[1];
   *(d++) = (double)(hull_verts[i].open);
   *(d++) = hull_verts[i].a_l;
   *(d++) = hull_verts[i].a_r;
 }
 RETVAL = NEWSV(545,0); /* 545 is arbitrary tag */
 PDL->SetSV_PDL(RETVAL,p);
OUTPUT:
 RETVAL

SV *
projmatrix(svrt)
SV *svrt
PREINIT:
 VERTEX *v;
 pdl *p;
 int dims[2];
 double *d,*d2;
 NUM x0[3],x1[3];
 NUM mat[9];
CODE: 
 v = SvVertex(svrt,"Flux::Vertex::projmatrix");
 
 if( ! v->next ) {
	RETVAL = &PL_sv_undef;
 } else {
	p = PDL->create(PDL_PERM);
	dims[0]=3;
	dims[1]=3;
	PDL->setdims(p,dims,2);
	p->datatype = PDL_D;
	PDL->allocdata(p);
	PDL->make_physical(p);
	x0[0] = v->x[0]; x0[1] = v->x[1]; x0[2] = v->x[2];
	x1[0] = v->next->x[0]; x1[1] = v->next->x[1]; x1[2] = v->next->x[2];
	projmatrix(mat,x0,x1);
	d2=p->data;
	d = mat;
	*(d2++)=*(d++);
	*(d2++)=*(d++);
	*(d2++)=*(d++);
	*(d2++)=*(d++);
	*(d2++)=*(d++);
	*(d2++)=*(d++);
	*(d2++)=*(d++);
	*(d2++)=*(d++);
	*(d2++)=*(d++);
	RETVAL = NEWSV(549,0); /* 549 is arbitrary tag */
	PDL->SetSV_PDL(RETVAL,p);
 }
OUTPUT:
 RETVAL

SV *
proj_neighbors(svrt,global=0)
SV *svrt
char global
PREINIT:
 VERTEX *v;
 pdl *p;
 DUMBLIST *dl;
 PDL_Long dims[2];
 int i;
 PDL_Double *d;
CODE:
 v = SvVertex(svrt,"Flux::Vertex::proj_neighbors");

 if(global) {
  dl = gather_neighbor_candidates(v,1);
 } else {
  vertex_update_neighbors(v,1);
  dl = &(v->neighbors);
 }
 project_n_fill(v,dl);
 
 /* Now stuff the dumblist of vertices into the return PDL. */
 p = PDL->create(PDL_PERM);
 dims[0] = 3;
 dims[1] = dl->n;
 PDL->setdims(p,dims,2);	
 p->datatype = PDL_D;
 PDL->allocdata(p);
 PDL->make_physical(p);
 d =  p->data;
 for(i=0;i<dl->n; i++) {
   *(d++) = ((VERTEX *)(dl->stuff[i]))->scr[0];
   *(d++) = ((VERTEX *)(dl->stuff[i]))->scr[1];
   *(d++) = ((VERTEX *)(dl->stuff[i]))->label;
 }
 RETVAL = NEWSV(545,0); /* 545 is arbitrary tag */
 PDL->SetSV_PDL(RETVAL,p);
OUTPUT:
 RETVAL

SV *
x(svrt)
SV *svrt
PREINIT:
 pdl *p;
 VERTEX *v;
 PDL_Long dims[1];
CODE:
 v = SvVertex(svrt,"Flux::Vertex::x");
 p = PDL->create(PDL_PERM);
 dims[0] = 3;
 PDL->setdims(p,dims,1);
 p->datatype = PDL_D;
 PDL->allocdata(p);
 ((PDL_Double *)p->data)[0] = v->x[0];
 ((PDL_Double *)p->data)[1] = v->x[1];
 ((PDL_Double *)p->data)[2] = v->x[2];
 RETVAL = NEWSV(550,0); /* 550 is arbitrary */
 PDL->SetSV_PDL(RETVAL,p);
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

