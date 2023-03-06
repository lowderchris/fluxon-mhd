/* Vertex.xs - glue code for the Flux::Vertex object
 * in perl.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) 2004-2008 Craig DeForest.  You may distribute this
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
 *  add_vertex_after   Constructs a VERTEX and adds it after the specified one.
 *  hull               Returns the projected hull of the vertex as a 7xN PDL
 *  photohull          Returns the real hull of the photospheric vertex as a 3xN PDL
 *  projmatrix         Returns the projection matrix used for hull, as a 3x3 PDL
 *  proj_neighbors     Calls vertex_update_neighbors.
 *  reconnect          Forces reconnection with another vertex
 *  x                  Returns the coordinates of the vertex - obviated by tied-hash
 * 
 * This is part of the FLUX 2.2 release (22-Nov-2008).
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


static FluxCore* FLUX; /* FLUX core functions (run-time linking) */
static SV *FluxCoreSV;


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
  v = SvVertex(vrt,"Flux::Vertex::_stringify",1);

  if(v->line) {
	sprintf(str,"vertex %5ld (fl %5ld): xyz=%7.3g,%7.3g,%7.3g, |B|=%7.3g, n=%5ld,p=%5ld, out/in:%2d/%2d\n",
	v->label,
	v->line->label,
	v->x[0],
	v->x[1],
	v->x[2],
	v->b_mag,
	v->next ? v->next->label : 0,
	v->prev ? v->prev->label : 0, 
	v->neighbors.n,
        v->nearby.n
        );
  } else {
	sprintf(str,"vertex %5ld (IMAGE; xyz may be invalid): xyz=%7.3g,%7.3g,%7.3g\n",v->label, v->x[0], v->x[1], v->x[2]);
  }

  RETVAL = str;
OUTPUT:
  RETVAL

void
delete(vrt)
 SV *vrt
PREINIT:
 VERTEX *v;
/**********************************************************************
 * delete - terminate with extreme prejudice
 */
CODE:
 v = SvVertex(vrt,"Flux::Vertex::delete",1);
 FLUX->delete_vertex(v);
 

long 
id(vrt)
 SV *vrt
PREINIT:
 VERTEX *v;
 SV *dv;
/**********************************************************************
 * id - return the longint ID as a perl scalar 
 * Kind of lame, as we could just pull the id out of the internal perl representation,
 * but calling SvVertex forces a consistency check.
 */
CODE:
 v = SvVertex(vrt,"Flux::Vertex::id",1);
 RETVAL = v->label;
OUTPUT:
 RETVAL

SV *
fluxon(vrt)
 SV *vrt
PREINIT:
 FLUXON *f;
 VERTEX *v;
 WORLD *w;
/**********************************************************************
 * fluxon - return the containing fluxon, as a Flux::Fluxon object
 */
CODE:
  v = SvVertex(vrt,"Flux::Vertex::fluxon",1);
  w = SvWorld(vrt, "Flux::Vertex::fluxon",1);
	
  RETVAL = (v->line ? 
	FLUX->new_sv_from_ptr(w, FT_FLUXON, v->line->label) :
	&PL_sv_undef
	);
OUTPUT:
 RETVAL

SV *
next(vrt)
 SV *vrt
PREINIT:
 VERTEX *v;
 WORLD *w;
/**********************************************************************
 * next - hop to the next vertex on the fluxon
 */
CODE:
 v = SvVertex(vrt,"Flux::Vertex::next",1);
 w = SvWorld(vrt, "Flux::Vertex::next",1);
 RETVAL = 
	(v->next ?
	FLUX->new_sv_from_ptr(w, FT_VERTEX, v->next->label) : 
	&PL_sv_undef
	);
OUTPUT:
 RETVAL

SV *
prev(vrt)
 SV *vrt
PREINIT:
 VERTEX *v;
 WORLD *w;
/**********************************************************************
 * prev - hop to the previous vertex on the fluxon
 */
CODE:
 v = SvVertex(vrt,"Flux::Vertex::prev",1);
 w = SvWorld(vrt,"Flux::Vertex::prev",1);
 RETVAL = 
	(v->prev ?
	FLUX->new_sv_from_ptr(w, FT_VERTEX, v->prev->label) :
	&PL_sv_undef
	);
OUTPUT:
 RETVAL

AV *
_adjacent(svrt,nearby)
SV *svrt
IV nearby
PREINIT:
 VERTEX *v;
 WORLD *w;
 AV *av;
 DUMBLIST *dl;
 int i;
 SV *sv;
/**********************************************************************
 * _adjacent - return a perl list of the neighbors or nearby list for a
 * vertex. 
 */
CODE:
 v = SvVertex(svrt,"Flux::Vertex::_adjacent",1);
 w = SvWorld(svrt, "Flux::Vertex::_adjacent",1);

 if(nearby)
  dl = &(v->nearby);
 else 
  dl = &(v->neighbors);
 

 av = newAV(); /* initialize array */
 av_clear(av);
 av_extend(av, dl->n);
 for(i=0; i<dl->n; i++) {
   VERTEX *v = (VERTEX *)(dl->stuff[i]);
   sv = FLUX->new_sv_from_ptr(w, FT_VERTEX, v->label);
   av_store(av, i, sv);
 }
 RETVAL = av;
OUTPUT:
 RETVAL

SV *
add_vertex_after(vsv, locsv)
SV *vsv
SV *locsv
PREINIT:
 VERTEX *v;
 VERTEX *nv;
 pdl *loc;
CODE:
 v = SvVertex(vsv,"Flux::Vertex::add_vertex_after",1);
	
 loc = PDL->SvPDLV(locsv);
 if(!loc || loc->ndims<1 || loc->dims[0] != 3) 
	croak("Flux::Vertex::add_vertex_after- requires a 3-PDL location");
 PDL->converttype(&loc, PDL_D, 1);
 PDL->make_physical(loc);
 nv = FLUX->new_vertex(0, 
	((PDL_Double *)(loc->data))[0],
	((PDL_Double *)(loc->data))[1],
	((PDL_Double *)(loc->data))[2],
	v->line
	);
 FLUX->add_vertex_after(v->line, v, nv);
 RETVAL = FLUX->new_sv_from_ptr(v->line->fc0->world, FT_VERTEX, nv->label);
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
 PDL_Indx dims[2];
 int i;
 PDL_Double *d;
CODE:
 /*********************************************
  * hull - return the Voronoi hull of a vertex
  */
 v = SvVertex(svrt,"Flux::Vertex::hull",1);
 hull_verts = FLUX->vertex_update_neighbors(v, global);
 
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
photohull(svrt,global=0)
SV *svrt
char global
PREINIT:
 VERTEX *v;
 pdl *p;
 HULL_VERTEX *hull_verts;
 PDL_Indx dims[2];
 int i;
 int n = 0;
 PDL_Double *d;
CODE:
 /*********************************************
  * hull - return the Voronoi hull of a vertex 
  * on the photosphere
  */
 v = SvVertex(svrt,"Flux::Vertex::hull",1); /**/
 hull_verts = FLUX->photosphere_vertex_update_neighbors(v, global, &n);
 
 /* Now stuff the hull vertices into the first PDL. */
 p = PDL->create(PDL_PERM);
 dims[0] = 4;
 dims[1] = n;
 PDL->setdims(p,dims,2);	
 p->datatype = PDL_D;
 PDL->allocdata(p);
 PDL->make_physical(p);
 d =  p->data;

 for(i=0;i<n; i++) {

   *(d++) = hull_verts[i].p[0];
   *(d++) = hull_verts[i].p[1];
   *(d++) = hull_verts[i].p[2];
   *(d++) = (double)(hull_verts[i].open);
 }
 RETVAL = NEWSV(545,0); /* 545 is arbitrary tag */ /**/
 PDL->SetSV_PDL(RETVAL,p);
OUTPUT:
 RETVAL

SV *
projmatrix(svrt)
SV *svrt
PREINIT:
 VERTEX *v;
 pdl *p;
 PDL_Indx dims[2];
 double *d,*d2;
 NUM x0[3],x1[3];
 NUM mat[9];
CODE: 
 v = SvVertex(svrt,"Flux::Vertex::projmatrix",1);
 
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
	FLUX->projmatrix(mat,x0,x1);
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
 PDL_Indx dims[2];
 int i;
 PDL_Double *d;
CODE:
 v = SvVertex(svrt,"Flux::Vertex::proj_neighbors",1);

 if(global) {
  dl = FLUX->gather_neighbor_candidates(v,1);
 } else {
  FLUX->vertex_update_neighbors(v,1);
  dl = &(v->neighbors);
 }
 FLUX->project_n_fill(v,dl);
 
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

void
reconnect(svv1, svv2, svpassno)
SV *svv1
SV *svv2
SV *svpassno
PREINIT:
 VERTEX *v1, *v2;
 long passno;
CODE:
 v1 = SvVertex(svv1, "Flux::Vertex::reconnect",1);
 v2 = SvVertex(svv2, "Flux::vertex::reconnect",1);
 if(SvROK(svpassno)) 
	croak("Flux::Vertex::reconnect: Can't take a PDL or reference as a passno");
if(svpassno == &PL_sv_undef) 
	passno = ++(v1->line->fc0->world->passno);
else
	passno = SvIV(svpassno);
 if(v1->line->fc0->world != v2->line->fc0->world) { 
	fprintf(stderr,"Error: vertices belong to different worlds!");
	return;
 }
 FLUX->reconnect_vertices(v1, v2, passno);

SV *
x(svrt)
SV *svrt
PREINIT:
 pdl *p;
 VERTEX *v;
 PDL_Indx dims[1];
CODE:
 v = SvVertex(svrt,"Flux::Vertex::x",1);
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
 require_pv("PDL/Core.pm");
 CoreSV = get_sv("PDL::SHARE",FALSE);
 if(CoreSV==NULL)     Perl_croak(aTHX_ "Can't load PDL::Core module (required by Flux::Vertex)");

 PDL = INT2PTR(Core*, SvIV( CoreSV ));  /* Core* value */
 if (PDL->Version != PDL_CORE_VERSION)
    Perl_croak(aTHX_ "Flux::Vertex needs to be recompiled against the newly installed PDL");


 require_pv("Flux/Core.pm");
 FluxCoreSV = get_sv("Flux::Core::FLUX",FALSE);
 if(FluxCoreSV == NULL)      Perl_croak(aTHX_ "Can't load Flux::Core module (required by Flux::Vertex)");
 
 FLUX = INT2PTR(FluxCore*, SvIV(FluxCoreSV));
 if(FLUX->CoreVersion != FLUX_CORE_VERSION) {
	Perl_croak(aTHX_ "Flux needs to be recompiled against the newly installed FLUX libraries");
}
