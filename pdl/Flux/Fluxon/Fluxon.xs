/* Fluxon.xs - glue code for the Flux::Fluxon object
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
 *  _stringify          <NA>
 *  vertex		<NA> - returns the corresponding vertex counting from start to finish.
 *  polyline            returns the locations of all the vertices in the fluxon as a 3xN PDL
 *  bfield              returns locations and B-field values of all vertices as a 6xN PDL
 *  dump_vecs           dumps a 17xN PDL containing a bunch of stuff (see Fluxon.pm).
 *  _new 		interface to new_fluxon in data.c, with vertex population
 *
 *  mutual_helicity	<NA> - calculates the mutual helicity of two fluxons.  (not working yet).
 * 
 * This file is part of the FLUX 2.2 release (22-Nov-2008).
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

static Core* PDL;  /* PDL core functions (run-time linking) */
static SV* CoreSV; /* gets perl var holding the core structures */

MODULE = Flux::Fluxon      PACKAGE = Flux::Fluxon


void
auto_open(fsv)
 SV *fsv
PREINIT:
 FLUXON *f;
CODE:
 f = SvFluxon(fsv,"Flux::Fluxon::Vertex",1);
 FLUX->fluxon_auto_open(f);



char *
_stringify(flx)
 SV *flx
PREINIT:
 FLUXON *f;
 char str[BUFSIZ];
/**********************************************************************
 * _stringify - generate a summary string about a fluxon. 
 */
CODE: 
  f = SvFluxon(flx,"Flux::Fluxon::_stringify",1);
  if(f) {
	  sprintf(str,"Fluxon %5.5ld: start-fc %ld, end-fc %ld   %ld vertices\n",f->label,f->fc0->label,f->fc1->label,f->v_ct);
  } else {
	  long label = FLUX->SvLabel(flx, "Flux::Fluxon::_stringify", "Flux::Fluxon");
	  sprintf(str,"Fluxon %ld not found - stale Perl link?\n",label);
  }
  RETVAL = str;
OUTPUT:
  RETVAL

void
delete(flx)
 SV *flx
PREINIT:
 FLUXON *f;
CODE:
 f = SvFluxon(flx,"Flux::Fluxon::delete",1);
 FLUX->delete_fluxon(f);

SV *
vertex(flx,vno)
 SV *flx
 IV vno
PREINIT:
 FLUXON *f;
 VERTEX *v;
 long i;
/**********************************************************************
 * vertex - given a vertex index location, return a Flux::Vertex object
 * pointing to it.
 */
CODE:
  f = SvFluxon(flx,"Flux::Fluxon::vertex",1);
  v = (VERTEX *)0;
  if(vno >= 0) 
    for(i=0, v=f->start; i<vno && v; i++, v=v->next)
      ;
  RETVAL = ( v ? 
             FLUX->new_sv_from_ptr(f->fc0->world, FT_VERTEX, v->label) :
	     &PL_sv_undef
	  );
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
 PDL_Indx dims[2];
/**********************************************************************
 * polyline - return a 3xn PDL containing the coordinates of each 
 * VERTEX in the fluxon.  Useful for rendering.
 */
CODE:
 f = SvFluxon(flx,"Flux::Fluxon::polyline",1);

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
 PDL_Indx dims[2];
/**********************************************************************
 * bfield - return a 6xN PDL containing the coordinates of each VERTEX
 * in the fluxon, together with the B-field components at that
 * location.
 */
CODE:
 f = SvFluxon(flx,"Flux::Fluxon::bfield",1);
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
 PDL_Indx dims[2];
/**********************************************************************
 * dump_vecs - return a 17xN PDL containing the coordinates of each VERTEX
 * in the fluxon, together with the B-field components at that
 * location and forces.
 */
CODE:
 f = SvFluxon(flx,"Flux::Fluxon::dump_vecs",1);
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



SV *
_new(wsv, fc0sv, fc1sv, flux, labelsv, vertssv)
 SV *wsv
 SV *fc0sv
 SV *fc1sv
 NV flux
 SV *labelsv
 SV *vertssv
PREINIT:
 FLUXON *f;
 VERTEX *v;
 WORLD *w;
 FLUX_CONCENTRATION *fc0;
 FLUX_CONCENTRATION *fc1;
 pdl *verts;
 PDL_Double *data;
 long label;
 /******************************
  * _new - generates a new fluxon.  For more options use the 
  * standard constructor Flux::Fluxon::new (in Fluxon.pm).
  * Two vertices are automagically created at the endpoints,
  * whether you specify intermediate vertices or no.
  * 
  * Returns a perl structure pointing to the new fluxon.
  */
CODE:
  w   = SvWorld(wsv,  "Flux::Fluxon::_new - world",1);
  fc0 = SvConc (fc0sv,"Flux::Fluxon::_new - fc0",1);
  fc1 = SvConc (fc1sv,"Flux::Fluxon::_new - fc1",1);

  if(!vertssv || vertssv == &PL_sv_undef || !(*(SvPV_nolen(vertssv)))) {
 	verts = 0;
  } else {
 	verts = PDL->SvPDLV(vertssv);

 	if(!verts) 
 		croak("Flux::Fluxon::_new - couldn't understand verts argument");
 	if(verts->ndims<1 || verts->ndims > 2) 
 		croak("Flux::Fluxon::_new - 3-PDL or 3xn-PDL required for verts argument");
 	if(verts->dims[0] != 3) 
 		croak("Flux::Fluxon::_new - 0th dim of verts argument must have size 3");
 	if(verts->ndims==1) {
 		verts->ndims=2;
 		verts->dims[1]=1;
  		verts->dimincs[1]=3*verts->dimincs[0];
 	}
 	PDL->converttype(&verts, PDL_D, 1);
         PDL->make_physical(verts);
 }

 label = SvIV(labelsv);
 f = FLUX->new_fluxon(1.0, fc0, fc1, label, 0);

 f->start = FLUX->new_vertex( -(f->label*2),   fc0->x[0], fc0->x[1], fc0->x[2], f );
 f->end   = FLUX->new_vertex( -(f->label*2)+1, fc1->x[0], fc1->x[1], fc1->x[2], f );
 f->start->next = f->end;
 f->end->prev = f->start;
 f->v_ct = 2;
 if(verts) {
 	long vno;
 	VERTEX *vlast = f->start;
 	long st = verts->dimincs[0];
 
 	for(vno=0;vno<verts->dims[1];vno++) {
 		long of = verts->dimincs[1] * vno;
 		FLUX->add_vertex_after(f, vlast, FLUX->new_vertex( 0,
 								   ((PDL_Double *)(verts->data))[ of           ],
 								   ((PDL_Double *)(verts->data))[ of + st      ],
 								   ((PDL_Double *)(verts->data))[ of + st + st ],
 								   f
 								 )
 					);
		vlast = vlast->next;
 	}
 }
 RETVAL = FLUX->new_sv_from_ptr(w, FT_FLUXON, f->label);
OUTPUT:
	 RETVAL		


void
reattach(flsv,fcsv)
 SV *flsv
 SV *fcsv
PREINIT:
 /**********************************************************************
  * reattach
  * 
  * Reattaches a fluxon to a new flux concentration, without changing 
  * the vertex geometry.
  */
  FLUX_CONCENTRATION *fc, *fc0;
  FLUXON *fl;
CODE:
  fl = SvFluxon(flsv, "Flux::Fluxon::reattach",1);
  fc = SvConc(fcsv,"Flux::Fluxon::reattach",1);
  
  if(fc->flux > 0) {
    /* source term */
    if(FLUX->tree_top(fl, fl_start_ln_of) != fl->fc0->lines) {
       croak("Flux::Fluxon::reattach - arena is corrupted! I give up! (start tree root != fc0->lines");
    }
    fl->fc0->lines = FLUX->tree_unlink( fl, fl_lab_of, fl_start_ln_of );
    fc->lines      = FLUX->tree_binsert( fc->lines, fl, fl_lab_of, fl_start_ln_of );
    fl->fc0 = fc;  
  } else if(fc->flux < 0) {
    if(FLUX->tree_top(fl, fl_end_ln_of) != fl->fc1->lines) {
       croak("Flux::Fluxon::reattach - arena is correupted! I give up! (end tree root != fc1->lines");
    }
    fl->fc1->lines = FLUX->tree_unlink( fl, fl_lab_of, fl_end_ln_of );
    fc->lines      = FLUX->tree_binsert( fc->lines, fl, fl_lab_of, fl_end_ln_of )    ;
    fl->fc1 = fc;
  } else {
    croak("Flux::Fluxon::reattach - destination concentration has zero flux, not sure what to do.  I give up.");
  }
    

double
mutual_helicity(fsv1, fsv2)
 SV *fsv1
 SV *fsv2
PREINIT:
 /**********************************************************************
  * mutual_helicity 
  * 
  * Calculate the mutual helicity of two fluxons.
  * (Should this be implemented as a fluxon method?)
  * Currently EXPERIMENTAL - spring 2008
  */
 WORLD *w;
 FLUXON *f1;
 FLUXON *f2;
 VERTEX *v1;
 VERTEX *v2;
 FLUX_CONCENTRATION *fc1b;
 FLUX_CONCENTRATION *fc1e;
 FLUX_CONCENTRATION *fc2b;
 FLUX_CONCENTRATION *fc2e;
 char open1b;
 char open1e;
 char open2b;
 char open2e;
 VERTEX *v;
 NUM mat[9];
 static NUM rotzx[9] = {0, 0, 1,
			0, 1, 0,
			-1,0, 0};
 static NUM origin[3] = {0,0,0};
 NUM scr[9],rmat[9];
CODE:
 f1 = SvFluxon(fsv1,"Flux::World::mutual_helicity",1);
 f2 = SvFluxon(fsv2,"Flux::World::mutual_helicity",1);
 fc1b = f1->fc0;
 fc1e = f1->fc1;
 fc2b = f2->fc0;
 fc2e = f2->fc1;
 w = fc1b->world;

 /* If we are using a spherical photosphere, convert everything to radius and
  * oblique Mercator projection.  Otherwise, just copy it.
  */
 if(w->photosphere.type == PHOT_PLANE) {
   for(v=f1->start; v; v=v->next)
     FLUX->cp_3d(v->scr, v->x);
   for(v=f2->start; v; v=v->next)
     FLUX->cp_3d(v->scr, v->x);
 } else if(w->photosphere.type == PHOT_SPHERE) {
    /* Finnd a matrix that places the origin at the intersection 
     * of the lines between the footpoints.
     */
   NUM x[3];
   FLUX->sum_3d(x, f1->start->x, f1->end->x);
   FLUX->sum_3d(x, x, f2->start->x);
   FLUX->sum_3d(x, x, f2->end->x);
   FLUX->scale_3d(x,x,0.25);
   FLUX->projmatrix(scr, origin, x);
   //FLUX->mat_mult_3x3(rmat, scr, rotzx);

   /* Now rmat rotates the "X-cross" to the X axis, which 
    * will become the spherical coordinate origin 
    */

   /* Loop over both fluxons.  The kludgey dual loop just
      keeps us from having to declare another subroutine.
    */
   for(v=f1->start; v; v= (v->next ? 
			   v->next : 
			   (v->line==f1 ? f2->start : 0))) {
     NUM vec;
     NUM scr[3];
     NUM *merc = v->scr;
 
     /* Perform the rotation */
     //FLUX->mat_vmult_3d(vec, rmat, v->x);

     /* Convert to spherical coordinates */
     //scr[2] = FLUX->norm_3d(vec);
     //scr[1] = asin(vec[2]/r);
     //scr[0] = atan2(vec[1],vec[0]);
     
     /* Convert to Mercator coordinates */
     merc[0] = scr[0];
   }
 }
printf("Not working yet.\n");
 RETVAL = 0;
OUTPUT:
 RETVAL

  
BOOT:
/**********************************************************************
 **********************************************************************
 **** bootstrap code -- load-time dynamic linking to pre-loaded PDL
 **** modules and core functions.   **/
 require_pv("PDL/Core.pm");
 CoreSV = perl_get_sv("PDL::SHARE",FALSE);
 if(CoreSV==NULL)     Perl_croak(aTHX_ "Can't load PDL::Core module (required by Flux::Fluxon)");

 PDL = INT2PTR(Core*, SvIV( CoreSV ));  /* Core* value */
 if (PDL->Version != PDL_CORE_VERSION)
    Perl_croak(aTHX_ "Flux::Fluxon needs to be recompiled against the newly installed PDL");


 require_pv("Flux/Core.pm");
 FluxCoreSV = perl_get_sv("Flux::Core::FLUX",FALSE);
 if(FluxCoreSV == NULL)      Perl_croak(aTHX_ "Can't load Flux::Core module (required by Flux::Fluxon)");
 
 FLUX = INT2PTR(FluxCore*, SvIV(FluxCoreSV));
 if(FLUX->CoreVersion != FLUX_CORE_VERSION) {
	Perl_croak(aTHX_ "Flux needs to be recompiled against the newly installed FLUX libraries");
}

