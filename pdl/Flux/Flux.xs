/* 
 * Flux.xs - glue code for the Flux shell class in perl.
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
 * This is version 1.0 of Flux.xs 
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
#include <flux/fluxperl.h>

#include <stdio.h>
#include <string.h>


#include "pdl.h"
#include "pdlcore.h"

static Core* PDL; /* PDL core functions (run-time linking)     */
static SV* CoreSV;/* gets perl var holding the core structures */

/******************************
 * helper routines...
 *
 * fieldptr is a kludge -- it is the opposite of the numeric structure 
 * at the top of Flux.pm.  Given a void * and a type and field code,
 * it returns an appropriate void * pointer.
 *
 */


void *fieldptr(void *foo, long typeno, long fieldno) {
 switch(typeno) {
  case 2: {
    VERTEX *v = (VERTEX *)foo;
    switch(fieldno) {
	case 1:  return (void *)&(v->line);       	break;
	case 2:  return (void *)&(v->prev);       	break;
	case 3:  return (void *)&(v->next);       	break;
	case 4:  return (void *)&(v->x[0]);       	break;
	case 5:  return (void *)&(v->neighbors);  	break;
   	case 6:  return (void *)&(v->nearby);     	break;
	case 7:  return (void *)&(v->scr[0]);     	break;
	case 8:  return (void *)&(v->r);          	break;
	case 9:  return (void *)&(v->a);          	break;
	case 10: return (void *)&(v->b_vec[0]);   	break;
	case 11: return (void *)&(v->b_mag);      	break;
	case 12: return (void *)&(v->f_v[0]);     	break;
	case 13: return (void *)&(v->f_s[0]);     	break;
	case 14: return (void *)&(v->f_t[0]);     	break;
	case 15: return (void *)&(v->f_s_tot);    	break;
	case 16: return (void *)&(v->f_v_tot);    	break;
	case 17: return (void *)&(v->r_v);        	break;
	case 18: return (void *)&(v->r_s);        	break;
	case 19: return (void *)&(v->r_cl);       	break;
	case 20: return (void *)&(v->label);      	break;
	case 21: return (void *)&(v->world_links.sum);  break;
	case 22: return (void *)&(v->world_links.n);    break;
	case 23: return (void *)&(v->world_links.up);   break;
	case 24: return (void *)&(v->world_links.left); break;
	case 25: return (void *)&(v->world_links.right);break;
	default: fprintf(stderr,"Unknown type,field (%d,%d) in Flux::World::fieldptr!\n",
	                 typeno,fieldno);
	         return (void *)0;
	         break;
       }
	}
     break;
   case 3: {
	FLUXON *f = (FLUXON *)foo;
	switch(fieldno) {
		case  1: return (void *)&(f->flux);		break;
		case  2: return (void *)&(f->label);		break;
		case  3: return (void *)&(f->start);		break;
		case  4: return (void *)&(f->end);		break;
		case  5: return (void *)&(f->v_ct);		break;
		case  6: return (void *)&(f->all_links);	break;
		case  7: return (void *)&(f->all_links.sum);	break;
		case  8: return (void *)&(f->all_links.n);	break;
		case  9: return (void *)&(f->all_links.up);	break;
		case 10: return (void *)&(f->all_links.left);	break;
		case 11: return (void *)&(f->all_links.right);	break;
		case 12: return (void *)&(f->start_links);	break;
		case 13: return (void *)&(f->start_links.sum);	break;
		case 14: return (void *)&(f->start_links.n);	break;
		case 15: return (void *)&(f->start_links.up);	break;
		case 16: return (void *)&(f->start_links.left);	break;
		case 17: return (void *)&(f->start_links.right);break;
		case 18: return (void *)&(f->end_links);	break;
		case 19: return (void *)&(f->end_links.sum);	break;
		case 20: return (void *)&(f->end_links.n);	break;
		case 21: return (void *)&(f->end_links.up);	break;
		case 22: return (void *)&(f->end_links.left);	break;
		case 23: return (void *)&(f->end_links.right);	break;
		case 24: return (void *)&(f->fc0);              break;
		case 25: return (void *)&(f->fc1);
		default: fprintf(stderr,"Unknown type,field (%d,%d) in Flux::World::fieldptr!\n",
				typeno,fieldno);
			return (void *)0;
			break;
	}
       }
       break;
     case 4: {
	WORLD *w = (WORLD *)foo;
	switch(fieldno) {
		case  1: return (void *)&(w->frame_number);           	break;
		case  2: return (void *)&(w->state);          		break;
		case  3: return (void *)&(w->concentrations);		break;
		case  4: return (void *)&(w->lines);			break;
		case  5: return (void *)&(w->vertices);			break;
		case  6: return (void *)&(w->photosphere);		break;
		case  7: return (void *)&(w->image);			break;
		case  8: return (void *)&(w->image2);			break;
		case  9: return (void *)&(w->locale_radius);		break;
		case 10: return (void *)&(w->fc_ob);			break;
		case 11: return (void *)&(w->fc_oe);			break;
		case 12: return (void *)&(w->fc_pb);			break;
		case 13: return (void *)&(w->fc_pe);			break;
		case 14: return (void *)&(w->verbosity);		break;
		case 15: return (void *)&(w->f_funcs);			break;
		case 16: return (void *)&(w->step_scale.b_power);	break;
		case 17: return (void *)&(w->step_scale.d_power);	break;
		case 18: return (void *)&(w->step_scale.s_power);	break;
		case 19: return (void *)&(w->step_scale.ds_power);	break;
		default: fprintf(stderr,"Unknown type,field (%d,%d) in Flux::World::fieldptr!\n",
				typeno, fieldno);		
			 return (void *)0;
		         break;
	 }
	}
	break;
     case 5: {
	FLUX_CONCENTRATION *fc = (FLUX_CONCENTRATION *)foo;
	switch(fieldno) {
		case  1: return (void *)&(fc->world);              	break;
		case  2: return (void *)&(fc->flux);		   	break;
		case  3: return (void *)&(fc->label);              	break;
		case  4: return (void *)&(fc->lines); 	           	break;
		case  5: return (void *)&(fc->links);              	break;
		case  6: return (void *)&(fc->links.sum);		break;
		case  7: return (void *)&(fc->links.n);			break;
		case  8: return (void *)&(fc->links.up);		break;
		case  9: return (void *)&(fc->links.left);		break;
		case 10: return (void *)&(fc->links.right);		break;
		case 11: return (void *)&(fc->x[0]);			break;
		case 12: return (void *)&(fc->locale_radius);		break;
		default: fprintf(stderr,"Unknown type,field (%d,%d) in Flux::World::fieldptr!\n",
				typeno, fieldno);
			return (void *)0;
			break;
	 }
        }
	break;
     default: fprintf(stderr,"Unknown type number %d in Flux::World::fieldptr!\n",typeno);
           return (void *)0;
           break;
  }
}

/**********************************************************************
 **********************************************************************
 ***
 *** XS definitions follow...
 ***/


MODULE = Flux    PACKAGE = Flux

SV *
_rnum(sv,typeno,fieldno)
  SV *sv
  long typeno
  long fieldno
PREINIT:
  void *ptr;
  NUM *np;
CODE:
	ptr = (void *)SvFluxPtr(sv,"_rnum","Flux");
	np = fieldptr(ptr,typeno,fieldno);
	if(np)
		RETVAL = newSVnv(*np);	
	else
		RETVAL = &PL_sv_undef;
OUTPUT:
	RETVAL

NV
_wnum(sv,typeno,fieldno,val)
 SV *sv
 long typeno
 long fieldno
 NV val
PREINIT:
 void *ptr;
 NUM *np;
CODE:
	ptr = (void *)SvFluxPtr(sv,"_wnum","Flux");
	np = fieldptr(ptr,typeno,fieldno);
	if(np)
		*np = val;
	else
		croak("Flux::_wnum: fieldptr failed");
	RETVAL = val;
OUTPUT:
	RETVAL

SV *
_rlong(sv,typeno,fieldno)
 SV *sv
 long typeno
 long fieldno
PREINIT:
 void *ptr;
 long *lp;
CODE:

	ptr = (void *)SvFluxPtr(sv,"_rlong","Flux");

	lp = fieldptr(ptr,typeno,fieldno);

	if(lp)
		RETVAL = newSViv(*lp);
	else
		RETVAL = &PL_sv_undef;
OUTPUT:
	RETVAL

IV
_wlong(sv,typeno,fieldno,val)
 SV *sv
 long typeno
 long fieldno
 IV val
PREINIT:
 void *ptr;
 long *lp;
CODE:
	ptr = (void *)SvFluxPtr(sv,"_wlong","Flux");
	lp = fieldptr(ptr,typeno,fieldno);
	if(lp)
		*lp = val;
	else
		croak("Flux::_wlong: fieldptr failed");
	RETVAL = val;
OUTPUT:
	RETVAL

SV *
_rvec(sv,typeno,fieldno)
 SV *sv
 long typeno
 long fieldno
PREINIT:
 void *ptr;
 NUM *np;
CODE:
	ptr = (void *)SvFluxPtr(sv,"_rvec","Flux");
	np = (NUM *)fieldptr(ptr,typeno,fieldno);
	if(np) {	
		pdl *p;
		PDL_Double *d;
		SV *psv;
		PDL_Long dim=3;

		p = PDL->create(PDL_PERM);
		PDL->setdims(p,&dim,1);
		p->datatype = PDL_D;
		PDL->allocdata(p);
		PDL->make_physical(p);
		d = p->data;
		d[0] = np[0];
		d[1] = np[1];
		d[2] = np[2];
		RETVAL = NEWSV(551,0); /* 551 is arbitrary */
		PDL->SetSV_PDL(RETVAL,p);
	} else
		RETVAL = &PL_sv_undef;
OUTPUT:
	RETVAL


SV *
_wvec(sv,typeno,fieldno,val)
 SV *sv
 long typeno
 long fieldno
 SV *val
PREINIT:
 void *ptr;
 NUM *np;
CODE: 
	ptr = (void *)SvFluxPtr(sv,"_wvec","Flux");
	np = (NUM *)fieldptr(ptr,typeno,fieldno);
	if(np) {
		pdl *p;
		PDL_Double *d;
	
		if(SvROK(val) && sv_derived_from(val,"PDL")) {
			/* printf("looks like a PDL...\n"); */
			p = PDL->SvPDLV(val);
		} else {
			/* printf("Got a non-PDL...\n"); */
			/* Hard way - got a not-PDL. Dive into perlspace to make one. */
			I32 foo;	
			SV *ret;
			ENTER;
			SAVETMPS;
			PUSHMARK(SP);
			XPUSHs(sv_2mortal(newSVpv("PDL",0)));
			XPUSHs(val);
			PUTBACK;
			foo = call_pv("PDL::new",G_SCALAR);
			SPAGAIN;
			if(foo != 1)
				croak("_wvec: Big trouble with PDL::new");
			ret = POPs;
			SvREFCNT_inc(ret);
			PUTBACK;
			FREETMPS;
			LEAVE;
			p = PDL->SvPDLV(sv_2mortal(ret));
		}
			
		PDL->converttype(&p,PDL_D,1); // Promote to double
		PDL->make_physical(p);
		if(p->nvals == 3) {
	/* printf("3 vals\n"); */
			d = p->data;
			np[0] = d[0];
			np[1] = d[1];
			np[2] = d[2];
		} else if(p->nvals == 1) {
	/* printf("1 val\n"); */
			d = p->data;
			np[0] = np[1] = np[2] = d[0];
		}
	} else 
		croak("Flux::_wvec: fieldptr failed");
	
	RETVAL = newSVsv(val);
OUTPUT:
	RETVAL
		

SV *
_rvertex(sv, typeno, fieldno)
 SV *sv
 long typeno
 long fieldno
PREINIT:
 void *ptr;
 VERTEX **vp;
CODE:
	ptr = (void *)SvFluxPtr(sv,"_rvertex","Flux");
	vp = (VERTEX **)fieldptr(ptr,typeno,fieldno);
	if(vp && *vp) {
		I32 foo;
		ENTER;
		SAVETMPS;
		PUSHMARK(SP);
		XPUSHs(sv_2mortal(newSVpv("Flux::Vertex",0)));
		XPUSHs(sv_2mortal(newSViv((IV)(*vp))));
		PUTBACK;
		foo = call_pv("Flux::Vertex::new_from_ptr",G_SCALAR);
		SPAGAIN;
		
		if(foo==1)
			RETVAL = POPs;
		else
			croak("Big trouble in _rvertex!");
		
		SvREFCNT_inc(RETVAL);
		PUTBACK;
		FREETMPS;
		LEAVE;
	} else {
		RETVAL = &PL_sv_undef;
	}
OUTPUT:
  	RETVAL

SV *
_rfluxon(sv, typeno, fieldno)
 SV *sv
 long typeno
 long fieldno
PREINIT:
 void *ptr;
 VERTEX **vp;
CODE:
	ptr = (void *)SvFluxPtr(sv,"_rfluxon","Flux");
	vp = (VERTEX **)fieldptr(ptr,typeno,fieldno);
	if(vp && *vp) {
		I32 foo;
		ENTER;
		SAVETMPS;
		PUSHMARK(SP);
		XPUSHs(sv_2mortal(newSVpv("Flux::Fluxon",0)));
		XPUSHs(sv_2mortal(newSViv((IV)(*vp))));
		PUTBACK;
		foo = call_pv("Flux::Fluxon::new_from_ptr",G_SCALAR);
		SPAGAIN;
		
		if(foo==1)
			RETVAL = POPs;
		else
			croak("Big trouble in _rfluxon!");
		
		SvREFCNT_inc(RETVAL);
		PUTBACK;
		FREETMPS;
		LEAVE;
	} else {
		RETVAL = &PL_sv_undef;
	}
OUTPUT:
  	RETVAL

SV *
_rconcentration(sv, typeno, fieldno)
 SV *sv
 long typeno
 long fieldno
PREINIT:
 void *ptr;
 FLUX_CONCENTRATION **vp;
CODE:
	ptr = (void *)SvFluxPtr(sv,"_rconcentration","Flux");
	vp = (FLUX_CONCENTRATION **)fieldptr(ptr,typeno,fieldno);
	if(vp && *vp) {
		I32 foo;
		ENTER;
		SAVETMPS;
		PUSHMARK(SP);
		XPUSHs(sv_2mortal(newSVpv("Flux::Concentration",0)));
		XPUSHs(sv_2mortal(newSViv((IV)(*vp))));
		PUTBACK;
		foo = call_pv("Flux::Concentration::new_from_ptr",G_SCALAR);
		SPAGAIN;
		
		if(foo==1)
			RETVAL = POPs;
		else
			croak("Big trouble in _rconcentration!");
		
		SvREFCNT_inc(RETVAL);
		PUTBACK;
		FREETMPS;
		LEAVE;
	} else {
		RETVAL = &PL_sv_undef;
	}
OUTPUT:
  	RETVAL


SV *
_rdumblist(class, sv, typeno, fieldno)
 SV *class
 SV *sv
 long typeno
 long fieldno
PREINIT:
 void *ptr;
 DUMBLIST *dl;
 void **array;
 char *classpv;
 AV *av;
CODE:
 I32 i;
 ptr = (void *)SvFluxPtr(sv,"_rdumblist","Flux");
 dl = (DUMBLIST *)fieldptr(ptr,typeno,fieldno);
 array = dl->stuff;
 static char constructorbuf[100];

 av = newAV();         /* new perl list */
 av_clear(av);         /* Make sure it's empty */
 av_extend(av, dl->n); /* pre-extend for efficiency */

 classpv = SvPVX(class);

 strncpy(constructorbuf,classpv,80);
 strcat(constructorbuf,"::new_from_ptr");

 /* printf("rdumblist: n=%d\n",dl->n);*/
 for(i=0;i<dl->n;i++) {
	SV *vert;
	I32 foo;
	ENTER;
	SAVETMPS;

	PUSHMARK(SP);
	XPUSHs(class);
	XPUSHs(sv_2mortal(newSViv((IV)(array[i]))));
	PUTBACK;
	foo = call_pv(constructorbuf,G_SCALAR);
	SPAGAIN;
	
	if(foo==1)
	 vert = POPs;
   	else 	
	  croak("Big trouble in _rdumblist - giving up");
	
	SvREFCNT_inc(vert);

	PUTBACK;
	FREETMPS;
	LEAVE;
 /*	printf("storing %d...",i);*/
	av_push( av, vert );
 }
 RETVAL = newRV_inc((SV*)av);
 /* printf("done\n");*/
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
    Perl_croak(aTHX_ "Flux needs to be recompiled against the newly installed PDL");
