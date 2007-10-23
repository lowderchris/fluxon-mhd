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

static FluxCore* FLUX; /* FLUX core functions (run-time linking) */
static SV *FluxCoreSV;

static Core* PDL; /* PDL core functions (run-time linking)     */
static SV* CoreSV;/* gets perl var holding the core structures */

/******************************
 * helper routines...
 *
 * fieldptr is a kludge -- it is the counterpart of the numeric structure 
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
	case 26: return (void *)&(v->energy);           break;
	case 27: return (void *)&(v->plan_step[0]);     break;
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
		case 25: return (void *)&(f->fc1);              break;
	        case 26: return (void *)&(f->plasmoid);         break;
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
	        case 20: return (void *)&(w->refct);                    break;
	        case 21: return (void *)&(w->rc_funcs);                 break;
	        case 22: return (void *)&(w->max_angle);                break;
	        case 23: return (void *)&(w->mean_angle);               break;
	        case 24: return (void *)&(w->dtau);                     break;
	        case 25: return (void *)&(w->rel_step);                 break;
	        case 26: return (void *)&(w->coeffs[0]);                break;
	        case 27: return (void *)&(w->n_coeffs);                 break;
	        case 28: return (void *)&(w->maxn_coeffs);              break;
		case 29: return (void *)&(w->handle_skew);              break;
		case 30: return (void *)&(w->auto_open);                break;
	        case 31: return (void *)&(w->default_bound);            break;
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
	        case 13: return (void *)&(fc->bound);                   break;
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
 FLUXON **fp;
CODE:
	ptr = (void *)SvFluxPtr(sv,"_rfluxon","Flux");
	fp = (FLUXON **)fieldptr(ptr,typeno,fieldno);
	if(fp && *fp) {
		I32 foo;
		ENTER;
		SAVETMPS;
		PUSHMARK(SP);
		XPUSHs(sv_2mortal(newSVpv("Flux::Fluxon",0)));
		XPUSHs(sv_2mortal(newSViv((IV)(*fp))));
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
_rfluxonlist(sv, typeno, fieldno)
 SV *sv
 long typeno
 long fieldno
PREINIT:
 void *ptr;
 FLUXON **fp;
CODE:
	ptr = (void *)SvFluxPtr(sv,"_rfluxonlist","Flux");
	fp = (FLUXON **)fieldptr(ptr,typeno,fieldno);
	

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


SV *
file_versions()
PREINIT:
 char buf[1024];
 char *newline="\n";
 SV *sv;
CODE:
 buf[0] = '\0';
 strcpy(buf,FLUX->code_info_data);
 strcat(buf,newline);
 strcat(buf,FLUX->code_info_io);
 strcat(buf,newline);
 strcat(buf,FLUX->code_info_geometry);
 strcat(buf,newline);
 strcat(buf,FLUX->code_info_model);
 strcat(buf,newline);
 strcat(buf,FLUX->code_info_physics);
 strcat(buf,newline);
 RETVAL = newSVpv(buf,0);
OUTPUT:
	RETVAL

NV
atan2_oct(y,x)
 NV y
 NV x
CODE:
	RETVAL = atan2_oct((NUM)y,(NUM)x);
OUTPUT:
	RETVAL

NV
cross_2d(x1,y1,x2,y2)
 NV x1
 NV y1
 NV x2
 NV y2
PREINIT:
	NUM p1[2];
	NUM p2[2];
CODE:
	p1[0]=x1; p1[1] = y1;
	p2[0]=x2; p2[1] = y2;
	RETVAL = cross_2d(p1,p2);
OUTPUT:
	RETVAL
			

SV *
_rcoeffs(sv, typeno, fieldno)
 SV *sv
 long typeno
 long fieldno
PREINIT:
 AV *av;
 WORLD *world;
CODE:
 I32 i;
 world = (void *)SvFluxPtr(sv, "_rcoeffs", "Flux::World");

 av = newAV();
 av_clear(av);
 av_extend(av, world->maxn_coeffs);

 for (i=0;i<MAXNUMCOEFFS && i<world->n_coeffs && world->coeffs[i];i++) {
	SV *vert;
        vert = newSVnv(world->coeffs[i]);
	av_push( av, vert );
 }

 RETVAL = newRV_inc((SV *) av);
OUTPUT:
 RETVAL


SV *
_wcoeffs(sv, typeno, fieldno, val)
 SV *sv
 long typeno
 long fieldno
 SV *val 
PREINIT:
 WORLD *world;
CODE:
 I32 i;
 pdl *p;
 PDL_Double *d;

 world = (void *)SvFluxPtr(sv,"_wcoeffs","Flux::World");

 if(SvROK(val) && sv_derived_from(val,"PDL")) {
	 // printf("looks like a PDL...\n"); 
	p = PDL->SvPDLV(val);
 } else {
	 // printf("Got a non-PDL...\n"); 
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
	croak("_wcoeff: Big trouble with PDL::new");
	ret = POPs;
	SvREFCNT_inc(ret);
	PUTBACK;
	FREETMPS;
	LEAVE;
	p = PDL->SvPDLV(sv_2mortal(ret));
 }

 PDL->converttype(&p,PDL_D,1); // Promote to double
 PDL->make_physical(p);

 for (i=0;i<p->nvals;i++) world->coeffs[i] = ((double *)(p->data))[i];
 world->n_coeffs = p->nvals;
/* ARD - Bizarre array problems were caused by val being returned mortal -
 *       too low a ref count led to it being freed.
 */

 SvREFCNT_inc(val);
 RETVAL = val;
OUTPUT:
 RETVAL		


SV *
_rforces(sv, typeno, fieldno)
 SV *sv
 long typeno
 long fieldno
PREINIT:
 AV *av;
 SV *rv;
 WORLD *world;
 void **f_funcs;
CODE:
 I32 i;
 world = (void *)SvFluxPtr(sv,"_rforces", "Flux::World");
 f_funcs = fieldptr(world, typeno, fieldno);

 av = newAV();
 for(i=0;i<N_FORCE_FUNCS && f_funcs[i]; i++) {
   int j;
   for(j=0; 
       FLUX->FLUX_FORCES[j].func && FLUX->FLUX_FORCES[j].func != f_funcs[i];
       j++
       ) 
	;

   if(FLUX->FLUX_FORCES[j].func) {
//printf("Match: %s (j=%d)\n",FLUX->FLUX_FORCES[j].name,j);
     sv = newSVpv(FLUX->FLUX_FORCES[j].name,strlen(FLUX->FLUX_FORCES[j].name));
   } else {
     char s[80];
     sprintf(s,"UNKNOWN (0x%x)",(unsigned long)(f_funcs[i]));
     sv = newSVpv(s,strlen(s));
   }
   av_store(av, av_len(av)+1, sv) || svREFCNT_dec(sv);
 }
 rv = newRV_noinc((SV *)av);
 SvREFCNT_inc(rv);
 RETVAL = rv;
OUTPUT:
 RETVAL


SV *
_wforces(sv, typeno, fieldno, val)
 SV *sv
 long typeno
 long fieldno
 SV *val
PREINIT:
 WORLD *world;
 void **f_funcs;
 SV *val_sv;
 AV *av;
CODE:
 I32 i,j,l;
 void ((*(new_f_funcs[N_FORCE_FUNCS]))());
 char *what;

 world = (void *)SvFluxPtr(sv, "_wforces","Flux::World");
 f_funcs = fieldptr(world, typeno, fieldno);

 if( !SvROK(val) ) {
	croak("Flux::_wforces: requires an array ref");
 }

 val_sv = SvRV(val);
 if(SvTYPE(val_sv) != SVt_PVAV) {
	croak("Flux::_wforces: non-array ref found (requires array ref)");
 }
 av = (AV *)val_sv;
 
 l = av_len(av)+1;
 for(i=0;i<l;i++) {
   SV **svp = av_fetch(av,i,0);
   what = (svp ? SvPV_nolen(*svp) : "");
   for(j=0;
     FLUX->FLUX_FORCES[j].func && strcmp(FLUX->FLUX_FORCES[j].name,what);
     j++)
       ;

   if(!(FLUX->FLUX_FORCES[j].func)) {
	char buf[10240];
	sprintf(buf, "Unknown force function '%s' (must be one of: %s", what, FLUX->FLUX_FORCES[0].name);
	for(i=1; FLUX->FLUX_FORCES[i].func; i++) {
		strcat(buf, ", ");
		strcat(buf, FLUX->FLUX_FORCES[i].name);
	}
	strcat(buf,")\n");
        croak(buf);
   } else {
      printf("Match: %s at %d\n",FLUX->FLUX_FORCES[j].name,i);
      new_f_funcs[i] = FLUX->FLUX_FORCES[j].func;
   }
  }

  for(j=0;j<l;j++) 
    f_funcs[j] = new_f_funcs[j];
  f_funcs[l] = 0;
  SvREFCNT_inc(val);
  RETVAL = val;
OUTPUT:
  RETVAL  


SV *
_rbound(sv, typeno, fieldno)
 SV *sv
 long typeno
 long fieldno
PREINIT:
 AV *av;
 SV *rv;
 void *ptr;
 void **field;
CODE:
 I32 i;
 ptr = (void *)SvFluxPtr(sv,"_rbound", "Flux");
 field = fieldptr(ptr,typeno, fieldno);

 av = newAV();
 for(i=0;FLUX->F_B_NAMES[i].func && 
	( (void *)(FLUX->F_B_NAMES[i].func) != *field);
	i++)
	;
 if(((void *)FLUX->F_B_NAMES[i].func) == *field) {
     sv = newSVpv(FLUX->F_B_NAMES[i].name,
	 strlen(  FLUX->F_B_NAMES[i].name));
 } else {
     char s[80];
     sprintf(s,"UNKNOWN (0x%x)",(unsigned long)(*field));
     sv = newSVpv(s,strlen(s));
 }
 RETVAL = sv;
OUTPUT:
 RETVAL



SV *
_wbound(sv, typeno, fieldno, val)
 SV *sv
 long typeno
 long fieldno
 SV *val
PREINIT:
 void *ptr;
 void **field;
 SV *val_sv;
 char *what;
CODE:
 I32 j;
 I32 l;

 ptr = (void *)SvFluxPtr(sv, "_wbound","Flux");
 field = fieldptr(ptr, typeno, fieldno);

 what = SvPV_nolen(val);

 if(!*what) {
	*field = 0;
 } else {
	 for(j=0;
	     FLUX->F_B_NAMES[j].func &&
	     strcmp(FLUX->F_B_NAMES[j].name, what);
	     j++)
	//printf("Trying 0x%x (%s); value is %x\n",FLUX->F_B_NAMES[j].func,FLUX->F_B_NAMES[j].name,*field);
	       ;
	 if(!(FLUX->F_B_NAMES[j].func)) {
		char buf[10240];
		sprintf(buf, "Unknown fluxon end condition '%s' (must be one of: %s", what, FLUX->F_B_NAMES[0].name);
		for(j=1;FLUX->F_B_NAMES[j].func;j++) {
			strcat(buf, ", ");
			strcat(buf, FLUX->F_B_NAMES[j].name);
	        }
		strcat(buf,")\n");
		croak(buf);
	 } else {
	     *field = (void*)(FLUX->F_B_NAMES[j].func);
	 }
  }
  RETVAL = newSVsv(val);
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

 perl_require_pv("Flux::Core");
 FluxCoreSV = perl_get_sv("Flux::Core::FLUX",FALSE);
 if(FluxCoreSV == NULL)      Perl_croak(aTHX_ "Can't load Flux::Core module (required b Flux)");
 
 FLUX = INT2PTR(FluxCore*, SvIV(FluxCoreSV));
 if(FLUX->CoreVersion != FLUX_CORE_VERSION) {
	Perl_croak(aTHX_ "Flux needs to be recompiled against the newly installed FLUX libraries");
}

  
