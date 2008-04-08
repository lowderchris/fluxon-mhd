/* 
 * Flux.xs - glue code for the Flux shell class in perl.
 *
 * The routines in Flux.xs focus on making the hash interface
 * work.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) 2004-2007 Craig DeForest.  You may distribute this
 * file under the terms of the Gnu Public License (GPL), version 2.
 * You should have received a copy of the GPL with this file.
 * If not, you may retrieve it from "http://www.gnu.org".
 *
 * This is version 2.0 of Flux.xs, released 31-Oct-2007.
 *
 *
 * WHAT IS WHAT: 
 *
 * The helper function "fieldptr" translates type and field number 
 * codes (defined here and in Flux.pm) into structure fields in the 
 * main FLUX structures defined in data.h. It is used to isolate 
 * knowledge about the C structures into one location, so that the 
 * individual access routines don't have to know about the FLUX data 
 * structures.
 * 
 * Most of the access routines also use the constructor defined in 
 * Core.xs and present in the FLUX core.
 * 
 * XS methods defined here are:
 * 
 * _new_from_ptr    Constructor for Perl interfaces to FLUX C
 *		    structures - this just trampolines into 
 *                  the generic constructor in Core.xs. 
 * destroy_sv	    Destructor for Perl interface objects -- this
 *		      just trampolines into the generic destructor
 *                    in Core.xs.
 *
 * _rnum	    Reads a NUM (floating-point) field as a Perl NV.
 * _wnum	    Writes a NUM (floating-point) field from a Perl scalar.
 *
 * _rlong	    Reads a long (integer) field as a Perl NV.
 * _wlong	    Writes a long (integer) field from a Perl scalar.
 *
 * _rvec	    Reads a vector (3-NUM) field as a PDL.
 * _wvec 	    Writes a vector (3-NUM) field from a PDL.
 *
 * _rvertex	    Reads a vertex field into a newly-constructed Flux::Vertex object.
 * _wvertex 	    Writes a vertex pointer from a Flux::Vertex object.
 * 
 * _rfluxon	    Reads a fluxon field into a newly-constructed Flux::Fluxon object.
 * _wfluxon	    Writes a fluxon pointer from a Flux::Fluxon object.
 *
 * _rconcentration  Reads a fc field into a newly-constructed Flux::Concentration.
 *
 * _rcoeffs	    Reads the current acceleration coefficients into a Perl list ref.
 * _wcoeffs	    Writes the current acceleration coefficients from a Perl list ref.
 *
 * _rforces	    Reads the current force laws into a Perl list of force-law names.
 * _wforces	    Writes the current force laws from a Perl list of force-law names;
 * 		     prints an informative message if a name can't be found.
 *
 * _rbound	    Reads a fluxon end condition function into a Perl scalar name.
 * _wbound	    Writes a fluxon end condition function from a Perl scalar name;
 * 		     prints an informative message if a name can't be found.
 *
 * _rrecon	    Reads the current reconnection conditions into a Perl list ref of
 * 			list refs, each of which contains (name,params) for a condition.
 * _wrecon 	    Writes the current reconnection conditions from a Perl list ref, each
 *                      elemtnt of which contains a test name and the parameters supplied.
 * 
 * _rworld	    Reads a world field into a newly-constructed Flux::World object.
 *  	  	    (Note that there is no _wworld -- it's neither easy nor desirable 
 *                  to change the WORLD to which objects belong.
 *
 * _rdumblist	    Reads a DUMBLIST into a Perl list of appropriate objects.  
 *                
 * file_versions    Returns a Perl string containing exact version info for each 
 * 		    C library file used to build the current libflux.a.
 *
 * atan2_oct  	    Access to atan2_oct (the sleazy octagonal atan2) in geometry.c
 *
 * cross_2d	    Access to cross_2d in geometry.c
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
  case FT_VERTEX: {
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
   case FT_FLUXON: {
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
     case FT_WORLD: {
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
		case 32: return (void *)&(w->photosphere2);		break;
		default: fprintf(stderr,"Unknown type,field (%d,%d) in Flux::World::fieldptr!\n",
				typeno, fieldno);		
			 return (void *)0;
		         break;
	 }
	}
	break;
     case FT_CONC: {
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
_new_from_ptr(wptr, type, label)
 IV type
 IV wptr
 IV label
CODE:
 RETVAL = FLUX->new_sv_from_ptr((WORLD *)wptr, type, label);
OUTPUT:
 RETVAL

void 
destroy_sv(me)
 SV *me
CODE:
 FLUX->destroy_sv(me);

SV *
_rnum(sv,typeno,fieldno)
  SV *sv
  long typeno
  long fieldno
PREINIT:
  void *ptr;
  NUM *np;
CODE:
	ptr = (void *)(FLUX->SvFluxPtr(sv,"_rnum","Flux",0,1));
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
	ptr = (void *)(FLUX->SvFluxPtr(sv,"_wnum","Flux",0,1));
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

	ptr = (void *)(FLUX->SvFluxPtr(sv,"_rlong","Flux",0,1));

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
	ptr = (void *)(FLUX->SvFluxPtr(sv,"_wlong","Flux",0,1));
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
	ptr = (void *)(FLUX->SvFluxPtr(sv,"_rvec","Flux",0,1));
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
	ptr = (void *)(FLUX->SvFluxPtr(sv,"_wvec","Flux",0,1));
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
	ptr = (void *)(FLUX->SvFluxPtr(sv,"_rvertex","Flux",0,1));
	vp = (VERTEX **)fieldptr(ptr,typeno,fieldno);
	RETVAL = (
		(vp && *vp) ? 
		 FLUX->new_sv_from_ptr((*vp)->line->fc0->world, FT_VERTEX, (*vp)->label) :
		 &PL_sv_undef
		);
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
	ptr = (void *)(FLUX->SvFluxPtr(sv,"_rfluxon","Flux",0,1));
	fp = (FLUXON **)fieldptr(ptr,typeno,fieldno);
	RETVAL = (	
		(fp && *fp) ?
		FLUX->new_sv_from_ptr((*fp)->fc0->world, FT_FLUXON, (*fp)->label) :
		&PL_sv_undef
		);
OUTPUT:
  	RETVAL

SV *
_rconcentration(sv, typeno, fieldno)
 SV *sv
 long typeno
 long fieldno
PREINIT:
 void *ptr;
 FLUX_CONCENTRATION **fcp;
CODE:
	ptr = (void *)(FLUX->SvFluxPtr(sv,"_rconcentration","Flux",0,1));
	if(!ptr) {
		RETVAL = &PL_sv_undef;
	} else {
		fcp = (FLUX_CONCENTRATION **)fieldptr(ptr,typeno,fieldno);
		RETVAL = (
		   (fcp && *fcp) ?
		   FLUX->new_sv_from_ptr( (*fcp)->world, FT_CONC, (*fcp)->label ) :
		   &PL_sv_undef
			);
	}
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
 world = (void *)(FLUX->SvFluxPtr(sv, "_rcoeffs", "Flux::World",0,1));

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

 world = (void *)(FLUX->SvFluxPtr(sv,"_wcoeffs","Flux::World",0,1));

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
 world = (void *)(FLUX->SvFluxPtr(sv,"_rforces", "Flux::World",0,1));
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

 world = (void *)(FLUX->SvFluxPtr(sv, "_wforces","Flux::World",0,1));
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
 SV *rv;
 void *ptr;
 void **field;
CODE:
 I32 i;
 ptr = (void *)(FLUX->SvFluxPtr(sv,"_rbound", "Flux",0,1));
 field = fieldptr(ptr,typeno, fieldno);

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

 ptr = (void *)(FLUX->SvFluxPtr(sv, "_wbound","Flux",0,1));
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



SV *
_rrecon(sv, typeno, fieldno)
 SV *sv
 long typeno
 long fieldno
PREINIT:
 AV *av;
 AV *params_av;
 SV *params;
 SV *rv;
 WORLD *w;
 void *ptr;
 void **field;
 I32 i;
 int ii;
 int jj;
CODE:
 w = SvWorld(sv,"_rrecon",1);

 av = newAV();

 for(ii=0; ii<N_RECON_FUNCS && w->rc_funcs[ii]; ii++) {
	 for(i=0;FLUX->FLUX_RECON[i].func && 
		( (void *)(FLUX->FLUX_RECON[i].func) != w->rc_funcs[ii]);
		i++)
		;
	 if(((void *)FLUX->FLUX_RECON[i].func) == w->rc_funcs[ii]) {
	     sv = newSVpv(FLUX->FLUX_RECON[i].name,
		 strlen(  FLUX->FLUX_RECON[i].name));

	     
	     params_av = newAV();
	     av_push(params_av, sv);

	     for(jj=0; jj < N_RECON_PARAMS; jj++) {
		av_push(params_av, newSVnv( w->rc_params[ii][jj] ));
	     }
	     av_push(av, (SV *)newRV_noinc((SV *)params_av)); 
	 } else {
	     char s[80];
	     sprintf(s,"UNKNOWN (0x%x)",(unsigned long)(w->rc_funcs[i]));
	     sv = newSVpv(s,strlen(s));

	     params_av = newAV();
	     av_push(params_av, sv);
	     av_push(av, (SV *)newRV_noinc((SV *)params_av));
         }
 }
 RETVAL = (SV *)newRV_noinc((SV *)av);
OUTPUT:
 RETVAL


SV *
_wrecon(sv, typeno, fieldno, val)
 SV *sv
 long typeno
 long fieldno
 SV *val
PREINIT:
 WORLD *w;
 int i;
 int j;
 RC_FUNC *rc_funcs[N_RECON_FUNCS];
 NUM rc_params[N_RECON_FUNCS][N_RECON_PARAMS];
CODE:
 w = SvWorld(sv, "_wrecon",1);

 for(i=0;i<N_RECON_FUNCS;i++) {
    rc_funcs[i] = 0;
    for(j=0;j<N_RECON_PARAMS;j++)
	rc_params[i][j] = 0;
 }
 if( val && val != &PL_sv_undef ) {   
   AV *val_av;
   static char errbuf[1024];
   if( ! SvROK ( val ) || SvTYPE(SvRV(val)) != SVt_PVAV ) 
	croak("_wrecon: error - reconnection field must be an array ref!\n");
   val_av = (AV *)SvRV(val);
   for(i=0; i<=av_len(val_av); i++) {
	SV **field = av_fetch(val_av, i, 0);
	AV *f_av;
	SV **svp;
	char *name;

	if(!field) {
		sprintf(errbuf,"_wrecon: error - null field in field %d of input\n",i);
		croak(errbuf);
	}
	
	if(!SvROK( *field ) || SvTYPE(SvRV(*field)) != SVt_PVAV) {
		sprintf(errbuf,"_wrecon: error - field %d of input isn't a list ref\n",i);
		croak(errbuf);
	}
	f_av= (AV *)SvRV(*field);

	svp = av_fetch(f_av, 0, 0);
	if(!svp) {
		sprintf(errbuf,"_wrecon: error - null name in field %d\n",i);
		croak(errbuf);
	}
	name = SvPV_nolen(*svp);
	
	// Do a linear search for the reconnection-function name in the table 
	for(j=0; FLUX->FLUX_RECON[j].func && strcmp(FLUX->FLUX_RECON[j].name,name); j++) 
		;

	//Not there -- generate a nice error message.
	if( !FLUX->FLUX_RECON[j].func ) {
		int k;
		sprintf(errbuf,"_wrecon: error - unknown name '%s' in field %d; valid names are: %s",
			name, i, FLUX->FLUX_RECON[0].name);
		for(k=1;FLUX->FLUX_RECON[k].func;k++) {
			strcat(errbuf, ", ");
			strcat(errbuf, FLUX->FLUX_RECON[k].name);
		}
		strcat(errbuf,".\n");
		croak(errbuf);
	}
	// There -- set the function in our stash.
	rc_funcs[i]=FLUX->FLUX_RECON[j].func;
	
	for(j=1; j <= av_len( f_av ) && j <= N_RECON_PARAMS; j++) {
		svp = av_fetch(f_av, j, 0);

		if(!svp || (*svp==&PL_sv_undef)) {
			printf("_wrecon: WARNING - null/undef param %d to %s (row %d) - proceeding.\n",j-1,name,i);
		} else {
			rc_params[i][j-1] = SvNV( *svp );
		}
	}
   }
 } 

 // Survived! Copy the cache into the world.
 for(i=0; i<N_RECON_FUNCS; i++) {
    w->rc_funcs[i] = rc_funcs[i];
    for(j=0; j<N_RECON_PARAMS; j++)
	w->rc_params[i][j] = rc_params[i][j];
 }
 
 RETVAL = val;
OUTPUT:
 RETVAL

SV *
_rworld(sv, typeno, fieldno)
 SV *sv
 long typeno
 long fieldno
PREINIT:
 void *ptr;
 WORLD **wp;
CODE:
	ptr = (void *)(FLUX->SvFluxPtr(sv,"_rfluxon","Flux",0,1));
	wp = (WORLD **)fieldptr(ptr,typeno,fieldno);
	RETVAL = ( 
		(wp && *wp) ? 
		FLUX->new_sv_from_ptr( (*wp), FT_WORLD, 0 ) :
		&PL_sv_undef
	      );
OUTPUT:
	RETVAL
	

SV *
_rdumblist(type, sv, typeno, fieldno)
 IV type
 SV *sv
 long typeno
 long fieldno
PREINIT:
 void *ptr;
 DUMBLIST *dl;
 void **array;
 AV *av;
 WORLD *w;
CODE:
 I32 i;
 ptr = (void *)(FLUX->SvFluxPtr(sv,"_rdumblist","Flux",0,1));
 w = (WORLD *)(FLUX->SvFluxPtr(sv, "_rdumblist","Flux",1,1));

 dl = (DUMBLIST *)fieldptr(ptr,typeno,fieldno);
 array = dl->stuff;

 av = newAV();         /* new perl list */
 av_clear(av);         /* Make sure it's empty */
 av_extend(av, dl->n); /* pre-extend for efficiency */

 /* printf("rdumblist: n=%d\n",dl->n);*/
 for(i=0;i<dl->n;i++) {
	SV *element;
	int label;
	switch(type) {
		case FT_VERTEX: label = ((VERTEX **)(dl->stuff))[i]->label; break;
		case FT_FLUXON: label = ((FLUXON **)(dl->stuff))[i]->label; break;
		case FT_WORLD:  label = 0; break;
		case FT_CONC:   label = ((FLUX_CONCENTRATION **)(dl->stuff))[i]->label; break;
		default: croak("Unknown type in _rdumblist!"); label=0; break;
	}
	element = FLUX->new_sv_from_ptr( w, type, label );
	av_push( av, element );
 }
 RETVAL = newRV_inc((SV*)av);
 /* printf("done\n");*/
OUTPUT:
 RETVAL		

SV *
_wdumblist(type, sv, typeno, fieldno, val)
 IV type
 SV *sv
 long typeno
 long fieldno
 SV *val
PREINIT:
 void *ptr;
 DUMBLIST *dl;
 void **array;
 AV *av;
 WORLD *w;
CODE:
 {
   I32 i,n;
   static DUMBLIST *dlcache = 0;

   if(type < MIN_FT || type > MAX_FT) 
     croak("_wdumblist: unknown type code!");

   ptr = (void *)( FLUX->SvFluxPtr(sv,"_wdumblist","Flux",0,1));
   w   = (WORLD *)(FLUX->SvFluxPtr(sv,"_wdumblist","Flux",1,1));
   dl = (DUMBLIST *)fieldptr(ptr, typeno, fieldno);
   array = dl->stuff;

   if( !val || val==&PL_sv_undef) {
     dl->n = 0;
     RETVAL = &PL_sv_undef;
   } else {
     if( !SvROK(val) || SvTYPE(SvRV(val)) != SVt_PVAV )
       croak("_wdumblist: requires an array ref or undef\n");
     av = (AV *)SvRV(val);
     n = av_len(av);

     if(!dlcache) {
       dlcache = FLUX->new_dumblist();
     }
     dlcache->n = 0;
     if(dlcache->size <= n) {
       FLUX->dumblist_grow(dlcache, n+1);
     }

     for(i=0;i<=n;i++) {
       SV **svp = av_fetch(av, i, 0);
       if(svp) {
	 ptr = FLUX->SvFluxPtr(*svp, "_wdumblist (array element)", classnames[type],0,0);
	 if(!ptr) 
	   croak("_wdumblist: one or more elements was not found in the World!");
	 dumblist_quickadd(dlcache, ptr);
       } else {
	 croak("_wdumblist: missing an element!");
       }
     }
     
     dl->n = 0;
     dumblist_snarf(dl, dlcache);
   }
   RETVAL = val;
 }
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

  
