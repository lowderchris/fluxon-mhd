/* 
 * World.xs - glue code for the Flux::World object
 * in perl.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) 2004-2008 Craig DeForest.  You may distribute this
 * file under the terms of the Gnu Public License (GPL), version 2.
 * You should have received a copy of the GPL with this file.
 * If not, you may retrieve it from "http://www.gnu.org".
 *
 * This file is part of the FLUX 2.2 release, 22-Nov-2008
 *
 * Routines in here:
 *   PERL INTERFACE        FLUX SUBROUTINE / FUNCTION
 *
 *   _read_world           read_world (io.c) (leading '_' for ref in World.pm)
 *   write_world       	   write_world (io.c)
 *
 *   fix_proximity     	   global_fix_proximity (model.c)
 *   _fix_curvature        global_fix_curvature (model.c); accessed mainly by
 *                         fix_curvature from World.pm.
 *  
 *   update_neighbors  	   world_update_neighbors             (model.c)
 *   f_length_check   	   world_fluxon_length_check          (model.c)
 *   update_force          world_update_mag                   (model.c)
 *   update_force_parallel world_update_mag_parallel          (model.c)
 *   relax_step            world_relax_step                   (model.c)
 *   relax_step_parallel   world_relax_step_parallel          (model.c)
 *
 *   verbosity         	NA (sets verbosity of the world; obviated by tied-hash interface)
 *   step_scales	Sets the scaling powers for step law.  Obviated by tied-hash interface.
 * 
 *   _fluxon_ids        Returns a perl list of all fluxon labels in the world
 *   _vertex_ids        Returns a perl list of all vertex labels in the world
 * 
 *   concentration      Returns a perl Flux::Concentration tied hash associated with the given id/label
 *   fluxon             Returns a perl Flux::Fluxon tied hash associated with the given id/label
 *   vertex  		Returns a perl Flux::Vertex tied hash associated with the given id/label
 *
 *   _new_concentration new_flux_concentration (data.c)
 *
 *   _forces		Returns the force list for the world as a set of identifier strings
 * 			  (obviated by the tied-hash interface)
 *   _set_force         Tries to interpret a string as a force and set the appropriate force pointer to it.
 *			  (obviated by the tied-hash interface)
 *
 *   _rcfuncs           Retrieve the reconnection criteria (obviated by the tied-hash interface)
 *   _set_rc            Set the reconnection criteria (obviated by the tied-hash interface)
 *   reconnect          global_recon_check (model.c)
 * 
 *   _b_flag            Legacy routine simply throws an error.
 *   _set_b_flag        Legacy routine simply throws an error.
 * 
 *   _photosphere       Sets the photospheric plane from a perl array
 *			  (obviated by the tied-hash interface, but don't delete: _wphot calls it!)
 *
 *   _stats             Collects force statistics over the whole world
 *
 *   _ls_dist_test      test harness for ls_dist in geometry.c
 *   _p_ls_dist_test    test harness for p_ls_dist in geometry.c
 *   _projmatrix        test harness for projmatrix in geometry.c
 *   _mat_vmult_3d      test harness for mat_vmult_3d in geometry.c
 *   _vec_mmult_3d      test harness for vec_mmult_3d in geometry.c
 *   _hull_points       test harness for hull_2d in geometry.c
 *
 *   _closest_vertex	returns the vertex that is closest to a particular location.  Uses elementary 
 *                      threading -- use the wrapper ("closest_vertex") in World.pm instead.
 *   _closest_simplex	Returns a list of VERTEXes describing a simplex close to the target point. 
 *                      Seems to not work very well -- needs replacement with a qhull-based 
 *                      triangulator.
 *   _interpolate_value Uses closest_simplex to find the surroundign simplex around a particular location,
 *                      and returns a value interpolated from the surrounding VERTEXes.  Requires the 
 *                      offset of the value from teh VERTEX structure -- it must be a NUM.
 *
 *   check              Runs basic tests on a World.
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
 * Some helper routines that are used internally
 */

/** fluxons_helper: callback to stuff a fluxon label into a list **/
static AV *fluxons_av;
static long fluxons_helper(void *tree, int lab_of, int ln_of, int depth) {
 SV *sv = newSViv(((FLUXON *)tree)->label);
 av_store(fluxons_av, av_len(fluxons_av)+1, sv) || svREFCNT_dec(sv);
 return 0;
}

static AV *vertices_av;
/** vertices_helper: callback to stuff a vertex label into a list **/
static long vertices_helper(void *tree, int lab_of, int ln_of, int depth) {
 SV *sv = newSViv(((VERTEX *)tree)->label);
 av_store(vertices_av, av_len(vertices_av)+1, sv) || svREFCNT_dec(sv);
 return 0;
}



/**********************************************************************
 **********************************************************************
 ***
 *** XS definitions follow...
 ***/
MODULE = Flux::World    PACKAGE = Flux::World

SV *
_read_world(s)
 char *s
PREINIT:
 FILE *f;
 WORLD *w;
 SV *sv;
/**********************************************************************/
/* read_world glue */
/**/
CODE:
	if(  !(f = fopen(s,"r"))  ) 
	    croak("Couldn't open file '%s' to read a Flux::World",s);

	w = FLUX->read_world(f,(WORLD *)0);
	fclose(f);


	if(w) {
		printf("Read '%s': %d vertices, %d fluxons, %d concentrations\n", 
			s,
			w->vertices->world_links.n,
			w->lines->all_links.n,
			w->concentrations->links.n
			);
	} else {
		printf("Error reading '%s'\n",s);
	}

	RETVAL = (w ? 
		FLUX->new_sv_from_ptr(w, FT_WORLD, 0) :
		&PL_sv_undef
	);
OUTPUT:
  RETVAL


void
_write_world(wsv,s)
  SV *wsv
  char *s
PREINIT:
  WORLD *w;
  FILE *f;
  SV *sv;
/**********************************************************************/
/* write_world glue */
/**/
CODE:
  w = SvWorld(wsv,"write_world",1);
  if( !(f=fopen(s,"w")) ) 
    croak("Couldn't open file '%s' to write a Flux::World",s);
  FLUX->fprint_world(f,w,"");
  fclose(f);


int
fix_proximity(wsv,alpha)
 SV *wsv
 NV alpha
PREINIT:
 WORLD *w;
/**********************************************************************/
/* fix_proximity glue */
/**/
CODE:
  w = SvWorld(wsv,"fix_proximity",1);
  RETVAL = FLUX->global_fix_proximity(w,alpha);
OUTPUT:
  RETVAL

int
_fix_curvature(wsv,alpha,beta)
 SV *wsv
 NV alpha
 NV beta
PREINIT:
 WORLD *w;
/**********************************************************************/
/* fix_curvature glue */
/**/
CODE:
  w = SvWorld(wsv,"fix_curvature",1);
  RETVAL = FLUX->global_fix_curvature(w,alpha,beta);
OUTPUT:
  RETVAL

void
update_neighbors(wsv, globalflag)
 SV *wsv
 IV globalflag
PREINIT:
 WORLD *w;
/**********************************************************************
 * update_neighbors glue
 */
CODE:
  w = SvWorld(wsv,"update_neighbors",1);
  FLUX->world_update_neighbors(w,globalflag);

void
f_length_check(wsv, globalflag)
 SV *wsv
 IV globalflag
PREINIT:
 WORLD *w;
/**********************************************************************
 * world_fluxon_length_check glue
 */
CODE:
  w = SvWorld(wsv,"f_length_check",1);
  FLUX->world_fluxon_length_check(w,globalflag);

NV
update_force(wsv,global=0)
SV *wsv
IV global
PREINIT:
 WORLD *w;
/**********************************************************************
 * update_force
 * Updates the forces everywhere in the model.
 */
CODE:
 w = SvWorld(wsv,"Flux::World::update_force",1);
 FLUX->world_update_mag(w,global);
 w->f_max;
OUTPUT:
 RETVAL

NV
update_force_parallel(wsv,global=0)
SV *wsv
IV global
PREINIT:
 WORLD *w;
/**********************************************************************
 * update_force
 * Updates the forces everywhere in the model, using the 
 * first-stage parallel entry point
 */
CODE:
 w = SvWorld(wsv,"Flux::World::update_force",1);
 FLUX->world_update_mag_parallel(w,global);
 RETVAL = w->f_max;
OUTPUT:
 RETVAL

void
relax_step(wsv,time)
SV *wsv
NV time
PREINIT:
 WORLD *w;
 /**********************************************************************
  * relax_step
  * glue for world_relax_step in model.c
  */
CODE:
 w = SvWorld(wsv,"Flux::World::relax_step",1);
 FLUX->world_relax_step(w,time);

void
relax_step_parallel(wsv,time)
SV *wsv
NV time
PREINIT:
 WORLD *w;
 /*********************************************************************
  * relax_step_parallel
  * glue for world_relax_step_parallel in model.c
  */
CODE:
 w = SvWorld(wsv,"Flux::World::relax_step",1);
 FLUX->world_relax_step_parallel(w,time);


IV
verbosity(wsv,verbosity=0)
 SV *wsv
 IV verbosity
PREINIT:
 WORLD *w;
 IV v;
/**********************************************************************
 * verbosity glue - set the verbosity of activity with this World
 */
CODE:
  w = SvWorld(wsv,"verbosity",1);
  if(items>1)
    w->verbosity = verbosity;
  RETVAL = w->verbosity;
OUTPUT:
  RETVAL


SV *
step_scales(wsv, href=&PL_sv_undef)
 SV *wsv
 SV *href
PREINIT:
 WORLD *w;
 HV *hash=0;
 SV *val, **valp;
/**********************************************************************
 * step_scales - import/export scaling laws
 */
CODE:
 w = SvWorld(wsv,"step_scales",1);
 if(items>1 && SvROK(href) && SvTYPE(SvRV(href))==SVt_PVHV) {
   /* Copy hash ref into WORLD */
   hash = (HV *)SvRV(href);

   /* Copy the hash values into the WORLD */
   valp = hv_fetch(hash, "b", 1, 0); 
   if(valp && *valp && *valp != &PL_sv_undef) 
	w->step_scale.b_power = SvNV(*valp);

   valp = hv_fetch(hash, "d", 1, 0); 
   if(valp && *valp && *valp != &PL_sv_undef) 
	w->step_scale.d_power = SvNV(*valp);

   valp = hv_fetch(hash, "s", 1, 0); 
   if(valp && *valp && *valp != &PL_sv_undef) 
	w->step_scale.s_power = SvNV(*valp);

   valp = hv_fetch(hash, "ds", 2, 0); 
   if(valp && *valp && *valp != &PL_sv_undef) 
	w->step_scale.ds_power = SvNV(*valp);
  }
  /* Generate an output value */
  hash=newHV();
  hv_clear(hash);
  hv_store(hash, "b", 1, newSVnv( w->step_scale.b_power), 0);
  hv_store(hash, "d", 1, newSVnv( w->step_scale.d_power), 0);
  hv_store(hash, "s", 1, newSVnv( w->step_scale.s_power), 0);
  hv_store(hash, "ds",2, newSVnv( w->step_scale.ds_power),0);
  RETVAL = newRV_noinc((SV *)hash);
OUTPUT:
  RETVAL 


AV *
_fluxon_ids(wsv)
 SV *wsv
PREINIT:
 WORLD *w;
/**********************************************************************
 * _fluxon_ids - generate a perl list of fluxon IDs associated with this
 * world 
 * 
 * Hands back an array ref instead of a list, so I don't have to hassle
 * with the perl stack...
 */
CODE:
 w = SvWorld(wsv,"Flux::World::_fluxon_ids",1);
 RETVAL = fluxons_av = newAV();  /* fluxons_av is static, at top */
 av_clear(RETVAL);
 if(w->lines) {
	av_extend(RETVAL, (w->lines->all_links).n);
 	FLUX->tree_walker(w->lines,fl_lab_of,fl_all_ln_of,fluxons_helper,0);
 }
OUTPUT:
 RETVAL


AV *
_vertex_ids(wsv)
 SV *wsv
PREINIT:
 WORLD *w;
/**********************************************************************
 * _vertex_ids - generate a perl list of vertex IDs associated with this world
 *
 * Hands back an array ref instead of a list, so I don't have to hassle
 * with the perl stack...
 */
CODE:
 w = SvWorld(wsv,"Flux::World::_vertex_ids",1);
 RETVAL = vertices_av = newAV(); /* vertices_av is static, at top */
 av_clear(RETVAL);
 if(w->vertices) {
	av_extend(RETVAL, (w->vertices->world_links).n);
	FLUX->tree_walker(w->vertices, v_lab_of, v_ln_of, vertices_helper,0);
 }
OUTPUT:
 RETVAL

SV *
concentration(wsv,id)
 SV *wsv
 IV id
PREINIT:
 WORLD *w;
 FLUX_CONCENTRATION *fc;
 SV *sv;
/**********************************************************************
 * concentration - generate and return a Flux::Concentration object
 * associated with the given id
 */
CODE:
 w = SvWorld(wsv,"Flux::World::concentration",1);
 fc = (FLUX_CONCENTRATION *)(FLUX->tree_find(w->concentrations,id,fc_lab_of, fc_ln_of));
 RETVAL = (fc ?
	FLUX->new_sv_from_ptr(w, FT_CONC, fc->label) :
	&PL_sv_undef
	);
OUTPUT:
 RETVAL

SV *
fluxon(wsv,id)
 SV *wsv
 IV id
PREINIT:
 WORLD *w;
 FLUXON *f;
 SV *sv;
/**********************************************************************
 * fluxon - generate and return a Flux::Fluxon object associated 
 * with the given id
 */
CODE:
  w = SvWorld(wsv,"Flux::World::fluxon",1);
  f = (FLUXON *)(FLUX->tree_find(w->lines,id,fl_lab_of,fl_all_ln_of));
  RETVAL = (f ? 
	FLUX->new_sv_from_ptr(w, FT_FLUXON, f->label) :
	&PL_sv_undef
	);
OUTPUT:
  RETVAL


SV *
vertex(wsv,id)
 SV *wsv
 IV id
PREINIT:
 WORLD *w;
 VERTEX *v;
 SV *sv;
/**********************************************************************
 * vertex - generate and return a Flux::Vertex object associated 
 * with the given id
 */
CODE:
  w = SvWorld(wsv,"Flux::World::vertex",1);
  v = (VERTEX *)(FLUX->tree_find((void *)(w->vertices),id,v_lab_of,v_ln_of));
  RETVAL = (v ?
	FLUX->new_sv_from_ptr(w, FT_VERTEX, v->label) :
	&PL_sv_undef
	);
OUTPUT:
  RETVAL

	
SV *
_new_concentration( wsv, xsv, flux, label, endsv );
 SV *wsv;
 SV *xsv;
 NV flux;
 IV label;
 SV *endsv;
PREINIT:
 WORLD *w;
 pdl *x;
 FLUX_CONCENTRATION *fc;
 /*********************************************************************
  * _new_concentration - glue for new_flux_concentration in data.c
  * internal - called only from new_concentration in World.pm...
  */
CODE:
 w = SvWorld(wsv, "Flux::World::new_concentration",1);
 x = PDL->SvPDLV(xsv);
 if(!x || x->ndims != 1 || x->dims[0] != 3) {
	croak("Flux::World::new_concentration requires a 3-PDL location!");
 }
 PDL->converttype(&x, PDL_D, 1);
 PDL->make_physical(x);
 fc = FLUX->new_flux_concentration(w, 
	((PDL_Double *)(x->data))[0], 
	((PDL_Double *)(x->data))[1],
	((PDL_Double *)(x->data))[2],
	flux,
	label
	);
 RETVAL = FLUX->new_sv_from_ptr(w, FT_CONC, fc->label);

 // Set the boundary, if specified...
 if(endsv && endsv != &PL_sv_undef && *(SvPV_nolen(endsv))) {
	I32 foo;
	ENTER;
	SAVETMPS;
	PUSHMARK(SP);
	XPUSHs(RETVAL);
	XPUSHs(sv_2mortal(newSViv( 5  ))); // magic number for flux concentrations
	XPUSHs(sv_2mortal(newSViv( 13 ))); // magic number for the bound field 
	XPUSHs(sv_2mortal(newSVsv(endsv)));
	PUTBACK;
	foo = call_pv("Flux::_wbound",G_SCALAR);
	SPAGAIN;
	PUTBACK;
	FREETMPS;
	LEAVE;
 }
OUTPUT:
 RETVAL


AV *
_forces(wsv)
 SV *wsv
PREINIT:
 WORLD *w;
 SV *sv;
 int i;
/**********************************************************************
 * forces
 * Retrieve the force list in the World, using strings and the 
 * conversion array defined at the top of physics.c
 */
CODE:
 w=SvWorld(wsv,"Flux::World::forces",1);

 av_clear(RETVAL = newAV());

 for(i=0;i<N_FORCE_FUNCS && w->f_funcs[i]; i++) {
   int j;
   for(j=0; 
       FLUX->FLUX_FORCES[j].func && FLUX->FLUX_FORCES[j].func != w->f_funcs[i];
       j++
       );
   if(FLUX->FLUX_FORCES[j].func) {
     sv = newSVpv(FLUX->FLUX_FORCES[j].name,strlen(FLUX->FLUX_FORCES[j].name));
   } else {
     char s[80];
     sprintf(s,"0x%x",(unsigned long)(w->f_funcs[i]));
     sv = newSVpv(s,strlen(s));
   }
   av_store(RETVAL, av_len(RETVAL)+1, sv) || svREFCNT_dec(sv);
   }
OUTPUT:
 RETVAL



IV
_set_force(wsv,where,what)
SV *wsv
IV where
char * what
PREINIT:
 WORLD *w;
 SV *sv;
 int j;
/**********************************************************************
 * _set_force
 * Try to interpret a string as a force name and set the 
 * force pointer to it.  Returns true on success, or false on failure
 */
CODE:
 w = SvWorld(wsv,"Flux::World::_set_force",1);
 if(where >= N_FORCE_FUNCS-1)
   croak("force index %d not allowed 0<i<%d\n",where,N_FORCE_FUNCS-1);

 if(*what==0) {
   w->f_funcs[where]=0;
   RETVAL=1;
 } else {
  for(j=0;
     FLUX->FLUX_FORCES[j].func && strcmp(FLUX->FLUX_FORCES[j].name,what);
     j++)
       ;
     if(!FLUX->FLUX_FORCES[j].func) {
       if(what[0]=='0' && what[1]=='x') {
	 unsigned long ul;
	 sscanf(what+2,"%x",&ul);
	 ((void **)(w->f_funcs))[where] = (void *)ul;
	 RETVAL = 3;
       } else {
         croak("Unknown force function '%s'\n",what);
       }
     } else {
        w->f_funcs[where] = FLUX->FLUX_FORCES[j].func;
	RETVAL = 2;
     }
  }
OUTPUT:
 RETVAL


AV *
_rcfuncs(wsv)
 SV *wsv
PREINIT:
 WORLD *w;
 SV *sv;
 SV *sv2;
 int i;
/**********************************************************************
 * _rcfuncs
 * Retrieve the reconnection criteria for the World, and return them in an
 * array ref.  Elements alternate between names of functions and lists of parameters.
 */
CODE:

  w=SvWorld(wsv,"Flux::World::_rcfuncs",1);
	
  if(w->verbosity) {
	printf("_rcfuncs...\n");
  }
  av_clear(RETVAL = newAV());
  for(i=0;i<N_RECON_FUNCS && w->rc_funcs[i]; i++) {
	int j;
	for( j=0;	
 	     FLUX->FLUX_RECON[j].func && FLUX->FLUX_RECON[j].func != w->rc_funcs[i];
	     j++
	);
	if(FLUX->FLUX_RECON[j].func) {
		sv = newSVpv(FLUX->FLUX_RECON[j].name,strlen(FLUX->FLUX_RECON[j].name));
	} else {
		char s[80];
		sprintf(s,"0x%x",(unsigned long)(w->rc_funcs[i]));
		sv=newSVpv(s,strlen(s));
	}
	if(w->verbosity)
		printf("x...");
	
	av_store(RETVAL, av_len(RETVAL)+1, sv) || svREFCNT_dec(sv);
	av_clear((AV *)(sv2 = (SV *)newAV()));
	for( j=0; j<N_RECON_PARAMS; j++ )  {
		sv = newSVnv(w->rc_params[i][j]);
		av_store((AV *)sv2, av_len((AV *)sv2)+1, sv) || svREFCNT_dec(sv);
	}
	sv2 = newRV_noinc(sv2);
	av_store(RETVAL, av_len(RETVAL)+1, sv2) || svREFCNTT_dec(sv2);
  }
OUTPUT:
 RETVAL


IV
_set_rc(wsv,where,what,params)
SV *wsv
IV where
char *what
SV *params
PREINIT:
	WORLD *w;
	SV *sv;
	AV *av;
	int i,j;
CODE:
	w = SvWorld(wsv, "Flux::World::_set_rc",1);
	/******************************
	* _set_rc 
	* Try to interpret a string as a reconnection criterion name 
	* and set the pointer appropriately.
	*/
	if(where >= N_RECON_FUNCS-1)
	  croak("reconnection criterion %d not allowed 0<i<%d\n",where,N_RECON_FUNCS-1);

	if(!what || !*what){
		w->rc_funcs[where]=0;
		RETVAL = 1;
	} else {
		for(j=0;
		FLUX->FLUX_RECON[j].func && strcmp(FLUX->FLUX_RECON[j].name,what);
		j++)
		;

		if(FLUX->FLUX_RECON[j].func) {
			w->rc_funcs[where] = FLUX->FLUX_RECON[j].func;
	
			if( params && 
			    params != &PL_sv_undef && 
			    (   !SvROK(params) || SvTYPE(SvRV(params)) != SVt_PVAV   )
		          )
				croak("_set_rc: requires an array ref params argument");

			else if( !params || params==&PL_sv_undef ) {
				for(i=0;i<N_RECON_PARAMS;i++) 
					w->rc_params[where][i] = 0;
			} else {
				av = (AV *)SvRV(params);
				for(i=0; i<N_RECON_PARAMS; i++) {
					SV **s1 = av_fetch( av, i, 0 );
					w->rc_params[where][i]= SvNV(s1?*s1:&PL_sv_undef);
				}
			}
			RETVAL = 2;
		} else {
			croak("Unknown reconnection criterion '%s'\n",what);
		}
	}
OUTPUT:
 RETVAL

IV 
reconnect(wsv)
SV *wsv
PREINIT:
	WORLD *w;
CODE:
	w = SvWorld(wsv, "Flux::World::reconnect",1);
	RETVAL = FLUX->global_recon_check(w);
OUTPUT:
 RETVAL

HV *
_stats(wsv)
SV *wsv
PREINIT:
 WORLD *w;
 struct VERTEX_STATS *st;
 /*********************************************************************
  * _stats - glue for world_collect_stats...
  */
CODE:
 w = SvWorld(wsv,"Flux::World::stats",1);
 st = FLUX->world_collect_stats(w);	
 RETVAL = newHV();
 hv_store(RETVAL, "n",    1, newSViv(st->n),0);
 if(st->n>0) {
   hv_store(RETVAL, "f_av",  4, newSVnv(st->f_acc/st->n),     0);	
   hv_store(RETVAL, "f_max", 4, newSVnv(st->f_max),           0);
   hv_store(RETVAL, "fs_av", 5, newSVnv(st->f_tot_acc/st->n), 0);
   hv_store(RETVAL, "fs_max",6, newSVnv(st->f_tot_max),       0);
   hv_store(RETVAL, "n_av",  4, newSVnv(st->n_acc/st->n),     0);
   hv_store(RETVAL, "n_max", 5, newSVnv(st->n_max),           0);
 }
OUTPUT:
 RETVAL
 

AV *
_ls_dist_test(a_0,a_1,a_2,b_0,b_1,b_2,c_0,c_1,c_2,d_0,d_1,d_2)
NV a_0
NV a_1
NV a_2
NV b_0
NV b_1
NV b_2
NV c_0
NV c_1
NV c_2
NV d_0
NV d_1
NV d_2
PREINIT:
  NUM a[3],b[3],c[3],d[3],P0[3],P1[3];
CODE:
/******************************
 *
 * _ls_dist_test: debuggging cruft.
 * This method should properly never be called, but it is present to 
 * enable testing of the ls_closest_approach function in the libraries.
 * You feed in two line-segments and get back the closest approach points.
 * The interface is rather primitive.
 */
  a[0] = a_0; a[1]=a_1; a[2]=a_2;
  b[0] = b_0; b[1]=b_1; b[2]=b_2;
  c[0] = c_0; c[1]=c_1; c[2]=c_2;
  d[0] = d_0; d[1]=d_1; d[2]=d_2;
  FLUX->ls_closest_approach(P0,P1,a,b,c,d);
  RETVAL = newAV();
  av_clear(RETVAL);
  av_extend(RETVAL,6);
  av_store(RETVAL, 0,newSVnv(P0[0]));
  av_store(RETVAL, 1,newSVnv(P0[1]));
  av_store(RETVAL, 2,newSVnv(P0[2]));
  av_store(RETVAL, 3,newSVnv(P1[0]));
  av_store(RETVAL, 4,newSVnv(P1[1]));
  av_store(RETVAL, 5,newSVnv(P1[2]));
OUTPUT:
 RETVAL

AV *
_p_ls_dist_test(a_0,a_1,a_2,b_0,b_1,b_2,c_0,c_1,c_2)
NV a_0
NV a_1
NV a_2
NV b_0
NV b_1
NV b_2
NV c_0
NV c_1
NV c_2
PREINIT:
  NUM a[3],b[3],c[3],P[3];
CODE:
/******************************
 *
 * _p_ls_dist_test: debugging cruft.
 * Test the p_ls_closest_approach function
 */
  a[0] = a_0; a[1]=a_1; a[2]=a_2;
  b[0] = b_0; b[1]=b_1; b[2]=b_2;
  c[0] = c_0; c[1]=c_1; c[2]=c_2;
  FLUX->p_ls_closest_approach(P,a,b,c);
  RETVAL = newAV();
  av_clear(RETVAL);
  av_extend(RETVAL,3);
  av_store(RETVAL,0,newSVnv(P[0]));
  av_store(RETVAL,1,newSVnv(P[1]));
  av_store(RETVAL,2,newSVnv(P[2]));
OUTPUT:
 RETVAL

SV *
_projmatrix(a_0,a_1,a_2,b_0,b_1,b_2)
NV a_0
NV a_1
NV a_2
NV b_0
NV b_1
NV b_2
PREINIT:
   NUM a[3],b[3];
   NUM matrix[9];
   PDL_Double *pd;
   NUM *d;
   pdl *p;
   PDL_Indx dims[2];
CODE:
/*******************************
 * _projmatrix: interface to the geometry.c projmatrix code
 */
 a[0] = a_0; a[1] = a_1; a[2] = a_2;
 b[0] = b_0; b[1] = b_1; b[2] = b_2;
 FLUX->projmatrix(matrix,a,b);
 p = PDL->create(PDL_PERM);
 dims[0] = 3;
 dims[1] = 3;
 PDL->setdims(p,dims,2);
 p->datatype = PDL_D;
 PDL->make_physical(p);
 pd = p->data;
 d = matrix;
 *(pd++) = *(d++);
 *(pd++) = *(d++);
 *(pd++) = *(d++);
 *(pd++) = *(d++);
 *(pd++) = *(d++);
 *(pd++) = *(d++);
 *(pd++) = *(d++);
 *(pd++) = *(d++);
 *(pd++) = *(d++);
 RETVAL = NEWSV(551,0); /* 551 is arbitrary */
 PDL->SetSV_PDL(RETVAL,p);
OUTPUT:
 RETVAL

SV *
_mat_vmult_3d(m_0,m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,a_0,a_1,a_2);
NV m_0
NV m_1
NV m_2
NV m_3
NV m_4
NV m_5
NV m_6
NV m_7
NV m_8
NV a_0
NV a_1
NV a_2
PREINIT:
	NUM a[3],b[3],m[9];
	pdl *p;
	PDL_Double *pd;
	PDL_Indx dims[2];
CODE:
 a[0] = a_0; a[1] = a_1; a[2] = a_2;
 m[0] = m_0; m[1] = m_1; m[2] = m_2;
 m[3] = m_3; m[4] = m_4; m[5] = m_5;
 m[6] = m_6; m[7] = m_7; m[8] = m_8;
 FLUX->mat_vmult_3d(b,m,a);
 p = PDL->create(PDL_PERM);
 dims[0] = 3;
 PDL->setdims(p,dims,1);
 p->datatype=PDL_D;
 PDL->make_physical(p);
 pd = (PDL_Double *)p->data;
 *(pd++) = b[0];
 *(pd++) = b[1];
 *(pd++) = b[2];
 RETVAL = NEWSV(552,0); /* 552 is arbitrary */
 PDL->SetSV_PDL(RETVAL,p);
OUTPUT:
 RETVAL

SV *
_vec_mmult_3d(m_0,m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,a_0,a_1,a_2);
NV m_0
NV m_1
NV m_2
NV m_3
NV m_4
NV m_5
NV m_6
NV m_7
NV m_8
NV a_0
NV a_1
NV a_2
PREINIT:
	NUM a[3],b[3],m[9];
	pdl *p;
	PDL_Double *pd;
	PDL_Indx dims[2];
CODE:
 a[0] = a_0; a[1] = a_1; a[2] = a_2;
 m[0] = m_0; m[1] = m_1; m[2] = m_2;
 m[3] = m_3; m[4] = m_4; m[5] = m_5;
 m[6] = m_6; m[7] = m_7; m[8] = m_8;
 FLUX->vec_mmult_3d(b,m,a);
 p = PDL->create(PDL_PERM);
 dims[0] = 3;
 PDL->setdims(p,dims,1);
 p->datatype=PDL_D;
 PDL->make_physical(p);
 pd = (PDL_Double *)p->data;
 *(pd++) = b[0];
 *(pd++) = b[1];
 *(pd++) = b[2];
 RETVAL = NEWSV(552,0); /* 552 is arbitrary */
 PDL->SetSV_PDL(RETVAL,p);
OUTPUT:
 RETVAL



SV *
_hull_points(pts)
SV *pts
PREINIT:
	pdl *p;
	pdl *o;
	DUMBLIST *horde,*rejects;
	HULL_VERTEX *hull;
	int i,j;
CODE:
/********************************************************************************
 * _hull_points - debugging routine for geometry.c code.  Accepts a 2xN PDL containing 2-D
 * points, and constructs a hull containing the origin and the points.  Returns the hull
 * in a 6xN PDL.  The (0,1) columns are the coordinates of the selected neighbors.  
 * (2,3) columns are the coordinates of each hull corner (if it is closed), or the 
 * angles of exit of the two branches of the hull (if they are open), and the (4) column
 * tells how to interpret the row (closed/open).  The (5) column contains the index number
 * of the corresponding original hull point. The (2,3) vertex goes to the LEFT of the 
 * corresponding (0,1) neighbor point.
 * (0,1) dimension is the hull point if it is closed, or 
 */
 printf("Parsing input pdl...\n");
 if(pts==NULL || pts==&PL_sv_undef) {
	fprintf(stderr,"_hull_points requires a 2xN PDL - got undef\n");
	RETVAL = &PL_sv_undef;
        goto _hull_exit;
 }
 p = PDL->SvPDLV(pts);
 if(!p) {
	fprintf(stderr,"_hull_points requires a 2xN PDL - but SvPDLV failed\n");
	RETVAL = &PL_sv_undef;
	goto _hull_exit;
}
if(p->ndims != 2 || p->dims[0] != 2) {
	fprintf(stderr,"_hull_points requires a 2xN PDL but got some other dimension\n");
	RETVAL = &PL_sv_undef;
	goto _hull_exit;
}
PDL->converttype(&p,PDL_D,1);
PDL->make_physdims(p);
PDL->make_physical(p);
printf("p is a %d-dim PDL (%d x %d)\n",p->ndims,p->dims[0],p->dims[1]);
/******************************
 * Stuff the PDL values into a collection of spankin'-new vertices
 */
horde = FLUX->new_dumblist();
FLUX->dumblist_grow(horde,p->dims[1]);
for(i=0;i<p->dims[1];i++) {
	double x,y;
	VERTEX *v = new_vertex(i,0,0,0,0);
	v->label = i;
	x = ((PDL_Double *)(p->data))[ i * (p->dimincs[1]) ];
	y = ((PDL_Double *)(p->data))[ i * (p->dimincs[1]) + p->dimincs[0] ];
	((VERTEX **)(horde->stuff))[i]  = v;
	v->scr[0] = x;
	v->scr[1] = y;
	v->scr[2] = 0;
	v->r = sqrt(x*x + y*y);
	v->a = atan2(y,x);
}
horde->n = p->dims[1];
/******************************
 * Allocate the hull data, and call hull_2d
 */
hull = (HULL_VERTEX *)malloc(sizeof(HULL_VERTEX) * p->dims[1]);
FLUX->hull_2d(hull,horde,rejects=new_dumblist());
/******************************
 * Create an output PDL
 */
o = PDL->create(PDL_PERM);
{
 PDL_Indx dims[2];
 dims[0]=6;
 dims[1]=horde->n;
 PDL->setdims(o, dims, 2);
}	
o->datatype = PDL_D;
PDL->resize_defaultincs(o);
PDL->make_physical(o);
/******************************
 * Copy the hull values into the output PDL
 */
for(i=0;i<horde->n;i++) {
  ((PDL_Double *)o->data)[ i * o->dimincs[1]                       ] = ((VERTEX *)(horde->stuff[i]))->scr[0];
  ((PDL_Double *)o->data)[ i * o->dimincs[1] +     o->dimincs[0] ] = ((VERTEX *)(horde->stuff[i]))->scr[1];
  ((PDL_Double *)o->data)[ i * o->dimincs[1] + 4 * o->dimincs[0] ] = hull[i].open;
  ((PDL_Double *)o->data)[ i * o->dimincs[1] + 5 * o->dimincs[0] ] = ((VERTEX *)(horde->stuff[i]))->label;
  if(hull[i].open){
   ((PDL_Double *)o->data)[ i * o->dimincs[1] + 2* o->dimincs[0] ] = hull[i].a_l;
   ((PDL_Double *)o->data)[ i * o->dimincs[1] + 3* o->dimincs[0] ] = hull[i].a_r;
  } else {
   ((PDL_Double *)o->data)[ i * o->dimincs[1] + 2 * o->dimincs[0] ] = hull[i].p[0];
   ((PDL_Double *)o->data)[ i * o->dimincs[1] + 3 * o->dimincs[0] ] = hull[i].p[1];
  }
}
/******************************
 * Clean up the data structures
 */
 for(i=0;i<horde->n;i++) {
	free(horde->stuff[i]);
 }
 free_dumblist(horde);
 for(i=0;i<rejects->n;i++) {
        free(rejects->stuff[i]);
 }
 free_dumblist(rejects);
 free(hull);
/******************************
 * Encapsulate the PDL in an SV for return
 */
 RETVAL = NEWSV(546,0); /* 546 is arbitrary tag */
 PDL->SetSV_PDL(RETVAL,o); 
 _hull_exit: ;
OUTPUT:
 RETVAL 

SV *
_closest_vertex(wsv, xsv, gsv=&PL_sv_undef, vsv=&PL_sv_undef)
SV *wsv
SV *xsv
SV *vsv
SV *gsv
PREINIT:
	/**************************************************
	 * closest_vertex - Flux::World interface to the 
	 * find_vertex_by_location routine in data.c.
	 * Just parse arguments and hand them off.  
         * This routine is internal because we want to be able
	 * to omit unused arguments.
	 */
	WORLD *w;
	VERTEX *v;
	POINT3D x;
	pdl *xpdl;
	int global;
CODE:
  w = SvWorld(wsv, "closest_vertex",1);
  printf("v is %d\n PL_sv_undef is %d\n",v,&PL_sv_undef);
  v = (!vsv || vsv==&PL_sv_undef) ? 0 : SvVertex(vsv, "closest_vertex",0);
  xpdl = PDL->SvPDLV(xsv);
  PDL->converttype(&xpdl,PDL_D,1);
  PDL->make_physical(xpdl);
  if(xpdl->ndims != 1 || xpdl->dims[0] != 3) {
	croak("closest_vertex: takes a 3-PDL only!\n");
  }
  x[0] = ((PDL_Double *)(xpdl->data))[0];
  x[1] = ((PDL_Double *)(xpdl->data))[1];
  x[2] = ((PDL_Double *)(xpdl->data))[2];
  global = SvIV(gsv);
  v = FLUX->find_vertex_by_location(x, w, v, global);
  if(v) {
    RETVAL = FLUX->new_sv_from_ptr(w, FT_VERTEX, v->label);
  } else {
    RETVAL = &PL_sv_undef;
  }
OUTPUT:
  RETVAL

AV *
_closest_simplex(wsv,xsv,gsv=&PL_sv_undef,vsv=&PL_sv_undef)
 SV *wsv
 SV *xsv
 SV *gsv
 SV *vsv
PREINIT:
 /**********************************************************************
  * closest_simplex -- return a list of VERTEXes describing a simplex 
  * close to the target point
  */
  WORLD *w;
  VERTEX *v;
  POINT3D x;
  pdl *xpdl;
  int global;
  DUMBLIST *dl;
  int i;
CODE:
  w = SvWorld(wsv,"closest_simplex",1);
  
  v = (!vsv || vsv==&PL_sv_undef)?0:SvVertex(vsv,"closest_vertex",0);
  xpdl = PDL->SvPDLV(xsv);
  PDL->converttype(&xpdl,PDL_D,1);
  PDL->make_physical(xpdl);
  if(xpdl->ndims != 1 || xpdl->dims[0] != 3) {
	croak("closest_vertex: takes a 3-PDL only!\n");
  }
  x[0] = ((PDL_Double *)(xpdl->data))[0];  
  x[1] = ((PDL_Double *)(xpdl->data))[1];  
  x[2] = ((PDL_Double *)(xpdl->data))[2];  
  
  global = SvIV(gsv);
  dl = FLUX->find_simplex_by_location(x, w, v, global);

  RETVAL = newAV();
  av_clear(RETVAL);
  av_extend(RETVAL,dl->n);
  for(i=0;i<dl->n;i++) {
    SV *sv = FLUX->new_sv_from_ptr(w, FT_VERTEX, ((VERTEX **)(dl->stuff))[i]->label);
    av_store(RETVAL, i, sv) || svREFCNT_dec(sv);
  }
OUTPUT:
 RETVAL  

NV
_interpolate_value( wsv, xsv, gsv, vsv, offset ) 
 SV *wsv
 SV *xsv
 SV *gsv
 SV *vsv
 IV offset
PREINIT:
 /***********************************************************************
 * _interpolate_value - takes a location and such as needed by _closest_simplex,
 * and a numerical offset to the NUM pointer you want to interpolate. 
 * It is normally called by interpolate_value, in World.pm, which calculated
 * the offset using the _ptr_offset routine in Flux.xs...
 */
 WORLD *w;
 VERTEX *v; 
 pdl *xpdl;
 POINT3D x;
 int global;
 DUMBLIST *dl;
 int i;
CODE:
 w = SvWorld(wsv, "_interpolate_value",1);
 v = (!vsv || vsv==&PL_sv_undef)?0:SvVertex(vsv, "closest_vertex",0);
 xpdl = PDL->SvPDLV(xsv);
 PDL->converttype(&xpdl, PDL_D, 1);
 PDL->make_physical(xpdl);
 if(xpdl->ndims != 1 || xpdl->dims[0] != 3) {
    croak("closest_vertex: takes a 3-PDL only!\n");
 }
 global = SvIV(gsv);

 x[0] = ((PDL_Double *)(xpdl->data))[0];
 x[1] = ((PDL_Double *)(xpdl->data))[1];
 x[2] = ((PDL_Double *)(xpdl->data))[2];

 RETVAL = FLUX->interpolate_value(x, w, v, global, offset);
OUTPUT:
 RETVAL

SV *
_interpolate_vector( wsv, xsv, gsv, vsv, offset ) 
 SV *wsv
 SV *xsv
 SV *gsv
 SV *vsv
 IV offset
PREINIT:
 /***********************************************************************
 * _interpolate_vector - takes a location and such as needed by _closest_simplex,
 * and a numerical offset to the vector you want to interpolate. 
 * It is normally called by interpolate_value, in World.pm, which calculated
 * the offset using the _ptr_offset routine in Flux.xs...
 * 
 * Returns a 3-PDL containing the interpolated vector.
 */
 WORLD *w;
 VERTEX *v; 
 pdl *xpdl;
 POINT3D x;
 int global;
 DUMBLIST *dl;
 pdl *opdl;
 NUM *pd;
 PDL_Indx dims[2];
 int i;
CODE:
 w = SvWorld(wsv, "_interpolate_value",1);
 v = (!vsv || vsv==&PL_sv_undef)?0:SvVertex(vsv, "closest_vertex",0);
 xpdl = PDL->SvPDLV(xsv);
 PDL->converttype(&xpdl, PDL_D, 1);
 PDL->make_physical(xpdl);
 if(xpdl->ndims != 1 || xpdl->dims[0] != 3) {
    croak("closest_vertex: takes a 3-PDL only!\n");
 }
 global = SvIV(gsv);

 x[0] = ((PDL_Double *)(xpdl->data))[0];
 x[1] = ((PDL_Double *)(xpdl->data))[1];
 x[2] = ((PDL_Double *)(xpdl->data))[2];

 dl = find_simplex_by_location(x,w,v,global);
 
 opdl = PDL->create(PDL_PERM);
 dims[0] = 3;
 PDL->setdims(opdl,dims,1);
 opdl->datatype = PDL_D;
 PDL->make_physical(opdl);
 pd = opdl->data;
 *(pd++) = FLUX->interpolate_value_simplex(x, dl, offset);
 *(pd++) = FLUX->interpolate_value_simplex(x, dl, offset + sizeof(NUM));
 *(pd++) = FLUX->interpolate_value_simplex(x, dl, offset + 2*sizeof(NUM));
 RETVAL = NEWSV(556,0); /* 556 is arbitrary */
 PDL->SetSV_PDL(RETVAL,opdl);
OUTPUT:
 RETVAL
 
IV
check( wsv )
 SV *wsv
PREINIT:
 WORLD *w;
CODE:
 /******************************
  * check - just springboard into world_check over in model.c..
  */
 w = SvWorld(wsv, "check",1);
 RETVAL = 0; // FLUX->world_check(w);
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
 if(FluxCoreSV == NULL)      Perl_croak(aTHX_ "Can't load Flux::Core module (required b Flux)");
 FLUX = INT2PTR(FluxCore*, SvIV(FluxCoreSV));
 if(FLUX->CoreVersion != FLUX_CORE_VERSION) {
	printf("FLUX->CoreVersion is %d; FLUX_%s is %d\n",FLUX->CoreVersion,"CORE_VERSION",FLUX_CORE_VERSION);
	Perl_croak(aTHX_ "Flux needs to be recompiled against the newly installed FLUX libraries");
}
