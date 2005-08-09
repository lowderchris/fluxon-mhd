/* 
 * World.xs - glue code for the Flux::World object
 * in perl.
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
 *   _read_world        read_world (io.c) (leading '_' for ref in World.pm)
 *   write_world       write_world (io.c)
 *
 *   fix_proximity     global_fix_proximity (model.c)
 *   fix_curvature     global_fix_curvature (model.c)
 *  
 *   update_neighbors  world_update_neighbors (model.c)
 *
 *   fluxons           ----  Return a perl list of fluxon IDs in the world
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
 * Some helper routines that are used internally
 */

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

	w = read_world(f,(WORLD *)0);

	sv = newSViv((IV)(w));
	RETVAL = newRV_noinc(sv);
	(void)sv_bless(RETVAL,gv_stashpv("Flux::World",TRUE));
	fclose(f);
OUTPUT:
	RETVAL


void
write_world(wsv,s)
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
  w = SvWorld(wsv,"write_world");
  if( !(f=fopen(s,"w")) ) 
    croak("Couldn't open file '%s' to write a Flux::World",s);
  fprint_world(f,w,"");
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
  w = SvWorld(wsv,"fix_proximity");
  RETVAL = global_fix_proximity(w,alpha);
OUTPUT:
  RETVAL

int
fix_curvature(wsv,alpha)
 SV *wsv
 NV alpha
PREINIT:
 WORLD *w;
/**********************************************************************/
/* fix_curvature glue */
/**/
CODE:
  w = SvWorld(wsv,"fix_curvature");
  RETVAL = global_fix_curvature(w,alpha);
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
  w = SvWorld(wsv,"update_neighbors");
  world_update_neighbors(w,globalflag);


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
  w = SvWorld(wsv,"verbosity");
  if(items>1)
    w->verbosity = verbosity;
  RETVAL = w->verbosity;
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
 w = SvWorld(wsv,"Flux::World::fluxon_ids");
 RETVAL = fluxons_av = newAV();  /* fluxons_av is static, at top */
 av_clear(RETVAL);
 if(w->lines) {
	av_extend(RETVAL, (w->lines->all_links).n);
 	tree_walker(w->lines,fl_lab_of,fl_all_ln_of,fluxons_helper);
 }
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
  w = SvWorld(wsv,"Flux::World::fluxon");
  f = (FLUXON *)tree_find(w->lines,id,fl_lab_of,fl_all_ln_of);
  if(f) {
     sv = newSViv((IV)(f));
     RETVAL = newRV_noinc(sv);
     (void)sv_bless(RETVAL,gv_stashpv("Flux::Fluxon",TRUE));
  } else {
     RETVAL = &PL_sv_undef;
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
 w=SvWorld(wsv,"Flux::World::forces");

 av_clear(RETVAL = newAV());

 for(i=0;i<N_FORCE_FUNCS && w->f_funcs[i]; i++) {
   int j;
   for(j=0; 
       FLUX_FORCES[j].func && FLUX_FORCES[j].func != w->f_funcs[i];
       j++
       );
   if(FLUX_FORCES[j].func) {
     sv = newSVpv(FLUX_FORCES[j].name,strlen(FLUX_FORCES[j].name));
   } else {
     char s[80];
     sprintf(s,"0x%x",(unsigned long)(w->f_funcs[i]));
     sv = newSVpv(s,strlen(s));
   }
   av_store(RETVAL, av_len(RETVAL)+1, sv) || svREFCNT_dec(sv);
   }
OUTPUT:
 RETVAL

void
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
 * force pointer to it.
 */
CODE:
 w = SvWorld(wsv,"Flux::World::_set_force");
 if(where >= N_FORCE_FUNCS-1)
   croak("force index %d not allowed 0<i<%d\n",where,N_FORCE_FUNCS-1);

 if(*what==0)
   w->f_funcs[where]=0;
 else {
  for(j=0;
     FLUX_FORCES[j].func && strcmp(FLUX_FORCES[j].name,what);
     j++)
       ;
     if(!FLUX_FORCES[j].func) {
       if(what[0]=='0' && what[1]=='x') {
	 unsigned long ul;
	 sscanf(what+2,"%x",&ul);
	 (void *)(w->f_funcs[where]) = (void *)ul;
       } else {
         croak("Unknown force function '%s'\n",what);
       }
     } else {
        w->f_funcs[where] = FLUX_FORCES[j].func;
     }
  }


IV
_b_flag(wsv)
 SV *wsv
PREINIT:
 WORLD *w;
/**********************************************************************
 * _b_flag
 * Retrieve the state of the B normalization flag
 */
CODE:
 w=SvWorld(wsv,"Flux::World::forces");
 RETVAL = w->f_over_b_flag;
OUTPUT:
 RETVAL

void
_set_b_flag(wsv,flag)
SV *wsv
IV flag
PREINIT:
 WORLD *w;
/**********************************************************************
 * _set_b_flag
 * Set the state of the B normalization flag
 */
CODE:
 w = SvWorld(wsv,"Flux::World::_set_force");
 w->f_over_b_flag = flag;


AV *
_photosphere(wsv,plane=0)
SV *wsv
SV *plane
PREINIT:
 WORLD *w;
 SV *sv;
 AV *av;
 int i;
/**********************************************************************
 * _photosphere
 * Set the location of the photospheric plane from a perl array
 */
CODE:
 w = SvWorld(wsv,"Flux::World::_set_plane");
 av_clear(RETVAL = newAV());
 if(plane==0 || plane==&PL_sv_undef) {
   /* dump */
   if(w->photosphere.plane) {
     av_extend(RETVAL,6);
     PLANE *p = w->photosphere.plane;
     sv=newSVnv(p->origin[0]); av_store(RETVAL,0,sv) || svREFCNT_dec(sv);
     sv=newSVnv(p->origin[1]); av_store(RETVAL,1,sv) || svREFCNT_dec(sv);
     sv=newSVnv(p->origin[2]); av_store(RETVAL,2,sv) || svREFCNT_dec(sv);
     sv=newSVnv(p->normal[0]); av_store(RETVAL,3,sv) || svREFCNT_dec(sv);
     sv=newSVnv(p->normal[1]); av_store(RETVAL,4,sv) || svREFCNT_dec(sv);
     sv=newSVnv(p->normal[2]); av_store(RETVAL,5,sv) || svREFCNT_dec(sv);
   }
 } else {
  PLANE *newp;
  /* set */
  if(!SvROK(plane) || SvTYPE(SvRV(plane))!= SVt_PVAV) 
    croak("Flux::World::_set_plane: 2nd arg must be undef or array ref");
  av = (AV *)SvRV(plane);
  if(av_len(av)<0) {
    newp = 0;
  } else {
    SV **s1;
    (newp = (PLANE *)malloc(sizeof(PLANE))) || (croak("malloc failed"),(PLANE *)0);
    s1 = av_fetch(av,0,0); newp->origin[0] = SvNV(s1?*s1:&PL_sv_undef);
    s1 = av_fetch(av,1,0); newp->origin[1] = SvNV(s1?*s1:&PL_sv_undef);
    s1 = av_fetch(av,2,0); newp->origin[2] = SvNV(s1?*s1:&PL_sv_undef);
    s1 = av_fetch(av,3,0); newp->normal[0] = SvNV(s1?*s1:&PL_sv_undef);
    s1 = av_fetch(av,4,0); newp->normal[1] = SvNV(s1?*s1:&PL_sv_undef);
    s1 = av_fetch(av,5,0); newp->normal[2] = SvNV(s1?*s1:&PL_sv_undef);
  }
  if(w->photosphere.plane) {
    free(w->photosphere.plane);
  }
  w->photosphere.plane = newp;
  w->photosphere.type = (newp ? PHOT_PLANE : PHOT_NONE);
}
OUTPUT:
 RETVAL

NV
update_force(wsv,global=0)
SV *wsv
IV global
PREINIT:
 WORLD *w;
 NUM *minmax;
/**********************************************************************
 * update_force
 * Updates the forces everywhere in the model.
 */
CODE:
 w = SvWorld(wsv,"Flux::World::update_force");
 minmax = world_update_mag(w,global);
 RETVAL = minmax[1];
OUTPUT:
 RETVAL

void
relax_step(wsv,time)
SV *wsv
NV time
PREINIT:
 WORLD *w;
CODE:
 w = SvWorld(wsv,"Flux::World::relax_step");
 world_relax_step(w,time);


HV *
_stats(wsv)
SV *wsv
PREINIT:
 WORLD *w;
 struct VERTEX_STATS *st;
CODE:
 w = SvWorld(wsv,"Flux::World::stats");
 st = world_collect_stats(w);	
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
  ls_closest_approach(P0,P1,a,b,c,d);
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
  p_ls_closest_approach(P,a,b,c);
  RETVAL = newAV();
  av_clear(RETVAL);
  av_extend(RETVAL,3);
  av_store(RETVAL,0,newSVnv(P[0]));
  av_store(RETVAL,1,newSVnv(P[1]));
  av_store(RETVAL,2,newSVnv(P[2]));
OUTPUT:
 RETVAL

AV *
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
CODE:
/*******************************
 * _projmatrix: interface to the geometry.c projmatrix code
 */
 a[0] = a_0; a[1] = a_1; a[2] = a_2;
 b[0] = b_0; b[1] = b_1; b[2] = b_2;
 projmatrix(matrix,a,b);
 RETVAL = newAV();
 av_clear(RETVAL);
 av_extend(RETVAL,9);
 av_store(RETVAL,0,newSVnv(matrix[0]));
 av_store(RETVAL,1,newSVnv(matrix[1]));
 av_store(RETVAL,2,newSVnv(matrix[2]));
 av_store(RETVAL,3,newSVnv(matrix[3]));
 av_store(RETVAL,4,newSVnv(matrix[4]));
 av_store(RETVAL,5,newSVnv(matrix[5]));
 av_store(RETVAL,6,newSVnv(matrix[6]));
 av_store(RETVAL,7,newSVnv(matrix[7]));
 av_store(RETVAL,8,newSVnv(matrix[8]));
OUTPUT:
 RETVAL

AV *
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
CODE:
 a[0] = a_0; a[1] = a_1; a[2] = a_2;
 m[0] = m_0; m[1] = m_1; m[2] = m_2;
 m[3] = m_3; m[4] = m_4; m[5] = m_5;
 m[6] = m_6; m[7] = m_7; m[8] = m_8;
 mat_vmult_3d(b,m,a);
 RETVAL = newAV();
 av_clear(RETVAL);
 av_extend(RETVAL,3);
 av_store(RETVAL,0,newSVnv(b[0]));
 av_store(RETVAL,1,newSVnv(b[1]));
 av_store(RETVAL,2,newSVnv(b[2]));
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
