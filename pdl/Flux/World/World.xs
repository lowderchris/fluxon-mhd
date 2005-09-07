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
	 ((void **)(w->f_funcs))[where] = (void *)ul;
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
_photosphere(wsv,plane=0,type=0)
SV *wsv
SV *plane
SV *type
PREINIT:
 WORLD *w;
 SV *sv;
 AV *av;
 int i;
/**********************************************************************
 * _photosphere
 * Set the location of the photospheric plane from a perl array, or
 * dump it into a perl array.
 */
CODE:
 w = SvWorld(wsv,"Flux::World::_set_plane");
 av_clear(RETVAL = newAV());
 if(!plane || plane==&PL_sv_undef) {
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
     sv=newSViv(w->photosphere.type); av_store(RETVAL,6,sv) || svREFCNT_dec(sv);
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
	
  if(!type || type==&PL_sv_undef) {
	  w->photosphere.type = (newp ? PHOT_PLANE : PHOT_NONE);
  } else {
          w->photosphere.type = SvIV(type);
  }
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
   int dims[2];
CODE:
/*******************************
 * _projmatrix: interface to the geometry.c projmatrix code
 */
 a[0] = a_0; a[1] = a_1; a[2] = a_2;
 b[0] = b_0; b[1] = b_1; b[2] = b_2;
 projmatrix(matrix,a,b);
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
	int dims[2];
CODE:
 a[0] = a_0; a[1] = a_1; a[2] = a_2;
 m[0] = m_0; m[1] = m_1; m[2] = m_2;
 m[3] = m_3; m[4] = m_4; m[5] = m_5;
 m[6] = m_6; m[7] = m_7; m[8] = m_8;
 mat_vmult_3d(b,m,a);
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
	int dims[2];
CODE:
 a[0] = a_0; a[1] = a_1; a[2] = a_2;
 m[0] = m_0; m[1] = m_1; m[2] = m_2;
 m[3] = m_3; m[4] = m_4; m[5] = m_5;
 m[6] = m_6; m[7] = m_7; m[8] = m_8;
 vec_mmult_3d(b,m,a);
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
horde = new_dumblist();
dumblist_grow(horde,p->dims[1]);
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
hull_2d(hull,horde,rejects=new_dumblist());
/******************************
 * Create an output PDL
 */
o = PDL->create(PDL_PERM);
{
 PDL_Long dims[2];
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
