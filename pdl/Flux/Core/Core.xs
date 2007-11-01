/* Core.xs - utility routines and dynamic linker for the FLUX
 * C - Perl interface.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) 2004-2007 Craig DeForest.  You may distribute this
 * file under the terms of the Gnu Public License (GPL), version 2.
 * You should have received a copy of the GPL with this file.
 * If not, you may retrieve it from "http://www.gnu.org".
 *
 * This file is part of FLUX 2.0 (31-Oct-2007).
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

static FluxCore FLUX_MASTER;
static FluxCore *FLUX;
static SV* FluxCoreSV;

/**********************************************************************
 * SvFluxPtr - given a FLUX object and a type code, return 
 * a pointer to the underlying object.  Involves doing a search 
 * through the WORLD for the object.  Returns NULL if no such
 * object exists, or croaks, depending on the value of the 
 * croak_on_null flag.
 */
void *SvFluxPtr( SV *sv, char *name, char *tstr, char wflag, char croak_on_null ) {
  char errbuf[160];
  MAGIC *m;
  WORLD *w;
  long type;
  long label;

  if( SvROK(sv) && sv_derived_from(sv, tstr)) {
    sv = SvRV(sv);
  } else {
    sprintf(errbuf,"%s requires a %s.\n",name,tstr);
    croak(errbuf);
  }

  // Unpack tied hashes
  if( SvTYPE(sv)==SVt_PVHV && (m=mg_find(sv,'P'))) // assignment in second term
    sv = m->mg_obj;

  // Underlying object should now be a list ref.  Dereference it.	
  if(SvROK(sv))
   sv = SvRV(sv);
 
  // Now we should have an AV.  
  if( SvTYPE(sv) == SVt_PVAV ) {
    SV **svp;
    void *out;

    w     = (WORLD *)( (svp = av_fetch((AV *)sv, 0, 0)) ? SvIV(*svp) : 0 );
 
    if( !w ) {
	if(croak_on_null) {
		sprintf(errbuf, "SvFluxPtr (%s): got a null WORLD!",name);
		croak(errbuf);
	} else {
		return (void *)0;
	}
    }

    if(wflag)
	return w;

    type  =          ( (svp = av_fetch((AV *)sv, 1, 0)) ? SvIV(*svp) : 0 );
    label =          ( (svp = av_fetch((AV *)sv, 2, 0)) ? SvIV(*svp) : 0 );

    switch(type) {
	case FT_VERTEX: 
		out = tree_find( w->vertices, label, v_lab_of,  v_ln_of ); break;
	case FT_FLUXON:
		out = tree_find( w->lines,    label, fl_lab_of, fl_all_ln_of ); break;
	case FT_WORLD:
		out = w; break;
	case FT_CONC:
		out = tree_find( w->concentrations, label, fc_lab_of, fc_ln_of ); break;
	default:
		sprintf(errbuf,"SvFluxPtr (%s): unknown type code %d!",name,type);
		croak(errbuf);
		break;
    }
  
    if(!out && croak_on_null) {
	sprintf(errbuf,"SvFluxPtr (%s): Failed to find a %s (type %d) with label %d!",name,classnames[type],type,label);
	croak(errbuf);
    }
    return out;
  } else {
    sprintf(errbuf,"SvFluxPtr (%s): big trouble - couldn't find the expected deep list ref!",name);
    croak(errbuf);
  }
}

/**********************************************************************
 * SvLabel - return the label field of a FLUX object
 */
long SvLabel( SV *sv, char *name, char *tstr ) {
  char errbuf[160];
  MAGIC *m;
  long label;

  if( SvROK(sv) && sv_derived_from(sv, tstr)) {
    sv = SvRV(sv);
  } else {
    sprintf(errbuf,"%s requires a %s.\n",name,tstr);
    croak(errbuf);
  }

  // Unpack tied hashes
  if( SvTYPE(sv)==SVt_PVHV && (m=mg_find(sv,'P'))) // assignment in second term
    sv = m->mg_obj;

  // Underlying object should now be a list ref.  Dereference it.	
  if(SvROK(sv))
   sv = SvRV(sv);
 
  // Now we should have an AV.  
  if( SvTYPE(sv) == SVt_PVAV ) {
    SV **svp;
    void *out;

    label =          ( (svp = av_fetch((AV *)sv, 2, 0)) ? SvIV(*svp) : 0 );
    return label;

  } else {
    sprintf(errbuf,"SvLabel (%s): big trouble - couldn't find the expected deep list ref!",name);
    croak(errbuf);
  }
}

/**********************************************************************
 * new_sv_from_ptr - construct a generic FLUX object from a world pointer, 
 * type code, and label.
 * 
 * This generates an AV containing the relevant information, blesses it into
 * the appropriate class, and then ties it to a hash.  The final return value
 * is a ref to the tied hash, blessed into the same class.
 *
 * The WORLD's refcount is incremented along the way.
 */
SV *new_sv_from_ptr( WORLD *wptr, int type, long label ) {
  SV *out;
  SV *rv;
  AV *av;
  HV *hv;
  HV *stash;
  char errbuf[80];

  if(type < MIN_FT || type > MAX_FT) {
	croak("new_sv_from_ptr: unknown type %d\n",type);
  }

  av = newAV();
  av_extend(av,2);
  av_store( av, 0, newSViv( (IV)wptr ));
  av_store( av, 1, newSViv(  type    ));
  av_store( av, 2, newSViv(  label   ));
  
  rv = newRV_noinc((SV *)av);
  stash = gv_stashpv( classnames[type], TRUE );
  sv_bless((SV *)rv, stash);
  
  hv = newHV();
  hv_magic(hv, (GV*)rv, PERL_MAGIC_tied);

  out = newRV_noinc(hv);
  sv_bless(out,stash);
  wptr->refct++;
  
  return out;
}

/**********************************************************************
 * destroy_sv - destroy an SV created by new_sv_from_ptr.
 * 
 * The WORLD's refcount is decremented, and if it reaches zero the 
 * WORLD gets deallocated.
 * 
 * Tying the inner list ref to the hash generates an additional reference,
 * so we decrement the reference count to that ref; that should make it
 * evaporate. 
 *
 * We get called twice - once with the hash and once with the inner list ref - 
 * so we ignore the list ref call.
 */
void destroy_sv(SV *sv) {
 MAGIC *m;
 WORLD *w;
 SV *av;
 SV *rv;

 // Only do anything if we get a valid FLUX object to unwrap: a hash ref
 // whose hash is tied to a list ref.  Then decrement the WORLD refcount,
 // unmagick the hash, and decrement the refcount on the underlying list ref.
 if( SvROK(sv)                          &&
     (sv = SvRV(sv))             	&& // assignment
     (SvTYPE(sv)==SVt_PVHV)       	&& // test
     (m = mg_find(sv, PERL_MAGIC_tied)) && // assignment
     (rv = m->mg_obj )           	&& // assignment
     (SvROK(rv))                   	&&	
     (av=SvRV(rv))                 	&& // assignment
     (SvTYPE(av)==SVt_PVAV)    	  ) 	   // test
    {
        SV **svp;

	// Fetch the world and decrement its refcount.
	w = (WORLD *)( (svp = av_fetch((AV *)av, 0, 0)) ? SvIV(*svp) : 0 );
	if( w ) {
	  w->refct--;
	  if(!w->refct)
	    free_world(w);
	}
	
	// Unmagick the tied hash and decrement the refcount 
	// on the underlying list ref.  (This will result in 
 	// us being invoked again with the list ref, but the 
  	// initial if() will fall through.)
	sv_unmagic(sv, PERL_MAGIC_tied);
	SvREFCNT_dec(rv);
   }
}


MODULE = Flux::Core      PACKAGE = Flux::Core

BOOT:
/**********************************************************************/
/**********************************************************************/
/**** bootstrap code -- load-time dynamic linking to pre-loaded PDL   */
/**** modules and core functions.   **/
 perl_require_pv("PDL::Core");
 CoreSV = perl_get_sv("PDL::SHARE",FALSE);
 if(CoreSV==NULL)     Perl_croak(aTHX_ "Can't load PDL::Core module (required by Flux::Fluxon)");
 PDL = INT2PTR(Core*, SvIV( CoreSV ));  /* Core* value */
 if (PDL->Version != PDL_CORE_VERSION)
    Perl_croak(aTHX_ "Flux::Fluxon needs to be recompiled against the newly installed PDL");
/**********************************************************************/
/**********************************************************************/
/**** FLUX bootstrap code -- allocate and populate the FluxCore struct. */
/****/
FLUX = &FLUX_MASTER;
FLUX->CoreVersion               = FLUX_CORE_VERSION;
FLUX->new_label  		= new_label;
FLUX->new_fluxon      		= new_fluxon;
FLUX->new_vertex      		= new_vertex;
FLUX->add_vertex_pos    	= add_vertex_pos;
FLUX->add_vertex_after		= add_vertex_after;
FLUX->new_flux_concentration	= new_flux_concentration;
FLUX->new_world			= new_world;
FLUX->free_world		= free_world;
FLUX->delete_flux_concentration	= delete_flux_concentration;
FLUX->delete_fluxon		= delete_fluxon;
FLUX->delete_vertex		= delete_vertex;
FLUX->clear_links		= clear_links;
FLUX->tree_top			= tree_top;
FLUX->tree_find			= tree_find;
FLUX->tree_insert		= tree_insert;
FLUX->tree_binsert		= tree_binsert;
FLUX->tree_unlink		= tree_unlink;
FLUX->tree_balance		= tree_balance;
FLUX->tree_balance_check	= tree_balance_check;
FLUX->tree_walk			= tree_walk;
FLUX->tree_walker		= tree_walker;
FLUX->new_dumblist		= new_dumblist;
FLUX->free_dumblist		= free_dumblist;
FLUX->dumblist_init		= dumblist_init;
FLUX->dumblist_add		= dumblist_add;
FLUX->dumblist_rm		= dumblist_rm;
FLUX->dumblist_sort		= dumblist_sort;
FLUX->dumblist_snarf		= dumblist_snarf;
FLUX->dumblist_grow		= dumblist_grow;
FLUX->dumblist_clear		= dumblist_clear;
FLUX->dumblist_crunch		= dumblist_crunch;
FLUX->fl_eq			= fl_eq;
FLUX->next_line			= next_line;
FLUX->footpoint_action		= footpoint_action;
FLUX->fprint_tree		= fprint_tree;
FLUX->fprint_node		= fprint_node;
FLUX->print_tree		= print_tree;
FLUX->fdump_fluxon		= fdump_fluxon;
FLUX->fprint_all_fluxon_node	= fprint_all_fluxon_node;
FLUX->fprint_all_fluxon_tree	= fprint_all_fluxon_tree;
FLUX->print_all_fluxon_tree	= print_all_fluxon_tree;
FLUX->fprint_v_nbors		= fprint_v_nbors;
FLUX->fprint_fc_line		= fprint_fc_line;
FLUX->fprint_fc_line_nonneg	= fprint_fc_line_nonneg;
FLUX->fprint_fls_by_fc		= fprint_fls_by_fc;
FLUX->fprint_fl_vertices	= fprint_fl_vertices;
FLUX->print_dumblist		= print_dumblist;
FLUX->print_world		= print_world;
FLUX->fprint_world		= fprint_world;
FLUX->read_world		= read_world;
FLUX->world_update_neighbors	= world_update_neighbors;
FLUX->world_update_mag		= world_update_mag;
FLUX->world_relax_step		= world_relax_step;
FLUX->fluxon_update_neighbors	= fluxon_update_neighbors;
FLUX->fluxon_update_mag		= fluxon_update_mag;
FLUX->fluxon_calc_step		= fluxon_calc_step;
FLUX->fluxon_relax_step		= fluxon_relax_step;
FLUX->vertex_update_neighbors	= vertex_update_neighbors;
FLUX->gather_neighbor_candidates= gather_neighbor_candidates;
FLUX->winnow_cmp_1		= winnow_cmp_1;
FLUX->winnow_neighbor_candidates= winnow_neighbor_candidates;
FLUX->hull_neighbors		= hull_neighbors;
FLUX->fix_proximity		= fix_proximity;
FLUX->fluxon_fix_proximity	= fluxon_fix_proximity;
FLUX->global_fix_proximity	= global_fix_proximity;
FLUX->fix_curvature		= fix_curvature;
FLUX->fluxon_fix_curvature	= fluxon_fix_curvature;
FLUX->global_fix_curvature	= global_fix_curvature;
FLUX->reconnect_vertices	= reconnect_vertices;
FLUX->vertex_recon_check	= vertex_recon_check;
FLUX->fluxon_recon_check	= fluxon_recon_check;
FLUX->global_recon_check	= global_recon_check;
FLUX->fc_cancel                 = fc_cancel;
FLUX->world_collect_stats	= world_collect_stats;
FLUX->fluxon_collect_stats	= fluxon_collect_stats;
FLUX->norm_2d			= norm_2d;
FLUX->norm_3d			= norm_3d;
FLUX->norm2_2d			= norm2_2d;
FLUX->norm2_3d			= norm2_3d;
FLUX->inner_2d			= inner_2d;
FLUX->inner_3d			= inner_3d;
FLUX->cross_2d			= cross_2d;
FLUX->cross_3d			= cross_3d;
FLUX->scale_2d			= scale_2d;
FLUX->scale_3d			= scale_3d;
FLUX->sum_3d			= sum_3d;
FLUX->diff_3d			= diff_3d;
FLUX->cp_3d			= cp_3d;
FLUX->rotmat_2d			= rotmat_2d;
FLUX->rotmat_2d_fr_slope	= rotmat_2d_fr_slope;
FLUX->mat_mult_2d		= mat_mult_2d;
FLUX->transpose_2x2		= transpose_2x2;
FLUX->transpose_3x3		= transpose_3x3;
FLUX->mat_mult_3d		= mat_mult_3d;
FLUX->mat_vmult_2d		= mat_vmult_2d;
FLUX->mat_vmult_3d		= mat_vmult_3d;
FLUX->vec_mmult_3d		= vec_mmult_3d;
FLUX->det_2d			= det_2d;
FLUX->det_3d			= det_3d;
FLUX->points2plane		= points2plane;
FLUX->p_l_intersection		= p_l_intersection;
FLUX->xy_l_intersection		= xy_l_intersection;
FLUX->p_inside_tri		= p_inside_tri;
FLUX->trivloop			= trivloop;
FLUX->projmatrix		= projmatrix;
FLUX->cart2_2d			= cart2_2d;
FLUX->cart2_3d			= cart2_3d;
FLUX->cart_2d			= cart_2d;
FLUX->cart_3d			= cart_3d;
FLUX->p_l_dist			= p_l_dist;
FLUX->p_ls_dist			= p_ls_dist;
FLUX->l_l_dist			= l_l_dist;
FLUX->p_ls_closest_approach	= p_ls_closest_approach;
FLUX->ls_closest_approach	= ls_closest_approach;
FLUX->ls_ls_dist		= ls_ls_dist;
FLUX->fl_segment_masked_dist		= fl_segment_masked_dist;
FLUX->fl_segment_masked_deluxe_dist	= fl_segment_masked_deluxe_dist;
FLUX->fl_segment_dist			= fl_segment_dist;
FLUX->fl_segment_deluxe_dist	= fl_segment_deluxe_dist;
FLUX->perp_bisector_2d		= perp_bisector_2d;
FLUX->intersection_2d		= intersection_2d;
FLUX->project_n_fill            = project_n_fill;
FLUX->hull_2d			= hull_2d;
FLUX->hull_2d_us		= hull_2d_us;
FLUX->code_info_model           = code_info_model;
FLUX->code_info_physics         = code_info_physics;
FLUX->code_info_io  	        = code_info_io;
FLUX->code_info_geometry        = code_info_geometry;
FLUX->code_info_data            = code_info_data;
FLUX->F_B_NAMES                 = F_B_NAMES;
FLUX->FLUX_RECON		= FLUX_RECON;
FLUX->FLUX_FORCES		= FLUX_FORCES;
FLUX->SvFluxPtr			= SvFluxPtr;
FLUX->SvLabel			= SvLabel;
FLUX->new_sv_from_ptr		= new_sv_from_ptr;
FLUX->destroy_sv		= destroy_sv;
sv_setiv(perl_get_sv("Flux::Core::FLUX",TRUE), PTR2IV(FLUX));


	
	