/* Convenient perl interface macros for the perl side of FLUX.
 *
 * You need to include the perl api and FLUX api stuff yourself.
 * This is just a couple of convenient macros for converting a perl
 * SV to a WORLD *, FLUXON *, FLUX_CONCENTRATION *, or VERTEX * as 
 * appropriate.  
 * 
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Southwest Research Institute, 2004
 * 
 * You may modify and/or distribute this software under the temrs of
 * the Gnu Public License, version 2.  You should have received a copy
 * of the license with this software, in the file COPYING to be found
 * in the top level directory of the distribution.  You may obtain
 * additional copies of the licesnse via http://www.gnu.org or by
 * writing to the Free Software Foundation, 59 Temple Place - Suite
 * 330, Boston, MA 02111-1307 USA.
 *
 * The software comes with NO WARRANTY.
 * 
 * You may direct questions, kudos, gripes, and/or patches to the
 * author, Craig DeForest, at "deforest@boulder.swri.edu".
 * 
 * This is version 1.0 of fluxperl.h - part of the FLUX 1.0 release.
 */



/******************************
 * SvFluxPtr(wsv, name, type) 
 *
 * wsv: an SV containing a reference to one of the FLUX types.  
 *   It should be a blessed reference pointing to either a scalar
 *   or a tied has that is tied to a scalar ref.
 *   The underlying scalar should contain the pointer.
 * name: a string containing caller info (in case of failure)
 * type: a string containing the class name needed for the dereference
 *
 * Returns: a (void *) with the same numeric value as the IV in the scalar
 */
static void *SvFluxPtr( SV *sv, char *name, char *tstr ) {
  MAGIC *m;
  
  if(SvROK(sv) && sv_derived_from(sv,tstr)) {
    sv = SvRV(sv);
  } else {
    croak("%s requires a %s.\n",(name,tstr));
  }

  /* Notice if this is a magical tied hash */
  if(  SvTYPE(sv)==SVt_PVHV  &&  (m=mg_find(sv,'P'))){  // assignment in second term
    /* It's a magical tied hash - dereference the underlying scalar ref instead */
    sv=m->mg_obj;
    if(SvROK(sv))
      sv=SvRV(sv);
  }

  return (void *)SvIV(sv);
}
  
/******************************
 * SvVertex
 * SvConc
 * SvFluxon
 * SvWorld
 *
 * Each macro calls SvFluxPtr to pull out the underlying pointer, and casts it correctly.
 */

#define SvVertex(sv,name) ( (VERTEX *)            (SvFluxPtr((sv),(name),"Flux::Vertex")       ) )
#define SvConc(sv,name)   ( (FLUX_CONCENTRATION *)(SvFluxPtr((sv),(name),"Flux::Concentration")) )
#define SvFluxon(sv,name) ( (FLUXON *)            (SvFluxPtr((sv),(name),"Flux::Fluxon")       ) )
#define SvWorld(sv,name)  ( (WORLD *)             (SvFluxPtr((sv),(name),"Flux::World")        ) )
