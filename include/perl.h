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
 */


/******************************
 * SvWorld(wsv, name) 
 *
 * wsv: an SV containing a reference to a WORLD.
 * name: the name of the calling routine, in case of failure
 *
 * Returns: the WORLD * (croaks on error).
 */
#define SvWorld(wsv,name) (   ( SvROK(wsv) && sv_derived_from((wsv),"Flux::World") )  ? (WORLD *)(SvIVx(SvRV(wsv))) : ( croak("%s requires a Flux::World.\n",(name)), (WORLD *)0 ) )


/******************************
 * SvFluxon(wsv, name) 
 *
 * wsv: an SV containing a reference to a FLUXON.
 * name: the name of the calling routine, in case of failure
 *
 * Returns: the FLUXON * (croaks on error).
 */
#define SvFluxon(wsv,name) (   ( SvROK(wsv) && sv_derived_from((wsv),"Flux::Fluxon") )  ? (FLUXON *)(SvIVx(SvRV(wsv))) : ( croak("%s requires a Flux::Fluxon.\n",(name)), (FLUXON *)0 ) )

/******************************
 * SvConc(wsv, name) 
 *
 * wsv: an SV containing a reference to a FLUX_CONCENTRATION.
 * name: the name of the calling routine, in case of failure
 *
 * Returns: the FLUX_CONCENTRATION * (croaks on error).
 */
#define SvConc(wsv,name) (   ( SvROK(wsv) && sv_derived_from((wsv),"Flux::Concentration") )  ? (FLUX_CONCENTRATION *)(SvIVx(SvRV(wsv))) : ( croak("%s requires a Flux::Concentration.\n",(name)), (FLUX_CONCENTRATION *)0 ) )


/******************************
 * SvVertex(wsv, name) 
 *
 * wsv: an SV containing a reference to a VERTEX.
 * name: the name of the calling routine, in case of failure
 *
 * Returns: the VERTEX * (croaks on error).
 */
#define SvVertex(wsv,name) (   ( SvROK(wsv) && sv_derived_from((wsv),"Flux::Vertex") )  ? (VERTEX *)(SvIVx(SvRV(wsv))) : ( croak("%s requires a Flux::Vertex.\n",(name)), (VERTEX *)0 ) )
