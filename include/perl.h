/* Convenient perl interface macros for the perl side of FLUX.
 *
 * You need to include the perl api and FLUX api stuff yourself.
 */


#define SvWorld(wsv,name) (   ( SvROK(wsv) && sv_derived_from((wsv),"Flux::World") )  ? (WORLD *)(SvIVx(SvRV(wsv))) : ( croak("%s requires a Flux::World.\n",(name)), (WORLD *)0 ) )

#define SvFluxon(wsv,name) (   ( SvROK(wsv) && sv_derived_from((wsv),"Flux::Fluxon") )  ? (FLUXON *)(SvIVx(SvRV(wsv))) : ( croak("%s requires a Flux::Fluxon.\n",(name)), (FLUXON *)0 ) )


#define SvConc(wsv,name) (   ( SvROK(wsv) && sv_derived_from((wsv),"Flux::Concentration") )  ? (FLUX_CONCENTRATION *)(SvIVx(SvRV(wsv))) : ( croak("%s requires a Flux::Concentration.\n",(name)), (FLUX_CONCENTRATION *)0 ) )

#define SvVertex(wsv,name) (   ( SvROK(wsv) && sv_derived_from((wsv),"Flux::Vertex") )  ? (VERTEX *)(SvIVx(SvRV(wsv))) : ( croak("%s requires a Flux::Vertex.\n",(name)), (VERTEX *)0 ) )


