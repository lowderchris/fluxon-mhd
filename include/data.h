/* Includes for FLUX -- field line and tie point structures 
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Southwest Research Institute, 2004
 * 
 * You may modify and/or distribute this software under the terms of
 * the Gnu Public License, version 2.  You should have received a copy
 * of the license with this software, in the file COPYING to be found
 * in the top level directory of the distribution.  You may obtain
 * additional copies of the license via http://www.gnu.org or by
 * writing to the Free Software Foundation, 59 Temple Place - Suite
 * 330, Boston, MA 02111-1307 USA.
 *
 * The software comes with NO WARRANTY.
 * 
 * You may direct questions, kudos, gripes, and/or patches to the
 * author, Craig DeForest, at "deforest@boulder.swri.edu".
 * 
 * This is version 1.1 of data.h - part of the FLUX 1.1 release
 */

#ifndef FLEM_DATA
#define FLEM_DATA 1

#define NUM double             /* Type of calculation.  Double is (in these  */
#define NUMDOUBLE 1            /* days of 64-bit processors) actually faster*/
                               /* than float, but uses twice as much memory!*/
/* NUM must be one of (float|double|longdouble).  Set NUMDOUBLE to (0|1|2),
 * and NUMCHAR to "", "l", or "L". [NUMCHAR modifies NUM parsing by printf --
 * the 'l' makes the parse a double.]
 */
#define NUMCHAR "l"

#include <stdlib.h>

/**********************************************************************
 * Data structures
 * 
 * Note that VERTEXs and FLUXONs both contain different sorts
 * of tree linking.  The trees are dealt with by a bunch of semi-generic
 * tree methods.  The tree methods keep track of a running sum of
 * both number of links and also magnetic flux.  The semi-genericity
 * comes from the variable offset from the struct address of the link
 * table used for the tree.  One of the running sums comes from the
 * magnetic flux (or any other NUM), which must be the very first item
 * in the struct.
 *
 * The VERTEXs and FLUXONs contain their own tree links to facilitate
 * moving between different data trees.  There is also a generic tree
 * construct which is used for passing around bunches of things.
 *
 */

typedef NUM POINT2D[2];
typedef NUM POINT3D[3];

typedef struct LINKS {
  void *left;
  void *right;
  void *up;
  long n;
  NUM sum;
  char balanced;
} LINKS;
  

typedef struct TREE {
  void *node;  /* node pointer serves as label! */
  LINKS links;
} TREE;
  

#ifndef NULL
#define NULL ((void *)0)
#endif

typedef struct DUMBLIST {
  void **stuff;
  unsigned int n;
  unsigned int size;
} DUMBLIST;

typedef struct VERTEX {
  NUM energy;                            /* calculated B energy for segment */
  struct FLUXON *line;                   /* Which line are we on? */
  struct VERTEX *prev,*next;  

     /* Vertex location in physical space */
  POINT3D x;           

     /* neighbors & nearby - infrastructure for geometrical work */
  DUMBLIST neighbors;
  DUMBLIST nearby;

     /* Scratch space used when this vertex is considered as a neighbor of
      * someone else; part of force calculation. */
  POINT3D scr; 		
  NUM r,a;      

     /* Magnetic field info is calculated by b_vec in physics.c.   */
     /* Remember that it is segment-centered, not vertex-centered. */
  POINT3D b_vec;                /* Stores the local B vector */
  NUM     b_mag;                /* Stores the magnitude of the B field */

     /* Force accumulators are filled by the force functions in physics.c. */
  POINT3D f_v;                  /* Stores force on the vertex */
  POINT3D f_s;                  /* Stores force on the following segment */
  POINT3D f_t;                  /* Stores total force */
  NUM f_s_tot;                  /* sum-of-magnitudes force on the segment (for help in auto-step-size)*/
  NUM f_v_tot;                  /* sum-of-magnitudes force on the vertex  (for help in auto-step-size)*/


  NUM r_v, r_s;                 /* Stores relevant lengthscale for vertex & segment forces */
  NUM r_cl;                     /* Stores closest-neighbor-approach radius */

  long label;                  
  LINKS world_links;      /* tree access to all the vertices in the world */
} VERTEX;

typedef struct FLUXON {
  NUM flux;                  /* Flux in Maxwells */
  long label;                /* Unique ID assigned at creation */
  

  struct VERTEX *start, *end; /* Endpoints of the field line (could be the

		                   same point for loops!!) */
  long v_ct;         /* Index counter -- should be the number
			 of VERTEXes in the fluxon.  */

  
  LINKS all_links;      /* structure containing all fluxons */
  LINKS start_links;    /* structure of fluxons sharing this one's start */
  LINKS end_links;      /* structure of fluxons sharing this one's end */
  struct FLUX_CONCENTRATION *fc0, *fc1;
  char plasmoid;

} FLUXON;

typedef struct FLUX_CONCENTRATION {
  struct WORLD *world;
  NUM flux;                  /* Total magnetic flux in Maxwells */
  long label;                /* Unique ID assigned at creation */
  struct FLUXON *lines;      /* Tree of lines coming out */
  struct LINKS links;        /* Tree for other concentrations */
  NUM x[3];                  /* Location of the concentration */
  NUM locale_radius;         /* Radius of the concentration's neighborhood */
  
  void (*bound)(VERTEX *v);  /* Boundary condition routine */
} FLUX_CONCENTRATION;

typedef struct PLANE {
  NUM origin[3];
  NUM normal[3];
} PLANE;

#define PHOT_NONE 0
#define PHOT_PLANE 1
#define PHOT_SPHERE 2
#define PHOT_CYL 3
extern char *BOUNDARY_NAMES[4];  /* in data.c */

typedef struct PHOTOSPHERE {
  PLANE *plane;
  int   type;
} PHOTOSPHERE;

#define WORLD_STATE_NEW     0
#define WORLD_STATE_LOADING 1
#define WORLD_STATE_LOADED  2
#define WORLD_STATE_WORKING 3
#define WORLD_STATE_RELAXED 4
#define WORLD_STATE_READY   5

#define N_FORCE_FUNCS 30    /* Number of functions allowed in the force list */
typedef struct WORLD {
  /* Globally memorable fields go here */
  long frame_number;        /* Identifies currently-being-worked-on frame */
  long state;               /* State of the world  */
  long refct;               /* Used for perl interface (deallocate world when this reaches zero) */
  FLUX_CONCENTRATION *concentrations;
  FLUXON *lines;
  VERTEX *vertices;

  PHOTOSPHERE photosphere;

  VERTEX *image, *image2;
  NUM locale_radius;        /* Default radius for concentrations' neighborhoods */
  
  FLUX_CONCENTRATION *fc_ob; /* bogus flux concentration to store open fluxons     */
  FLUX_CONCENTRATION *fc_oe; /* bogus flux concentration to store open fluxons     */
  FLUX_CONCENTRATION *fc_pb; /* bogus flux concentration to store plasmoid fluxons */
  FLUX_CONCENTRATION *fc_pe; /* bogus flux concentration to store plasmoid fluxons */

  long verbosity;           /* Verbose flag turns on/off debugging lines */
  
  void ((*(f_funcs[N_FORCE_FUNCS]))());
  
  struct { 
    NUM b_power; 
    NUM d_power;
    NUM s_power;
    NUM ds_power;
  } step_scale;

} WORLD;

const char *world_state_name(WORLD *a);

/**********************************************************************
 * Barf opcodes
 */

#define BARF_MALLOC 1

/**********************************************************************
 * Routines in data.c
 */

void barf(int barf_on_op, char *where);

long new_label(long label);

/* Offsets defined here... */
static const FLUXON f__samp;
static const fl_lab_of  = (long)&(f__samp.label) - (long)&f__samp;
static const fl_all_ln_of = (long)&(f__samp.all_links) - (long)&f__samp;
static const fl_start_ln_of = (long)&(f__samp.start_links) - (long)&f__samp;
static const fl_end_ln_of = (long)&(f__samp.end_links) - (long)&f__samp;

static const FLUX_CONCENTRATION fc__samp;
static const fc_lab_of = (long)&(fc__samp.label) - (long)&fc__samp;
static const fc_ln_of = (long)&(fc__samp.links) - (long)&fc__samp;

static const TREE tree__samp;
static const tree_lab_of = (long)&(tree__samp.node) - (long)&tree__samp;
static const tree_ln_of  = (long)&(tree__samp.links) - (long)&tree__samp;

static const VERTEX v__samp;
static const v_neighbor_of = (long)&(v__samp.neighbors) - (long)&(v__samp);
static const v_nearby_of = (long)&(v__samp.nearby) - (long)&(v__samp);
static const v_ln_of = (long)&(v__samp.world_links) - (long)&(v__samp);
static const v_lab_of = (long)&(v__samp.label) - (long)&(v__samp);
inline FLUXON *new_fluxon(NUM flux,  
			 FLUX_CONCENTRATION *c0,
			 FLUX_CONCENTRATION *c1,
			 long label,
			 char plasmoid
			 );

inline VERTEX *new_vertex(long label, NUM x, NUM y, NUM z,
			      FLUXON *fluxon);

inline VERTEX *spawn_new_vertex(NUM x, NUM y, NUM z, VERTEX *parent);

       int add_vertex_pos     (FLUXON *fl, long pos,  VERTEX *v);
inline int add_vertex_after   (FLUXON *fl, VERTEX *nbr,  VERTEX *v);



inline FLUX_CONCENTRATION *new_flux_concentration(WORLD *, NUM x, NUM y, NUM z,
						  NUM flux,
						  long label
						  );

WORLD *new_world();

void free_world();
void delete_flux_concentration();
void delete_fluxon();
void delete_vertex();

#define MAX_TREE_DEPTH 1000000
void clear_links(LINKS *links);
inline void *tree_top(void *tree, int link_offset);
inline void *tree_find(void *tree, long label, int label_offset, int link_offset);
inline void *tree_insert(void *root, void *item, int label_offset, int link_offset);
inline void *tree_binsert(void *root, void *item, int label_offset, int link_offset);
inline void *tree_unlink(void *data, int label_offset, int link_offset);
inline void *tree_bunlink(void *data, int label_offset, int link_offset);
inline void *tree_balance(void *tree, int label_offset, int link_offset);
inline char tree_balance_check(void *tree, int link_offset);

inline long tree_walk(void *tree, int label_ofset, int link_offset, long ((*func)()));
long tree_walker(void *tree, int label_ofset, int link_offset, long ((*func)()),int depth);

DUMBLIST *new_dumblist();
void free_dumblist(DUMBLIST *dl);
inline void dumblist_add(DUMBLIST *dl, void *a);
inline void dumblist_delete(DUMBLIST *dl, void *a);
inline void dumblist_rm(DUMBLIST *dl, int i);
inline void dumblist_sort(DUMBLIST *dl, int((*compare)(void *a, void *b)));
inline void dumblist_snarf(DUMBLIST *dl, DUMBLIST *source);
inline void dumblist_grow(DUMBLIST *dl, int size);
inline void dumblist_clear(DUMBLIST *dl);

inline char fl_eq(NUM a, NUM b);


/* defines for debugging malloc wrapper in data.c */
#define MALLOC_WORLD 1
#define MALLOC_FLUXON 2
#define MALLOC_VERTEX 3
#define MALLOC_FC 4
#define MALLOC_DL 5
#define MALLOC_DL_L 6
#define MALLOC_VL 7
#define MALLOC_DLS_ARENA 8
#define MALLOC_PLANE 9
#define MALLOC_MISC 10
#define MALLOC_MAXTYPENO 10


char *malloc_types[MALLOC_MAXTYPENO+1];


#define valid_malloc_type(x) ((x)>0 && (x)<=MALLOC_MAXTYPENO)


/****************************************
 * Memory handling options
 * We have the ability to use either vanilla malloc,
 * a debugging/fencing malloc, or a perl-SV-based malloc.
 * If libflux is linked to the PDL front end, then 
 * you should use the perl malloc.
 */

#ifdef USE_DEBUGGING_MALLOC

#define localmalloc(x,y) flux_malloc(x,y)
#define localfree(x) flux_free(x)

void flux_free(void *p);
char *flux_malloc(long size, int what_for);
void flux_memcheck();

#else
#ifdef USE_PERL_MALLOC

#define localmalloc(x,y) flux_perl_malloc(x)
#define localfree(x) flux_perl_free(x)

 char *flux_perl_malloc(long size);
 void flux_perl_free(void *where);

#include "EXTERN.h"
#include "perl.h"
#include "perlapi.h"

#else
#ifdef USE_PADDED_MALLOC

#define localmalloc(x,y) flux_padded_malloc(x)
#define localfree(x) flux_padded_free(x)

 char *flux_padded_malloc(long size);
 void flux_padded_free(void *ptr);
#else

/*** Default case: just use normal malloc ***/
 #define localmalloc(x,y) malloc(x)
 #define localfree(x) free(x)

#endif /**** use_padded_malloc test ***/
#endif /**** use_perl_malloc test ****/
#endif /**** malloc switch ****/



#endif /* overall file include */


/******************************
 *
 * version-tracking: informational strings about each file in use
 */
extern char *code_info_model;
extern char *code_info_physics;
extern char *code_info_io;
extern char *code_info_geometry;
extern char *code_info_data;

