/* Includes for FLUX -- field line and tie point structures 
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Craig DeForest, 2004-2008
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
 * This file is part of FLUX 2.2 (22-Nov-2008).
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
#include <float.h>
#define FP_EPSILON DBL_EPSILON /*change to FLT_EPSILON or LDBL_EPSILON if NUM is changed */

#include <stdlib.h>
#include <stddef.h>

/* MD5 hashing is now used for new VERTEX ID's CED 10-Feb-2009 */
#include <time.h>
#include <unistd.h>

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
  long passno;  /* see geometry.h, geometry.c (anti-repeat gate) */
  NUM r,a;      

     /* Magnetic field info is calculated by b_vec in physics.c.   */
     /* Remember that it is segment-centered, not vertex-centered. */
  POINT3D b_vec;                /* Stores the local B vector */
  NUM     b_mag;                /* Stores the magnitude of the B field */
  NUM energy;                            /* calculated B energy for segment */

     /* Force accumulators are filled by the force functions in physics.c. */
     /* As of v3, the vertex pseudoforce can be tracked indepenently of the */
     /* plasma force that is localized at the vertex.  For the purpose of relaxation */
     /* step sizing, the f_v_ps magnitude is added to f_v_tot, not to a separate variable. */
  POINT3D f_v;                  /* Stores force on the vertex */
  POINT3D f_v_ps;               /* vertex pseudoforce - used for vertex distribution */
  POINT3D f_s;                  /* Stores force on the following segment */
  POINT3D f_t;                  /* Stores total force */
  NUM f_s_tot;                  /* sum-of-magnitudes force on the segment (for help in auto-step-size)*/
  NUM f_v_tot;                  /* sum-of-magnitudes force on the vertex  (for help in auto-step-size)*/
  NUM f_n_tot;                  /* norm of the total force */


  NUM r_v, r_s;                 /* Stores relevant lengthscale for vertex & segment forces */
  NUM r_cl;                     /* Stores closest-approach radius of any vertex */
  NUM r_ncl;                    /* Stores closest-approach radius of neighbors */

  long label;                  
  LINKS world_links;      /* tree access to all the vertices in the world */
  POINT3D plan_step;       /* planned step */


  /* Single-fluid data -- these are valid halfway down the following segment. */
  NUM rho;         // Density (measured in particles (electrons) per unit volume)
  NUM T;           // Internal energy (sort-of-Kelvin; see conversion in World object)
  POINT3D p;       // Momentum (particle-distances per unit time)
  NUM A;           // Cross-sectional area (for conservation of stuff)

} VERTEX;

#define V_ISDUMMY(v) ( ((v)->label <= -1 && (v)->label >=-4) || (!(v)->line) || !((v)->line->label) )

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

#define FC_OB_ID  -1
#define FC_OE_ID  -2
#define FC_PB_ID  -3
#define FC_PE_ID  -4
#define FC_IM10_ID -8
#define FC_IM11_ID -9
#define FC_IM20_ID -10
#define FC_IM21_ID -11
#define FC_ISDUMMY(fc) ( (fc)->label < 0 && (fc)->label >= -11 )

#define FL_IM1_ID -9
#define FL_IM2_ID -10
#define FL_ISDUMMY(fl) ( (fl)->label == -9 || (fl)->label == -10 )

typedef struct FLUX_CONCENTRATION {
  struct WORLD *world;
  NUM flux;                  /* Total magnetic flux in Maxwells */
  long label;                /* Unique ID assigned at creation */
  struct FLUXON *lines;      /* Tree of lines coming out */
  struct LINKS links;        /* Tree for other concentrations */
  NUM x[3];                  /* Location of the concentration */
  NUM locale_radius;         /* Radius of the concentration's neighborhood */
  
  long passno;               /* Multiple-pass detection for hull calculation */
  void (*bound)(VERTEX *v);  /* Boundary condition routine */
  void (*pbound)(VERTEX *v); /* Plasma boundary condition routine */
  
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

#define MAXNUMCOEFFS 1024 /* Number of coefficients to create */

#define N_FORCE_FUNCS 30    /* Number of functions allowed in the force list */
#define N_RECON_FUNCS 10    /* Number of reconnection condition functions in the recon list */
#define N_RECON_PARAMS 8    /* Number of world-global parameters to be stored for each reconnection function */

#define N_M_PARAMS 5        /* Number of world-global parameters for a generic mass scaling function */

#define N_NEIGHBORS_PREALLOC 10 /* Number of neighbors to pre-allocate for each vertex.  Saves time on file read. */

typedef VERTEX *((RC_FUNC)(VERTEX *v, NUM *params));

typedef struct WORLD {
  /* Globally memorable fields go here */
  long frame_number;        /* Identifies currently-being-worked-on time step */
  long state;               /* State of the world  */
  long refct;               /* Used for perl interface (deallocate world when this reaches zero) */
  FLUX_CONCENTRATION *concentrations;
  FLUXON *lines;
  VERTEX *vertices;

  PHOTOSPHERE photosphere;
  PHOTOSPHERE photosphere2;

  FLUX_CONCENTRATION *fc_im10,*fc_im11,*fc_im20,*fc_im21;
  FLUXON *fl_im,*fl_im2;
  VERTEX *image, *image2,*image3 ,*image4;

  NUM locale_radius;        /* Default radius for concentrations' neighborhoods */

  long auto_open;             // Flag indicating whether a full-auto open boundary is being maintained
  
  FLUX_CONCENTRATION *fc_ob; /* bogus flux conc. to store open fluxons     */
  FLUX_CONCENTRATION *fc_oe; /* bogus flux conc. to store open fluxons     */
  FLUX_CONCENTRATION *fc_pb; /* bogus flux conc. to store plasmoid fluxons */
  FLUX_CONCENTRATION *fc_pe; /* bogus flux conc. to store plasmoid fluxons */

  void (*default_bound)(VERTEX *v);  /* Boundary condition routine */

  long verbosity;           /* Verbose flag turns on/off debugging lines */
  
  void ((*(f_funcs[N_FORCE_FUNCS]))());         /* force functions (incl. B field etc.) */
  void ((*(m_funcs[N_FORCE_FUNCS]))());         /* mass calculation functions */
  RC_FUNC *rc_funcs[N_RECON_FUNCS];             /* reconnection threshold detectors */
  NUM rc_params[N_RECON_FUNCS][N_RECON_PARAMS]; /* reconnection threshold parameters */

  NUM m_params[N_M_PARAMS];
  
  struct { 
    NUM b_power; 
    NUM d_power;
    NUM s_power;
    NUM ds_power;
  } step_scale;

  int passno;
  int handle_skew; // flag used to tweak Voronoi cell wall projection; see geometry.c

  // ARD added code to keep largest angle between fluxels and mean angle
  // between fluxels for WORLD 

  NUM max_angle;
  NUM mean_angle;
  // ARD - world should keep step and dt 
  NUM dtau;
  long rel_step;
  // ARD - array of coefficients 
  NUM coeffs[MAXNUMCOEFFS];
  long n_coeffs;
  long maxn_coeffs;

  // CED - how many separate processes should be run at once when calculating forces?
  long concurrency; 

  // Force and angle accumulators for global statistics;
  NUM f_min;
  NUM f_max;
  NUM fr_min;
  NUM fr_max;
  NUM ca_min;
  NUM ca_max;
  NUM ca_acc;
  NUM ca_ct;

  // Switch use of mass
  long use_fluid;
  NUM k_b; /* Boltzmann's constant in (n_e dist^2 / time^2) / K (settable but initialized) */
  long gravity_type;
  POINT3D gravity_origin;
  NUM g;   /* space / time^2  OR  space^3/time^2 */

  // temporary dumblist used for I/O (vertex number translation)
  DUMBLIST *tmplist;


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

long new_label(long request);
long new_vertex_label(long request);
long hash_vertex_label(long request, WORLD *w);

/* Offsets defined here... */
static const long fl_lab_of      = offsetof(FLUXON,label);
static const long fl_all_ln_of   = offsetof(FLUXON, all_links);
static const long fl_start_ln_of = offsetof(FLUXON, start_links);
static const long fl_end_ln_of   = offsetof(FLUXON, end_links);

static const long fc_lab_of      = offsetof(FLUX_CONCENTRATION, label);
static const long fc_ln_of       = offsetof(FLUX_CONCENTRATION, links);

static const long tree_lab_of = offsetof(TREE, node);
static const long tree_ln_of  = offsetof(TREE, links);

static const long v_neighbor_of = offsetof(VERTEX, neighbors);
static const long v_nearby_of   = offsetof(VERTEX, nearby);
static const long v_ln_of       = offsetof(VERTEX, world_links);
static const long v_lab_of      = offsetof(VERTEX, label);

 FLUXON *new_fluxon(NUM flux,  
			 FLUX_CONCENTRATION *c0,
			 FLUX_CONCENTRATION *c1,
			 long label,
			 char plasmoid
			 );

 VERTEX *new_vertex(long label, NUM x, NUM y, NUM z,
			      FLUXON *fluxon);

 VERTEX *spawn_new_vertex(NUM x, NUM y, NUM z, VERTEX *parent);

       int add_vertex_pos     (FLUXON *fl, long pos,  VERTEX *v);
 int add_vertex_after   (FLUXON *fl, VERTEX *nbr,  VERTEX *v);



 FLUX_CONCENTRATION *new_flux_concentration(WORLD *, NUM x, NUM y, NUM z,
						  NUM flux,
						  long label
						  );
void check_special_concentration(FLUX_CONCENTRATION *f);

WORLD *new_world();

void free_world( WORLD *w );
void delete_flux_concentration(FLUX_CONCENTRATION *fc);
void delete_fluxon(FLUXON *f);
void delete_vertex(VERTEX *v);
void unlink_vertex(VERTEX *v);
void vertex_clear_neighbors(VERTEX *v);
void vertex_add_neighbor(VERTEX *v, VERTEX *neighbor);
int vertex_renumber(VERTEX *v, long newlab);
int fluxon_renumber(FLUXON *f, long newlab);
int concentration_renumber(FLUX_CONCENTRATION *fc, long newlab);

#define MAX_TREE_DEPTH 1000000
void clear_links(LINKS *links);
 void *tree_top(void *tree, int link_offset);
 void *tree_find(void *tree, long label, int label_offset, int link_offset);
 void *tree_insert(void *root, void *item, int label_offset, int link_offset);
 void *tree_binsert(void *root, void *item, int label_offset, int link_offset);
 void *tree_unlink(void *data, int label_offset, int link_offset);
 void *tree_bunlink(void *data, int label_offset, int link_offset);
 void *tree_balance(void *tree, int label_offset, int link_offset);
 char tree_balance_check(void *tree, int link_offset);

 long tree_walk(void *tree, int label_ofset, int link_offset, long ((*func)()));
long tree_walker(void *tree, int label_ofset, int link_offset, long ((*func)()),int depth);

long stw_helper(void *node, int lab, int link, int depth);
long safe_tree_walker(void *tree, int label_offset, int link_offset, long((*func)()), int depth);


DUMBLIST *new_dumblist();
void free_dumblist(DUMBLIST *dl);
void dumblist_init(DUMBLIST *foo);
void dumblist_clean(DUMBLIST *foo);
void dumblist_add(DUMBLIST *dl, void *a);
void dumblist_quickadd(DUMBLIST *dl, void *a);
void dumblist_delete(DUMBLIST *dl, void *a);
void dumblist_rm(DUMBLIST *dl, int i);
void dumblist_sort(DUMBLIST *dl, int((*compare)(void *a, void *b)));
void dumblist_snarf(DUMBLIST *dl, DUMBLIST *source);
void dumblist_grow(DUMBLIST *dl, int size);
void dumblist_clear(DUMBLIST *dl);
void dumblist_crunch(DUMBLIST *dl);

char fl_eq(NUM a, NUM b);


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



extern char *malloc_types[MALLOC_MAXTYPENO+1];

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
