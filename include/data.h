/* Data structures for flr -- field line and tie point structures */
/* */
/* See accompanying index.html file */
#ifndef FLEM_DATA
#define FLEM_DATA 1

#define NUM double             /* Type of calculation.  Double is (in these  */
#define NUMDOUBLE 1            /* days of 64-bit processors) actually faster*/
                               /* than float, but uses twice as much memory!*/
/* NUM must be one of (float|double|longdouble).  Set NUMDOUBLE to (0|1|2),
   and NUMCHAR to "", "l", or "L". */
#define NUMCHAR "l"

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
  

#define NULL ((void *)0)

typedef struct DUMBLIST {
  void **stuff;
  unsigned int n;
  unsigned int size;
} DUMBLIST;

typedef struct VERTEX {
  struct FLUXON *line;                      /* Which line are we on? */
  struct VERTEX *prev,*next;  /* List links        */
  POINT3D x;                                /* x,y,z location */
  POINT3D scr;                /* Scratch vector -- comes in handy for 
			       * generating neighbor lists and the like */
  DUMBLIST neighbors;
  DUMBLIST nearby;

  NUM r,a;                      /* Scratch radius and angle for neighbor stuff */

  POINT3D f_v;                  /* Stores force on the vertex */
  POINT3D f_s;                  /* Stores force on the following segment */
  POINT3D f_t;                  /* Stores total force */
  NUM f_s_tot;                  /* sum-of-magnitudes force on the segment (for help in auto-step-size)*/
  NUM f_v_tot;                  /* sum-of-magnitudes force on the vertex  (for help in auto-step-size)*/


  NUM r_v, r_s;                 /* Stores relevant lengthscale for vertex & segment forces */
  NUM r_cl;                     /* Stores closest-neighbor-approach radius */

  long label;                  
} VERTEX;

typedef struct FLUXON {
  NUM flux;                  /* Flux in Maxwells */
  long label;                /* Unique ID assigned at creation */
  

  struct VERTEX *start, *end; /* Endpoints of the field line (could be the

		                   same point for loops!!) */
  void (*start_b)(),(*end_b)(); /* Boundary condition routines */

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
} FLUX_CONCENTRATION;

typedef struct PLANE {
  NUM origin[3];
  NUM normal[3];
} PLANE;

#define WORLD_STATE_NEW     0
#define WORLD_STATE_LOADING 1
#define WORLD_STATE_LOADED  2
#define WORLD_STATE_WORKING 3
#define WORLD_STATE_RELAXED 4
#define WORLD_STATE_READY   5


typedef struct WORLD {
  /* Globally memorable fields go here */
  long frame_number;        /* Identifies currently-being-worked-on frame */
  int  state;               /* State of the world  */
  FLUX_CONCENTRATION *concentrations;
  FLUXON *lines;

  PLANE *photosphere;
  VERTEX *image, *image2;
  NUM locale_radius;        /* Default radius for concentrations' neighborhoods */
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
#endif /* overall file include */

