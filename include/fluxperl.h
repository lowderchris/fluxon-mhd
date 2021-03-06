/* Convenient perl interface macros for the perl side of FLUX.
 *
 * You need to include the perl api and FLUX api stuff yourself.
 * This is just a couple of convenient macros for converting a perl
 * SV to a WORLD *, FLUXON *, FLUX_CONCENTRATION *, or VERTEX * as 
 * appropriate.  
 * 
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Craig DeForest, 2004-2008
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
 * This file is part of FLUX 2.2 (22-Nov-2008).
 */

#define FC FLUX_CONCENTRATION
#define FLUX_CORE_VERSION 18

typedef struct FluxCore {

  long CoreVersion;

  /************  functions from data.h */
  long (*new_label)(long label);
  FLUXON *(*new_fluxon)     (NUM flux, FC *c0, FC *c1, long label, char plasmoid);
  VERTEX *(*new_vertex)     (long labl, NUM x, NUM y, NUM z, FLUXON *fluxon);
  int    (*add_vertex_pos)  (FLUXON *fl, long pos, VERTEX *v);
  int    (*add_vertex_after)(FLUXON *fl, VERTEX *nbr, VERTEX *v);
  FC     *(*new_flux_concentration)(WORLD *w, NUM x, NUM y, NUM z, NUM flux, long label);
  WORLD  *(*new_world)();
  void   (*free_world)();
  void   (*delete_flux_concentration)();
  void   (*delete_fluxon)();
  void   (*delete_vertex)();
  int    (*vertex_renumber)(VERTEX *v, long newlab);
  int    (*fluxon_renumber)(FLUXON *f, long newlab);
  int    (*concentration_renumber)(FLUX_CONCENTRATION *fc, long newlab);


  void (*clear_links)(LINKS *links);
  void *(*tree_top)          (void *tree, int link_offset);
  void *(*tree_find)         (void *tree, long label, int label_offset, int link_offset);
  void *(*tree_insert)       (void *root, void *item, int label_offset, int link_offset);
  void *(*tree_binsert)      (void *root, void *item, int label_offset, int link_offset);
  void *(*tree_unlink)       (void *data, int label_offset, int link_offset);
  void *(*tree_balance)      (void *tree, int label_offset, int link_offset);
  char  (*tree_balance_check)(void *tree, int link_offset);

  long (*tree_walk)  (void *tree, int lab_of, int ln_of, long ((*func)()));
  long (*tree_walker)(void *tree, int lab_of, int ln_of, long ((*func)()),int depth);

  long (*stw_helper)(void *node, int lab, int link, int depth);
  long (*safe_tree_walker)(void *tree, int lab_of, int ln_of, long ((*func)()),int depth);

  DUMBLIST *(*new_dumblist)();
  void (*free_dumblist)     (DUMBLIST *dl);
  void (*dumblist_init)     (DUMBLIST *dl);
  void (*dumblist_add)      (DUMBLIST *dl, void *a);
  void (*dumblist_delete)   (DUMBLIST *dl, void *a);
  void (*dumblist_rm)       (DUMBLIST *dl, int i);
  void (*dumblist_sort)     (DUMBLIST *dl, int ((*compare)(void *a, void *b)));
  void (*dumblist_snarf)    (DUMBLIST *dl, DUMBLIST *source);
  void (*dumblist_grow)     (DUMBLIST *dl, int size);
  void (*dumblist_clear)    (DUMBLIST *dl);
  void (*dumblist_crunch)   (DUMBLIST *dl);

  char (*fl_eq)(NUM a, NUM b);

  

  /************  functions from io.h  */
  char *(*next_line)(FILE *file);
  
  int (*footpoint_action)(WORLD *a, char *s);
  
  void (*fprint_tree)(FILE *f, void *t, int lb, int lk, int idt, void ((*prntr)()));
  void (*fprint_node)(FILE *f, void *foo, int indent, int lab_o, int lk_o);
  void (*print_tree)(void *t, int lab_o, int lk_o, int indent, void ((*printer)()));
  void (*fdump_fluxon)(FILE *f, FLUXON *foo, int indent);
  void (*fprint_all_fluxon_node)(FILE *f, FLUXON *foo, int indent);
  void (*fprint_all_fluxon_tree)(FILE *f, FLUXON *foo);
  void (*print_all_fluxon_tree)(FLUXON *foo);
  void (*fprint_v_nbors)(FILE *f, void *foo, int indnt, int lab_o, int ln_o);
  void (*fprint_fc_line)(FILE *f, void *foo, int indnt, int lab_o, int ln_o);
  void (*fprint_fc_line_nonneg)(FILE *f, void *foo, int indnt, int lab_o, int ln_o);
  void (*fprint_fls_by_fc)(FILE *f, void *foo, int indnt, int lab_o, int ln_o);
  void (*fprint_fl_vertices)(FILE *f, void *foo, int indnt, int lab_o, int ln_o);
  
  void (*print_dumblist)(DUMBLIST *foo, void ((*item_printer)()));

  int (*print_world)(WORLD *a, char*header);
  int (*fprint_world)(FILE *file, WORLD *world, char *header);

  WORLD *(*read_world)(FILE *file, WORLD *a);

  int (*binary_dump_header)(int fd);
  int (*binary_dump_WORLD)(int fd, WORLD *w);
  int (*binary_dump_CONCENTRATION)(int fd, FLUX_CONCENTRATION *fc);
  int (*binary_dump_FLUXON)(int fd, FLUXON *fl);
  int (*binary_dump_neighbors)(int fd, FLUXON *fl);
  int (*binary_dump_end)(int fd);
  int (*binary_dump_fluxon_pipe)(int fd, FLUXON *f);

  WORLD *(*binary_read_dumpfile)(int fd, WORLD *w);
  int (*binary_read_fluxon_pipe)(long size, char *buf, WORLD *w);
  
  /***********  functions from model.h */
  
  
  void (*world_update_neighbors) (WORLD *a, char global);
  int (*world_update_mag)       (WORLD *a, char global);
  int (*world_update_mag_parallel)(WORLD *a, char global);
  void (*world_fluxon_length_check) (WORLD *a, char global);  
  void (*world_relax_step)       (WORLD *a, NUM t);
  void (*world_relax_step_parallel) (WORLD *a, NUM t);

  int (*fluxon_update_neighbors)(FLUXON *fl, char gl);
  NUM (*fluxon_update_mag)      (FLUXON *fl, char gl, void ((**f_funcs)()));

  void (*fluxon_calc_step) (FLUXON *fl, NUM t);
  void (*fluxon_relax_step)(FLUXON *fl, NUM t);

  HULL_VERTEX *(*vertex_update_neighbors)(VERTEX *v, char global);

  DUMBLIST *(*gather_neighbor_candidates)(VERTEX *v,char global);
  void (*image_find)(PHOTOSPHERE *phot, VERTEX *image, VERTEX *v);

  int  (*winnow_cmp_1)(void *a, void *b);
  void (*winnow_neighbor_candidates)(VERTEX *v, DUMBLIST *horde);

  HULL_VERTEX *(*hull_neighbors)(VERTEX *v, DUMBLIST *horde);

  int (*fix_proximity)(VERTEX *V, NUM scale_thresh);
  int (*fluxon_fix_proximity)(FLUXON *F, NUM scale_thresh);
  int (*global_fix_proximity)(WORLD *w, NUM scale_thresh);

  int (*fix_curvature)(VERTEX *V, NUM curve_thresh_high, NUM curve_thresh_low);
  int (*fluxon_fix_curvature)(FLUXON *F, NUM curve_thresh_high, NUM curve_thresh_low);
  int (*global_fix_curvature)(WORLD *w, NUM curve_thresh_high, NUM curve_thresh_low);

  void (*reconnect_vertices)( VERTEX *v1, VERTEX *v2, long passno );
  int (*vertex_recon_check)( VERTEX *v1, long passno );
  long (*fluxon_recon_check)( FLUXON *f, long passno );
  long (*global_recon_check)( WORLD *w );

  int (*fc_cancel)( FLUX_CONCENTRATION *fc0, FLUX_CONCENTRATION *fc1 );

  VERTEX_STATS *(*world_collect_stats)(WORLD *a);
  void (*fluxon_collect_stats)(FLUXON *fl, VERTEX_STATS *st);

  void (*fluxon_auto_open)(FLUXON *f);
  DUMBLIST *(*gather_photosphere_neighbor_candidates)(VERTEX *v, char global);
  HULL_VERTEX *(*photosphere_hull_neighbors)(VERTEX *v, DUMBLIST *horde);
  HULL_VERTEX *(*photosphere_vertex_update_neighbors)(VERTEX *v, char global, int *n_ptr);

  /************  routines from geometry.h */
  NUM (*norm_2d)(NUM *x);
  NUM (*norm_3d)(NUM *x);
  NUM (*norm2_2d)(NUM *x);
  NUM (*norm2_3d)(NUM *x);
  NUM (*inner_2d)(NUM *p0, NUM *p1);
  NUM (*inner_3d)(NUM *p0, NUM *p1);
  NUM (*cross_2d)(NUM *p0, NUM *p1);
  void (*cross_3d)(NUM *out, NUM *p0, NUM *p1);
  void (*scale_2d)(NUM *out, NUM *a, NUM alpha);
  void (*scale_3d)(NUM *out, NUM *a, NUM alpha);
  void (*sum_3d)(NUM *out, NUM *a, NUM *b);
  void (*diff_3d)(NUM *out, NUM *a, NUM *b);
  void (*cp_3d)(NUM *out, NUM *a);
  
  void (*rotmat_2d)(NUM *out, NUM alpha);
  void (*rotmat_2d_fr_slope)(NUM *out, NUM dy, NUM dx);
  void (*mat_mult_2d)(NUM *out, NUM *a, NUM *b);
  void (*transpose_2x2)(NUM *mat);
  void (*transpose_3x3)(NUM *mat);
  void (*mat_mult_3d)(NUM *out, NUM *a, NUM *b);
  void (*mat_vmult_2d)(NUM *out, NUM *mat, NUM *v);
  void (*mat_vmult_3d)(NUM *out, NUM *mat, NUM *v);
  void (*vec_mmult_3d)(NUM *out, NUM *mat, NUM *v);
  NUM (*det_2d)(NUM *mat);
  NUM (*det_3d)(NUM *mat);

  void (*points2plane)(PLANE *plane, NUM *p0, NUM *p1, NUM *p2);
  int (*p_l_intersection)(NUM *out, PLANE *plane, NUM *p0, NUM *p1);
  int (*xy_l_intersection)(NUM *out, NUM *p0, NUM *p1);
  int (*p_inside_tri)(NUM *tri0, NUM *tri1, NUM *tri2, NUM *p);
  int (*trivloop)(FLUXON *f);

  void (*projmatrix)(NUM *out, NUM *x0_3, NUM *x1_3);
  NUM (*cart2_2d)(NUM *x1, NUM *x2);
  NUM (*cart2_3d)(NUM *x1, NUM *x2);
  NUM (*cart_2d)(NUM *x1, NUM *x2);
  NUM (*cart_3d)(NUM *x1, NUM *x2);
  NUM (*p_l_dist)(NUM *p0, NUM *x0, NUM *x1); 
  NUM (*p_ls_dist)(NUM *p0, NUM *x0, NUM *x1);
  NUM (*l_l_dist)( NUM a0[3], NUM b0[3], NUM c0[3], NUM d0[3]);
  void (*p_ls_closest_approach)(NUM p0[3], NUM a0[3], NUM b0[3], NUM c0[3]);
  void (*ls_closest_approach)(NUM p0[3], NUM p1[3], NUM a0[3], NUM b0[3], NUM c0[3], NUM d0[3]);
  NUM (*ls_ls_dist)(NUM a0[3], NUM b0[3], NUM c0[3], NUM d0[3]);
  NUM (*fl_segment_masked_dist)(VERTEX *v0, VERTEX *v1); 
  NUM (*fl_segment_masked_deluxe_dist)(NUM P0[3], NUM P1[3], VERTEX *v0, VERTEX *v1);
  NUM (*fl_segment_dist)(VERTEX *v1, VERTEX *v2);        
  NUM (*fl_segment_deluxe_dist)(NUM P0[3], NUM P1[3], VERTEX *v0, VERTEX *v1);

  int (*perp_bisector_2d)(NUM *out, NUM *P, NUM *Q); 
  int (*intersection_2d)(NUM *out, NUM *L1, NUM *L2);
  void (*project_n_fill)(VERTEX *v, DUMBLIST *horde);
  void (*hull_2d)(HULL_VERTEX *out, DUMBLIST *horde, DUMBLIST *rejects);
  void (*hull_2d_us)(HULL_VERTEX *hull, DUMBLIST *horde, VERTEX *central_v);

  int (*in_simplex)(POINT3D P0, POINT3D P1, POINT3D P2, POINT3D P3, POINT3D X);
  int (*above_plane)(POINT3D A, POINT3D B, POINT3D C, POINT3D X);
  DUMBLIST *(*find_simplex_by_location)(POINT3D x, WORLD *w, VERTEX *v, int global);
  DUMBLIST *(*find_nsimplex_by_location)(POINT3D x, WORLD *w, VERTEX *v, int global);
  VERTEX *(*find_vertex_by_location)(POINT3D x, WORLD *w, VERTEX *v, int global);

  NUM (*interpolate_lin_3d)( POINT3D x, NUM p[12], NUM val[4], int n, int tint);
  NUM (*interpolate_value_simplex)( POINT3D x, DUMBLIST *dl, int val_offset, int tint);
  NUM (*interpolate_value)( POINT3D x, WORLD *w, VERTEX *v, int global, int val_offset, int tint);
  void (*project_n_fill_photosphere)(VERTEX *v, DUMBLIST *horde);
  void (*hull_2d_us_photosphere)(HULL_VERTEX *hull, DUMBLIST *horde, VERTEX *central_v);


  /************  globals from data.h */
  char *code_info_model;
  char *code_info_physics;
  char *code_info_io;
  char *code_info_geometry;
  char *code_info_data;

  /*********** globals from model.c */
  struct F_B_NAMES *F_B_NAMES;

  /***********  globals from physics.h */
  struct FLUX_RECON *FLUX_RECON;
  struct FLUX_FORCES *FLUX_FORCES;


  /***********  functions in Core.xs */
  void *(*SvFluxPtr) (SV *sv, char *name, char *tstr, char wlflag, char croak_on_null);
  long (*SvLabel)    (SV *sv, char *name, char *tstr);
  void (*SvChangeLabel) (SV *sv, long newlabel, char *tstr);
  SV *(*new_sv_from_ptr) (WORLD *wptr, int type, long label);
  void (*destroy_sv)  (SV *sv);
} FluxCore;

  
/******************************
 * SvVertex
 * SvConc
 * SvFluxon
 * SvWorld
 *
 * Each macro calls SvFluxPtr to pull out the underlying pointer, and
 * casts it correctly.  Requires that the FLUX core has been
 * initialized in the current source file.
 */

#define SvVertex(sv,name,croak) ( (VERTEX *)            (FLUX->SvFluxPtr((sv),(name),"Flux::Vertex",0,(croak))))
#define SvConc(sv,name,croak)   ( (FLUX_CONCENTRATION *)(FLUX->SvFluxPtr((sv),(name),"Flux::Concentration",0,(croak))))
#define SvFluxon(sv,name,croak) ( (FLUXON *)            (FLUX->SvFluxPtr((sv),(name),"Flux::Fluxon",0,(croak))))
#define SvWorld(sv,name,croak)  ( (WORLD *)             (FLUX->SvFluxPtr((sv),(name),"Flux",1,(croak))))

/******************************
 * Type definitions...
 * These should be synchronized with the $typecodes hash in Flux.pm.
 * Also check the class names in Flux.xs.
 */
#define FT_LINKS  1
#define FT_VERTEX 2
#define FT_FLUXON 3
#define FT_WORLD  4
#define FT_CONC   5

#define MIN_FT 2
#define MAX_FT 5

static char *classnames[] = {
  "",
  "",
  "Flux::Vertex",
  "Flux::Fluxon",
  "Flux::World",
  "Flux::Concentration"
};

			 
