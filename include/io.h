/**********************************************************************
 * io.h -- I/O routine headers for FLUX
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Southwest Research Institute, 2004-2007
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
 * This file is part of FLUX 2.0 (31-Oct-2007).
 */

#ifndef FLEM_IO
#define FLEM_IO 1

#include "data.h"
#include <stdio.h>

#include <unistd.h>
#include <string.h> 
#include <sys/types.h>  
#include <sys/stat.h>
#include <sys/wait.h>  
#include <signal.h>
#include <stdio.h>

typedef struct FOOTPOINT_SPEC {
  long label;
  NUM x[3];
  NUM flux;
} FOOTPOINT_SPEC;

char *next_line(FILE *file);  /* Read line and skip comments */

/* Parse and act on a single non-comment line from a footpoint file */
int footpoint_action(WORLD *a, char *s);

/* Output routines */

/* Generic tree output to any file */
void fprint_tree(FILE *f, void *t, int lb, int lk, int idt, void ((*prntr)()));
void fprint_node(FILE *f, void *foo, int indent, int lab_o, int lk_o);

/* Generic tree output to stdout (hooks to above) */
void print_tree(void *t, int lab_o, int lk_o, int indent, void ((*printer)()));
void print_node(void *foo, int indent, int lab_o, int lk_o);

/* Fieldline output functions for print_tree */
void fdump_fluxon(FILE *f, FLUXON *foo, int indent);
void fprint_all_fluxon_node(FILE *f, FLUXON *foo, int indent);
void fprint_all_fluxon_tree(FILE *f, FLUXON *foo);
void print_all_fluxon_tree(FLUXON *foo);

/* Line output functions for state files */
void fprint_v_nbors(FILE *f, void *foo, int indnt, int lab_o, int ln_o);
void fprint_fc_line(FILE *f, void *foo, int indnt, int lab_o, int ln_o);
void fprint_fc_line_nonneg(FILE *f, void *foo, int indnt, int lab_o, int ln_o);
void fprint_fls_by_fc(FILE *f, void *foo, int indnt, int lab_o, int ln_o);
void fprint_fl_vertices(FILE *f, void *foo, int indnt, int lab_o, int ln_o);

/* Generic dumblist output */
void print_dumblist(DUMBLIST *foo, void ((*item_printer)()));

/* Output a whole world to a file */
int print_world(WORLD *a, char*header);
int fprint_world(FILE *file, WORLD *world, char *header);

/* Read a whole world from a file */
WORLD *read_world(FILE *file, WORLD *a);

/**********************************************************************
 **********************************************************************
 * 2-D Graphics output routines 
 */

int gl_2d_start(char *fname, float l, float r, float b, float t);
char *gl_2d_finish(float display_time, int display_mode);

void gl_2d_point(float x, float y, float radius, float colors[3]); 
void gl_2d_line (float x0, float y0, float x1, float y1, float colors[3]); 

void gl_2d_scr_poly(DUMBLIST *horde, float colors[3]);
void gl_2d_scr_list(DUMBLIST *horde,float colors[3]);


/**********************************************************************
 **********************************************************************
 * binary dump/restore stuff 
 */

// Some fences - unusual numbers with funny names 
#define BD_FENCE         0xAB5C155A
#define NEIGHBOR_FENCE   0xBED51DE0
#define NEARBY_FENCE     0xBED51DE1
#define VERTEX_END_FENCE 0x5AFEB055

// Assign type codes sequentially, and set BD_MAX_TYPENO to the value of the highest one.
#define BD_END                       1
#define BD_FLUXON_PIPE               2
#define BD_MAX_TYPENO                2

// The read/write routines.  Note that some (the _pipe routines) are intended ONLY 
// for returning state to a parent via a pipe -- they pass pointers directly!
int binary_dump_field(int fd, long code, long len, char *buf);
int binary_dump_end(int fd);
int binary_dump_fluxon_pipe(int fd, FLUXON *f);

int binary_read_dumpfile(int fd, WORLD *w);
int binary_read_fluxon_pipe(long size, char *buf, WORLD *w);
int binary_read_fluxon(int fd, FLUXON *f);


// This structure, containing the fixed elements of a VERTEX,
// is used for the fluxon_pipe routines (for parallelization)

typedef struct VERTEX_PIPE_FIXED {
  long label;
  POINT3D x;
  POINT3D scr;
  POINT3D b_vec;
  POINT3D f_v;
  POINT3D f_s;
  POINT3D f_t;
  POINT3D plan_step;
  NUM r;
  NUM a;
  NUM b_mag;
  NUM energy;
  NUM f_s_tot;
  NUM f_v_tot;
  NUM r_v, r_s;
  NUM r_cl;
  long neighbors_n;
  long nearby_n;
  NUM f_n_tot;
} VERTEX_PIPE_FIXED;


#endif /* overall file include */


