/**********************************************************************
 * geomview.c -- GeomView interface routines for FLEM
 * 
 * For now, this just consists of modules for dumping 
 * portions of fluxons into OOGL format.  Eventually, it will
 * include stuff for interacting with GeomView rather than just dumping
 * stuff to it.
 * 
 * This code is mostly for interfacing with GeomView via OOGL.  
 * It is mostly deprecated in favor of the newer PDL interface.
 *
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
 */
#include "geomview.h"
#include "model.h"
#include "hvectext.h"
#include <math.h>

#define CUBE_VERTICES
/******************************
 * expand_vertex_rim:  helper routine that expands a vertex's rim into a 
 * polygon normal to the local direction of the line. 
 * It returns a (3xn_faces) array of NUMs in static workspace.
 * 
 * You also have to feed in a radius.
 * 
 * expand_vertex_rim currently caches its trigonometric results --
 * this might be a waste of time, given that trig is getting cheaper 
 * every month.
 *
 * There's some twisting code that turns out not to be necessary most
 * of the time.  Might be worth tweaking to work right if a significant
 * number of field lines get bungied up:  GL runs MUCH faster if no polygons
 * intersect.
 *
 */

NUM *expand_vertex_rim(VERTEX *v, NUM r, int n_faces) {
  static NUM *buf = 0;
  static int bufsiz = 0;

  static NUM *pointcache = 0;
  static int cached_n_faces = 0;

  if(!v) {
    fprintf(stderr,"Warning: expand_vertex_rim received a null vertex (geomview.c)..\n");
    return 0;
  }

  /* Make sure the output buffer is there */
  if(bufsiz < n_faces) {
    if(buf)
      free(buf);
    buf = (NUM *)malloc(sizeof(NUM) * n_faces * 3 * 2);

    if(pointcache) 
      free(pointcache);
    pointcache = (NUM *)malloc(sizeof(NUM) * n_faces * 2 * 2);
    cached_n_faces = 0;

    bufsiz = n_faces*2;
    if(!buf || !pointcache) {
      fprintf(stderr,"Malloc failed in expand_vertex_rim (geomview.c)!\n");
      exit(-99);
    }
  }

  /* Do the actual work */
  { 
    NUM m0[9];
    NUM m1[9];
    NUM *mp;
    NUM m[9];
    NUM p[3];
    NUM a[3],ax[3];
    static NUM prev_x[3] = {1,0,0};
    static NUM prev_y[3] = {0,1,0};
    NUM nx;
    NUM *pc;
    int i;

    /* Project the line into the Z axis (and its perp. plane into 2-space) */
    projmatrix(m,
	       v->prev ? v->prev->x : v->x,
	       v->next ? v->next->x : v->x
	       );

    /* Generate points */
    if(cached_n_faces != n_faces) {
      pc = pointcache;
      for(i=0;i<n_faces;i++) {
	*(pc++) = cos(i * 2 * M_PI / n_faces);
	*(pc++) = sin(i * 2 * M_PI / n_faces);
      }
    }

    /* Stick the points into the buffer */
    pc = pointcache;
    p[2] = 0;
    for(i=0;i<n_faces;i++) {
      p[0] = *(pc++) * r;
      p[1] = *(pc++) * r;

      vec_mmult_3d(a, m, p);
      sum_3d(buf+i*3,a,v->x);

#ifdef never
      if(i==0) {
	prev_x[0] = a[0];
	prev_x[1] = a[1];
	prev_x[2] = a[2];
      }
#endif
    }

  } /* end of convenience block */

  return buf;
} /* end of routine */



/**********************************************************************
 * OOGL_vertex
 * vertices are implemented as OFF icosahedra at the moment. 
 * The OFF structure is to allow coloration of the vertices.
 * There's no reason for the icosahedron, except that it's cool.
 * (It may be sort of expensive to use for simple coolness value...)
 */
#ifdef ICOS_VERTICES
static NUM cache[36] = {
  0.0, 0.0, 2.0,
  1.788854, 0.000000, 0.894427,
  0.552786, 1.701302, 0.894427,
  -1.447214, 1.051462, 0.894427,
  -1.447214, -1.051462, 0.894427,
  0.552786, -1.701302, 0.894427,
  1.447214, 1.051462, -0.894427,
  -0.552786, 1.701302, -0.894427,
  -1.788854, 0.000000, -0.894427,
  -0.552786, -1.701302, -0.894427,
  1.447214, -1.051462, -0.894427,
  0.0, 0.0, -2.0
};

static int vertices[60] = {
  2,0,1,     3,0,2,     4,0,3,
  5,0,4,     1,0,5,     2,1,6,
  7,2,6,     3,2,7,     8,3,7,
  4,3,8,     9,4,8,     5,4,9,
  10,5,9,    6,1,10,    1,5,10,
  6,11,7,    7,11,8,    8,11,9,
  9,11,10,   10,11,6
};

static int n_verts=12;
static int n_gon=3;
static int n_faces=20;
static NUM r_scale=2;
#endif

#ifdef TETRA_VERTICES
static NUM cache[12] = {
  0.7746, -0.4472,-0.4472
  ,-0.7746,-0.4472,-0.4472
  ,0      , 0.7746,-0.4472
  ,0      , 0     ,0.7746
};

static int vertices[12] = {
  0,1,2,  0,2,3,  0,1,3,  1,2,3
};

static int n_verts=4;
static int n_gon=3;
static int n_faces=4;
static NUM r_scale=0.5;
#endif

#ifdef CUBE_VERTICES
static NUM cache[24]= {
  -1,-1,-1,   -1,-1,1,    -1,1,1,   -1,1,-1,
  1, -1,-1,    1,-1,1,     1,1,1,    1,1,-1
};

static int vertices[24] = {
  0,1,2,3,   0,1,5,4,   0,3,7,4,
  6,7,4,5,   6,7,3,2,   6,5,1,2};

static int n_verts=8;
static int n_gon=4;
static int n_faces=6;
static NUM r_scale=1.0;
#endif

void OOGL_vertex(FILE *f, NUM *x, NUM radius, NUM *color, char redefine) {
  NUM col_def[4] = {1.0,0.5,0.5,1.0};
  int i,j;
  char *s;

  if(!color) 
    color = col_def;
  
  
  fprintf(f,"{ appearance { material { ambient %g %g %g diffuse %g %g %g } } INST transform { \n\t%g 0 0 0\n\t0 %g 0 0\n\t0 0 %g 0\n\t%g %g %g 1 }\n geom { "
	  ,color[0],color[1],color[2],color[0],color[1],color[2]
	  ,radius,radius,radius,x[0],x[1],x[2]);
  
  
  if(redefine) {
    
    fprintf(f,"define vertex = {OFF\n");
    
    fprintf(f, "%d\t%d\t%d \n"
	    ,n_verts,n_faces,0
	    ,x[0],x[1],x[2]
	    );
    
    
    /* Dump vertices */
    for(i=0;i<n_verts;i++) {
      fprintf(f,"%8g\t%8g\t%8g\n"
	      ,cache[i*3+0]/r_scale
	      ,cache[i*3+1]/r_scale
	      ,cache[i*3+2]/r_scale
	      );
    }
    
    s = "3\t%2d %2d %2d\t%g %g %g %g\n";
    /* Dump vertex indices */
    for(i=0;i<n_faces;i++) {
      fprintf(f,"%d ",n_gon);
      for(j=0;j<n_gon;j++) {
	fprintf(f,"%d ",vertices[i*n_gon+j]);
      }
      fprintf(f,"\n");
      //      fprintf(f,"%g %g %g %g\n",color[0],color[1],color[2],color[3]);
    }
    fprintf(f,"}}\n");
  }
  else
    fprintf(f," : vertex }");
  
  fprintf(f,"}\n");
}

/**********************************************************************
 * OOGL_label
 * 
 * OOGL_label generates a label object using the hvectext module.  You feed in 
 * a location, normal vector, text direction vector, height
 * (in normalized coordinates), and string to print.
 * The command passes through the shell (for now) so be careful about
 * printing special characters -- especially "'"!
 *
 * Alignment is via the hvectext alignment options -- you just pass in 
 * the string to hand to hvectext.  The options are 'n','s','e','w',
 * 'ne','nw','se','sw','c', corresponding to the eight compass points 
 * looking down onto the text on an imaginary map.  'c' is centered.
 * If you pass in 0 instead of a string, you get 'w' alignment.
 */
void OOGL_label(FILE *f, char *string, NUM x[3], NUM dir[3], NUM normal[3], NUM height, char *align, NUM *color) {
  static int lno = 0;
  NUM dir_hat[3]= {1,0,0};
  NUM normal_hat[3] = {0,0,1};
  NUM third[3];
  NUM third_hat[3] = {0,1,0};
  int i;
  NUM a;

  /* Figure direction and normal */
  /* Cross direction with normal vectors to get textual "up", then */
  /* cross direction and textual "up" to get an honest-to-god normal. */

  if(dir && (a = norm_3d(dir))) 
    scale_3d(dir_hat,dir,1.0/a);

  if(normal && (a = norm_3d(normal)))
    cp_3d(normal_hat,normal);  /* Not really a hat vector (yet) */

  cross(third,normal_hat,dir_hat);
  if(a = norm_3d(third))
    scale_3d(third_hat,third,1.0/a);

  cross(normal_hat,dir_hat,third_hat);/* No need to rescale (unit by const.) */
  

  /* If a color is supplied, give one - else leave it as the default color */
  if(color) 
    fprintf(f,"{ appearance material {edgecolor %g %g %g}\n ",color[0],color[1],color[2]);
  else 
    fprintf(f,"{");
    


  /* Generate the transformation matrix.  Since this is the computer- */
  /* graphics world, matrices are really their own inverses.  Argh.   */
  fprintf(f,"### text: '%s'\n",string);
  
  fprintf(f,"INST transform {\n\t%g %g %g 0\n\t%g %g %g 0\n\t%g %g %g 0\n\t%g %g %g 1}\n\tgeom "
	  ,dir_hat[0],dir_hat[1],dir_hat[2]
	  ,third_hat[0],third_hat[1],third_hat[2]
	  ,normal_hat[0],normal_hat[1],normal_hat[2]
	  ,x[0],x[1],x[2]
	  );

  hvectext(f, 0, 0, height ? height : 0.1, string);

  fprintf(f,"}\n\n");

}


/**********************************************************************
 **********************************************************************
 *** 
 *** Not-so-general-purpose stuff starts here...
 ***
 */
 
  
/******************************
 * OOGL_dump_fluxon
 * Dump a whole fluxon as a giant OFF.
 * (WARNING: This may be expanded soon to include color and sliding-scale
 * colors....)
 */
#define OOGL_N_FACES 4
#define OOGL_RADIUS 0.05
void OOGL_dump_fluxon(FILE *f, FLUXON *fl, NUM radius, int n_faces, char label, NUM *color, NUM *vertex_color) {
  NUM col_defaults[4]={0.5,1.0,0.5,1.0};
  NUM vcol_defaults[4]={1.0,0.5,0.5,1.0};
  int i,j,vct;
  VERTEX *v;


  if(!color)
    color=col_defaults;

  if(!vertex_color)
    vertex_color = vcol_defaults;
    
  if(!f || !fl) {
    fprintf(stderr,"OOGL_dump_fluxon: Null file or fluxon -- ignoring.\n");
    return ;
  }

  fprintf(f,"\n{LIST\n");

  if(fl->v_ct < 2) {
    fprintf(f,"\n# Skipped trivial fluxon with %d vertices...\n",fl->v_ct);
    return;
  }

  /* Start the output */
  fprintf(f,"{OFF\n");
  
  /* File size */
  fprintf(f,"%d\t%d\t%d # NVerts, NFaces, NEdges\n",
	  n_faces * fl->v_ct,
	  n_faces * (fl->v_ct-1) + 2,
	  n_faces * (fl->v_ct * 2 - 1)
	  );

  /* Dump vertices */
  i=0;
  for(v = fl->start; v && i<fl->v_ct; v=v->next) {
    NUM *rim = expand_vertex_rim(v,radius*2.0/3.0,n_faces);
    for(j=0;j<n_faces;j++) {
      fprintf(f,"%8g\t%8g\t%8g # vtx %d  no. %d\n"
	      ,finite(rim[0])?rim[0]:0.0
	      ,finite(rim[1])?rim[1]:0.0
	      ,finite(rim[2])?rim[2]:0.0
	      ,i,j
	      );
      rim += 3;
    }
    i++;
  }
  
  /* Complain if necessary -- v should be zero AND i should == fl->v_ct. */
  if(v || i < fl->v_ct) {
    if(v) {
      fprintf(f,"# WARNING: line truncated due to v_ct too high (%d; i=%d)\n",fl->v_ct,i);
      fprintf(stderr,"WARNING: v_ct too high in OOGL_dump_fluxon.\n");
    } else {
      fprintf(f,"# WARNING: incorrect COFF output due to v_ct too low\n");
      fprintf(stderr,"WARNING: v_ct too low in OOGL_dump_fluxon.\n");
    }
  }

  /* Dump edge faces */
  for(i=0;i<fl->v_ct-1;i++) 
    for(j=0;j<n_faces;j++)
      fprintf(f,"4\t%d\t%d\t%d\t%d\t%g %g %g %g\t# vtx %d no. %d\n",
	      i     * n_faces + j,
	      i     * n_faces + (j + 1) % n_faces,
	      (i+1) * n_faces + (j + 1) % n_faces,
	      (i+1) * n_faces + j,
	      color[0],color[1],color[2],color[3],
	      i,j
	      );

  /* Dump end faces */
  fprintf(f,"%d\t",n_faces);
  for(i=0;i<n_faces;i++) 
    fprintf(f,"%d\t",i);
  fprintf(f,"%g %g %g %g\n",color[0],color[1],color[2],color[3]);

  /* reverse order to point the normal outward! */
  fprintf(f,"%d\t",n_faces);
  for(i=0;i<n_faces;i++)
    fprintf(f,"%d\t", n_faces * fl->v_ct - i - 1); 
  fprintf(f,"%g %g %g %g\n",color[0],color[1],color[2],color[3]);

  fprintf(f,"}\n"); /* Close out the OFF object */

  /* Dump spheres where all the vertices *really* lie */
  for(v=fl->start;v;v=v->next) {
    OOGL_vertex(f,v->x,radius*2,0,(v==fl->start));
  }

  /* Label all the segments, at both ends */
  if(label) {
    for(v = fl->start;v->next; v=v->next) {
      NUM x[3], seg[3],seg1[3], seg_hat[3], offset[3];
      NUM dir[3];
      NUM norm[3] = {0,0,1};
      NUM norm_hat[3];
      char buf[1024];
      
      
      diff_3d(seg,v->next->x,v->x);
      
      cross_3d(dir,seg,norm);
      sum_3d(dir,dir,norm);
      scale_3d(offset,dir,2 * radius / norm_3d(dir));
      
      cross_3d(norm_hat,dir,seg);
      scale_3d(norm_hat,norm_hat,1.0/norm_3d(norm_hat));
      scale_3d(seg_hat,seg,1.0/norm_3d(seg));
      
      sprintf(buf,"L%d: V%d",v->line->label,v->label);
      
      
      scale_3d(seg1,seg,0.1);
      sum_3d(x,v->x,seg1);
      sum_3d(x,x,offset);
      sum_3d(norm,norm_hat,seg_hat); /* tilt 45 degrees */
      OOGL_label(f,buf, x, dir, norm, radius*2, 0,color);
            
    }
  }
  fprintf(f,"}\n\n"); /* Close out the LIST object */
}


/**********************************************************************
 * OOGL_neighbor_lines
 * 
 * Slightly less general-purpose than the other OOGL routines, 
 * neighbor_lines takes a VERTEX and draws lines to all of that
 * VERTEX's neighbors (according to the VERTEX's neighbor structure).
 */
static int line_no = 0;
void OOGL_neighbor_lines(FILE *f, VERTEX *v, DUMBLIST *neighbors, NUM *c1, NUM *c2) {
  NUM default_c1[4] = {1.0,1.0,0.0,1.0}; /* Default near source is yellow */
  NUM default_c2[4] = {0.5,0.75,0.75,1.0}; /* Default near neighbor is blue*/
  int i;
  
  if(!c1)
    c1 = default_c1;
  if(!c2)
    c2 = default_c2;

  if(!v) return;
  if(!(neighbors->n)) return;

  fprintf(f,"{LIST\nappearance {+vect shading smooth} VECT # %d\n",line_no++);
  fprintf(f,"%d %d %d # polylines vertices colors\n",neighbors->n, 2*neighbors->n,2*neighbors->n);

  /* Two vertices per polyline; n polylines */
  for(i=0;i<neighbors->n;i++)
    fprintf(f,"2 ");
  fprintf(f,"\n");

  /* Two colors per polyline; n polylines */
  for(i=0;i<neighbors->n;i++)
    fprintf(f,"2 ");
  fprintf(f,"\n");
  
  
  /* Dump vertex locations */
  for(i=0;i<neighbors->n;i++) {
    VERTEX *a = ((VERTEX *)(neighbors->stuff[i]));
    if(a->next && v->next) {
      fprintf(f,"%8g\t%8g\t%8g # %d\n"
	      ,v->x[0]*0.52+v->next->x[0]*0.48
	      ,v->x[1]*0.52+v->next->x[1]*0.48
	      ,v->x[2]*0.52+v->next->x[2]*0.48  /* Central point --offset to avoid superposition */
	      ,line_no++);
      fprintf(f,"%8g\t%8g\t%8g\n"
	      ,a->x[0]*0.48+a->next->x[0]*0.52
	      ,a->x[1]*0.48+a->next->x[1]*0.52
	      ,a->x[2]*0.48+a->next->x[2]*0.52);  /* Neighbor */
    } else {
      fprintf(f,"%8g\t%8g\t%8g\n"
	      ,v->x[0],v->x[1],v->x[2]);
      fprintf(f,"%8g\t%8g\t%8g\n"
	      ,a->x[0],a->x[1],a->x[2]);
    }
  }

  /* Dump colors */
  for(i=0;i<neighbors->n;i++) 
    fprintf(f,"%8g %8g %8g %8g %8g %8g %8g %8g # %d\n"
	    ,c1[0],c1[1],c1[2],c1[3],c2[0],c2[1],c2[2],c2[3],line_no++);

  fprintf(f,"}\n");
}

/**********************************************************************
 * OOGL_neigborhood_cell
 * 
 * OOGL_neighborhood_cell takes a VERTEX and draws the Voronoi cell
 * around it at the midpoint.  The distances used are one quarter of
 * the actual distances.
 */
static inline void unproject(NUM *dest, NUM *source, NUM *matrix, NUM scale, NUM *offset) {
  NUM a[3];
  NUM b[3];
  a[2] = 0;
  a[1] = source[1];
  a[0] = source[0];

  vec_mmult_3d(b,matrix,a);
  scale_3d(b,b,scale);
  sum_3d(dest,offset,b);
}
  
void OOGL_neighborhood_cell(FILE *f, VERTEX *v, NUM *color,NUM *c2) {
  HULL_VERTEX *hull;
  NUM default_color[4]={1.0,0.2,1.0,1.0}; /* purple */
  NUM default_c2[4]={1.0,0,0.3,1.0}; /* red */
  NUM scale=0.1;
  NUM pm[9];
  NUM a[3],b[3];
  NUM centroid[3];
  NUM *foo;
  int i,j,k;
  int n_closed,n_cl_2;

  if(v->neighbors.n == 0)
    return;

  if(!color)
    color=default_color;

  if(!c2)
    c2 = default_c2;
    
  if(!v->next)
    return;

  /* Get the perpendicular-plane projection matrix */
  projmatrix(pm, v->x, v->next->x);

  /* Generate line segment centroid */
  sum_3d(centroid,v->x,v->next->x);
  scale_3d(centroid,centroid,1.0/2);

  /* Generate the hull */
  hull = hull_neighbors(v, &(v->neighbors));
  foo = (NUM *)malloc(sizeof(NUM)*3*v->neighbors.n);
  
  fprintf(f,"{LIST #Voronoi-diagram # %d\n",line_no++);


  /******************************
   ** Draw the hull
   */
  
  /* Count how many actual lines we'll draw */
  /* n_closed gets number of line segments with two closed hull_vertices; */
  /* n_cl_2 gets number of closed hull_vertices */
  n_closed = v->neighbors.n;
  for(n_cl_2=j=i=0;i<v->neighbors.n; i++) {
      if(hull[i].open) {
	n_closed--;
	if(!hull[(i+1)%v->neighbors.n].open)
	  n_closed--;
      }
      else {
	n_cl_2++;
      }
  }
  if(n_closed > 0) {

    fprintf(f,"{VECT #hull\n");
    fprintf(f,"%d %d %d # polylines vertices colors (hull)\n",n_closed ,n_closed*2,n_closed);
    
    /* 2 vertices per line */
    for(i=0;i<n_closed;i++) {
      fprintf(f,"%d ",2);
    }
    fprintf(f,"\n");
    for(i=0;i<n_closed;i++) {
      fprintf(f,"%d ",1);
    }
    fprintf(f,"\n");
    
    for(j=k=i=0;i<v->neighbors.n; i++) {
      if(!hull[i].open) {
	if(i==0 || hull[i-1].open) {
	  unproject(a,hull[i].p,pm,scale,centroid);
	}
	cp_3d(foo+3*(k++),a);
	
	if(!hull[(i+1)%v->neighbors.n].open) {
	  fprintf(f,"%8g\t%8g\t%8g\n",a[0],a[1],a[2]);
	  
	  unproject(a,hull[(i+1)%v->neighbors.n].p,pm,scale,centroid);
	  fprintf(f,"%8g\t%8g\t%8g\n",a[0],a[1],a[2]);
	  j++;
	}
      }
    }
    
    if(j!=n_closed) {
      fprintf(stderr,"ASSERTION FAILED:  j(%d) != n_closed(%d)\n",j,n_closed);
      exit(125);
    }
    if(k!=n_cl_2) {
      fprintf(stderr,"ASSERTION FAILED:  k(%d) != n_cl_2(%d)\n",k,n_cl_2);
      exit(124);
    }
    
    for(i=0;i<n_closed;i++) {
      NUM alpha = (i / (k+(k==0)));
      
      fprintf(f,"%8g %8g %8g %8g # color(line %d)\n"
	      ,color[0]*(1 - alpha)
	      ,color[1]*(1 - alpha)
	      ,color[2]*(1 - alpha)
	      ,color[3]
	      ,line_no++
	      );
    }
    
    fprintf(f,"}\n");


    /******************************
     ** Draw lines from each closed vertex to its associated neighbors.
     */
    if(n_cl_2 > 0) {
      fprintf(f,"{VECT \n");
      fprintf(f,"%d %d %d \n",n_cl_2,n_cl_2*3,n_cl_2);
      for(i=0;i<n_cl_2;i++)
	fprintf(f,"%d ",3);
      fprintf(f,"\n");
      for(i=0;i<n_cl_2; i++)
	fprintf(f,"1 ");
      fprintf(f,"\n");
      
      for(k=i=0;i<v->neighbors.n; i++) {
	if(!hull[i].open) {
	  unproject(a,((VERTEX **)(v->neighbors.stuff))[i]->scr,pm,scale,centroid);
	  fprintf(f,"%8g\t%8g\t%8g\n",a[0],a[1],a[2]);
	  
	  unproject(a,hull[i].p,pm,scale,centroid);
	  fprintf(f,"%8g\t%8g\t%8g\n",a[0],a[1],a[2]);
	  
	  unproject(a,((VERTEX **)(v->neighbors.stuff))[(i+1)%v->neighbors.n]->scr,pm,scale,centroid);
	  fprintf(f,"%8g\t%8g\t%8g\n",a[0],a[1],a[2]);
	  k++;
	}
      }
      if(k!=n_cl_2) {
	fprintf(stderr,"ASSERTION FAILED:  k(%d) != n_cl_2(%d)\n",k,n_cl_2);
	exit(126);
      }
      for(i=0;i<n_cl_2;i++)
	fprintf(f,"0 0 0 1\n");
      fprintf(f,"}\n");
    }
    
    /**************
     ** Draw the hull's vertices, and label them
     **/
    for(i=0;i<v->neighbors.n;i++) {
      char buf[80];
      if(!hull[i].open){
	OOGL_vertex(f,foo+3*i,0.001,color,(i==0));
	sprintf(buf," %d - %d"
		,((VERTEX *)(v->neighbors.stuff[i]))->label
		,((VERTEX *)(v->neighbors.stuff[(i+1)%v->neighbors.n]))->label);
	OOGL_label(f,buf,foo+3*i,0,0,0.003,0,color);
      }
    }
  }    
 /*******************
  ** Draw the Delaunay lines 
  **/
 fprintf(f,"{VECT # neighbors\n");
 fprintf(f,"%d %d %d #polylines vertices colors (neighbors)\n",v->neighbors.n, v->neighbors.n*2,v->neighbors.n);
 for(i=0;i<v->neighbors.n; i++) 
   fprintf(f,"%d ",2);
   fprintf(f,"\n");
 for(i=0;i<v->neighbors.n; i++)
   fprintf(f,"1 ");
 fprintf(f,"\n");


 for(i=0;i<v->neighbors.n;i++) {
   fprintf(f,"%8g %8g %8g\n",centroid[0],centroid[1],centroid[2]);

   cp_3d(b,((VERTEX *)(v->neighbors.stuff[i]))->scr);
   b[2] = 0;
   scale_3d(b,b,scale);

   /* Now transform it back out of the plane of the hull diagram 
      and into world coordinates  */
   vec_mmult_3d(a,pm,b);
   sum_3d(a,a,centroid);

   /* Draw the line and save it for later vertexification */
   fprintf(f,"%8g %8g %8g\n",a[0],a[1],a[2]);
   cp_3d(foo+3*i,a);
 }
 for(i=0;i<v->neighbors.n;i++) {
   NUM alpha = i/(float)((v->neighbors.n - 1) + (v->neighbors.n == 1));
   
   fprintf(f,"%8g %8g %8g %8g # color\n"
	   ,c2[0]*(1-alpha)
	   ,c2[1]*(1-alpha)
	   ,c2[2]*(1-alpha)
	   ,c2[3]
	   );
 }

 fprintf(f,"}\n");


 /******************************
  ** Draw the neighbors' projected locations, and label them.
  **/
 for(i=0;i<v->neighbors.n;i++) {
   char buf[80];
   NUM direction[3] = {0,1,1};

   OOGL_vertex(f,foo+3*i,0.0015,c2,(i==0));

   sprintf(buf," %d",((VERTEX *)(v->neighbors.stuff[i]))->label);
   OOGL_label(f,buf,foo+3*i,direction,0,0.005,0,c2);
 }
 
 fprintf(f,"}\n");

}      


/**********************************************************************
 * OOGL_forces
 * 
 * OOGL_forces takes a VERTEX and draws a vector representing the
 * total magnetic force on it.  
 */
static void scale_force(NUM out[3], NUM in[3]) {
  NUM n1;
  n1 = norm_3d(in);
  scale_3d(out,in,log10(n1+1) * 0.1 / (n1+1));
}


void OOGL_draw_force(FILE *f, VERTEX *v, NUM *color) {
  NUM default_color[4] = {1.0,1.0,0.0,1.0}; /* yellow */
  NUM a[3];
  NUM f_[3];
  int i;

  if(!color)
    color = default_color;

  if(!v->next)
    return;

  fprintf(f,"{VECT #force on vertex %d\n",v->label);
  fprintf(f,"3 6 3\n");
  fprintf(f,"2 2 2\n1 1 1\n"); /* 2 vertices per line, one color per line */

  /* Total force on vertex */
  cp_3d(a,v->x);
  fprintf(f,"%8g\t%8g\t%8g\n",a[0],a[1],a[2]);
  scale_force(f_,v->f_t);
  sum_3d(a,a,f_);
  fprintf(f,"%8g\t%8g\t%8g\n",a[0],a[1],a[2]);
  
  /* "Vertex" force */
  cp_3d(a,v->x);
  fprintf(f,"%8g\t%8g\t%8g\n",a[0],a[1],a[2]);
  scale_force(f_,v->f_v);
  sum_3d(a,a,f_);
  fprintf(f,"%8g\t%8g\t%8g\n",a[0],a[1],a[2]);
  
  /* "Segment" force in the middle of the segment */
  sum_3d(a,v->x,v->next->x);
  scale_3d(a,a,0.5);
  fprintf(f,"%8g\t%8g\t%8g\n",a[0],a[1],a[2]);
  scale_force(f_,v->f_s);
  sum_3d(a,a,f_);
  fprintf(f,"%8g\t%8g\t%8g\n",a[0],a[1],a[2]);
  
  fprintf(f,"1 1 1 1\n"); /* white for total force */
  fprintf(f,"%8g %8g %8g %8g \n",color[0],color[1],color[2],color[3]);
  fprintf(f,"%8g %8g %8g %8g \n",color[0],color[1],color[2],color[3]);
  fprintf(f,"}\n");
}

/**********************************************************************
 * OOGL_vect
 * 
 * Draw an arbitrary vect
 * 
 **********************************************************************/
void OOGL_vect(FILE *f, int vct, NUM **points, NUM *color) {
  NUM def_clr[] = {0,0,0,0};
  NUM *clr = color;
  int i;

  fprintf(f,"{VECT\n");

  fprintf(f,"1 %d 1\n %d\n1\n",vct,vct);
  
  for(i=0;i<vct;i++,points++)
    fprintf(f,"%g %g %g\n",(*points)[0],(*points)[1],(*points)[2]);
  
  if(!clr)
    clr = def_clr;

  fprintf(f,"%g %g %g %g\n}\n",clr[0],clr[1],clr[2],clr[3]);
  
}

void OOGL_segment(FILE *f, NUM *a, NUM *b, NUM *color) {
  NUM *(p[2]);
  p[0] = a;
  p[1] = b;
  OOGL_vect(f,2,p,color);
}

