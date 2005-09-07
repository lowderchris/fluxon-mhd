/**********************************************************************
 * io.c -- I/O routines for FLEM
 *
 * Mostly text I/O; now has graphics output for extra cheesy flavor!
 * OpenGL output is implemented, well, oddly.  The model is one of
 * building a display list and then, well, displaying it.  You send
 * stuff to an output file and then ask for it to be displayed.  The 
 * current implementation just writes a perl script that uses the OpenGL
 * glue to display the list you want.  Then when you close/display the
 * visual, it closes the file and executes the script. 
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
 * This is io.c version 1.0 - part of the FLUX 1.0 release.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <regex.h>
#include <math.h>

#include "data.h"
#include "io.h"
#include "physics.h" /* for FORCES translation */

/**********************************************************************
 * 
 * next_line -- read a line from an open footpoint action file
 * (open it with fopen()).  Skip over comments.  Return either 
 * the next line, or 0 on EOF or other error.
 * Obliterate whatever came last (there's only a one line buffer
 * here).
 */

static char buf[BUFSIZ];

char *next_line(FILE *file) {
  char *s;

  if(!file) return NULL;
  
  for(buf[0]='#',buf[1]='\0';buf[0]=='#';){
    char *s;
    if(!(fgets(buf,BUFSIZ,file))) return NULL;

    for(s=buf;*s && isspace(*s);s++) /* Check for all-whitespace */
      ;
    if(!*s)                          /* All-whitespace: */
      buf[0] = '#',buf[1]='\0';      /* treat this buffer as a comment */
  }
  

  /* Get rid of the newline character */
  for(s=buf; *s && *s != '\n'; s++)
    ;
  *s='\0';
    
  return buf;
}


/**********************************************************************
 *
 * footpoint_action -- take action on a WORLD based on the string 
 * contained in the parameter (which is presumably a line from next_line).
 *
 * Return zero if normal completion, -1 if the line contained a 
 * terminator.  Return 1 if there was an unrecognized directive, 2
 * if there was a recognized directive that failed.
 *
 * Initial whitespace is ignored.
 * Comments are silently ignored. 
 * (Comment lines have '#' as the first non-whitespace character)
 *
 */
int footpoint_action(WORLD *world, char *s) {
  char *ds = 0;
  char *badstr = 0;
  char c;
  int i,n;
  long vertex_label;
  long l0,l1,fl0;
  NUM x0[3],x1[3],flux0,flux1;
  static char mscan[80]="";
  static char lscan[80]="";
  static char vscan[80]="";
  static char gscan[80]="";

  /* Skip over initial whitespace and comment characters */
  while(*s && isspace(*s))
    *s++;

  /* Return on blank line or comment */
  if(!*s || *s=='#') return 0;

  c = toupper(*s);

  if(world==NULL) {
    fprintf(stderr,"footpoint_action got a null world!\n");
    exit(2);
  }

  switch(*s) {
  case 'F': /* FRAME <n>
	     * FINAL 
	     */
    if(toupper(s[1])=='R') {
      n=sscanf(s,"%*s %ld",&i);
      if(n != 1) 
	badstr = "Bad parse in FRAME line";
      else {
	if(world->state == WORLD_STATE_NEW) {
	  world->frame_number = i;
	  world->state = WORLD_STATE_LOADING;
	} else if(world->state == WORLD_STATE_READY) {
	  if(world->frame_number < i) {
	    world->state = WORLD_STATE_LOADING;
	    world->frame_number = i;
	  } else {
	    badstr = "New frame number is less than old one";
	  }
	} else {
	  badstr = "World state is not good";
	}
      }
    } 
    else if(toupper(s[1]) == 'I') {
      if(world->state == WORLD_STATE_LOADING) {
	world->state = WORLD_STATE_LOADED;
	return -1;
      } else {
	badstr = "World is not ready for this state";
      }
    } else {
      badstr = "Unknown 'F' directive";
    }
    break;

  case 'M': /* MOVE <label> <x> <y> <z> <flux>
	     */
  case 'N': /* NEW  <label> <x> <y> <z> <flux>
	     */

    /* Initialize the scan string, if necessary */
    if(!*mscan) {
      sprintf(mscan,"%%*s %%ld %%%sf %%%sf %%%sf %%%sf",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
    }

    n = sscanf(s,mscan, &l0,x0,x0+1,x0+2,&flux0);

    if(n != 5) {
      badstr = "Couldn't parse NEW line";
    } else {
      if(toupper(*s)=='N') {

	/* New stuff handler */

	FLUX_CONCENTRATION *fc;
	
	fc = new_flux_concentration( world,x0[0],x0[1],x0[2],flux0,l0 );
	world->concentrations = 
	  tree_binsert( world->concentrations
			,  fc 
			, fc_lab_of, fc_ln_of );
      } else {
	
	/* Motion handler */
	FLUX_CONCENTRATION *fc;
	fc = tree_find(world->concentrations
		       , l0
		       , fc_lab_of, fc_ln_of);
	if(fc==NULL) {
	  badstr = "Couldn't find labelled flux concentration in the tree";
	} else {
	  /* Copy the new flux and location into the concentration. */
	  NUM *n0 = fc->x;
	  NUM *n1 = x0;
	  fc->flux = flux0;
	  *n0++ = *n1++;
	  *n0++ = *n1++;
	  *n0 = *n1;
	}
      }
    }
    break;
    
  case 'D': /* DELETE <label>
	     */
    badstr = "Delete not implemented";
    break;
    
  case 'S': /* SPLIT <l1> <l2> <x1> <y1> <z1> <x2> <y2> <z2> <fl1> <fl2> 
	     * (Split <l1> into two parts with two locations and fluxes 
	     */
    badstr = "Split not implemented";
    break;

  case 'J': /* JOIN <l1> <l2> <x> <y> <z> <flux>
	     *  (Join <l1> and <l2> into <l1> at <x>, <y>, <z>)
	     */
    badstr = "Join not implemented";
    break;

  case 'E': /* EMERGE */
    badstr = "Emerge not implemented";
    break;

  case 'C': /* CANCEL */
    badstr = "Cancel not implemented";
    break;

  /******************************
   * End of footpoint motion codes; start of fluxon and vertex 
   * specification codes
   */
  case 'L': /* LINE <fl_lab> <lab1> <lab2> <flux> 
	     *  Create a field line from 1 to 2 with label and flux given 
	     */
    
    if(!*lscan) {
      sprintf(lscan,"%%*s %%ld %%ld %%ld %%%sf",NUMCHAR,NUMCHAR,NUMCHAR);
    }

    n = sscanf(s,lscan, &fl0, &l0, &l1, &flux0);
    
    if(n != 4) {
      badstr = "Couldn't parse LINE line";
    } else {
      FLUX_CONCENTRATION *fc0, *fc1;
      fc0 = tree_find(world->concentrations, l0, fc_lab_of, fc_ln_of);
      fc1 = tree_find(world->concentrations, l1, fc_lab_of, fc_ln_of);
      if(!fc0 || !fc1) {
	char *badbuf = (char *)localmalloc(BUFSIZ,MALLOC_MISC);
	sprintf(badbuf,"Found a fluxon specifier between concentrations %ld and %ld, but they \ncame up %ld and %ld in tree_find (one doesn't exist)!  \nThis error message leaked %d bytes (don't let it happen again!)\n",l0,l1,fc0,fc1,BUFSIZ);
      } else if(fc0->flux * fc1->flux >= 0) {
	badstr = (char *)localmalloc(BUFSIZ,MALLOC_MISC);
	sprintf(badstr,"This fluxon connects two flux concentrations of the same sign flux, or\none of its flux tubes has zero flux. Line %ld; concentrations %ld (%g) and %ld (%g)\n",fl0, flux0, l0, fc0->flux, l1, fc1->flux);
      }	else {
	/* Check if the field line exists in either of the two 
	   concentrations' local lists */
	FLUXON *f;
	
	if((f = tree_find(world->lines,fl0,fl_lab_of,fl_all_ln_of)) &&
	   (f->fc0 != fc0 || f->fc1 != fc1)
	   ) {
	  badstr = (char *)localmalloc(BUFSIZ,MALLOC_MISC);
	  sprintf(badstr,"Hey!  fluxon %ld already exists (from conc. %ld - %ld)\n"); 
	} else {
	  FLUXON *f0;

	  if(f) {
	    f->flux = fabs(flux0);
	  } else {
	    /* Regularize fluxon order */
	    if(fc0->flux < 0) {  /* Yes, I know about the triple-^= trick. */
	      FLUX_CONCENTRATION *a;    
	      a=fc1;
	      fc1=fc0;
	      fc0=a;
	    }
	    
	    /* Manufacture the new fluxon */
	    f0 = new_fluxon(fabs(flux0),fc0,fc1,fl0,0);
	    
	    /* Link it in (three ways!) */
	    /* This works because each flux concentration can only be the */
	    /* beginning of all its lines, or the end of all its lines! */
	    fc0->lines = tree_binsert(fc0->lines, f0, fl_lab_of, fl_start_ln_of);
	    fc1->lines = tree_binsert(fc1->lines, f0, fl_lab_of, fl_end_ln_of);
	    world->lines   = tree_binsert(world->lines,   f0, fl_lab_of, fl_all_ln_of);
	    
	    /* Make sure it has at least the two end vertices */
	    {
	      VERTEX *v0, *v1;
	      v0 = new_vertex(0,fc0->x[0],fc0->x[1],fc0->x[2],f0);
	      v1 = new_vertex(0,fc1->x[0],fc1->x[1],fc1->x[2],f0);
	      if(!v0 || !v1) {
		badstr = "Couldn't make trivial vertices for fluxon!\n";
	      } else {
		if( add_vertex_pos(f0,0,v0) || add_vertex_pos(f0,-1,v1) ) 
		  badstr = "Problems with trivial vertex addition\n";
	      }
	    }
	  }
	}
      }
    }

    break;
    
  case 'V': /* VERTEX <fl_lab> <vert_lab> <pos> <x> <y> <z> 
             *  Create a new vertex on field line <fl_lab> at position
	     *  <pos> (0 means start; -1 means end) along the existing one.
	     */
    if(!*vscan) {
      sprintf(vscan,"%%*s %%ld %%ld %%ld %%%sf %%%sf %%%sf",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
    }
    
    n = sscanf(s,vscan,&l0,&vertex_label,&l1,x0,x0+1,x0+2);
    if(n != 6) {
      badstr = "Couldn't parse VERTEX line";
    } else {
      VERTEX *v0;
      FLUXON *f;
      f = tree_find(world->lines, l0, fl_lab_of, fl_all_ln_of);
      if(!f) {
	badstr = "Tried to add a vertex to a nonexistent field line";
      } else {
	v0 = new_vertex(vertex_label, x0[0], x0[1], x0[2], f);
	if(world->verbosity > 2) printf("vertex_label = %d; v0->label = %d\n",vertex_label, v0->label);
	if(add_vertex_pos(f, l1, v0))
	  badstr = "VERTEX: add_vertex_pos returned error";
      }
    }      
    break;
  
  case 'G': /* GLOBAL <cmd> <arg1> <arg2> ... */
    {
      char gcmd[81];
      int off,n;
      n = sscanf(s,"%*s %80s %n",gcmd,&off);

      if(n != 1) {
	static char bstr[80];
	sprintf(bstr,"Couldn't parse GLOBAL line (n=%d)\n",n);
	badstr = bstr;
      } else {
	switch(*gcmd) {
	case 'F': /* GLOBAL FORCE ... */
	  {
	    int i,j,o2;
	    char fname[81];
	    void ((*(f_funcs[N_FORCE_FUNCS]))());
	    
	    for(i=0; 
		s[off] && i<N_FORCE_FUNCS && sscanf(s+off,"%80s %n",fname,&o2);
		i++) {

	      for(j=0; 
		  FLUX_FORCES[j].func && strcmp(FLUX_FORCES[j].name,fname);
		  j++)
		;
	      if(FLUX_FORCES[j].func) {
		f_funcs[i] = FLUX_FORCES[j].func;
	      } else {
		long foo;
		if( sscanf(fname,"0x%x",&foo) )
		  (void *)(f_funcs[i]) = (void *)foo;
		else {
		  static char bstr[160];
		  sprintf(bstr,"Couldn't parse GLOBAL FORCE LINE (bad force at '%s')",s+off);
		  badstr = bstr;
		  goto global_escape;
		}
	      } /* end of scanning branch */

	      off += o2;
	    } /* end of force specifier loop */

	    for(;i<N_FORCE_FUNCS;i++) /* pad with zeroes */
	      (void *)(f_funcs[i]) = (void *)0;

	    /* We didn't barf by now -- all forces must be OK.  Copy 'em in. */
	    for(j=0;j<N_FORCE_FUNCS;j++) 
	      world->f_funcs[j] = f_funcs[j];
	  } /* end of GLOBAL FORCE branch */
	  break;

	case 'P': /* GLOBAL PHOTOSPHERE ... */
	  {
	    PLANE *p = (PLANE *)0;
	    int typ;
	    int nval;

	    if(s[off]!='<') {  /* do nothing for <NONE> */
	      p = (PLANE *)localmalloc(sizeof(PLANE),MALLOC_PLANE);
	      char phscan[80];
	      sprintf(phscan,"%%%sf %%%sf %%%sf %%%sf %%%sf %%%sf %%%sf",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
	      if( 6 > (nval = sscanf(s+off,phscan
			      ,&(p->origin[0])
			      ,&(p->origin[1])
			      ,&(p->origin[2])
			      ,&(p->normal[0])
			      ,&(p->normal[1])
			      ,&(p->normal[2])
			      ,&typ
				     ) ) )
		badstr = "couldn't parse six values from GLOBAL PHOTOSPHERE ";
	      
	      if(nval==6)
		typ = PHOT_PLANE;
	      
	      if(!badstr) {
		world->photosphere.type = typ;
		if(world->photosphere.plane)
		  localfree(world->photosphere.plane);
		world->photosphere.plane = p;
	      }
	      world->photosphere.type = PHOT_PLANE;

	    } else { /* <NONE> */
	      world->photosphere.plane = 0;
	      world->photosphere.type = 0;
	    }
	  } /* end of GLOBAL PHOTOSPHERE convenience block */
          break;

	case 'B':  /* BFIELD_NORMALIZATION */
	  if(gcmd[1] == 'O' || gcmd[1]=='o') {
	    /* GLOBAL BOUNDARY ... */
	    int n;
	    int type_code;
	    PLANE *p = (PLANE *)localmalloc(sizeof(PLANE),MALLOC_PLANE);
	    char phscan[80];


	    switch(s[off]) {

	    case PHOT_PLANE+'0': 
	    case 'P':
	      /* Plane */
	      sprintf(phscan,"%%*s %%%sf %%%sf %%%sf %%%sf %%%sf %%%sf",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
	      if(6 > sscanf(s+off,phscan,
			    &(p->origin[0]),
			    &(p->origin[1]),
			    &(p->origin[2]),
			    &(p->normal[0]),
			    &(p->normal[1]),
			    &(p->normal[2])
			    )
		 )
		badstr = "Couldn't parse six values from GLOBAL BOUNDARY PLANE";
	      type_code = PHOT_PLANE;
	      
	      break;
	      
	    case PHOT_SPHERE+'0':
	    case 'S':
	      /* Sphere */
	      sprintf(phscan,"%%*s %%%sf %%%sf %%%sf %%%sf",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
	      if(4 > sscanf(s+off,phscan,
			    &(p->origin[0]),
			    &(p->origin[1]),
			    &(p->origin[2]),
			    &(p->normal[0])
			    )
		 )
		badstr = "Couldn't parse four values from GLOBAL BOUNDARY SPHERE";
	      p->normal[1] = 0;
	      p->normal[2] = 0;
	      type_code = PHOT_SPHERE;
	      break;
	      
	    case PHOT_CYL+'0':
	    case 'C':
	      /* Cylinder */
	      sprintf(phscan,"%%*s %%%sf",NUMCHAR);
	      if(1 > sscanf(s+off,phscan,
			    &(p->origin[0])
			    )
		 ) {
		badstr = "Couldn't parse cylindrical radius for GLOBAL BOUNDARY CYLINDER";
		printf("s+off: %s\nphscan: %s\n",s+off,phscan);
	      }
	      p->origin[1] = p->origin[2] = p->normal[0] = p->normal[1] = p->normal[2] = 0;
	      type_code = PHOT_CYL;
	      break;

	    case '0':
	    case 'N':
	      /* None */
	      free(p);
	      p=0;
	      type_code = 0;

	    default:

	      badstr = "Unknown boundary condition in GLOBAL BOUNDARY declaration\n";
	      break;
	    }

	    if(!badstr) { /* Parse was okay */
	      world->photosphere.type = type_code;
	      if(world->photosphere.plane) 
		localfree(world->photosphere.plane);
	      world->photosphere.plane = p;

	    } else { /* Clean up so we don't leak too much memory */
	      if(p)
		free(p);
	    }

            /* end of GLOBAL BOUNDARY case */

 	  } else {

	    /** GLOBAL B_NORMALIZATION ... **/
	    int i;
	    sscanf(s+off,"%d",&i) || (badstr="couldn't parse F_OVER_B state");
	    world->f_over_b_flag = i;
	  }
	  break;
	case 'S': /* GLOBAL STATE */
	  sscanf(s+off,"%d",&(world->state)) || 
	    (badstr="couldn't parse GLOBAL STATE");
	  break;

	default: /* GLOBAL <FOO> unrecognized */
	  fprintf(stderr,"Unrecognized GLOBAL command in line '%s'\n",s);
	  return 1;
	  global_escape: break;
	} /* end of gcmd switch */
      } /* end of ok-GLOBAL-line branch */
    } /* end of GLOBAL convenience block */
    break;

  default: /* Oops */
    fprintf(stderr,"Unrecognized command in line '%s'\n",s);
    return 1;
    break;
  }

  if(badstr) {
    fprintf(stderr, "%s; line: '%s'\n",badstr,s);
    return 2;
  } else {
    /*    printf("OK: '%s'\n",s);*/
    return 0;
  }
}
   
 
/**********************************************************************
 * Read_world
 *
 * Reads in a complete world from the specified FILE.   If you give
 * it an older WORLD, then the older WORLD is used as context.
 * Only one frame is read from the FILE, so that you can 
 * make multiple calls to read_world to update the state of a model.
 *
 * To determine whether the file is finished, check feof(file) after
 * exit.
 */

WORLD *read_world(FILE *file, WORLD *a) {
  char *s;
  int error = 0;
  
  do{
    if( a == NULL ) {
      a = new_world();
    }

    if(s = next_line(file)) {
      error = footpoint_action(a, s);
    }
  } while (s && !error);
  return a;
}    

/**********************************************************************
 * Print_world
 * Prints out a complete world in a format suitable for restoral by
 * read_world, to the specified FILE.  Returns an error status.
 * The third argument contains some stuff to stick into the top of the 
 * printout, as a comment line (or an infinite number!).  You're responsible
 * for putting '#' characters at the start of the header lines!
 */

int print_world(WORLD *a, char *header) {
  return fprint_world(stdout, a, header);
}

int fprint_world(FILE *file, WORLD *world, char *header) {
  char *s;
  int i;

  if(file==NULL) {
    fprintf(stderr,"print_world: null file supplied. Returning.\n");
    return 1;
  }

  if(world==NULL) {
    fprintf(stderr,"print_world:  got asked to print a null world.\n");
    return 1;
  }

  /* Output the header */
  fputs(header, file);
  
  /* Output the frame number that's in the world */
  fprintf(file, "\n\n######################################\n\n");
  fprintf(file, "FRAME %ld\n",world->frame_number);


  /***** Write out global configuration */
    /* force list */
  fprintf(file,"GLOBAL FORCES ");
  for(i=0;i<N_FORCE_FUNCS && world->f_funcs[i];i++) {
    int j;
    for(j=0; FLUX_FORCES[j].func && FLUX_FORCES[j].func != world->f_funcs[i]; j++) 
      ;
    if(FLUX_FORCES[j].func)
      fprintf(file,"%s ",FLUX_FORCES[j].name);
    else 
      fprintf(file,"0x%x ",world->f_funcs[i]);
  }
  fprintf(file,"\n");

  fprintf(file, "GLOBAL STATE %d  (%s)\n",world->state,world_state_name(world));

  /* global boundary type and location */
  {
    PLANE *p = world->photosphere.plane;
    char fmt[80];
    sprintf(fmt,"%%s %%s %%%sg %%%sg %%%sg %%%sg %%%sg %%%sg\n",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
    fprintf(file,fmt
	    , "GLOBAL BOUNDARY "
	    , BOUNDARY_NAMES[ world->photosphere.type ]
	    , p?p->origin[0]:0, p?p->origin[1]:0, p?p->origin[2]:0
	    , p?p->normal[0]:0, p?p->normal[1]:0, p?p->normal[2]:0
	    );
  }

  fprintf(file,"GLOBAL BFIELD_NORMALIZATION %d",world->f_over_b_flag);
  fprintf(file,"\n"); /* leave an extra space after the globals */
    
  /* Write out all flux concentration locations */
  fprint_tree(file, world->concentrations, fc_lab_of, fc_ln_of, 0, fprint_fc_line);

  /* Output the field lines */
  fputc('\n',file);
  fprint_tree(file, world->concentrations, fc_lab_of, fc_ln_of, 0, fprint_fls_by_fc);

}    
    

    
/**********************************************************************
 **********************************************************************
 ** 
 ** Tree printing routines -- useful for debugging, probably not
 ** used much in production code.  (Would go in i/o but really are
 ** meant for debugging trees and not as the main i/o method).
 **/

/**********************************************************************
 * print_tree: print out a whole tree
 */
void fprint_tree(FILE *file, void *tree, int label_offset, int link_offset, int indent, void ((*printer)())){ 
  LINKS *node_links = (LINKS *)(tree+link_offset);

  if(tree == NULL) {
    printf("<<NO TREE>>\n");fflush(stdout);
    return;
  }

  if(node_links->left != NULL) 
    fprint_tree(file, node_links->left, label_offset, link_offset, indent+1, printer);

  (*printer)(file, tree, indent, label_offset, link_offset);

  if(node_links->right != NULL) 
    fprint_tree(file, node_links->right, label_offset, link_offset, indent+1, printer);

}

void print_tree(void *tree, int label_offset, int link_offset, int indent, void ((*printer)())) {
  fprint_tree(stdout,tree,label_offset,link_offset,indent,printer);
}

/**********************************************************************
 * print_node -- dump a node to the output.  Gives minimal, generic 
 * information.  Use more specific dumping routines for more info.
 * Or write your own... :-)
 */

void fprint_node(FILE *f, void *foo,int indent, int label_offset, int link_offset) {
  int i;
  char buf[10240];
  char *bend = buf;
  
  *bend = 0;
  for(i=0;i<indent;i++) {
    *(bend++) = ' ';
  }
  *(bend) = '\0';
  
  fprintf(f 
	 , "%s%s links: %3.3d  label:%8.8d  total:%8.2g  flux:%8.2g\n"
	 , buf
	 , buf
	 , ((LINKS *)(foo + link_offset))->n
	 , *((long *)(foo + label_offset))
	 , ((LINKS *)(foo + link_offset))->sum
	 , *((NUM *)foo)
	 );
}

/**********************************************************************
 * print_all_fluxon_node
 * Prints out a node in the all_links tree for a particular fluxon.
 * (specific because the sum information depends on which set of links 
 * is being used). This is really just a convenience call to print_node.
 */
void fprint_all_fluxon_node(FILE *f, FLUXON *foo, int indent) {
  fprint_node(f, foo
	     , indent
	     , (long)(&(foo->label)) - (long)(foo)
	     , (long)(&(foo->all_links)) - (long)(foo)
	     );
}


/**********************************************************************
 * fprint_all_fluxon_tree
 * Prints out a fluxon tree linked through the all_links structure
 * (intended to print out all fluxons in the current world).  
 * Convenience call to print_tree.
 */

void fprint_all_fluxon_tree(FILE *f, FLUXON *foo) {
  int lab, lin;
  if(f==0) 
    f=stdout;
  lab = (long)&(foo->label) - (long)foo;
  lin = (long)&(foo->all_links) - (long)foo;
  fprint_tree(f, foo,lab,lin,0,fprint_node);
}  

void print_all_fluxon_tree(FLUXON *foo) {
  fprint_all_fluxon_tree(stdout,foo);
}


/**********************************************************************
 * fdump_fluxon
 * Prints out all the vertices in a fluxon.
 */
void fdump_fluxon(FILE *f, FLUXON *foo, int indent) {
  char buf[10240];
  char *bend = buf;
  VERTEX *v;
  int i;
  int j;

  if(f==0)
    f=stdout;

  *bend = 0;
  for(i=0;i<indent;i++) {
    *(bend++) = ' ';
  }
  *(bend) = '\0';

  fprint_all_fluxon_node(f,foo,indent);

  for((i=0),(v = foo->start); v ; v=v->next) {
    fprintf(f,"%s v %3d: lab %11d, loc %x x:(%6.3g, %6.3g, %6.3g) neigh: %3d, near: %3d, next: %x, prev: %x\n"
	    ,buf
	    , i
	    , v->label
	    , v
	    , v->x[0]
	    , v->x[1]
	    , v->x[2]
	    , v->neighbors.n
	    , v->nearby.n
	    , v->next
	    , v->prev
	    );
    fprintf(f,"\tNeighbors: ");

    for(j=0;j<v->neighbors.n;j++) 
      fprintf(f,"V%4d ",((VERTEX *)v->neighbors.stuff[j])->label);
    fprintf(f,"\n");

    if(v==foo->end) 
      fprintf(f,"%s%s   end\n",buf,buf);
    i++;
  }
}
	  
/**********************************************************************
 * fprint_all_fluxon_tree
 * Prints out a fluxon tree linked through the all_links structure
 * (intended to print out all fluxons in the current world).  
 * Convenience call to print_tree.
 */

void fdump_all_fluxon_tree(FILE *f, FLUXON *foo) {
  int lab, lin;
  lab = (long)&(foo->label) - (long)foo;
  lin = (long)&(foo->all_links) - (long)foo;
  if(f==0)
    f=stdout;
  fprint_tree(f, foo,lab,lin,0,fdump_fluxon);
}  

void dump_all_fluxon_tree(FLUXON *foo) {
  fdump_all_fluxon_tree(stdout,foo);
}



/**********************************************************************
 **********************************************************************
 ***** Line output printers -- 
 ***** Auxiliary routines for print_world.
 */

/****************************************
 * fprint_fc_line 
 * Flux concentration state output
 */
void fprint_fc_line(FILE *f, 
		    void *foo, int indent, int label_offset, int link_offset) {
  FLUX_CONCENTRATION *fc=foo;
  static char fc_line_format[80]="";
  if(!*fc_line_format) {
    sprintf(fc_line_format,"NEW\t%%ld\t%%%sg\t%%%sg\t%%%sg\t%%%sg\n",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
  }
  fprintf(f, fc_line_format,fc->label, fc->x[0], fc->x[1], fc->x[2], fc->flux);
}


/****************************************
 * fprint_fls_by_fc
 * Prints all the fluxons that START at a particular flux concentration.
 * Since positive flux concentrations start lines and negative ones end lines,
 * the distincion is easy to make.
 * Uses print_tree ...
 */
void fprint_fls_by_fc(FILE *f,
		      void *fc0,  int i, int label_offset, int link_offset) {
  FLUX_CONCENTRATION *fc = fc0;
  if(fc->flux >= 0) {
    fprintf(f,"# FC %ld\n",fc->label);
    fprint_tree(f, fc->lines, fl_lab_of, fl_start_ln_of, 0, fprint_fl_vertices);
  }
}

/****************************************
 * fprint_fl_vertices
 * Prints a fluxon definition line, and definition lines for all
 * the vertices in the fluxon.
 */
void fprint_fl_vertices(FILE *f,
			void *fl0, int i, int lab_of, int lin_of) {
  FLUXON *fl = fl0;
  static char fl_format[80]="";
  static char v_format[80]="";
  if(!*fl_format) {
    sprintf(fl_format,"LINE\t%%ld\t%%ld\t%%ld\t%%%sg\n",NUMCHAR);
    sprintf( v_format,"\tVERTEX\t%%ld\t%%ld\t%%ld\t%%%sg\t%%%sg\t%%%sg\n",
	     NUMCHAR,NUMCHAR,NUMCHAR);
  }


  fprintf(f, fl_format
	  , fl->label
	  , fl->fc0 ? fl->fc0->label : -1
	  , fl->fc1 ? fl->fc1->label : -1
	  , fl->flux);
  
  /* Print all the vertices (leave out the first and last) */
  {
    VERTEX *v = fl->start;
    int i = 1;

    for(; v && (v->next && (v->next != fl->end)); i++){
      v = v->next;
      fprintf(f, v_format
	      , v->line ? v->line->label : -1
	      , v->label
	      , i
	      , v->x[0]
	      , v->x[1]
	      , v->x[2]
	      );
    }
    fprintf(f,"\n");
  }
}
    
    
  
/**********************************************************************
 * print_dumblist
 */
void print_dumblist(DUMBLIST *foo, void ((*item_printer)())) {
  int i;

  if(!foo) 
    return;

  printf("Dumblist has %d items (size is %d):\n",foo->n, foo->size);
  for(i=0;i<foo->n;i++) {
    if(!item_printer) {
      printf("\t%d: %d\n",i,foo->stuff[i]);
    }
    else {
      (*item_printer)(i,foo);
    }
  }
}

/**********************************************************************
 **********************************************************************
 *****  GRAPHICS STUFF STARTS HERE
 *****  For now, this is a pretty stoopid diagnostic/debugging tool --
 *****  but eventually it will have several different output options.
 *****  You call gl_2d_start() to start writing a display script,
 *****  then you call the various gl_2d routines to dump display features
 *****  into the script.  Finally, you call gl_2d_finish() to close and
 *****  execute the display script.  Check individual routines for calling
 *****  conventions.

 *****  This is all pretty cheesy just now -- it should be handled
 *****  better and could almost certainly be handled with local gl calls.
 */ 


/**********************************************************************
 * gl_2d_start(FILENAME, l,r,b,t )
 * FILENAME is optional -- pass in 0 to make the routine
 * generate one for you.  
 * 
 * The l, r, b, and t parameters set the size of the viewport in 
 * scientific coordinates.  You need them.
 * 
 */
static FILE *gl_2d_file = 0;
static char gl_2d_fname[BUFSIZ+1];
static NUM gl_current_depth = 0;
static NUM gl_depth_step=1e-6;
static int gl_2d_fnum = 0;

int gl_2d_start(char *fname, float l, float r, float b, float t) {
  gl_current_depth = 0;

  if(gl_2d_file) {
    fprintf(stderr,"Warning: gl_2d_start abandoning unfinished file `%s'.\n",gl_2d_fname);
  }

  if(fname) {
    strncpy(gl_2d_fname,fname,BUFSIZ);
  } else {
    sprintf(gl_2d_fname,"/tmp/%d.graph%d",getpid(),gl_2d_fnum++);
  }

  if (!(gl_2d_file = fopen(gl_2d_fname,"w"))) /* assignment */ {
    fprintf(stderr,"Warning: gl_2d_start couldn't open file `%s' for writing\n",gl_2d_fname);
    return -1;
  }
  

  fprintf(gl_2d_file,"%s","#!/usr/bin/perl\n" \
"#\n" \
"# 2-D graphics dump file from FLEM \n" \
"#\n" \
"use OpenGL;\n" \
"glpOpenWindow(attributes=>[GLX_RGBA,GLX_DEPTH_SIZE,1],width=>600,height=>600,mask=>0x1ffffff);\n" \
"glEnable(GL_DEPTH_TEST);\n" \
"glClearColor(0.5,0.5,0.5,1);\n" \
"glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);\n" \
"glLoadIdentity;\n");

  fprintf(gl_2d_file,"glOrtho(%g,%g,%g,%g,-1,1);\n",l,r,b,t);


  fprintf(gl_2d_file,"sub draw {\n" \
"glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);\n");

  return 0;
}

/**********************************************************************
 * gl_2d_finish( display_time, display_mode )
 * Finish up a 2d gl dump file.  Returns the file name in a 
 * temporary variable.  The display_time parameter sets how long the
 * window is displayed, in seconds.  (Give it a negative number to wait
 * forever or until killed with the X manager).  The display_mode 
 * parameter tells what FLEM does with the output:
 *      0  -   do nothing
 *      1  -   run the display dump and wait for it to exit
 *      2  -   run the display dump and don't wait (killing it the next
 *             time gl_2d_finishes its run).
 */

static int gl_2d_pending = 0;
static char gl_2d_last[BUFSIZ+1];

char *gl_2d_finish( float display_time, int display_mode ) {
  struct stat stat_buf;

  if(!gl_2d_file) {
    fprintf(stderr,"gl_2d_finish:  Tried to finish a file without starting one!\n");
    *gl_2d_last = 0;
    return gl_2d_last;
  }

  /* Put the filename in a safe place */
  strncpy(gl_2d_last,gl_2d_fname,BUFSIZ);
  
  /* Write out the end condition for the window */
  fprintf(gl_2d_file, "glFlush();\n");
  fprintf(gl_2d_file,"\n}\n");

  fprintf(gl_2d_file,"\n\n\ndraw();\n");

  fprintf(gl_2d_file,"\n" \
"##Set up event handlers\n" \
"$cb{&ConfigureNotify} = sub{ my($e,$w,$h)=@_; };\n" \
"$cb{4}=$cb{&KeyPress}=sub{exit 0;}; # 4 is mouse down\n" \
"$cb{12}=$cb{&GraphicsExpose}=sub{ draw(); }; \n");

  if(display_time < 0) {
    fprintf(gl_2d_file,"\n" \
"while(1) {\n" \
"  while($p=XPending) {\n" \
"    @e=&glpXNextEvent;\n" \
"    if($s=$cb{$e[0]}) {\n" \
" #     print \"handling event @e\\n\";\n" \
"      &$s(@e);\n" \
"    } else {\n" \
"#      print \"event was @e\\n\";\n" \
"    }\n" \
"  }\n" \
"  sleep 0.1;\n" \
"} \n");

  } else {
    fprintf(gl_2d_file,"sleep %g\n",display_time);
  }

  /* Close the file and set permissions */
  fclose(gl_2d_file);
  gl_2d_file = 0;

  if(stat(gl_2d_fname,&stat_buf)) {
    fprintf(stderr,"gl_2d_finish: couldn't stat the file being finished (`%s')\n",gl_2d_fname);
    perror("");
    *gl_2d_last = 0;
    return gl_2d_last;
  }
  chmod(gl_2d_fname,  stat_buf.st_mode | S_IXUSR);

  /* Figure out what sort of subprocess to run */
  if(gl_2d_pending) {
    kill(gl_2d_pending,SIGKILL);
    gl_2d_pending = 0;
  }

  {
    int pid;
    int foo;

    switch(display_mode) {
    case 1:
      if(pid = fork())
	waitpid(pid,0,0);
      else {
	execl(gl_2d_fname,gl_2d_fname,0);
	fprintf(stderr,"exec failed (`%s')\n",gl_2d_fname);
	perror("");
	exit(-1);
      }
      break;

    case 2:
      if(pid=fork())
	gl_2d_pending = pid;
      else {
	execl(gl_2d_fname,gl_2d_fname,0);
	fprintf(stderr,"exec failed (`%s')\n",gl_2d_fname);
	perror("");
	exit(-1);
      }
      break;

    default:
      break;
    }
  }

  strncpy(gl_2d_last,gl_2d_fname,BUFSIZ);
  *gl_2d_fname=0;
  return gl_2d_last;
}

/**********************************************************************
 * gl_2d_point -- draw a "point" (actually, a little octagon) in the 
 * gl buffer
 */
void gl_2d_point(float x, float y, float radius, float colors[3]) {
  float rt2= 1.4142135623731;
  float x_off[9] = {-1, 1, 1+rt2 , 1+rt2, 1, -1, -1-rt2, -1-rt2, -1};
  float y_off[9] = {1+rt2,1+rt2,1,-1,-1-rt2,-1-rt2,-1,1,1+rt2};
  int i;

  if(!gl_2d_file) {
    fprintf(stderr,"gl_2d_point: no 2d output file is open!\n");
    return;
  }

  radius /= (1+rt2);  /* Make up for scale in the offsets */

  fprintf(gl_2d_file,"\n## gl_2d_point(%g,%g,%g)\n",x,y,radius);
  fprintf(gl_2d_file,"glColor3f(%g,%g,%g);\n",colors[0],colors[1],colors[2]);
  fprintf(gl_2d_file,"glBegin(GL_POLYGON);\n");
  for(i=0;i<9;i++) 
    fprintf(gl_2d_file,"\tglVertex3f(%g,%g,%g);\n",
	    x + radius * x_off[i],
	    y + radius * y_off[i],
	    (gl_current_depth += gl_depth_step)
	    );
  fprintf(gl_2d_file,"glEnd();\n");

  return;
}    
  
void gl_2d_scr_poly(DUMBLIST *horde, float colors[3]) {
  int i;
  if(!gl_2d_file) {
    fprintf(stderr,"gl_2d_scr_poly: no 2d output file!\n");
    return;
  }

  fprintf(gl_2d_file,"\n\n##gl_2d_vertex_scr\n");
  fprintf(gl_2d_file,"glColor3f(%g,%g,%g);\n",colors[0],colors[1],colors[2]);
  fprintf(gl_2d_file,"glBegin(GL_POLYGON);\n");
  for(i=0;i<horde->n;i++) {
    fprintf(gl_2d_file,"\tglVertex3f(%g,%g,%g);\n",
	    ((VERTEX *)((horde->stuff)[i]))->scr[0],
	    ((VERTEX *)((horde->stuff)[i]))->scr[1],
	    (gl_current_depth += gl_depth_step));
  }
  fprintf(gl_2d_file,"\tglVertex3f(%g,%g,0);\n",
	  ((VERTEX *)((horde->stuff)[0]))->scr[0],
	  ((VERTEX *)((horde->stuff)[0]))->scr[1]);
	  
  fprintf(gl_2d_file,"glEnd();\n");
}

void gl_2d_line(float x0, float y0, float x1, float y1, float colors[3]) {
  if(!gl_2d_file) {
    fprintf(stderr,"gl_2d_line: no 2d output file is open!\n");
    return;
  }

  fprintf(gl_2d_file,"\n\n## gl_2d_line(%g,%g,%g,%g)\n",x0,y0,x1,y1);
  fprintf(gl_2d_file,"glColor3f(%g,%g,%g);\n",colors[0],colors[1],colors[2]);
  fprintf(gl_2d_file,"glBegin(GL_LINE_STRIP);\n");
  fprintf(gl_2d_file,"\tglVertex3f(%g,%g,%g);\n",x0,y0,(gl_current_depth += gl_depth_step));
  fprintf(gl_2d_file,"\tglVertex3f(%g,%g,%g);\n",x1,y1,(gl_current_depth += gl_depth_step));
  fprintf(gl_2d_file,"glEnd();\n");

  return;
}
  

void gl_2d_scr_list(DUMBLIST *horde,float colors[3]) {
  int i;
  if(!gl_2d_file) {
    fprintf(stderr,"gl_2d_scr_list: no 2d output file!\n");
    return;
  }

  for(i=0;i<horde->n;i++) {
    gl_2d_point(
		((VERTEX *)((horde->stuff)[i]))->scr[0],
		((VERTEX *)((horde->stuff)[i]))->scr[1],
		1,
		colors
		);
  }
}
    

