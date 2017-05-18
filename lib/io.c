/**********************************************************************
 * io.c -- I/O routines for FLUX
 *
 * Mostly text I/O.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Southwest Research Institute, 2004-2008
 * 
 * You may modify and/or distribute this software under the terms of
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
 * This file is part of the FLUX 2.2 release (22-Nov-2008).
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <regex.h>
#include <math.h>

#include "data.h"
#include "io.h"
#include "model.h" 
#include "physics.h" /* for FORCES translation */

char *code_info_io="%%%FILE%%%";

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
static char mscan[80]="";
static char lscan[80]="";
static char lscan2[80]="";
static char vscan[80]="";
static char gscan[80]="";
static char vnscan[80]="";

int footpoint_action(WORLD *world, char *s) {
  char *ds = 0;
  char *badstr = 0;
  char c;
  int i,n;
  long vl0, vl1; // scratch vertex labels
  long vertex_label;
  long l0,l1,fl0;
  NUM x0[3],x1[3],flux0,flux1;

  /* Skip over initial whitespace */
  while(*s && isspace(*s))
    s++;

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
      n=sscanf(s,"%*s %d",&i);
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
             * NEW  <label> <x> <y> <z> <flux> <boundary> <locale_radius>
	     */
    /* Initialize the scan string, if necessary */
    if(!*mscan) {
      sprintf(mscan,"%%*s %%ld %%%sf %%%sf %%%sf %%%sf %%s %%%sf",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
    }
    {
      char boundary_name[BUFSIZ];
      NUM radius;

      n = sscanf(s,mscan, &l0, x0, x0+1, x0+2, &flux0, boundary_name, &radius);
      
      if(n != 5 && !(n==7 && toupper(*s)=='N')) {
	badstr = "Couldn't parse NEW line";
      } else {
	if(toupper(*s)=='N') {
	  void *bptr;
	  
	  if(n==7) {
	    bptr = boundary_name_to_ptr(boundary_name);
	    //	    if(!bptr) {
	    //  badstr = "invalid boundary name in NEW line";
	    //  goto global_escape;
	    //  }
	  }
	  
	  if( l0<0 ) {
	    printf("WARNING: automagic FC %ld is given in file -- ignoring this gaffe. (Probably OK)\n",l0);
	    break;
	  }
	  
	  
	  /* New stuff handler */
	  
	  FLUX_CONCENTRATION *fc;
	  
	  fc = new_flux_concentration( world,x0[0],x0[1],x0[2],flux0, l0 );

	  if(n==7) {
	    fc->bound = bptr;
	    fc->locale_radius = radius;
	  }

	} else if(toupper(*s)=='M') {
	  
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
	} else {
	  badstr = "NEW/MOVE line is neither a NEW nor a MOVE!  This should never happen.\n";
	}
      }
    } // end of convenience block 

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
  case 'L': /* 
	     * LINE <fl_lab> <vlab1> <vlab2> <lab1> <lab2> <flux> <sx> <sy> <sz> <ex> <ey> <ez> (12-parameter form)
	     * LINE <fl_lab> <vlab1> <vlab2> <lab1> <lab2> <flux>  (6-parameter form)
	     * LINE <fl_lab> <lab1> <lab2> <flux> (4 parameter form)
	     *  Create a field line from 1 to 2 with label and flux given 
	     */
    if(!*lscan) {
      sprintf(lscan,"%%*s %%ld %%ld %%ld %%ld %%ld %%%sf %%%sf %%%sf %%%sf %%%sf %%%sf %%%sf",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
    }
    
    {
      NUM sx, sy, sz, ex, ey, ez;
      
      n = sscanf(s,lscan, &fl0, &vl0, &vl1, &l0, &l1, &flux0, &sx, &sy, &sz, &ex, &ey, &ez);
      
      if(n < 6) {
      // Fall back to older form 
	
	if(!*lscan2) {
	  sprintf(lscan2,"%%*s %%ld %%ld %%ld %%%sf",NUMCHAR);
	}
	
	n = sscanf(s,lscan2, &fl0, &l0, &l1, &flux0);
	
	if(n != 4) {
	  badstr = "Couldn't parse LINE line";
	  n=0;
	} else {
	  // n==4
	  vl0 = vl1 = 0;
	}
      }
      
      if(n != 0) {
      //fprintf(stderr,"We are at another line now %ld\n",fl0); //debugging
	
	if(fl0 <= -11 || fl0 >0) {
	  FLUX_CONCENTRATION *fc0, *fc1;
	  fc0 = tree_find(world->concentrations, l0, fc_lab_of, fc_ln_of);
	  fc1 = tree_find(world->concentrations, l1, fc_lab_of, fc_ln_of);
	  
	  if(!fc0 || !fc1) {
	    char *badbuf = (char *)localmalloc(BUFSIZ,MALLOC_MISC);
	    sprintf(badbuf,"Found a fluxon specifier between concentrations %ld and %ld, but they \ncame up %ld and %ld in tree_find (one doesn't exist)!  \nThis error message leaked %d bytes (don't let it happen again!)\n",l0,l1,(long)fc0,(long)fc1,BUFSIZ);
	    badstr = badbuf;
	    
	  } else if(fc0->flux * fc1->flux >= 0) {
	    badstr = (char *)localmalloc(BUFSIZ,MALLOC_MISC);
	    sprintf(badstr,"This fluxon connects two flux concentrations of the same sign flux, or\none of its flux tubes has zero flux. Line %ld (flux %g); concentrations %ld (%g) and %ld (%g)\n",fl0, flux0, l0, fc0->flux, l1, fc1->flux);
	    break;
	  }	else {
	    /* Check if the field line exists in either of the two 
	       concentrations' local lists */
	    FLUXON *f;
	    
	    if((f = tree_find(world->lines,fl0,fl_lab_of,fl_all_ln_of)) &&
	       (f->fc0 != fc0 || f->fc1 != fc1)
	       ) {
	      badstr = (char *)localmalloc(BUFSIZ,MALLOC_MISC);
	      sprintf(badstr,"Hey!  fluxon %ld already exists (from conc. %ld - %ld)\n",f->label,f->fc0->label, f->fc1->label); 
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
		
		if(fc0->label==-3 && fc1->label==-4) {
		  f0->plasmoid = 1;
		} else if(fc0->label==-3 || fc1->label==-4) {
		  badstr = (char *)localmalloc(BUFSIZ,MALLOC_MISC);
		  sprintf(badstr,"Fluxon %ld uses one (but not both) of the reserved plasmoid FC's (-3 and -4)! Not allowed.",f->label);
		}
		
		if(!badstr) {
		  /* Make sure it has at least the two end vertices */
		  VERTEX *v0, *v1;
		  if(vl0==0) 
		    vl0 = -(f0->label*2 + 100);
		  if(vl1==0) 
		    vl1 = -(f0->label*2 + 100)+1;
		  
		  if( n==6 ) {
		    sx = fc0->x[0]; sy = fc0->x[1]; sz = fc0->x[2];
		    ex = fc1->x[0]; ey = fc1->x[1]; ez = fc1->x[2];
		  }
		  
		  v0 = new_vertex( vl0, sx, sy, sz, f0);
		  v1 = new_vertex( vl1, ex, ey, ez, f0);
		  
		  if(!v0 || !v1) {
		    badstr = "Couldn't make trivial vertices for fluxon!\n";
		  } else {
		    if( add_vertex_pos(f0,0,v0) || add_vertex_pos(f0,-1,v1) ) 
		      badstr = "Problems with trivial vertex addition\n";
		  } // end of vertex-creation check
		}// end of error check
	      } // end of fluxon existence & flux concentration consistency check
	    } // end of fluxon existence check
	  } // end of sign consistency check
	} // end of label check
      } // end of parsing check
    } // end of 'L' case
    
    break;
    
  case 'V': 

    if(s[1]=='E' || s[1]=='e') {

      /* VERTEX <fl_lab> <vert_lab> <pos> <x> <y> <z> 
       *  Create a new vertex on field line <fl_lab> at position
       *  <pos> (0 means start; -1 means end) along the existing one.
       */
      if(!*vscan) {
	sprintf(vscan,"%%*s %%ld %%ld %%ld %%%sf %%%sf %%%sf",NUMCHAR,NUMCHAR,NUMCHAR);
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
	  if(world->verbosity > 2) printf("vertex_label = %ld; v0->label = %ld\n",vertex_label, v0->label);
	  if(add_vertex_pos(f, l1, v0))
	    badstr = "VERTEX: add_vertex_pos returned error";
	}
      }      
    
    } else if(s[1]=='_' || s[1]=='N' || s[1]=='n') {
      VERTEX *v;
      VERTEX *n;
      long vlab;
      long nlab;
      int n_scan;

      /* V_NEIGHBOR <vert_lab> <vn_lab> */
      //      printf("(V_NEIGHBOR not yet supported...)\n");
      n_scan = sscanf(s,"%*s %ld %ld",&vlab,&nlab);
      if(n_scan!=2) {
	badstr = "Couldn't parse V_NEIGHBOR line";
      } else if( nlab != 1 && nlab != 2 ) {
	v = tree_find(world->vertices, vlab, v_lab_of, v_ln_of);
	n = tree_find(world->vertices, nlab, v_lab_of, v_ln_of);
	if(!v || !n) {
	  fprintf(stderr,"Warning: V_NEIGHBOR line asked for targets %ld and %ld; at least one of 'em isn't defined...\n",vlab,nlab);
	} else {
	  dumblist_add( &(v->neighbors), n );
	  dumblist_add( &(n->nearby),    v );
	}
      }
    } else {
      badstr = "V-huh?";
    }

    break;
  
  case 'G': /* GLOBAL <cmd> <arg1> <arg2> ... */
    {
      char gcmd[81];
      int off,n;
      int fred;
      NUM bert;

      n = sscanf(s,"%*s %80s %n",gcmd,&off);
      if(n != 1) {
	static char bstr[80];
	sprintf(bstr,"Couldn't parse GLOBAL line (n=%d)\n",n);
	badstr = bstr;
      } else {
	switch(*gcmd) {
	case 'C':  /* GLOBAL CONCURRENCY or GLOBAL COEFFICIENTS */
	  if(gcmd[1]=='O' && gcmd[2]=='E') {/* GLOBAL COEFFICIENTS */
	    int num_to_read, i, off, off2;
	    double coeff;
	    if (sscanf(s, "%*s %*s %d %n", &num_to_read, &off) != 1) 
	      badstr = "Couldn't parse number of coefficients to read";
	    for (i=0;i<num_to_read;i++) {
	      if (sscanf(s+off, "%lf %n", &coeff, &off2) !=1) {
		badstr = "Couldn't parse world->coeffs";
	      }
	      off += off2;
	      world->coeffs[i] = coeff;
	    }
	    world->n_coeffs = num_to_read;
	  }
	  else if(gcmd[1]=='O' && gcmd[2]=='N') { /* GLOBAL CONCURRENCY */
	    int nval;
	    int concurrency;
	    if( 1 > (nval = sscanf(s,"%*s %*s %d",&concurrency))) {
	      badstr = "couldn't parse GLOBAL CONCURRENCY line";
	    } else {
	      world->concurrency = concurrency;
	    }
	  }
	  else  {
	    badstr = "Couldn't parse GLOBAL C... (CONCURRENCY? COEFFICIENTS?)";
	  }
	  break;
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
		unsigned long foo;
		if( sscanf(fname,"0x%lx",&foo) )
		  f_funcs[i] = (void *)foo;
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

	      /* f_funcs[i] = (void *)0; */

	      (f_funcs[i]) = (void *)0;


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

	case 'B':  /* BFIELD_NORMALIZATION or BOUNDARY */
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
	      /* Cylinder, radius is length of normal */
	      sprintf(phscan,"%%*s %%%sf %%%sf %%%sf %%%sf %%%sf %%%sf",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
	      if(6 > sscanf(s+off,phscan,
			    &(p->origin[0]),
			    &(p->origin[1]),
			    &(p->origin[2]),
			    &(p->normal[0]),
			    &(p->normal[1]),
			    &(p->normal[2])
			    )
		 ) {
		badstr = "Couldn't parse cylindrical radius for GLOBAL B2 CYLINDER";
	      }
	      type_code = PHOT_CYL;
	      break;

	    case '0':
	    case 'N':
	      /* None */
	      free(p);
	      p=0;
	      type_code = 0;
	      break;

	    default:
	      {
		static char badbuf[250];
		sprintf(badbuf,"Unknown boundary condition in GLOBAL BOUNDARY - off is %d, s[off] is '%c'; line is '%s'\n",off,s[off],s);
		badstr = badbuf;
	      }
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
 	  } else if(gcmd[1] == '2') { 

	    /* GLOBAL B2, second boundary */
	    fprintf(stderr,"second photosphere detected\n");
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
		badstr = "Couldn't parse six values from GLOBAL B2 PLANE";
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
		badstr = "Couldn't parse four values from GLOBAL B2 SPHERE";
	      p->normal[1] = 0;
	      p->normal[2] = 0;
	      type_code = PHOT_SPHERE;
	      break;
	      
	    case PHOT_CYL+'0':
	    case 'C':
	      /* Cylinder */
	      sprintf(phscan,"%%*s %%%sf %%%sf %%%sf %%%sf %%%sf %%%sf",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
	      if(6 > sscanf(s+off,phscan,
			    &(p->origin[0]),
			    &(p->origin[1]),
			    &(p->origin[2]),
			    &(p->normal[0]),
			    &(p->normal[1]),
			    &(p->normal[2])
			    )
		 ) {
		badstr = "Couldn't parse cylindrical radius for GLOBAL B2 CYLINDER";
	      }
	      type_code = PHOT_CYL;
	      break;

	    case '0':
	    case 'N':
	      /* None */
	      free(p);
	      p=0;
	      type_code = 0;
	      break;

	    default:
	      {
		static char badbuf[250];
		sprintf(badbuf,"Unknown boundary condition in GLOBAL B2 - off is %d, s[off] is '%c'; line is '%s'\n",off,s[off],s);
		badstr = badbuf;
	      }
	      break;
	    }

	    if(!badstr) { /* Parse was okay */
	      world->photosphere2.type = type_code;
	      if(world->photosphere2.plane) 
		localfree(world->photosphere2.plane);
	      world->photosphere2.plane = p;

	    } else { /* Clean up so we don't leak too much memory */
	      if(p)
		free(p);
	    }

	    /* end of GLOBAL B2*/
	  } else {

	    fprintf(stderr,"WARNING: deprecated 'GLOBAL B_FLAG' line found in file...\n");
	    /** GLOBAL B_NORMALIZATION ... **/
	    int i;
	    sscanf(s+off,"%d",&i) || (badstr="couldn't parse F_OVER_B state");
	    
	    if(i==0) {
	      fprintf(stderr,"\t(setting b_power=0,d_power=2,s_power=0, ds_power=0)\n");
	      world->step_scale.b_power = 0;
	      world->step_scale.d_power = 2;
	      world->step_scale.s_power = 0;
	      world->step_scale.ds_power = 0;
	    } else {
	      fprintf(stderr,"\t(setting b_power=1,d_power=2,s_power=0, ds_power=0)\n");
	      world->step_scale.b_power = 1;
	      world->step_scale.d_power = 2;
	      world->step_scale.s_power = 0;
	      world->step_scale.ds_power = 0;
	    }

	  }
	  break;
	case 'O': /* OPEN <x> <y> <z> <r> <auto> - where is the open sphere, and should open lines be handled automagically? */
	  {
	    char phscan[80];
	    sprintf(phscan, "%%%sf %%%sf %%%sf %%%sf %%d",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);

	    if(5 > sscanf(s+off,phscan,
			  &(world->fc_oe->x[0]),
			  &(world->fc_oe->x[1]),
			  &(world->fc_oe->x[2]),
			  &(world->fc_oe->locale_radius),
			  &(world->auto_open)
			  )
	       ) {
	      badstr = "Couldn't parse five values from GLOBAL OPEN ";
	    } else {
	      world->fc_ob->x[0] = world->fc_oe->x[0];
	      world->fc_ob->x[1] = world->fc_oe->x[1];
	      world->fc_ob->x[2] = world->fc_oe->x[2];
	      world->fc_ob->locale_radius = world->fc_oe->locale_radius;
	    }
	  }
	  break;
	case 'R':  /* ARD Added - read rel_step */
	  if (sscanf(s, "%*s %*s %ld", &(world->rel_step)) != 1) {
	    badstr = "Couldn't parse world->rel_step";
	  }
	  break;
	case 'D':  /* ARD Added - read_dtau */
	  if (sscanf(s, "%*s %*s %lf", &(world->dtau)) != 1) {
	    badstr = "Couldn't parse world->dtau";
	  }
	  break;
	case 'S': /* GLOBAL STATE  or GLOBAL SCALING or GLOBAL SKEW */

	  switch(*(gcmd+1)) {
	  case 'C': case 'c': /* GLOBAL SCALING */
	    {
	      int ok_to_scan = 1;
	      void *scanvar;

	      switch(*(s+off)) {
	      case 'B': case 'b':
		scanvar = &(world->step_scale.b_power);
		break;
	      case 'D': case 'd':
		switch(*(s+off+1)) {
		case 'S': case 's':
		  scanvar = &(world->step_scale.ds_power);
		  break;
		case 0: case ' ': case '_': case '\t':
		  scanvar = &(world->step_scale.d_power);
		  break;
		default:
		  ok_to_scan = 0;
		  break;
		}
		break;
	      case 'S': case 's':
		scanvar = &(world->step_scale.s_power);
		break;
	      default:
		ok_to_scan = 0;
		break;
	      }

	      if(ok_to_scan) {
		//		char mscan[100];
		//		sprintf(mscan,"%%*s %%%sf",NUMCHAR);
		//		sscanf(mscan,scanvar);
		sscanf(s+off,"%*s %"NUMCHAR"f", (double *)scanvar);
	      } else {
		badstr = s;
	      }

	    }
	    break;

	  case 'T': case 't': /* GLOBAL STATE */
	    sscanf(s+off,"%ld",&(world->state)) || 
	      (badstr="couldn't parse GLOBAL STATE");
	    break;

	  case 'K': case 'k': /* GLOBAL SKEW */
	    sscanf(s+off,"%d",&(world->handle_skew)) ||
	      (badstr = "couldn't parse GLOBAL SKEW");
	    break;

	  }
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
  int seekable = 0;
  int reporting = 0;
  long startpos;
  long endpos = 0;
  long last_report = 0;
  long pos;

  // Check if it's seekable
  startpos = ftell(file);
  if(startpos >= 0) {
    if(fseek(file, startpos, SEEK_SET) >= 0) {
      seekable = 1;
    }
  }

  // If seekable, scan ahead to find the end of the file
  if(seekable) {
    fseek(file, 0, SEEK_END);
    endpos = ftell(file);
    fseek(file, startpos, SEEK_SET);
  }

  // Report progress on files greater than a megabyte
  if(seekable && (endpos - startpos > 1024*1024)) {
    printf("FLUX: Reading %g MB...\n",(double)(endpos-startpos)/1024/1024);
    reporting = 1;
    last_report = startpos;
  }
  
  do{
    if( a == NULL ) {
      a = new_world();
    }

    if( (s = next_line(file)) ) {       // assignment
      error = footpoint_action(a, s);
    }

    if(seekable && reporting) {
      pos = ftell(file);
      if(last_report==startpos || (pos - last_report >= 1024 * 100)) {
	printf("\r%4.2f%% completed (%5.2f MB), %ld fluxons, %ld vertices  ",100 * (double)(pos - startpos)/(double)(endpos - startpos), (double)(pos - startpos) / 1024 / 1024, a->lines->all_links.n, a->vertices->world_links.n);
	fflush(stdout);
	last_report = pos;
      }
    }
    
  } while (s && !error);

  if(seekable && reporting) {
    printf("done\n");
  }


  world_update_ends(a);

  a->vertices = tree_balance(a->vertices, v_lab_of, v_ln_of);

  world_check(a);

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
  int i, n;

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
      fprintf(file,"0x%lx ",(unsigned long)(world->f_funcs[i]));
  }
  fprintf(file,"\n");

  fprintf(file, "GLOBAL STATE %ld  (%s)\n",world->state,world_state_name(world));

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

  {
    PLANE *p2 = world->photosphere2.plane;
    /* photosphere2 */
    char fmt[80];
    sprintf(fmt,"%%s %%s %%%sg %%%sg %%%sg %%%sg %%%sg %%%sg\n",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
    fprintf(file,fmt
	    , "GLOBAL B2 "
	    , BOUNDARY_NAMES[ world->photosphere2.type ]
	    , p2?p2->origin[0]:0, p2?p2->origin[1]:0, p2?p2->origin[2]:0
	    , p2?p2->normal[0]:0, p2?p2->normal[1]:0, p2?p2->normal[2]:0
	    );
  }

  fprintf(file,"GLOBAL SCALING B %g\n", world->step_scale.b_power);
  fprintf(file,"GLOBAL SCALING D %g\n", world->step_scale.d_power);
  fprintf(file,"GLOBAL SCALING S %g\n", world->step_scale.s_power);
  fprintf(file,"GLOBAL SCALING DS %g\n", world->step_scale.ds_power);
  
  /*  ARD 020207 Added writing code fgor RSTEP and DTAU */

  fprintf(file,"GLOBAL RSTEP %ld\n", world->rel_step);
  fprintf(file,"GLOBAL DTAU %lf\n", world->dtau);

  /* ARD Added code for writing out global coefficients */
 
  fprintf(file,"GLOBAL COEFFICIENTS %ld", world->n_coeffs);
  for (i=0;i<world->n_coeffs;i++) {
    fprintf(file, " %f",  world->coeffs[i]); 
  }
  fprintf(file,"\n");
  fprintf(file,"GLOBAL OPEN %g %g %g %g %ld\n",
	  world->fc_ob->x[0],
	  world->fc_ob->x[1],
	  world->fc_ob->x[2],
	  world->fc_ob->locale_radius,
	  world->auto_open
	  );

  fprintf(file,"GLOBAL SKEW_HANDLING %d\n", world->handle_skew);

  fprintf(file,"GLOBAL CONCURRENCY %ld\n",world->concurrency);
  
  fprintf(file,"\n\n"); /* leave an extra space after the globals */
    
  /* Write out all flux concentration locations */
  fprint_tree(file, world->concentrations, fc_lab_of, fc_ln_of, 0, fprint_fc_line_nonneg);

  /* Output the field lines */
  fputc('\n',file);
  fprint_tree(file, world->concentrations, fc_lab_of, fc_ln_of, 0, fprint_fls_by_fc);

  /* Output neighbor relations */
  fputc('\n',file);
  fprint_tree(file, world->vertices, v_lab_of, v_ln_of, 0, fprint_v_nbors);

  return 1;
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
    //    printf("<<NO TREE>>\n");fflush(stdout);
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
	 , "%s%s links: %3.3ld  label:%8.8ld  total:%8.2g  flux:%8.2g\n"
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
    fprintf(f,"%s v %3d: lab %11ld, loc %lx x:(%6.3g, %6.3g, %6.3g) neigh: %3d, near: %3d, next: %lx, prev: %lx\n"
	    ,buf
	    , i
	    , v->label
	    , (unsigned long)v
	    , v->x[0]
	    , v->x[1]
	    , v->x[2]
	    , v->neighbors.n
	    , v->nearby.n
	    , (unsigned long)v->next
	    , (unsigned long)v->prev
	    );
    fprintf(f,"\tNeighbors: ");

    for(j=0;j<v->neighbors.n;j++) 
      fprintf(f,"V%4ld ",((VERTEX *)v->neighbors.stuff[j])->label);
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


/******************************
 * fprint_v_nbors 
 * Print the neighbors of a vertex in a parsable format
 */
void fprint_v_nbors(FILE *f,
		    void *foo, int indent, int label_offset, int link_offset) {
  VERTEX *v = foo;
  int i;
  for(i=0;i<v->neighbors.n;i++) {
    fprintf(f, 
	    "VNEIGHBOR %ld %ld\n", 
	    v->label, 
	    (((VERTEX **)(v->neighbors.stuff))[i])->label);
  }
}
 
/****************************************
 * fprint_fc_line 
 * Flux concentration state output
 */
void fprint_fc_line(FILE *f, 
		    void *foo, int indent, int label_offset, int link_offset) {
  FLUX_CONCENTRATION *fc=foo;
  static char fc_line_format[80]="";
  if(!*fc_line_format) {
    sprintf(fc_line_format,"NEW\t%%ld\t%%%sf\t%%%sf\t%%%sf\t%%%sg\t%%s\t%%%sg\n",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
  }
  fprintf(f, fc_line_format,fc->label, fc->x[0], fc->x[1], fc->x[2], fc->flux, boundary_ptr_to_name(fc->bound), fc->locale_radius );
}

void fprint_fc_line_nonneg(FILE *f, 
		    void *foo, int indent, int label_offset, int link_offset) {
  FLUX_CONCENTRATION *fc=foo;
  if( fc->label >= 0 )
    fprint_fc_line(f,foo,indent,label_offset,link_offset);
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
  if(fc->lines && (fc->lines->fc0 == fc)) {
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
    sprintf(fl_format,"LINE\t%%ld\t%%ld\t%%ld\t%%ld\t%%ld\t%%%sg            %%%sg %%%sg %%%sg   %%%sg %%%sg %%%sg\n",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
    sprintf( v_format,"\tVERTEX\t%%ld\t%%ld\t%%ld\t%%%sf\t%%%sf\t%%%sf\n",
	     NUMCHAR,NUMCHAR,NUMCHAR);
  }


  fprintf(f, fl_format
	  , fl->label
	  , fl->start ? fl->start->label : 0
	  , fl->end ? fl->end->label : 0
	  , fl->fc0 ? fl->fc0->label : -999
	  , fl->fc1 ? fl->fc1->label : -999
	  , fl->flux
	  , fl->start->x[0]
	  , fl->start->x[1]
	  , fl->start->x[2]
	  , fl->end->x[0]
	  , fl->end->x[1]
	  , fl->end->x[2]
	  );
  
  /* Print all the vertices (leave out the first and last) */
  {
    VERTEX *v = fl->start;
    int i = 1;

    for(; v && (v->next && (v->next != fl->end)); i++){
      v = v->next;
      fprintf(f, v_format
	      , v->line ? v->line->label : -999
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
      printf("\t%d: %lx\n",i,(unsigned long)(foo->stuff[i]));
    }
    else {
      (*item_printer)(i,foo);
    }
  }
}


/**********************************************************************
 **********************************************************************
 * Binary serialization / activation
 *
 * These routines implement a binary file format for FLUX, but it is not
 * complete in the sense of replacing the ASCII .flux files generated by
 * print_world and read by read_world.  It is more compact, hence more 
 * efficient, and uses machine-dependent binary architecture.  It is 
 * primarily used for parallelization.
 *
 * The file format is tagged for easy extension.  
 * It is a collection of tags: Fence, ID, offset to next ID.  
 * The tags repeat until an end tag is encountered.
 * 
 * To generate a file, just start sending packets of data into an open
 * file descriptor, and terminate with the end field (using binary_dump_end).
 *
 * If you want to store a portable file, you should start by calling
 * binary_dump_header(), which will send a check header to ensure 
 * binary architecture compatibility.  
 *
 * To read a file, call binary_read_dumpfile with an open descriptor and the 
 * WORLD to which it should apply.  It will read and parse all the elements 
 * in the file, modifying the arena as they are processed, and close the 
 * file descriptor on exit.  If you send NULL instead of a valid WORLD, 
 * then binary_read_dumpfile will allocate a WORLD for you.
 *
 *
 */

/**********************************************************************
 * Dump modules - you deliver a file descriptor and whatever else you're
 * dumping, and it gets sent to the file descriptor, along witha  pointer
 * to the next tag.n
 */


/******************************
 * binary_dump_field -- sends the header and serialized data out to the 
 * dump file (or pipe).
 */
int binary_dump_field(int fd, long code, long len, char *buf) {
  long hdrbuf[3];
  long i;

  hdrbuf[0] = BD_FENCE;
  hdrbuf[1] = code;
  hdrbuf[2] = len;
  
  i = write(fd, hdrbuf, 3*sizeof(long));
  
  if( i != 3*sizeof(long) )  {
    perror("binary_dump_field");
    return 1;
  }
  
  if(len > 0) {
    i = write( fd, buf, len );
    if( i != len )
      return 2;
  }

  return 0;
}
    
/******************************
 * dump routines - each one sends one packet out to the dump file, containing 
 * its own type of data. Currently there are fluxon_forces (dumping all the vertex
 * force data for an entire fluxon)
 */


/**********
 * Local utility: a single I/O buffer that grows to a high-water mark
 */
static char *binary_buffer = 0;
static long binary_buflen = 0;
static int check_binary_buf( long desired_len ) {
  if(!binary_buffer) {
    binary_buffer = (char *)malloc(10240);
    binary_buflen = 10240;
    if(!binary_buffer) {
      fprintf(stderr,"Failed to get an I/O buffer (tried for %d bytes) in FLUX library - check_binary_buf",10240);
      exit(1);
    }
  }
  if(desired_len > binary_buflen) {
    if(binary_buffer)
      free(binary_buffer);
    binary_buffer = (char *)malloc(desired_len * 2);
    binary_buflen = desired_len * 2;
    if(!binary_buffer) {
      fprintf(stderr,"Failed to grow the binary I/O buffer - tried for %ld bytes",binary_buflen);
      exit(1);
    }
  }
  return 1;
}

/******************************
 * binary_dump_header  & binary_read_header - send/check  a start marker in a file.  
 * Not strictly necessary for local dump/restore, but allows endianness and data size checking.
 * 
 * The header format is several unsigned longs:
 *    0:  0xAABBCCDD - check endianness
 *    1:  version number (currently 1)
 *    2:  sizeof(NUM)
 *    3:  3.1416 in floating point
 */
int binary_dump_header(int fd) {
  long *foo;
  check_binary_buf( 4 * sizeof(long) );
  foo = (long *)binary_buffer;
  foo[0] = 0xAABBCCDD;
  foo[1] = 1;
  foo[2] = sizeof(NUM);
  *(float *)(&(foo[3])) = 3.1416;
  binary_dump_field(fd, BD_HDR, 4 * sizeof(long), binary_buffer);
  return 1;
}

int binary_read_header(long size, char *buf, WORLD *w) {
  char *me = "binary_read_header";
  long *foo = (long *)buf;
  float f;

  if(foo[0] != 0xAABBCCDD) {
    fprintf(stderr,"%s, endian check: expected %lx, got %lx (NUXI problem?)\n",me, (unsigned long)(0xAABBCCDD),(unsigned long)(foo[0]));
    return 1;
  }
  if(foo[1] != 1) {
    fprintf(stderr,"%s, version check: expected %d, got %ld (oops)\n",me, 1,foo[1]);
    return 2;
  }

  if(foo[2] != sizeof(NUM)) {
    fprintf(stderr,"%s, precision check: sizeof(NUM) is %ld, but file claims %ld (oops)\n",me, sizeof(NUM), foo[2]);
    return 3;
  }
  
  f = *(float *)(&(foo[3]));
  if(f != 3.1416) {
    fprintf(stderr,"%s, binary format check: expected 3.1416, got %g",me,f);
    return 4;
  }
  return 0;
}


/******************************
 * binary_dump_WORLD
 * 
 * Doesn't dump the whole shebang, only the global stuff 
 * (and none of the major data structures such as the
 * concentration, vertex, or fluxon trees).
 */
int binary_dump_WORLD(int fd, WORLD *w) {
  char *ptr;
  char *s;
  int i,j;

  // Allocate buffer space.   As always, overkill is the norm.
  check_binary_buf( sizeof(WORLD) +        // Room for structure fields 
		    sizeof(long)*8  +      // Room for a few extras
		    80 * (1 + N_FORCE_FUNCS*2 + N_RECON_FUNCS + N_RECON_PARAMS + N_M_PARAMS));
  ptr = binary_buffer;

  *(long *)ptr = 1;  // WORLD dump version number
  ptr += sizeof(long);
  
  *(long *)ptr = sizeof(WORLD);         //sizeof(WORLD) structure as an additional check
  ptr += sizeof(long);
  
  // Copy the WORLD structure, in its entirety
  for(i=0, s=(char *)w; i<sizeof(WORLD); i++)
    *(ptr++) = *(s++);
  for(; i%16; i++) 
    *(ptr++) = 0;

  
  // Now patch it up.


  // Leave a fence before the variable part
  *(long *)ptr = WORLD_VAR_FENCE; 
  ptr += sizeof(long);

  // skip the concentration, line, and vertex trees...

  // fill in photospheric planes...
  if(w->photosphere.plane) {
    *(NUM *)ptr = w->photosphere.plane->origin[0];   ptr += sizeof(NUM);
    *(NUM *)ptr = w->photosphere.plane->origin[1];   ptr += sizeof(NUM);
    *(NUM *)ptr = w->photosphere.plane->origin[2];   ptr += sizeof(NUM);
    *(NUM *)ptr = w->photosphere.plane->normal[0];   ptr += sizeof(NUM);
    *(NUM *)ptr = w->photosphere.plane->normal[1];   ptr += sizeof(NUM);
    *(NUM *)ptr = w->photosphere.plane->normal[2];   ptr += sizeof(NUM);
  }
  if(w->photosphere2.plane) {
    *(NUM *)ptr = w->photosphere2.plane->origin[0];   ptr += sizeof(NUM);
    *(NUM *)ptr = w->photosphere2.plane->origin[1];   ptr += sizeof(NUM);
    *(NUM *)ptr = w->photosphere2.plane->origin[2];   ptr += sizeof(NUM);
    *(NUM *)ptr = w->photosphere2.plane->normal[0];   ptr += sizeof(NUM);
    *(NUM *)ptr = w->photosphere2.plane->normal[1];   ptr += sizeof(NUM);
    *(NUM *)ptr = w->photosphere2.plane->normal[2];   ptr += sizeof(NUM);
  }

  // skip all the magic-concentration and image stuff...


  // Store the default-boundary name...
  s = boundary_ptr_to_name( w->default_bound );
  if(!s) { s = ""; }
  for(i=0;i<80;i++) {
    *ptr++ = *s;
    if(*s)
      s++;
  }

    // f_funcs as names
  for(i=0; i<N_FORCE_FUNCS; i++) {
    s = force_ptr_to_str(w->f_funcs[i]);
    if(!s) { s = ""; }
    for(j=0; j<80; j++) {
      *ptr++ = *s;
      if(*s)
	s++;
    }
  }
  
  // m_funcs as names
  for(i=0;i<N_FORCE_FUNCS;i++) {
    s = force_ptr_to_str(w->m_funcs[i]);
    if(!s) { s = ""; }
    for(j=0; j<80; j++) {
      *ptr++ = *s;
      if(*s)
	s++;
    }
  }

  // rc_funcs as names
  for(i=0; i<N_RECON_FUNCS; i++) {
    s = recon_ptr_to_str(w->rc_funcs[i]);
    if(!s) { s="";}
    for(j=0; j<80; j++) {
      *ptr++ = *s;
      if(*s)
	s++;
    }
  }
  
  // final fence
  *(long *)ptr = WORLD_END_FENCE;
  ptr += sizeof(long);

  return binary_dump_field(fd, BD_WORLD, ptr - binary_buffer, binary_buffer);
}

int binary_read_WORLD(long size, char *buf, WORLD *w) {
  char *me = "binary_read_world";
  char *ptr = buf;
  char *s;
  WORLD w0;
  long phot_ok;
  int i,j;
  
  if( *(long *)ptr != 1 ) {
    fprintf(stderr,"%s: expected WORLD dump version %d, got %ld\n",me,1,*(long *)ptr);
    return 1;
  }
  ptr += sizeof(long);


  if( *(long *)ptr != sizeof(WORLD) ) {
    fprintf(stderr,"%s: expected WORLD with a size of %ld; found %ld (oops)\n",me, sizeof(WORLD), *(long *)ptr);
    return 2;
  }
  ptr += sizeof(long);

  // Copy the WORLD structure, in its entirety
  for(i=0, s = (char *)(&w0); i<sizeof(WORLD); i++)
    *(s++) = *(ptr++);
  for(; i%16; i++) 
    ptr++;

  // Copy all the non-pointer fields into the WORLD we were given.
  w->frame_number = w0.frame_number;
  w->state        = w0.state;
  // skip refct
  // skip concentrations
  // skip lines
  // skip vertices

  // do photosphere types, but skip planes for now
  // (handle them after the main copy)
  w->photosphere.type = w0.photosphere.type;
  w->photosphere2.type = w0.photosphere2.type;

  // Skip all the magic-concentration and image stuff...
  
  w->locale_radius = w0.locale_radius;
  w->auto_open = w0.auto_open;
  
  // Skip the open & plasmoid concentrations
  // Skip the default boundary condition -- will fill in later from string.
  
  w->verbosity = w0.verbosity;
  
  // Skip force, m, and rc funcs and  rc params.
  
  // Copy m params
  for(i=0; i<N_M_PARAMS; i++)
    w->m_params[i] = w0.m_params[i];
  
  w->step_scale.b_power = w0.step_scale.b_power;
  w->step_scale.d_power = w0.step_scale.d_power;
  w->step_scale.s_power = w0.step_scale.s_power;
  w->step_scale.ds_power = w0.step_scale.ds_power;
  
  w->passno = w0.passno;
  w->handle_skew = w0.handle_skew;

  w->max_angle = w0.max_angle;
  w->mean_angle = w0.mean_angle;
  w->dtau = w0.dtau;
  w->rel_step = w0.rel_step;

  for(i=0;i<MAXNUMCOEFFS;i++) 
    w->coeffs[i] = w0.coeffs[i];
  
  w->n_coeffs = w0.n_coeffs;
  w->maxn_coeffs = w0.maxn_coeffs;
  
  w->concurrency = w0.concurrency;
  
  w->f_min = w0.f_min;
  w->f_max = w0.f_max;
  w->fr_min = w0.fr_min;
  w->fr_max = w0.fr_max;
  w->ca_min = w0.ca_min;
  w->ca_max = w0.ca_max;
  w->ca_acc = w0.ca_acc;
  w->ca_ct = w0.ca_ct;

  // RC params
  for(i=0; i<N_RECON_FUNCS;i++) {
    for(j=0; j<N_RECON_PARAMS; j++) {
      w->rc_params[i][j] = w0.rc_params[i][j];
    }
  }

  // M params
  for(i=0; i<N_M_PARAMS; i++) {
    w->m_params[i] = w0.m_params[i];
  }


  // Check fence
  if( *(long *)ptr != WORLD_VAR_FENCE ) {
    fprintf(stderr,"%s: missed variable fence! (expected 0x%lx, got 0x%lx)\n",me,WORLD_VAR_FENCE, *(long *)ptr);
    for(i=-16; i<=16; i++) {
      fprintf(stderr, "rel. pos %d: 0x%lx\n",i, *((long *)ptr + i));
    }
    return 3;
  }
  ptr += sizeof(long);
  
  // Finished copying the actual structure.  Now interpret the ancillary stuff - pointers and strings.

  if(w0.photosphere.plane) {
    if(!w->photosphere.plane) {
      w->photosphere.plane = (PLANE *)localmalloc(sizeof(PLANE),MALLOC_PLANE);
    }
    w->photosphere.plane->origin[0] = *(NUM *)ptr;  ptr += sizeof(NUM);
    w->photosphere.plane->origin[1] = *(NUM *)ptr;  ptr += sizeof(NUM);
    w->photosphere.plane->origin[2] = *(NUM *)ptr;  ptr += sizeof(NUM);
    w->photosphere.plane->normal[0] = *(NUM *)ptr;  ptr += sizeof(NUM);
    w->photosphere.plane->normal[1] = *(NUM *)ptr;  ptr += sizeof(NUM);
    w->photosphere.plane->normal[2] = *(NUM *)ptr;  ptr += sizeof(NUM);
  } else {
    if(w->photosphere.plane)
      localfree(w->photosphere.plane);
    w->photosphere.plane = 0;
  }

  if(w0.photosphere2.plane) {
    if(!w->photosphere2.plane) {
      w->photosphere2.plane = (PLANE *)localmalloc(sizeof(PLANE),MALLOC_PLANE);
    }
    w->photosphere2.plane->origin[0] = *(NUM *)ptr;  ptr += sizeof(NUM);
    w->photosphere2.plane->origin[1] = *(NUM *)ptr;  ptr += sizeof(NUM);
    w->photosphere2.plane->origin[2] = *(NUM *)ptr;  ptr += sizeof(NUM);
    w->photosphere2.plane->normal[0] = *(NUM *)ptr;  ptr += sizeof(NUM);
    w->photosphere2.plane->normal[1] = *(NUM *)ptr;  ptr += sizeof(NUM);
    w->photosphere2.plane->normal[2] = *(NUM *)ptr;  ptr += sizeof(NUM);
  } else {
    if(w->photosphere2.plane)
      localfree(w->photosphere2.plane);
    w->photosphere2.plane = 0;
  }

  
  // Start interpreting strings...

  // Default boundary condition
  w->default_bound = boundary_name_to_ptr(ptr);
  ptr += 80;

  // Force funcs
  for(i=0; i<N_FORCE_FUNCS; i++)  {
    w->f_funcs[i] = force_str_to_ptr(ptr);
    ptr += 80;
  }

  // M funcs
  for(i=0; i<N_FORCE_FUNCS; i++) {
    w->m_funcs[i] = force_str_to_ptr(ptr);
    ptr += 80;
  }
  
  // RC funcs
  for(i=0; i<N_RECON_FUNCS; i++) {
    w->rc_funcs[i] = recon_str_to_ptr(ptr);
    ptr += 80;
  }
  
  // final fence
  if( *(long *)ptr != WORLD_END_FENCE ) {
    fprintf(stderr,"%s: missed end-of-WORLD fence! (expected 0x%lx, got 0x%lx)\n",me, WORLD_END_FENCE, *(long *)ptr);
    for(i=-80; i<=80; i++) {
      fprintf(stderr, "rel. pos %d: 0x%lx\n",i, *((long *)ptr + i));
    }
    return 4;
  }

  return 0;
}

/******************************
 * binary_dump_CONCENTRATION
 * 
 * Dumps a FLUX_CONCENTRATION, minus fluxons and vertices.
 */
int binary_dump_CONCENTRATION(int fd, FLUX_CONCENTRATION *fc) {
  char *ptr;
  char *s;
  int i,j;
  
  check_binary_buf( sizeof(FLUX_CONCENTRATION) +
		    sizeof(long)*8 +
		    80 * 2
		    );
  ptr = binary_buffer;

  *(long *)ptr = 1; // CONCENTRATION dump version number
  ptr += sizeof(long);

  *(long *)ptr = sizeof(FLUX_CONCENTRATION); 
  ptr += sizeof(long);

  // Copy the FLUX_CONCENTRATION structure, in its entirety
  for(i=0, s=(char *)fc; i<sizeof(FLUX_CONCENTRATION); i++)
    *(ptr++) = *(s++);
  for(; i%16; i++)
    *(ptr++) = 0;

  // Write a fence (re-use the WORLD fence)
  *(long *)ptr = WORLD_VAR_FENCE;
  ptr += sizeof(long);


  // Now patch up the concentration...

  // skip world
  // skip lines tree
  // skip concentrations links
  
  // store the boundary condition pointer as a string
  s = boundary_ptr_to_name( fc->bound );
  if(!s) { s=""; }
  for(i=0;i<80;i++) {
    *ptr++ = *s;
    if(*s)
      s++;
  }

  // all done - write another fence (re-use the WORLD fence);
  *(long *)ptr = WORLD_END_FENCE;
  ptr += sizeof(long);
  
  return binary_dump_field(fd, BD_CONCENTRATION, ptr - binary_buffer, binary_buffer);
}
 
int binary_read_CONCENTRATION(long size, char *buf, WORLD *w) {
  char *me = "binary_read_concentration";
  char *ptr = buf;
  char *s;
  FLUX_CONCENTRATION fc;
  FLUX_CONCENTRATION *f;
  int i;
  char new_conc;

  if( *(long *)ptr != 1 ) { // Check for version
    fprintf(stderr,"%s: version number is wrong (expected %d, got %ld)\n",me, 1, *(long *)ptr);
    return 1;
  }
  ptr += sizeof(long);

  if( *(long *)ptr != sizeof(FLUX_CONCENTRATION)) {
    fprintf(stderr,"%s: size of FLUX_CONCENTRATION is wrong (expected %ld, got %ld)\n",me, sizeof(FLUX_CONCENTRATION), *(long *)ptr);
    return 2;
  }
  ptr += sizeof(long);

  for(i=0, s = (char *)&fc; i<sizeof(FLUX_CONCENTRATION); i++)
    *(s++) = *(ptr++);
  for(; i%16; i++)
    *(ptr++) = 0;

  if(*(long *)ptr != WORLD_VAR_FENCE) {
    fprintf(stderr,"%s: Missed variable field fence! (expected 0x%lx, got 0x%lx)\n",me,WORLD_VAR_FENCE, *(long *)ptr);
    for(i=-10; i<=10; i++) {
      fprintf(stderr, "rel. pos %d: 0x%lx\n",i, *((long *)ptr + i));
    }
    return 3;
  }
  ptr += sizeof(long);

  fc.bound = boundary_name_to_ptr(ptr);
  ptr += 80;
  
  if(*(long *)ptr != WORLD_END_FENCE) {
    fprintf(stderr,"%s: Missed end-of-concentration fence! (expected 0x%lx, got 0x%lx)\n",me,WORLD_VAR_FENCE, *(long *)ptr);
    for(i=-10; i<=10; i++) {
      fprintf(stderr, "rel. pos %d: 0x%lx\n",i, *((long *)ptr + i));
    }
    return 4;
  }

  // Now copy the fc into a new (or old!) flux concentration in the WORLD.
  f = tree_find(w->concentrations, fc.label, fc_lab_of, fc_ln_of);
  new_conc = ( !f );
  if(new_conc) {
    f = new_flux_concentration(w, fc.x[0], fc.x[1], fc.x[2], fc.flux, fc.label);
    check_special_concentration(f);
  }
  
  f->world = w;
  f->flux = fc.flux;
  f->label = fc.label;
  f->x[0] = fc.x[0];
  f->x[1] = fc.x[1];
  f->x[2] = fc.x[2];
  f->locale_radius = fc.locale_radius;
  f->passno = fc.passno;
  f->bound = fc.bound;

  return 0;
}
    
/******************************
 * binary_dump_FLUXON
 * 
 * Dumps a FLUXON (including concentration label), including vertices but not neighbors.
 */

int binary_dump_FLUXON(int fd, FLUXON *f) {
  char *me = "binary_dump_FLUXON";
  char *ptr;
  char *s;
  int i,j;
  VERTEX *v;

  check_binary_buf( sizeof(FLUXON) + 
		    sizeof(long)*8 + 
		    (sizeof(VERTEX)+16) * f->v_ct
		    );
  
  ptr = binary_buffer;
  
  *(long *)ptr = 1; // FLUXON dump version number
  ptr += sizeof(long);

  *(long *)ptr = sizeof(FLUXON);  // size of a fluxon
  ptr += sizeof(long);

  *(long *)ptr = sizeof(VERTEX);  // size of a vertex
  ptr += sizeof(long);

  *(long *)ptr = f->fc0 ? f->fc0->label : 0;  // start FC label
  ptr += sizeof(long);
		   
  *(long *)ptr = f->fc1 ? f->fc1->label : 0;  // end FC label
  ptr += sizeof(long);

  // Dump the FLUXON structure in its entirety.
  for( i=0, s=(char *)f; i<sizeof(FLUXON); i++) 
    *(ptr++) = *(s++);
  for(; i%16; i++)
    *(ptr++) = 0;

  // leave a fence (re-use the WORLD fence)
  *(long *)ptr = WORLD_VAR_FENCE;
  ptr += sizeof(long);

  // dump the vertices as a unit...

  for(i=0, v=f->start; i<f->v_ct && v; i++,v=v->next) {
    s = (char *)v;
    for(j=0;j<sizeof(VERTEX); j++)
      *(ptr++) = *(s++);
    for(; j%16; j++)
      *(ptr++) = 0;
  }
  if(i != f->v_ct || v) {
    fprintf(stderr,"%s: inconsistent vertex count in fluxon %ld; giving up (i=%d)\n",me, f->label,i);
    return 3;
  }

  // leave another fence 
  *(long *)ptr = WORLD_END_FENCE;
  ptr += sizeof(long);

  return binary_dump_field(fd, BD_FLUXON, ptr - binary_buffer, binary_buffer);
}

int binary_read_FLUXON(long size, char *buf, WORLD *w) {
  char *me = "binary_read_FLUXON";
  char *ptr= buf;
  char *s;
  FLUXON f0;
  FLUXON *f;
  long fc0lab, fc1lab, lab;
  FLUX_CONCENTRATION *fc0, *fc1;
  VERTEX v0, *v, *vp;
  int i,j;
  char fluxon_was_new;

  // Check for version
  if( *(long *)ptr != 1) {
    fprintf(stderr,"%s: version number is wrong (expected %d, got %ld)\n",me,1,*(long *)ptr);
    return 1;
  }
  ptr += sizeof(long);

  // Check for size of a fluxon
  if( *(long *)ptr != sizeof(FLUXON)) {
    fprintf(stderr,"%s: FLUXON structure has wrong size (expected %ld, got %ld)\n",me, sizeof(FLUXON), *(long *)ptr);
    return 2;
  }
  ptr += sizeof(long);

  // Check for size of a vertex
  if( *(long *)ptr != sizeof(VERTEX)) {
    fprintf(stderr,"%s: VERTEX structure has wrong size (expected %ld, got %ld)\n",me, sizeof(VERTEX), *(long *)ptr);
    return 3;
  }
  ptr += sizeof(long);

  // start concentration
  fc0lab = *(long *)ptr;
  fc0 = tree_find(w->concentrations, fc0lab, fc_lab_of, fc_ln_of);
  if(!fc0) {
    fprintf(stderr,"%s: couldn't find start concentration %ld in existing WORLD!\n",me, fc0lab);
    return 4;
  }
  ptr += sizeof(long);

  // end concentration
  fc1lab = *(long *)ptr;
  fc1 = tree_find(w->concentrations, fc1lab, fc_lab_of, fc_ln_of);
  if(!fc1) {
    fprintf(stderr,"%s: couldn't find end concentration %ld in existing WORLD!\n",me, fc1lab);
    return 5;
  }
  ptr += sizeof(long);

  // Read the FLUXON into buffer...
  for(i=0, s=(char *)&f0; i<sizeof(FLUXON); i++) 
    *(s++) = *(ptr++);
  for(; i%16; i++)
    ptr++;

  lab = f0.label;
  f = tree_find(w->lines, lab, fl_lab_of, fl_all_ln_of);
  if(!f) {
    f = new_fluxon( f0.flux, fc0, fc1, lab, f0.plasmoid );
    fluxon_was_new = 1;
  }

  // Copy the data into the output FLUXON.
  f->flux = f0.flux;
  f->label = f0.label;
  // f->start = 0; // redundant -- 0 on creation or after our purge above
  // f->end = 0;   // redundant -- 0 on creation or after our purge above
  f->v_ct = f0.v_ct;

  // Handle linking in the old-fluxon case
  if( !fluxon_was_new ) {
    if( f->fc0 != fc0 ) {
      f->fc0->lines = tree_unlink( f, fl_lab_of, fl_start_ln_of );
      fc0->lines = tree_binsert( fc0->lines, f, fl_lab_of, fl_start_ln_of );
    }
    if( f->fc1 != fc1 ) {
      f->fc1->lines = tree_unlink( f, fl_lab_of, fl_end_ln_of );
      fc1->lines = tree_binsert( fc1->lines, f, fl_lab_of, fl_end_ln_of );
    }
    f->fc0 = fc0;
    f->fc1 = fc1;
  }

  // Get rid of old vertices if this isn't a new fluxon...
  if(!fluxon_was_new) {
    for(v=f->start; v; v=f->start )
      delete_vertex(v);
  }
  
  // Check the fence
  if( *(long *)ptr != WORLD_VAR_FENCE ) {
    fprintf(stderr,"%s: missed variable fence!\n",me);
    if(fluxon_was_new) 
      delete_fluxon(f);
    return 6;
  }
  ptr += sizeof(long);
      
  // Grab vertices...
  vp = 0;

  for(i=0; i<f0.v_ct; i++) {
    void *l,*r,*up;
    long n;
    NUM sum;

    // Snarf the new vertex into the vertex buffer
    for(j=0, s=(char *)&v0; j<sizeof(VERTEX); j++) 
      *(s++) = *(ptr++);
    for(; j%16; j++)
      ptr++;

    //    printf("f label is %d, f->fc0 label is %d, world is %d; ",f->label, f->fc0->label, f->fc0->world);
    v = new_vertex( v0.label, v0.x[0], v0.x[1], v0.x[2], f );

    // Copy the vertex buffer fields into the new vertex.
    v->line = f;
    v->prev = vp;
    if(vp)
      vp->next = v;
    else
      f->start = v;
    v->next = 0; // will be overwritten on next i iteration.

    // Skip x as new_vertex handles that.

    v->neighbors.n = 0;
    v->neighbors.size = 0;
    v->neighbors.stuff = 0;
    v->nearby.n = 0;
    v->nearby.size = 0;
    v->nearby.stuff = 0;

    cp_3d(v->scr, v0.scr);
    
    v->passno = v0.passno;
    
    v->r = v0.r;
    v->a = v0.a;

    cp_3d(v->b_vec, v0.b_vec);

    v->b_mag = v0.b_mag;
    v->energy = v0.energy;

    cp_3d(v->f_v, v0.f_v);
    cp_3d(v->f_s, v0.f_s);
    cp_3d(v->f_t, v0.f_t);
    v->f_s_tot = v0.f_s_tot;
    v->f_v_tot = v0.f_v_tot;
    v->f_n_tot = v0.f_n_tot;
    v->r_v    = v0.r_v;
    v->r_s    = v0.r_s;
    v->r_cl   = v0.r_cl;
    v->r_ncl  = v0.r_ncl;
    v->label  = v0.label;
    cp_3d(v->plan_step, v0.plan_step);

    vp = v;
  }

  f->end = v;

  if( *(long *)ptr  != WORLD_END_FENCE ) {
    fprintf(stderr,"%s: missed final fence in fluxon %ld, arena may be corrupted.\n",me, f->label);
    for(i=-10; i<=10; i++) {
      fprintf(stderr, "rel. pos %d: 0x%lx\n",i, *((long *)ptr + i));
    }

    return 7; 
  }
 
  return 0;
}
		    
/******************************
 * binary_dump_neighbors - write the neighbor list of all vertices in a fluxon.
 */
  
int binary_dump_neighbors(int fd, FLUXON *f) {
  char *me = "binary_dump_neighbors";
  char *ptr;
  char *s;
  int i,j;
  VERTEX *v;

  // Accumulate total number of neighbors on this fluxon
  for(v=f->start, j=i=0; v; v=v->next, j++) {
    i += v->neighbors.n;
  }

  if(j != f->v_ct)  {
    fprintf(stderr,"%s: inconsistent vertex count in fluxn %ld (v_ct is %ld, counted %d neighbors); giving up\n", me, f->label, f->v_ct, j);
    return 1;
  }
  
  check_binary_buf( sizeof(long) * i  +
		    sizeof(long) * 5 +
		    (3 * sizeof(long)) * f->v_ct
		    );
  ptr = binary_buffer;

  *(long *)ptr = 1;  // Version number
  ptr += sizeof(long);

  *(long *)ptr = f->label; // Which fluxon's vertices to dump
  ptr += sizeof(long);
  
  *(long *)ptr = f->v_ct; // Total number of vertices to dump
  ptr += sizeof(long);

  // Loop over vertices within the fluxon...
  for(v=f->start; v; v=v->next) {
    *(long *)ptr = NEIGHBOR_FENCE;  // Fence for safety
    ptr += sizeof(long);

    *(long *)ptr = v->label;        // Label of this vertex
    ptr += sizeof(long);
    
    *(long *)ptr = v->neighbors.n;  // Neighbor count
    ptr += sizeof(long);
    
    for(j=0;j<v->neighbors.n;j++)  { // The neighbors (as labels)
      *(long *)ptr = ((VERTEX **)(v->neighbors.stuff))[j]->label;
      ptr += sizeof(long);
    }
  }
  
  *(long *)ptr = WORLD_END_FENCE;
  ptr+= sizeof(long);

  return binary_dump_field(fd, BD_NEIGHBORS, ptr - binary_buffer, binary_buffer);
}

int binary_read_neighbors(long size, char *buf, WORLD *w) {
  char *me = "binary_read_neighbors";
  char *ptr = buf;
  char *s;
  long i,j, ct;
  FLUXON *f;
  VERTEX *v;
 
  if( *(long *)ptr != 1 ) { // Check version number
    fprintf(stderr,"%s: wrong version (expected %d, got %ld)\n",me, 1, *(long *)ptr);
    return 1;
  }
  ptr += sizeof(long);

  // Read fluxon label, and find it in the tree...
  f = tree_find( w->lines, *(long *)ptr, fl_lab_of, fl_all_ln_of );
  if(!f) {
    fprintf(stderr, "%s: couldn't find fluxon %ld! I give up.\n",me, *(long *)ptr);
    return 2;
  }
  ptr += sizeof(long);

  // Check number of vertices
  if( (*(long *)ptr) != f->v_ct ) {
    fprintf(stderr, "%s: fluxon %ld thinks it has %ld vertices, not %ld as advertised in the file...\n",me, f->label, f->v_ct, *(long *)ptr);
    
    {
      VERTEX *vv;
      int i;
      for(vv=f->start, i=0; vv; vv=vv->next, i++)
	;
      fprintf(stderr,"    ... and it seems to actually have %d vertices.\n",i);
      if(i == *(long *)ptr) {
	f->v_ct = i;
	fprintf(stderr,"  Fixed it.\n");
      } else {
	fprintf(stderr,"  I give up.\n");
	return 3;
      }
    }
  }
  ptr += sizeof(long);

  for( v=f->start, i=0; v; v=v->next, i++ ) {
    if( *(long *)ptr != NEIGHBOR_FENCE ) {
      fprintf(stderr,"%s: missed neighbor fence on vertex #%ld of fluxon %ld; giving up\n",me, i, f->label);
      return 4;
    }
    ptr += sizeof(long);

    if( *(long *)ptr != v->label ) {
      fprintf(stderr,"%s: file has vertex %ld at position #%ld of fluxon %ld, we have %ld instead.  Giving up.\n",me, *(long *)ptr, i, f->label, v->label);
      return 5;
    }
    ptr += sizeof(long);
    
    // Get count.
    ct = *(long *)ptr;
    ptr += sizeof(long);

    // Clear the dumblist, then add the neighbors from the file...
    vertex_clear_neighbors(v);
    for(j=0; j<ct; j++) {
      VERTEX *vn = tree_find(w->vertices, *(long *)ptr, v_lab_of, v_ln_of );
      if(!vn) {
	fprintf(stderr,"%s: vertex %ld is called out as a neighbor to %ld (fluxon %ld, pos. #%ld) but doesn't exist!\n",me, *(long *)ptr, v->label, f->label, j);
	return 6;
      }
      vertex_add_neighbor( v, vn);
      ptr += sizeof(long);
    }
  }

  
  if(*(long *)ptr != WORLD_END_FENCE) {
    fprintf(stderr, "%s: missed end fence\n",me);
    return 7;
  }
  
  return 0;
}
    
  
/******************************
 * binary_dump_end - send an end marker to the file.  Does not close the file descriptor.
 */
int binary_dump_end(int fd) {
  binary_dump_field( fd, BD_END, 0, 0 );
  close(fd);
  return 1;
}


/******************************
 * binary_dump_fluxon_pipe - send a fluxon's complete data out to the file, including 
 * neighbor/nearby information as pointers (so it is not suitable for different instantiations
 * of the model -- only for communicating state from a daughter process back to the parent).
 *
 * Format: a long with version number, a long indicating the fluxon ID, a long indicating 
 * v_ct, and then the VERTEX_FORCE_DUMPs in order.
 *
 *
 * Because the vertex data is variable length (we include the dumblists) we use a 
 * fence, for belt-and-suspenders error checking goodness.
 * 
 * WARNING:
 * dump_fluxon_forces is intended for communicating back from a daughter process for
 * first-stage parallelization.  Hence we actually send pointers (whose value should
 * be preserved across forking), rather than labels.  
 * 
 * Format: 
 *  [long]               fence (version number) 
 *  [long]               fluxon label
 *  [FLUXON *]           fluxon pointer
 *  [long]               number of vertices
 *
 *  [VERTEX_FORCE_DUMP]  force dump structure
 *  [void * x n]         neighbor pointers
 *  [long]               trailing fence
 * 
 *  - fixed-length block (the VERTEX_FORCE_DUMP structure)
 *
 */
#define MAX_AV_NEIGHBORS 100
int binary_dump_fluxon_pipe( int fd, FLUXON *f) {
  int v_ct = f->v_ct;
  int i,j;
  long v_ct_max;
  VERTEX_PIPE_FIXED *vd;
  VERTEX *v;
  char *dex;
  long neighbors_found = 0;
  long neighbors_allowed_for = MAX_AV_NEIGHBORS * f->v_ct;
  long len = 
    3 * sizeof(long)  + sizeof( FLUXON * ) +   // Header length
    f->v_ct * ( sizeof(long) +                 // Fence
		sizeof(VERTEX_PIPE_FIXED) +    // fixed-length structure
		MAX_AV_NEIGHBORS * sizeof(VERTEX *) +       // Allow 100 neighbors and 100 nearby - way overkill
		sizeof(long)                   // Trailing fence
		);

  // Make sure we've got a large enough buffer
  check_binary_buf( len );
  dex = binary_buffer;
  
  *(long *)dex = 3;         // Version number of the fluxon force dump;
  dex += sizeof(long);
  
  *(long *)dex = f->label;  // Label of the fluxon
  dex += sizeof(long);

  *(FLUXON **)dex = f;      // fluxon pointer itself
  dex += sizeof(FLUXON *);
  
  *(long *)dex = f->v_ct;   // Number of vertices to expect
  dex += sizeof(long);
  
  /*** Now copy vertex info into the buffer... ***/
  v_ct_max = 0;
  for(v=f->start, i=0; v; v=v->next, i++) {

    if(v->neighbors.n > v_ct_max) 
      v_ct_max = v->neighbors.n;

    if(v->neighbors.n > MAX_AV_NEIGHBORS) {
      fprintf(stderr,"binary_dump_fluxon_pipe: WARNING - vertex %ld is wonky (%d neighbors!)",v->label, v->neighbors.n);
    }

    // i is just for error checking - complain here
    if( i >= f->v_ct ) {
      fprintf(stderr,"binary_dump_fluxon_pipe: Wonky fluxon %ld - has more vertices than v_ct (%ld)! I quit\n",f->label, f->v_ct);
      return 1;
    }

    // Now copy the fixed part into the buffer...
    vd = (VERTEX_PIPE_FIXED *)dex;
    dex += sizeof(VERTEX_PIPE_FIXED);

    vd->label  = v->label;
    cp_3d( vd->x,         v->x);
    cp_3d( vd->scr,       v->scr);
    cp_3d( vd->b_vec,     v->b_vec);
    cp_3d( vd->f_v,       v->f_v);
    cp_3d( vd->f_s,       v->f_s);
    cp_3d( vd->f_t,       v->f_t);
    cp_3d( vd->plan_step, v->plan_step);
    vd->r       = v->r;
    vd->a       = v->a;
    vd->b_mag   = v->b_mag;
    vd->energy  = v->energy;
    vd->f_s_tot = v->f_s_tot;
    vd->f_v_tot = v->f_v_tot;
    vd->r_v     = v->r_v;
    vd->r_s     = v->r_s;
    vd->r_cl    = v->r_cl;
    vd->r_ncl   = v->r_ncl;
    vd->neighbors_n = v->neighbors.n;
    vd->f_n_tot = v->f_n_tot;
    vd->A = v->A;

    neighbors_found += v->neighbors.n;
    if(neighbors_found > neighbors_allowed_for){
      fprintf(stderr,"binary_dump_fluxon_pipe: Whoa!  Never expected to see this.  Found %ld neighbors along fluxon %ld (which has %ld vertices)- that's over %d per vertex! I give up.\n",neighbors_found, f->label, f->v_ct, MAX_AV_NEIGHBORS);
      fflush(stderr);
      return 2;
    }

    // Now copy the neighbors...
    for(j=0; j<v->neighbors.n; j++) {
      *(void **)dex = v->neighbors.stuff[j];
      dex += sizeof(void *);
    }
    
    // Finally, add a fence.
    *(long *)dex = VERTEX_END_FENCE;
    dex += sizeof(long);
  }
  
  /** Dump the buffer to the file **/
  binary_dump_field( fd, BD_FLUXON_PIPE, dex - binary_buffer, binary_buffer );
  return 1;
}

/**********************************************************************
 * binary_read_dumpfile - read in a complete binary dump, returning the
 * WORLD object you've been working on (or 0 on failure).
 *
 */
WORLD *binary_read_dumpfile ( int fd, WORLD *w ) {
  long hdrbuf[3];
  int ct;
  int pos = 0;
  long type;
  long len;
  char *me = "FLUX: binary_read_dumpfile";
  char allocated_world = 0;

  if(!w) {
    w = new_world();
    allocated_world = 1;
  }

  do {

    // Read the header (try twice if necessary)
    ct = read( fd, (char *)hdrbuf, sizeof(long)*3 );
    if(ct != sizeof(long)*3) {
      int ct2 = -10; // any ol' large negative number
      
      if(ct>=0) 
	ct2 = read(fd, (char *)hdrbuf + ct, sizeof(long)*3 - ct);
      
      if(ct2+ct != sizeof(long)*3) {
	fprintf(stderr,"%s: failed to read field header, position %d",me,pos);
	perror("read returned the error");
	if(allocated_world)
	  free_world(w);
	return 0;
      }
    }

    // Parse header: check for fence, and get data type and length
    pos += ct;
    if(hdrbuf[0] != BD_FENCE) {
      fprintf(stderr,"%s: failed to find fence (expected %lx, got %lx), position %d", me, BD_FENCE, hdrbuf[0],pos);
      if(allocated_world)
	free_world(w);
      return 0;
    }

    type = hdrbuf[1];
    len = hdrbuf[2];
    if(type <=0 || type > BD_MAX_TYPENO) {
      fprintf(stderr,"%s: failed to find valid type (expected 1-%d, got %ld), position %d",me,BD_MAX_TYPENO,type,pos);
      if(allocated_world)
	free_world(w);
      return 0;
    }
    check_binary_buf( len );

    // Read the packet (try twice if necessary)
    ct = read( fd, binary_buffer, len );
    if(ct != len) {
      int ct2 = -10; // any ol' large negative number

      if(ct>=0)
	ct2= read(fd, binary_buffer+ct, len-ct);
      
      if(ct2+ct != len)  {
	fprintf(stderr,"%s: failed to read %ld bytes (got %d bytes from position %d, type %ld; second try yielded %d bytes for %d total)\n",me,len,ct, pos, type,ct2,ct+ct2);
	perror("read returned the error");
	if(allocated_world)
	  free_world(w);
	return 0;
      }
    }

    pos += len;

    if(w->verbosity>1) {
      printf("read_dumpfile: found fence; type=%ld, len=%ld\n",hdrbuf[1],hdrbuf[2]);
    }

    ///// Dispatch read to the correct reader...
    {
      int (*reader)(long size, char *buf, WORLD *w) = 0;

      switch(type) {
      case BD_END:
	close(fd);
	return w;
	break;
      case BD_FLUXON_PIPE:
	reader = binary_read_fluxon_pipe;
	break;
      case BD_WORLD:
	reader = binary_read_WORLD;
	break;
      case BD_CONCENTRATION:
	reader = binary_read_CONCENTRATION;
	break;
      case BD_FLUXON:
	reader = binary_read_FLUXON;
	break;
      case BD_NEIGHBORS:
	reader = binary_read_neighbors;
	break;
      case BD_POSITION:
	reader = binary_read_flpos;
	break;
      case BD_STEP:
	reader = binary_read_flstep;
	break;
      default:
	fprintf(stderr,"WARNING: typecode %ld not implemented - you should never see this message from binary_read_dumpfile, in io.c...\n",type);
	return 0;
	break;
      }
	
      if(reader) {
	if( (*reader)(len, binary_buffer, w) ) {
	  fprintf(stderr,"%s",buf);
	  if(allocated_world)
	    free_world(w);
	  return 0;
	}
      } else {
	fprintf(stderr,"binary_read_dumpfile: you should not ever see this message...\n");
      }

    }// end of dispatch convenience block

  } while(1);
}    

/**********************************************************************
 * read routines - these are dispatched from the binary dump reader.
 * Each gets a buffer with the binary data from the packet.
 * 
 * Each read routine returns 0 on success, nonzero on failure.
 * On failure, the I/O buffer has an error message in it.
 */


/******************************
 * binary_read_fluxon_pipe
 * 
 * Reads a fluxon field produced by binary_dump_fluxon_pipe.
 * 
 * WARNING: the fluxon must have been dumped by this process or a daughter, 
 * as actual pointers are used rather than labels. 
 *
 * Version number must match...
 *
 * Also, we accumulate min/max force and angle info per vertex -
 * this uses the vertex_accumulate_f_minmax routine in model.c.
 *
 */
int binary_read_fluxon_pipe( long size, char *buf, WORLD *w ) {
  long f_lab;
  long v_ct;
  FLUXON *f;
  VERTEX_PIPE_FIXED *vd;
  VERTEX *v;
  int i,j,k;
  char *me = "FLUX I/O library: binary_read_fluxon_forces";
  char *dex = buf;
  static DUMBLIST *dl = 0;
  
  if(dl==0) {
    dl = new_dumblist();
  }

  if(w->verbosity>1) {
    printf("read_fluxon_pipe...");
    fflush(stdout);
  }



  if(size< 3*sizeof(long)) {
    fprintf(stderr, "%s: inconceivable packet size %ld is too small!\n",me,size);
    return 1;
  }

  if( (*(long *)dex) != 3) { // Check version number
    fprintf(stderr, "%s: packet is the wrong version (%ld, I am version %d)\n",me,*(long *)dex,3);
    return 2;
  }
  dex+= sizeof(long);

  f_lab = *(long *)dex; // Label of the fluxon
  dex += sizeof(long);

  f = *(FLUXON **)dex;  // Fluxon pointer itself
  dex += sizeof(FLUXON *);
  
  if(f->label != f_lab) {
    fprintf(stderr,"%s: found bogus fluxon pointer in dump! (expected fluxon %ld, got %ld)",me,f_lab, f->label);
    return 3;
  }

  v_ct = *(long *)dex; // Vertex count
  dex += sizeof(long);
  
  if(f->v_ct != v_ct) {
    fprintf(stderr,"%s: vertex count doesn't agree in fluxon %ld (expected %ld from file, found %ld in fluxon)\n",me,f_lab,v_ct, f->v_ct);
    return 4;
  }

  if(w->verbosity>1) {
    printf(" r: f=%ld,vct=%ld    ",f->label,f->v_ct);
    fflush(stdout);
  }

  for(v=f->start, i=0; i<v_ct; v=v->next, i++) {
    vd = (VERTEX_PIPE_FIXED *)dex;
    dex += sizeof(VERTEX_PIPE_FIXED);

    if(!vd || !v || vd->label != v->label) {
      fprintf(stderr,"%s: vertex label mismatch in fluxon %ld, pos %d - file expected %ld, found %ld\n",me,f_lab,i,vd?vd->label:-1,v?v->label:-1);
      return 5;
    }

    cp_3d( v->x,          vd->x );
    cp_3d( v->scr,        vd->scr );
    cp_3d( v->b_vec,      vd->b_vec );
    cp_3d( v->f_v,        vd->f_v );
    cp_3d( v->f_s,        vd->f_s );
    cp_3d( v->f_t,        vd->f_t );
    cp_3d( v->plan_step,  vd->plan_step );
    v->r        = vd->r;
    v->a        = vd->a;
    v->b_mag    = vd->b_mag;
    v->energy   = vd->energy;
    v->f_s_tot  = vd->f_s_tot;
    v->f_v_tot  = vd->f_v_tot;
    v->r_v      = vd->r_v;
    v->r_s      = vd->r_s;
    v->r_cl     = vd->r_cl;
    v->r_ncl    = vd->r_ncl;
    v->f_n_tot  = vd->f_n_tot;
    v->A        = vd->A;
    
    // Note: passno not copied!  (the daughter may have incremented it a lot; we want to make sure
    // nothing has a passno higher than the current world passno....)

    if(w->verbosity>1) {
      printf("vertex %ld: found %ld neighbors...\n",v->label,vd->neighbors_n);
    }

    /* Accumulate force and angle.  We lag by one vertex, to get the curvature right. */
    if(v->prev && v->prev->prev) {
      vertex_accumulate_f_minmax(v->prev,w); // in model.c
    }

    /////Handle neighbors...
    // First - collect the neighbors in a local dumblist.
    dl->n = 0;
    for(j=0; j<vd->neighbors_n; j++) {
      dumblist_add( dl, *(void **)dex );
      dex += sizeof(void *);
    }

    // Next - check the fence!
    if( *(long *)dex != VERTEX_END_FENCE ) {
      fprintf(stderr,"%s: missed end fence while copying vertex %d (label %ld) into fluxon %ld - giving up!\n",me, i, v->label, f->label);
      return 6;
    }
    dex += sizeof(long);


#if DEBUG_DUPE_CHECK    
    /* Delete this check in production code */
    {
      int ii,jj;
      for(ii=1; ii<dl->n;ii++) {
	for(jj=0;jj<ii;jj++) {
	  if(dl->stuff[ii]==dl->stuff[jj]) {
	    printf("YOIKS! daughter gave us neighbor dupes ( vertex %d, pos %d and %d)\n",((VERTEX **)(dl->stuff))[ii]->label,ii,jj);
	  }
	}
      }

      for(ii=1;ii<v->neighbors.n; ii++) {
	for(jj=0;jj<ii;jj++) {
	  if(v->neighbors.stuff[ii]==v->neighbors.stuff[jj]) {
	    printf("YOW! parent had neighbor dupes (vertex %d, pos %d and %d)\n",((VERTEX **)(v->neighbors.stuff))[ii]->label,ii,jj);
	  }
	}
      }
    }
#endif

    ////// Now merge the neighbors with the local neighbor list. 


    // Step 1: delete anything that is in the current list but not in the daughter-delivered neighbor list.
    {
      long passno = ++(w->passno); 

      for(j=0; j<v->neighbors.n; j++) {

	// Check for dupes -- if passno matches our current passno, then we've been here before.
	if( ((VERTEX *)(v->neighbors.stuff[j]))->passno == passno) {

	  // Delete it!  (Don't delete "nearby" back-pointer entries -- that was handled the 
	  // first time we touched this VERTEX.)
	  dumblist_rm( &(v->neighbors), j );
	  // Decrement counter - last element moved to this slot, so we have to re-do the slot.
	  j--;

	} else { 
	  int keep;
	  // not a dupe -- check to make sure this neigbor is in the daughter-delivered list.

	  for(k=keep=0; k<dl->n && !keep; k++) 
	    keep = ( dl->stuff[k] == v->neighbors.stuff[j] );
	  
	  // Mark that we've been here...
	  ((VERTEX *)(v->neighbors.stuff[j]))->passno = passno; 

	  if(!keep) {
	    // Delete the back link first, then remove the item from the neighbor list.
	    dumblist_delete(  &( ((VERTEX *)(v->neighbors.stuff[j]))->nearby ), v );
	    dumblist_rm( &(v->neighbors), j );
	    // Decrement counter - last element moved to this slot, so we have to re-do the slot.
	    j--;
	  } // end of don't-keep contingency
	} // end of non-dupe contingency
      } // end of j-loop over the neighbors
    } // end of convenience block for merging

    // Step 2:  add any new neighbors reported by the daughter, who are not in the current list.
    for(k=0; k<dl->n; k++) {
      int found;
      for(j=found=0; j<v->neighbors.n && !found; j++) 
	found = ( dl->stuff[k] == v->neighbors.stuff[j] );

      if(!found) {
	dumblist_add( &(v->neighbors), dl->stuff[k] );
	dumblist_add( &( ((VERTEX *)(dl->stuff[k]))->nearby ), v );
      }
    }

    if(w->verbosity>1) {
      printf("after transfer: %d neighbors.      ",v->neighbors.n);
    }

#if DEBUG_DUPE_CHECK
    /* Delete this check in production code */
    {
      int ii,jj;
      for(ii=1;ii<v->neighbors.n;ii++) {
	for(jj=0;jj<ii;jj++) {
	  if(v->neighbors.stuff[ii]==v->neighbors.stuff[jj]) {
	    printf("YAWN: after merge, we had neighbor dupes (vertex %d, pos %d and %d)\n",((VERTEX **)(v->neighbors.stuff))[ii]->label,ii,jj);
	  }
	}
      }
    }
#endif

  }

  return 0;
}

    
  
/******************************
 * binary_dump_flpos
 * 
 * Dumps just the X vectors of the associated vertices.
 * Used for parallelization.
 */

int binary_dump_flpos( int fd, FLUXON *f) {
  long len = 
    sizeof(NUM) * 3 * f->v_ct              // vectors
    + 3*sizeof(long)                       // version, fluxon label, vertex count
    + sizeof(FLUXON *)                     // fluxon pointer
    + sizeof(long);                        // fence
  char *dex;
  VERTEX *v;
  WORLD *w = f->fc0->world;
  long verbosity = w->verbosity;
  long i;

  check_binary_buf(len);
  dex = binary_buffer;
  *(long *)dex = 1; // Version number
  dex += sizeof(long);

  *(long *)dex = f->label; // Label of the fluxon
  dex += sizeof(long) ;

  *(FLUXON **)dex = f; // fluxon pointer itself
  dex += sizeof(FLUXON *);

  *(long *)dex = f->v_ct; // Number of vertices to expect
  dex += sizeof(long);

  *(long *)dex = 0xAABBCCDD; // lay down a fence
  dex += sizeof(long);

  for(v=f->start, i=0; v; v=v->next,i++) {
    if(i >= f->v_ct) {
      fprintf(stderr,"binary_dump_positions: fluxon %ld v_ct is incorrect!\n", f->label);
      return 1;
    }

    if(verbosity >= 2) {
      printf("  dumping position for fluxon %ld, vertex %ld: (%g,%g,%g)\n",f->label, v->label, v->x[0],v->x[1],v->x[2]);
    }

    *(NUM *)dex = v->x[0];
    dex += sizeof(NUM);

    *(NUM *)dex = v->x[1];
    dex += sizeof(NUM);

    *(NUM *)dex = v->x[2];
    dex += sizeof(NUM);
  }

  if(verbosity >= 1)
    printf("Dumping %ld bytes of position data for fluxon %ld...\n",len,f->label);

  if(dex - binary_buffer != len) {
    printf("binary_dump_flpos: expected %ld bytes, but offset is %ld\n",len, dex - binary_buffer);
  }

  binary_dump_field(fd, BD_POSITION, len, binary_buffer);

  return 1;
}

int binary_read_flpos( long size, char *buf, WORLD *w) {
  char *ptr = buf;
  char *me = "binary_read_flpos";
  long v_ct;
  long label;
  FLUXON *f;
  VERTEX *v;
  long i;
  long verbosity = w->verbosity;

  if ( *(long *)ptr != 1 ) { // Check version number
    fprintf(stderr,"%s: wrong version  (expected %d, got %ld)\n", me, 1, *(long *)ptr);
    return 1;
  }
  ptr += sizeof(long);

  label = *(long *)ptr; // Read fluxon label
  ptr += sizeof(long);

  f = *(FLUXON **)ptr;   // Read fluxon pointer
  ptr += sizeof(FLUXON *);


  if( f->label != label ) { 
    fprintf(stderr, "%s: fluxon mismatch! Followed pointer, didn't get the label %ld\n",me, label);
    return 1;
  }

  v_ct = *(long *)ptr;  // Read expected number of vertices
  ptr += sizeof(long);
  if(v_ct != f->v_ct) {
    fprintf(stderr,"%s: fluxon %ld has %ld vertices, not %ld as advertised...\n",me, label, f->v_ct, v_ct);
    return 1;
  }

  if(*((long *)ptr) != 0xAABBCCDD) {
    printf("Problem with read_binary_flpos -- fence not found!\n");
    return 1;
  }
  ptr += sizeof(long);


  if(w->verbosity > 1)
    printf("Reading data for %ld vertices in fluxon %ld\n",f->v_ct,f->label);

  for(v=f->start, i=0; i<v_ct; v=v->next,i++) {


    v->x[0] = *(NUM *)ptr;
    ptr += sizeof(NUM);
    
    v->x[1] = *(NUM *)ptr;
    ptr += sizeof(NUM);
    
    v->x[2] = *(NUM *)ptr;
    ptr += sizeof(NUM);

    if(verbosity >= 2) {
      printf("  reading fluxon %ld, vertex %ld: (%g,%g,%g)\n",f->label, v->label, v->x[0],v->x[1],v->x[2]);
    }

    
  }

  return 0;
}


  
/******************************
 * binary_dump_flstep
 * 
 * Dumps just the planned step vector for each vertex in a fluxon.
 * Used for parallelization.
 */

int binary_dump_flstep( int fd, FLUXON *f) {
  long len = 
    sizeof(NUM) * 3 * f->v_ct              //  vectors
    + 3*sizeof(long)                       //  three longs: version, fluxon label, vertex count
    + sizeof(FLUXON *)                     //  fluxon pointer
    + 8*sizeof(NUM)                        // accumulators
    + sizeof(long)                         // fence
    ;
  char *dex;
  VERTEX *v;
  WORLD *w = f->fc0->world;
  int verbosity = w->verbosity;
  long i;

  check_binary_buf(len);
  dex = binary_buffer;
  *(long *)dex = 1; // Version number
  dex += sizeof(long);

  *(long *)dex = f->label; // Label of the fluxon
  dex += sizeof(long) ;

  *(FLUXON **)dex = f; // fluxon pointer itself
  dex += sizeof(FLUXON *);

  *(long *)dex = f->v_ct; // Number of vertices to expect
  dex += sizeof(long);

  *(NUM *)dex = w->f_min;
  dex += sizeof(NUM);
  
  *(NUM *)dex = w->f_max;
  dex += sizeof(NUM);
  
  *(NUM *)dex = w->fr_min;
  dex += sizeof(NUM);
  
  *(NUM *)dex = w->fr_max;
  dex += sizeof(NUM);

  *(NUM *)dex = w->ca_min;
  dex += sizeof(NUM);

  *(NUM *)dex = w->ca_max;
  dex += sizeof(NUM);

  *(NUM *)dex = w->ca_acc;
  dex += sizeof(NUM);

  *(NUM *)dex = w->ca_ct;
  dex += sizeof(NUM);

  *(long *)dex = 0x01234567;
  dex += sizeof(long);

  if(verbosity >= 2) {
    printf("binary_dump_flstep: dumping shifts for fluxon %ld:", f->label);
  }

  for(v=f->start, i=0; v; v=v->next,i++) {
    if(i >= f->v_ct) {
      fprintf(stderr,"binary_dump_flstep: fluxon %ld v_ct is incorrect!\n", f->label);
      return 1;
    }

    if(verbosity >= 2) {
      printf(" (%g,%g,%g)\t",v->plan_step[0], v->plan_step[1], v->plan_step[2]);
    }

    *(NUM *)dex = v->plan_step[0];
    dex += sizeof(NUM);

    *(NUM *)dex = v->plan_step[1];
    dex += sizeof(NUM);

    *(NUM *)dex = v->plan_step[2];
    dex += sizeof(NUM);
  }
  if(verbosity >= 2) {
    printf("\n");
  }
  if(verbosity >= 1)
    printf("binary_dump_flstep: dumping %ld bytes of position data for fluxon %ld...\n",len,f->label);

  if( dex - binary_buffer  != len) {
    printf("binary_dump_flstep: len was %ld, offset is %ld\n",len, dex - binary_buffer);
  }

  binary_dump_field(fd, BD_STEP, len, binary_buffer);
  return 1;
}

int binary_read_flstep( long size, char *buf, WORLD *w) {
  char *ptr = buf;
  char *me = "binary_read_flstep";
  long v_ct;
  long label;
  FLUXON *f;
  VERTEX *v;
  int verbosity = w->verbosity;
  long i;

  if ( *(long *)ptr != 1 ) { // Check version number
    fprintf(stderr,"%s: wrong version  (expected %d, got %ld)\n", me, 1, *(long *)ptr);
    return 1;
  }
  ptr += sizeof(long);

  label = *(long *)ptr; // Read fluxon label
  ptr += sizeof(long);

  f = *(FLUXON **)ptr;   // Read fluxon pointer
  ptr += sizeof(FLUXON *);


  if( f->label != label ) { 
    fprintf(stderr, "%s: fluxon mismatch! Followed pointer, didn't get the label %ld\n",me, label);
    return 1;
  }

  v_ct = *(long *)ptr;  // Read expected number of vertices
  ptr += sizeof(long);
  if(v_ct != f->v_ct) {
    fprintf(stderr,"%s: fluxon %ld has %ld vertices, not %ld as advertised...\n",me, label, f->v_ct, v_ct);
    return 1;
  }


  // Accumulate minmax stuf finto the current accumulator
  if(w->f_min > *(NUM *)ptr || w->f_min < 0)
    w->f_min = *(NUM *)ptr;
  ptr += sizeof(NUM);

  if(w->f_max < *(NUM *)ptr || w->f_max < 0)
    w->f_max = *(NUM *)ptr;
  ptr += sizeof(NUM);

  if(w->fr_min > *(NUM *)ptr || w->fr_min < 0) 
    w->fr_min = *(NUM *)ptr;
  ptr += sizeof(NUM);

  if(w->fr_max < *(NUM *)ptr || w->fr_max < 0)
    w->fr_max = *(NUM *)ptr;
  ptr += sizeof(NUM);

  if(w->ca_min > *(NUM *)ptr || w->ca_min < 0) 
    w->ca_min = *(NUM *)ptr;
  ptr += sizeof(NUM);

  if(w->ca_max < *(NUM *)ptr || w->ca_max < 0) 
    w->ca_max = *(NUM *)ptr;
  ptr += sizeof(NUM);

  w->ca_acc += *(NUM *)ptr;
  ptr += sizeof(NUM);

  w->ca_ct += *(NUM *)ptr;
  ptr += sizeof(NUM);
						
  if(*(long *)ptr != 0x01234567) {
    int i;
    printf("Problem with read_binary_flstep -- fence not found!\n");
    printf("looking for 01 23 45 67 (bigendian) or 67 45 23 01 (littleendian); dumping from ptr-16 to ptr+15:\n");
    for(i=-16; i<15;i++) {
      printf(" %x",*((unsigned char *)ptr+i));
    }
    printf("\n");
	   
    return 1;
  }
  ptr += sizeof(long);

  if(verbosity > 1)
    printf("binary_read_flstep: Reading data for %ld vertices in fluxon %ld\n",f->v_ct,f->label);

  for(v=f->start, i=0; i<v_ct; v=v->next,i++) {

    v->plan_step[0] = *(NUM *)ptr;
    ptr += sizeof(NUM);

    v->plan_step[1] = *(NUM *)ptr;
    ptr += sizeof(NUM);
    
    v->plan_step[2] = *(NUM *)ptr;
    ptr += sizeof(NUM);
    if(verbosity > 1) {
      printf(" (%g,%g,%g)\t",v->plan_step[0],v->plan_step[1],v->plan_step[2]);
    }
  }

  if(verbosity > 1) {
    printf("\n");
  }
  return 0;
}


