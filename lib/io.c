/**********************************************************************
 * io.c -- I/O routines for FLUX
 *
 * Mostly text I/O.
 *
 * This file is part of FLUX, the Field Line Universal relaXer.
 * Copyright (c) Southwest Research Institute, 2004-2007
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
 * This file is part of the FLUX 2.0 release (31-Oct-2007).
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <regex.h>
#include <math.h>

#include "data.h"
#include "io.h"
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

	if( l0<0 ) {
	  printf("WARNING: automagic FC %ld is given in file -- ignoring this gaffe. (Probably OK)\n",l0);
	  break;
	}
	

	/* New stuff handler */

	FLUX_CONCENTRATION *fc;

	fc = new_flux_concentration( world,x0[0],x0[1],x0[2],flux0, l0 );

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
      sprintf(lscan,"%%*s %%ld %%ld %%ld %%ld %%ld %%%sf",NUMCHAR,NUMCHAR,NUMCHAR);
    }

    n = sscanf(s,lscan, &fl0, &vl0, &vl1, &l0, &l1, &flux0);
    if(n != 6) {
      // Fall back to older form 

      if(!*lscan2) {
	sprintf(lscan2,"%%*s %%ld %%ld %%ld %%%sf",NUMCHAR);
      }
      
      n = sscanf(s,lscan2, &fl0, &l0, &l1, &flux0);
      
      if(n != 4) {
	badstr = "Couldn't parse LINE line";
	n=0;
      } else {
	vl0 = vl1 = 0;
      }
    }

    if(n != 0) {

      if(fl0 <= -11 || fl0 >0) {
	FLUX_CONCENTRATION *fc0, *fc1;
	fc0 = tree_find(world->concentrations, l0, fc_lab_of, fc_ln_of);
	fc1 = tree_find(world->concentrations, l1, fc_lab_of, fc_ln_of);
	//printf("line %d: fc %d (%x) - fc %d (%x)\n",fl0,l0,fc0,l1,fc1);
	if(!fc0 || !fc1) {
	  char *badbuf = (char *)localmalloc(BUFSIZ,MALLOC_MISC);
	  sprintf(badbuf,"Found a fluxon specifier between concentrations %ld and %ld, but they \ncame up %ld and %ld in tree_find (one doesn't exist)!  \nThis error message leaked %d bytes (don't let it happen again!)\n",l0,l1,fc0,fc1,BUFSIZ);
	  badstr = badbuf;
	  
	} else if(fc0->flux * fc1->flux >= 0) {
	  badstr = (char *)localmalloc(BUFSIZ,MALLOC_MISC);
	  sprintf(badstr,"This fluxon connects two flux concentrations of the same sign flux, or\none of its flux tubes has zero flux. Line %ld; concentrations %ld (%g) and %ld (%g)\n",fl0, flux0, l0, fc0->flux, l1, fc1->flux);
	  break;
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
		
		v0 = new_vertex( vl0, fc0->x[0],fc0->x[1],fc0->x[2],f0);
		v1 = new_vertex( vl1, fc1->x[0],fc1->x[1],fc1->x[2],f0);
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
    }
    
    break;
    
  case 'V': 

    if(s[1]=='E' || s[1]=='e') {

      /* VERTEX <fl_lab> <vert_lab> <pos> <x> <y> <z> 
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
	  (stderr,"Warning: V_NEIGHBOR line asked for targets %d and %d; at least one of 'em isn't defined...\n",vlab,nlab);
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
	    long concurrency;
	    if( 1 > (nval = sscanf(s+off,"%D",&concurrency))) {
	      badstr = "couldn't parse GLOBAL CONCURRENCY line\n";
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
		long foo;
		if( sscanf(fname,"0x%x",&foo) )
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
		printf("s+off: %s\nphscan: %s\n",s+off,phscan);
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
		printf("s+off: %s\nphscan: %s\n",s+off,phscan);
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
	      badstr = "Couldn't parse four values from GLOBAL OPENcd ";
	    } else {
	      world->fc_ob->x[0] = world->fc_oe->x[0];
	      world->fc_ob->x[1] = world->fc_oe->x[1];
	      world->fc_ob->x[2] = world->fc_oe->x[2];
	      world->fc_ob->locale_radius = world->fc_oe->locale_radius;
	    }
	  }
	  break;
	case 'R':  /* ARD Added - read rel_step */
	  if (sscanf(s, "%*s %*s %d", &(world->rel_step)) != 1) {
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
		char mscan[100];
		sprintf(mscan,"%%*s %%%sf",NUMCHAR);
		sscanf(mscan,scanvar);
	      } else {
		badstr = s;
	      }

	    }
	    break;

	  case 'T': case 't': /* GLOBAL STATE */
	    sscanf(s+off,"%d",&(world->state)) || 
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
  
  do{
    if( a == NULL ) {
      a = new_world();
    }

    if(s = next_line(file)) {
      error = footpoint_action(a, s);
    }
  } while (s && !error);
  
  world_update_ends(a);

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

  fprintf(file,"GLOBAL RSTEP %d\n", world->rel_step);
  fprintf(file,"GLOBAL DTAU %lf\n", world->dtau);

  /* ARD Added code for writing out global coefficients */
 
  fprintf(file,"GLOBAL COEFFICIENTS %d", world->n_coeffs);
  for (i=0;i<world->n_coeffs;i++) {
    fprintf(file, " %f",  world->coeffs[i]); 
  }
  fprintf(file,"\n");
  fprintf(file,"GLOBAL OPEN %g %g %g %g %d\n",
	  world->fc_ob->x[0],
	  world->fc_ob->x[1],
	  world->fc_ob->x[2],
	  world->fc_ob->locale_radius,
	  world->auto_open
	  );

  fprintf(file,"GLOBAL SKEW_HANDLING %d\n", world->handle_skew);

  fprintf(file,"GLOBAL CONCURRENCY %d\n",world->concurrency);
  
  fprintf(file,"\n\n"); /* leave an extra space after the globals */
    
  /* Write out all flux concentration locations */
  fprint_tree(file, world->concentrations, fc_lab_of, fc_ln_of, 0, fprint_fc_line_nonneg);

  /* Output the field lines */
  fputc('\n',file);
  fprint_tree(file, world->concentrations, fc_lab_of, fc_ln_of, 0, fprint_fls_by_fc);

  /* Output neighbor relations */
  fputc('\n',file);
  fprint_tree(file, world->vertices, v_lab_of, v_ln_of, 0, fprint_v_nbors);
}    
    
/**********************************************************************
 **********************************************************************
 ** FREEZING ROUTINES --
 ** 
 ** These generate snippets of Perl code that reconstitute the 
 ** World in its current state.
 ** (FIXME: write this stuff!)
 **    char *freeze_vertex
 **    char *freeze_fluxon
 **    char *freeze_fc
 **    char *freeze_world
 **/

    
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
    sprintf(fc_line_format,"NEW\t%%ld\t%%%sg\t%%%sg\t%%%sg\t%%%sg\n",NUMCHAR,NUMCHAR,NUMCHAR,NUMCHAR);
  }
  fprintf(f, fc_line_format,fc->label, fc->x[0], fc->x[1], fc->x[2], fc->flux);
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
    sprintf(fl_format,"LINE\t%%ld\t%%ld\t%%ld\t%%ld\t%%ld\t%%%sg\n",NUMCHAR);
    sprintf( v_format,"\tVERTEX\t%%ld\t%%ld\t%%ld\t%%%sg\t%%%sg\t%%%sg\n",
	     NUMCHAR,NUMCHAR,NUMCHAR);
  }


  fprintf(f, fl_format
	  , fl->label
	  , fl->start ? fl->start->label : 0
	  , fl->end ? fl->end->label : 0
	  , fl->fc0 ? fl->fc0->label : -999
	  , fl->fc1 ? fl->fc1->label : -999
	  , fl->flux);
  
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
    

/**********************************************************************
 **********************************************************************
 */
/**********************************************************************
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
 * To read a file, call binary_dump_read with an open descriptor and the 
 * WORLD to which it should apply.  It will read and parse all the elements 
 * in the file, modifying the arena as they are processed, and close the 
 * file descriptor on exit.
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
  
  if( i != 3*sizeof(long) ) 
    return 1;
  
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
      fprintf(stderr,"Failed to grow the binary I/O buffer - tried for %d bytes",binary_buflen);
      exit(1);
    }
  }
}

/******************************
 * dump_end - send an end marker to the file.  Does not close the file descriptor.
 */
int binary_dump_end(int fd) {
  binary_dump_field( fd, BD_END, 0, 0 );
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
  for(v=f->start, i=0; v; v=v->next, i++) {

    // i is just for error checking - complain here
    if( i >= f->v_ct ) {
      fprintf(stderr,"binary_dump_fluxon_pipe: Wonky fluxon %d - has more vertices than v_ct (%d)! I quit\n",f->label, f->v_ct);
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
    vd->neighbors_n = v->neighbors.n;

    neighbors_found += v->neighbors.n;
    if(neighbors_found > neighbors_allowed_for){
      fprintf(stderr,"binary_dump_fluxon_pipe: Whoa!  Never expected to see this.  Found %d neighbors along fluxon %d - that's over %d per vertex! I give up.\n",f->label, neighbors_found, MAX_AV_NEIGHBORS);
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
}



/**********************************************************************
 * binary_read_dumpfile - read in a complete binary dump
 *
 */
int binary_read_dumpfile ( int fd, WORLD *w ) {
  long hdrbuf[3];
  int ct;
  int pos = 0;
  long type;
  long len;
  char *me = "FLUX: binary_read_dumpfile";
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
	return 1;
      }
    }


    // Parse header: check for fence, and get data type and length
    pos += ct;
    if(hdrbuf[0] != BD_FENCE) {
      fprintf(stderr,"%s: failed to find fence (expected %x, got %x), position %d", me, BD_FENCE, hdrbuf[0]);
      return 2;
    }
    type = hdrbuf[1];
    len = hdrbuf[2];
    if(type <=0 || type > BD_MAX_TYPENO) {
      fprintf(stderr,"%s: failed to find valid type (expected 1-%d, got %d), position %d",me,BD_MAX_TYPENO,type);
      return 3;
    }
    check_binary_buf( len );

    // Read the packet (try twice if necessary)
    ct = read( fd, binary_buffer, len );
    if(ct != len) {
      int ct2 = -10; // any ol' large negative number

      if(ct>=0)
	ct2= read(fd, binary_buffer+ct, len-ct);
      
      if(ct2+ct != len)  {
	fprintf(stderr,"%s: failed to read %d bytes (got %d bytes from position %d, type %d; second try yielded %d bytes for %d total)\n",me,len,ct, pos, type,ct2,ct+ct2);
	perror("read returned the error");
	return 4;
      }
    }

    pos += len;

    if(w->verbosity) {
      printf("read_dumpfile: found fence; type=%d, len=%d\n",hdrbuf[1],hdrbuf[2]);
    }
    
    switch(type) {
    case BD_END:
      close(fd);
      return 0;
      break;
    case BD_FLUXON_PIPE:
      if(binary_read_fluxon_pipe(len, binary_buffer, w)) {
	fprintf(stderr,buf);
	return 1;
      }
      break;
    default:
      fprintf(stderr,"WARNING: typecode %d not implemented - you should never see this message from binary_read_dumpfile, in io.c...\n",type);
      break;
    }
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

  if(w->verbosity) {
    printf("read_fluxon_pipe...");
    fflush(stdout);
  }



  if(size< 3*sizeof(long)) {
    fprintf(stderr, "%s: inconceivable packet size %d is too small!\n",me,size);
    return 1;
  }

  if( (*(long *)dex) != 3) { // Check version number
    fprintf(stderr, "%s: packet is the wrong version (%d, I am version %d)\n",me,*(long *)dex,3);
    return 2;
  }
  dex+= sizeof(long);

  f_lab = *(long *)dex; // Label of the fluxon
  dex += sizeof(long);

  f = *(FLUXON **)dex;  // Fluxon pointer itself
  dex += sizeof(FLUXON *);
  
  if(f->label != f_lab) {
    fprintf(stderr,"%s: found bogus fluxon pointer in dump! (expected fluxon %d, got %d)",me,f_lab, f->label);
    return 3;
  }

  v_ct = *(long *)dex; // Vertex count
  dex += sizeof(long);
  
  if(f->v_ct != v_ct) {
    fprintf(stderr,"%s: vertex count doesn't agree in fluxon %d (expected %d from file, found %d in fluxon)\n",me,f_lab,v_ct, f->v_ct);
    return 4;
  }

  if(w->verbosity) {
    printf(" r: f=%d,vct=%d    ",f->label,f->v_ct);
    fflush(stdout);
  }

  for(v=f->start, i=0; i<v_ct; v=v->next, i++) {
    vd = (VERTEX_PIPE_FIXED *)dex;
    dex += sizeof(VERTEX_PIPE_FIXED);

    if(vd->label != v->label) {
      fprintf(stderr,"%s: vertex label mismatch in fluxon %d, pos %d - file expected %d, found %d\n",me,f_lab,i,vd->label,v->label);
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

    if(w->verbosity) {
      printf("vertex %d: found %d neighbors...\n",v->label,vd->neighbors_n);
    }

    // First - collect the neighbors in a local dumblist.
    dl->n = 0;
    for(j=0; j<vd->neighbors_n; j++) {
      dumblist_add( dl, *(void **)dex );
      dex += sizeof(void *);
    }

    // Next - check the fence!
    if( *(long *)dex != VERTEX_END_FENCE ) {
      fprintf(stderr,"%s: missed end fence while copying vertex %d (label %d) into fluxon %d - giving up!\n",me, i, v->label, f->label);
      return 6;
    }
    dex += sizeof(long);

    ////// Now - merge the neighbors with the local neighbor list. 

    // Delete current vertex neighbors that aren't in the new list.
    // This is complicated by the need to maintain back-links.
    for(j=0; j<v->neighbors.n; j++) {
      int found = 0;
      for(k=0; k<dl->n && !found; k++) 
	found = ( dl->stuff[k] == v->neighbors.stuff[j] );
      if(!found) {
	// Delete the back link first, then remove the item from the neighbor list.
	dumblist_delete(  &( ((VERTEX *)(v->neighbors.stuff[j]))->nearby ), v );
	dumblist_rm( &(v->neighbors), j );
	// Decrement counter - last element moved to this slot, so we have to re-do the slot.
	j--;
      }
    }
    
    // Finally - add any items from the new list that aren't already neighbors.
    for(k=0; k<dl->n; k++) {
      int found = 0;
      for(j=0; j<v->neighbors.n && !found; j++) 
	found = ( dl->stuff[k] == v->neighbors.stuff[j] );
      if(!found) {
	dumblist_add( &(v->neighbors), dl->stuff[k] );
	dumblist_add( &( ((VERTEX *)(dl->stuff[k]))->nearby ), v );
      }
    }

    if(w->verbosity) {
      printf("after transfer: %d neighbors.      ",v->neighbors.n);
    }

  }

  return 0;
}

    
  
    
  
    
