/* data.c 
 *
 * Basic data manipulation routines and definitions for FLUX -- how to 
 * handle VERTEXs and FLUXONs etc.
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
 *
 *
 * Functions:
 *   new_fluxon - Generate a new field line structure (empty)
 *   new_vertex  - Generate a vertex from X,Y,Z, and fluxon references.
 *   add_vertex - Given a vertex and an X,Y, and Z location, create
 *             a new vertex and link it into the list at that location.
 *   find_neighbor_candidates - Find the neighbor-candidates of a vertex.
 *             The vertex gets the union of its neighbor list, its 
 *             fluxon partners' lists, and the neighbors of everybody
 *             on that list.  If the vertex is the first on its fluxon,
 *             then its flux concentrations are searched for neighbors.  If
 *             it is the first fluxon on the flux concentration, then 
 *             the neighbors its Delauney neighbors.  (These are currently
 *             found by exhaustive Voronoi cell construction; this is slow
 *             but should work for now.)
 *   delete_vertex  - Unlink and delete a vertex. (unlinks neighbors too)
 *   (NOT YET) kill_fluxon - Unlink and delete all of a fluxons vertexs.
 *
 *  clear_links - Initialize a tree-link data structure.
 *  tree_find   - Find a node in a tree (there are several kinds!).
 *  tree_insert - insert a node into a tree.
 *  tree_binsert - (the call of choice) - Insert a node into a tree, 
 *                   and rebalance the tree if necessary.
 *  tree_unlink - Remove a node from a tree.
 *  tree_bunlink - (the call of choice) - Remove a node from a tree,
 *                   and rebalance the tree if necessary.
 *  tree_balance - Balance a tree.
 *  tree_top    - Return a pointer to the top of a tree.
 *
 *  dumblist_add    - Add an item to a dumblist (no duplicates allowed)
 *  dumblist_delete - Delete an item from a dumblist.
 *  dumblist_snarf  - Copy a dumblist into another one.
 *  
 */
#include "data.h"
#include "physics.h" /* for declaration of force subroutines */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/******************************
 * Boundary condition names, by number
 */
char *BOUNDARY_NAMES[] = {
  "NONE",
  "PLANE",
  "SPHERE",
  "CYL"
};




/**********************************************************************
 * barf
 * 
 * Calls perror and exits, with some indication of where you started from.
 */
static char program_name[80];
void init_barf(char *whoami) {
  strncpy(program_name,whoami,79);
  program_name[79]=0;
}

void barf(int barf_on_op, char *where) {
  char buf[10240];
  switch(barf_on_op) {
  case BARF_MALLOC:
    fprintf(stderr,"Problems with Malloc!\n");
    break;
  default:
    break;
  }
  sprintf(buf,"%s(%s)",program_name,where);
  perror(buf);
  exit(barf_on_op);
}


/**********************************************************************
 * new_label
 * Generates a new long-int label for a fluxon.  With 2^64 labels
 * available, we don't bother trying to conserve -- just plough through
 * them in sequence.  Fieldlines and vertexs share a label numbering
 * system.   If you pass in 0, you get the next sequentially available
 * label.  If you pass in a positive label number, you get that number back but
 * also the last_label cache is set to at least that number.  There are 
 * potential collision problems here, but the solution is not to 
 * allocate any default-numbered nodes (e.g. from simulation) until after
 * all the set-numbered nodes (e.g. from a file) have been allocated.
 * 
 * If you pass in -1, you get back the current maximum-label cache
 * and it is not incremented.
 *  
 * Labels seem mainly useful for I/O (storing stuff to files), but
 * they come in handy for debugging too.  It might be possible 
 * eventually to remove them from the code.
 */
long new_label(long request){
  static long max_label=0;
  if(request <= -1) return max_label;
  if(!request) return ++max_label;

  if(max_label < request)
    max_label=request;
  return request;
}

/**********************************************************************
 * new_vertex_label -- same as new_label, except for VERTEXes.
 */
static long max_vertex_label=0;
long new_vertex_label(long request) {
  long out;
  switch(request) {
  case -1: 
    out = max_vertex_label;
    break;
  case -2: 
    out = 0;
    break;
  case 0:
    out = ++max_vertex_label;
    break;
  default:
    out = request;
    if(max_vertex_label < request)
      max_vertex_label = request;
    break;
  }
  return out;
}

/***********************************************************************
 * new_fluxon 
 */
FLUXON *new_fluxon(      NUM flux, 
			 FLUX_CONCENTRATION *c1, 
			 FLUX_CONCENTRATION *c2, 
			 long label,
			 char plasmoid)
{
  int f_no = 0;
  struct FLUXON *nf;
  
  nf = (struct FLUXON *)localmalloc(sizeof(struct FLUXON),MALLOC_FLUXON);

  if(!nf) barf(BARF_MALLOC,"new_fluxon");
  
  nf->flux            = flux;
  nf->label           = new_label(label);
  nf->start           = 0;
  nf->end             = 0;
  nf->v_ct            = 0;

  nf->start_b         = 0;
  nf->end_b           = 0;

  clear_links(&(nf->all_links));
  nf->all_links.sum = nf->flux;
  
  if(c1 == NULL) {
    fprintf(stderr,"new_fluxon:  no starting flux conc.!  Proceeding anyway...\n");

  } 
  clear_links(&(nf->start_links));
  nf->start_links.sum = nf->flux;
  nf->fc0 = c1;


  if(c2 == NULL) {
    fprintf(stderr,"new_fluxon: no ending flux conc.!  Proceeding anyway...\n");
  }

  clear_links(&(nf->end_links));
  nf->end_links.sum = nf->flux;
  nf->fc1 = c2;

  nf->plasmoid = 0;

  return nf;
}

/***********************************************************************
 * new_vertex
 *
 * Returns a vertex given an X,Y,Z, and fluxon. 
 * The neighbor list is empty -- call set_neighbors to get that.
 *
 */

inline VERTEX *new_vertex(long label, NUM x, NUM y, NUM z, FLUXON *fluxon) { 
  VERTEX *tp;
  int i;

  tp = (VERTEX *)localmalloc(sizeof(VERTEX),MALLOC_VERTEX);

  if(!tp) barf(BARF_MALLOC,"new_vertex");
  tp->line = fluxon;
  tp->prev = 0;
  tp->next = 0;

  tp->x[0] = x;
  tp->x[1] = y;
  tp->x[2] = z;

  dumblist_clear( &(tp->neighbors) );
  dumblist_clear( &(tp->nearby) );


  tp->b_mag = 0;
  tp->b_vec[0] = tp->b_vec[1] = tp->b_vec[2] = 0;

  tp->label = new_vertex_label(label);

  return tp;
}  

/**********************************************************************
 * new_flux_concentration
 * 
 * Returns a flux concentration at the given location
 *
 */

FLUX_CONCENTRATION *new_flux_concentration(
					  WORLD *world,
					   NUM x, NUM y, NUM z,
					   NUM flux,
					   long label
					   ) {
  FLUX_CONCENTRATION *fc;
  
  fc = (FLUX_CONCENTRATION *)localmalloc(sizeof(FLUX_CONCENTRATION),MALLOC_VERTEX);
  if(!fc) barf(BARF_MALLOC,"new_flux_concentration");
  fc->world = world;
  fc->flux = flux;
  fc->label = new_label(label);
  fc->lines = (FLUXON *)0;
  fc->x[0] = x;
  fc->x[1] = y;
  fc->x[2] = z;
  fc->locale_radius = 0; /* locale_radius isn't used just now */
  clear_links(&(fc->links));
  fc->links.sum = fc->flux;
  return fc;
}

/**********************************************************************
 * new_world
 *
 */
WORLD *new_world() {
  WORLD *a = (WORLD *)localmalloc(sizeof(WORLD),MALLOC_WORLD);
  if(!a) barf(BARF_MALLOC,"new_world");
  
  a->frame_number = 0;
  a->state = WORLD_STATE_NEW;
  a->concentrations = NULL;
  a->lines = NULL;

  a->photosphere.type = 0;
  a->photosphere.plane = NULL;

  /* Initialize the two dummy vertices for the mirroring */
  a->image = new_vertex(1,0,0,0,NULL);
  a->image2 = new_vertex(2,0,0,0,NULL);
  a->image->next = a->image2;
  a->image2->prev = a->image;

  a->verbosity = 0;  /* No verbose printouts by default */

  /* Put in a sensible default force list */
  a->f_funcs[0] = b_eqa;
  a->f_funcs[1] = f_curv_m;
  a->f_funcs[2] = f_p_eqa_radial; 
  a->f_funcs[3] = f_vert;
  a->f_funcs[4] = 0;

  a->f_over_b_flag = 1; 

  return a;
}


/**********************************************************************
 **********************************************************************
 ***** 
 ***** List routines: handle the singly linked list structures in the
 ***** fluxons themselves
 *****
 */

/**********************************************************************
 * unlink_vertex: Remove the given vertex from its list and its
 * neighborhood.  Neighborhood information gets sloshed around.
 */
void unlink_vertex(VERTEX *v) {
  int i;

  if(!v)
    return;
  if(!v->prev || !v->next) {
    fprintf(stderr,"unlink_vertex ignoring an end condition...\n");
    return;
  }

  v->prev->next = v->next;
  v->next->prev = v->prev;
  
  for(i=0;i<v->neighbors.n;i++) {
    VERTEX *a = ((VERTEX **)(v->neighbors.stuff))[i];
    dumblist_delete( &(a->nearby), v);
    dumblist_add(    &(a->nearby), v->next);
    dumblist_add(    &(a->nearby), v->prev);
  }
  dumblist_snarf(&(v->prev->neighbors), &(v->neighbors));
  dumblist_snarf(&(v->next->neighbors), &(v->neighbors));

  for(i=0;i<v->nearby.n;i++) {
    VERTEX *a = ((VERTEX **)(v->nearby.stuff))[i];
    dumblist_delete(&(a->neighbors), v);
    dumblist_add(   &(a->neighbors), v->next);
    dumblist_add(   &(a->neighbors), v->prev);
  }
  dumblist_snarf(&(v->prev->nearby), &(v->nearby));
  dumblist_snarf(&(v->next->nearby), &(v->nearby));
}

void delete_vertex(VERTEX *v) {
  if(!v) 
    return;
  unlink_vertex(v);
  if(v->neighbors.stuff) {
    localfree(v->neighbors.stuff);
    v->neighbors.stuff=0;
  }
  if(v->nearby.stuff) {
    localfree(v->nearby.stuff);
    v->nearby.stuff=0;
  }
  localfree(v);
}

/**********************************************************************
 * add_vertex_pos: stick vertex <v> into fluxon <f> at position <pos>.
 * Return status 0=success, 1=error.  
 *
 * Input:  Node counting starts at 0 (which is a trivial node); 
 * -1 implies the last node (which is the other trivial node). 
 * -2 adds to the end of the list of middle nodes.  1 adds to the
 * beginning of the list of middle nodes.  Numbers are truncated such
 * that it's impossible to scrozzle the final node by feeding in a 
 * positive number -- you have to use -1 to get that.
 *
 * High numbers 
 */

int add_vertex_pos(FLUXON *f, long pos, VERTEX *v) {
  
  if(!f || !v) 
    return 1;

  if(v->next != 0 || v->prev != 0) {
    fprintf(stderr,"Hey! add_vertex_pos got a *list*, not a vertex.  Not supported yet.\n");
    return 1;
  }

  /* Handle the empty-fluxon case */
  if(f->start == 0 || f->end == 0){
    if(f->start == 0 && f->end == 0) {
      f->start = v;
      f->end = v;
      f->v_ct = 1;
      v->next = 0;
      v->prev = 0;
      return 0;
    } else {
      fprintf(stderr, "Encoutered a fluxon with at START but no END, or vice versa!\n");
      return 1;
    }
  }

  /* Handle counting-forwards */
  if(pos >= 0) {
    int i;
    VERTEX *v0;

    if( f->v_ct <= pos ) {
      if(!f->end->prev) {
	fprintf(stderr,"add_vertex_pos got a 1-node fl with a large pos. position.  Not kosher!\n");
	return 1;
      }
      return add_vertex_after(f, f->end->prev, v);
    }

    /* If it's faster to count backwards, count backwards */
    if(pos > f->v_ct/2) {
      for(i=f->v_ct-1, v0 = f->end; i>pos; i++)
	v0=v0->prev;
    } else {
      for(i = 0, v0 = f->start; v0->next && i<pos; i++)
	v0 = v0->next;
    }
    return add_vertex_after(f, v0->prev, v);
  }

  /* handle counting-backwards (no check because everything else returns!) */
  /* if(pos<0) */
  {
    int i;
    VERTEX *v0;

    if( f->v_ct < -pos ) {
      return add_vertex_after(f, NULL, v);
    }
    
    /* Add optional counting-forwards here for minor speed improvement */

    for(i = -1, v0 = f->end; v0->prev && i>pos; i--) 
      v0 = v0->prev;
    
    return add_vertex_after(f, v0, v);
  }
}

/**********************************************************************
 * add_vertex_after: Stick vertex <v> into fluxon <f>, just
 * AFTER <neighbor>.  If <neighbor> is 0, then <v> goes at the beginning
 * of the list.  Return an error flag.  The newly inserted vertex
 * gets the prior and next vertices' neighbor lists, too.
 */
int add_vertex_after(FLUXON *f, VERTEX *neighbor, VERTEX *v) {

  /* Insert at beginning of list */
  if(!neighbor) {
    v->prev = 0;
    v->next = f->start;
    if(f->start) {
      f->start->prev = v;  /* There was at least one other vertex */
      f->v_ct++;
    }
    else {
      f->end = v;          /* This is the first vertex */
      f->v_ct=1;           /* Enforce the v_ct to be correct */
    }
    f->start = v;
  }

  else {
    if(neighbor->line != f) {
      fprintf(stderr,"add_vertex_after: neighbor is not on its fluxon!\n");
      return 1;
    }
    
    /* Insert at end of list */
    if(!neighbor->next) {
      if(f->end != neighbor) { /* Print an error and abort if inconsistent */
	fprintf(stderr,"Strange vertex list found by add_vertex_after\n");
	exit(1021);
      }
      v->prev = neighbor;
      v->next = 0;
      neighbor->next = v;
      f->end = v;
      f->v_ct++;
    }
    
    else {
      /* Insert in middle of list */
      v->next = neighbor->next;
      v->prev = neighbor;
      neighbor->next->prev = v;
      neighbor->next = v;
      f->v_ct++;
    }  
  }

  /* Fix up neighbor lists */
  if(v->prev) 
    dumblist_snarf(&(v->neighbors),&(v->prev->neighbors));
  if(v->next)
    dumblist_snarf(&(v->neighbors),&(v->next->neighbors));
  {
    int i;
    for(i=0;i<v->neighbors.n;i++) {
      dumblist_add( &(((VERTEX *)(v->neighbors.stuff[i]))->nearby), (void *)v);
      dumblist_add( &(((VERTEX *)(v->neighbors.stuff[i]))->neighbors), (void *)v);
      dumblist_add( &(v->nearby), v->neighbors.stuff[i]);
    }
  }

  return 0;
}

/**********************************************************************
 **********************************************************************
 *****
 *****  The tree routines:  handle multiply linked tree structures
 *****  with long-int labels and [left,right,up] link structures.
 *****
 *****  If you use a different struct type with these routines,
 *****  make sure that the first element is a NUM that can be
 *****  running-summed over.
 *****
 *****
 */

/**********************************************************************
 * clear_links
 * Clears out the links for a new data structure.
 */
void clear_links(LINKS *links) {
  links->left  = NULL;
  links->right = NULL;
  links->up    = NULL;
  links->n     = 1;
  links->sum   = 0;
  links->balanced =1;
}

/**********************************************************************
 * tree_top
 * Return the top node of a tree
 */
void *tree_top(void *tree, int link_offset) {
  long i;
  for(i=0; ((LINKS *)(tree+link_offset))->up && i<MAX_TREE_DEPTH; i++)
    tree= ((LINKS *)(tree+link_offset))->up;
  if(i==MAX_TREE_DEPTH) {
    fprintf(stderr,"Exceeded maximum tree depth!  Something's wrong here...\n");
    exit(-1);
  }
  return tree;
}

/**********************************************************************
 * tree_walk
 * Walk a tree, calling a function once per node in label order.
 * The function must return an int with 0=ok, nonzero=error.
 * It gets the node pointer and a depth meter.  It should know the
 * offsets!
 * Strong typing is for people with weak minds!
 * Returns 0 on normal completion, the failure value on error.
 *
 * tree_walk is a convenience function for tree_walker, which does the
 * work and also carts along a depth meter.
 */
long tree_walker(void *tree
		, int label_offset
		, int link_offset
		, long((*func)())
		, long depth){
  long err=0;

  LINKS *l = (LINKS *)(tree+link_offset);

  if(!tree)
    return 1;

  if(l->left) 
    err = tree_walker(l->left,label_offset, link_offset, func, depth+1);
  
  if(!err)
    err = (*func)(tree, label_offset, link_offset,depth);

  if(!err && l->right)
    err = tree_walker(l->right, label_offset, link_offset, func, depth+1);

  return err;
}
  
long tree_walk(void *tree
	       , int label_offset, int link_offset
	       , long ((*func)())){
  return tree_walker(tree,label_offset,link_offset,func,0);
}

/**********************************************************************
 * tree_find
 * Given a data type that includes tree links and an offset to the label
 * field and the tree link array, find the lebelled item in the tree,
 * or NULL.  
 *
 * Returns a pointer to the labelled field.
 */
void *tree_find(void *tree, long label, int label_offset, int link_offset){
  long node_label;
  LINKS *links;   

  if(tree==NULL) return NULL;
  node_label = *((long *)(tree+label_offset));
  if(node_label == label) return tree;

  links = (LINKS *)((void *)tree + link_offset);

  if( node_label > label) {
    if( links->left == NULL ) return NULL;
    return tree_find( links->left,
		      label,
		      label_offset, 
		      link_offset 
		      );
  } 

  if( links->right == NULL ) return NULL;
  return tree_find( links->right,
		    label,
		    label_offset,
		    link_offset
		    );
}


/**********************************************************************
 * tree_insert
 * Given a tree and a datum, insert the datum into the tree.  
 * Balancing is given proper thought.  The root of the tree is returned
 * (because it might change!)
 *
 * Because this is a single-node insertion only, the tree 
 * links from them item should both be NULL before insertion.  If that's
 * not true, the item is inserted anyway (possibly breaking whatever
 * tree it came from), and a warning message is printed.
 *
 * As an efficiency measure (in part because tree_balance uses tree_insert
 * to shuffle things around), tree_insert does not call tree_balance on 
 * branches that may become imbalanced.  It's best to call tree_b_insert
 * if that's what you want.
 */
void *tree_binsert(void *root, void *item, int label_offset, int link_offset) {
  return tree_balance(tree_insert(root, item, label_offset, link_offset),
		      label_offset, link_offset);
}
  
void *tree_insert(void *root, void *item, int label_offset, int link_offset) {
  long node_label;
  LINKS *item_links = (LINKS *)(item+link_offset);
  void *foo;
  long foo_label;

  if(!item) {
    fprintf(stderr,"HEY! tree_insert can't insert null items! I'm giving this tree back unchanged.\n");
    return root;
  }

  node_label = *((long *)(item+label_offset));
  
  /* Check assertion that this is a simple item and not a tree in its
   * own right */
  if( (item_links->left != NULL) ||
      (item_links->right != NULL) ||
      (item_links->up != NULL) ||
      (item_links->n > 1)) {
    fprintf(stderr,"HEY!  tree_insert got item #%ld, which isn't clean! I refuse to insert this.\n");
    return root;
  }
  
  /* Handle the no-root-at-all case */
  if(root == NULL) {
    (((LINKS *)(item+link_offset))->n)=1;
    return item;
  }

  /* Find the proper parent location for the new item */
  for( foo=root; foo; ) {
    long foo_label = *( (long *)(foo+label_offset) );
    LINKS *foo_links = (LINKS *)(foo+link_offset);

    if( foo_label < node_label ) {
      if(foo_links->right == NULL) {
	foo_links->right = item;
	break;
      } else {
	foo = foo_links->right;
      }
    }

    else if( foo_label > node_label) {
      if(foo_links->left == NULL) {
	foo_links->left = item;
	break;
      } else {
	foo = foo_links->left;
      } 
    }

    else if(foo == item) 
      break; /* Item is already in tree... */
    else {
      /* Item is not already in tree but label is not unique. */
      fprintf(stderr,"Assertion failed!  label %ld is not unique! Proceeding anyway!\n(Casey Jones, you're in trouble!)\n");
      if(foo_links->left == NULL) {
	foo_links->left = item;
	break;
      } else {
	foo = foo_links->left;
      }
    }
  }

  if(foo == NULL) {
    fprintf(stderr,"tree_insert should never get here.  You lose.\n");
    exit(2);
  }

  /* Link it up! (and notate the item as a leaf) */
  item_links->up = foo;
  item_links->n = 1;
  item_links->sum = *((NUM *)item);
  item_links->balanced = 1;

  /* Increment tree-size figures and set imbalance flags for the branch */
  for(; foo; foo = ((LINKS *)(foo+link_offset))->up) {
    LINKS *fool = (LINKS *)(foo+link_offset);
    (fool->n)++;
    (fool->sum) += *((NUM *)item);
    tree_balance_check(foo,link_offset);
  }

  return root;
}

/**********************************************************************
 * tree_balance_check
 * Given a node in a tree, check whether it is balanced (non-recursively)
 * using its node counts and its descendents' balanced flags.  Return
 * the balanced flag, but also put it into the node.  Allows 30% slop
 * in the size of the two branches before it complains.
 */
char tree_balance_check(void *foo, int link_offset) {
  LINKS *l;
  long ln, rn;
  char lb, rb;
  if(!foo) return 1;
  l = (LINKS *)(foo+link_offset);

  ln = l->left  ? ((LINKS *)(l->left  + link_offset))->n : 0;
  rn = l->right ? ((LINKS *)(l->right + link_offset))->n : 0;

  lb = l->left  ? ((LINKS *)(l->left  + link_offset))->balanced : 1;
  rb = l->right ? ((LINKS *)(l->right + link_offset))->balanced : 1;

  return 
    l->balanced = 
    (lb && rb) 
    ? ( ( ((ln - rn) <= 1) && ((rn - ln) <= 1) )
	|| ( fabs((float)ln - (float)rn)/((float)ln+(float)rn) < 0.3) )
    : 0;
}
  
/**********************************************************************
 * tree_unlink
 * Given an item in a tree, unlink the item from that tree.  The item
 * is not freed -- only unlinked from the tree!
 *
 * As always, you have to supply the offsets to the label and 
 * link fields of the tree data structure.
 * 
 * Returns the new root, which might chagne.
 */

void *tree_unlink(void *data, int label_offset, int link_offset) {
  LINKS *data_links = (LINKS *)(data+link_offset);
  void *foo, *root;

  /* If the data have no branches, then it's easy. */
  if(data_links->left == NULL && data_links->right == NULL) {
    for(foo = data_links->up; foo; foo = ((LINKS *)(foo+link_offset))->up) {
      if( (--(((LINKS *)(foo+link_offset))->n))  < 0) {
	fprintf(stderr,"You hoser -- negative node count found!\n");
      }
      ((LINKS *)(foo+link_offset))->sum -= *((NUM *)data);
      tree_balance_check(foo, link_offset);
      root = foo;
    }

    /* Unlink from parent, if possible. */
    if(data_links->up) {
      LINKS *foolinks = (LINKS *)(((void *)data_links->up)+link_offset);
      if(foolinks->left == data) {
	foolinks->left = NULL;
      } else if(foolinks->right == data) {
	foolinks->right = NULL;
      } else {
	fprintf(stderr,"tree_unlink:  parent doesn't point to this node!\n");
      }
    }

    clear_links(data+link_offset);
    ((LINKS *)(data + link_offset))->sum = *(NUM *)data;
    return root;
  }


  /* If the data have only one branch, then it's almost as easy. 
   * In this case, we remove the node and bump up the next node down.
   */
  if(data_links->left == NULL || data_links->right == NULL) {
   
    if(data_links->up) {
      LINKS *foolinks = (LINKS *)(((void *)data_links->up)+link_offset);
      if(foolinks->left == data) {
	foolinks->left = 
	  data_links->left ? 
	  data_links->left : 
	  data_links->right;
      } else if(foolinks->right == data) {
	foolinks->right =
	  data_links->left ? 
	  data_links->left :
	  data_links->right;
      } else {
	fprintf(stderr,"tree_delete: HEY!  The tree seems to have come unglued...\n");
      }
    }
    else {
      if(data_links->left == NULL) 
	root = data_links->right;
      else
	root = data_links->left;
    }

    /* Link its branches back up to its parent */
    if(data_links->left != NULL)
      ((LINKS *)(data_links->left + link_offset))->up = data_links->up;
    else 
      ((LINKS *)(data_links->right + link_offset))->up = data_links->up;


    for(foo = data_links->up; foo; foo = ((LINKS *)(foo+link_offset))->up) {
      (((LINKS *)(foo+link_offset))->n)--;
      ((LINKS *)(foo+link_offset))->sum -= *((NUM *)data);
      tree_balance_check(foo,link_offset);
      root = foo;
    }
    
    clear_links(data+link_offset);
    ((LINKS *)(data + link_offset))->sum = *(NUM *)data;
    return root;
  }
      
  /* If the data have two branches, then we promote the leftmost node
   * of the right hand branch to the top dog position.  That doesn't
   * necessarily preserve nice balance, but we can always fix that later
   * if necessary.
   */
  if( (data_links->left != NULL) && (data_links->right != NULL) ) {
    void *foo;

    for( foo = data_links->right; 
	 ((LINKS *)(foo + link_offset))->left; 
	 foo = ((LINKS *)(foo + link_offset))->left
	 )
      ;
    tree_unlink(foo,label_offset, link_offset);
    
    if(data_links->up == NULL) {
      root = foo;
    } else {
      LINKS *uplink = (LINKS *)(data_links->up + link_offset);

      /* Find root and also adjust the running sum to reflect the proper
       *  flux (the 'n' values are OK because all entries count for '1'.) 
       */
      for( root = data_links->up; 
	   ((LINKS *)(root + link_offset))->up;
	   (root= ((LINKS *)(root + link_offset))->up,
	    ((LINKS *)(root + link_offset))->sum += (*(NUM *)foo - *(NUM *)data))
	   )
	tree_balance_check(root,link_offset)
	  ;


      if (uplink->left == data)
	uplink->left = foo;
      else if(uplink->right == data)
	uplink->right = foo;
      else 
	fprintf(stderr,"tree_delete:  Wow!  The tree seems to have come unglued.\n");
      
      /* Link old neighbors to new place */
      ((LINKS *)(foo+link_offset))->left  = data_links->left;
      ((LINKS *)(foo+link_offset))->right = data_links->right;

      /* Backlink this place to the neighbors */
      if( data_links->left )
	((LINKS *)((LINKS *)(foo+link_offset))->left  + link_offset)->up = foo;
      if( data_links->right )
	((LINKS *)((LINKS *)(foo+link_offset))->right + link_offset)->up = foo;

	
      ((LINKS *)(foo+link_offset))->n 
	= ( (data_links->left ? ((LINKS *)(data_links->left + link_offset))->n : 0) + 
	    (data_links->right ? ((LINKS *)(data_links->right + link_offset))->n : 0) +
	    1 );

      ((LINKS *)(foo+link_offset))->sum
	= ( (data_links->left ? ((LINKS *)(data_links->left + link_offset))->sum : 0) +
	    (data_links->right? ((LINKS *)(data_links->right + link_offset))->sum : 0) +
	    *((NUM *)(foo)) );


      tree_balance_check(foo,link_offset);
      tree_balance_check(root,link_offset);

      clear_links(data+link_offset);
      ((LINKS *)(data + link_offset))->sum = *(NUM *)data;
      
      return root;
    }
  }
}

/**********************************************************************
 * tree_balance
 *
 * If there's an imbalance between left and right sides, 
 * unlink the rightmost <n> leaves from the left side or the leftmost
 * <n> leaves from the right side, and stick 'em on the opposite side.
 * Then rebalance both subtrees.  Return the new tree.
 *
 * This could be made more efficient by optimizing the tree_unlink
 * call a bit (and keeping it local):  rather than finding the rightmost
 * leaf every time, it's possible to "hold one's place".  
 * 
 * Note to self:  copy a real tree algorithm someday.  (not that this
 * bit of code will be the slowest in the bunch -- the "hot spot" is
 * in the geometry section.)  2-oct-2000
 *
 */
void *tree_balance(void *tree, int label_offset, int link_offset) {
  int i;
  void *l_tree, *r_tree;
  void *upnode;
  long l_n, r_n;
  LINKS *l_link, *r_link, *t_link;
  
  if(tree == NULL) 
    return tree;


  t_link = (LINKS *)(tree+link_offset);
  l_tree = t_link->left;
  r_tree = t_link->right;

  l_link = (LINKS *)(l_tree + link_offset);
  r_link = (LINKS *)(r_tree + link_offset);

  l_n = l_tree ? l_link->n : 0;
  r_n = r_tree ? r_link->n : 0;

  /** Some shuffling has to be done:  sever the head of the tree from the branches. **/
  if( ! t_link->balanced ) {
    char l_tweak=0;
    char r_tweak=0;

    if(l_tree) l_link->up = NULL;
    if(r_tree) r_link->up = NULL;
    
    t_link->left = NULL;
    t_link->right = NULL;
    t_link->up = NULL;
    t_link->n = 1;
    t_link->sum = *((NUM *)tree);
    upnode = t_link->up;

    
    /** Shuffle from left to right as necessary **/
    for(i=0; i < (l_n - r_n)/2 ; i++) {
      void *foo, *f;
      l_tweak=1;
      
      r_tree = tree_insert(r_tree, tree, label_offset, link_offset);


      for(foo = l_tree; f = ((LINKS *)(foo + link_offset))->right; foo = f) ;
      l_tree = tree_unlink(foo, label_offset, link_offset);

      tree = foo;
    }

    /** Shuffle from right to left as necessary **/
    for(i=0; i< (r_n - l_n)/2 ; i++) {
      void *foo, *f;
      r_tweak = 1;
      
      l_tree = tree_insert(l_tree, tree, label_offset, link_offset);
      
      for(foo = r_tree; f = ((LINKS *)(foo + link_offset))->left; foo = f) ;
      r_tree = tree_unlink(foo, label_offset, link_offset);

      tree = foo;
    }


    if(l_tweak && r_tweak) {
      fflush(stdout);
      fprintf(stderr,"balance_tree:  Got an impossible dual-tweak!\n");
      fflush(stderr);
    }
  }
   
  l_link = (LINKS *)(l_tree+link_offset);
  r_link = (LINKS *)(r_tree + link_offset);

  if(l_tree && (!l_link->balanced)) 
    l_tree = tree_balance(l_tree, label_offset, link_offset);

  if(r_tree && (!r_link->balanced))
    r_tree = tree_balance(r_tree, label_offset, link_offset);
  
  /** Mend the links **/
  t_link = (LINKS *)(tree+link_offset);
  t_link->up = NULL;
  t_link->left = l_tree;
  t_link->right = r_tree;
  t_link->n = 1;
  t_link->sum = *((NUM *)tree);
  t_link->balanced = 1;
  
  if(l_tree) {
    l_link = (LINKS *)(l_tree + link_offset);
    l_link->up = tree;
    t_link->n += l_link->n;
    t_link->sum += l_link->sum;
  }
  
  if(r_tree) {
    r_link = (LINKS *)(r_tree + link_offset);
    r_link->up = tree;
    t_link->n += r_link->n;
    t_link->sum += r_link->sum;
  }
  
  return tree;
  
}


/**********************************************************************
 * world_state_name
 * dumb -- just produces a string with the name of the state the world
 * (the given one, not the big one) is in.  
 */
const char *world_state_name(WORLD *a){
  if(a==NULL) {
    return "NULL WORLD";
  }
  switch(a->state) {
  case WORLD_STATE_NEW:     return "NEW";
  case WORLD_STATE_LOADING: return "LOADING";
  case WORLD_STATE_LOADED:  return "LOADED";
  case WORLD_STATE_WORKING: return "WORKING";
  case WORLD_STATE_RELAXED: return "RELAXED";
  case WORLD_STATE_READY:   return "READY";
  default: return "UNKNOWN";
  }
}

/**********************************************************************
 **********************************************************************
 ***** Dumblist routines
 ***** dumb lists keep track of vertices' neighbors.  They're pretty
 ***** dumb, but here's a generic implementation.  Each is an 
 ***** array of unsorted pointers to objects, and knows how many things
 ***** are in the list and how big the allocated space is.  
 *****
 ***** They grow to the power of 1.5 closest to the high-water mark and 
 ***** stay there till you get rid of 'em.  That's useful for short lists 
 ***** of things (say under 50) where you don't want lots of malloc 
 ***** overhead.  
 *****
 ***** For sanity reasons, these things crash if you try to make a dumblist
 ***** with more than 2^15 things in it.
 *****/
DUMBLIST *new_dumblist() {
  DUMBLIST *foo = (DUMBLIST *)localmalloc(sizeof(DUMBLIST),MALLOC_DL);
  foo->n = 0;
  foo->size = 0;
  foo->stuff = 0;
}
 
void free_dumblist(DUMBLIST *foo) {
  if(foo->stuff) 
    localfree (foo->stuff);
  localfree (foo);
}

void dumblist_quickadd(DUMBLIST *dl, void *a) {
  int i;
  void **foo;
  if(!dl || !a) return;
  if(!(dl->stuff)) {
    dl->stuff = (void **)localmalloc( 16*sizeof(void *),MALLOC_DL_L);
    dl->size = 16;
    dl->n = 1;
    dl->stuff[0] = a;
    return;
  }
  if(dl->n >= dl->size) 
    dumblist_grow(dl, dl->size +1);
  
  dl->stuff[dl->n] = a;
  dl->n++;
}
    
void dumblist_add(DUMBLIST *dl, void *a) {
  int i;
  void **foo;
  if(!dl || !a) return;

  /* If it's a null list, allocate some space (this could 
     actually fit in the extend-list conditions, but what the heck)
   */
  if(!(dl->stuff)) {
    dl->stuff = (void **)localmalloc( 16 * sizeof(void *),MALLOC_DL_L);
    dl->size = 16;
    dl->n = 1;
    dl->stuff[0] = a;
    return;
  }

  /* Traverse the list to see if a is already in it! */
  /* (this is what makes it dumb...) */
  for(foo=dl->stuff,i=0;i<dl->n && *foo != a;foo++,i++)
    ;
  if(i<dl->n) 
    return;
  
  /* OK, a isn't in the list -- add it on the end. */
  if(dl->n >= dl->size)
    dumblist_grow(dl,dl->size+1);

  dl->stuff[(dl->n)++] = a;
  return;
}

/**********************************************************************
 * dumblist_delete -- attempt to remove an item from a dumb list.
 * Finds each occurence of the item and blasts 'em.  If no such 
 * item is found, then nothing happens.  The variant,
 * sorted_dumblist_delete, knocks off as soon as a larger number than 
 * the requested one is found.  It could be reworked to run in 
 * log time rather than linear time, but who has the time?
 */

inline void dumblist_delete(DUMBLIST *dl, void *a) {
  int i;
  void **foo;

  if(!dl || !a) return;

  for(foo=dl->stuff,i=0; i<dl->n; foo++, i++) 
    if(*foo == a) {
      dumblist_rm(dl,i);
      i--;
    }
}
    
/**********************************************************************
 * dumblist_rm -- remove the ith item from a dumblist.  Move the last
 * item in the dumblist to the vanished position.
 */

inline void dumblist_rm(DUMBLIST *dl, int i) {
  if(i < dl->n) {
    dl->n--;
    if( i < dl->n )
      dl->stuff[i] = dl->stuff[dl->n];
  }
}


/**********************************************************************
 * dumblist_grow
 * Grow a dumblist to the required size.  Leave some margin.
 */
void dumblist_grow(DUMBLIST *dl, int size) {
  void **foo;
  void **a,**b;
  int i;
  long newsize;

  newsize = (3 * size)/2;

  if(newsize < 16)         /* Minimum size */
    newsize = 16;

  if(newsize <= dl->size)  /* No extra work please */
    return;
  
  if(newsize > 32768) {
    if(newsize > 1073741820) {
      fprintf(stderr,"dumblist_grow:  size is too big!\n");
      exit(-30);
    }
    else {
      int a,b=0;
      fprintf(stderr,"dumblist_grow: Warning -- size is %d!\n",newsize);
      //j      a /= b; // throw an arithmetic exception
    }
  }
  
  /* This seems to corrupt the arena on my linux box -- replaced 
   * with more pedestrian version below.  --CED 07-Oct-2004
   *
   * dl->stuff = (void **)realloc(dl->stuff,newsize * sizeof(void *));
   */
  {
    void **oldstuff, **ost, **newstuff;
    int i;
    oldstuff = ost = dl->stuff;
    dl->stuff = newstuff = (void **)localmalloc(newsize * sizeof(void *),MALLOC_DL_L);
    if(ost) {
      for(i=0;i<dl->size; i++)
      *(newstuff++) = *(oldstuff++);
      localfree(ost);
    }
  }
    
    dl->size = newsize;
}
    

static void **dls_wk = 0;
static long dls_size = 0;         /* Size of dumblist we can currently qsort */
static long dls_sz = 0;           /* Size of allocated space for qsorting */
static int ((*dls_cmp)(void *a, void *b));

/* Sort the given dumblist -- just the n elements starting at l --
 * using the workspace starting at element wk_start.  Return the 
 * start of the sorted list in the workspace.  Eat up workspace linearly --
 * there should be plenty (never call this routine from outside; call the 
 * interface routine dumblist_sort instead!)
 * 
 * This is just your basic qsort routine -- the main reasons to 
 * reimplement it here rather than using the libraries are:
 *    (A) it's kind of fun to mess with qsort after all these years;
 *    (B) I wanted one that would munch up a known workspace with
 *        known-max-size lists, to avoid lots of mallocation.
 * 
 * This routine is doomed to be part of the 'hot spot' for the model --
 * it should probably be retweaked for maximum speed, or rewritten 
 * to omit the internal function call for each compare (pre-process the
 * list to use floating point size indices?)
 * 
 */
void **dumblist_qsort(void **l, int n, void **wk_start) {
  if(n<=1) {
    if(!n)        /* Paranoia */
      return l;

    *wk_start = *l;  /* n=1: return the input list (it's sorted!) */
    return wk_start;
  }

  {  /* Block avoids variable initialization for the trivial case */
    int nl=n/2;
    int nr=n-nl;
    int n_copied; 
    void **out,**out1;
    void **left, **left0;
    void **right, **right0;
    
    /* Divide the list in half (nl and nr), and sort each half. */
    if(n>2) {
      left0  = left  = dumblist_qsort(l    , nl, wk_start      );
      right0 = right = dumblist_qsort(l+nl , nr, left + nl + 1 );
    } else {  /* n==2 */
      /* No need to iterate to the hairy end for 2 trivial lists. */
      left0 = left = wk_start;
      right0 = right = left+1;
      *left = *l;
      *right = l[1];
    }
    
    out1 = out = right + nr + 1;
    
    /* Merge the lists, skipping duplicates and copying back into place. */
    for(n_copied = 0; 
	n_copied < n;
	n_copied++) {
      int l_done = (left  - left0 >= nl)   || (*left == NULL)  ;
      int r_done = (right - right0 >= nr)  || (*right == NULL) ;
      int foo;

      if(!l_done) {
	if(!r_done) {
	  switch (foo=((*dls_cmp)(*left,*right))) {
	    
	  case -1:
	    *(out++) = *(left++);
	    break;
	    
	  case 1:
	    *(out++) = *(right++);
	    break;
	    
	  case 0:
	    *(out++) = *(right++);
	    left++;
	    break;
	    
	  default:
	    fprintf(stderr,"Unknown response %d from sort function in dumblist_qsort.  Giving up.\n",foo);
	    exit(-10);
	    break;
	    
	  }

	} else  /* r_done && !l_done */
	  *(out++) = *(left++);
      }
      else if(!r_done) /* !r_done && l_done */
	*(out++) = *(right++);
      else /* l_done and r_done are true */ {
	break; /* Don't waste any more iterations! */
      }
    }

    *out = 0; /* Lay down a fence */

    return out1;
  } /* end of convenience block */
} /* end of dumblist_qsort */


/******************************
 * dumblist_shellsort is a basic comb/shell sorter.
 * It's useful for lists too small to need quicksort (up to, say, 50 elements).
 * @@@ DIRTY - FIXME
 */
void dumblist_shellsort( DUMBLIST *dl, int ((*cmp)(void *a, void *b)) ) {
  int i,j,increment;
  void *temp;
  int n = dl->n;
  int iter=0;

#ifdef testing_with_bubblesort
  // Bubblesort to make sure shellsort isn't scrozzling memory
  char done;
  int ret;
  do {
    done = 1;
    for(i=0;i<n-1;i++) {
      if(  (ret = ((*cmp)(dl->stuff[i],dl->stuff[i+1]))) > 0) {
	void *foo = dl->stuff[i+1];
	dl->stuff[i+1] = dl->stuff[i];
	dl->stuff[i] = foo;
	done = 0;
      }
      if(((*cmp)(dl->stuff[n-2-i],dl->stuff[n-1-i]))>0) {
	void *foo = dl->stuff[n-2-i];
	dl->stuff[n-2-i] = dl->stuff[n-1-i];
	dl->stuff[n-1-i] = foo;
	done = 0;
      }
    }
  } while(!done);
#endif

  increment = (dl->n / 3) || 1;

  while(increment>0)
    { 
      for(i=0; i < dl->n; i++) {

	j = i;
	temp = dl->stuff[i];
	while((j >= increment) && ( (*cmp)(dl->stuff[ j - increment ],temp) > 0) ) {
	  dl->stuff[j] = dl->stuff[j-increment];
	  j = j - increment;
	}
	dl->stuff[j] = temp;
      }
      increment >>= 1;
    }
}

/******************************
 * dumblist_crunch crunches a sorted dumblist, 
 * removing duplicate elements.
 */
void dumblist_crunch(DUMBLIST *dl,int((*cmp)(void *a, void *b))) {
  void **a, **b;
  int i,j;

  b = &(dl->stuff[0]);
  a = &(dl->stuff[1]);
  j =0;
  for(i=1;i<dl->n;i++) {
    if((*cmp)(*b,*a)) {
      b++; j++;
      if(b != a) 
	*b = *a;
    }
    a++;
  }
  dl->n = j+1;
}

/**********************************************************************
 * dumblist_sort
 * Sort a dumblist by some external criterion.  You feed in 
 * a function that returns -1, 0, or 1 depending on whether the
 * first argument is less than, equal to, or greater than the 
 * second, and the list is sorted for you. 
 *
 * The sorting takes place in a workspace that is static, so that
 * it isn't always being remalloc'ed. 
 *
 * See notes in dumblist_qsort -- this is doomed to be a hot-spot
 * routine, but isn't particularly well optimized.
 *
 * If the number of elements is small, then shellsort is used as
 * it's faster than qsort for the small stuff.
 *
 */

void dumblist_sort(DUMBLIST *dl, int ((*cmp)(void *a, void *b))) {
  void ** sorted;
  void **dl_a;
  int i;

  /* use shellsort all the time for now */
  if( dl->n <= 100 ) {             /* use shellsort for the small cases */
    dumblist_shellsort(dl,cmp);
    dumblist_crunch(dl,cmp);
  } else {
    int odls_size;
    /* Make sure that there's enough scratch space */
    if(dls_size < dl->size) {
      long levels = 2;
      unsigned long a;
      /* This is kind of wasteful but should be tiny compared to the
	 main data structure -- so don't be makin' global gather_neighbor
	 calls if you can possibly avoid it.  

	 Reserve enough space for two copies at each level, and enough
	 levels to cover the whole depth.

      */
      odls_size = dls_size;
      dls_size = a = dl->size * 2;
      while(a >>= 1)   /* shift and test for nonzero */
	levels++;
      dls_sz = (2*dls_size)*levels;
      
      /* realloc seems to bust the arena -- so manually reallocate...
       * dls_wk = (void **)realloc(dls_wk, dls_sz * sizeof(void *));
       */
      if(dls_wk)
	localfree(dls_wk);
      dls_wk = (void **)localmalloc(dls_sz * sizeof(void *),MALLOC_DLS_ARENA);
    }
    
    /* Launch the sort */
    dls_cmp = cmp;
    sorted = dumblist_qsort(dl->stuff,dl->n,dls_wk);
    
    /* Copy the sorted list back into the vertex. */
    dl_a = dl->stuff;
    for(i=0;i<dl->n && *sorted;i++) {
      *(dl_a++) = *(sorted++);
    }
    dl->n = i;
  }
}
  

int ptr_cmp(void *a, void *b) { /* Helper routine for sorting by pointer */
  if(a<b) return -1;
  else if(a>b) return 1;
  return 0;
}

/**********************************************************************
 * dumblist_snarf 
 * Copy every item from the source dumblist into the destination.
 * They're quicksorted to avoid duplication.  This is a bit slow,
 * because in the normal course of things the model snarfs up vertices
 * and *then* sorts 'em by angle.  But in pathological cases this 
 * doesn't eat as much memory as not sorting each time.
 */
void dumblist_snarf(DUMBLIST *dest, DUMBLIST *source) {
  int i;
  void **b;
  void **b1;
  void **c;
  
  if( !dest || !source || dest==source) {
    fprintf(stderr,"dumblist_snarf: Got null src or dest, or dest==src! (d=%d, s=%d)\n",dest,source);
    fflush(stderr);
    return;
  }
  if(dest->size <= dest->n + source->n) 
    dumblist_grow(dest,dest->n + source->n);
  
  for(b=source->stuff, c=&(dest->stuff[dest->n]), i=0; i<source->n; i++) {
    *(c++) = *(b++);
  }
  
  dest->n += source->n;

}

/**********************************************************************
 * dumblist_clear
 * Zeroes out a dumblist.  Especially useful for the dumblists in VERTEXes,
 * because they're included in the base struct instead of being pointers out.
 */
inline void dumblist_clear(DUMBLIST *foo) {
  foo->n = 0;
  foo->size = 0;
  foo->stuff= 0;
}

void dd_vertex_printer(void *a) {
  printf("vertex #%d",((VERTEX *)a)->label);
}

inline void dumblist_dump(DUMBLIST *foo, void ((*printer)(void *a))) {
  int i;

  if(!foo) {
    printf("dumblist_dump got a null list!\n");
  } else {
      printf("list has %D elements\n");
      for (i=0;i<foo->n;i++) {
	printf("\t%d: ",i);
	if(printer)
	  (*printer)(foo->stuff[i]);
        else 
	  printf("element ptr=%d",foo->stuff[i]);
	printf("\n");
      }
  }
}
      
      

/**********************************************************************
 * fl_eq
 * Test whether two NUMs are reasonably close.  
 * Currently, 100ppm is the slop.
 */
inline char fl_eq(NUM a, NUM b) {
  return ( (a==0 && b==0) ||
	   ( !(a==0 || b==0) &&
	     (a/b > 0) &&
	     ( (a/b < 1.0001) && (b/a < 1.0001) )
	     )
	   );
}
  

/**********************************************************************
 *
 * flux_malloc
 *
 * A wrapper for malloc that includes optional fencing.
 * Only gets compiled if you asked for debugging malloc with
 * -DUSE_DEBUGGING_MALLOC at compile time.
 *
 * flux_malloc allocates a fence of size fencecount*sizeof(long) 
 * on either side of each newly allocated bit of memory, then fills
 * the fences with 0xDEADBEEF, and monitors them on every memory operation.
 * If the fences get touched, a warning gets printed.
 * 
 * Of course, there is a rather severe hit on performance! 
 * 
 * With small fences (4 longs), memory contention issues still
 * appear in the perl side -- notably, the malloc arena gets corrupted
 * and the code crashes.  But none of the fences ever sustains a hit.
 * With large fences (100 longs), the crashes go away -- but none
 * of the fences ever seem to get touched. 
 * 
 */
#ifdef USE_DEBUGGING_MALLOC

static char *flux_malloc_types[MALLOC_MAXTYPENO] = 
  {"World","Fluxon","Vertex","Flux Concentration",
   "Dumblist","Dumblist stuff field","Vertex list","DLS arena",
   "Photospheric plane", "Misc"
  };

#define MALLOC_MAXNO (10000)
static struct flux_malloc_alloc {
    void *where;
    long size;
    int type;
  } flux_malloc_alloc[MALLOC_MAXNO];
  static int flux_malloc_next = 0;

#define fencecount (100)
#define fencesize  (4*fencecount)

char *flux_malloc(long size, int what_for) {
  char *p;
  long *p2;

  flux_memcheck();
  
  if(flux_malloc_next >= MALLOC_MAXNO){
    fprintf(stderr,"Filled malloc table; ending test\n");
    exit(1);
  }

  if(!valid_malloc_type(what_for)) {
    fprintf(stderr,"Invalid type received by flux_malloc\n");
    // Throw an exception
    {
    int i = 0;
    int j = 2;
    j /= i;
    }
    return (void *)0; 
  }
  
  // call malloc and pad for the fence.
  p = (char *)malloc(size + 2*fencesize);
  p2 = (long *)p;
  { int i;
  for(i=0;i<fencecount;i++) *(p2++) = 0xDEADBEEF;
  p = (char *)p2;
  p2 = (long *)(p+size);
  for(i=0;i<fencecount;i++) *(p2++) = 0xDEADBEEF;
  }

  flux_malloc_alloc[flux_malloc_next].where = p;
  flux_malloc_alloc[flux_malloc_next].size = size;
  flux_malloc_alloc[flux_malloc_next].type = what_for;
  flux_malloc_next++;
  return (void *)p;
}
  
void flux_free(void *p) {
  int i;
  int ok = 0;
  flux_memcheck();

  for(i=0;!ok && i<flux_malloc_next;i++) 
    if (p==flux_malloc_alloc[i].where) 
      ok = i ;
  
  if(ok) {
    flux_malloc_alloc[ok].where = 0;
    flux_malloc_alloc[ok].size = 0;
    flux_malloc_alloc[ok].type = 0;
    free((char *)p-fencesize);
    return;
  } else {
    fprintf(stderr,"%%%%%%flux_free: got an unknown block to free");
    // Throw an exception
    {
      int i=0,j=2;
      j/=i;
    }
    free((char *)p-fencesize);
    return;
  }    
}

void flux_memcheck() {
  long badct = 0;
  long freed = 0;
  int i;
  char flag =0;
  static FILE *foo=0;

  if(!foo) {
    foo = fopen("foo.txt","w");
  }
  fprintf(foo,"%c======memcheck table:\n",(char)12);
  for(i=0;i<flux_malloc_next;i++) {
    char type_bad=0;
    char fence_bad=0;
    char *p;
    long *p2;
    fprintf(foo,"%4.4d: 0x%8.8x size=%6.6d (%4.4x) type=%d (%s)\n",
	   i,flux_malloc_alloc[i].where,
	   flux_malloc_alloc[i].size,flux_malloc_alloc[i].size,
	    flux_malloc_alloc[i].type,
	    malloc_types[flux_malloc_alloc[i].type]
	   );
    if(!flux_malloc_alloc[i].where) {
      freed++;
    } else {
      int ii;
      p = flux_malloc_alloc[i].where;
      p2 = (long *)(p-fencesize);
      
      for(ii=0; ii<fencecount;ii++) 
	fence_bad += (*(p2++) != 0xDEADBEEF);
      p2 = (long *)(p+flux_malloc_alloc[i].size);
      if(!valid_malloc_type(flux_malloc_alloc[i].type)) 
	type_bad = 1;
      for(ii=0; ii<fencecount; ii++)
	fence_bad += (*(p2++) != 0xDEADBEEF);
    }
    if(fence_bad) {
      int ii;
      char buf[fencecount*2+10];
      char *s=buf;
      // Announce the problem
      fprintf(stderr,"%%%%%%flux_memcheck: slot %d, 0x%x (%s) has a bad fence! (%d bad longs)\n",i,p,type_bad?"UNKNOWN_TYPE":flux_malloc_types[flux_malloc_alloc[i].type],fence_bad);
      for(ii=0;ii<fencecount;ii++) {
	*(s++) = ( ((long *)(flux_malloc_alloc[i].where - fencesize))[ii]==0xDEADBEEF ) ? 'X' : '.';
      }
      *(s++) = ' '; *(s++) = '|'; *(s++) = ' ';
      for(ii=0;ii<fencecount;ii++) {
	*(s++) = ( ((long *)(flux_malloc_alloc[i].where + flux_malloc_alloc[i].size))[ii]==0xDEADBEEF ) ? 'X' : '.';
      }
      fprintf(stderr,"  MAP: %s\n",buf);
      // Should mend the fence here...
    }
    
    if(fence_bad || type_bad)
      badct++;
  }
  fprintf(foo, "%sflux_memcheck (both-sides fencing): %d table entries (%d freed) OK, %d good, %d bad\n",(badct==0?"   ":"%%%"),i,freed,i-freed-badct,badct);
  fflush(foo);
}

char *malloc_types[MALLOC_MAXTYPENO+1] = {
  "Nothing",
  "World","Fluxon","Vertex","Flux Conc.",
  "Dumblist","Dumblist_stuff","Vertex List","DLS Arena",
  "Plane"};

 /* end of debugging-malloc code */

#else
#ifdef USE_PERL_MALLOC

/******************************
 * perl malloc uses the perl allocator; this was an idea 
 * to avoid memory contention issues that cropped up with the 
 * perl interface.  It doesn't help.
 * --CED 27-Oct-2004
 */
char *flux_perl_malloc(long size) {
  char *out;
  New(0,out,size,char);
  if(!out) {
    int i,j;
    fprintf(stderr,"Out of memory!\n");
    i=2; j=0;
    i /= j; /* throw arithmetic exception */
    exit(3498);
  }
  printf("Alloc'ed %d bytes (%x)\n",size,out);
  fflush(stdout);
  return out;
}

void flux_perl_free(void *foo) {
  //  Safefree(foo);   /* do nothing for now */
}

#else
#ifdef USE_PADDED_MALLOC

/********** 
 * Padded malloc allocates more space but doesn't keep track of it.  That
 * seems to make the perl allocator happy, I'm not sure why. (allocation pattern
 * is similar to debugging malloc, but we don't bother with keeping track of the
 * fencing).  
 *  --CED 27-Oct-2004
 */
#define padding (500)
char *flux_padded_malloc(long size) {
  return malloc(size+2*padding)+padding;
}

void flux_padded_free( void *foo ) {
  free(foo-padding);
}
#endif  
#endif
#endif

  

    
