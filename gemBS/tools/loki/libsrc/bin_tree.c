/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - CNG, Evry                                *
 *                                                                          *
 *                       October 2002                                       *
 *                                                                          *
 * bin_tree.c:                                                              *
 *                                                                          *
 * Routines for binary tree stuff                                           *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>

#include "lk_malloc.h"
#include "bin_tree.h"

struct bin_node *rotate_left(struct bin_node *node)
{
  struct bin_node *l,*gc;
	
  l=node->left;
  if(l->balance==-1) {
    node->left=l->right;
    l->right=node;
    node->balance=l->balance=0;
    node=l;
  } else {
    gc=l->right;
    l->right=gc->left;
    gc->left=l;
    node->left=gc->right;
    gc->right=node;
    switch(gc->balance) {
    case -1:
      node->balance=1;
      l->balance=0;
      break;
    case 0:
      node->balance=l->balance=0;
      break;
    case 1:
      node->balance=0;
      l->balance=-1;
    }
    gc->balance=0;
    node=gc;
  }
  return node;
}

struct bin_node *rotate_right(struct bin_node *node)
{
  struct bin_node *r,*gc;
	
  r=node->right;
  if(r->balance==1) {
    node->right=r->left;
    r->left=node;
    node->balance=r->balance=0;
    node=r;
  } else {
    gc=r->left;
    r->left=gc->right;
    gc->right=r;
    node->right=gc->left;
    gc->left=node;
    switch(gc->balance) {
    case -1:
      node->balance=0;
      r->balance=1;
      break;
    case 0:
      node->balance=r->balance=0;
      break;
    case 1:
      node->balance=-1;
      r->balance=0;
    }
    gc->balance=0;
    node=gc;
  }
  return node;
}

static struct bin_node *_insert_bin_node(struct bin_node *node,const void *s,struct bin_node **target,int *balanced,
					 int (*compar)(const void *,const void *),struct bin_node *(*alloc)(const void *))
{
  int i;
	
  i=compar(node->data,s);
  if(i<0) {
    if(node->left) node->left=_insert_bin_node(node->left,s,target,balanced,compar,alloc);
    else {
      node->left=alloc(s);
      if(target) *target=node->left;
      *balanced=0;
    }
    if(!(*balanced)) {
      switch(node->balance) {
      case -1:
	node=rotate_left(node);
	*balanced=1;
	break;
      case 0:
	node->balance=-1;
	break;
      case 1:
	node->balance=0;
	*balanced=1;
      }
    }
  } else if(i>0) {
    if(node->right) node->right=_insert_bin_node(node->right,s,target,balanced,compar,alloc);
    else {
      node->right=alloc(s);
      if(target) *target=node->right;
      *balanced=0;
    }
    if(!(*balanced)) {
      switch(node->balance) {
      case -1:
	node->balance=0;
	*balanced=1;
	break;
      case 0:
	node->balance=1;
	break;
      case 1:
	node=rotate_right(node);
	*balanced=1;
      }
    }
  } else {
    if(target) *target=node;
    *balanced=1;
  }
  return node;
}

struct bin_node *alloc_bin_node(const void *s)
{
  struct bin_node *nd;
	
  nd=lk_malloc(sizeof(struct bin_node));
  nd->left=nd->right=0;
  nd->balance=0;
  nd->data=(void *)s;
  return nd;
}

struct bin_node *insert_bin_node(struct bin_node *node,const void *s,struct bin_node **target,
				 int (*compar)(const void *,const void *),struct bin_node *(*alloc)(const void *))
{
  struct bin_node *nd;
  int balanced;
	
  if(node) nd=_insert_bin_node(node,s,target,&balanced,compar,alloc);
  else {
    nd=alloc(s);
    if(target) *target=nd;
  }
  return nd;
}

struct bin_node *find_bin_node(struct bin_node *node,const void *s,int (*compar)(const void *,const void *))
{
  int i;
	
  i=compar(node->data,s);
  if(i<0) {
    if(node->left) return find_bin_node(node->left,s,compar);
    return 0;
  } else if(i>0) {
    if(node->right) return find_bin_node(node->right,s,compar);
    return 0;
  }
  return node;
}

void traverse_bin_tree(struct bin_node *node,void (*op)(struct bin_node *))
{
  struct bin_node *node1;
  if(node->left) traverse_bin_tree(node->left,op);
  node1=node->right;
  op(node);
  if(node1) traverse_bin_tree(node1,op);
}

void traverse_bin_tree1(struct bin_node *node,void (*op)(struct bin_node *,void *data),void *data)
{
  if(node->left) traverse_bin_tree1(node->left,op,data);
  op(node,data);
  if(node->right) traverse_bin_tree1(node->right,op,data);
}

void free_bin_tree(struct bin_node *node,void (*free_node)(void *))
{
  if(node->left) free_bin_tree(node->left,free_node);
  if(node->right) free_bin_tree(node->right,free_node);
  if(free_node && node->data) free_node(node->data);
  free(node);
}
