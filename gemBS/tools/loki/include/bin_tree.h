#ifndef _BIN_TREE_H_
#define _BIN_TREE_H_

struct bin_node {
	struct bin_node *left,*right;
	void *data;
	int balance;
};

struct bin_node *rotate_left(struct bin_node *);
struct bin_node *rotate_right(struct bin_node *);
struct bin_node *insert_bin_node(struct bin_node *,const void *,struct bin_node **,
		 int (*)(const void *,const void *),struct bin_node *(*)(const void *));
struct bin_node *find_bin_node(struct bin_node *,const void *,int (*)(const void *,const void *));
struct bin_node *alloc_bin_node(const void *);
void free_bin_tree(struct bin_node *,void (*)(void *));
void traverse_bin_tree(struct bin_node *,void (*)(struct bin_node *));
void traverse_bin_tree1(struct bin_node *,void (*)(struct bin_node *,void *),void *);

#endif
