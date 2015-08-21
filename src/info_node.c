#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "appc.h"
#include "info_node.h"
#include "structs.h"
#include "estate.h"

info_node
get_inode (char * id
           ){

  info_node curr_inode = root_inode;
  info_node next_inode;

    /* locate the info_node corresponding to the file name input argument */
  while(strstr((curr_inode -> str_id),id) == NULL) {

    next_inode = curr_inode -> next;
    curr_inode = next_inode;

    if (curr_inode == NULL) {
      fprintf(stderr, "parse_input.c: get_inode, no info node can be found\
 with a str_id == %s\n", id);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  return curr_inode;
}

info_node
init_inode (char * s,
            int ns,
            int nt
            ){

  int j,k;

  static info_node last_inode;
  int str_sz = strlen(s);

  int * el_idxs;
  char * si = malloc(strlen(s)+1); /* info node id string */

  info_node new_inode;

  if((new_inode = malloc(sizeof(struct info_node_s))) == NULL ){
    fprintf(stderr, "info_node.c:function init_inode, malloc: failed \
to allocate memory for \"new_inode\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((el_idxs = malloc((ns*2)*sizeof(int))) == NULL ){
    fprintf(stderr, "info_node.c:function init_inode, malloc: failed \
to allocate memory for \"el_idxs\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }


  for (j=0; j<str_sz; j++) {
    si[j] = s[j];
  }

  new_inode -> str_id = si;

  if (n_inodes == 0) { /* there is no root info node defined  */
    n_inodes = 1;
    new_inode -> last = NULL;
    root_inode = new_inode;
    last_inode = root_inode;
  }
  else { /* the root info node is already defined */
    n_inodes++;
    new_inode -> last = last_inode;
    last_inode -> next = new_inode;
  }
  new_inode -> next = NULL;
  new_inode -> idx = n_inodes-1;

  return new_inode;
}
