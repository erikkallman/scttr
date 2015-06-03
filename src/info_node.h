#ifndef INFO_NODE_H
#define INFO_NODE_H
#include "structs.h"

extern info_node root_inode;

info_node
get_inode (char * idx);

info_node
init_inode (char * s,
            int * mom,
            int ns,
            int nt
                );

/* function set_symtrans

   * synopsis:
   set_symtrans finds all symmetric transitions between initial and
   intermediate states, and adds them to their corresponding
   intermediate states.

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
set_symtrans (info_node inode_root);


#endif /* INFO_NODE_H */
