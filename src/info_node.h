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

#endif /* INFO_NODE_H */
