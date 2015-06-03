#ifndef STATE_H
#define STATE_H
/* contains some hacky global variables to take care of the scope of the data used in the 2p->1s program */

extern int ns;
extern int nt;
extern n_is, n_fs;

extern double tmax_d;
extern double tmax_q;
extern double e0;

extern int * idxs_map;
extern double ** parsed_input;

#endif /* STATE_H */
