#ifndef STATE_H
#define STATE_H
/* contains some hacky global variables to take care of the scope of the data used in the 2p->1s program */

extern int ns;
extern int nt;

/* number of states in each energy range provided in the input */
extern int n_is;
extern int n_tmax;/* maximum number of intermediate state transitions */
extern int n_gfs;

extern double tmax_d;
extern double tmax_q;
extern double e0;

extern double ** parsed_input;
extern int * idx_map;
extern int sz_inp;
#endif /* STATE_H */
