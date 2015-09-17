#ifndef SPEC_INFO_H
#define SPEC_INFO_H
#include "structs.h"

extern spec_info root_sinfo;
extern int n_sinfo;

/* suffixes for the output and input files */
extern const char * dat_sfx;

extern const char * dat_sfx;
extern const char * plot_sfx;
extern const char * log_sfx;
extern const char * bin_sfx;
extern const char * tmp_sfx;
extern const char * log_sfx;

metadata
init_md ();

int
free_md (metadata md);

spec_info
get_sinfo (char * idx);

spec_info
init_sinfo (metadata md);

#endif /* SPEC_INFO_H */
