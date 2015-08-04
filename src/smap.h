#ifndef SMAP_H
#define SMAP_H

/* function calc_smap

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
void
calc_smap_m (char * fn_infile,
             char * dat_fpstr,
             char * plot_fpstr,
             double * state_er,
             double * state_t,
             double * res,
             double * fwhm_inp
             );

void
write_log (double * state_er,
           double * state_t,
           double * res,
           double * fwhm_inp,
           char * fn_relpath,
           char * log_fpstr,
           int n_max
           );

#endif /* SMAP_H */
