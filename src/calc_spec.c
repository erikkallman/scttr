#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <complex.h>
#include "calc_spec.h"
#include "dyn_array.h"
#include "parse_input.h" /* struct types */
#include "sci_const.h"
#include "std_num_ops.h" /* power */
#include "std_char_ops.h" /* power */
#include "spectrum_info.h"

void
calc_spec (spec_info s
             ) {

  metadata md = s -> md;

  char * dat_fpstr = concs(3,md->outpath,md->inp_fn,dat_sfx);
  char * plot_fpstr = concs(3,md->outpath,md-> inp_fn,plot_sfx);

  FILE * fp;
  FILE * fp_plot_in;
  FILE * fp_plot_out;

  size_t len;
  ssize_t read;

  char * line;

  unsigned long progress;

  unsigned long ksum;

  int j,k; /* control loop indices */

  /* index variables for looping over the rixs map */
  int max_x;
  int max_y;
  int x,y;

  int pstep; /* determines at what values progression is printed */

  int is_num;
  int is_idx;

  /* create local copies and pointers to variables whos name would
   otherwise clutter the rixsmap loop */
  int n_sfs = s->scr->is2fs->n_el;
  int is_pos;
  int * is2fs = s -> scr -> is2fs -> a;
  int * gs2is = s -> scr -> gs2is -> a;
  int * is_idxs = s -> scr -> is_idxs -> a;
  int * ii_start = s -> scr -> ii_start -> a;

  double rmax = -0.1;

  double de_gi,de_if;

  double o_x, o_y;
  double bw;

  /* energies and transition moments along the axis of .. */
  /* .. excitiation (x-axis) */
  double ediff_x,tmom_gi;

  /* .. energy transfer (y-axis)*/
  double ediff_y,tmom_if;

  double emin_x = s -> scr -> emin_x;
  double emax_x = s -> scr -> emax_x;

  double emin_y = s -> scr -> emin_y;
  double emax_y = s -> scr -> emax_y;

  double de_x,de_y,fwhm_in,fwhm_tr;

  /* variables used in the Kramers-Heisenberg formula */
  double complex tmp;

  /* excitation energy (x-axis) broadening parameters */
  double grms_in; /* gaussian RMS value */
  double gvar_in; /* square root of the variance of the gaussian */

  /* transfer energy (y-axis) broadening parameters */
  double grms_tr;
  double gvar_tr;

  /* get a pointer to the 0th column of iis */
  double ** omega_x;
  double ** omega_y;
  double ** rixsmap;

  /* open the placeholder file */
  if((fp=fopen(dat_fpstr, "w"))==NULL) {
        fprintf(stderr,"calc_spec.c, function calc_spec: unable to open the output file %s.\n",dat_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  /* open the placeholder file */
  if((fp_plot_in=fopen("../src/plot_template", "r"))==NULL) {
    fprintf(stderr,"calc_spec.c, function calc_spec: unable to open the output file %s.\n","../src/plot_template");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((fp_plot_out=fopen(plot_fpstr, "w"))==NULL) {
    fprintf(stderr,"calc_spec.c, function calc_spec: unable to open the output file %s.\n",plot_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  fwhm_in = md->fwhm[0]/AUTOEV;
  fwhm_tr = md->fwhm[1]/AUTOEV;

  de_x = md->res[0]/AUTOEV;
  de_y = md->res[1]/AUTOEV;

  max_x = (int)((emax_x-emin_x)/de_x);
  max_y = (int)((emax_y-emin_y)/de_y);

  if((omega_x = malloc(max_x*sizeof(double*))) == NULL ){
    fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed \
to allocate memory for \"omega_x\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((omega_y = malloc(max_x*sizeof(double*))) == NULL ){
    fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed \
to allocate memory for \"omega_y\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((rixsmap = malloc(max_x*sizeof(double*))) == NULL ){
    fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed \
to allocate memory for \"rixsmap\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<max_x; j++) {
    if((rixsmap[j] = malloc(max_y*sizeof(double))) == NULL ){
      fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed \
to allocate memory for \"rixsmap[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }

    if((omega_x[j] = malloc(max_y*sizeof(double))) == NULL ){
      fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed \
to allocate memory for \"omega_x[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }

    if((omega_y[j] = malloc(max_y*sizeof(double))) == NULL ){
      fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed \
to allocate memory for \"omega_y[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  progress = max_x*max_y;
  pstep = progress/1000;

  grms_in = 2.0*powerl((fwhm_in/(2*sqrt(2*log(2)))),2);
  gvar_in = fwhm_in/(2.0*sqrt(2.0*log(2)))*sqrt(2.0*3.1415927);

  grms_tr = 2.0*powerl((fwhm_tr/(2*sqrt(2*log(2)))),2);
  gvar_tr = fwhm_tr/(2.0*sqrt(2.0*log(2)))*sqrt(2.0*3.1415927);

  printf( "  - calculating the RIXS map .. \n");

  ksum = 0;

  for (x=0; x<max_x; x++) {
    o_x = emin_x+(x*de_x);

    for (y=0; y<max_y; y++) {
      ksum++;
      rixsmap[x][y] = 0;
      o_y = emin_y+(y*de_y);
      omega_x[x][y] = o_x;
      omega_y[x][y] = o_y;

      k = 0;

      for (j=0; j<n_sfs; j++) {
        is_pos = ii_start[j];
        tmom_if = s->trs[4][is2fs[j]];
        de_if = s->trs[3][is2fs[j]] - s->trs[2][is2fs[j]];

        is_idx = is_idxs[j];
        is_num = (int)s->trs[0][is2fs[j]];
        tmp = 0 + 0*I;

        for (k=is_pos; ((is_idx = is_idxs[k]) != -1); k++) {

          tmom_gi = s->trs[4][gs2is[is_idx]];
          de_gi = s->trs[3][gs2is[is_idx]] - s->trs[2][gs2is[is_idx]];

          bw = get_rbdist(s->e0,s->trs[2][gs2is[is_idx]]);

          ediff_x = omega_x[x][y] - de_gi;
          ediff_y = omega_y[x][y] - de_gi - de_if;

          tmp += tmom_gi*tmom_if*bw/(-de_gi + omega_x[x][y]\
                                     - (fwhm_in/2)*I);

          tmp*=(exp(-(powerl(ediff_x,2))/grms_in)/gvar_in*de_x) * (exp(-(powerl(ediff_y,2))/grms_tr)/gvar_tr*de_y);

        }
        tmp = fabsc(tmp);
        tmp *= tmp;
        rixsmap[x][y] += creal(tmp);
      }

      if ((ksum%pstep) == 0) {
        printf( "    progress: %.2f%%\r", (((float)ksum/(float)progress)*100));
        fflush(stdout);
      }
    }
  }

  printf( "      progress: 100%%\n");
  fflush(stdout);
  printf( "    .. done.\n");

  printf( "  - normalizing the calculated RIXS map ..\n");
  for (x=0; x<max_x; x++) {
    for (y=0; y<max_y; y++) {
      if (rixsmap[x][y] > rmax) {
        rmax = rixsmap[x][y];
      }
    }
  }

  s->sfac = rmax;

  printf( "  - writing RIXS map to file: %s ..",dat_fpstr);
  for (x=0; x<max_x; x++) {
    for (y=0; y<max_y; y++) {
      rixsmap[x][y] = rixsmap[x][y]/rmax;
      fprintf(fp,"%le %le %le\n", (omega_x[x][y])*AUTOEV, omega_y[x][y]*AUTOEV, rixsmap[x][y]);
      fflush(fp);
    }
    fprintf(fp,"\n");
    fflush(fp);
  }
  printf( " done.\n");

  line = NULL;
  len = 0;

  /* construct a gnuplot script from the energy ranges */
  while((read = getline(&line, &len, fp_plot_in)) != -1){
    fprintf(fp_plot_out, "%s",line);
  }

  fprintf(fp_plot_out, "set output \"./%s.png\"\n",md -> inp_fn);
  fprintf(fp_plot_out, "set title \"%s\" font \"Helvetica,40\" offset 0,1\n",md -> inp_fn);
  fprintf(fp_plot_out, "splot [%le:%le][%le:%le] \"./%s.dat\" u (($1+xshift)/1):2:($3*sc) with pm3d title \"\"",omega_x[0][0]*AUTOEV,omega_x[max_x-1][0]*AUTOEV,omega_y[0][0]*AUTOEV,omega_y[0][max_y-1]*AUTOEV,md -> inp_fn);

  if (fclose(fp_plot_in) != 0) {
    fprintf(stderr, "calc_spec.c, function calc_spec: unable to close some of the files:\n%s\n", "../src/plot_template");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  if (fclose(fp_plot_out) != 0) {
    fprintf(stderr, "calc_spec.c, function calc_spec: unable to close some of the files:\n%s\n", plot_fpstr);
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  if (fclose(fp) != 0) {
    fprintf(stderr, "calc_spec.c, function calc_spec: unable to close some of the files:\n%s\n", dat_fpstr);
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  free(line);

  for (j=0; j<max_x; j++) {
    free(omega_x[j]);
    free(omega_y[j]);
    free(rixsmap[j]);
  }
  free(omega_x);
  free(omega_y);

}

void
write_log (spec_info s
           ) {

  int j,k,l; /* looping variables */
  metadata md = s -> md;

  char * log_fpstr = concs(3,md->outpath,md->inp_fn,log_sfx);
  printf( "  - writing parameter log to file: %s ..",log_fpstr);

  fflush(stdout);

  FILE * fp;

  /* open the placeholder file */
  if((fp=fopen(log_fpstr, "w"))==NULL) {
    printf("Cannot open file %s.\n",log_fpstr);
  }
  char date_string[100];
  time_t now = time(NULL);
  struct tm *t = localtime(&now);
  strftime(date_string, sizeof(date_string)-1, "%d/%m/%Y %H:%M", t);

  /* number of screened states in each respective category */
  int n_gs, n_sgs;
  int n_is, n_sis;
  int n_fs, n_sfs;

  int tot_st, tot_t; /* total number of screened transitions and total total
                        (screened + non screened) for a given energy range */

  int gs_idx,fs_idx;

  int flag; /* switch between gathering intensity values and printing them  */

  float t_idx;

  double pt; /* final state printin threshold value */
  double tot_sint, tot_tint; /* total screened and unscreened intensities for
                                all states in a range*/

  double tmp_sint, tmp_tint;
  double sint_tot, tint_tot;

  double e_gs,e_is,e_fs;
  double bw;
  double tmom_jk,tmom_kl;

  double tmp_tmax1,tmp_tmax2;
  double max_ts3; /* maximum intensity for transitions from state s of type 3 */

  n_gs = n_sgs = 0;
  n_is = n_sis = 0;
  n_fs = n_sfs = 0;

  pt = 0.05;

  fprintf(fp, "parameter log for calc_spec calculation on the input file \"%s\", \
date %s.\n\n",md->inpath,date_string);

  fprintf(fp, "note: for all transitions from a certain state in a certain\
energy range, the \nstrongest transitions is marked \"<-max<f>([state index])\" in the log.\n\n");



  fprintf(fp,    "\n====================== START ======================\n" );
  fprintf(fp,      "================ general parameters ===============\n\n");

  fprintf(fp, "screening threshold, %% of maximum screening parameter (if in\
 ground state energy\nrange, maximum boltzmann weight, else, maximum\
transition moment for a given\nmomentum type) values for energy ranges used \
in the calculation: \n\n" );
  for (j=1,t_idx=1; j<=md->state_er[0]; j+=2, t_idx +=0.5) {
    fprintf(fp, "[%f,%f] = %.2f%%\n", md->state_er[j], md->state_er[j+1],md->state_t[(int)t_idx]*100);
  }

  fprintf(fp, "\nbroadening values (eV, fwhm) for energy ranges (eV) used in the calculation: \n" );
  fprintf(fp, "[%f,%f] = %f , Gaussian\n", md->state_er[3], md->state_er[4],md->fwhm[0]);
  fprintf(fp, "[%f,%f] = %f, Gaussian\n", md->state_er[5], md->state_er[6],md->fwhm[1]);

  fprintf(fp, "\ntemperature = %.2f C, %.2f K \n", (float)TEXP, (float)TEXP-273.15);

  fprintf(fp, "\nresolution (eV) in incident and energy transfer direction: %le %le \n",md->res[0],md->res[1]);

  fprintf(fp, "\nscaling factor for rixs map normalization: %le \n",s->sfac);

  fprintf(fp,     "\n================ general parameters ===============\n");
  fprintf(fp,    "======================= END =======================\n\n\n" );

  /* print data for the ground state screening */
  fprintf(fp,    "\n====================== START ======================\n" );
  fprintf(fp,"=== states in group 1, range [%d,%d] ===\n\n",(int)md->state_er[1],(int)md->state_er[2]);

  gs_idx = 0;
  while((gs_idx < s->n_trans) && (gs_idx >= 0)){

    /* screen ground state */
    e_gs = s->trs[2][gs_idx];
    bw = get_rbdist(s->e0,e_gs);

    if (inrange((e_gs-s->e0)*AUTOEV,md->state_er[1],md->state_er[2])) {
      n_gs++;
      if (bw >= md->state_t[1]) {
        n_sgs++;
        k = gs_idx;

        fprintf(fp,"I(1) = %d, E = %f, BW = %f\n\n",(int)s->trs[0][gs_idx], (e_gs-s->e0)*AUTOEV, bw);
      }
    } else {
      break;
    }
    gs_idx = get_inext(s,(int)s->trs[0][gs_idx]);
  }
  fprintf(fp,"\n%d of %d , %d%% of the (1) states were screened out.)\n",\
          n_gs-n_sgs, n_gs,100-(int)((n_sgs/n_gs)*100));

  fprintf(fp,  "\n=== states in group 1, range [%d,%d] ===\n",(int)md->state_er[1],(int)md->state_er[2]);
    fprintf(fp,    "======================= END =======================\n\n\n" );

  /* print data for the intermediate state transitions */
  fprintf(fp,    "====================== START ======================\n" );
  fprintf(fp,"=== states in group 2, range [%d,%d] ===\n\n",(int)md->state_er[3],(int)md->state_er[4]);

  tint_tot = sint_tot = gs_idx = 0;
  tot_st = tot_t = 0;
  tot_sint = tot_tint = 0;

  while((gs_idx < s->n_trans) && (gs_idx >= 0)){

    /* screen ground state */
    e_gs = s->trs[2][gs_idx];
    bw = get_rbdist(s->e0,e_gs);
    flag = 0;
    n_is = n_sis = 0;
    if (inrange((e_gs-s->e0)*AUTOEV,md->state_er[1],md->state_er[2]) &&\
        (bw >= md->state_t[1])) {

      k = gs_idx;

      /* loop over all intermediate transitions from this state */
      tmp_sint = tmp_tint = sint_tot = tint_tot = 0;
      flag = 0;
      n_sis = n_is = 0;
      fprintf(fp,"I(1) = %d, E = %f, BW = %f to.. \n",gs_idx, (e_gs-s->e0)*AUTOEV, bw);
      while((int)s->trs[0][k] == (int)s->trs[0][gs_idx]){

        e_is = s->trs[3][k];
        if (inrange((e_is-s->e0)*AUTOEV,md->state_er[3],md->state_er[4])) {

          tmom_jk = s->trs[4][k];
          /* screen intermediate transition */
          if (s->trs[5][k] == 1) {
            tmp_tmax1 = s -> tmax_d;
          } else {
            tmp_tmax1 = s->tmax_q;
          }
          if ((tmom_jk/tmp_tmax1) > md->state_t[2]){

            if (flag == 0) {
              n_sis++;
              tot_st++;
              tmp_sint += tmom_jk;
              sint_tot += tmom_jk;
              tot_sint += tmom_jk;
            } else {

              fprintf(fp, "  I(2) = %d, E(2) = %f, type %d, <f> = %.10e ,\
<f>/scr = %le%%, <f>/tot = %le%% \n",(int)s->trs[1][k],\
(float)((e_is-s->e0)*AUTOEV),(int)s->trs[5][k], tmom_jk, \
                      (tmom_jk/tmp_sint)*100, (tmom_jk/tmp_tint)*100);
            }
          }
          if (flag == 0) {
            n_is++;
            tot_t++;
            tmp_tint += tmom_jk;
            tint_tot += tmom_jk;
            tot_tint += tmom_jk;
          }
        }
        k++;
        if (((int)s->trs[0][k] != (int)s->trs[0][gs_idx]) && (flag == 0)) {
          k = gs_idx;
          flag = 1;
        } else if (((int)s->trs[0][k] != (int)s->trs[0][gs_idx]) && (flag == 1)) {
          fprintf(fp,"\n  summary, state %d:\n    after screeing out all states %.3f%% below the maximum intensity,\n    - %d of %d (%.2f%%) of the (3) states,\n    - %.2f%% of the total intensity for transitions from this state,\n    was screened out.\n",(int)s->trs[0][gs_idx], (float)(md->state_t[2]*100),n_is-n_sis, n_is,100-(((float)n_sis/(float)n_is)*100), 100-(sint_tot/tint_tot)*100);
          flag = 0;
          break;
        }
      }
    } else {
      break;
    }
    gs_idx = get_inext(s,(int)s->trs[0][gs_idx]);
  }
  fprintf(fp,"\n\nsummary, energy range [%d,%d]:\n  after screening states in this range,\n  - %d out of %d states,\n  - %.3f%% of the total intensity for all transitions from group 1 ([%d,%d]),\n  were screened out.\n", (int)md->state_er[3],(int)md->state_er[4], tot_t - tot_st, tot_t, 100-(tot_sint/tot_tint)*100, (int)md->state_er[1],(int)md->state_er[2]);
  fprintf(fp,"\n=== states in group 2, range [%d,%d] ===\n",(int)md->state_er[3],(int)md->state_er[4]);
  fprintf(fp,    "======================= END =======================\n\n\n" );


  gs_idx = 0;

  /* screen ground state */
  e_gs = s->trs[2][gs_idx];
  bw   = get_rbdist(s->e0,e_gs);
  flag = 0;
  n_is = n_sis = 0;

  k = gs_idx;

  sint_tot = tint_tot = 0;
  tot_sint = tot_tint = 0;
  tot_st   = tot_t = 0;

  /* print data for the final state transitions */
  fprintf(fp,    "====================== START ======================\n" );
  fprintf(fp,"=== states in group 3, range [%d,%d] ===\n\n",(int)md->state_er[5],(int)md->state_er[6]);

  fprintf(fp,"\nprinting states %.2f%% as intense as the transition maximum for a given initial state.\n\n",pt*100);

  while((int)s->trs[0][k] == (int)s->trs[0][gs_idx]){

    e_is = s->trs[3][k];
    if (inrange((e_is-s->e0)*AUTOEV,md->state_er[3],md->state_er[4])) {
      if (flag == 0) {
        n_is++;
      }

      tmom_jk     = s->trs[4][k];
      /* screen intermediate transition */
      if (s->trs[5][k] == 1) {
        tmp_tmax1 = s -> tmax_d;
      } else {
        tmp_tmax1 = s->tmax_q;
      }
      if ((tmom_jk/tmp_tmax1) > md->state_t[2]){
        fprintf(fp,"\nI(2) = %d, E = %f to.. \n\n",(int)s->trs[1][k], (e_is-s->e0)*AUTOEV);
        if ((fs_idx = get_i(s,s->trs[1][k])) != -1) {

          l = fs_idx;
          tmp_sint                       = tmp_tint = 0;
          n_sfs                          = n_fs = 0;
          flag                           = 0;
          max_ts3                        = -1;
          while((int)s->trs[0][l] == (int)s->trs[0][fs_idx]){
            e_fs = s->trs[3][l];
            if (inrange((e_fs-s->e0)*AUTOEV,md->state_er[5],md->state_er[6])) {

              tmom_kl = s->trs[4][l];
              /* screen final state transition */
              if (s->trs[5][l] == 1) {
                tmp_tmax2 = s -> tmax_d;
              } else {
                tmp_tmax2 = s->tmax_q;
              }
              if ((tmom_kl/tmp_tmax2) > md->state_t[3]){
                if (flag   == 0) {
                  n_sfs++;
                  tot_st++;
                  tmp_sint += tmom_kl;
                  sint_tot += tmom_kl;
                  tot_sint += tmom_kl;

                  if (max_ts3 < tmom_kl) {
                    max_ts3 = tmom_kl;
                  }
                } else {
                  if (tmom_kl           >= (max_ts3 * pt)) {
                    fprintf(fp, "  I(3)  = %d, E(3) = %f, type %d, <f> = %le ,<f>/scr = %.2e%%, <f>/tot = %.2e%% ",(int)s->trs[1][l], (e_fs-s->e0)*AUTOEV,(int)s->trs[5][l], tmom_kl, ((tmom_kl/tmp_sint)*100), (tmom_kl/tmp_tint)*100);
                    if(tmom_kl == max_ts3){
                      fprintf(fp,"<- max <f>(%d)",(int)s->trs[0][l]);

                    }
                    fprintf(fp,"\n");
                  }
                }
              }
              if (flag   == 0) {
                tot_t++;
                n_fs++;
                tmp_tint += tmom_kl;
                tint_tot += tmom_kl;
                tot_tint += tmom_kl;
              }
            }
            l++;
            if (((int)s->trs[0][l] != (int)s->trs[0][fs_idx]) && (flag == 0)) {
              l = fs_idx;
              flag = 1;

            } else if (((int)s->trs[0][l] != (int)s->trs[0][fs_idx]) && (flag == 1)) {
              fprintf(fp,"\n  summary, state %d:\n    after screeing out all states %.3f%% below the maximum intensity,\n    - %d of %d (%.2f%%) of the (3) states,\n    - %.2f%% of the total intensity for transitions from this state,\n    was screened out.\n",(int)s->trs[0][fs_idx], md->state_t[2]*100,n_fs-n_sfs, n_fs,100-(((float)n_sfs/(float)n_fs)*100), 100-(sint_tot/tint_tot)*100);
              flag = 0;
            }
          }
        }
      }
    }
    k++;
  }

  fprintf(fp,"\n\nsummary, energy range [%d,%d]:\n  after screening states in this range,\n  - %d out of %d states,\n  - %.3f%% of the total intensity for all transitions from group 2 ([%d,%d]),\n  were screened out.\n", (int)md->state_er[5],(int)md->state_er[6], tot_t - tot_st, tot_t, 100-(tot_sint/tot_tint)*100, (int)md->state_er[3],(int)md->state_er[4]);
  printf( " done.\n" );
  fprintf(fp,"\n=== states in group 3, range [%d,%d] ===\n",(int)md->state_er[5],(int)md->state_er[6]);
  fprintf(fp,    "======================= END =======================\n\n\n" );

  if (fclose(fp) != 0) {
    fprintf(stderr, "calc_spec.c, function write_log: unable to close file:\n%s\n", log_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  free(log_fpstr);
}
