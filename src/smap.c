#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include "state.h"
#include "dynarray.h"
#include "smap.h"
#include "parse_input.h" /* struct types */
#include "structs.h"
#include "sci_const.h"
#include "std_num_ops.h" /* power */
#include "spec.h"

static double xshift;

extern int nt;
extern double tmax_d, tmax_q, e0;
extern int * idxs_map;
extern double ** parsed_input;

/* const double xshift = -19/AUTOEV; /\* Fe2pCN 1s -> 3d *\/ */

void
calc_smap_m (char * fn_infile,
             char * dat_fpstr,
             char * plot_fpstr,
             double * state_er,
             double * state_t,
             double * res,
             double * fwhm_inp
             ) {

  FILE * fp;
  FILE * fp_plot_in;
  FILE * fp_plot_out;

  size_t len;
  ssize_t read;

  char * line;

  unsigned long progress;

  unsigned long ksum;

  int j,k,l,m,n; /* control loop indices */
  int maxgridj;
  int maxgridk;
  int jgrid,kgrid;

  int gs_idx,is_idx,fs_idx;
  int n_gs,n_is, n_fs;

  double rmax = -0.1;

  double e_gs,e_is,e_fs, de_jk, de_kl;

  double tmp_tmax1,tmp_tmax2;

  double omegain, omegaut, omega_gs;
  double bw;
  double ediffj,ediff_jk,tmom_jk,ediffk,ediff_km,tmom_kl;
  double eminj,emaxj,emink,emaxk,dej,dek,fwhm_in,fwhm_tr;
  double emin_gs,emax_gs,de_gs;
  double tmp_e;
  /* variables used in the Kramers-Heisenberg formula */
  double c1,c2,tmp;

  /* excitation energy broadening parameters */
  double grms_in; /* gaussian RMS value */
  double gvar_in; /* square root of the variance of the gaussian */
  /* transfer energy broadening parameters */
  double grms_tr;
  double gvar_tr;

  /* get a pointer to the 0th column of iis */
  double ** omega_x;
  double ** omega_y;
  double ** rixsmap;
  double ** tmp_evals = malloc(2*sizeof(double*));

  /* open the placeholder file */
  if((fp=fopen(dat_fpstr, "w"))==NULL) {
        fprintf(stderr,"smap.c, function calc_smap_m: unable to open the output file %s.\n",dat_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  /* open the placeholder file */
  if((fp_plot_in=fopen("../src/plot_template", "r"))==NULL) {
    fprintf(stderr,"smap.c, function calc_smap_m: unable to open the output file %s.\n","../src/plot_template");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((fp_plot_out=fopen(plot_fpstr, "w"))==NULL) {
    fprintf(stderr,"smap.c, function calc_smap_m: unable to open the output file %s.\n",plot_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<2; j++) {
    if((tmp_evals[j] = malloc(nt*sizeof(double))) == NULL ){
      fprintf(stderr, "smap.c:function calc_smap_m, malloc: failed \
to allocate memory for \"tmp_evals[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  /* figure out the energy ranges by analyzing the energy ranges in the input */
  n_gs = n_is = n_fs = 0;
  gs_idx = 0;
  while(gs_idx < nt){

    /* screen ground state */
    e_gs = parsed_input[2][gs_idx];
    if (ISINSIDE((e_gs-e0)*AUTOEV,state_er[1],state_er[2])){

      bw = get_rbdist(e0,e_gs);
      if(bw >= state_t[1]) {
        n_gs++;
        k = gs_idx;
        /* loop over all intermediate transitions from this state */
        while((int)parsed_input[0][k] == (int)parsed_input[0][gs_idx]){

          e_is = parsed_input[3][k];
          if (ISINSIDE((e_is-e0)*AUTOEV,state_er[3],state_er[4])) {
            tmom_jk = parsed_input[4][k];

            /* screen intermediate transition */
            if (parsed_input[5][k] == 1) {
              tmp_tmax1 = tmax_d;
            } else {
              tmp_tmax1 = tmax_q;
            }
            if ((tmom_jk/tmp_tmax1) > state_t[2]){

              /* set the energy range value */
              /* if (pyth_distl(e_is-e0*AUTOEV,state_er[4]) <=\ */
              /*      pyth_distl(e_is-e0*AUTOEV,state_er[3])) { */
              /*   if ((e_is-e0)*AUTOEV > emaxj) { */
              /*     emaxj = (e_is-e0)*AUTOEV; */
              /*   } */
              /* } */
              /* else if ((e_is-e0)*AUTOEV < eminj) { */
              /*   eminj = (e_is-e0)*AUTOEV; */
              /* } */
              tmp_evals[0][n_is++] = (e_is-e0)*AUTOEV;
 /*              else{ */
 /*                fprintf(stderr, "smap.c, fuction calc_smap_m found intermediate\ */
 /* state of energy %le that was outside of the range %le to %le\n", (e_is-e0)*AUTOEV, state_er[3], state_er[4]); */
 /*                printf( "program terminating due to the previous error.\n"); */
 /*                printf( "eminj = %le, emaxj = %le\n", eminj,emaxj); */
 /*                printf( "emink = %le, emaxk = %le\n", emink,emaxk); */
 /*                printf( "p1 = %le, p2 = %le\n", pyth_distl(e_fs-e0*AUTOEV,state_er[4]),pyth_distl(e_fs-e0*AUTOEV,state_er[5])); */
 /*                exit(EXIT_FAILURE); */
 /*              } */

              is_idx = k;
              if ((fs_idx = get_i(parsed_input[1][k])) != -1) {
                de_jk = parsed_input[3][k] - parsed_input[2][gs_idx];

                l = fs_idx;

                while((int)parsed_input[0][l] == (int)parsed_input[0][fs_idx]){
                  e_fs = parsed_input[3][l];
                  if (ISINSIDE((e_fs-e0)*AUTOEV,state_er[5],state_er[6])) {
                    tmom_kl = parsed_input[4][l];

                    /* screen final state transition */
                    if (parsed_input[5][l] == 1) {
                      tmp_tmax2 = tmax_d;
                    } else {
                      tmp_tmax2 = tmax_q;
                    }
                    if ((tmom_kl/tmp_tmax2) > state_t[3]){
                      /* set the energy range value */
                      /* if (pyth_distl(e_fs-e0*AUTOEV,state_er[6]) <=\ */
                      /*     pyth_distl(e_fs-e0*AUTOEV,state_er[5])) { */
                      /*   if ((e_fs-e0)*AUTOEV > emaxk) { */
                      /*     emaxk = (e_fs-e0)*AUTOEV; */
                      /*   } */
                      /* } */
                      /* else if ((e_fs-e0)*AUTOEV < emink) { */
                      /*   emink = (e_fs-e0)*AUTOEV; */
                      /* } */
                      tmp_evals[1][n_fs++] = (e_fs-e0)*AUTOEV;

 /*                      else{ */
 /*                        fprintf(stderr, "smap.c, fuction calc_smap_m found final\ */
 /* state of energy %le that was outside of the range %le to %le\n", (e_fs-e0)*AUTOEV, state_er[5], state_er[6]); */
 /*                        printf( "program terminating due to the previous error.\n"); */
 /*                          printf( "eminj = %le, emaxj = %le\n", eminj,emaxj); */
 /*                          printf( "emink = %le, emaxk = %le\n", emink,emaxk); */
 /*                          printf( "p1 = %le, p2 = %le\n", pyth_distl(e_fs-e0*AUTOEV,state_er[5]),pyth_distl(e_fs-e0*AUTOEV,state_er[6])); */
 /*                        exit(EXIT_FAILURE); */
 /*                      } */
                    }
                  }
                  l++;
                }
              }
            }
          }
          k++;
        }
      }
    } else {
      break;
    }

    gs_idx = get_inext((int)parsed_input[0][gs_idx]);

  }

  eminj = floor((get_min(tmp_evals[0],n_is)-2))/AUTOEV;
  emaxj = ceil((get_maxl(tmp_evals[0],n_is) + 2))/AUTOEV;
  emink = floor((get_min(tmp_evals[1],n_fs)-2))/AUTOEV;
  emaxk = ceil((get_maxl(tmp_evals[1],n_fs) + 2))/AUTOEV;

  fwhm_in = fwhm_inp[0]/AUTOEV;
  fwhm_tr = fwhm_inp[1]/AUTOEV;


  dej = res[0]/AUTOEV;
  dek = res[1]/AUTOEV;

  maxgridj = (int)((emaxj-eminj)/dej);
  maxgridk = (int)((emaxk-emink)/dek);

  if((omega_x = malloc(maxgridj*sizeof(double*))) == NULL ){
    fprintf(stderr, "smap.c:function calc_smap, malloc: failed \
to allocate memory for \"omega_x\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((omega_y = malloc(maxgridj*sizeof(double*))) == NULL ){
    fprintf(stderr, "smap.c:function calc_smap, malloc: failed \
to allocate memory for \"omega_y\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((rixsmap = malloc(maxgridj*sizeof(double*))) == NULL ){
    fprintf(stderr, "smap.c:function calc_smap, malloc: failed \
to allocate memory for \"rixsmap\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<maxgridj; j++) {
    if((rixsmap[j] = malloc(maxgridk*sizeof(double))) == NULL ){
      fprintf(stderr, "smap.c:function calc_smap, malloc: failed \
to allocate memory for \"rixsmap[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }

    if((omega_x[j] = malloc(maxgridk*sizeof(double))) == NULL ){
      fprintf(stderr, "smap.c:function calc_smap, malloc: failed \
to allocate memory for \"omega_x[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }

    if((omega_y[j] = malloc(maxgridk*sizeof(double))) == NULL ){
      fprintf(stderr, "smap.c:function calc_smap, malloc: failed \
to allocate memory for \"omega_y[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  progress = maxgridj*maxgridk;

  /* x */
  grms_in = 2.0*powerl((fwhm_in/(2*sqrt(2*log(2)))),2);
  gvar_in = fwhm_in/(2.0*sqrt(2.0*log(2)))*sqrt(2.0*3.1415927);

  /* y */
  grms_tr = 2.0*powerl((fwhm_tr/(2*sqrt(2*log(2)))),2);
  gvar_tr = fwhm_tr/(2.0*sqrt(2.0*log(2)))*sqrt(2.0*3.1415927);

  printf( "  - calculating the RIXS map .. \n");

  ksum = 0;

  for (jgrid=0; jgrid<maxgridj; jgrid++) {
    omegain = eminj+(jgrid*dej);

    for (kgrid=0; kgrid<maxgridk; kgrid++) {
      ksum++;
      rixsmap[jgrid][kgrid] = 0;
      omegaut = emink+(kgrid*dek);
      omega_x[jgrid][kgrid] = omegain;
      omega_y[jgrid][kgrid] = omegaut;

      gs_idx = get_i(1);

      while(gs_idx < nt){
        /* for (j=1; j<nt+1; j++) { /\* loop over all states *\/ */

        /* screen ground state */
        e_gs = parsed_input[2][gs_idx];
        bw = get_rbdist(e0,e_gs);

        if (ISINSIDE((e_gs-e0)*AUTOEV,state_er[1],state_er[2]) &&\
            (bw >= state_t[1])) {
          k = gs_idx;

          /* printf( "gs[%d] = %d\n",k, (int)parsed_input[0][k]); */
          /* loop over all intermediate transitions from this state */
          while((int)parsed_input[0][k] == (int)parsed_input[0][gs_idx]){

            e_is = parsed_input[3][k];
            if (ISINSIDE((e_is-e0)*AUTOEV,state_er[3],state_er[4])) {
              tmom_jk = parsed_input[4][k];

              /* screen intermediate transition */
              if (parsed_input[5][k] == 1) {
                tmp_tmax1 = tmax_d;
              } else {
                tmp_tmax1 = tmax_q;
              }
              if ((tmom_jk/tmp_tmax1) > state_t[2]){

                is_idx = k;
                /* printf( "  is[%d] = %d, %le\n",k, (int)parsed_input[1][is_idx],parsed_input[3][is_idx]); */

                if ((fs_idx = get_i(parsed_input[1][k])) != -1) {
                  de_jk = parsed_input[3][k] - parsed_input[2][gs_idx];

                  l = fs_idx;

                  while((int)parsed_input[0][l] == (int)parsed_input[0][fs_idx]){
                    e_fs = parsed_input[3][l];
                    /* printf( "is_idx = %d TO fs_idx = %d, e_is = %le e_fs = %le %d\n", (int)parsed_input[0][fs_idx], (int)parsed_input[1][l],(parsed_input[2][fs_idx]-e0)*AUTOEV, (e_fs-e0)*AUTOEV, ISINSIDE((e_fs-e0)*AUTOEV,state_er[5],state_er[6])); */
                    /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
                    /* exit(1); */

                    if (ISINSIDE((e_fs-e0)*AUTOEV,state_er[5],state_er[6])) {
                      /* printf( "YEP!\n" ); */
                      /*   printf( "    fs[%d] = %d\n",l, (int)parsed_input[1][l]); */
                      /*   sleep(1); */

                      tmom_kl = parsed_input[4][l];

                      /* screen final state transition */
                      if (parsed_input[5][l] == 1) {
                        tmp_tmax2 = tmax_d;
                      } else {
                        tmp_tmax2 = tmax_q;
                      }
                      if ((tmom_kl/tmp_tmax2) > state_t[3]){


                        de_kl = (parsed_input[3][l] - parsed_input[2][is_idx]);
                        /* .. excitation energy */
                        ediffj = omega_x[jgrid][kgrid] - de_jk;

                        /* .. energy transfer */
                        /* ediffk = omega_y[jgrid][kgrid] - (de_jk + de_kl); */
                        ediffk = omega_y[jgrid][kgrid] - de_kl;
                        tmp = tmom_jk*tmom_kl*bw*exp(-(powerl(ediffj,2))/grms_in)/gvar_in*dej;
                        tmp *= exp(-(powerl(ediffk,2))/grms_tr)/gvar_tr*dek;
                        rixsmap[jgrid][kgrid] += tmp;
                        /* printf( "map[%d][%d] = %le, ej %le, ek %le, mom %le, gj %le, gk %le\n", jgrid, kgrid, rixsmap[jgrid][kgrid], ediffj,ediffk,tmom_jk*tmom_kl*bw,exp(-(powerl(ediffj,2))/grms_in)/gvar_tr*dej,exp(-(powerl(ediffk,2))/grms_tr)/gvar_tr*dek); */
                      }/* else { */
                      /*   printf( " NOPE! %le %le\n", (tmom_kl/tmp_tmax2), tmom_kl); */
                      /* } */

                    }
                    l++;
                  }
                }
              }
            }
            k++;
          }
        } else {
          break;
        }

        gs_idx = get_inext((int)parsed_input[0][gs_idx]);
        /* printf( "%d %d %d\n",gs_idx , nt, j); */
        /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
        /* exit(1); */
      }

      printf( "    progress: %.2f%%\r", (((float)ksum/(float)progress)*100));
      fflush(stdout);
    }

    /* printf( "    progress: %.2f%%\r", (((float)ksum/(float)progress)*100)); */
    /* fflush(stdout); */
  }

  printf( "      progress: 100%%\n");
  fflush(stdout);
  printf( "    .. done.\n");

  printf( "  - normalizing the calculated RIXS map ..");
  /* normalize the map */
  for (jgrid=0; jgrid<maxgridj; jgrid++) {
    for (kgrid=0; kgrid<maxgridk; kgrid++) {
      if (rixsmap[jgrid][kgrid] > rmax) {
        rmax = rixsmap[jgrid][kgrid];
      }
    }
  }
  printf( " done.\n");

  printf( "  - writing RIXS map to file: %s ..",dat_fpstr);
  for (jgrid=0; jgrid<maxgridj; jgrid++) {
    for (kgrid=0; kgrid<maxgridk; kgrid++) {
      rixsmap[jgrid][kgrid] = rixsmap[jgrid][kgrid]/rmax;
      /* write the map */
      fprintf(fp,"%le %le %le\n", (omega_x[jgrid][kgrid])*AUTOEV, omega_y[jgrid][kgrid]*AUTOEV, rixsmap[jgrid][kgrid]);
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

  fprintf(fp_plot_out, "set output \"./%s.png\"\n",fn_infile);
  fprintf(fp_plot_out, "set title \"%s\" font \"Helvetica,40\" offset 0,1\n",fn_infile);
  fprintf(fp_plot_out, "splot [%le:%le][%le:%le] \"./%s.dat\" u (($1+xshift)/1):2:($3*sc) with pm3d title \"\"",omega_x[0][0]*AUTOEV,omega_x[maxgridj-1][0]*AUTOEV,omega_y[0][0]*AUTOEV,omega_y[0][maxgridk-1]*AUTOEV,fn_infile);

  if (fclose(fp_plot_in) != 0) {
    fprintf(stderr, "smap.c, function calc_smap_m: unable to close some of the files:\n%s\n", "../src/plot_template");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  if (fclose(fp_plot_out) != 0) {
    fprintf(stderr, "smap.c, function calc_smap_m: unable to close some of the files:\n%s\n", plot_fpstr);
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  if (fclose(fp) != 0) {
    fprintf(stderr, "smap.c, function calc_smap_m: unable to close some of the files:\n%s\n", dat_fpstr);
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  free(line);

  for (j=0; j<maxgridj; j++) {
    free(omega_x[j]);
    free(omega_y[j]);
  }
  free(omega_x);
  free(omega_y);

}

void
write_log (double * state_er,
           double * state_t,
           double * res,
           double * fwhm_inp,
           char * fn_relpath,
           char * log_fpstr,
           int n_max
           ) {
  printf( "  - writing parameter log to file: %s ..",log_fpstr);

  fflush(stdout);

  int j,k,l; /* looping variables */

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

  int tmp_ln; /* temporary pointer to line number */

  int gs_idx,is_idx,fs_idx;

  int flag; /* switch between gathering intensity values and printing them  */
  int n_maxem = 5;

  double pt = 0.6; /* final state printin threshold value */
  double s_int, t_int; /* total screened and unscreened intensities for a
                          given state*/
  double tot_sint, tot_tint; /* total screened and unscreened intensities for
                                all states in a range*/

  double tmp_sint, tmp_tint;
  double sint_tot, tint_tot;

  double e_gs,e_is,e_fs;
  double bw;
  double tmom_jk,tmom_kl;

  double tmp_tmax1,tmp_tmax2;
  double max_ts3; /* maximum intensity for transitions from state s of type 3 */
  double * max_em;

  int fst[5] = {0}; /* n_max maximum final state transitions */

  n_gs = n_sgs = 0;
  n_is = n_sis = 0;
  n_fs = n_sfs = 0;

  fprintf(fp, "parameter log for smap calculation on the input file \"%s\", \
date %s.\n\n",fn_relpath,date_string);

  fprintf(fp, "note: for all transitions from a certain state in a certain\
energy range, the \nstrongest transitions is marked \"<-max<f>([state index])\" in the log.\n\n");

  float t_idx = 0;

  fprintf(fp,    "\n====================== START ======================\n" );
  fprintf(fp,      "================ general parameters ===============\n\n");

  fprintf(fp, "screening threshold, %% of maximum screening parameter (if in\
 ground state energy\nrange, maximum boltzmann weight, else, maximum\
transition moment for a given\nmomentum type) values for energy ranges used \
in the calculation: \n\n" );
  for (j=1,t_idx=1; j<=state_er[0]; j+=2, t_idx +=0.5) {
    fprintf(fp, "[%f,%f] = %.2f%%\n", state_er[j], state_er[j+1],state_t[(int)t_idx]*100);
  }

  fprintf(fp, "\nbroadening values (eV, fwhm) for energy ranges (eV) used in the calculation: \n" );
  fprintf(fp, "[%f,%f] = %f , Gaussian\n", state_er[3], state_er[4],fwhm_inp[0]);
  fprintf(fp, "[%f,%f] = %f, Gaussian\n", state_er[5], state_er[6],fwhm_inp[1]);

  fprintf(fp, "\ntemperature = %.2f C, %.2f K \n", (float)TEXP, (float)TEXP-273.15);

  fprintf(fp, "\nresolution (eV) in incident and energy transfer direction: %le %le \n",res[0],res[1]);

  fprintf(fp,     "\n================ general parameters ===============\n");
  fprintf(fp,    "======================= END =======================\n\n\n" );

  /* print data for the ground state screening */
  /* printf( "print data for the ground state transitions\n" ); */
  fprintf(fp,    "\n====================== START ======================\n" );
  fprintf(fp,"=== states in group 1, range [%d,%d] ===\n\n",(int)state_er[1],(int)state_er[2]);

  gs_idx = 0;
  while(gs_idx < nt){

    /* screen ground state */
    e_gs = parsed_input[2][gs_idx];
    bw = get_rbdist(e0,e_gs);

    if (ISINSIDE((e_gs-e0)*AUTOEV,state_er[1],state_er[2])) {
      n_gs++;
      if (bw >= state_t[1]) {
        n_sgs++;
        k = gs_idx;

        fprintf(fp,"I(1) = %d, E = %f, BW = %f\n\n",(int)parsed_input[0][gs_idx], (e_gs-e0)*AUTOEV, bw);
      }
    } else {
      break;
    }
    gs_idx = get_inext((int)parsed_input[0][gs_idx]);
  }
  fprintf(fp,"\n%d of %d , %d%% of the (1) states were screened out.)\n",\
          n_gs-n_sgs, n_gs,100-(int)((n_sgs/n_gs)*100));

  fprintf(fp,  "\n=== states in group 1, range [%d,%d] ===\n",(int)state_er[1],(int)state_er[2]);
    fprintf(fp,    "======================= END =======================\n\n\n" );
  /* append, to the beginning, the data on number of states etc. */

  /* print data for the intermediate state transitions */
  /* printf( "print data for the intermediate state transitions\n" ); */
  fprintf(fp,    "====================== START ======================\n" );
  fprintf(fp,"=== states in group 2, range [%d,%d] ===\n\n",(int)state_er[3],(int)state_er[4]);

  tint_tot = sint_tot = gs_idx = 0;
  tot_st = tot_t = 0;
  tot_sint = tot_tint = 0;

  while(gs_idx < nt){
    /* for (j=1; j<nt+1; j++) { /\* loop over all states *\/ */

    /* screen ground state */
    e_gs = parsed_input[2][gs_idx];
    bw = get_rbdist(e0,e_gs);
    flag = 0;
    n_is = n_sis = 0;
    if (ISINSIDE((e_gs-e0)*AUTOEV,state_er[1],state_er[2]) &&\
        (bw >= state_t[1])) {

      k = gs_idx;
      /* printf( "gs[%d] = %d\n",k, (int)parsed_input[0][k]); */
      /* loop over all intermediate transitions from this state */
      tmp_sint = tmp_tint = sint_tot = tint_tot = 0;
      flag = 0;
      n_sis = n_is = 0;
      fprintf(fp,"I(1) = %d, E = %f, BW = %f to.. \n",gs_idx, (e_gs-e0)*AUTOEV, bw);
      while((int)parsed_input[0][k] == (int)parsed_input[0][gs_idx]){

        e_is = parsed_input[3][k];
        if (ISINSIDE((e_is-e0)*AUTOEV,state_er[3],state_er[4])) {

          tmom_jk = parsed_input[4][k];
          /* screen intermediate transition */
          if (parsed_input[5][k] == 1) {
            tmp_tmax1 = tmax_d;
          } else {
            tmp_tmax1 = tmax_q;
          }
          if ((tmom_jk/tmp_tmax1) > state_t[2]){
            is_idx = k;

            if (flag == 0) {
              n_sis++;
              tot_st++;
              tmp_sint += tmom_jk;
              sint_tot += tmom_jk;
              tot_sint += tmom_jk;
            } else {

              fprintf(fp, "  I(2) = %d, E(3) = %f, type %d, <f> = %.10e ,\
<f>/scr = %le%%, <f>/tot = %le%% \n",(int)parsed_input[1][k],\
(float)((e_is-e0)*AUTOEV),(int)parsed_input[5][k], tmom_jk, \
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
        if (((int)parsed_input[0][k] != (int)parsed_input[0][gs_idx]) && (flag == 0)) {
          k = gs_idx;
          flag = 1;
        } else if (((int)parsed_input[0][k] != (int)parsed_input[0][gs_idx]) && (flag == 1)) {
          /* fprintf(fp,"\n%d of %d (%.2f%%) of the (2) states screened out, along with (%.2f%%) of the total intensity for transitions from state %d.\n",n_is-n_sis, n_is,100-(((float)n_sis/(float)n_is)*100), 100-((float)sint_tot/(float)tint_tot)*100,(int)parsed_input[0][gs_idx]); */
          fprintf(fp,"\n  summary, state %d:\n    after screeing out all states %.3f%% below the maximum intensity,\n    - %d of %d (%.2f%%) of the (3) states,\n    - %.2f%% of the total intensity for transitions from this state,\n    was screened out.\n",(int)parsed_input[0][gs_idx], (float)(state_t[2]*100),n_is-n_sis, n_is,100-(((float)n_sis/(float)n_is)*100), 100-(sint_tot/tint_tot)*100);
          flag = 0;
          break;
        }
      }
    } else {
      break;
    }
    gs_idx = get_inext((int)parsed_input[0][gs_idx]);
  }
  fprintf(fp,"\n\nsummary, energy range [%d,%d]:\n  after screening states in this range,\n  - %d out of %d states,\n  - %.3f%% of the total intensity for all transitions from group 1 ([%d,%d]),\n  were screened out.\n", (int)state_er[3],(int)state_er[4], tot_t - tot_st, tot_t, 100-(tot_sint/tot_tint)*100, (int)state_er[1],(int)state_er[2]);
  fprintf(fp,"\n=== states in group 2, range [%d,%d] ===\n",(int)state_er[3],(int)state_er[4]);
  fprintf(fp,    "======================= END =======================\n\n\n" );


  gs_idx = 0;

  /* screen ground state */
  e_gs = parsed_input[2][gs_idx];
  bw   = get_rbdist(e0,e_gs);
  flag = 0;
  n_is = n_sis = 0;

  k = gs_idx;
  /* printf( "gs[%d] = %d\n",k, (int)parsed_input[0][k]); */
  /* loop over all intermediate transitions from this state */

  sint_tot = tint_tot = 0;
  tot_sint = tot_tint = 0;
  tot_st   = tot_t = 0;
  /* print data for the final state transitions */
  /* printf( "print data for the final state transitions\n" ); */
  fprintf(fp,    "====================== START ======================\n" );
  fprintf(fp,"=== states in group 3, range [%d,%d] ===\n\n",(int)state_er[5],(int)state_er[6]);

  fprintf(fp,"\nprinting states %.2f%% as intense as the transition maximum for a given initial state.\n\n",pt*100);

  while((int)parsed_input[0][k] == (int)parsed_input[0][gs_idx]){

    e_is = parsed_input[3][k];
    if (ISINSIDE((e_is-e0)*AUTOEV,state_er[3],state_er[4])) {
      if (flag == 0) {
        n_is++;
      }

      tmom_jk     = parsed_input[4][k];
      /* screen intermediate transition */
      if (parsed_input[5][k] == 1) {
        tmp_tmax1 = tmax_d;
      } else {
        tmp_tmax1 = tmax_q;
      }
      if ((tmom_jk/tmp_tmax1) > state_t[2]){
        is_idx                                    = k;
        fprintf(fp,"\nI(2) = %d, E = %f to.. \n\n",(int)parsed_input[1][k], (e_is-e0)*AUTOEV);
        if ((fs_idx = get_i(parsed_input[1][k])) != -1) {

          l = fs_idx;
          tmp_sint                       = tmp_tint = 0;
          n_sfs                          = n_fs = 0;
          flag                           = 0;
          max_ts3                        = -1;
          while((int)parsed_input[0][l] == (int)parsed_input[0][fs_idx]){
            e_fs = parsed_input[3][l];
            if (ISINSIDE((e_fs-e0)*AUTOEV,state_er[5],state_er[6])) {

              tmom_kl = parsed_input[4][l];
              /* screen final state transition */
              if (parsed_input[5][l] == 1) {
                tmp_tmax2 = tmax_d;
              } else {
                tmp_tmax2 = tmax_q;
              }
              if ((tmom_kl/tmp_tmax2) > state_t[3]){
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
                    fprintf(fp, "  I(3)  = %d, E(3) = %f, type %d, <f> = %le ,<f>/scr = %.2e%%, <f>/tot = %.2e%% ",(int)parsed_input[1][l], (e_fs-e0)*AUTOEV,(int)parsed_input[5][l], tmom_kl, ((tmom_kl/tmp_sint)*100), (tmom_kl/tmp_tint)*100);
                    if(tmom_kl == max_ts3){
                      fprintf(fp,"<- max <f>(%d)",(int)parsed_input[0][l]);

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
            if (((int)parsed_input[0][l] != (int)parsed_input[0][fs_idx]) && (flag == 0)) {
              l = fs_idx;
              flag = 1;

            } else if (((int)parsed_input[0][l] != (int)parsed_input[0][fs_idx]) && (flag == 1)) {
              fprintf(fp,"\n  summary, state %d:\n    after screeing out all states %.3f%% below the maximum intensity,\n    - %d of %d (%.2f%%) of the (3) states,\n    - %.2f%% of the total intensity for transitions from this state,\n    was screened out.\n",(int)parsed_input[0][fs_idx], state_t[2]*100,n_fs-n_sfs, n_fs,100-(((float)n_sfs/(float)n_fs)*100), 100-(sint_tot/tint_tot)*100);
              flag = 0;
            }
          }
        }
      }
    }
    k++;
  }

  /* fprintf(fp,"\n%d out of %d (%.2f%%) of the (2) states kept after screening. (%.2f%%) of the total intensity for transitions from state %d was screened. \n",n_tsfs, n_tfs,((n_sfs/n_fs)*100), (sint_tot/tint_tot)*100,(int)parsed_input[0][fs_idx]); */
  fprintf(fp,"\n\nsummary, energy range [%d,%d]:\n  after screening states in this range,\n  - %d out of %d states,\n  - %.3f%% of the total intensity for all transitions from group 2 ([%d,%d]),\n  were screened out.\n", (int)state_er[5],(int)state_er[6], tot_t - tot_st, tot_t, 100-(tot_sint/tot_tint)*100, (int)state_er[3],(int)state_er[4]);
  printf( " done.\n" );
  fprintf(fp,"\n=== states in group 3, range [%d,%d] ===\n",(int)state_er[5],(int)state_er[6]);
  fprintf(fp,    "======================= END =======================\n\n\n" );

  if (fclose(fp) != 0) {
    fprintf(stderr, "smap.c, function write_log: unable to close file:\n%s\n", log_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
}
