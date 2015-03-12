#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dynarray.h"
#include "smap.h"
#include "parse_input.h" /* struct types */
#include "e_state_ll.h" /* for operations on e_state structs */
#include "info_ll.h" /* the ll of input data */
#include "rmap_structs.h"
#include "sci_const.h"
#include "signal.h"
#include "std_num_ops.h" /* power */

void
calc_smap (char * method,
           char * inode_id,
           mdda_s * mdda /* screened indices */
           ) {
  printf( "  -calculating RIXS map.\n\n" );
  /* printf( "calc_smap got this method: %s \n", method); */
  FILE * fp;
  /* open the placeholder file */
  if((fp=fopen("/home/kimchi/dev/rmap/output/map.dat", "w"))==NULL) {
    printf("Cannot open file %s.\n","/home/kimchi/dev/rmap/output/map.dat");
  }

  int j,k,l,m; /* control loop indices */
  int maxgridj = 100;
  int maxgridk = 100;
  int jgrid,kgrid;

  int n_gs = mdda_get(mdda, 0, 0);

  int n_is,n_fs;
  int gs_idx,is_idx,fs_idx;

  info_node iroot = get_inode(inode_id);

  double rmax = -0.1;
  double xshift = 7000.0/AUTOEV;
  double omegain, omegaut;
  double ediffj,ediff_jk,tmom_jk,ediffk,ediff_km,tmom_km;
  double eminj,emaxj,emink,emaxk,dej,dek,fwhm;

  /* variables used in the Kramers-Heisenberg formula */
  double c1,c2,tmp;
  mdda_s * root_mdda = (mdda -> root);
  mdda_s * curr_mdda = mdda;
  mdda_s * next_mdda = (mdda -> next);

  /* get a pointer to the 0th column of iis */
  mdda_s * iis = root_mdda -> branch;
  curr_mdda = next_mdda;

  double ** omega_x;
  double ** omega_y;
  double ** rixsmap;
  double ** deltae;
  e_state tmp_gs, tmp_is;

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

  if((deltae = malloc(maxgridj*sizeof(double*))) == NULL ){
    fprintf(stderr, "smap.c:function calc_smap, malloc: failed \
to allocate memory for \"deltae\"\n");
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

    if((deltae[j] = malloc(maxgridk*sizeof(double))) == NULL ){
      fprintf(stderr, "smap.c:function calc_smap, malloc: failed \
to allocate memory for \"deltae[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }
  mdda2s(root_mdda);
  for (j=1; j<n_gs+1; j++) { /* loop over ground states */

    gs_idx= mdda_get(mdda, 0, j);
    tmp_gs = get_state(iroot,gs_idx);
    n_is = mdda_get(mdda, j, 0);
    for (k=1; k<n_is+1; k++) { /* loop over intermediate states */
      is_idx = mdda_get(mdda, j, k);
      deltae[j][k] = get_ediff(iroot, gs_idx, is_idx) + xshift;
      deltae[j][k] = deltae[k][j];
    }
  }

  fwhm = (double)0.9/AUTOEV;

  eminj = (double)6900/AUTOEV;
  emaxj = eminj + (double)200/AUTOEV;
  dej = (emaxj-eminj)/(double)maxgridj;

  emink = -(double)2/AUTOEV;
  emaxk = emink + (double)15/AUTOEV;
  dek = (emaxk-emink)/(double)maxgridk;

  c1 = 2*powerl((fwhm/(2*sqrt(2*log(2)))),2);
  c2 = fwhm/(2*sqrt(2*log(2)))*sqrt(2*3.1415927);

  tmp = 0;
  for (jgrid=0; jgrid<maxgridj; jgrid++) {
    omegain = eminj+(jgrid*dej);

    for (kgrid=0; kgrid<maxgridk; kgrid++) {
      rixsmap[jgrid][kgrid] = 0;
      omegaut = emink+(kgrid*dek);

      omega_x[jgrid][kgrid] = omegain;
      omega_y[jgrid][kgrid] = omegaut;
      for (j=1; j<n_gs+1; j++) { /* loop over ground states */
        /* printf( "ut = %le, in = %le", omegaut, omegain); */
        gs_idx= mdda_get(mdda, 0, j);

        tmp_gs = get_state(iroot,gs_idx);
        n_is = mdda_get(mdda, j, 0);
        for (k=1; k<n_is+1; k++) { /* loop over intermediate states */
          is_idx = mdda_get(mdda, j, k);
          tmom_jk = get_trans(tmp_gs, is_idx);

          /* ediff_jk = get_ediff(iroot, gs_idx, is_idx); */

            /* the FS idices are not necessarily in order, so we have to look up
             the right offset for the right FS index first */
          for (l=1; l<n_is+1; l++) {
            if (mdda_get(iis, 0, l) == is_idx) {
              break;
            }
          }

          tmp_is = get_state(iroot,is_idx);
          n_fs = mdda_get(iis, l, 0);

          for (m=1; m<n_fs+1; m++) {/* loop over final states */
            fs_idx = mdda_get(iis, l, m);

            tmom_km = get_trans(tmp_is, fs_idx);
            /* ediff_km = get_ediff(iroot, is_idx, fs_idx); */

            /* ediffj = omega_x[jgrid][kgrid] - ediff_jk + xshift; */
            /* ediffk = omega_y[jgrid][kgrid] - (ediff_jk - ediff_km + 2*xshift); */
            ediffj = omega_x[jgrid][kgrid] - deltae[j][k];
            ediffk = omega_y[jgrid][kgrid] - (deltae[j][k] - deltae[k][m]);
            /* printf( "ej = %le, ejk = %le , ek = %le, ekm = %le\n", ediffj, \ */
                    /* ediff_jk, ediffk, ediff_km); */
            /* printf( "tmom_jk = %le\n", tmom_jk); */
            /* printf( "tmom_km = %le\n", tmom_km); */

            /* printf( "tmp0=%le\n", tmp); */
            /* printf( "tmoms=%le, exp=%le, dej=%le\n", tmom_jk*tmom_km, exp(-(powerl(ediffj,2))/c1), c2*dej); */
            /* printf( "power=%le, total=%le\n", (powerl(ediffj,2)), -powerl(ediffj,2)/c1); */
            tmp = tmom_jk*tmom_km*exp(-(powerl(ediffj,2))/c1)/c2*dej;
            /* printf( "tmp1=%le\n", tmp); */
            tmp *= exp(-(powerl(ediffk,2))/c1)/c2*dek;
            /* printf( "tmp2=%le\n", tmp); */
            rixsmap[jgrid][kgrid] += tmp;

            /* sleep(1); */
          }
        }
      }
    }
  }
  int jgmx,kgmx;
  /* normalize the map */
  for (jgrid=0; jgrid<maxgridj; jgrid++) {
    for (kgrid=0; kgrid<maxgridk; kgrid++) {
      if (rixsmap[jgrid][kgrid] > rmax) {
        rmax = rixsmap[jgrid][kgrid];
        jgmx = jgrid;
        kgmx = kgrid;
      }
    }
  }
  printf( "  -writing RIXS map to file:\n    %s\n\n","/home/kimchi/dev/rmap/output/map.dat");
  for (jgrid=0; jgrid<maxgridj; jgrid++) {
    for (kgrid=0; kgrid<maxgridk; kgrid++) {
      rixsmap[jgrid][kgrid] = rixsmap[jgrid][kgrid]/rmax;
      /* write the map */
      fprintf(fp,"%le %le %le\n", omega_x[jgrid][kgrid]*AUTOEV, omega_y[jgrid][kgrid]*AUTOEV, rixsmap[jgrid][kgrid]);
    }
    fprintf(fp,"\n");
  }

  for (j=0; j<maxgridj; j++) {
    free(omega_x[j]);
     free(omega_y[j]);
  }
  free(omega_x);
  free(omega_y);
  fclose(fp);
  mdda_free(mdda);

}
