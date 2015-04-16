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

static double xshift;

/* const double xshift = -19/AUTOEV; /\* Fe2pCN 1s -> 3d *\/ */

void
calc_smap_m (char * method,
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

  int gs_idx,is_idx,fs_idx;
  int n_gfs = mdda_get(mdda, 0, 0);

  int n_is,n_fs;

  info_node iroot = get_inode(inode_id);

  double e_gs = iroot -> root_e_state -> e_val;
  double rmax = -0.1;

  double omegain, omegaut;
  double bw;
  double bw_sum = iroot -> bw_sum;
  double ediffj,ediff_jk,tmom_jk,ediffk,ediff_km,tmom_km;
  double eminj,emaxj,emink,emaxk,dej,dek,fwhm,fwhm_trs;
  double tmp_e;
  /* variables used in the Kramers-Heisenberg formula */
  double c1,c2,tmp;
  mdda_s * root_mdda = (mdda -> root);
  mdda_s * curr_mdda = mdda;
  mdda_s * next_mdda = (mdda -> next);

  /* get a pointer to the 0th column of iis */
  mdda_s * igs = root_mdda;
  mdda_s * iis = root_mdda -> branch;
  int n_is_proc = mdda_get(iis,0,0);
  curr_mdda = next_mdda;

  double ** omega_x;
  double ** omega_y;
  double ** rixsmap;
  e_state tmp_gs, tmp_is, tmp_fs;

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
  printf( "  -printing screening parameters:\n");
  printf( "    maximum IS transition intensity = %le\n",  (iroot -> mt_is));
  printf( "    maximum FS transition intensity = %le\n\n",  (iroot -> mt_fs));
  /* mdda_show(root_mdda); */
  /* mdda2s(root_mdda); */
  /* e_statelist2s(iroot,1); */

  /* fwhm = (double)0.9/AUTOEV; */
  fwhm = (double)1.06/AUTOEV;
  /* for Fe2p 2p->3d transitions */
  /* eminj = (double)(716/AUTOEV); */
  /* emaxj = (double)(744/AUTOEV); */
  /* dej = (emaxj-eminj)/(double)maxgridj; */

  /* for Fe2p 1s->3d transitions */
  xshift = -40/AUTOEV;
  eminj = (double)(7148/AUTOEV);
  emaxj = (double)(7160/AUTOEV);
  dej = (emaxj-eminj)/(double)maxgridj;

  /* for Fe3p 2p->3d transitions */
  /* eminj = (double)(726/AUTOEV); */
  /* emaxj = (double)(750/AUTOEV); */
  /* dej = (emaxj-eminj)/(double)maxgridj; */

  /* for Fe3p 1s->p3d transitions */
  /* xshift = -49.5/AUTOEV; */
  /* eminj = (double)(7158/AUTOEV); */
  /* emaxj = (double)(7168/AUTOEV); */
  /* dej = (emaxj-eminj)/(double)maxgridj; */

  /* Fe2pCN 1s -> 3d */
  /* eminj = (double)(7125/AUTOEV); */
  /* emaxj = (double)(7140/AUTOEV); */
  /* dej = (emaxj-eminj)/(double)maxgridj; */

  emink = -(double)(2/AUTOEV);
  emaxk = emink + (double)(14/AUTOEV);
  dek = (emaxk-emink)/(double)maxgridk;

  /* emink = (double)(1390/AUTOEV); */
  /* emaxk = emink + (double)(25/AUTOEV); */
  /* dek = (emaxk-emink)/(double)maxgridk; */

  c1 = 2.0*powerl((fwhm/(2*sqrt(2*log(2)))),2);
  c2 = fwhm/(2.0*sqrt(2.0*log(2)))*sqrt(2.0*3.1415927);

  for (jgrid=0; jgrid<maxgridj; jgrid++) {
    omegain = eminj+(jgrid*dej);

    for (kgrid=0; kgrid<maxgridk; kgrid++) {
      rixsmap[jgrid][kgrid] = 0;
      omegaut = emink+(kgrid*dek);

      omega_x[jgrid][kgrid] = omegain;
      omega_y[jgrid][kgrid] = omegaut;
      for (j=1; j<n_gfs+1; j++) { /* loop over ground states */

        /* printf( "ut = %le, in = %le", omegaut, omegain); */
        gs_idx= mdda_get(igs, 0, j);
        tmp_gs = get_state_si(iroot,gs_idx);
        /* printf( "[%d]/%d, %le\n", gs_idx,n_gfs, tmp_gs -> e_val); */
        /* sleep(1); */
        /* printf( "%le\n", tmp_gs -> e_val); */
        /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
        /* exit(1); */

        /* bw = get_rbdist(e_gs, tmp_gs -> e_val); */
        bw = (tmp_gs -> bw);
        n_is = mdda_get(igs, j, 0);
        for (k=1; k<n_is+1; k++) { /* loop over intermediate states */
          is_idx = mdda_get(igs, j, k);

          tmom_jk = get_trans(tmp_gs, is_idx);

          /* ediff_jk = get_ediff(iroot, gs_idx, is_idx); */

          /* the FS idices are not necessarily in order, so we have to look up
             the right offset for the right FS index first */
          for (l=1; l<n_is_proc+1; l++) {
            if (mdda_get(iis, 0, l) == is_idx) {
              break;
            }
          }
          /* printf( "l=%d/%d", l, n_is); */
          tmp_is = get_state_si(iroot,is_idx);
          /* printf( "  [%d,%d]/%d, %le\n", gs_idx,is_idx,n_is, tmom_jk); */

          n_fs = mdda_get(iis, l, 0);

          for (m=1; m<n_fs+1; m++) {/* loop over final states */

            fs_idx = mdda_get(iis, l, m);

            tmp_fs = get_state_si(iroot,fs_idx);
            tmom_km = get_trans(tmp_is, fs_idx);

            tmp_e = get_ediff(iroot, gs_idx, is_idx);

            /* differences for.. */
            ediffj = omega_x[jgrid][kgrid] - tmp_e; /* .. excitation energy */

            /* .. energy transfer */
            ediffk = omega_y[jgrid][kgrid] - (tmp_e + get_ediff(iroot, is_idx, fs_idx));

            tmp = tmom_jk*tmom_km*bw*exp(-(powerl(ediffj,2))/c1)/c2*dej;

            /* add lorentzian broadening w.r.t the excitation energy and
             gaussian broadening w.r.t the energy transfer */
            tmp *= lorz(ediffj,1.25);
            tmp *= exp(-(powerl(ediffk,2))/c1)/c2*dek;
            rixsmap[jgrid][kgrid] += tmp;
          }
        }
      }
    }
  }

  /* normalize the map */
  for (jgrid=0; jgrid<maxgridj; jgrid++) {
    for (kgrid=0; kgrid<maxgridk; kgrid++) {
      if (rixsmap[jgrid][kgrid] > rmax) {
        rmax = rixsmap[jgrid][kgrid];
      }
    }
  }
  printf( "  -writing RIXS map to file:\n    %s\n\n","/home/kimchi/dev/rmap/output/map.dat");
  for (jgrid=0; jgrid<maxgridj; jgrid++) {
    for (kgrid=0; kgrid<maxgridk; kgrid++) {
      rixsmap[jgrid][kgrid] = rixsmap[jgrid][kgrid]/rmax;
      /* write the map */
      fprintf(fp,"%le %le %le\n", (omega_x[jgrid][kgrid] + xshift)*AUTOEV, omega_y[jgrid][kgrid]*AUTOEV, rixsmap[jgrid][kgrid]);
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

void
calc_smap_dbg (char * method,
               char * inode_id,
               mdda_s * mdda /* screened indices */
               ) {
  printf( "  -calculating RIXS map for debugging purposes.\n\n" );
  /* printf( "calc_smap got this method: %s \n", method); */
  FILE * fp;
  /* open the placeholder file */
  if((fp=fopen("/home/kimchi/dev/rmap/output/map_dbg.dat", "w"))==NULL) {
    printf("Cannot open file %s.\n","/home/kimchi/dev/rmap/output/map_dbg.dat");
  }

  int j,k,l,m; /* control loop indices */
  int maxgridj = 100;
  int maxgridk = 100;
  int jgrid,kgrid;

  int gs_idx,is_idx,fs_idx;
  int n_gfs = mdda_get(mdda, 0, 0);

  int n_is,n_fs;

  info_node iroot = get_inode(inode_id);

  double e_gs = iroot -> root_e_state -> e_val;
  double rmax = -0.1;

  double omegain, omegaut;
  double bw;
  double bw_sum = iroot -> bw_sum;
  double ediffj,ediff_jk,tmom_jk,ediffk,ediff_km,tmom_km;
  double eminj,emaxj,emink,emaxk,dej,dek,fwhm;

  /* variables used in the Kramers-Heisenberg formula */
  double c1,c2,tmp;
  mdda_s * root_mdda = (mdda -> root);
  mdda_s * curr_mdda = mdda;
  mdda_s * next_mdda = (mdda -> next);

  /* get a pointer to the 0th column of iis */
  mdda_s * igs = root_mdda;
  mdda_s * iis = root_mdda -> branch;
  int n_is_proc = mdda_get(iis,0,0);
  curr_mdda = next_mdda;

  double ** omega_x;
  double ** omega_y;
  double ** rixsmap;

  e_state tmp_gs, tmp_is, tmp_fs;

  printf( "  -printing screening parameters:\n");
  printf( "    maximum IS transition intensity = %le\n",  (iroot -> mt_is));
  printf( "    maximum FS transition intensity = %le\n\n",  (iroot -> mt_fs));
  mdda_show(root_mdda);
  /* mdda2s(root_mdda); */
  /* e_statelist2s(iroot,1); */

  for (j=1; j<n_gfs+1; j++) { /* loop over ground states */

    gs_idx= mdda_get(igs, 0, j);
    tmp_gs = get_state_si(iroot,gs_idx);

    n_is = mdda_get(igs, j, 0);
    for (k=1; k<n_is+1; k++) { /* loop over intermediate states */
      is_idx = mdda_get(igs, j, k);

      tmom_jk = get_trans(tmp_gs, is_idx);

      /* the FS idices are not necessarily in order, so we have to look up
         the right offset for the right FS index first */
      for (l=1; l<n_is_proc+1; l++) {
        if (mdda_get(iis, 0, l) == is_idx) {
          break;
        }
      }

      tmp_is = get_state_si(iroot,is_idx);
      n_fs = mdda_get(iis, l, 0);

      for (m=1; m<n_fs+1; m++) {/* loop over final states */

        fs_idx = mdda_get(iis, l, m);
        tmp_fs = get_state_si(iroot,fs_idx);

        tmom_km = get_trans(tmp_is, fs_idx);

        fprintf(fp,"%d %d %d %le %le %le\n", gs_idx, is_idx, fs_idx, (get_ediff(iroot, gs_idx, is_idx)+xshift)*AUTOEV, (get_ediff(iroot, gs_idx, is_idx) + get_ediff(iroot, is_idx, fs_idx))*AUTOEV, tmom_jk*tmom_km);

      }
    }
  }

  fclose(fp);
  mdda_free(mdda);

}
