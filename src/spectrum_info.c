#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sci_const.h"
#include "formats.h"
#include "spectrum_info.h"

spec_info root_sinfo;
int n_sinfo = 0;

/* suffixes for the output and input files */
const char * dat_sfx = ".dat";
const char * plot_sfx = ".gp";
const char * log_sfx = ".txt";
const char * bin_sfx = ".bin";
const char * tmp_sfx = ".tmp";

metadata
init_md () {
  metadata new_md;

  if((new_md = malloc(sizeof(struct metadata_s))) == NULL ){
    fprintf(stderr, "spectrum_info.c:function init_md, malloc: failed \
to allocate memory for \"new_md\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  return new_md;
}

int
free_md (metadata md) {

  free(md -> outpath);
  free(md -> inpath);
  free(md -> inp_fn);
  free(md);
  return EXIT_SUCCESS;
}

spec_info
get_sinfo (char * id
           ){

  spec_info curr_sinfo = root_sinfo;
  spec_info next_sinfo;

    /* locate the spectrum_info node corresponding to the file name input argument */
  while(strstr((curr_sinfo -> md -> inp_fn),id) == NULL) {

    next_sinfo = curr_sinfo -> next;
    curr_sinfo = next_sinfo;

    if (curr_sinfo == NULL) {
      fprintf(stderr, "parse_input.c: get_sinfo, no spec_info node can be found\
 with a str_id == %s\n", id);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  return curr_sinfo;
}


screen
init_screen (spec_info s){

  screen scr = malloc(sizeof(struct screen_s));

  scr -> emin_x = (s->md->state_er[4]-2)/AUTOEV;
  scr -> emax_x = (s->md->state_er[3]+2)/AUTOEV;
  scr -> emin_y = (s->md->state_er[6]-2)/AUTOEV;
  scr -> emax_y = (s->md->state_er[5]+2)/AUTOEV;

  return scr;
}

int
free_screen (screen scr) {

  free(scr -> is2fs);
  free(scr -> is_idxs);
  free(scr -> gs2is);

  return EXIT_SUCCESS;
}

spec_info
init_sinfo (metadata md
            ){

  static spec_info last_sinfo;

  spec_info new_sinfo;

  if((new_sinfo = malloc(sizeof(struct spec_info_s))) == NULL ){
    fprintf(stderr, "spectrum_info.c:function init_sinfo, malloc: failed \
to allocate memory for \"new_sinfo\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  new_sinfo -> md = md;

  if (n_sinfo == 0) { /* there is no root info node defined  */
    n_sinfo = 1;
    new_sinfo -> last = NULL;
    root_sinfo = new_sinfo;
    last_sinfo = root_sinfo;
  }
  else { /* the root info node is already defined */
    n_sinfo++;
    new_sinfo -> last = last_sinfo;
    last_sinfo -> next = new_sinfo;
  }
  new_sinfo -> next = NULL;
  new_sinfo -> idx = n_sinfo-1;

  return new_sinfo;
}

int
free_sinfo (spec_info s) {
  int j;

  for (j=0; j<6; j++) {
    free(s->trs[j]);
  }
  free(s->trs);
  free(s->idx_map);

  free_md(s->md);
  free_screen(s->scr);

  free(s);

  return EXIT_SUCCESS;
}
