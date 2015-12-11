#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <regex.h>
#include <stdarg.h>
#include <unistd.h>
#include <ctype.h>
#include "std_char_ops.h"
#include "cache_opt.h"
#define MAX_ERROR_MSG 0x1000

const char *base_dir = "/sys/devices/system/cpu/";
struct cache_lvl *root_lvl;
int n_cl;

struct cache
{
  unsigned int lvl; /**< Cache level */
  unsigned int l_sz; /**< Cache line sz (byte) */
  unsigned int tot_sz; /**< Total cache sz (byte)*/
  char type; /**< Cache type. NULL = Unknown type,  D = Data, I = Instruction, U = Unified */
  unsigned int assoc; /**< Cache associativity (number of) */
  unsigned int nl; /**< Number of cache lines */

  unsigned int n_fpl; /**< Number single precision floating point words that can be stored on one line of the cache */
  unsigned int n_dpl; /**< Number double precision floating point words that can be stored on one line of the cache */

  unsigned int n_ccpu; /**< Number of connected CPUs (index) */
  unsigned int *ccpus; /**< This cache is connected to these CPUs */

  char *datdir; /**< The directory from which the data in this truct was read */
  struct cache *next_cache;
  struct cache *last_cache;

};

struct cache_lvl
{
  unsigned int n_caches;
  unsigned int lvl;

  struct cache *root_cache;

  struct cache_lvl *next_lvl;
  struct cache_lvl *last_lvl;

};

FILE *
openfile( const char *dirname, const char *filename, const char *mode )
{
  char pathname[1024];   /* should alwys be big enough */
  FILE *fp;

  sprintf(pathname, "%s/%s", dirname, filename);
  fp = fopen( pathname, mode );

  return fp;
}


/* Compile the regular expression described by "regex_text" into
   "r". */
static int compile_regex (regex_t * r, const char * regex_text)
{
    int status = regcomp (r, regex_text, REG_EXTENDED|REG_NEWLINE);
    if (status != 0) {
        char error_message[MAX_ERROR_MSG];
        regerror (status, r, error_message, MAX_ERROR_MSG);
        printf ("Regex error compiling '%s': %s\n",
                 regex_text, error_message);
        return 1;
    }
    return 0;
}

static int match_regex (regex_t * r, const char * to_match)
{
  int nomatch;
  /* "P" is a pointer into the string which points to the end of the
     previous match. */
  const char * p = to_match;
  /* "N_matches" is the maximum number of matches allowed. */
  const int n_matches = 10;
  /* "M" contains the matches found. */
  regmatch_t m[n_matches];

  while (1) {
    nomatch = regexec (r, p, n_matches, m, 0);
    if (!nomatch) {
      return 1;
    }
    else {
      return 0;
    }

    p += m[0].rm_eo;
  }
  return 0;
}

int
set_ccpus (struct cache *c)
{
  int j;

  long val;
  FILE * p = openfile(c -> datdir, "shared_cpu_list", "r");

  char s1[1024], *s2;
  char int_buf[1024];

  if (p) {
    fscanf(p, "%s", &s1);
    s2 = s1;
    while(*s2) {
      if (isdigit(*s2)) {
        val = strtol(s2, &s2, 10);
        int_buf[c -> n_ccpu++] = val;
      } else {
        s2++;
      }
    }
  } else {
    printf("function set_ccpus: no cachce fae2\n");
  }

  c -> ccpus = malloc(c -> n_ccpu * sizeof(int));

  for (j = 0; j < c -> n_ccpu; j++) {
    c -> ccpus[j] = int_buf[j];
  }

  fclose(p);

  return 0;

}

int
set_size (struct cache *c)
{
  FILE * p = openfile(c->datdir, "size", "r");
  /* printf("%s",c->dat ); */
  unsigned int i = 0;
  if (p) {
    fscanf(p, "%d", &i);
    fclose(p);
  } else {
    printf("function set_size: no cachce fae\n");
  }

  c -> tot_sz = (int)i*100;

  return 0;
}

int
set_lvl (struct cache *c)
{
  FILE * p = openfile(c->datdir, "level", "r");
  /* printf("%s",c->dat ); */
  unsigned int i = 0;
  if (p) {
    fscanf(p, "%d", &i);
    fclose(p);
  } else {
    printf("function set_lvl: no cachce fae\n");
  }

  c -> lvl = (int)i;

  return 0;
}

int
set_lsz (struct cache *c)
{
  FILE * p = openfile(c->datdir, "coherency_line_size", "r");
  /* printf("%s",c->dat ); */
  unsigned int i = 0;
  if (p) {
    fscanf(p, "%d", &i);
    fclose(p);
  } else {
    printf("function set_lsz: no cachce fae\n");
  }

  c -> l_sz = (int)i;

  return 0;
}

int
set_assoc (struct cache *c)
{
  FILE * p = openfile(c->datdir, "ways_of_associativity", "r");
  /* printf("%s",c->dat ); */
  unsigned int i = 0;
  if (p) {
    fscanf(p, "%d", &i);
    fclose(p);
  } else {
    printf("function set_assoc: no cachce fae\n");
  }

  c -> assoc = (int)i;

  return 0;
}

int
set_type (struct cache *c)
{
  FILE * p = openfile(c->datdir, "type", "r");

  char type;
  if (p) {
    fscanf(p, "%c", &type);
    fclose(p);
  } else {
    printf("function set_lsz: no cachce fae\n");
  }

  c -> type = (char)type;

  return 0;
}

struct cache_lvl *
init_cache_lvl (int l)
{
  struct cache_lvl *cl = malloc(sizeof(struct cache_lvl));
  struct cache_lvl *last_cl;

  cl -> n_caches = 0;
  cl -> lvl = l;
  cl -> root_cache = NULL;

  if (n_cl == 0) {
    root_lvl = cl;
    root_lvl -> next_lvl = NULL;
    root_lvl -> last_lvl = NULL;

  } else {

    last_cl = root_lvl;
    while(last_cl -> next_lvl != NULL){
      last_cl = last_cl -> next_lvl;
    }

    last_cl -> next_lvl = cl;
    cl -> last_lvl = last_cl;
    cl -> next_lvl = NULL;
  }

  n_cl++;
  return cl;
}

struct cache_lvl *
get_cache_lvl (struct cache *c)
{
  struct cache_lvl *last_cl = NULL;

  if (n_cl != 0) {

    last_cl = root_lvl;
    while((last_cl != NULL)
          && (last_cl -> lvl != c -> lvl)) {
      last_cl = last_cl -> next_lvl;
    }
  }

  return last_cl;
}

int
add_cache (struct cache *c)
{
  struct cache_lvl *cl = get_cache_lvl(c);
  struct cache * next_c;
  struct cache * last_c;

  if (cl == NULL) { /* lvl not found */
    return 1;
  } else {
    if (cl -> n_caches == 0) {
      cl -> root_cache = c;
      last_c = NULL;
    } else {
      /* add the cache */
      next_c = cl -> root_cache;
      last_c = next_c;

      while( next_c != NULL){
        last_c = next_c;
        next_c = next_c -> next_cache;
      }

      last_c -> next_cache = c;
    }

    c -> next_cache = NULL;
    c -> last_cache = last_c;

    cl -> n_caches++;
  }

  return 0;
}

struct cache *
init_cache (char *c_dir)
{
  struct cache *c = malloc(sizeof(struct cache));
  c -> datdir = c_dir;
  c -> n_ccpu = 0;

  set_lvl(c);
  set_lsz(c);
  set_size(c);
  set_type(c);
  set_assoc(c);
  set_ccpus(c);

  if (add_cache(c) == 1) { /* no cache lvl found for the lvl of c */
    init_cache_lvl(c -> lvl);
    add_cache(c);
  }

  c -> nl = c -> tot_sz / c -> l_sz;
  c -> n_fpl = c -> l_sz/sizeof(float);
  c -> n_dpl = c -> l_sz/sizeof(double);

  return c;
}

int
free_cache (struct cache *c)
{

  free(c -> ccpus);
  free(c -> datdir);
  free(c);

  return 0;
}

int
free_all_caches (void)
{

  struct cache *c;
  struct cache *nc;

  struct cache_lvl *cl;
  struct cache_lvl *ncl = root_lvl;

  while (ncl != NULL) {
    cl = ncl;
    nc = cl -> root_cache;
    while (nc != NULL) {
      c = nc;
      nc = c -> next_cache;
      free_cache(c);
    }
    ncl = cl -> next_lvl;
    free(cl);
  }

  return 0;
}

int
cache2str ()
{
  int j;
  struct cache_lvl *last_cl = root_lvl;
  struct cache *c;

  while(last_cl != NULL){

    printf("Cache lvl %d\n", last_cl -> lvl );
    c = last_cl -> root_cache;

    while(c != NULL){
      if (c -> type != 'I' ) {
        printf("\tlvl = %d\n", c -> lvl);
        printf("\tl_sz = %d\n", c -> l_sz);
        printf("\ttot_sz = %d\n", c -> tot_sz);
        printf("\ttype = %c\n", c -> type);
        printf("\tassoc = %d\n", c -> assoc);
        printf("\tnl = %d\n", c -> nl);
        printf("\tn_ccpu = %d\n", c -> n_ccpu);
        printf("\tn_fpl = %d\n", c -> n_fpl);
        printf("\tn_dpl = %d\n", c -> n_dpl);
        for (j = 0; j < c -> n_ccpu; j++) {
          printf("\t  ccpus[%d] = %d\n", j, c -> ccpus[j]);
        }
        printf("\n" );
      }

      c = c -> next_cache;
    }
    printf("\n" );
    printf("\n" );
    last_cl = last_cl -> next_lvl;
  }

  return 0;
}

int
set_cachetopo ()
{
  struct dirent *curr_dir;
  struct dirent *last_dir;
  struct stat st;

  regex_t r;
  const char *regex_text;

  char *base_subdir;
  char *cache_dir;

  regex_text = "[^[:digit:]]+[[:digit:]]";
  compile_regex(&r, regex_text);

  DIR* d = opendir(base_dir);
  DIR* d_sub;

  if(d == NULL)     {

    printf("list_dir : %s : %s \n", base_dir, strerror(errno));

    return 0;
  }

  while((curr_dir = readdir(d))) {

    /* exclude the . and .. dirs */
    if (stat(curr_dir->d_name, &st) != 0) {

      if (match_regex(&r, curr_dir->d_name)) {
        /* a cpu# directory was found */
        last_dir = curr_dir;

        base_subdir = concs(3,base_dir,curr_dir->d_name,"/cache/\0");
        d_sub = opendir(base_subdir);

        /* loop over the cache directories */
        while((curr_dir = readdir(d_sub))) {
          if (match_regex(&r, curr_dir->d_name)) {

            /* initialize a new cache struct with the information
             contained in this subdir */
            cache_dir = concs(3,base_subdir,curr_dir->d_name,"/\0");
            init_cache(cache_dir);
            free(cache_dir);

          }
        }
        free(base_subdir);
        curr_dir = last_dir;
      }
    }
  }

  closedir(d);
  regfree (& r);

  return 0;
}

struct ccfg *
set_ccfg (float usage)
{
  /* Assume usage% of cache is used for the looping data set. */
  /* Loop over the cache, starting from the lowest cache. */

  set_cachetopo();
  struct cache_lvl *last_cl = root_lvl;
  struct cache *c;
  struct ccfg *cfg = malloc(sizeof(struct ccfg));

  cfg -> l_sz = 0;
  cfg -> usage = usage;

  while(last_cl != NULL){

    c = last_cl -> root_cache;
    while(c != NULL){
      if (c -> type != 'I' ) {
        if((c -> l_sz < cfg -> l_sz) || (cfg -> l_sz == 0)){
          cfg -> l_sz = c -> l_sz;
          cfg -> nl = (int)(c -> nl * usage);
          cfg -> tot_sz = cfg -> l_sz * cfg -> nl;
        } else {
          if ((c -> n_ccpu > 1)
              && ((int)((c -> tot_sz) * usage)
                  <= (c -> n_ccpu) * cfg -> tot_sz)) {
            /* if it is, does n_shared*m fit inside of it? */
            /* if not, reduce the size of m such that it will
             by changing the number of lines in m */
            while(--(cfg -> nl) > 0) {
              cfg -> tot_sz = cfg -> l_sz * cfg -> nl;
              if (((int)((c -> tot_sz) * usage)
                  <= (c -> n_ccpu) * cfg -> tot_sz)) {
                /* the number of lines allcfgated to m have been reduced
                 such that n_ccpu numbers of m will fit in c */
                break;
              }
            }
          }
        }
      }
      c = c -> next_cache;
    }
    last_cl = last_cl -> next_lvl;
  }

  return cfg;
}

int
get_row_chunk (int nw, int wsz, struct ccfg *c) {
  return (int)c -> tot_sz / (nw*wsz);
}
