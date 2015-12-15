#ifndef CACHE_OPT_H
#define CACHE_OPT_H

extern struct ccfg *cache_cfg;

struct ccfg
{
  int l_sz; /**< Number of bytes per line. */
  int nl; /**< Number of cache lines used by this chunk. */
  int tot_sz; /**< Total size (in bytes) usable by this chunk */

  int n_f; /**< Number single precision floating point words that can be stored on one line of the cache */
  int n_d; /**< Number double precision floating point words that can be stored on one line of the cache */

  float usage; /**< How much (in percent) of the total cache size that is used by this chunk. */
};

struct ccfg *
set_ccfg (float usage);

int
get_row_chunk (int nw, int wsz, struct ccfg *c);

int
cache2file (char *log_fn);

int
cache2str ();


#endif /* CACHE_OPT_H */
