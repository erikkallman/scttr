#ifndef DYN_ARRAY_H
#define DYN_ARRAY_H
#define DEFAULT_CAP 10
#define DEFAULT_INC 10

struct da_s;
typedef struct da_s * da;

struct da_s{
  int n_el; /* number of elements */
  int cap; /* data capacity in number of elements */
  int inc; /* increment upon extension */
  int * a;
};

int
da_extend (da src);

int
da_shrink (da src);

int
da_del_us (da ar,
           int el);

int
da_set (da ar,
        int el,
        int val);

int
da_append (da ar,
           int val);

int
da_get (da a,
        int el);

int
da_getlast (da ar);

da
da_init (int cap,
         int inc);

#endif /* DYN_ARRAY_H */
