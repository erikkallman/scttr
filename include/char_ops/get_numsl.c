#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "std_char_ops.h"
#include "get_numsl.h"

int
get_numinstr (char * s,
              char * buf,
              int idx,
              int str_len
              ){

  int j,k,l;

  char num_key[5] = {'-','.','E','+','\0'};
  /* mode = 1 corresponds to reading digits, = 0, to reading anything else */
  int mode = 0;
  char c = NULL;
  char last_c = NULL;

  /* if we're searching this string more in the future, we can store the index
   of the last digit so that we wont have to loop over the entire string again */
  l = 0;
  k = 0;

  for (j=0; j<=str_len; c = s[j], j++) {
    /* check if we're still reading a number and avoid dashed lines */
    if (((charinstr(num_key,c) != 0) || (isdigit(c) != 0)) \
        && ((isalpha(last_c) == 0 ) || (last_c == 'E'))) {

      /* check if we're reading the right digit, else ignore the result */
      if (k == idx) {
        buf[l++] = c;
      }

      if (mode == 0) {
        /* now proceeding to read a digit */
        mode = bin_flip(mode);
      }
    }
    /* if we read a non-number character and we're in read mode, we have
     read the end of a digit. */
    else if(mode == 1) {
      if (++k > idx) { /* increase the counter for read numbers */
        break;
      }
      mode = bin_flip(mode); /* go back to reading non numbers */
    }

    /* if in reading number mode and found something non number
       the digit we were reading has now ended. increment k to
       note that we have read a number.
    */
    last_c = c;
  }
  buf[l] = '\0';
  return l;
}

int
get_numsl (char * str,
           int * idxs_out, /* idexes of numbers in string that we want */
           int str_len, /* length of the string in str */
           int n_idxs,
           ...){
  int j,k,l;
  va_list argv;

  double * tmp_num;
  char * num_buf;
  /* store the memory address to the variadic input arguments so that we can
     assign them locally */
  va_start(argv, n_idxs);

  /* a number in the string can maximally be as long as the string */
  num_buf = malloc(str_len);

  for (j=0; j<n_idxs; j++) { /* loop over indexes */

    tmp_num = va_arg(argv, double*); /* grab the next vararg */

    /* find the correctly indexed number in the string and store it. */
    l = get_numinstr(str, num_buf, idxs_out[j],str_len);

    for (k=0; k<l; k++) {}
    *tmp_num = satof(num_buf, k+1); /* extract the next memory location */
  }

  free(num_buf);

  va_end(argv);
  return 0;
}
