#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include "std_char_ops.h"
#include "get_numsl.h"

char*
get_numinstr (char * s,
              int idx,
              int str_len
              ){

  int j,k,l;
  j = k = l = 0;

  int n_digits_found = 0;

  static char num_key[] = {'-','.','E',};

  /* mode = 1 corresponds to reading digits, = 0, to reading anything else */
  int mode = 0;
  char c;
  char last_c;
  char * numstr; /* a string containing the extracted number */
  char num_buf[BUF_SIZE*2];

  /* if we're searching this string more in the future, we can store the index
   of the last digit so that we wont have to loop over the entire string again */

  for (; j<str_len; c = s[j]) {

    /* check if we're still reading a number and avoid dashed lines */
    if (((strchr(num_key,c) != NULL) || (isdigit(c) != 0))      \
        && ((isalpha(last_c) == 0 ) || (last_c == 'E'))) {

      /* check if we're reading the right digit, else ignore the result */
      if (k == idx) {
        num_buf[l++] = c;
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
    j++;
  }

  numstr = malloc(j);
  printf( "num_buf=%s,len=%d\n", num_buf,str_len);
  sleep(1);
  for (k=0; k<l; k++) {
    /* store the number and return a pointer to it.
       this gets freed up by the caller. */
    numstr[k] = num_buf[k];
  }
  return numstr;
}

int
get_numsl (char * str,
           int * idxs_out, /* idexes of numbers in string that we want */
           int str_len,
           int n_idxs,
           ...){
  int j,k;
  va_list argv;
  char * numstr;

  double * tmp_num;
  /* store the memory address to the variadic input arguments so that we can
     assign them locally */
  va_start(argv, n_idxs);

  for (j=0; j<n_idxs; j++) { /* loop over indexes */

    tmp_num = va_arg(argv, double*); /* grab the next vararg */

    /* find the correctly indexed number in the string and store it. */
    numstr = get_numinstr(str,idxs_out[j],str_len);
    printf( "string=%s", str);
    *tmp_num = sci_atof(numstr); /* extract the next memory location */

    free(numstr);
  }

  va_end(argv);
  return 0;
}
