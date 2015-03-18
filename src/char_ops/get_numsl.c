#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include "std_char_ops.h"
#include "get_numsl.h"

int
get_numinstr (char * s,
              char * buf,
              int idx,
              int str_len
              ){

  int j,k,l;

  int n_digits_found = 0;
  char num_key[5] = {'-','.','E','+'};

  /* mode = 1 corresponds to reading digits, = 0, to reading anything else */
  int mode = 0;
  char c;
  char last_c;
  char * numstr; /* a string containing the extracted number */

  /* if we're searching this string more in the future, we can store the index
   of the last digit so that we wont have to loop over the entire string again */
  l = 0;
  k = 0;
  for (j=0; j<str_len; c = s[j], j++) {

    /* check if we're still reading a number and avoid dashed lines */
    if (((strchr(num_key,c) != NULL) || (isdigit(c) != 0)) \
        && ((isalpha(last_c) == 0 ) || (last_c == 'E'))) {

      /* check if we're reading the right digit, else ignore the result */
      if (k == idx) {
        /* printf( "reading @k=%d\n", k); */
        /* printf( "string = %s\n",s ); */
        /* printf( "char = %c\n",c ); */
        /* printf( "\n\n\n" ); */
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
        /* printf( "got all chars!\n"); */
        /* sleep(1); */
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

  return l+1;
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
  char * numstr;
  char * num_buf;
  /* store the memory address to the variadic input arguments so that we can
     assign them locally */
  va_start(argv, n_idxs);

  num_buf = malloc(BUF_SIZE*2);

  for (j=0; j<n_idxs; j++) { /* loop over indexes */

    tmp_num = va_arg(argv, double*); /* grab the next vararg */

    /* find the correctly indexed number in the string and store it. */
    l = get_numinstr(str, num_buf, idxs_out[j],str_len);

    numstr = malloc(l);
    /* printf( "\n\nget_numsl got this in return:\n" ); */
    for (k=0; k<l-1; k++) {
      /* store the number and return a pointer to it.
         this gets freed up by the caller. */
      numstr[k] = num_buf[k];
      /* printf( "%c", numstr[k]); */
    }

    /* *tmp_num = sci_atof(numstr); /\* extract the next memory location *\/ */
    *tmp_num = satof(numstr, k); /* extract the next memory location */

    free(numstr);
    /* numstr = NULL; */
  }

  free(num_buf);
  /* num_buf = NULL; */
  va_end(argv);
  return 0;
}
