#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "limits.h"
#include "std_num_ops.h"
#include "std_char_ops.h"

/* check if both x and y are dashes */
#define isddash(x,y) ((((x) == '-') && ((y) == '-')) ? 1 : 0)

size_t
fpread(void *buffer,
       size_t size,
       size_t mitems,
       size_t offset,
       FILE *fp){

     if (fseek(fp, offset, SEEK_SET) != 0)
         return 0;
     return fread(buffer, size, mitems, fp);
}

/* function strscmp

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
strscmp (const char * str1,
         const char ** str2,
         int n_str) {
  int j;

  for (j=0; j<n_str; j++) {
    if (strstr(str1, str2[j])){
      /* printf( "found %s = %s\n", str1, str2[j]); */
      return j;
    }
  }

  return -1;
}

int
isanyalpha (char * s,
           int len) {
  int j;

  char num_key[6] = {'-','.','E','+','\0'};

  for (j=0; j<len; j++) {
    if ((isalpha(s[j]) &&
        (charinstr(num_key,s[j]) == 0)) ||
        (s[j] == '#') ||
        (s[j] == '*')
        ){
      return 1;
    }
  }

  return 0;
}

int
send_range_qmsg (double * state_er,
                 double eval
                ){
  int j; /* counter for number of values printed */
  int k; /* counter for number of lines of values printed */
  int r; /* response from the query */


  printf("\n\nWhich of the following energy ranges do you want to extend with %le\
?\n",  eval);

  for (j=2,k=1,r=1; r<state_er[0]+1; r++,j++) {
    if (j == 2) {
      if (r > 1) {
        printf( "\n" );
      }
      printf( "%d: ",k);
      k++;
      j = 0;
    }
    printf( "%le ", state_er[r]);
  }

  printf( "\nPick one of the above indices : ");
  scanf("%d",&r);
  return r;
}

int
charinstr (char * str,
           char c){
  int j;
  for (j=0; str[j] != '\0' ; j++) {
    if (str[j] == c) {
      return 1;
    }
  }
  return 0;
}

int
satopow(char * s,
        int len
        ){

  int p=-1;
  int j = len-1;
  int psign;
  double pval = 0;

  /* extract the power of the number from the back of the string*/
  while (j>0) {
    if (s[j] == '+') {
      psign = 1;
      break;
    }
    else if(s[j] == '-'){
      psign = -1;
      break;
    }
    else if(s[j] == '.'){
      psign = 1;
      p = 0;
      break;
    }
    p++;
    j--;
  }
  if (j == 0) {
    return 0;
  }
  /* check so that we're not trying to store a number that is larger than
     what can be contained in a double */
  if (p > MAX_POWERL) {
    fprintf(stderr, "std_char_ops.c, function satof: attempted to store a number of power greater than what can be held in a double (MAX_POWERL = %d < power = %d).\n",MAX_POWERL, p);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  while((p--)>0)
    pval = pval*10 + (s[++j] - '0');
  /* if (psign<0) { */
  /*   pval = 1/pval; */
  /* } */

  return pval*psign;
}

double
satof(char * s,
      int len){

  int j;
  int p = 1; /* every number is atleast to the 0th power of ten */

  int psign; /* sign of the power */
  int nsign; /* sign of the number */

  double v = 0;
  double pval = satopow(s, len);
  /* printf( "==pval origin = %le, len %d\n", pval, len); */
  if (s[0] == '-') {
    nsign = -1;
    j = 1;
  } else {
    nsign = 1;
    j = 0;
  }

  /* locate the decimal point, if there is one */
  for (; isdigit(s[j]) && (j < len); j++){
    v = v*10 + (double)(s[j] - '0');
  }

  /* printf( "v.1 = %le\n", v); */
  /* printf( "==pval pre = %le\n", pval); */
  if (s[j++] == '.') {
    /* account for all powers below the decimal */
    while(isdigit(s[j]) && (s[j] != '\0')){
      /* printf( "char = %d\n", (s[j] - '0')); */
      v = v*10 + (double)(s[j] - '0');
      /* printf( "val = %le\n", v); */
      pval = pval-1;
      j++;
    }
  }

  /* printf( "==pval post1 = %le\n", pval); */

  if (pval<0) {
    pval = 1/powerl(10,-pval);
  } else {
    pval = powerl(10,pval);
  }

  /* printf( "==pval post2 = %le\n", pval); */
  /* printf( "v.2 = %le\n", v); */
  /* printf( "v.3 = %le\n", ((double)nsign)*v*pval); */

  return ((double)nsign)*v*pval;
}

double
sci_atof(char * s){
  double val,pow;
  int sign,i,esign,exp;
  int power(int base,int exp);

  for(i=0;isspace(s[i]);i++)
    ;

  sign=(s[i]=='-')?-1:1;

  if(s[i]=='+' || s[i] == '-')
    i++;

  for(val=0.0;isdigit(s[i]);i++)
    val = (double)10.0 * val + (s[i] - '0');

  if(s[i]=='.')
    i++;

  for(pow=1.0;isdigit(s[i]);i++)
    {
      val = 10.0 * val + (s[i] - '0');
      pow *= 10.0;
    }

  if(s[i]=='e' || s[i] =='E')
    i++;
  if(s[i]=='+' || s[i] =='-')
    {
      esign=s[i];
      i++;
    }

  for(exp=0;isdigit(s[i]);i++)
    exp=10.0 * exp + (s[i] - '0');

  if( esign == '-')
    return (double)sign * (val / pow) / (double)power(10,exp);
  else
    return (double)sign * (val / pow) * (double)power(10,exp);
}

int
isdashes (char * s,
         int len) {
  int j,k;
  char c,last_c;
  for (j=0; j<len; j++) {
    c = s[j];
    if (isddash(c,last_c)) {
      return 1;/* found at least one double dash */
    }
    last_c = c;
  }
  /* we got through the string without finding a double dash" */
  return 0;
}

int
isempty (char * s,
         int len) {
  int j,k;

  static char key[5] = { ' ', '\n', '\t','-', '\0' }; /* key for checking empty strings. */
  for (j=0; j<len; j++) {
    if (strchr(key, s[j]) == NULL) {
      /* as soon as strchr returns something that isnt a match,
         of the key, the string is not empty */
      return 0;/* match found, it isnt empty */
    }
  }
  /* we got through the string without finding a match to key.
     the string is "empty" */
  return 1;
}

int
isdigitin (char * s,
           int len) {
  int j;
  for (j=0; j<len; j++) {
    if (isdigit(s[j])){
      return 1;
    }
  }

  return 0;
}
