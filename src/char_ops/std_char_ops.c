#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "limits.h"
#include "std_num_ops.h"
#include "std_char_ops.h"

/* check if both x and y are dashes */
#define isddash(x,y) ((((x) == '-') && ((y) == '-')) ? 1 : 0)

int
satopow(char * s,
        int len
        ){

  int p=-1;
  int j = len;
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
  printf( "pval = %le\n", pval);

  if (s[0] == '-') {
    nsign = -1;
    j = 1;
  } else {
    nsign = 1;
    j = 0;
  }

  /* locate the decimal point, if there is one */
  for (; isdigit(s[j]); j++){
    printf( "char = %d\n", (s[j] - '0'));

    /* v += (double)((s[j] - '0') * p++); */
    v = v*10 + (s[j] - '0');
  }

  printf( "v.1 = %le\n", v);
  printf( "==pval pre = %le\n", pval);
  if (s[j++] == '.') {
    /* account for all powers below the decimal */
    while(isdigit(s[j]) || (s[j] != '\0')){
      printf( "char = %d\n", (s[j] - '0'));
      v = v*10 + (s[j] - '0');
      printf( "val = %le\n", v);
      pval = pval-1;
      j++;
    }
  }

  printf( "==pval post1 = %le\n", pval);

  if (pval<0) {
    pval = 1/powerl(10,-pval);
  } else {
    pval = powerl(10,pval);
  }

  printf( "==pval post2 = %le\n", pval);
  printf( "v.2 = %le\n", v);
  printf( "v.3 = %le\n", ((double)nsign)*v*pval);
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
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
