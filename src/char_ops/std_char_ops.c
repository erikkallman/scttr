#include <string.h>
#include "char_ops.h"

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
    val = 10.0 * val + (s[i] - '0');

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
    return sign * (val / pow) / power(10,exp);
  else

    return sign * (val / pow) * power(10,exp);
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
