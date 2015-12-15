#include <stdio.h>
#include <time.h>
#include "timing.h"

char *
get_loctime (char *time_s)
{
  int j;
  time_t rawtime;

  time ( &rawtime );
  char *s = asctime(localtime(&rawtime));

  for (j = 0; j < 19; j++) {
    time_s[j] = s[j];
  }
  time_s[19] = '\0';

  return time_s;
}
