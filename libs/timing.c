#include <stdio.h>
#include <time.h>
#include <string.h>
#include "timing.h"
#include "std_char_ops.h"
#include "std_num_ops.h"

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

/* calculate the shift in date between months or dates. assumes d2 to be
futher, or equal to, in time than d1. */
int
get_datediff(char *d1, char *d2, const char **date, int ndates)
{
  int j;
  int i_d1, i_d2 = -1;
  int dd = 0; /* index difference in dates */

  for (j = 0; j < ndates; j++) {

    /* avoid the strcmp if the d1 has been found */
    if ((i_d1 == -1) && (strcmp(d1, date[j]) == 0)) {
      i_d1 = j;
    }
    if ((i_d2 == -1) && (strcmp(d2, date[j]) == 0)) {
      i_d2 = j;
    }
  }

  dd = i_d2 - i_d1;

  if ((i_d1 == -1) || (i_d2 == -1)) {
    return -1;
  }
  else {
    return dd;
  }
}

/* Extract the time signature from a date timestamp and return the difference (in time) between them in seconds */
int
get_timediff (char *stamp1, char *stamp2)
{
  int j = 0;
  int ndays, nmonths, nmins, nsek, nhours, ntotsec = 0;
  int time1, time2 = 0;

  /* arrays of hours, minutes, seconds*/
  int time1_a[3];
  int time2_a[3];

  char *time1s = malloc(7*sizeof(char));
  char *time2s = malloc(7*sizeof(char));
  time1s[6] = '\0';
  time2s[6] = '\0';

  char *day1 = malloc(4*sizeof(char));
  char *day2 = malloc(4*sizeof(char));

  char *month1 = malloc(4*sizeof(char));
  char *month2 = malloc(4*sizeof(char));

  /* extract the time digits */
  get_nchars(stamp1, time1s, 11, 19, 2);
  get_nchars(stamp2, time2s, 11, 19, 2);

  /* timestamps without comma seperators */
  time1 = atoi(time1s);
  time2 = atoi(time2s);

  for (j = 0; j < 3; j++) {
    time1_a[j] = get_digit(time1, 2, j*2);
    time2_a[j] = get_digit(time2, 2, j*2);
  }

  /* extract the day */
  get_nchars(stamp1, day1, 4, 7, 1);
  get_nchars(stamp2, day2, 4, 7, 1);

  /* extract the month */
  get_nchars(stamp1, month1, 0, 4, 1);
  get_nchars(stamp2, month2, 0, 4, 1);

  nmonths = get_datediff(month1, month2, MONTHS, NMONTHS);
  ndays = get_datediff(day1, day2, DAYS, NDAYS);

  /* add up all contributions to the total seconds */
  nhours = time2_a[0] - time1_a[0];
  nmins = time2_a[1] - time1_a[1];
  nsek = time2_a[0] - time1_a[0];

  ntotsec = nmonths*SPMONTH + ndays*SPDAY + nhours*SPHOUR + nmins*SPMIN
    + nsek;

  return ntotsec;
}
