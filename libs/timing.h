#ifndef TIMING_H
#define TIMING_H

static const int NDAYS = 7;
static const int NMONTHS = 12;
static const int SPDAY = 86400; /* seconds per day */
static const int SPMONTH = 2592000; /* seconds per month */
static const int SPMIN = 60;
static const int SPHOUR = 3600;

static const char *DAYS[] = {"Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"};

static const char *MONTHS[] = {"Jan",  "Feb",  "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

char *
get_loctime (char *time_s);

int
get_datediff(char *d1, char *d2, const char **date, int ndates);

int
get_timediff (char *stamp1, char *stamp2);
#endif /* TIMING_H */
