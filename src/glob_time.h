#ifndef GLOB_TIME_H
#define GLOB_TIME_H
/**
   * @file glob_time.h
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains global variables storing timing values
   (in seconds) obtained using the omp_get_wtime function.
   */

extern double para_t; /* <** Parallel walltime */
extern double serial_t; /* <** Serial walltime */
extern double total_t; /* <** Total walltime */

#endif /* GLOB_TIME_H */
