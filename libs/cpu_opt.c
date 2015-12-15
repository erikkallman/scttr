#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "scttr_cfg.h"
#include "cpu_opt.h"
#include "libcpuid.h"
#include "std_char_ops.h"

int
set_ccnuma_affinity (void)
{
  int j;
  int n_cpus;
  int *cpu_idx;
  /* strings used for creating the environment variable  */
  char *cpu_map;
  char numstr[15];

  if (!cpuid_present()) {
    printf("Sorry, your CPU doesn't support CPUID!\n");
    return -1;
  }

  struct cpu_raw_data_t raw;
  struct cpu_id_t data;

  n_cpus = env2int("OMP_NUM_THREADS");

  if (!n_cpus) {
    if (cpuid_get_raw_data(&raw) < 0) {
      printf("Sorry, cannot get the CPUID raw data.\n");
      printf("Error: %s\n", cpuid_error());
      return -2;
    }
    if (cpu_identify(&raw, &data) < 0) {
      printf("Sorrry, CPU identification failed.\n");
      printf("Error: %s\n", cpuid_error());
      return -3;
    }
    n_cpus = data.num_cores;
  }

  cpu_idx = malloc(n_cpus * sizeof(int));
  for (j = 0; j < n_cpus; j++) {
    cpu_idx[j] = j;
  }

  cpu_map = a2str(cpu_idx, n_cpus,' ');

  if (COMPILER_TYPE == 1){
    setenv ("GOMP_CPU_AFFINITY", cpu_map, 1);
    sprintf(numstr, "%d", n_cpus);
    setenv ("OMP_NUM_THREADS", numstr, 1);
    sprintf(numstr, "%d", 0);
    setenv ("OMP_DYNAMIC", numstr, 1);
  }

  free(cpu_map);
  return 1;
}
