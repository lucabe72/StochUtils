#ifndef __DISTRIB__
#define __DISTRIB__

#define PMF_SIZE 2000

struct distribution {
  struct distribution *next;
  struct sys_distr *group;
  void *window;
  double *values;
  char *name;
};

struct sys_distr {
  int max_distr_dim;
  int n_distr;
  struct distribution *first;
};

void init_distribs(struct sys_distr *sd);
#endif
