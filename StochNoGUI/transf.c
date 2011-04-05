#include <stdio.h>

int main(int argc, char *argv[])
{
  while(!feof(stdin)) {
    int val, res;
    double p;

    res = fscanf(stdin, "%d %lf\n", &val, &p);
    if (res == 2) {
      printf("%d %16.14f\n", val / 100, p);
    }
  }

  return 0;
}
