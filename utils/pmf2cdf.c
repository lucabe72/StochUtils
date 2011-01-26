#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define MAX 5000

int main(int argc, char *argv[])
{
    double sum = 0;

    while(!feof(stdin)) {
        int res, val;
        double p;

        res = scanf("%d %lf\n", &val, &p);
        if (res == 2) {
            sum += p;
            printf("%d %16.14f\n", val, sum);
        }
    }

    return 0;
}
