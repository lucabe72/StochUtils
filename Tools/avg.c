#include <stdio.h>
#include <math.h>

static const double t_distribution[30][2] = {
        {6.314, 12.706},
        {2.920, 4.303},
        {2.353, 3.182},
        {2.132, 2.776},
        {2.015, 2.571},
        {1.943, 2.447},
        {1.895, 2.365},
        {1.860, 2.306},
        {1.833, 2.262},
        {1.812, 2.228},
        {1.796, 2.201},
        {1.782, 2.179},
        {1.771, 2.160},
        {1.768, 2.145},
        {1.753, 2.131},
        {1.746, 2.120},
        {1.740, 2.110},
        {1.734, 2.101},
        {1.729, 2.093},
        {1.725, 2.086},
        {1.721, 2.080},
        {1.717, 2.074},
        {1.714, 2.069},
        {1.711, 2.064},
        {1.708, 2.060},
        {1.706, 2.056},
        {1.703, 2.052},
        {1.701, 2.048},
        {1.699, 2.045},
        {1.697, 2.042}
};

int main(int argc, char *argv[])
{
    double val, sum = 0, sum2 = 0;
    int res, n = 0;

    while(!feof(stdin)) {
        res = scanf("%lf\n", &val);
        if (res == 1) {
          n++;
          sum += val;
          sum2 += val * val;
        }
    }

    printf("%f %f\n", sum / n, sqrt(sum2 / n - (sum / n) * (sum / n)));
    if (n < 31) {
      printf("%f\n", t_distribution[n - 2][1] * sqrt(sum2 / n - (sum / n) * (sum / n)));
    }
    return 0;
}
