#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int planet, epoch;
    double time, rsky, vsky;
} TimesRow;

typedef struct {
    double time, rv;
} RVRow;

/* Copy of numpy.isclose */
int double_close(double x, double y) {
    const double rtol = 1E-5;
    const double atol = 1E-8;

    if (fabs(x - y) <= (atol + rtol * fabs(y))) {
        return 1;
    } else {
        return 0;
    }
}

void check_times_rows(const TimesRow *received, const TimesRow *expected, int line_no) {
    if (received->planet != expected->planet) {
        fprintf(stderr, "Planets do not match: %d != %d; line %d\n",
                received->planet, expected->planet, line_no);
        exit(1);
    }

    if (received->epoch != expected->epoch) {
        fprintf(stderr, "Epochs do not match: %d != %d; line %d\n",
                received->epoch, expected->epoch, line_no);
        exit(1);
    }
    
    if (!double_close(received->time, expected->time)) {
        fprintf(stderr, "Times do not match: %lf != %lf; line %d\n",
                received->time, expected->time, line_no);
        exit(1);
    }

    if (!double_close(received->rsky, expected->rsky)) {
        fprintf(stderr, "Rsky does not match: %lf != %lf; line %d\n",
                received->rsky, expected->rsky, line_no);
        exit(1);
    }

    if (!double_close(received->vsky, expected->vsky)) {
        fprintf(stderr, "Vsky does not match: %lf != %lf; line %d\n",
                received->vsky, expected->vsky, line_no);
        exit(1);
    }
}

void check_rv_rows(const RVRow *received, const RVRow *expected, int line_no) {
    if (!double_close(received->time, expected->time)) {
        fprintf(stderr, "Time does not match: %lf != %lf; line %d\n",
                received->time, expected->time, line_no);
        exit(1);
    }

    if (!double_close(received->rv, expected->rv)) {
        fprintf(stderr, "RV does not match: %lf != %lf; line %d\n",
                received->rv, expected->rv, line_no);
        exit(1);
    }
}

int main() {
    const char *times_file = "Times";
    const char *times_expected_file = "Times.expected";
    const char *rv_file = "RV_out";
    const char *rv_expected_file = "RV_out.expected";

    FILE *times, *times_expected;
    FILE *rv, *rv_expected;

    {
        printf("Checking %s\n", times_file);
        times = fopen(times_file, "r");
        times_expected = fopen(times_expected_file, "r");

        int line_no = 0;
        TimesRow received, expected;
        while (fscanf(times, "%d %d %lf %lf %lf\n",
                    &(received.planet), &(received.epoch),
                    &(received.time), &(received.rsky), &(received.vsky)) != EOF) {
            fscanf(times_expected, "%d %d %lf %lf %lf\n",
                    &(expected.planet), &(expected.epoch),
                    &(expected.time), &(expected.rsky), &(expected.vsky));

            check_times_rows(&received, &expected, line_no);

            line_no++;

        }
        fclose(times);
        fclose(times_expected);
    }

    {
        printf("Checking %s\n", rv_file);
        rv = fopen(rv_file, "r");
        rv_expected = fopen(rv_expected_file, "r");

        int line_no = 0;
        RVRow received, expected;
        while (fscanf(rv, "%lf %lf\n", &(received.time), &(received.rv)) != EOF) {
            fscanf(rv_expected, "%lf %lf\n", &(expected.time), &(expected.rv));

            check_rv_rows(&received, &expected, line_no);
            line_no++;
        }

        fclose(rv);
        fclose(rv_expected);
    }

    printf("OK\n");


    return 0;
}
