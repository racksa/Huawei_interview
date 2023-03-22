#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define PI 3.14159265358979323846

// Recursive function to compute the DFT using SRT radix-2 algorithm
void dft(double complex *x, int N, double complex *X) {
    if (N == 1) {
        X[0] = x[0];
    } else {
        int k;
        double complex W_N = cexp(-2 * PI * I / N);
        double complex W = 1.0;
        double complex *X0 = (double complex*)malloc(N/2 * sizeof(double complex));
        double complex *X1 = (double complex*)malloc(N/2 * sizeof(double complex));
        double complex *x0 = (double complex*)malloc(N/2 * sizeof(double complex));
        double complex *x1 = (double complex*)malloc(N/2 * sizeof(double complex));
        for (k = 0; k < N/2; k++) {
            x0[k] = x[2*k];
            x1[k] = x[2*k+1];
        }
        dft(x0, N/2, X0);
        dft(x1, N/2, X1);
        for (k = 0; k < N/2; k++) {
            X[k] = X0[k] + W * X1[k];
            X[k + N/2] = X0[k] - W * X1[k];
            W = W * W_N;
        }
        free(X0);
        free(X1);
        free(x0);
        free(x1);
    }
}

int main() {
    int N = 8;  // length of input sequence
    double complex x[N];  // input sequence
    double complex X[N];  // DFT of input sequence
    int k;
    for (k = 0; k < N; k++) {
        x[k] = sin(2 * PI * k / N) + cos(4 * PI * k / N) * I;
    }
    dft(x, N, X);
    printf("DFT of input sequence:\n");
    for (k = 0; k < N; k++) {
        printf("%f + %fi\n", creal(X[k]), cimag(X[k]));
    }
    return 0;
}