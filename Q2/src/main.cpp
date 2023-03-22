#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;

// Recursive function to compute the DFT using SRT radix-2 algorithm
void dft(vector<complex<double>>& x, vector<complex<double>>& X) {
    int N = x.size();
    if (N == 1) {
        X[0] = x[0];
    } else {
        complex<double> W_N = exp(-2.0 * M_PI * complex<double>(0, 1) / double(N));
        complex<double> W = 1.0;
        vector<complex<double>> X0(N/2);
        vector<complex<double>> X1(N/2);
        vector<complex<double>> x0(N/2);
        vector<complex<double>> x1(N/2);
        for (int k = 0; k < N/2; k++) {
            x0[k] = x[2*k];
            x1[k] = x[2*k+1];
        }
        dft(x0, X0);
        dft(x1, X1);
        for (int k = 0; k < N/2; k++) {
            X[k] = X0[k] + W * X1[k];
            X[k + N/2] = X0[k] - W * X1[k];
            W = W * W_N;
        }
    }
}

int main() {
    int N = 7;  // length of input sequence
    int N_padded = pow(2, ceil(log2(N)));  // nearest power of 2 to N
    vector<complex<double>> x(N_padded);  // input sequence padded with zeros
    vector<complex<double>> X(N_padded);  // DFT of input sequence
    // Initialize input sequence
    for (int k = 0; k < N; k++) {
        x[k] = sin(2 * M_PI * k / N) + cos(4 * M_PI * k / N) * complex<double>(0, 1);
    }
    // Compute DFT of input sequence
    dft(x, X);
    // Print DFT of input sequence
    cout << "DFT of input sequence:" << endl;
    for (int k = 0; k < N_padded; k++) {
        cout << X[k] << endl;
    }
    return 0;
}