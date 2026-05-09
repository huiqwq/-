#include <iostream>
using namespace std;

const int MAX_N = 2050;
double s[MAX_N][MAX_N];
double v[MAX_N];
double sum[MAX_N];
__attribute__((noinline))
void f_normal(int n) {
    for (int i = 0; i < n; i++) sum[i] = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            sum[i] += s[j][i] * v[j];
}

__attribute__((noinline))
void f_cache(int n) {
    for (int i = 0; i < n; i++) sum[i] = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            sum[j] += s[i][j] * v[i];
}

__attribute__((noinline))
void f_unroll(int n) {
    for (int i = 0; i < n; i++) sum[i] = 0;
    for (int i = 0; i < n; i++) {
        int j;
        for (j = 0; j <= n - 4; j += 4) {
            sum[j]   += s[i][j]   * v[i];
            sum[j+1] += s[i][j+1] * v[i];
            sum[j+2] += s[i][j+2] * v[i];
            sum[j+3] += s[i][j+3] * v[i];
        }
        for (; j < n; j++)
            sum[j] += s[i][j] * v[i];
    }
}

void init() {
    for (int i = 0; i < MAX_N; i++) {
        for (int j = 0; j < MAX_N; j++)
            s[i][j] = i + j;
        v[i] = i;
    }
}

int main() {
    init();

    int sizes[] = {128,512,1024,2048};

    for (int n : sizes) {
        for (int t = 0; t < 1000; t++) f_normal(n);
        for (int t = 0; t < 1000; t++) f_cache(n);
        for (int t = 0; t < 1000; t++) f_unroll(n);

        cout << sum[0] << endl;
    }

    return 0;
}