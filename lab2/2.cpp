#include <bits/stdc++.h>
using namespace std;

const int N = 2048 * 2048;
long long arr[N];
long long c[N]; 
long long sum[N];
void init() {
    for (int i = 0; i < N; i++) {
        arr[i] = rand() % 100;
    }
}

long long sum_normal(int n) {
    long long sum = 0;
    for (int i = 0; i < n; i++) {
        sum += arr[i];
    }
    return sum;
}

long long sum_ilp(int n) {
    long long sum1 = 0, sum2 = 0;
    int i;
    for (i = 0; i + 1 < n; i += 2) {
        sum1 += arr[i];
        sum2 += arr[i + 1];
    }
    for (; i < n; i++) {
        sum1 += arr[i];
    }
    return sum1 + sum2;
}

void sum__recursive1(int n) {
    if (n == 1) return;
    for (int i = 0; i < n / 2; i++) {
        c[i] += c[n - i - 1];
    }
    sum__recursive1(n / 2);
}

long long sum__recursive2(int n) {
    for (long long m = n; m > 1; m /= 2) {
        for (long long i = 0; i < m / 2; i++) {
            c[i] = c[i * 2] + c[i * 2 + 1];
        }
    }
    return c[0];
}

int main() {
    init();
    int sizes[] = {2048*2048};
    int repeat = 100;
    for (int n : sizes) {
        for (int t = 0; t < repeat; t++)
            sum_normal(n);

        for (int t = 0; t < repeat; t++)
            sum_ilp(n);

        for (int t = 0; t < repeat; t++){
        	for (int i = 0; i < n; i++) c[i] = arr[i];
            sum__recursive1(n);
       	}
		cout<<c[0]<<endl;
        for (int t = 0; t < repeat; t++){
        	for (int i = 0; i < n; i++) c[i] = arr[i];
            sum__recursive2(n);
       	}
		cout<<c[0]<<endl;
        cout << sum_normal(n) << endl;
        cout << sum_ilp(n) << endl;
    }

    return 0;
}