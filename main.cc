#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <cstdint>
// 可以自行添加需要的头文件

void fRead(int *a, int *b, int *n, int *p, int input_id){
    // 数据输入函数
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin>>*n>>*p;
    for (int i = 0; i < *n; i++){
        fin>>a[i];
    }
    for (int i = 0; i < *n; i++){   
        fin>>b[i];
    }
}

void fCheck(int *ab, int n, int input_id){
    // 判断多项式乘法结果是否正确
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++){
        int x;
        fin>>x;
        if(x != ab[i]){
            std::cout<<"多项式乘法结果错误"<<std::endl;
            return;
        }
    }
    std::cout<<"多项式乘法结果正确"<<std::endl;
    return;
}

void fWrite(int *ab, int n, int input_id){
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++){
        fout<<ab[i]<<'\n';
    }
}

void poly_multiply(int *a, int *b, int *ab, int n, int p){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            ab[i+j]=(1LL * a[i] * b[j] % p + ab[i+j]) % p;
        }
    }
}

long long qpow(long long a, long long b, int p) {
    long long res = 1;
    while (b) {
        if (b & 1) res = res * a % p;
        a = a * a % p;
        b >>= 1;
    }
    return res;
}

struct Montgomery32{
    using u32 = uint32_t;
    using u64 = uint64_t;

    u32 mod;
    u32 inv;
    u32 r2;
    explicit Montgomery32(u32 p) : mod(p) {
        u32 x = 1;
        for (int i = 0; i < 5; ++i) {
            x *= 2 - mod * x;
        }
        inv = 0 - x;

        u64 r = (u64(1) << 32) % mod;
        r2 = u32(r * r % mod);
    }
    u32 reduce(u64 x) const {
    //    u32 q = u32(x) * u64(inv);
        u32 q = (u32)((x & 0xffffffffULL) * inv);
        u64 y = (x + u64(q) * mod) >> 32;
        if (y >= mod) y -= mod;
        return u32(y);
    }

    u32 init(u32 x) const {
        return reduce(u64(x) * r2);
    }
        u32 value(u32 x) const {
        return reduce(x);
    }

    u32 one() const {
        return init(1);
    }

    u32 add(u32 a, u32 b) const {
        u32 c = a + b;
        if (c >= mod) c -= mod;
        return c;
    }

    u32 sub(u32 a, u32 b) const {
        return a >= b ? a - b : a + mod - b;
    }

    u32 mul(u32 a, u32 b) const {
        return reduce(u64(a) * b);
    }
};

void ntt(std::vector<int>& f, int p, bool inverse) {
    int n = f.size();
    Montgomery32 mont((uint32_t)p);

    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) std::swap(f[i], f[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        int wn_normal = (int)qpow(3, (p - 1) / len, p);
        if (inverse) {
            wn_normal = (int)qpow(wn_normal, p - 2, p);
        }
        uint32_t wn = mont.init((uint32_t)wn_normal);

        for (int i = 0; i < n; i += len) {
            uint32_t w = mont.one();
            for (int j = 0; j < len / 2; ++j) {
                uint32_t u = (uint32_t)f[i + j];
                uint32_t v = mont.mul(w, (uint32_t)f[i + j + len / 2]);

                f[i + j] = (int)mont.add(u, v);
                f[i + j + len / 2] = (int)mont.sub(u, v);

                w = mont.mul(w, wn);
            }
        }
    }

    if (inverse) {
        uint32_t inv_n = mont.init((uint32_t)qpow(n, p - 2, p));
        for (int i = 0; i < n; ++i) {
            f[i] = (int)mont.mul((uint32_t)f[i], inv_n);
        }
    }
}

void poly_multiply_ntt(int *a, int *b, int *ab, int n, int p) {
    int lim = 1;
    while (lim < 2 * n - 1) lim <<= 1;
    
    Montgomery32 mont((uint32_t)p);
    std::vector<int> A(lim), B(lim);

    for (int i = 0; i < n; ++i) {
        A[i] = (int)mont.init((uint32_t)a[i]);
        B[i] = (int)mont.init((uint32_t)b[i]);
    }

    ntt(A, p, false);
    ntt(B, p, false);

    for (int i = 0; i < lim; ++i) {
        A[i] = (int)mont.mul((uint32_t)A[i], (uint32_t)B[i]);
    }

    ntt(A, p, true);

    for (int i = 0; i < 2 * n - 1; ++i) {
        ab[i] = (int)mont.value((uint32_t)A[i]);
    }
}


int a[300000], b[300000], ab[300000];
int main(int argc, char *argv[])
{
    
    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 263882790666241
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    int test_begin = 0;
    int test_end = 3;
    for(int i = test_begin; i <= test_end; ++i){
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));
        auto Start = std::chrono::high_resolution_clock::now();
        // TODO : 将 poly_multiply 函数替换成你写的 ntt
        //poly_multiply(a, b, ab, n_, p_);
        poly_multiply_ntt(a,b,ab,n_,p_);
        auto End = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
        ans += elapsed.count();
        fCheck(ab, n_, i);
        std::cout<<"average latency for n = "<<n_<<" p = "<<p_<<" : "<<ans<<" (us) "<<std::endl;
        // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        // 禁止使用 cout 一次性输出大量文件内容
        fWrite(ab, n_, i);
    }
    return 0;
}
 
