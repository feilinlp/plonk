#include <mcl/bn.hpp>
#include <vector>
#include <cmath>

using namespace mcl;
using namespace bn;
using namespace std;

// Hardcoded generator of Fr
const int GENERATOR = 5;

Fr findPrimitiveRoot(size_t N) {
    Fr order = -1;

    size_t i = 1;
    while (i < N / 2) {
        i *= 2;
        Fr::squareRoot(order, order);
    }

    return order;
}

Fr getRootOfUnity(size_t n) {
    mpz_class r("21888242871839275222246405745257275088548364400416034343698204186575808495617");
    mpz_class exp = (r - 1) / n;

    Fr omega;
    Fr::pow(omega, Fr(GENERATOR), exp);
    return omega;
}

size_t bitReverse(size_t x, size_t logN) {
    size_t res = 0;
    for (size_t i = 0; i < logN; ++i) {
        res <<= 1;
        res |= (x >> i) & 1;
    }
    return res;
}

void ntt_transform(vector<Fr> &A, Fr omega) {
    size_t n = A.size();
    size_t logN = log2(n);
    for (size_t i = 0; i < n; ++i) {
        size_t j = bitReverse(i, logN);
        if (i < j) swap(A[i], A[j]);
    }

    for (size_t len = 2; len <= n; len <<= 1) {
        Fr wlen;
        Fr::pow(wlen, omega, n / len);
        for (size_t i = 0; i < n; i += len) {
            Fr w = 1;
            for (size_t j = 0; j < len / 2; ++j) {
                Fr u = A[i + j];
                Fr v = A[i + j + len / 2] * w;
                A[i + j] = u + v;
                A[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

void ntt_inverse(vector<Fr> &A, Fr omega) {
    size_t n = A.size();
    Fr omega_inv;
    Fr::inv(omega_inv, omega);
    ntt_transform(A, omega_inv);

    Fr n_inv;
    Fr::inv(n_inv, Fr(n));
    for (Fr &x : A) x *= n_inv;
}

vector<Fr> polynomial_multiply(vector<Fr> A, vector<Fr> B) {
    size_t n = 1;
    while (n < A.size() + B.size() - 1) n *= 2;

    A.resize(n, Fr(0));
    B.resize(n, Fr(0));

    Fr omega = getRootOfUnity(n);

    ntt_transform(A, omega);
    ntt_transform(B, omega);

    vector<Fr> C(n);
    for (size_t i = 0; i < n; i++) {
        C[i] = A[i] * B[i];
    }

    ntt_inverse(C, omega);

    while (!C.empty() && C.back().isZero()) {
        C.pop_back();
    }

    return C;
}
