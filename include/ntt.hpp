#ifndef NTT_HPP
#define NTT_HPP

#include <mcl/bn.hpp>
#include <vector>
#include <cassert>
#include <cmath>

using namespace mcl;
using namespace bn;
using namespace std;

/**
 * @brief Performs bit reversal for NTT
 * @param x The number to reverse bits for
 * @param logN Log base 2 of the array size
 * @return Bit-reversed number
 */
size_t bitReverse(size_t x, size_t logN);

/**
 * @brief Finds a primitive N-th root of unity in the finite field
 * @param N The order of the root (must be power of 2)
 * @return A primitive N-th root of unity
 */
Fr findPrimitiveRoot(size_t N);

/**
 * @brief Performs Number Theoretic Transform (NTT) on the input array
 * @param A Input/output vector of field elements (size must be power of 2)
 * @param omega Primitive N-th root of unity where N = A.size()
 * 
 * The function performs in-place NTT transformation.
 * After calling this function, A[i] will contain the i-th coefficient
 * of the NTT of the original polynomial.
 */
void ntt_transform(vector<Fr> &A, Fr omega);

/**
 * @brief Performs inverse Number Theoretic Transform (INTT)
 * @param A Input/output vector of field elements (size must be power of 2)
 * @param omega Primitive N-th root of unity where N = A.size()
 * 
 * This function undoes the NTT transformation by using the inverse
 * of omega and scaling by 1/N.
 */
void ntt_inverse(vector<Fr> &A, Fr omega);

/**
 * @brief Multiplies two polynomials using NTT
 * @param A First polynomial coefficients
 * @param B Second polynomial coefficients
 * @return Product polynomial coefficients
 * 
 * The result size will be a.size() + b.size() - 1.
 * Input vectors will be zero-padded to the next power of 2 if needed.
 */
vector<Fr> polynomial_multiply(vector<Fr> A, vector<Fr> B);

Fr getRootOfUnity(size_t n);

// Direct polynomial multiplication in NTT domain
void polynomial_multiply_ntt(vector<Fr> &A_ntt, const vector<Fr> &B_ntt);

// Add polynomials in NTT domain
void polynomial_add_ntt(vector<Fr> &A_ntt, const vector<Fr> &B_ntt);

// Scale polynomial in NTT domain
void polynomial_scale_ntt(vector<Fr> &A_ntt, const Fr &scalar);

// Convert polynomial to NTT domain with proper padding
vector<Fr> to_ntt_domain(const vector<Fr> &poly, size_t target_size, Fr omega);

#endif // NTT_HPP
