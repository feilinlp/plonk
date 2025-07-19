#include "../include/plonk.hpp"
#include "../include/ntt.hpp"
#include "../include/kzg.hpp"
#include "../include/debug.hpp"
#include <openssl/evp.h>
#include <sstream>
#include <mcl/bn.hpp>
#include <vector>
#include <map>
#include <cmath>

using namespace std;
using namespace mcl;
using namespace bn;

vector<uint8_t> toBytes(const Fr &x) {
    string s = x.getStr(16); 
    return vector<uint8_t>(s.begin(), s.end());
}

vector<uint8_t> toBytes(const G1 &g) {
    string s;
    s = g.getStr();
    return vector<uint8_t>(s.begin(), s.end());
}

Fr hashToFr(const vector<vector<uint8_t>>& inputs) {
    // Create and initialize the digest context
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) {
        throw runtime_error("Failed to create EVP_MD_CTX");
    }
    
    // Initialize the digest operation with SHA256
    if (EVP_DigestInit_ex(ctx, EVP_sha256(), nullptr) != 1) {
        EVP_MD_CTX_free(ctx);
        throw runtime_error("Failed to initialize SHA256 digest");
    }
    
    // Update the digest with each input
    for (const auto& input : inputs) {
        if (EVP_DigestUpdate(ctx, input.data(), input.size()) != 1) {
            EVP_MD_CTX_free(ctx);
            throw runtime_error("Failed to update digest");
        }
    }
    
    // Finalize the digest
    uint8_t hash[EVP_MAX_MD_SIZE];
    unsigned int hash_len;
    if (EVP_DigestFinal_ex(ctx, hash, &hash_len) != 1) {
        EVP_MD_CTX_free(ctx);
        throw runtime_error("Failed to finalize digest");
    }
    
    // Clean up the context
    EVP_MD_CTX_free(ctx);
    
    // Create and return the Fr result
    Fr result;
    result.setArrayMask(hash, hash_len); 
    return result;
}

bool isQuadraticResidue(Fr w) {
    Fr power = -1;
    Fr::pow(w, w, power / 2);

    return w == 1;
}

Fr quadraticNonResidue() {
    Fr ret;
    ret.setByCSPRNG();

    while (isQuadraticResidue(ret)) { ret.setByCSPRNG(); }
    return ret;
}

bool inCoset(vector<Fr> h, Fr k2) {
    size_t n = h.size() / 2;
    for (size_t i = 0; i < n; i++) {
        if (h[n + 1 + i] == k2) return true;
    }
    return false;
} 

vector<vector<Fr>> lagrangeBasis(size_t n, Fr omega) {
    vector<vector<Fr>> lagrange(n + 1);
    for (size_t i = 0; i < n; ++i) {
        vector<Fr> e(n, Fr(0));
        e[i] = Fr(1);  
        
        ntt_inverse(e, omega); 
        lagrange[i] = e;
    }
    lagrange[n] = lagrange[0];
    return lagrange;
}

Fr evaluatePolynomial(vector<Fr> p, vector<Fr> x) {
    Fr ret = 0;
    for (size_t i = 0; i < p.size(); i++) { ret += p[i] * x[i]; }
    return ret;
}

vector<Fr> polynomialDivision(vector<Fr> a, size_t n) {
    // Remove leading zeros
    while (!a.empty() && a.back().isZero()) {
        a.pop_back();
    }
    
    if (a.size() <= n) {
        return vector<Fr>(1, Fr(0));
    }
    
    vector<Fr> quotient(a.size() - n, Fr(0));
    
    // Copy high-degree coefficients to quotient and add them to low-degree terms
    for (size_t i = a.size() - 1; i >= n; i--) {
        quotient[i - n] = a[i];
        a[i - n] += a[i]; 
        a[i] = 0;
    }
    
    return quotient;
}

vector<Fr> addPolynomials(vector<Fr> a, vector<Fr> b) {
    vector<Fr> ret(max(a.size(), b.size()), 0);

    for (size_t i = 0; i < a.size(); i++) ret[i] += a[i];
    for (size_t i = 0; i < b.size(); i++) ret[i] += b[i];

    return ret;
}
