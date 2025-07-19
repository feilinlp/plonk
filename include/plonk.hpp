#ifndef PLONK_HPP
#define PLONK_HPP

#include "kzg.hpp"
#include <mcl/bn.hpp>
#include <vector>
#include <map>

using namespace std;
using namespace mcl;
using namespace bn;

G1 evaluatePolynomial(vector<Fr> &p, vector<G1> &x);
Fr evaluatePolynomial(vector<Fr> p, vector<Fr> x);
vector<vector<Fr>> lagrangeBasis(size_t n, Fr omega);

class Plonk {
public:
    struct Circuit {
        size_t n;
        vector<Fr> qm, ql, qr, qo, qc;
        vector<Fr> h;
        vector<Fr> permutation = {0}; // o*
        Fr k1, k2;
    };

    struct Preprocess {
        KZG::PublicKey pk;
        Circuit circuit;
        vector<vector<Fr>> lagrange;

        vector<Fr> qm, ql, qr, qo, qc;
        vector<Fr> so1, so2, so3;
    };

    struct Transcript {
        // ROUND 1
        vector<Fr> a, b, c;
        KZG::Commitment ax, bx, cx;

        // ROUND 2
        vector<Fr> z;
        KZG::Commitment zx;

        // ROUND 3
        vector<Fr> t, lo, mid, hi, pi;
        KZG::Commitment lox, midx, hix;

        // ROUND 4
        KZG::Witness av, bv, cv, so1v, so2v, zwv;

        // ROUND 5
        vector<Fr> r, w, ww;
        KZG::Commitment wx, wwx;
    };

    struct Challenge {
        Fr beta, gamma, alpha, v, u, e;
    };

    struct Witness {
        KZG::Commitment ax, bx, cx;
        KZG::Commitment zx;
        KZG::Commitment lox, midx, hix;
        KZG::Commitment wx, wwx;
        KZG::Witness av, bv, cv, so1v, so2v, zwv;
    };

    struct Verifier {
        KZG::Commitment qmx, qlx, qrx, qox, qcx;
        KZG::Commitment so1x, so2x, so3x;
    };

    static Circuit initialize(size_t n, vector<Fr> qm, vector<Fr> ql, vector<Fr> qr, vector<Fr> qo, vector<Fr> qc, vector<int> permutation);
    static Preprocess setup(Fr x, Circuit circuit);
    static Witness prove(Preprocess preprocess, size_t l, vector<Fr> w);
    static Verifier preprocess(Preprocess preprocess);
    static bool verify(size_t n, Preprocess preprocess, Witness witness, size_t l, vector<Fr> w);
};

#endif // PLONK_HPP
