#ifndef PLONK_HPP
#define PLONK_HPP

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
        vector<G1> srs;
        vector<G2> srs2;
        Circuit circuit;
        vector<vector<Fr>> lagrange;

        vector<Fr> qm, ql, qr, qo, qc;
        vector<Fr> so1, so2, so3;
    };

    struct Transcript {
        // ROUND 1
        vector<Fr> a, b, c;
        G1 ax, bx, cx;

        // ROUND 2
        vector<Fr> z;
        G1 zx;

        // ROUND 3
        vector<Fr> t, lo, mid, hi, pi;
        G1 lox, midx, hix;

        // ROUND 4
        Fr av, bv, cv, so1v, so2v, zwv;

        // ROUND 5
        vector<Fr> r, w, ww;
        G1 wx, wwx;

    };

    struct Challenge {
        Fr beta, gamma, alpha, v, u, e;
    };

    struct Witness {
        G1 ax, bx, cx;
        G1 zx;
        G1 lox, midx, hix;
        G1 wx, wwx;
        Fr av, bv, cv, so1v, so2v, zwv;
    };

    struct Verifier {
        G1 qmx, qlx, qrx, qox, qcx;
        G1 so1x, so2x, so3x, x;
    };

    static Circuit initialize(size_t n, vector<Fr> qm, vector<Fr> ql, vector<Fr> qr, vector<Fr> qo, vector<Fr> qc, vector<int> permutation);
    static Preprocess setup(Fr x, Circuit circuit);
    static Witness prove(Preprocess preprocess, size_t l, vector<Fr> w);
    static Verifier preprocess(Preprocess preprocess);
    static bool verify(size_t n, Preprocess preprocess, Witness witness, size_t l, vector<Fr> w);
};

#endif // PLONK_HPP
