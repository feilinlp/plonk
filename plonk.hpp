#ifndef PLONK_HPP
#define PLONK_HPP

#include <mcl/bn.hpp>
#include <vector>
#include <map>

using namespace std;
using namespace mcl;
using namespace bn;

class Plonk {
public:
    struct Circuit {
        int n;
        vector<Fr> qm, ql, qr, qo, qc;
        vector<Fr> h = {1};
        vector<Fr> permutation = {0}; // o*
        Fr k1, k2;
    };

    struct Preprocess {
        vector<G1> srs;
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

    struct Witness {
        G1 ax, bx, cx;
        G1 zx;
        G1 lox, midx, hix;
        G1 wx, wwx;
        Fr av, bv, cv, so1v, so2v, zwv;
    };

    static Circuit initialize(int n, vector<Fr> qm, vector<Fr> ql, vector<Fr> qr, vector<Fr> qo, vector<Fr> qc, vector<int> permutation);
    static Preprocess setup(Fr x, Circuit circuit);
    static Witness prove(Preprocess preprocess, int l, vector<Fr> w);
    // static bool verify();
};

#endif // PLONK_HPP
