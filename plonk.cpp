#include "plonk.hpp"
#include "ntt.hpp"
#include "kzg.hpp"
#include <mcl/bn.hpp>
#include <vector>
#include <map>
#include <cmath>

using namespace std;
using namespace mcl;
using namespace bn;

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
    int n = h.size() / 2;
    for (int i = 0; i < n; i++) {
        if (h[n + 1 + i] == k2) return true;
    }
    return false;
} 

vector<vector<Fr>> lagrangeBasis(int n, Fr omega) {
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

G1 evaluatePolynomial(vector<Fr> &p, vector<G1> &x) {
    G1 ret;

    G1::mulVec(ret, x.data(), p.data(), p.size());

    return ret;
}

Fr evaluatePolynomial(vector<Fr> p, vector<Fr> x) {
    Fr ret = 0;
    for (int i = 0; i < p.size(); i++) { ret += p[i] * x[i]; }
    return ret;
}


// All vectors are 1-indexed
Plonk::Circuit Plonk::initialize(int n, vector<Fr> qm, vector<Fr> ql, vector<Fr> qr, vector<Fr> qo, vector<Fr> qc, vector<int> permutation) {
    Plonk::Circuit ret;

    ret.n = n;
    ret.qm = qm;
    ret.ql = ql;
    ret.qr = qr;
    ret.qo = qo;
    ret.qc = qc;

    Fr w = findPrimitiveRoot(n);
    if (!isQuadraticResidue(w)) throw runtime_error("Not quadratic residue!");

    Fr curr = w;
    for (int i = 0; i < n; i++) {
        ret.h.push_back(curr);
        Fr::mul(curr, curr, w);
    }

    Fr k1 = quadraticNonResidue();
    ret.k1 = k1;
    k1 *= w;
    for (int i = 0; i < n; i++) {
        ret.h.push_back(k1);
        Fr::mul(k1, k1, w);
    }

    Fr k2 = quadraticNonResidue();
    while (inCoset(ret.h, k2)) { k2 = quadraticNonResidue(); }
    ret.k2 = k2;
    Fr::mul(k2, k2, w);
    for (int i = 0; i < n; i++) {
        ret.h.push_back(k2);
        Fr::mul(k2, k2, w);
    }

    for (int i = 1; i <= 3 * n; i++) {
        ret.permutation.push_back(ret.h[permutation[i]]);
    }

    return ret;
}

Plonk::Preprocess Plonk::setup(Fr x, Plonk::Circuit circuit) {
    Plonk::Preprocess ret;
    ret.lagrange = lagrangeBasis(circuit.n, circuit.h[1]);
    ret.circuit = circuit;

    ret.qm.resize(circuit.n, 0);
    ret.ql.resize(circuit.n, 0);
    ret.qr.resize(circuit.n, 0);
    ret.qo.resize(circuit.n, 0);
    ret.qc.resize(circuit.n, 0);
    ret.so1.resize(circuit.n, 0);
    ret.so2.resize(circuit.n, 0);
    ret.so3.resize(circuit.n, 0);

    // Create the SRS
    G1 one;
    mapToG1(one, 1);

    Fr curr = 1;
    for (int i = 0; i <= circuit.n + 5; i++) { 
        G1 ins;
        G1::mul(ins, one, curr);  
        
        ret.srs.push_back(ins);
        
        Fr::mul(curr, curr, x);
    }

    for (int i = 1; i <= circuit.n; i++) {
        for (int j = 0; j < ret.lagrange[i].size(); j++) {
            ret.qm[j] += circuit.qm[i] * ret.lagrange[i][j];
            ret.ql[j] += circuit.ql[i] * ret.lagrange[i][j];
            ret.qr[j] += circuit.qr[i] * ret.lagrange[i][j];
            ret.qo[j] += circuit.qo[i] * ret.lagrange[i][j];
            ret.qc[j] += circuit.qc[i] * ret.lagrange[i][j];

            ret.so1[j] += circuit.permutation[i] * ret.lagrange[i][j];
            ret.so2[j] += circuit.permutation[circuit.n + i] * ret.lagrange[i][j];
            ret.so3[j] += circuit.permutation[2 * circuit.n + i] * ret.lagrange[i][j];
        }
    }

    return ret;
}

vector<Fr> abc(int n, Fr b1, Fr b2, Fr* w, vector<vector<Fr>> lagrange) {
    vector<Fr> ret(n + 2, 0);
    
    // b1 * X^(n+1) - b1 * X + b2 * X^n - b2
    ret[n+1] = b1;
    ret[1] = -b1;
    
    ret[n] = b2;
    ret[0] = -b2;
    
    for (int i = 1; i <= n; i++) {
        for (int j = 0; j <= lagrange[i].size(); j++) {
            ret[j] += w[i] * lagrange[i][j];
        }
    }

    return ret;
}

void roundOne(Plonk::Transcript &ret, Plonk::Preprocess preprocess, vector<Fr> w, vector<Fr> b) {
    // Alternatively, pass original vector using data() and add offset.
    int n = preprocess.circuit.n;
    ret.a = abc(n, b[1], b[2], &w[0], preprocess.lagrange);
    ret.b = abc(n, b[3], b[4], &w[n], preprocess.lagrange);
    ret.c = abc(n, b[5], b[6], &w[2 * n], preprocess.lagrange);

    ret.ax = evaluatePolynomial(ret.a, preprocess.srs);
    ret.bx = evaluatePolynomial(ret.b, preprocess.srs);
    ret.cx = evaluatePolynomial(ret.c, preprocess.srs);
}

void roundTwo(Plonk::Transcript &ret, Plonk::Preprocess preprocess, vector<Fr> w, vector<Fr> b, Fr beta, Fr gamma) { 
    int n = preprocess.circuit.n;
    ret.z.resize(n + 3, 0);
    ret.z[n+2] = b[7];
    ret.z[n+1] = b[8];
    ret.z[n] = b[9];

    ret.z[2] = -b[7];
    ret.z[1] = -b[8];
    ret.z[0] = -b[9];

    for (int i = 0; i < preprocess.lagrange[1].size(); i++) {
        ret.z[i] += preprocess.lagrange[1][i];
    }

    Fr ratio = 1;
    for (int i = 1; i < n; i++) {
        for (int j = i; j < 3 * n; j += n) {
            Fr temp = gamma, _;
            Fr::add(temp, temp, w[j]);

            Fr::mul(_, beta, preprocess.circuit.h[j]);
            Fr::add(temp, temp, _);

            Fr::mul(ratio, ratio, temp);

            temp = gamma;
            Fr::add(temp, temp, w[j]);

            Fr::mul(_, beta, preprocess.circuit.permutation[j]);
            Fr::add(temp, temp, _);

            Fr::div(ratio, ratio, temp);
        }

        for (int j = 0; j < preprocess.lagrange[i+1].size(); j++) {
            ret.z[j] += preprocess.lagrange[i+1][j] * ratio;
        }
    }

    ret.zx = evaluatePolynomial(ret.z, preprocess.srs);
}

// Divide by ZH
vector<Fr> polynomialDivision(vector<Fr> &a, size_t n) {
    // Remove leading zeros
    while (!a.empty() && a.back().isZero()) {
        a.pop_back();
    }
    
    if (a.size() <= n) {
        return vector<Fr>(1, Fr(0));
    }
    
    vector<Fr> quotient(a.size() - n, Fr(0));
    
    // Copy high-degree coefficients to quotient and add them to low-degree terms
    for (int i = a.size() - 1; i >= n; i--) {
        quotient[i - n] = a[i];
        a[i - n] += a[i]; 
        a[i] = 0;
    }
    
    return quotient;
}

vector<Fr> addPolynomials(vector<Fr> a, vector<Fr> b) {
    vector<Fr> ret(max(a.size(), b.size()), 0);

    for (int i = 0; i < a.size(); i++) ret[i] += a[i];
    for (int i = 0; i < b.size(); i++) ret[i] += b[i];

    return ret;
}

void split(Plonk::Transcript &ret, int n, vector<Fr> b) {
    ret.lo.resize(n + 1, 0);
    ret.mid.resize(n + 1, 0);
    ret.hi.resize(n + 6, 0);

    for (int i = 0; i < n; i++) {
        ret.lo[i] = ret.t[i];
        ret.mid[i] = ret.t[i + n];
    }

    for (int i = 2 * n; i < ret.t.size(); i++) { ret.hi[i - 2 * n] = ret.t[i]; }

    ret.lo[n] += b[10];
    ret.mid[0] -= b[10];
    ret.mid[n] += b[11];
    ret.hi[0] -= b[11];
}

void roundThree(Plonk::Transcript &ret, Plonk::Preprocess preprocess, int l, vector<Fr> w, vector<Fr> b, Fr alpha, Fr beta, Fr gamma) {
    int n = preprocess.circuit.n;
    Fr omega = preprocess.circuit.h[1];
    
    ret.t.resize(3 * n, 0);

    // Constraint
    vector<Fr> temp, _;
    temp = polynomial_multiply(ret.a, ret.b, omega);
    temp = polynomial_multiply(temp, preprocess.qm, omega);

    _ = polynomial_multiply(ret.a, preprocess.ql, omega);
    temp = addPolynomials(temp, _);

    _.clear();
    _ = polynomial_multiply(ret.b, preprocess.qr, omega);
    temp = addPolynomials(temp, _);

    _.clear();
    _ = polynomial_multiply(ret.c, preprocess.qo, omega);
    temp = addPolynomials(temp, _);

    temp = addPolynomials(temp, preprocess.qc);

    vector<Fr> pi(n + 1, 0);

    for (int i = 1; i <= l; i++) {
        for (int j = 0; j < preprocess.lagrange[i].size(); j++) {
            pi[j] += preprocess.lagrange[i][j] * (-1) * w[i];
        }
    }

    ret.pi = pi;

    temp = addPolynomials(temp, pi);

    ret.t = addPolynomials(ret.t, temp);

    // Checks for Round 2 polynomials
    // First part
    temp.clear();
    temp.assign(ret.a.begin(), ret.a.end());
    temp[1] += beta; temp[0] += gamma;

    _.clear();
    _.assign(b.begin(), b.end());
    _[1] += beta * preprocess.circuit.k1; _[0] += gamma;
    temp = polynomial_multiply(temp, _, omega);

    _.clear();
    _.assign(ret.c.begin(), ret.c.end());
    _[1] += beta * preprocess.circuit.k2; _[0] += gamma;
    temp = polynomial_multiply(temp, _, omega);

    temp = polynomial_multiply(temp, ret.z, omega);

    for (int i = 0; i < temp.size(); i++) { temp[i] *= alpha; }

    ret.t = addPolynomials(ret.t, temp);

    // Second part
    vector<Fr> zw;
    zw.assign(ret.z.begin(), ret.z.end());

    vector<Fr> so1, so2, so3;
    so1.assign(preprocess.so1.begin(), preprocess.so1.end());
    so2.assign(preprocess.so2.begin(), preprocess.so2.end());
    so3.assign(preprocess.so3.begin(), preprocess.so3.end());

    for (int i = 0; i < so1.size(); i++) { so1[i] *= beta; }
    for (int i = 0; i < so2.size(); i++) { so2[i] *= beta; }
    for (int i = 0; i < so3.size(); i++) { so3[i] *= beta; }

    // Calculate z(Xw)
    for (int i = 0; i < zw.size(); i++) { zw[i] *= preprocess.circuit.h[i]; }

    temp.clear();
    temp.assign(ret.a.begin(), ret.a.end());
    for (int i = 0; i < so1.size(); i++) { temp[i] += so1[i]; }
    temp[0] += gamma;

    _.clear();
    _.assign(b.begin(), b.end());
    for (int i = 0; i < so2.size(); i++) { _[i] += so2[i]; }
    _[0] += gamma;
    temp = polynomial_multiply(temp, _, omega);

    _.clear();
    _.assign(b.begin(), b.end());
    for (int i = 0; i < so3.size(); i++) { _[i] += so3[i]; }
    _[0] += gamma;
    temp = polynomial_multiply(temp, _, omega);

    temp = polynomial_multiply(temp, zw, omega);

    for (int i = 0; i < temp.size(); i++) { temp[i] *= -alpha; }

    ret.t = addPolynomials(ret.t, temp);

    // Induction base case check
    ret.z[0] -= 1;
    temp.clear();
    temp = polynomial_multiply(ret.z, preprocess.lagrange[1], omega);
    for (int i = 0; i < temp.size(); i++) { temp[i] *= alpha * alpha; }
    ret.z[0] += 1;

    ret.t = addPolynomials(ret.t, temp);
    ret.t = polynomialDivision(ret.t, n); // Divide by ZH

    split(ret, n, b);

    ret.lox = evaluatePolynomial(ret.lo, preprocess.srs);
    ret.midx = evaluatePolynomial(ret.mid, preprocess.srs);
    ret.hix = evaluatePolynomial(ret.hi, preprocess.srs);

    vector<Fr>().swap(zw);
    vector<Fr>().swap(temp);
    vector<Fr>().swap(_);
    vector<Fr>().swap(so1);
    vector<Fr>().swap(so2);
    vector<Fr>().swap(so3);
}

void roundFour(Plonk::Transcript &ret, Plonk::Preprocess preprocess, vector<Fr> vs) {
    ret.av = evaluatePolynomial(ret.a, vs);
    ret.bv = evaluatePolynomial(ret.b, vs);
    ret.cv = evaluatePolynomial(ret.c, vs);

    ret.so1v = evaluatePolynomial(preprocess.so1, vs);
    ret.so2v = evaluatePolynomial(preprocess.so2, vs);

    Fr curr = 1;
    for (int i = 0; i < vs.size(); i++) {
        vs[i] *= curr;
        curr *= preprocess.circuit.h[1];
    }

    int count = 0;
    for (int i = 0; i < ret.z.size(); i++) { 
        if (ret.z[i] != 0) count++;
    }
    ret.zwv = evaluatePolynomial(ret.z, vs);
}

void roundFive(Plonk::Transcript &ret, Plonk::Preprocess preprocess, Fr alpha, Fr beta, Fr gamma, vector<Fr> vs) {
    // Evaluating r(X)
    int n = preprocess.circuit.n;
    ret.r.resize(3 * n, 0); 

    for (int i = 0; i < preprocess.qm.size(); i++) { ret.r[i] = ret.av * ret.bv * preprocess.qm[i]; }
    for (int i = 0; i < preprocess.ql.size(); i++) { ret.r[i] += ret.av * preprocess.ql[i]; }
    for (int i = 0; i < preprocess.qr.size(); i++) { ret.r[i] += ret.bv * preprocess.qr[i]; }
    for (int i = 0; i < preprocess.qo.size(); i++) { ret.r[i] += ret.cv * preprocess.qo[i]; }
    for (int i = 0; i < preprocess.qc.size(); i++) { ret.r[i] += preprocess.qc[i]; }

    Fr piv = evaluatePolynomial(ret.pi, vs);
    ret.r[0] += piv;

    for (int i = 0; i < ret.z.size(); i++) { 
        ret.r[i] += ret.z[i] * alpha * (ret.av + beta * vs[1] + gamma) * 
            (ret.bv + beta * vs[1] * preprocess.circuit.k1 + gamma) * 
            (ret.cv + beta * vs[1] * preprocess.circuit.k2 + gamma); 
    }

    vector<Fr> so3;
    so3.assign(preprocess.so3.begin(), preprocess.so3.end());
    for (int i = 0; i < so3.size(); i++) { so3[i] *= beta; }
    so3[0] += ret.cv + gamma; 

    for (int i = 0; i < so3.size(); i++) {
        ret.r[i] -= so3[i] * alpha * ret.zwv * (ret.av + beta * ret.so1v + gamma) * 
            (ret.bv + beta * ret.so2v + gamma);
    }

    Fr l1v = evaluatePolynomial(preprocess.lagrange[1], vs);
    for (int i = 0; i < ret.z.size(); i++) { ret.r[i] += alpha * alpha * ret.z[i] * l1v; }
    ret.r[0] -= alpha * alpha * l1v;

    Fr zhv = vs[n] - 1;
    for (int i = 0; i < ret.lo.size(); i++) { ret.r[i] -= ret.lo[i] * zhv; }
    for (int i = 0; i < ret.mid.size(); i++) { ret.r[i] -= ret.mid[i] * zhv * vs[n]; }
    for (int i = 0; i < ret.hi.size(); i++) { ret.r[i] -= ret.hi[i] * zhv * vs[n] * vs[n]; }

    // Evaluating W(X)
    Fr u;
    u.setByCSPRNG(); // Should be H(transcript)
    ret.w.assign(ret.r.begin(), ret.r.end());

    for (int i = 0; i < ret.a.size(); i++) { ret.w[i] += ret.a[i] * vs[1]; }
    for (int i = 0; i < ret.b.size(); i++) { ret.w[i] += ret.b[i] * vs[2]; }
    for (int i = 0; i < ret.c.size(); i++) { ret.w[i] += ret.c[i] * vs[3]; }
    for (int i = 0; i < preprocess.so1.size(); i++) { ret.w[i] += preprocess.so1[i] * vs[4]; }
    for (int i = 0; i < preprocess.so2.size(); i++) { ret.w[i] += preprocess.so2[i] * vs[5]; }

    ret.w[0] -= ret.av * vs[1] + ret.bv * vs[2] + ret.cv * vs[3] + ret.so1v * vs[4] + ret.so2v * vs[5];
    ret.w = divideByLinear(ret.w, vs[1]);
    ret.wx = evaluatePolynomial(ret.w, preprocess.srs);

    // Evaluating Ww(X)
    ret.ww.assign(ret.z.begin(), ret.z.end());
    ret.ww[0] -= ret.zwv;

    ret.ww = divideByLinear(ret.ww, vs[1] * preprocess.circuit.h[1]);
    ret.wwx = evaluatePolynomial(ret.ww, preprocess.srs);
}

// Can be optimized by not converting back-and-forth polynomials from lagrange basis to monomial basis
// Missing KZG commitments
Plonk::Witness Plonk::prove(Plonk::Preprocess preprocess, int l, vector<Fr> w) {
    Plonk::Transcript transcript;
    Plonk::Witness ret;

    vector<Fr> b(22); b[0] = 0;
    for (int i = 1; i <= 11; i++) {
        b[i].setByCSPRNG();
    }
    
    // ROUND 1
    roundOne(transcript, preprocess, w, b);

    // ROUND 2
    Fr beta, gamma; // Should be H(transcript, 0) and H(transcript, 1) respectively
    beta.setByCSPRNG();
    gamma.setByCSPRNG();
    roundTwo(transcript, preprocess, w, b, beta, gamma);

    // ROUND 3
    Fr alpha;
    alpha.setByCSPRNG(); // Should be H(transcript)
    roundThree(transcript, preprocess, l, w, b, alpha, beta, gamma);

    // ROUND 4
    Fr v; 
    v.setByCSPRNG(); // Should be H(transcript)
    vector<Fr> vs(preprocess.circuit.n + 5, 1);
    Fr curr = 1;
    for (int i = 0; i < vs.size(); i++) {
        vs[i] = curr;
        curr *= v;
    }
    roundFour(transcript, preprocess, vs);

    // ROUND 5
    roundFive(transcript, preprocess, alpha, beta, gamma, vs);

    // Copy values from transcript to witness
    ret.ax = transcript.ax;
    ret.bx = transcript.bx;
    ret.cx = transcript.cx;
    ret.zx = transcript.zx;

    ret.lox = transcript.lox;
    ret.midx = transcript.midx;
    ret.hix = transcript.hix;

    ret.wx = transcript.wx;
    ret.wwx = transcript.wwx;

    ret.av = transcript.av;
    ret.bv = transcript.bv;
    ret.cv = transcript.cv;

    ret.so1v = transcript.so1v;
    ret.so2v = transcript.so2v;
    ret.zwv = transcript.zwv;

    return ret;
}

bool verify() {
    // 
}
