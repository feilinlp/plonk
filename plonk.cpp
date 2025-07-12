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

G1 evaluatePolynomialG1(vector<Fr> p, vector<G1> x) {
    G1 ret;

    G1::mulVec(ret, x.data(), p.data(), p.size());

    return ret;
}

Fr evaluatePolynomial(vector<Fr> p, vector<Fr> x) {
    Fr ret = 0;
    for (size_t i = 0; i < p.size(); i++) { ret += p[i] * x[i]; }
    return ret;
}

// All vectors are 1-indexed
Plonk::Circuit Plonk::initialize(size_t n, vector<Fr> qm, vector<Fr> ql, vector<Fr> qr, vector<Fr> qo, vector<Fr> qc, vector<int> permutation) {
    Plonk::Circuit ret;

    ret.n = n;
    ret.qm = qm;
    ret.ql = ql;
    ret.qr = qr;
    ret.qo = qo;
    ret.qc = qc;

    Fr w = findPrimitiveRoot(n);
    // if (!isQuadraticResidue(w)) cout << "Not quadratic residue!\n";

    Fr curr = w;
    ret.h = {1};
    for (size_t i = 0; i < n; i++) {
        ret.h.push_back(curr);
        Fr::mul(curr, curr, w);
    }

    Fr k1 = quadraticNonResidue();
    ret.k1 = k1;
    k1 *= w;
    for (size_t i = 0; i < n; i++) {
        ret.h.push_back(k1);
        Fr::mul(k1, k1, w);
    }

    Fr k2 = quadraticNonResidue();
    while (inCoset(ret.h, k2)) { k2 = quadraticNonResidue(); }
    ret.k2 = k2;
    Fr::mul(k2, k2, w);
    for (size_t i = 0; i < n; i++) {
        ret.h.push_back(k2);
        Fr::mul(k2, k2, w);
    }

    for (size_t i = 1; i <= 3 * n; i++) {
        ret.permutation.push_back(ret.h[permutation[i]]);
    }

    return ret;
}

void debugSetup(Plonk::Preprocess ret, Plonk::Circuit circuit) {
    cout << "Checking qm, ql, qr, qo, qc values...\n";
    for (size_t i = 1; i <= circuit.n; i++) {
        cout << "At " << i << "\n";
        cout << evaluatePoly(ret.qm, circuit.h[i]) << endl;
        cout << evaluatePoly(ret.ql, circuit.h[i]) << endl;
        cout << evaluatePoly(ret.qr, circuit.h[i]) << endl;
        cout << evaluatePoly(ret.qo, circuit.h[i]) << endl;
        cout << evaluatePoly(ret.qc, circuit.h[i]) << endl;

        Fr o1 = evaluatePoly(ret.so1, circuit.h[i]);
        Fr o2 = evaluatePoly(ret.so2, circuit.h[i]);
        Fr o3 = evaluatePoly(ret.so3, circuit.h[i]);
        if (o1 == circuit.permutation[i] && 
            o2 == circuit.permutation[circuit.n + i] && 
            o3 == circuit.permutation[circuit.n * 2 + i]) cout << "SO verified\n";
        else cout << "SO failed\n";
    }
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
    G1::mul(one, one, 1);

    Fr curr = 1;
    for (size_t i = 0; i <= circuit.n + 5; i++) { 
        G1 ins;
        G1::mul(ins, one, curr);  
        
        ret.srs.push_back(ins);
        
        Fr::mul(curr, curr, x);
    }

    G2 two;
    mapToG2(two, 1);

    ret.srs2.push_back(two);
    G2::mul(two, two, x);
    ret.srs2.push_back(two);
    
    for (size_t i = 1; i <= circuit.n; i++) {
        for (size_t j = 0; j < ret.lagrange[i].size(); j++) {
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

    // debugSetup(ret, circuit);

    return ret;
}

vector<Fr> abc(size_t n, Fr b1, Fr b2, Fr* w, vector<vector<Fr>> lagrange) {
    vector<Fr> ret(n + 2, 0);
    
    // b1 * X^(n+1) - b1 * X + b2 * X^n - b2
    ret[n+1] = b1;
    ret[1] = -b1;
    
    ret[n] = b2;
    ret[0] = -b2;
    
    for (size_t i = 1; i <= n; i++) {  
        for (size_t j = 0; j < lagrange[i].size(); j++) { 
            ret[j] += w[i] * lagrange[i][j];
        }
    }

    return ret;
}

void debugRoundOne(Plonk::Transcript ret, Plonk::Preprocess prep, vector<Fr> w) {
    cout << "Checking values of abc...\n";
    for (size_t i = 1; i <= prep.circuit.n; i++) {
        Fr a = evaluatePoly(ret.a, prep.circuit.h[i]);
        Fr b = evaluatePoly(ret.b, prep.circuit.h[i]);
        Fr c = evaluatePoly(ret.c, prep.circuit.h[i]);

        cout << "At " << i << "\n";
        if (a == w[i]) cout << "a passed " << a << "\n";
        else cout << "a failed\n";
        if (b == w[prep.circuit.n + i]) cout << "b passed " << b << "\n";
        else cout << "b failed\n";
        if (c == w[prep.circuit.n * 2 + i]) cout << "c passed " << c << "\n";
        else cout << "c failed\n";
    } 
}

void roundOne(Plonk::Transcript &ret, Plonk::Preprocess preprocess, vector<Fr> w, vector<Fr> b) {
    // Alternatively, pass original vector using data() and add offset.
    size_t n = preprocess.circuit.n;
    ret.a = abc(n, b[1], b[2], &w[0], preprocess.lagrange);
    ret.b = abc(n, b[3], b[4], &w[n], preprocess.lagrange);
    ret.c = abc(n, b[5], b[6], &w[2 * n], preprocess.lagrange);

    ret.ax = evaluatePolynomialG1(ret.a, preprocess.srs);
    ret.bx = evaluatePolynomialG1(ret.b, preprocess.srs);
    ret.cx = evaluatePolynomialG1(ret.c, preprocess.srs);

    // debugRoundOne(ret, preprocess, w);
}

void debugRoundTwo(Plonk::Transcript ret, Plonk::Preprocess prep, vector<Fr> w, vector<Fr> b, Plonk::Challenge challs) {
    cout << "Checking values of z...\n";
    Fr ratio = 1;
    for (size_t i = 1; i <= prep.circuit.n; i++) {
        cout << "At " << i << " ";
        Fr zx;
        Fr::pow(zx, prep.circuit.h[i], prep.circuit.n);
        zx -= 1;

        zx *= (b[7] * prep.circuit.h[2 * i] + b[8] * prep.circuit.h[i] + b[9]);

        zx += ratio;

        size_t n = prep.circuit.n;
        ratio *= w[i] + challs.beta * prep.circuit.h[i] + challs.gamma;
        ratio *= w[n + i] + challs.beta * prep.circuit.h[n + i] + challs.gamma;
        ratio *= w[2 * n + i] + challs.beta * prep.circuit.h[2 * n + i] + challs.gamma;

        ratio /= w[i] + challs.beta * prep.circuit.permutation[i] + challs.gamma;
        ratio /= w[n + i] + challs.beta * prep.circuit.permutation[n + i] + challs.gamma;
        ratio /= w[2 * n + i] + challs.beta * prep.circuit.permutation[2 * n + i] + challs.gamma;
        

        if (zx == evaluatePoly(ret.z, prep.circuit.h[i])) cout << "z passed\n";
        else cout << "z failed\n";
    }
}

void roundTwo(Plonk::Transcript &ret, Plonk::Preprocess preprocess, vector<Fr> w, vector<Fr> b, Plonk::Challenge challs) { 
    size_t n = preprocess.circuit.n;
    ret.z.resize(n + 3, 0);
    ret.z[n+2] = b[7];
    ret.z[n+1] = b[8];
    ret.z[n] = b[9];

    ret.z[2] = -b[7];
    ret.z[1] = -b[8];
    ret.z[0] = -b[9];

    for (size_t i = 0; i < preprocess.lagrange[1].size(); i++) {
        ret.z[i] += preprocess.lagrange[1][i];
    }

    Fr ratio = 1;
    for (size_t i = 1; i < n; i++) {
        for (size_t j = i; j < 3 * n; j += n) {
            Fr temp = challs.gamma, _;
            Fr::add(temp, temp, w[j]);

            Fr::mul(_, challs.beta, preprocess.circuit.h[j]);
            Fr::add(temp, temp, _);

            Fr::mul(ratio, ratio, temp);

            temp = challs.gamma;
            Fr::add(temp, temp, w[j]);

            Fr::mul(_, challs.beta, preprocess.circuit.permutation[j]);
            Fr::add(temp, temp, _);

            Fr::div(ratio, ratio, temp);
        }

        for (size_t j = 0; j < preprocess.lagrange[i+1].size(); j++) {
            ret.z[j] += preprocess.lagrange[i+1][j] * ratio;
        }
    }

    ret.zx = evaluatePolynomialG1(ret.z, preprocess.srs);

    // debugRoundTwo(ret, preprocess, w, b, challs);
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

void split(Plonk::Transcript &ret, size_t n, vector<Fr> b) {
    ret.lo.resize(n + 1, 0);
    ret.mid.resize(n + 1, 0);
    ret.hi.resize(n + 6, 0);

    for (size_t i = 0; i < n; i++) {
        ret.lo[i] = ret.t[i];
        ret.mid[i] = ret.t[i + n];
    }

    for (size_t i = 2 * n; i < ret.t.size(); i++) { ret.hi[i - 2 * n] = ret.t[i]; }

    ret.lo[n] += b[10];
    ret.mid[0] -= b[10];
    ret.mid[n] += b[11];
    ret.hi[0] -= b[11];
}

void debugSplit(Plonk::Transcript ret, Plonk::Preprocess prep) {
    size_t n = prep.circuit.n;

    cout << "Checking split...\n";

    for (size_t i = 1; i <= n; i++) {
        cout << "At " << i << " ";
        Fr wn, w2n;
        Fr::pow(wn, prep.circuit.h[i], n);
        Fr::pow(w2n, prep.circuit.h[i], 2 * n);
        if (evaluatePoly(ret.lo, prep.circuit.h[i]) + 
            wn * evaluatePoly(ret.mid, prep.circuit.h[i]) + 
            w2n * evaluatePoly(ret.hi, prep.circuit.h[i]) == evaluatePoly(ret.t, prep.circuit.h[i])) cout << "Split passed\n";
        else cout << "Split failed\n";
    }

    cout << "Checking split polynomial...\n";
    // cout << ret.t.size() << endl;
    // cout << ret.lo.size() << endl;
    // cout << ret.mid.size() << endl;
    // cout << ret.hi.size() << endl;
    bool match = true;
    for (size_t i = 0; i < n; i++) {
        if (ret.lo[i] != ret.t[i]) match = false;
    }
    if (ret.t[n] != ret.lo[n] + ret.mid[0]) match = false;
    for (size_t i = n + 1; i < 2 * n; i++) {
        if (ret.mid[i - n] != ret.t[i]) match = false;
    }
    if (ret.t[2 * n] != ret.mid[n] + ret.hi[0]) match = false;
    for (size_t i = 2 * n + 1; i < ret.t.size(); i++) {
        if (ret.hi[i - n - n] != ret.t[i]) match = false;
    }

    if (match) cout << "Polynomials verified.\n";
    else cout << "Polynomials failed.\n";
}

void debugRoundThree(Plonk::Transcript ret, Plonk::Preprocess prep) {
    size_t n = prep.circuit.n;

    cout << "Checking t values...\n";

    for (size_t i = 1; i <= n; i++) {
        cout << "At " << i << " ";
        cout << evaluatePoly(ret.t, prep.circuit.h[i]) << endl;
    }
}

// Have yet to include public inputs
void roundThree(Plonk::Transcript &ret, Plonk::Preprocess preprocess, size_t l, vector<Fr> w, vector<Fr> b, Plonk::Challenge challs) {
    size_t n = preprocess.circuit.n;
    Fr omega = preprocess.circuit.h[1];
    
    ret.t.resize(3 * n + 5, 0);

    // Constraint
    vector<Fr> temp, _;

    temp = polynomial_multiply(ret.a, ret.b);
    temp = polynomial_multiply(temp, preprocess.qm);

    _ = polynomial_multiply(ret.a, preprocess.ql);
    temp = addPolynomials(temp, _);
    
    _.clear();
    _ = polynomial_multiply(ret.b, preprocess.qr);
    temp = addPolynomials(temp, _);
    
    _.clear();
    _ = polynomial_multiply(ret.c, preprocess.qo);
    temp = addPolynomials(temp, _);

    temp = addPolynomials(temp, preprocess.qc);

    // vector<Fr> pi(n + 1, 0);

    // for (size_t i = 1; i <= l; i++) {
    //     for (size_t j = 0; j < preprocess.lagrange[i].size(); j++) {
    //         pi[j] -= preprocess.lagrange[i][j] * w[i];
    //     }
    // }

    // ret.pi = pi;

    // temp = addPolynomials(temp, pi);

    ret.t = addPolynomials(ret.t, temp);

    // Checks for Round 2 polynomials
    // First part
    temp.clear();
    temp.assign(ret.a.begin(), ret.a.end());
    temp[1] += challs.beta; temp[0] += challs.gamma;

    _.clear();
    _.assign(ret.b.begin(), ret.b.end());
    _[1] += challs.beta * preprocess.circuit.k1; _[0] += challs.gamma;
    temp = polynomial_multiply(temp, _);

    _.clear();
    _.assign(ret.c.begin(), ret.c.end());
    _[1] += challs.beta * preprocess.circuit.k2; _[0] += challs.gamma;
    temp = polynomial_multiply(temp, _);

    temp = polynomial_multiply(temp, ret.z);

    for (size_t i = 0; i < temp.size(); i++) { temp[i] *= challs.alpha; }

    ret.t = addPolynomials(ret.t, temp);

    // Second part
    vector<Fr> zw;
    zw.assign(ret.z.begin(), ret.z.end());

    vector<Fr> so1, so2, so3;
    so1.assign(preprocess.so1.begin(), preprocess.so1.end());
    so2.assign(preprocess.so2.begin(), preprocess.so2.end());
    so3.assign(preprocess.so3.begin(), preprocess.so3.end());

    for (size_t i = 0; i < so1.size(); i++) { so1[i] *= challs.beta; }
    for (size_t i = 0; i < so2.size(); i++) { so2[i] *= challs.beta; }
    for (size_t i = 0; i < so3.size(); i++) { so3[i] *= challs.beta; }

    // Calculate z(Xw)
    for (size_t i = 1; i < zw.size(); i++) { 
        Fr hp;
        Fr::pow(hp, preprocess.circuit.h[1], i);
        zw[i] *= hp; 
    }

    temp.clear();
    temp.assign(ret.a.begin(), ret.a.end());
    for (size_t i = 0; i < so1.size(); i++) { temp[i] += so1[i]; }
    temp[0] += challs.gamma;

    _.clear();
    _.assign(ret.b.begin(), ret.b.end());
    for (size_t i = 0; i < so2.size(); i++) { _[i] += so2[i]; }
    _[0] += challs.gamma;
    temp = polynomial_multiply(temp, _);

    _.clear();
    _.assign(ret.c.begin(), ret.c.end());
    for (size_t i = 0; i < so3.size(); i++) { _[i] += so3[i]; }
    _[0] += challs.gamma;
    temp = polynomial_multiply(temp, _);

    temp = polynomial_multiply(temp, zw);

    for (size_t i = 0; i < temp.size(); i++) { temp[i] *= -challs.alpha; }

    ret.t = addPolynomials(ret.t, temp);

    // Induction base case check
    ret.z[0] -= 1;
    temp.clear();
    temp = polynomial_multiply(ret.z, preprocess.lagrange[1]);
    for (size_t i = 0; i < temp.size(); i++) { temp[i] *= challs.alpha * challs.alpha; }
    ret.z[0] += 1;

    ret.t = addPolynomials(ret.t, temp);
    
    // debugRoundThree(ret, preprocess);
    
    vector<Fr> temp_t;
    temp_t = polynomialDivision(ret.t, n); // Divide by ZH

    vector<Fr> zh(n+1, 0); zh[n] = 1; zh[0] = -1;

    temp.clear();
    temp = polynomial_multiply(temp_t, zh);
    bool correct = true;
    for (size_t i = 0; i < temp.size(); i++) {
        if (temp[i] != ret.t[i]) {
            correct = false;
            break;
        }
    }
    if (correct) cout << "Correct polynomial division!\n";
    else cout << "Wrong polynomial division!\n";

    ret.t.clear();
    ret.t.assign(temp_t.begin(), temp_t.end());

    split(ret, n, b);

    ret.lox = evaluatePolynomialG1(ret.lo, preprocess.srs);
    ret.midx = evaluatePolynomialG1(ret.mid, preprocess.srs);
    ret.hix = evaluatePolynomialG1(ret.hi, preprocess.srs);

    // debugSplit(ret, preprocess);
}

void roundFour(Plonk::Transcript &ret, Plonk::Preprocess preprocess, vector<Fr> vs) {
    ret.av = evaluatePoly(ret.a, vs[1]);
    ret.bv = evaluatePoly(ret.b, vs[1]);
    ret.cv = evaluatePoly(ret.c, vs[1]);

    ret.so1v = evaluatePoly(preprocess.so1, vs[1]);
    ret.so2v = evaluatePoly(preprocess.so2, vs[1]);

    Fr curr = 1;
    for (size_t i = 0; i < vs.size(); i++) {
        vs[i] *= curr;
        curr *= preprocess.circuit.h[1];
    }

    ret.zwv = evaluatePolynomial(ret.z, vs);
}

void debugRoundFive(Plonk::Transcript transcript, Plonk::Preprocess preprocess, Fr v, bool first) {
    cout << "Checking r values...\n";
    
    if (first) {
        Fr t = evaluatePoly(transcript.t, v);
        Fr r = evaluatePoly(transcript.r, v);
        Fr vn;
        Fr::pow(vn, v, preprocess.circuit.n);
        vn -= 1;

        if (r == t * vn) cout << "First half of r verified.\n";
        else cout << "First half of r failed.\n";
        return;
    }

    for (size_t i = 1; i <= preprocess.circuit.n; i++) {
        cout << evaluatePoly(transcript.r, preprocess.circuit.h[i]) << endl;
    }
    cout << "V: " << evaluatePoly(transcript.r, v) << endl;
}

void debugRoundFive(Plonk::Transcript transcript, Plonk::Preprocess preprocess, Fr v, size_t f) {
    if (f == 1) {
        if (evaluatePoly(transcript.w, v) == 0) cout << "w verified.\n";
        else cout << "w failed.\n";
        return;
    }

    if (evaluatePoly(transcript.ww, v * preprocess.circuit.h[1]) == 0) cout << "ww verified.\n";
    else cout << "ww failed.\n";
}

void roundFive(Plonk::Transcript &ret, Plonk::Preprocess preprocess, Plonk::Challenge challs, vector<Fr> vs) {
    // Evaluating r(X)
    size_t n = preprocess.circuit.n;
    ret.r.resize(3 * n + 5, 0); 

    for (size_t i = 0; i < preprocess.qm.size(); i++) { ret.r[i] += ret.av * ret.bv * preprocess.qm[i]; }
    for (size_t i = 0; i < preprocess.ql.size(); i++) { ret.r[i] += ret.av * preprocess.ql[i]; }
    for (size_t i = 0; i < preprocess.qr.size(); i++) { ret.r[i] += ret.bv * preprocess.qr[i]; }
    for (size_t i = 0; i < preprocess.qo.size(); i++) { ret.r[i] += ret.cv * preprocess.qo[i]; }
    for (size_t i = 0; i < preprocess.qc.size(); i++) { ret.r[i] += preprocess.qc[i]; }

    // Fr piv = evaluatePolynomial(ret.pi, vs);
    // ret.r[0] += piv;

    for (size_t i = 0; i < ret.z.size(); i++) { 
        ret.r[i] += ret.z[i] * challs.alpha * (ret.av + challs.beta * vs[1] + challs.gamma) * 
            (ret.bv + challs.beta * vs[1] * preprocess.circuit.k1 + challs.gamma) * 
            (ret.cv + challs.beta * vs[1] * preprocess.circuit.k2 + challs.gamma); 
    }

    vector<Fr> so3;
    so3.assign(preprocess.so3.begin(), preprocess.so3.end());
    for (size_t i = 0; i < so3.size(); i++) { so3[i] *= challs.beta; }
    so3[0] += ret.cv + challs.gamma; 

    for (size_t i = 0; i < so3.size(); i++) {
        ret.r[i] -= so3[i] * challs.alpha * ret.zwv * (ret.av + challs.beta * ret.so1v + challs.gamma) * 
            (ret.bv + challs.beta * ret.so2v + challs.gamma);
    }

    Fr l1v = evaluatePoly(preprocess.lagrange[1], vs[1]);
    for (size_t i = 0; i < ret.z.size(); i++) { ret.r[i] += challs.alpha * challs.alpha * ret.z[i] * l1v; }
    ret.r[0] -= challs.alpha * challs.alpha * l1v;

    // debugRoundFive(ret, preprocess, vs[1], true);

    Fr zhv = vs[n] - 1;
    for (size_t i = 0; i < ret.lo.size(); i++) { ret.r[i] -= ret.lo[i] * zhv; }
    for (size_t i = 0; i < ret.mid.size(); i++) { ret.r[i] -= ret.mid[i] * zhv * vs[n]; }
    for (size_t i = 0; i < ret.hi.size(); i++) { ret.r[i] -= ret.hi[i] * zhv * vs[n] * vs[n]; }

    // debugRoundFive(ret, preprocess, vs[1], false);

    // Evaluating W(X)
    vector<Fr> us(n + 5, 1);
    Fr curr = challs.u;
    for (size_t i = 1; i <= n + 4; i++) {
        us[i] = curr;
        curr *= challs.u;
    }

    ret.w.assign(ret.r.begin(), ret.r.end());

    for (size_t i = 0; i < ret.a.size(); i++) { ret.w[i] += ret.a[i] * us[1]; }
    for (size_t i = 0; i < ret.b.size(); i++) { ret.w[i] += ret.b[i] * us[2]; }
    for (size_t i = 0; i < ret.c.size(); i++) { ret.w[i] += ret.c[i] * us[3]; }
    for (size_t i = 0; i < preprocess.so1.size(); i++) { ret.w[i] += preprocess.so1[i] * us[4]; }
    for (size_t i = 0; i < preprocess.so2.size(); i++) { ret.w[i] += preprocess.so2[i] * us[5]; }

    ret.w[0] -= ret.av * us[1] + ret.bv * us[2] + ret.cv * us[3] + ret.so1v * us[4] + ret.so2v * us[5];

    // debugRoundFive(ret, preprocess, vs[1], 1);

    ret.w = divideByLinear(ret.w, vs[1]);

    while (!ret.w.empty() && ret.w.back().isZero()) {
        ret.w.pop_back();
    }

    ret.wx = evaluatePolynomialG1(ret.w, preprocess.srs);

    // Evaluating Ww(X)
    ret.ww.clear();
    ret.ww.assign(ret.z.begin(), ret.z.end());
    ret.ww[0] -= ret.zwv;

    // debugRoundFive(ret, preprocess, vs[1], 2);

    ret.ww = divideByLinear(ret.ww, vs[1] * preprocess.circuit.h[1]);
    ret.wwx = evaluatePolynomialG1(ret.ww, preprocess.srs);
}

// Can be optimized by not converting back-and-forth polynomials from lagrange basis to monomial basis
// Missing KZG commitments
Plonk::Witness Plonk::prove(Plonk::Preprocess preprocess, size_t l, vector<Fr> w) {
    Plonk::Transcript transcript;
    Plonk::Witness ret;
    ret.challs.k1 = preprocess.circuit.k1;
    ret.challs.k2 = preprocess.circuit.k2;

    vector<Fr> b(22); b[0] = 0;
    for (size_t i = 1; i <= 11; i++) {
        b[i].setByCSPRNG();
    }
    
    // ROUND 1
    roundOne(transcript, preprocess, w, b);

    // ROUND 2
    ret.challs.beta.setByCSPRNG(); // Should be H(transcript, 0)
    ret.challs.gamma.setByCSPRNG(); // Should be H(transcript, 1)
    roundTwo(transcript, preprocess, w, b, ret.challs);

    // ROUND 3
    ret.challs.alpha.setByCSPRNG(); // Should be H(transcript)
    roundThree(transcript, preprocess, l, w, b, ret.challs);

    // ROUND 4
    Fr v; 
    v.setByCSPRNG(); // Should be H(transcript)
    ret.challs.v = v;

    vector<Fr> vs(preprocess.circuit.n + 5, 1);
    Fr curr = 1;
    for (size_t i = 0; i < vs.size(); i++) {
        vs[i] = curr;
        curr *= v;
    }
    roundFour(transcript, preprocess, vs);

    // ROUND 5
    ret.challs.u.setByCSPRNG(); // Should be H(transcript)
    roundFive(transcript, preprocess, ret.challs, vs);

    ret.challs.e.setByCSPRNG(); // Should be H(transcript)

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

Plonk::Verifier Plonk::preprocess(Plonk::Preprocess preprocess) {
    Plonk::Verifier ret;

    ret.qmx = evaluatePolynomialG1(preprocess.qm, preprocess.srs);
    ret.qlx = evaluatePolynomialG1(preprocess.ql, preprocess.srs);
    ret.qrx = evaluatePolynomialG1(preprocess.qr, preprocess.srs);
    ret.qox = evaluatePolynomialG1(preprocess.qo, preprocess.srs);
    ret.qcx = evaluatePolynomialG1(preprocess.qc, preprocess.srs);

    ret.so1x = evaluatePolynomialG1(preprocess.so1, preprocess.srs);
    ret.so2x = evaluatePolynomialG1(preprocess.so2, preprocess.srs);
    ret.so3x = evaluatePolynomialG1(preprocess.so3, preprocess.srs);
    
    ret.x = preprocess.srs[1];

    return ret;
}

// Verifier can only access up until index l of w
bool Plonk::verify(size_t n, Plonk::Preprocess prep, Plonk::Witness witness, size_t l, vector<Fr> w) {
    Plonk::Verifier verifier = preprocess(prep);
    Fr omega = prep.circuit.h[1];

    Plonk::Challenge challs = witness.challs;

    // Supposed to include KZG Validation
    if (!witness.ax.isValid() || !witness.bx.isValid() || !witness.cx.isValid() || 
    !witness.zx.isValid() || !witness.lox.isValid() || !witness.midx.isValid() ||
    !witness.hix.isValid() || !witness.wx.isValid() || !witness.wwx.isValid()) return false;

    if (!witness.av.isValid() || !witness.bv.isValid() || !witness.cv.isValid() ||
    !witness.so1v.isValid() || !witness.so2v.isValid() || !witness.zwv.isValid()) return false;

    for (size_t i = 1; i <= l; i++) if (!w[i].isValid()) return false;

    cout << "✓ Data check passed\n";

    // Compute challenges (skipped as it is generated randomly)

    // Compute zero polynomial evaluation
    Fr zhv;
    Fr::pow(zhv, challs.v, n);
    zhv -= 1;

    // Compute lagrange polynomial evaluation
    Fr l1v = zhv * omega / n;
    l1v /= (challs.v - omega);

    // Compute public input polynomial evaluation
    // vector<Fr> pi(n + 1, 0);

    // for (size_t i = 1; i <= l; i++) {
    //     for (size_t j = 0; j < prep.lagrange[i].size(); j++) {
    //         pi[j] -= prep.lagrange[i][j] * w[i];
    //     }
    // }
    // Fr piv = evaluatePoly(pi, challs.v);

    // Constant term of r 
    Fr r0 = - l1v * challs.alpha * challs.alpha - 
        challs.alpha * (witness.av + challs.beta * witness.so1v + challs.gamma) *
        (witness.bv + challs.beta * witness.so2v + challs.gamma) *
        (witness.cv + challs.gamma) * witness.zwv;

    // r0 += piv;
    
    // Batched polynomial commitment D = r - r0 + uz
    G1 D;
    G1 qmx, qlx, qrx, qox;
    G1::mul(qmx, verifier.qmx, witness.av * witness.bv);
    G1::mul(qlx, verifier.qlx, witness.av);
    G1::mul(qrx, verifier.qrx, witness.bv);
    G1::mul(qox, verifier.qox, witness.cv);
    D = qmx + qlx + qrx + qox + verifier.qcx;

    Fr temp = (witness.av + challs.beta * challs.v + challs.gamma) * 
        (witness.bv + challs.beta * challs.k1 * challs.v + challs.gamma) *
        (witness.cv + challs.beta * challs.k2 * challs.v + challs.gamma) * challs.alpha +
        l1v * challs.alpha * challs.alpha + challs.e;
    
    G1 zx;
    G1::mul(zx, witness.zx, temp);
    D += zx;

    temp = (witness.av + challs.beta * witness.so1v + challs.gamma) * 
        (witness.bv + challs.beta * witness.so2v + challs.gamma) *
        challs.alpha * challs.beta * witness.zwv;
    G1 so3x;
    G1::mul(so3x, verifier.so3x, temp);
    D -= so3x;

    Fr vn, v2n;
    Fr::pow(vn, challs.v, n);
    Fr::pow(v2n, challs.v, 2 * n);

    G1 midx, hix;
    G1::mul(midx, witness.midx, vn);
    G1::mul(hix, witness.hix, v2n);

    G1 _ = witness.lox + midx + hix;
    G1::mul(_, _, zhv);
    D -= _;

    G1 minus; G1::mul(minus, witness.zx, challs.e);
    G1 add; G1::mul(add, prep.srs[0], r0);

    vector<Fr> us(6, 1);
    Fr curr = challs.u;
    for (size_t i = 1; i <= 5; i++) {
        us[i] = curr;
        curr *= challs.u;
    }

    G1 ax, bx, cx, so1x, so2x;
    G1::mul(ax, witness.ax, us[1]);
    G1::mul(bx, witness.bx, us[2]);
    G1::mul(cx, witness.cx, us[3]);
    G1::mul(so1x, verifier.so1x, us[4]);
    G1::mul(so2x, verifier.so2x, us[5]);
    G1 F = D + ax + bx + cx + so1x + so2x;

    G1 E;
    temp = us[1] * witness.av + us[2] * witness.bv + us[3] * witness.cv +
        us[4] * witness.so1v + us[5] * witness.so2v - r0 + challs.e * witness.zwv;
    G1::mul(E, prep.srs[0], temp);

    GT left, right;
    G1 wwx;
    G1::mul(wwx, witness.wwx, challs.e);
    pairing(left, witness.wx + wwx, prep.srs2[1]);

    G1 wx;
    G1::mul(wx, witness.wx, challs.v);
    G1::mul(wwx, witness.wwx, challs.v * omega * challs.e);
    pairing(right, wx + wwx + F - E, prep.srs2[0]);

    return left == right;
}
