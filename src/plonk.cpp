#include "../include/plonk.hpp"
#include "../include/ntt.hpp"
#include "../include/kzg.hpp"
#include "../include/debug.hpp"
#include "../include/helper.hpp"
#include <openssl/evp.h>
#include <sstream>
#include <mcl/bn.hpp>
#include <vector>
#include <map>
#include <cmath>

using namespace std;
using namespace mcl;
using namespace bn;

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
    ret.pk = KZGSetup(circuit.n, x);
    
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

void roundOne(Plonk::Transcript &ret, Plonk::Preprocess preprocess, vector<Fr> w, vector<Fr> b) {
    // Alternatively, pass original vector using data() and add offset.
    size_t n = preprocess.circuit.n;
    ret.a = abc(n, b[1], b[2], &w[0], preprocess.lagrange);
    ret.b = abc(n, b[3], b[4], &w[n], preprocess.lagrange);
    ret.c = abc(n, b[5], b[6], &w[2 * n], preprocess.lagrange);

    ret.ax = commit(preprocess.pk, ret.a);
    ret.bx = commit(preprocess.pk, ret.b);
    ret.cx = commit(preprocess.pk, ret.c);
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

    ret.zx = commit(preprocess.pk, ret.z);
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

void roundThree(Plonk::Transcript &ret, Plonk::Preprocess preprocess, size_t l, vector<Fr> w, vector<Fr> b, Plonk::Challenge challs) {
    size_t n = preprocess.circuit.n;
    Fr omega = preprocess.circuit.h[1];
    
    // Determine the target size for NTT operations (must be power of 2)
    size_t ntt_size = 1;
    while (ntt_size < 3 * n + 5) ntt_size *= 2;

    Fr ntt_omega = getRootOfUnity(ntt_size);
    
    // Convert all polynomials to NTT domain once
    vector<Fr> a_ntt = to_ntt_domain(ret.a, ntt_size, ntt_omega);
    vector<Fr> b_ntt = to_ntt_domain(ret.b, ntt_size, ntt_omega);
    vector<Fr> c_ntt = to_ntt_domain(ret.c, ntt_size, ntt_omega);
    vector<Fr> z_ntt = to_ntt_domain(ret.z, ntt_size, ntt_omega);
    
    vector<Fr> qm_ntt = to_ntt_domain(preprocess.qm, ntt_size, ntt_omega);
    vector<Fr> ql_ntt = to_ntt_domain(preprocess.ql, ntt_size, ntt_omega);
    vector<Fr> qr_ntt = to_ntt_domain(preprocess.qr, ntt_size, ntt_omega);
    vector<Fr> qo_ntt = to_ntt_domain(preprocess.qo, ntt_size, ntt_omega);
    vector<Fr> qc_ntt = to_ntt_domain(preprocess.qc, ntt_size, ntt_omega);
    
    vector<Fr> so1_ntt = to_ntt_domain(preprocess.so1, ntt_size, ntt_omega);
    vector<Fr> so2_ntt = to_ntt_domain(preprocess.so2, ntt_size, ntt_omega);
    vector<Fr> so3_ntt = to_ntt_domain(preprocess.so3, ntt_size, ntt_omega);
    
    // Evaluate t
    vector<Fr> t_ntt(ntt_size, Fr(0)), temp, _;

    temp = a_ntt;
    polynomial_multiply_ntt(temp, b_ntt);
    polynomial_multiply_ntt(temp, qm_ntt);
    polynomial_add_ntt(t_ntt, temp);

    temp.clear(); temp = a_ntt;
    polynomial_multiply_ntt(temp, ql_ntt);
    polynomial_add_ntt(t_ntt, temp);

    temp.clear(); temp = b_ntt;
    polynomial_multiply_ntt(temp, qr_ntt);
    polynomial_add_ntt(t_ntt, temp);

    temp.clear(); temp = c_ntt;
    polynomial_multiply_ntt(temp, qo_ntt);
    polynomial_add_ntt(t_ntt, temp);

    polynomial_add_ntt(t_ntt, qc_ntt);

    // Evaluate PI
    vector<Fr> pi(n + 1, 0);
    for (size_t i = 1; i <= l; i++) {
        for (size_t j = 0; j < preprocess.lagrange[i].size(); j++) {
            pi[j] -= preprocess.lagrange[i][j] * w[i];
        }
    }
    ret.pi = pi;
    
    // Convert PI to NTT domain properly
    vector<Fr> pi_ntt = to_ntt_domain(pi, ntt_size, ntt_omega);
    polynomial_add_ntt(t_ntt, pi_ntt);

    // Checks for Round 2 polynomials
    // First part
    temp.clear(); temp = a_ntt;
    _.clear(); _.resize(ntt_size, Fr(0));
    _[1] = challs.beta; _[0] = challs.gamma;
    ntt_transform(_, ntt_omega); 
    polynomial_add_ntt(temp, _);

    _.clear(); _.resize(ntt_size, Fr(0));
    _[1] = challs.beta * preprocess.circuit.k1; _[0] = challs.gamma;
    ntt_transform(_, ntt_omega); 
    polynomial_add_ntt(_, b_ntt);
    polynomial_multiply_ntt(temp, _);

    _.clear(); _.resize(ntt_size, Fr(0));
    _[1] = challs.beta * preprocess.circuit.k2; _[0] = challs.gamma;
    ntt_transform(_, ntt_omega); 
    polynomial_add_ntt(_, c_ntt);
    polynomial_multiply_ntt(temp, _);

    polynomial_multiply_ntt(temp, z_ntt);
    polynomial_scale_ntt(temp, challs.alpha);
    polynomial_add_ntt(t_ntt, temp);

    // Second part
    vector<Fr> gamma_ntt(ntt_size, Fr(0));
    gamma_ntt[0] = challs.gamma; 
    ntt_transform(gamma_ntt, ntt_omega); 

    temp.clear(); temp = so1_ntt;
    polynomial_scale_ntt(temp, challs.beta);
    polynomial_add_ntt(temp, gamma_ntt);
    polynomial_add_ntt(temp, a_ntt);

    _.clear(); _ = so2_ntt;
    polynomial_scale_ntt(_, challs.beta);
    polynomial_add_ntt(_, gamma_ntt);
    polynomial_add_ntt(_, b_ntt);
    polynomial_multiply_ntt(temp, _);

    _.clear(); _ = so3_ntt;
    polynomial_scale_ntt(_, challs.beta);
    polynomial_add_ntt(_, gamma_ntt);
    polynomial_add_ntt(_, c_ntt);
    polynomial_multiply_ntt(temp, _);
    
    vector<Fr> zw_ntt = ret.z;
    for (size_t i = 1; i < ret.z.size(); i++) {
        Fr temp;
        Fr::pow(temp, omega, i);
        zw_ntt[i] *= temp; 
    }
    zw_ntt.resize(ntt_size, Fr(0));
    zw_ntt = to_ntt_domain(zw_ntt, ntt_size, ntt_omega);

    polynomial_multiply_ntt(temp, zw_ntt);
    polynomial_scale_ntt(temp, -challs.alpha);
    polynomial_add_ntt(t_ntt, temp);

    // Induction base case
    temp.clear(); temp = z_ntt;
    _.clear(); _.resize(ntt_size, Fr(0));
    _[0] = -1;
    ntt_transform(_, ntt_omega); 
    polynomial_add_ntt(temp, _);
    
    vector<Fr> l1_ntt = to_ntt_domain(preprocess.lagrange[1], ntt_size, ntt_omega);
    polynomial_multiply_ntt(temp, l1_ntt);
    polynomial_scale_ntt(temp, challs.alpha * challs.alpha);
    
    polynomial_add_ntt(t_ntt, temp);
    
    // Convert t back to coefficient form
    ntt_inverse(t_ntt, ntt_omega);
    ret.t = t_ntt;

    // Divide by ZH
    ret.t = polynomialDivision(ret.t, n);
    
    // Split t polynomial
    split(ret, n, b);
    
    // Commit to the pieces
    ret.lox = commit(preprocess.pk, ret.lo);
    ret.midx = commit(preprocess.pk, ret.mid);
    ret.hix = commit(preprocess.pk, ret.hi);
}

void roundFour(Plonk::Transcript &ret, Plonk::Preprocess preprocess, vector<Fr> vs) {
    ret.av = createWitness(preprocess.pk, ret.a, vs[1]);
    ret.bv = createWitness(preprocess.pk, ret.b, vs[1]);
    ret.cv = createWitness(preprocess.pk, ret.c, vs[1]);

    ret.so1v = createWitness(preprocess.pk, preprocess.so1, vs[1]);
    ret.so2v = createWitness(preprocess.pk, preprocess.so2, vs[1]);

    ret.zwv = createWitness(preprocess.pk, ret.z, vs[1] * preprocess.circuit.h[1]);
}

void roundFive(Plonk::Transcript &ret, Plonk::Preprocess preprocess, Plonk::Challenge challs, vector<Fr> vs) {
    // Evaluating r(X)
    size_t n = preprocess.circuit.n;
    ret.r.resize(3 * n + 5, 0); 

    for (size_t i = 0; i < preprocess.qm.size(); i++) { ret.r[i] += ret.av.qi * ret.bv.qi * preprocess.qm[i]; }
    for (size_t i = 0; i < preprocess.ql.size(); i++) { ret.r[i] += ret.av.qi * preprocess.ql[i]; }
    for (size_t i = 0; i < preprocess.qr.size(); i++) { ret.r[i] += ret.bv.qi * preprocess.qr[i]; }
    for (size_t i = 0; i < preprocess.qo.size(); i++) { ret.r[i] += ret.cv.qi * preprocess.qo[i]; }
    for (size_t i = 0; i < preprocess.qc.size(); i++) { ret.r[i] += preprocess.qc[i]; }

    Fr piv = evaluatePolynomial(ret.pi, vs);
    ret.r[0] += piv;

    for (size_t i = 0; i < ret.z.size(); i++) { 
        ret.r[i] += ret.z[i] * challs.alpha * (ret.av.qi + challs.beta * vs[1] + challs.gamma) * 
            (ret.bv.qi + challs.beta * vs[1] * preprocess.circuit.k1 + challs.gamma) * 
            (ret.cv.qi + challs.beta * vs[1] * preprocess.circuit.k2 + challs.gamma); 
    }

    vector<Fr> so3;
    so3.assign(preprocess.so3.begin(), preprocess.so3.end());
    for (size_t i = 0; i < so3.size(); i++) { so3[i] *= challs.beta; }
    so3[0] += ret.cv.qi + challs.gamma; 

    for (size_t i = 0; i < so3.size(); i++) {
        ret.r[i] -= so3[i] * challs.alpha * ret.zwv.qi * (ret.av.qi + challs.beta * ret.so1v.qi + challs.gamma) * 
            (ret.bv.qi + challs.beta * ret.so2v.qi + challs.gamma);
    }

    Fr l1v = evaluatePoly(preprocess.lagrange[1], vs[1]);
    for (size_t i = 0; i < ret.z.size(); i++) { ret.r[i] += challs.alpha * challs.alpha * ret.z[i] * l1v; }
    ret.r[0] -= challs.alpha * challs.alpha * l1v;

    Fr zhv = vs[n] - 1;
    for (size_t i = 0; i < ret.lo.size(); i++) { ret.r[i] -= ret.lo[i] * zhv; }
    for (size_t i = 0; i < ret.mid.size(); i++) { ret.r[i] -= ret.mid[i] * zhv * vs[n]; }
    for (size_t i = 0; i < ret.hi.size(); i++) { ret.r[i] -= ret.hi[i] * zhv * vs[n] * vs[n]; }

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

    ret.w[0] -= ret.av.qi * us[1] + ret.bv.qi * us[2] + ret.cv.qi * us[3] + ret.so1v.qi * us[4] + ret.so2v.qi * us[5];

    ret.w = divideByLinear(ret.w, vs[1]);

    while (!ret.w.empty() && ret.w.back().isZero()) {
        ret.w.pop_back();
    }

    ret.wx = commit(preprocess.pk, ret.w);

    // Evaluating Ww(X)
    ret.ww.clear();
    ret.ww.assign(ret.z.begin(), ret.z.end());
    ret.ww[0] -= ret.zwv.qi;

    ret.ww = divideByLinear(ret.ww, vs[1] * preprocess.circuit.h[1]);
    ret.wwx = commit(preprocess.pk, ret.ww);
}

Plonk::Witness Plonk::prove(Plonk::Preprocess preprocess, size_t l, vector<Fr> w) {
    Plonk::Transcript transcript;
    Plonk::Witness ret;
    Plonk::Challenge challs;

    // Transcript to be hashed to generate challenges.
    vector<vector<uint8_t>> hashed;

    vector<Fr> b(22); b[0] = 0;
    for (size_t i = 1; i <= 11; i++) {
        b[i].setByCSPRNG();
    }
    
    // ROUND 1
    roundOne(transcript, preprocess, w, b);
    hashed.push_back(toBytes(transcript.ax.c));
    hashed.push_back(toBytes(transcript.bx.c));
    hashed.push_back(toBytes(transcript.cx.c));
    hashed.push_back(toBytes(Fr(0)));

    // ROUND 2
    challs.beta = hashToFr(hashed);
    hashed[hashed.size() - 1] = toBytes(Fr(1));
    challs.gamma = hashToFr(hashed);
    roundTwo(transcript, preprocess, w, b, challs);
    hashed[hashed.size() - 1] = toBytes(transcript.zx.c);

    // ROUND 3
    challs.alpha = hashToFr(hashed);
    roundThree(transcript, preprocess, l, w, b, challs);
    hashed.push_back(toBytes(transcript.lox.c));
    hashed.push_back(toBytes(transcript.midx.c));
    hashed.push_back(toBytes(transcript.hix.c));

    // ROUND 4
    challs.v = hashToFr(hashed);

    vector<Fr> vs(preprocess.circuit.n + 5, 1);
    Fr curr = 1;
    for (size_t i = 0; i < vs.size(); i++) {
        vs[i] = curr;
        curr *= challs.v;
    }
    roundFour(transcript, preprocess, vs);

    hashed.push_back(toBytes(transcript.av.qi));
    hashed.push_back(toBytes(transcript.bv.qi));
    hashed.push_back(toBytes(transcript.cv.qi));
    hashed.push_back(toBytes(transcript.so1v.qi));
    hashed.push_back(toBytes(transcript.so2v.qi));
    hashed.push_back(toBytes(transcript.zwv.qi));

    // ROUND 5
    challs.u = hashToFr(hashed);
    roundFive(transcript, preprocess, challs, vs);
    hashed.push_back(toBytes(transcript.wx.c));
    hashed.push_back(toBytes(transcript.wwx.c));

    challs.e = hashToFr(hashed);

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

    ret.qmx = commit(preprocess.pk, preprocess.qm);
    ret.qlx = commit(preprocess.pk, preprocess.ql);
    ret.qrx = commit(preprocess.pk, preprocess.qr);
    ret.qox = commit(preprocess.pk, preprocess.qo);
    ret.qcx = commit(preprocess.pk, preprocess.qc);

    ret.so1x = commit(preprocess.pk, preprocess.so1);
    ret.so2x = commit(preprocess.pk, preprocess.so2);
    ret.so3x = commit(preprocess.pk, preprocess.so3);

    return ret;
}

Plonk::Challenge evaluateChalls(Plonk::Witness witness) {
    Plonk::Challenge challs;
    vector<vector<uint8_t>> hashed;

    hashed.push_back(toBytes(witness.ax.c));
    hashed.push_back(toBytes(witness.bx.c));
    hashed.push_back(toBytes(witness.cx.c));
    hashed.push_back(toBytes(Fr(0)));

    challs.beta = hashToFr(hashed);
    hashed[hashed.size() - 1] = toBytes(Fr(1));
    challs.gamma = hashToFr(hashed);
    hashed[hashed.size() - 1] = toBytes(witness.zx.c);

    challs.alpha = hashToFr(hashed);
    hashed.push_back(toBytes(witness.lox.c));
    hashed.push_back(toBytes(witness.midx.c));
    hashed.push_back(toBytes(witness.hix.c));

    challs.v = hashToFr(hashed);
    hashed.push_back(toBytes(witness.av.qi));
    hashed.push_back(toBytes(witness.bv.qi));
    hashed.push_back(toBytes(witness.cv.qi));
    hashed.push_back(toBytes(witness.so1v.qi));
    hashed.push_back(toBytes(witness.so2v.qi));
    hashed.push_back(toBytes(witness.zwv.qi));

    challs.u = hashToFr(hashed);
    hashed.push_back(toBytes(witness.wx.c));
    hashed.push_back(toBytes(witness.wwx.c));

    challs.e = hashToFr(hashed);

    return challs;
}

// Verifier can only access up until index l of w
bool Plonk::verify(size_t n, Plonk::Preprocess prep, Plonk::Witness witness, size_t l, vector<Fr> w) {
    Plonk::Verifier verifier = preprocess(prep);
    Fr omega = prep.circuit.h[1];

    // Compute challenge based on witness / transcript
    Plonk::Challenge challs = evaluateChalls(witness);

    // Basic validity checks
    if (!witness.ax.c.isValid() || !witness.bx.c.isValid() || !witness.cx.c.isValid() || 
    !witness.zx.c.isValid() || !witness.lox.c.isValid() || !witness.midx.c.isValid() ||
    !witness.hix.c.isValid() || !witness.wx.c.isValid() || !witness.wwx.c.isValid()) return false;

    for (size_t i = 1; i <= l; i++) if (!w[i].isValid()) return false;

    // Prepare batch verification for KZG proofs
    vector<KZG::BatchItem> items;
    
    // Add all KZG evaluations to batch
    items = addItems(prep, witness, verifier, challs);

    // Generate random coefficients for batch verification
    vector<Fr> random_coeffs = generateRandomCoeffs(items.size());

    // Perform batch verification of all KZG proofs
    if (!verifyBatch(prep.pk, items, random_coeffs)) {
        cout << "✗ KZG batch verification failed\n";
        return false;
    }

    cout << "✓ KZG batch verification passed\n";

    // Compute zero polynomial evaluation
    Fr zhv;
    Fr::pow(zhv, challs.v, n);
    zhv -= 1;

    // Compute lagrange polynomial evaluation
    Fr l1v = zhv * omega / n;
    l1v /= (challs.v - omega);

    // Compute public input polynomial evaluation
    vector<Fr> pi(n + 1, 0);

    for (size_t i = 1; i <= l; i++) {
        for (size_t j = 0; j < prep.lagrange[i].size(); j++) {
            pi[j] -= prep.lagrange[i][j] * w[i];
        }
    }
    Fr piv = evaluatePoly(pi, challs.v);

    // Constant term of r 
    Fr r0 = - l1v * challs.alpha * challs.alpha - 
        challs.alpha * (witness.av.qi + challs.beta * witness.so1v.qi + challs.gamma) *
        (witness.bv.qi + challs.beta * witness.so2v.qi + challs.gamma) *
        (witness.cv.qi + challs.gamma) * witness.zwv.qi;

    r0 += piv;
    
    // Batched polynomial commitment D = r - r0 + uz
    G1 D;
    G1 qmx, qlx, qrx, qox;
    G1::mul(qmx, verifier.qmx.c, witness.av.qi * witness.bv.qi);
    G1::mul(qlx, verifier.qlx.c, witness.av.qi);
    G1::mul(qrx, verifier.qrx.c, witness.bv.qi);
    G1::mul(qox, verifier.qox.c, witness.cv.qi);
    D = qmx + qlx + qrx + qox + verifier.qcx.c;

    Fr temp = (witness.av.qi + challs.beta * challs.v + challs.gamma) * 
        (witness.bv.qi + challs.beta * prep.circuit.k1 * challs.v + challs.gamma) *
        (witness.cv.qi + challs.beta * prep.circuit.k2 * challs.v + challs.gamma) * challs.alpha +
        l1v * challs.alpha * challs.alpha + challs.e;
    
    G1 zx;
    G1::mul(zx, witness.zx.c, temp);
    D += zx;

    temp = (witness.av.qi + challs.beta * witness.so1v.qi + challs.gamma) * 
        (witness.bv.qi + challs.beta * witness.so2v.qi + challs.gamma) *
        challs.alpha * challs.beta * witness.zwv.qi;

    G1 so3x;
    G1::mul(so3x, verifier.so3x.c, temp);
    D -= so3x;

    Fr vn, v2n;
    Fr::pow(vn, challs.v, n);
    Fr::pow(v2n, challs.v, 2 * n);

    G1 midx, hix;
    G1::mul(midx, witness.midx.c, vn);
    G1::mul(hix, witness.hix.c, v2n);

    G1 _ = witness.lox.c + midx + hix;
    G1::mul(_, _, zhv);
    D -= _;

    vector<Fr> us(6, 1);
    Fr curr = challs.u;
    for (size_t i = 1; i <= 5; i++) {
        us[i] = curr;
        curr *= challs.u;
    }

    G1 ax, bx, cx, so1x, so2x;
    G1::mul(ax, witness.ax.c, us[1]);
    G1::mul(bx, witness.bx.c, us[2]);
    G1::mul(cx, witness.cx.c, us[3]);
    G1::mul(so1x, verifier.so1x.c, us[4]);
    G1::mul(so2x, verifier.so2x.c, us[5]);
    G1 F = D + ax + bx + cx + so1x + so2x;

    G1 E;
    temp = us[1] * witness.av.qi + us[2] * witness.bv.qi + us[3] * witness.cv.qi +
        us[4] * witness.so1v.qi + us[5] * witness.so2v.qi - r0 + challs.e * witness.zwv.qi;
    G1::mul(E, prep.pk.g1[0], temp);

    GT left, right;
    G1 wwx;
    G1::mul(wwx, witness.wwx.c, challs.e);
    pairing(left, witness.wx.c + wwx, prep.pk.g2[1]);

    G1 wx;
    G1::mul(wx, witness.wx.c, challs.v);
    G1::mul(wwx, witness.wwx.c, challs.v * omega * challs.e);
    pairing(right, wx + wwx + F - E, prep.pk.g2[0]);

    cout << "✓ Final pairing check: " << (left == right ? "PASSED" : "FAILED") << "\n";

    return left == right;
}
