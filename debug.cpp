// This debug is used before KZG is implemented.
// Some implementations will not work as expected.

#include "plonk.hpp"
#include "ntt.hpp"
#include "kzg.hpp"
#include "debug.hpp"
#include <mcl/bn.hpp>
#include <vector>
#include <map>
#include <cmath>

using namespace std;
using namespace mcl;
using namespace bn;

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

void debugRoundThree(Plonk::Transcript ret, Plonk::Preprocess prep) {
    size_t n = prep.circuit.n;

    cout << "Checking t values...\n";

    for (size_t i = 1; i <= n; i++) {
        cout << "At " << i << " ";
        cout << evaluatePoly(ret.t, prep.circuit.h[i]) << endl;
    }
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
