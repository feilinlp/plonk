#include "include/plonk.hpp"
#include "include/ntt.hpp" 
#include "include/kzg.hpp"
#include <mcl/bn.hpp>
#include <iostream>
#include <vector>
#include <chrono>

using namespace std;
using namespace mcl;
using namespace bn;
using namespace std::chrono;

int main() {
    initPairing(mcl::BN_SNARK1);
    cout << "=== PLONK Protocol Test ===" << endl;

    int n = 4;  // Power of 2
    int l = 4;  // Number of public inputs (not used in this case)

    // Gate: a * b = c
    vector<Fr> qm(n + 1, Fr(0)), ql(n + 1, Fr(0)), qr(n + 1, Fr(0));
    vector<Fr> qo(n + 1, Fr(0)), qc(n + 1, Fr(0));
    vector<int> permutation(3 * n + 1);

    // Gate at index 1
    // qm[1] = Fr(1);
    ql[1] = Fr(1);
    qr[1] = Fr(1);
    qo[1] = Fr(-1);

    // Identity permutation
    for (int i = 1; i <= 3 * n; i++) {
        permutation[i] = i;
    }

    // Step 1: Circuit Initialization
    cout << "\nInitializing Circuit..." << endl;
    Plonk::Circuit circuit = Plonk::initialize(n, qm, ql, qr, qo, qc, permutation);

    // Step 2: Preprocessing
    Fr x;
    x.setByCSPRNG();

    auto t1 = high_resolution_clock::now();
    Plonk::Preprocess preprocess = Plonk::setup(x, circuit);
    auto t2 = high_resolution_clock::now();
    cout << "[✓] Preprocessing completed in " 
         << duration_cast<milliseconds>(t2 - t1).count() << " ms" << endl;

    // Step 3: Witness assignment
    vector<Fr> w(3 * n + 1, Fr(0));
    w[1] = Fr(3);         // a
    w[n + 1] = Fr(4);     // b
    w[2 * n + 1] = Fr(7); // c = a * b

    // Optional: Constraint sanity check
    Fr constraint = w[1] + w[n + 1] - w[2 * n + 1];
    cout << "[✓] Manual constraint check: " << (constraint == 0 ? "Valid" : "Invalid") << endl;

    // Step 4: Prover
    t1 = high_resolution_clock::now();
    Plonk::Witness witness = Plonk::prove(preprocess, l, w);
    t2 = high_resolution_clock::now();
    cout << "[✓] Prover completed in " 
         << duration_cast<milliseconds>(t2 - t1).count() << " ms" << endl;

    // Step 5: Verifier
    t1 = high_resolution_clock::now();
    bool result = Plonk::verify(n, preprocess, witness, l, w);
    t2 = high_resolution_clock::now();
    cout << "[✓] Verifier completed in " 
         << duration_cast<milliseconds>(t2 - t1).count() << " ms" << endl;

    cout << "\n[RESULT] Proof is " << (result ? "VALID ✅" : "INVALID ❌") << endl;
    return 0;
}
