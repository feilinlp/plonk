#include "kzg.hpp"
#include <mcl/bn.hpp>
#include <mcl/lagrange.hpp>

using namespace std;
using namespace mcl;
using namespace bn;

KZG::PublicKey setup(size_t t) {
    KZG::PublicKey pk;
    pk.t = t;
    pk.g1.resize(t+1);
    pk.g2.resize(t+1);

    Fr a;
    a.setByCSPRNG();
    
    // Get generators for G1 and G2
    G1 g1;
    hashAndMapToG1(g1, "G1_Generator", strlen("G1_Generator"));
    G2 g2;
    hashAndMapToG2(g2, "G2_Generator", strlen("G2_Generator"));

    Fr power = 1;
    for (size_t i = 0; i <= t; i++) {
        G1::mul(pk.g1[i], g1, power);
        G2::mul(pk.g2[i], g2, power);
        Fr::mul(power, power, a);
    }

    return pk;
}

KZG::Commitment commit(KZG::PublicKey pk, vector<Fr> q) {
    KZG::Commitment comm; 
    comm.c.clear(); 

    for (size_t i = 0; i < q.size(); i++) {
        G1 temp;
        G1::mul(temp, pk.g1[i], q[i]);
        G1::add(comm.c, comm.c, temp);
    }

    return comm;
}

// Evaluate the value of q(i)
Fr evaluatePoly(vector<Fr> q, Fr i) {
    if (q.empty()) return Fr(0);
    
    Fr result;
    evaluatePolynomial(result, q.data(), q.size(), i);
    return result;
}

vector<Fr> divideByLinear(vector<Fr> q, Fr i) {
    vector<Fr> result(q.size() - 1);

    // Synthetic division from highest degree down
    Fr carry = 0;
    for (int j = q.size() - 1; j >= 1; --j) {
        carry = carry * i + q[j];
        result[j - 1] = carry;
    }

    return result;
}

KZG::Witness createWitness(KZG::PublicKey pk, vector<Fr> q, Fr i) {
    KZG::Witness witness;
    witness.i = i;
    witness.w.clear();

    witness.qi = evaluatePoly(q, i);

    q[0] -= witness.qi;
    vector<Fr> result = divideByLinear(q, i);

    for (size_t j = 0; j < result.size(); j++) {
        G1 temp;
        G1::mul(temp, pk.g1[j], result[j]);
        G1::add(witness.w, witness.w, temp);
    }

    return witness;
}

bool verifyEval(KZG::PublicKey pk, KZG::Commitment comm, Fr i, KZG::Witness witness) {
    GT left, right1, right2;
    pairing(left, comm.c, pk.g2[0]); // e(C, g)

    G2 gi;
    G2::mul(gi, pk.g2[0], i); // g^i

    G2 temp;
    G2::sub(temp, pk.g2[1], gi); // g^a / g^i
    pairing(right1, witness.w, temp); // e(w_i, g^a / g^i)

    pairing(right2, pk.g1[0], pk.g2[0]); // e(g,g)
    
    GT right;
    GT::pow(right, right2, witness.qi);
    GT::mul(right, right1, right);

    return left == right;
}
