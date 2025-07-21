#include "../include/kzg.hpp"
#include <mcl/bn.hpp>
#include <mcl/lagrange.hpp>

using namespace std;
using namespace mcl;
using namespace bn;

KZG::PublicKey KZGSetup(size_t t, Fr x) {
    KZG::PublicKey pk;
    pk.t = t;

    G1 one;
    mapToG1(one, 1);
    G1::mul(one, one, 1);

    Fr curr = 1;
    for (size_t i = 0; i <= t + 5; i++) { 
        G1 ins;
        G1::mul(ins, one, curr);  
        
        pk.g1.push_back(ins);
        
        Fr::mul(curr, curr, x);
    }

    G2 two;
    mapToG2(two, 1);

    pk.g2.push_back(two);
    G2::mul(two, two, x);
    pk.g2.push_back(two);

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
    int n = q.size();
    vector<Fr> result(n - 1);
    Fr carry = 0;

    for (int j = n - 1; j >= 0; --j) {
        carry = carry * i + q[j];
        if (j > 0)
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

vector<Fr> generateRandomCoeffs(size_t count) {
    vector<Fr> coeffs(count);
    
    for (size_t i = 0; i < count; i++) {
        coeffs[i].setByCSPRNG();
    }
    
    return coeffs;
}

// Batch verification with provided random coefficients
bool verifyBatch(KZG::PublicKey pk, const vector<KZG::BatchItem>& batch_items, const vector<Fr>& coeffs) {
    if (batch_items.empty()) return true;
    if (batch_items.size() != coeffs.size()) return false;
    
    size_t n = batch_items.size();
    
    G1 commitments;
    commitments.clear();
    
    for (size_t i = 0; i < n; i++) {
        G1 temp;
        G1::mul(temp, batch_items[i].commitment.c, coeffs[i]);
        G1::add(commitments, commitments, temp);
    }
    
    GT left;
    pairing(left, commitments, pk.g2[0]);
    
    // Compute the evaluation part: e(g, g)^(∑_i γ_i * v_i)
    Fr evaluations;
    evaluations.clear();
    
    for (size_t i = 0; i < n; i++) {
        Fr temp;
        Fr::mul(temp, coeffs[i], batch_items[i].value);
        Fr::add(evaluations, evaluations, temp);
    }
    
    GT evaluation_pairing;
    pairing(evaluation_pairing, pk.g1[0], pk.g2[0]);
    GT::pow(evaluation_pairing, evaluation_pairing, evaluations);
    
    // Divide left side by evaluation part
    GT::div(left, left, evaluation_pairing);
    
    // Compute right side: e(∑_i γ_i * w_i, g^α) / e(∑_i γ_i * w_i * z_i, g)
    G1 combined_witness;
    combined_witness.clear();
    
    G1 combined_witness_points;
    combined_witness_points.clear();
    
    for (size_t i = 0; i < n; i++) {
        // Add γ_i * w_i to combined witness
        G1 temp_witness;
        G1::mul(temp_witness, batch_items[i].witness.w, coeffs[i]);
        G1::add(combined_witness, combined_witness, temp_witness);
        
        // Add γ_i * w_i * z_i to combined witness points
        Fr temp_coeff;
        Fr::mul(temp_coeff, coeffs[i], batch_items[i].point);
        G1 temp_point;
        G1::mul(temp_point, batch_items[i].witness.w, temp_coeff);
        G1::add(combined_witness_points, combined_witness_points, temp_point);
    }
    
    // e(∑_i γ_i * w_i, g^α)
    GT right1;
    pairing(right1, combined_witness, pk.g2[1]);
    
    // e(∑_i γ_i * w_i * z_i, g)
    GT right2;
    pairing(right2, combined_witness_points, pk.g2[0]);
    
    // Compute right side
    GT right;
    GT::div(right, right1, right2);
    
    return left == right;
}
