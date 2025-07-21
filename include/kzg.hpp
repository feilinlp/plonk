#ifndef KZG_HPP
#define KZG_HPP

#include <mcl/bn.hpp>

using namespace mcl;
using namespace bn;
using namespace std;

class KZG {
public:
    struct PublicKey {
        vector<G1> g1; 
        vector<G2> g2;
        size_t t; 
    };

    struct Commitment {
        G1 c;
    };

    struct Witness {
        Fr i;
        G1 w; // Witness
        Fr qi; // Evaluated value
    };

    struct BatchItem {
        Commitment commitment;
        Fr point;
        Fr value;
        Witness witness;
    };
};

KZG::PublicKey KZGSetup(size_t t, Fr x);

KZG::Commitment commit(KZG::PublicKey pk, vector<Fr> q);

Fr evaluatePoly(vector<Fr> q, Fr i);

vector<Fr> divideByLinear(vector<Fr> q, Fr i);

KZG::Witness createWitness(KZG::PublicKey pk, vector<Fr> q, Fr i);

bool verifyEval(KZG::PublicKey pk, KZG::Commitment comm, Fr i, KZG::Witness witness);

vector<Fr> generateRandomCoeffs(size_t count);

bool verifyBatch(KZG::PublicKey pk, const vector<KZG::BatchItem>& batch_items, const vector<Fr>& random_coeffs);

#endif // KZG_HPP
