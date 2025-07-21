#include "plonk.hpp"
#include "ntt.hpp"
#include "kzg.hpp"
#include "debug.hpp"
#include <openssl/evp.h>
#include <sstream>
#include <mcl/bn.hpp>
#include <vector>
#include <map>
#include <cmath>

using namespace std;
using namespace mcl;
using namespace bn;

vector<uint8_t> toBytes(const Fr &x);

vector<uint8_t> toBytes(const G1 &g);

Fr hashToFr(const vector<vector<uint8_t>>& inputs);

bool isQuadraticResidue(Fr w);

Fr quadraticNonResidue();

bool inCoset(vector<Fr> h, Fr k2);

vector<vector<Fr>> lagrangeBasis(size_t n, Fr omega);

Fr evaluatePolynomial(vector<Fr> p, vector<Fr> x);

vector<Fr> polynomialDivision(vector<Fr> a, size_t n);

vector<Fr> addPolynomials(vector<Fr> a, vector<Fr> b);

vector<KZG::BatchItem> addItems(Plonk::Preprocess prep, Plonk::Witness witness, Plonk::Verifier verifier, Plonk::Challenge challs);
