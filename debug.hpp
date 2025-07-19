// This debug is used before KZG is implemented.
// Some implementations will not work as expected.

#include "ntt.hpp"
#include "kzg.hpp"
#include <mcl/bn.hpp>
#include <vector>
#include <map>
#include <cmath>

using namespace std;
using namespace mcl;
using namespace bn;

void debugSetup(Plonk::Preprocess ret, Plonk::Circuit circuit);

void debugRoundOne(Plonk::Transcript ret, Plonk::Preprocess prep, vector<Fr> w);

void debugRoundTwo(Plonk::Transcript ret, Plonk::Preprocess prep, vector<Fr> w, vector<Fr> b, Plonk::Challenge challs);

void debugRoundThree(Plonk::Transcript ret, Plonk::Preprocess prep);

void debugSplit(Plonk::Transcript ret, Plonk::Preprocess prep);

void debugRoundFive(Plonk::Transcript transcript, Plonk::Preprocess preprocess, Fr v, bool first);

void debugRoundFive(Plonk::Transcript transcript, Plonk::Preprocess preprocess, Fr v, size_t f);
