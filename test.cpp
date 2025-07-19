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
     cout << "=== PLONK Protocol Test with Public Inputs ===" << endl;

     int n = 4;  // Circuit size (power of 2)
     int l = 1;  // Number of public inputs

     cout << "\nCircuit: Proving knowledge of x such that x + 5 = public_result" << endl;
     cout << "Public input: result = 8" << endl;
     cout << "Secret: x = 3" << endl;

     // Set up gate polynomials
     vector<Fr> qm(n + 1, Fr(0)), ql(n + 1, Fr(0)), qr(n + 1, Fr(0));
     vector<Fr> qo(n + 1, Fr(0)), qc(n + 1, Fr(0));
     vector<int> permutation(3 * n + 1);

     // Input constraint
     ql[1] = Fr(1);

     // a + b - c = 0
     ql[2] = Fr(1);   // coefficient of left input (x)
     qr[2] = Fr(1);   // coefficient of right input (5) 
     qo[2] = Fr(-1);  // coefficient of output (result)
     qc[2] = Fr(0);   // no constant term
     qm[2] = Fr(0);   // no multiplication

     // Identity permutation (no copy constraints for this simple example)
     for (int i = 1; i <= 3 * n; i++) {
          permutation[i] = i;
     }

     cout << "\nStep 1: Circuit Initialization..." << endl;
     Plonk::Circuit circuit = Plonk::initialize(n, qm, ql, qr, qo, qc, permutation);

     cout << "Step 2: Preprocessing (Setup)..." << endl;
     Fr setup_param;
     setup_param.setByCSPRNG();

     auto t1 = high_resolution_clock::now();
     Plonk::Preprocess preprocess = Plonk::setup(setup_param, circuit);
     auto t2 = high_resolution_clock::now();
     cout << "[✓] Preprocessing completed in " 
          << duration_cast<milliseconds>(t2 - t1).count() << " ms" << endl;

     cout << "\nStep 3: Witness Assignment..." << endl;

     // Witness vector
     vector<Fr> w(3 * n + 1, Fr(0));

     // PUBLIC INPUT
     w[1] = Fr(8); 

     // PRIVATE WITNESS
     Fr secret_x = Fr(3);  
     w[2] = secret_x;      
     w[3] = Fr(0);         
     w[4] = Fr(0);         

     // SECOND COPY 
     w[n + 1] = Fr(0);
     w[n + 2] = Fr(5);
     w[n + 3] = Fr(0);
     w[n + 4] = Fr(0);

     // THIRD COPY 
     w[2*n + 1] = Fr(0);
     w[2*n + 2] = Fr(8);
     w[2*n + 3] = Fr(0);
     w[2*n + 4] = Fr(0);

     cout << "Public input: w[1] = " << w[1].getStr() << endl;
     cout << "Secret value: w[2] = " << w[2].getStr() << endl;
     cout << "Constant: " << w[n+2].getStr() << endl;
     cout << "Result: " << w[2*n+2].getStr() << endl;

     // Constraint verification 
     Fr constraint_check = ql[2] * w[2] + qr[2] * w[n+1] + qo[2] * w[2*n+1] + qc[2];
     cout << "[✓] Gate constraint check: " << (constraint_check == 0 ? "Valid" : "Invalid") << endl;

     cout << "\nStep 4: Generating Proof..." << endl;
     t1 = high_resolution_clock::now();
     Plonk::Witness witness = Plonk::prove(preprocess, l, w);
     t2 = high_resolution_clock::now();
     cout << "[✓] Prover completed in " 
          << duration_cast<milliseconds>(t2 - t1).count() << " ms" << endl;

     cout << "\nStep 5: Verification..." << endl;

     vector<Fr> public_inputs(l + 1);  
     public_inputs[0] = Fr(0);         
     for (int i = 1; i <= l; i++) {
          public_inputs[i] = w[i];      
     }

     cout << "Verifier knows public inputs: ";
     for (int i = 1; i <= l; i++) {
          cout << "w[" << i << "] = " << public_inputs[i].getStr() << " ";
     }
     cout << endl;

     t1 = high_resolution_clock::now();
     bool result = Plonk::verify(n, preprocess, witness, l, public_inputs);
     t2 = high_resolution_clock::now();
     cout << "[✓] Verifier completed in " 
          << duration_cast<milliseconds>(t2 - t1).count() << " ms" << endl;

     cout << "\n=== RESULT ===" << endl;
     cout << "Proof verification: " << (result ? "VALID ✅" : "INVALID ❌") << endl;

     if (result) {
          cout << "\n🎉 Successfully proved knowledge of x = " << secret_x.getStr() 
               << " such that x + 5 = " << w[1].getStr() 
               << " without revealing x to the verifier!" << endl;
     }

     // Test with wrong public input to verify security
     cout << "\n=== Security Test: Wrong Public Input ===" << endl;
     vector<Fr> wrong_public(l + 1);
     wrong_public[0] = Fr(0);
     wrong_public[1] = Fr(10);  // Wrong public input

     bool wrong_result = Plonk::verify(n, preprocess, witness, l, wrong_public);
     cout << "Verification with wrong public input: " << (wrong_result ? "VALID ❌ (ERROR!)" : "INVALID ✅ (Correct)") << endl;

     return 0;
}
