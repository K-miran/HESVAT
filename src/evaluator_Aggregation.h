/*
 * @file       evaluator_Aggregation.h, header file
 * @brief      functions of encryption and decryption for tasks 2 and 3
 *
 * @copyright  MIT License
 */


#ifndef __EVALUATOR__
#define __EVALUATOR__

#include <vector>
#include "seal/seal.h"

using namespace std;
using namespace seal;


class HEenc{
public:
    Encryptor& encryptor;
    BatchEncoder& encoder;
  
    HEenc(Encryptor& encryptor, BatchEncoder& encoder):  encryptor(encryptor), encoder(encoder) {}
    
    // encryption of 0-encoding
    void encrypt_data(vector<vector<vector<Ciphertext>>>& ct, vector<vector<vector<bool>>> geno, int n_posns, int n_allels, const int nctxts = 0);
    
    // encryption of 1-encoding
    void encrypt_data(vector<vector<vector<Ciphertext>>>& ct, vector<vector<vector<vector<bool>>>> geno, int n_posns, int n_allels, const int nctxts = 0);
};


class HEeval{
public:
    BatchEncoder& encoder;
    Decryptor& decryptor;
    Evaluator& evaluator;
     
    HEeval(BatchEncoder& encoder, Decryptor& decryptor, Evaluator& evaluator):  encoder(encoder), decryptor(decryptor), evaluator(evaluator) {}
    
    void decrypt_result(vector<vector<int>>& res,  vector<vector<Ciphertext>> ct, int n_posns, int n_alleles,  const int nctxts = 0);
};

#endif
