#include <vector>
#include "utils.h"
#include "thread.h"
#include "seal/seal.h"
 
using namespace std;
using namespace seal;

// task 1
class HEutils{
public:
    Encryptor& encryptor;
    Decryptor& decryptor;
    BatchEncoder& encoder;
    
    HEutils(Encryptor& encryptor, Decryptor& decryptor, BatchEncoder& encoder):  encryptor(encryptor), decryptor(decryptor), encoder(encoder) {}

    // Encryption of input data
    void encrypt_data(vector<Ciphertext> &res, int *variant_vector, long nctxts, long last_slot_count, int n_threads);
    
    void encode_variants(vector<vector<Plaintext>> &impact_pt, vector<vector<bool>> &indicator,
                       uint64_t *impact_vector, size_t nbitsize, size_t nwords, long nctxts, long last_slot_count);
    
    // Decryption of result
    void decrypt_result(vector<uint64_t> &res, vector<vector<Ciphertext>> res_ct,
                        vector<vector<bool>> indicator, int n_posns, size_t nbitsize, size_t nwords, long nctxts, long last_slot_count);
};


// functions of encryption and decryption for tasks 2 and 3
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
