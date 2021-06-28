/*
 * @file       evaluator_Aggregation.cpp, cpp file
 * @brief      functions of encryption and decryption for tasks 2 and 3
 *
 * @author     Miran Kim
 * @copyright  GNU Pub License
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <bitset>
#include <sys/resource.h>   // check the memory
#include <iomanip>          // std::setprecision
#include <unistd.h>         // sleep
 
#include "utils.h"
#include "thread.h"
#include "evaluator_Aggregation.h"

using namespace std;
using namespace seal;


/*
 Encryption of genotype matrics for task2
 
 @param[in] geno, genotype matrix encoded by 0-method (bool[n_individuals][n_allels][n_posns])
 @param[in] n_posns, the number of positions
 @param[in] n_alleles, the number of alleles
 @param[in] nctxts, the number of ciphertexts
 @param[out] ct, ciphertexts that represent geno (ct[n_individuals][n_alleles][nctxts])
*/

void HEenc::encrypt_data(vector<vector<vector<Ciphertext>>>& ct,  vector<vector<vector<bool>>> geno, int n_posns, int n_allels, int nctxts)
{
    int nctxts_total = ct[0][0].size();
    if(nctxts == 0){
        nctxts = nctxts_total;
    }
    int n_individuals = geno.size();
    size_t slot_count = encoder.slot_count();
    long last_slot_count = (n_posns - slot_count * (nctxts_total - 1));
 
    struct rusage usage;
    long memoryscale = (1 << 20);
   
    for(int allele_i = 0; allele_i < n_allels; allele_i++){
        for(int i = 0; i < n_individuals; ++i){
            
            MT_EXEC_RANGE(nctxts, first, last);
            for(int j = first; j < last; j++){
                vector<uint64_t> input_vector(slot_count, 0ULL);
                int jst = j * slot_count;
                int kend = (j == (nctxts_total - 1)? last_slot_count: slot_count);
                
                for(int k = 0; k < kend; ++k){
                    int posn = jst + k;
                    input_vector[k] = geno[i][allele_i][posn];
                }

                Plaintext plaintext;
                encoder.encode(input_vector, plaintext);
                encryptor.encrypt(plaintext, ct[i][allele_i][j]);
            }
            MT_EXEC_RANGE_END
            
            if((i + 1) % 10 == 0){
                cout << "(i,allele) = (" << i << "," << allele_i << "): " ;
                getrusage(RUSAGE_SELF, &usage);
                cout << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
            }
        }
    }
}

/*
 Encryption of genotype matrics for task2
 
 @param[in] geno, genotype matrix encoded by 1-method (bool[2][n_individuals][n_allels][n_posns])
 @param[in] n_posns, the number of positions
 @param[in] n_alleles, the number of alleles
 @param[in] nctxts, the number of ciphertexts
 @param[out] ct, ciphertexts that represent geno (ct[n_individuals][n_alleles][nctxts])
*/

void HEenc::encrypt_data(vector<vector<vector<Ciphertext>>>& ct, vector<vector<vector<vector<bool>>>> geno, int n_posns, int n_allels, int nctxts)
{
    int nctxts_total = ct[0][0].size();
    if(nctxts == 0){
        nctxts = nctxts_total;
    }
    
    int n_individuals = geno[0].size();
    size_t slot_count = encoder.slot_count();
    long last_slot_count = (n_posns - slot_count * (nctxts_total - 1));
  
    struct rusage usage;
    long memoryscale = (1 << 20); // 2^20: linux, 2^30: mac
    
    for(int allele_i = 0; allele_i < n_allels; allele_i++){ // n_allels
        for(int i = 0; i < n_individuals; ++i){
            cout << "(i,allele) = (" << i << "," << allele_i << "): " ;
            
            //vector<vector<uint64_t>> input_temp(nctxts, vector<uint64_t>(slot_count, 0ULL));
            MT_EXEC_RANGE(nctxts, first, last);
            for(int j = first; j < last; j++){
                vector<uint64_t> input_vector(slot_count, 0ULL);
                int jst = j * slot_count;
                int kend = (j == (nctxts_total - 1)? last_slot_count: slot_count);
                
                for(int k = 0; k < kend; ++k){
                    int posn = jst + k;
                    input_vector[k] = ((geno[1][i][allele_i][posn] << 1) | (geno[0][i][allele_i][posn]));// first merge two bits and then put!
                }

                Plaintext plaintext;
                encoder.encode(input_vector, plaintext);
                encryptor.encrypt(plaintext, ct[i][allele_i][j]);
            }
            MT_EXEC_RANGE_END
          
            getrusage(RUSAGE_SELF, &usage);
            cout << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
        }
    }
}

/*
 Encryption of genotype matrics for task2
 
 @param[out] ct, ciphertexts that represent geno (ct[n_alleles][nctxts])
 @param[in] n_posns, the number of positions
 @param[in] n_alleles, the number of alleles
 @param[in] nctxts, the number of ciphertexts
 @param[out] res, plain result (res[n_alleles][n_posns])
*/
 
void HEeval::decrypt_result(vector<vector<int>>& res, vector<vector<Ciphertext>> ct, int n_posns, int n_alleles, int nctxts)
{
    int nctxts_total = ct[0].size();
    if(nctxts == 0){
        nctxts = nctxts_total;
    }
    
    size_t slot_count = encoder.slot_count();
    long last_slot_count = (n_posns - slot_count * (nctxts_total - 1));
    
    for(int allele_i = 0; allele_i < n_alleles; allele_i++){
        MT_EXEC_RANGE(nctxts, first, last);
        for(int k = first; k < last; k++){
            vector<uint64_t> buf_result;
            Plaintext plain_result;
            decryptor.decrypt(ct[allele_i][k], plain_result);
            encoder.decode(plain_result, buf_result);

            int jend = (k == (nctxts_total - 1)? last_slot_count: slot_count);
            for(int j = 0; j < jend; ++j){
                int l = k * slot_count + j;
                res[allele_i][l] = (buf_result[j]);
            }
        }
        MT_EXEC_RANGE_END
    }
}
