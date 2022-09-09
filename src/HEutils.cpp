#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "math.h"
#include "time.h"
#include <chrono>
#include <sys/resource.h>
#include <omp.h>

#include "utils.h"
#include "thread.h"
#include "HEutils.h"

using namespace std;

void HEutils::encrypt_data(vector<Ciphertext> &res, 
                           int *variant_vector, long nctxts, long last_slot_count, int n_threads)
{
    size_t slot_count = encoder.slot_count();
    res.resize(nctxts);
    
    // openmp version (0.8sec)
//#pragma omp parallel num_threads(n_threads) // 0.8 sec
//    {
//#pragma omp for
//        for(int i = 0; i < nctxts; i++){
//            //vector<uint64_t> input_vector(slot_count, 0ULL);
//            vector<uint64_t> input_vector;
//
//            int jst = i * slot_count;
//            int jend = (i == (nctxts - 1)? jst + last_slot_count: jst + slot_count);
//
//            input_vector.assign(variant_vector + jst, variant_vector + jend);
//
//            Plaintext plaintext;
//            encoder.encode(input_vector, plaintext);
//            encryptor.encrypt(plaintext, res[i]);
//        }
//    }
    
    
    // MT version (0.08sec)
    MT_EXEC_RANGE(nctxts, first, last);
    for(int i = first; i < last; i++){
        vector<uint64_t> input_vector;

        int jst = i * slot_count;
        int jend = (i == (nctxts - 1)? jst + last_slot_count: jst + slot_count);

        input_vector.insert(input_vector.end(), variant_vector + jst, variant_vector + jend);

        Plaintext plaintext;
        encoder.encode(input_vector, plaintext);
        encryptor.encrypt(plaintext, res[i]);
    }
    MT_EXEC_RANGE_END
    
}

void HEutils::encode_variants(vector<vector<Plaintext>> &impact_pt, vector<vector<bool>> &indicator,
                            uint64_t *impact_vector, size_t nbitsize, size_t nwords, long nctxts, long last_slot_count)
{
    size_t slot_count = encoder.slot_count();
    uint64_t unit = (1 << nbitsize) - 1;
    impact_pt.resize(nctxts, vector<Plaintext>(nwords));
    indicator.resize(nctxts, vector<bool>(nwords, true));
    
    
    MT_EXEC_RANGE(nctxts, first, last);
    for(int i = first; i < last; i++){
        int unit_slot_cnt = (i == (nctxts - 1)? last_slot_count: slot_count);
        vector<vector<uint64_t>> buf_vector(nwords, vector<uint64_t> (unit_slot_cnt));
        
        int jst = i * slot_count;
        int jend = jst + unit_slot_cnt;
        
        long l = 0;
        for(long j = jst; j < jend; ++j){
            uint64_t val = impact_vector[j];
            // w-bit representation
            for(int k = 0; k < nwords; ++k){
                buf_vector[k][l] = ((val >> (k * nbitsize)) & unit);
            }
            l++;
        }
        
        for(int k = 0; k < nwords; ++k){
            encoder.encode(buf_vector[k], impact_pt[i][k]);
            
            // check if the entries of buf_vector are all zero
            int j = 0;
            while(1){
                if(buf_vector[k][j] != 0){
                    break;
                }
                else{ // buf_vector[k][j] = 0
                    j++;
                    if(j == unit_slot_cnt){
                        indicator[i][k] = false;
                        break;
                    }
                }
            }
        }
    }
    MT_EXEC_RANGE_END
}


void HEutils::decrypt_result(vector<uint64_t> &res, vector<vector<Ciphertext>> res_ct,
                   vector<vector<bool>> indicator, int n_posns, size_t nbitsize, size_t nwords, long nctxts, long last_slot_count)
{
    size_t slot_count = encoder.slot_count();
    res.resize(n_posns, 0ULL);
    vector<vector<vector<uint64_t>>> buf_result(nctxts, vector<vector<uint64_t>> (nwords, vector<uint64_t> (slot_count)));
    
    MT_EXEC_RANGE(nctxts, first, last);
    for(int i = first; i < last; i++){
        for(int k = 0; k < nwords; ++k){
            if(indicator[i][k]){
                Plaintext plain_result;
                decryptor.decrypt(res_ct[i][k], plain_result);
                encoder.decode(plain_result, buf_result[i][k]);
            }
        }
        
        // buf_result[i][k]: vector
        int jend = (i == (nctxts - 1)? last_slot_count: slot_count);
        for(int j = 0; j < jend; ++j){
            int l = i * slot_count + j;
            for(int k = 0; k < nwords; ++k){
                res[l] |= (buf_result[i][k][j] << (nbitsize * k));
            }
        }
    }
    MT_EXEC_RANGE_END
}


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
         
#ifdef DEBUG
            if((i + 1) % 10 == 0){
                cout << "(i,allele) = (" << i << "," << allele_i << "): " ;
                getrusage(RUSAGE_SELF, &usage);
                cout << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
            }
#endif
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
#ifdef DEBUG
            cout << "(i,allele) = (" << i << "," << allele_i << "): " ;
#endif
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
