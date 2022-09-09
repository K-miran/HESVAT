/*
 * @file       testSVAT.cpp, cpp file
 * @brief      function for running secure outsourcing of targeted variant annotation and genotype aggregation
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
#include <omp.h>            // openmp

#include "seal/seal.h"
#include "utils.h"
#include "thread.h"
#include "HEutils.h"
#include "testPVAT.h"

#define DEBUG true
#define omp false
#define precompute true

using namespace std;
using namespace seal;


/*
 @Task1: Perform the secure variant annotation
 @param[in] variant_fp, filename of an annotation vector
 @param[in] impact_fp, filename of a variant loci vector
 @param[in] int n_posns, the number of positions
*/
void run_secure_variant_annotation(vector<uint64_t> &res, const char *variant_fp, const char *impact_fp, int n_posns, long nvariant)
{
    size_t poly_modulus_degree = (1 << 11);
    size_t plaintext_modulus_size = 19;
    
    size_t nwords = 3;
    size_t nbitsize = 19;   // 2^w-bit repn
    
    // condition1: nwords * nbitsize < 56
    // condition2: nbitsize < plaintet_modulus_size

    chrono::high_resolution_clock::time_point time_start;
    chrono::high_resolution_clock::time_point time_end;
     
    long memoryscale = (1 << 20); // 2^20: linux, 2^30: mac
    struct rusage usage;
    
    /* Load data (annotation vector, variant loci vector ) */
    
    int *variant_vector = new int[n_posns];
    FILE *ptr;
    ptr = fopen(variant_fp, "rb");
    fread(variant_vector, sizeof(int), n_posns, ptr);
    
    uint64_t *impact_vector = new uint64_t[n_posns];
    FILE *ptr1;
    ptr1 = fopen(impact_fp, "rb");
    fread(impact_vector, sizeof(uint64_t), n_posns, ptr1);
    
    getrusage(RUSAGE_SELF, &usage);
    cout << "Load Data : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
   
    cout << "+----------------------------------------+" << endl;
    cout << "|             Key Generation             |" << endl;
    cout << "+----------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    /*  Key Setup for HE computation  */
    EncryptionParameters parms(scheme_type::bfv);
    
    parms.set_poly_modulus_degree(poly_modulus_degree);
    
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, plaintext_modulus_size));

    SEALContext context(parms);
    
    auto qualifiers = context.first_context_data()->qualifiers();
     
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    
    auto &context_data = *context.key_context_data();
   
   /* Print the size of the true (product) coefficient modulus */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << endl;
    std::cout << "|   plain_modulus size: " << log2(context_data.parms().plain_modulus().value()) << std::endl;
    
    time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Scheme generation (milliseconds) : " << time_diff.count()/1000.0 << " ms] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    cout << "+----------------------------------------+" << endl;
    cout << "|     Encryption (annotation vector)     |" << endl;
    cout << "+----------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    BatchEncoder batch_encoder(context);
    HEutils heutils(encryptor, decryptor, batch_encoder);
    
    size_t slot_count = batch_encoder.slot_count();
    long nctxts = (ceil)((n_posns/(double)slot_count));              // the number of generated new ciphertexts (=1665 in this dataset)
    long last_slot_count = (n_posns - slot_count * (nctxts - 1));   // the number of plaintext slots in the last ciphertext
    
    vector<Ciphertext> variant_ct; // [nctxts]
    int n_threads = Thread::availableThreads(); // omp_get_max_threads();
    
    heutils.encrypt_data(variant_ct, variant_vector, nctxts, last_slot_count, n_threads);
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Encryption (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

    cout << "+----------------------------------------+" << endl;
    cout << "|          Encode (variant loci)         |" << endl;
    cout << "+----------------------------------------+" << endl;

    time_start = chrono::high_resolution_clock::now();

    vector<vector<Plaintext>> impact_pt;    // the generated impact plaintext polynomial
    vector<vector<bool>> indicator;     // indicator vector where the plaintext polynomial is zero or not (if zero, we can skip plaintext-ciphertext multiplication)
     
    heutils.encode_variants(impact_pt, indicator, impact_vector, nbitsize, nwords, nctxts, last_slot_count);
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Encode (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    cout << "+----------------------------------------+" << endl;
    cout << "|               Evaluation               |" << endl;
    cout << "+----------------------------------------+" << endl;

    time_start = chrono::high_resolution_clock::now();
    
    // multiplication with variant_ct and impact_pt
    vector<vector<Ciphertext>> res_ct(nctxts, vector<Ciphertext>(nwords));
    MT_EXEC_RANGE(nctxts, first, last);
    for(int i = first; i < last; i++){
        for(int k = 0; k < nwords; ++k){
            if(indicator[i][k]){
                evaluator.multiply_plain(variant_ct[i], impact_pt[i][k], res_ct[i][k]);
            }
        }
    }
    MT_EXEC_RANGE_END
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Evaluation (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "Noise budget in result: " << decryptor.invariant_noise_budget(res_ct[0][0]) << " bits" << endl;
    
    cout << "+----------------------------------------+" << endl;
    cout << "|              Decryption                |" << endl;
    cout << "+----------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    heutils.decrypt_result(res, res_ct, indicator, n_posns, nbitsize, nwords, nctxts, last_slot_count);
       
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Decryption (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
}



// DEL or INS
void run_secure_variant_annotation_large(vector<vector<uint64_t>> &res, string dataset_str, vector<vector<string>> chr_info, string target)
{
    chrono::high_resolution_clock::time_point time_start;
    chrono::high_resolution_clock::time_point time_end;
    
    long memoryscale = (1 << 20); // 2^20: linux, 2^30: mac
    struct rusage usage;
    
    // Step 0. Read the data
    int nchr = chr_info.size(); // number of chromosome
    //nchr = 2;
    
    string variant_folder = dataset_str + "VARIANT_SIGNAL/";      // the folder of annotation vectors
    string impact_folder = dataset_str + "ANNOTATION_SIGNAL/";     // the folder of annotation signals
    
    // vector<vector<int>> variant_vector1 (nchr);
    res.resize(nchr);
    int **variant_vector = new int*[nchr];               // variant_vector[nchr][n_posns]
    uint64_t **impact_vector = new uint64_t*[nchr];     // impact_vector[nchr][n_posns]
    string *chr_id = new string[nchr];
    long *n_posns = new long[nchr];
    
    for(int ch = 0; ch < nchr; ch++){
        chr_id[ch] = chr_info[ch][0];// chromosome id
        n_posns[ch] = (long) atoi(chr_info[ch][1].c_str());  // number of variants
        if(target == "SNV"){
            n_posns[ch] *= 4; 
        }
        string variant_str = variant_folder + chr_id[ch] + ".bin"; // file for the annotation vector
        string impact_str = impact_folder + chr_id[ch] + ".bin";   // file for the impact vector
        
        const char *variant_fp = variant_str.c_str();
        const char *impact_fp = impact_str.c_str();
        
        // if you want plain computation ...
        //run_secure_variant_annotation(res[ch], variant_fp, impact_fp, n_posns[ch], 10000);
        
        /* Load data (annotation vector, variant loci vector ) */
        variant_vector[ch] = new int[n_posns[ch]];              // to be encrypted
        FILE *ptr;
        ptr = fopen(variant_fp, "rb");
        fread(variant_vector[ch], sizeof(int), n_posns[ch], ptr);

        impact_vector[ch] = new uint64_t[n_posns[ch]];// not to be encrypted
        FILE *ptr1;
        ptr1 = fopen(impact_fp, "rb");
        fread(impact_vector[ch], sizeof(uint64_t), n_posns[ch], ptr1);

        getrusage(RUSAGE_SELF, &usage);
        cout << "Load Data " + chr_id[ch] << ".bin (" << n_posns[ch] << "): " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    }
    

    cout << "+----------------------------------------+" << endl;
    cout << "|             Key Generation             |" << endl;
    cout << "+----------------------------------------+" << endl;
    time_start = chrono::high_resolution_clock::now();

    size_t poly_modulus_degree = (1 << 11);
    size_t plaintext_modulus_size = 19;

    size_t nwords = 3;
    size_t nbitsize = 19;   // 2^w-bit repn

    // condition1: nwords * nbitsize < 56
    // condition2: nbitsize < plaintet_modulus_size

    /*  Key Setup for HE computation  */
    EncryptionParameters parms(scheme_type::bfv);

    parms.set_poly_modulus_degree(poly_modulus_degree);

    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, plaintext_modulus_size));

    SEALContext context(parms);

    auto qualifiers = context.first_context_data()->qualifiers();

    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    auto &context_data = *context.key_context_data();

   /* Print the size of the true (product) coefficient modulus */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << endl;
    std::cout << "|   plain_modulus size: " << log2(context_data.parms().plain_modulus().value()) << std::endl;

    time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Scheme generation (milliseconds) : " << time_diff.count()/1000.0 << " ms] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

    cout << "+----------------------------------------+" << endl;
    cout << "|     Encryption (annotation vector)     |" << endl;
    cout << "+----------------------------------------+" << endl;

    time_start = chrono::high_resolution_clock::now();

    BatchEncoder batch_encoder(context);
    HEutils heutils(encryptor, decryptor, batch_encoder);

    size_t slot_count = batch_encoder.slot_count();
    int n_threads = Thread::availableThreads(); // omp_get_max_threads();
    long *nctxts = new long[nchr];
    vector<vector<Ciphertext>> variant_ct(nchr); // [nchr][nctxts]
    int jst, jend;
    long *last_slot_count = new long[nchr];
    long total_nctxts = 0;
    
    for(int ch = 0; ch < nchr; ch++){
        nctxts[ch] = (ceil)((n_posns[ch]/(double)slot_count));
        cout << chr_id[ch] << ".id (" << nctxts[ch] << ") ... " ;
        last_slot_count[ch] = (n_posns[ch] - slot_count * (nctxts[ch] - 1));
        heutils.encrypt_data(variant_ct[ch], variant_vector[ch], nctxts[ch], last_slot_count[ch], n_threads);
        total_nctxts += nctxts[ch];
    }
    cout  << "encryption done" << endl;
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Encryption "  << nchr << " chrs/"<< total_nctxts << " cts (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

    cout << "+----------------------------------------+" << endl;
    cout << "|          Encode (variant loci)         |" << endl;
    cout << "+----------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();

    vector<vector<vector<Plaintext>>> impact_pt(nchr); // [nchr][nctxts][nwords]
    vector<vector<vector<bool>>> indicator(nchr);       // [nchr][nctxts][nwords]
    
    for(int ch = 0; ch < nchr; ch++){
        cout << chr_id[ch] << ".id ... " ;
        heutils.encode_variants(impact_pt[ch], indicator[ch],
                              impact_vector[ch], nbitsize, nwords, nctxts[ch], last_slot_count[ch]);
    }
    cout << "encode done " << endl;

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Encode (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;


    cout << "+----------------------------------------+" << endl;
    cout << "|               Evaluation               |" << endl;
    cout << "+----------------------------------------+" << endl;

    time_start = chrono::high_resolution_clock::now();

    vector<vector<vector<Ciphertext>>> res_ct (nchr); //[nchr][ctxts][nwords]

    for(int ch = 0; ch < nchr; ch++){
        res_ct[ch].resize(nctxts[ch]);

        MT_EXEC_RANGE(nctxts[ch], first, last);
        for(int i = first; i < last; i++){
            res_ct[ch][i].resize(nwords);
            for(int k = 0; k < nwords; ++k){
                if(indicator[ch][i][k]){
                    evaluator.multiply_plain(variant_ct[ch][i], impact_pt[ch][i][k], res_ct[ch][i][k]);
                }
            }
        }
        MT_EXEC_RANGE_END
    }


    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Evaluation (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "Noise budget in result: " << decryptor.invariant_noise_budget(res_ct[0][0][0]) << " bits" << endl;

    cout << "+----------------------------------------+" << endl;
    cout << "|              Decryption                |" << endl;
    cout << "+----------------------------------------+" << endl;

    time_start = chrono::high_resolution_clock::now();

    res.resize(nchr); // [nchr][n_posns]

    for(int ch = 0; ch < nchr; ch++){
        heutils.decrypt_result(res[ch], res_ct[ch], indicator[ch], n_posns[ch], nbitsize, nwords, nctxts[ch], last_slot_count[ch]);
    }

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Decryption (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
}

/*
 @Task2: Perform the secure variant aggregation
 @param[in] n_posns, the number of positions
 @param[in] encoding, genotype encoding method
 "0" indicates the "variant existence" is stored in the matrix as 0/1 values.  This file contains genotypes encoded with 1 bit per individual.
 "1" indicates "actual genotype" is stored as 0/1/2 values. This file is encoded with 2 bits per individual.
 @param[in] n_individuals, the number of individuals
 @param[in] n_alleles_trials, the number of alleles
*/

void run_secure_variant_aggregation(int n_posns, int encoding, int n_individuals, int n_alleles_trials)
{
    chrono::high_resolution_clock::time_point time_start;
    chrono::high_resolution_clock::time_point time_end;
    
    long memoryscale = (1 << 20); // 2^20: linux, 2^30: mac
    struct rusage usage;
    
    int n_alleles = 5;
   
    time_start = chrono::high_resolution_clock::now();
    
    /* Load data */
    vector<vector<vector<bool>>> geno0;
    vector<vector<vector<vector<bool>>>> geno1;
    
    if(encoding == 0){
        read_genotype_matrix(geno0, n_posns, encoding, n_individuals, n_alleles);
    } else if(encoding == 1){
        read_genotype_matrix(geno1, n_posns, encoding, n_individuals, n_alleles);
    }
   
    time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Load data (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|           Key Generation           |" << endl;
    cout << "+------------------------------------+" << endl;
    
    size_t poly_modulus_degree = (1 << 11);
    size_t plaintext_modulus_size = 14; // 14, 16
    
    /*  Key Setup for HE computation */
    EncryptionParameters parms(scheme_type::bfv);
    
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, plaintext_modulus_size));

    SEALContext context(parms);
    
    auto qualifiers = context.first_context_data()->qualifiers();
     
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    
    auto &context_data = *context.key_context_data();
   
    /* Print the size of the true (product) coefficient modulus */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << endl;
    std::cout << "|   plain_modulus size: " << log2(context_data.parms().plain_modulus().value()) << std::endl;
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Scheme generation (milliseconds) : " << time_diff.count()/1000.0 << " ms] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|         Encryption (data)          |" << endl;
    cout << "+------------------------------------+" << endl;
    
    BatchEncoder batch_encoder(context);
    HEenc heenc(encryptor, batch_encoder);
    
    size_t slot_count = batch_encoder.slot_count();
    long nctxts = (ceil)((n_posns/(double)slot_count)); // 3224 in this experiment
    long last_slot_count = (n_posns - slot_count * (nctxts - 1));
    int nctxts_trial = nctxts; // nctxts = 3224
    
    std::cout << "|   nindividuals=" << n_individuals << ", nallels=" << n_alleles_trials << ", nposns=" << n_posns << endl;
    std::cout << "|   nslots=" << slot_count << ", nlastslot=" << last_slot_count << ", nctxts=(nvariants/slot)=" << nctxts << ", Total ctxts(5 alleles) = " << n_individuals * n_alleles * nctxts << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    vector<vector<vector<Ciphertext>>> variant_ct(n_individuals, vector<vector<Ciphertext>> (n_alleles, vector<Ciphertext> (nctxts)));
    std::cout << "|   ct.size = [" << n_individuals << "," << n_alleles << "," << nctxts << "]" << endl;
    
    if(encoding == 0){
        heenc.encrypt_data(variant_ct, geno0, n_posns, n_alleles_trials, nctxts_trial);
    } else if (encoding == 1){
        heenc.encrypt_data(variant_ct, geno1, n_posns, n_alleles_trials, nctxts_trial);
    }

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Encryption (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "Modulus chain index for encryption: q[" << context.get_context_data(variant_ct[0][0][0].parms_id()) ->chain_index() << "]" << endl;
    cout << "Noise budget in fresh encryption: " << decryptor.invariant_noise_budget(variant_ct[0][0][0]) << " bits" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;

    HEeval heeval(batch_encoder, decryptor, evaluator);
    
    time_start = chrono::high_resolution_clock::now();

    vector<vector<Ciphertext>> res_ct(n_alleles_trials, vector<Ciphertext> (nctxts));
    for(int allele_i = 0; allele_i < n_alleles_trials; allele_i++){
        MT_EXEC_RANGE(nctxts_trial, first, last);
        for(int k = first; k < last; k++){
            res_ct[allele_i][k] = variant_ct[0][allele_i][k];
            for(int i = 1; i < n_individuals; ++i){
                evaluator.add_inplace(res_ct[allele_i][k], variant_ct[i][allele_i][k]);
            }
        }
        MT_EXEC_RANGE_END
    }
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Evaluation (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "Modulus chain index for encryption: q[" << context.get_context_data(res_ct[0][0].parms_id()) ->chain_index() << "]" << endl;
    cout << "Noise budget in result: " << decryptor.invariant_noise_budget(res_ct[0][0]) << " bits" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|            Decryption              |" << endl;
    cout << "+------------------------------------+" << endl;

    time_start = chrono::high_resolution_clock::now();

    vector<vector<int>> res (n_alleles, vector<int>(n_posns, 0ULL));
    heeval.decrypt_result(res, res_ct, n_posns, n_alleles_trials, nctxts_trial);
  
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Decryption (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
 
#if (DEBUG)
    cout << "+----------------------------------------+" << endl;
    cout << "|              Correctness               |" << endl;
    cout << "+----------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    string filename = "data/Aggregation_Matrix_Data/task2_res_plain_"  + to_string(encoding) + "_" + to_string(n_individuals) + ".txt";
    vector<vector<int>> plain_geno_res;   //[n_posns][n_alleles]
    
    ifstream openFile(filename.data());
    if(openFile.is_open()) {
        read_data(plain_geno_res, filename, n_alleles_trials);
        
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        getrusage(RUSAGE_SELF, &usage);
        cout << "Read Precomputed results (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    } else{
        if(encoding == 0){
            colSums(plain_geno_res, geno0, n_individuals, n_alleles_trials, n_posns);// n_posns
        } else if(encoding == 1){
            colSums(plain_geno_res, geno1, n_individuals, n_alleles_trials, n_posns);// n_posns
        }
        write_data(plain_geno_res, filename);
        cout << "Calculate results and store them (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    }

    int nfail = 0;
    for(int allele_i = 0; allele_i < n_alleles_trials; ++allele_i){
        int jend = ((nctxts > nctxts_trial)? nctxts_trial * slot_count: n_posns);
        for(int j = 0; j < jend; ++j){
            if(plain_geno_res[j][allele_i] != res[allele_i][j]){
                cout << j << "-pos: incorrect (" << plain_geno_res[j][allele_i] << ", " << res[allele_i][j] << endl;
                nfail++;
            }
        }
    }
    
    if(nfail != 0){
        cout << "> " << nfail << " fails" << endl;
    } else{
        cout << "> All Passed" << endl;
    }
#endif
}


/*
 @Task3: Perform the secure variant aggregation with re-encryption

 @param[in] n_posns, the number of positions
 @param[in] encoding, genotype encoding method
 "0" indicates the "variant existence" is stored in the matrix as 0/1 values.  This file contains genotypes encoded with 1 bit per individual.
 "1" indicates "actual genotype" is stored as 0/1/2 values. This file is encoded with 2 bits per individual.
 @param[in] n_alleles_trials, the number of alleles
*/

void run_secure_variant_reencrypt_aggregation(int n_posns, int encoding, int n_alleles_trials)
{
    chrono::high_resolution_clock::time_point time_start;
    chrono::high_resolution_clock::time_point time_end;
    
    long memoryscale = (1 << 20); // 2^20: linux, 2^30: mac
    struct rusage usage;
    
    int n_alleles = 5;
    int n_individuals_0 = 500;
    int n_individuals_1 = 1000;
     
    time_start = chrono::high_resolution_clock::now();
    
    /* Load data */
    vector<vector<vector<bool>>> geno0_500;     // 0-enc, 500 samples
    vector<vector<vector<bool>>> geno0_1000;    // 0-enc, 1000 samples
    
    vector<vector<vector<vector<bool>>>> geno1_500;     // 1-enc, 500 samples
    vector<vector<vector<vector<bool>>>> geno1_1000;    // 1-enc, 1000 samples
    
    if(encoding == 0){
        read_genotype_matrix(geno0_500, n_posns, encoding, n_individuals_0, n_alleles);
        read_genotype_matrix(geno0_1000, n_posns, encoding, n_individuals_1, n_alleles);
    } else if(encoding == 1){
        read_genotype_matrix(geno1_500, n_posns, encoding, n_individuals_0, n_alleles);
        read_genotype_matrix(geno1_1000, n_posns, encoding, n_individuals_1, n_alleles);
    }
   
    time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Load data (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    cout << "+-------------------------------------------+" << endl;
    cout << "|               Key Generation              |" << endl;
    cout << "+-------------------------------------------+" << endl;
    
    size_t poly_modulus_degree = (1 << 12);
    size_t plaintext_modulus_size = 17;
    
    /*  Key Setup for HE computation */
    EncryptionParameters parms(scheme_type::bfv);

    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {35, 40}));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, plaintext_modulus_size));

    SEALContext context(parms);

    auto qualifiers = context.first_context_data()->qualifiers();
    
    /* Generate a key pair (sk, pk) */
    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
    Evaluator evaluator(context);
    
    /* Generate a key pair (sk0, pk0) */
    KeyGenerator keygen_0(context);
    auto secret_key_0 = keygen_0.secret_key();
    PublicKey public_key_0;
    keygen_0.create_public_key(public_key_0);
    Encryptor encryptor_0(context, public_key_0);
    Decryptor decryptor_0(context, secret_key_0);
    
    /* Generate a key pair (sk1, pk1) */
    KeyGenerator keygen_1(context);
    auto secret_key_1 = keygen_1.secret_key();
    PublicKey public_key_1;
    keygen_1.create_public_key(public_key_1);
    Encryptor encryptor_1(context, public_key_1);
    Decryptor decryptor_1(context, secret_key_1);
    
    /* Key-switchinkg keys */
    // switch key from secret_key_* to secret_key
    RelinKeys switch_keys_0;
    keygen.create_kswitch_key(secret_key_0, switch_keys_0);

    RelinKeys switch_keys_1;
    keygen.create_kswitch_key(secret_key_1, switch_keys_1);

    /* Print the size of the true (product) coefficient modulus */
    auto &context_data = *context.key_context_data();
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << endl;
    std::cout << "|   plain_modulus size: " << log2(context_data.parms().plain_modulus().value()) << std::endl;
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Scheme generation (milliseconds) : " << time_diff.count()/1000.0 << " ms] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
 
    cout << "+-------------------------------------------+" << endl;
    cout << "|             Encryption (data)             |" << endl;
    cout << "+-------------------------------------------+" << endl;
    
    BatchEncoder batch_encoder(context);
    HEenc heenc0(encryptor_0, batch_encoder);
    HEenc heenc1(encryptor_1, batch_encoder);
    
    size_t slot_count = batch_encoder.slot_count();
    long nctxts = (ceil)((n_posns/(double)slot_count)); // 3224
    long last_slot_count = (n_posns - slot_count * (nctxts - 1));
    
    int nctxts_trial = nctxts; // 1662

    
    std::cout << "|   nindividuals=" << n_individuals_0 << ", nallels=" << n_alleles_trials << ", nposns=" << n_posns << endl;
    std::cout << "|   nslots=" << slot_count << ", nlastslot=" << last_slot_count << ", nctxts=(nvariants/slot)=" << nctxts << ", Total ctxts-0(5 alleles) = " << n_individuals_0 * n_alleles * nctxts <<  ", Total ctxts-1(5 alleles) = " << n_individuals_1 * n_alleles * nctxts << endl;
    
    // Input vector<vector<vector<uint8_t>>> geno[n_individuals][5][n_posns]
    /* Encryption of data0 */
    time_start = chrono::high_resolution_clock::now();
    
    vector<vector<vector<Ciphertext>>> variant_ct_0(n_individuals_0,
                                                    vector<vector<Ciphertext>> (n_alleles, vector<Ciphertext> (nctxts)));
    std::cout << "|   ct0.size = [" << n_individuals_0 << "][" << n_alleles << "][" << nctxts << "]" << endl;
    
    if(encoding == 0){
        heenc0.encrypt_data(variant_ct_0, geno0_500, n_posns, n_alleles_trials, nctxts_trial);
    } else if (encoding == 1){
        heenc0.encrypt_data(variant_ct_0, geno1_500, n_posns, n_alleles_trials, nctxts_trial);
    }
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Encryption-0: " << time_diff.count()/(1000000.0) << " sec] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

    /* Encryption of data1 */
    time_start = chrono::high_resolution_clock::now();

    vector<vector<vector<Ciphertext>>> variant_ct_1(n_individuals_1, vector<vector<Ciphertext>> (n_alleles, vector<Ciphertext> (nctxts)));
    std::cout << "|   ct1.size = [" << n_individuals_1 << "][" << n_alleles << "][" << nctxts << "]" << endl;

    if(encoding == 0){
        heenc1.encrypt_data(variant_ct_1, geno0_1000, n_posns, n_alleles_trials, nctxts_trial);
    } else if (encoding == 1){
        heenc1.encrypt_data(variant_ct_1, geno1_1000, n_posns, n_alleles_trials, nctxts_trial);
    }

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "Encryption-1: " << time_diff.count()/(1000000.0) << " sec] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

    cout << "Modulus chain index for encryption: q[" << context.get_context_data(variant_ct_0[0][0][0].parms_id()) ->chain_index() << "]" << endl;
    cout << "Noise budget in fresh encryption (c0): " << decryptor_0.invariant_noise_budget(variant_ct_0[0][0][0]) << " bits" << endl;
    cout << "Noise budget in fresh encryption (c1): " << decryptor_1.invariant_noise_budget(variant_ct_1[0][0][0]) << " bits" << endl;

    cout << "+-------------------------------------------+" << endl;
    cout << "|                 Evaluation                |" << endl;
    cout << "+-------------------------------------------+" << endl;

    HEeval heeval(batch_encoder, decryptor, evaluator);
    int n_threads = Thread::availableThreads(); // omp_get_max_threads();

    cout << "Step 1: re-encryption" << endl;

    time_start = chrono::high_resolution_clock::now();
    
    for(int allele_i = 0; allele_i < n_alleles_trials; allele_i++){
        cout << allele_i << ":" ;
        for(int k = 0; k < nctxts_trial; ++k){
            //if(k%10==0){cout << k << ",";}
            MT_EXEC_RANGE(n_individuals_0, first, last);
            for(int i = first; i < last; i++){
                evaluator.switch_keys_inplace(variant_ct_0[i][allele_i][k], switch_keys_0);
            }
            MT_EXEC_RANGE_END
        }
        cout << endl;
    }

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "KeySwitch-0: " << time_diff.count()/(1000000.0) << " sec] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

    time_start = chrono::high_resolution_clock::now();
    
    for(int allele_i = 0; allele_i < n_alleles_trials; allele_i++){
        cout << allele_i << ":" ;
        for(int k = 0; k < nctxts_trial; ++k){
            //if(k%10==0){cout << k << ",";}
            MT_EXEC_RANGE(n_individuals_1, first, last);
            for(int i = first; i < last; i++){
                evaluator.switch_keys_inplace(variant_ct_1[i][allele_i][k], switch_keys_1);
            }
            MT_EXEC_RANGE_END
        }
        cout << endl;
    }

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    getrusage(RUSAGE_SELF, &usage);
    cout << "KeySwitch-1: " << time_diff.count()/(1000000.0) << " sec] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

    cout << "Modulus chain index for encryption (c0): q[" << context.get_context_data(variant_ct_0[0][0][0].parms_id()) ->chain_index() << "]" << endl;
    cout << "Noise budget in fresh encryption (c0): " << decryptor.invariant_noise_budget(variant_ct_0[0][0][0]) << " bits" << endl;
    cout << "Modulus chain index for encryption (c1): q[" << context.get_context_data(variant_ct_1[0][0][0].parms_id()) ->chain_index() << "]" << endl;
    cout << "Noise budget in fresh encryption (c1): " << decryptor.invariant_noise_budget(variant_ct_1[0][0][0]) << " bits" << endl;
    cout << endl;

    // Step 2: Aggregate
    cout << "Step 2: aggregation" << endl;
    time_start = chrono::high_resolution_clock::now();

    vector<vector<Ciphertext>> res_ct(n_alleles_trials, vector<Ciphertext> (nctxts));

    for(int allele_i = 0; allele_i < n_alleles_trials; allele_i++){
        MT_EXEC_RANGE(nctxts_trial, first, last);
        for(int k = first; k < last; k++){
            res_ct[allele_i][k] = variant_ct_0[0][allele_i][k];
            for(int i = 1; i < n_individuals_0; ++i){
                evaluator.add_inplace(res_ct[allele_i][k], variant_ct_0[i][allele_i][k]);
            }
        }
        MT_EXEC_RANGE_END
        
        MT_EXEC_RANGE(nctxts_trial, first, last);
        for(int k = first; k < last; k++){
            for(int i = 1; i < n_individuals_1; ++i){
                evaluator.add_inplace(variant_ct_1[0][allele_i][k], variant_ct_1[i][allele_i][k]);
            }
        }
        MT_EXEC_RANGE_END

        MT_EXEC_RANGE(nctxts_trial, first, last);
        for(int k = first; k < last; k++){
            evaluator.add_inplace(res_ct[allele_i][k], variant_ct_1[0][allele_i][k]);
        }
        MT_EXEC_RANGE_END
    }

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Aggregation (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "Modulus chain index for encryption: q[" << context.get_context_data(res_ct[0][0].parms_id()) ->chain_index() << "]" << endl;
    cout << "Noise budget in result: " << decryptor.invariant_noise_budget(res_ct[0][0]) << " bits" << endl;

    cout << "+-------------------------------------------+" << endl;
    cout << "|                Decryption                 |" << endl;
    cout << "+-------------------------------------------+" << endl;

    time_start = chrono::high_resolution_clock::now();

    vector<vector<int>> res (n_alleles, vector<int>(n_posns, 0ULL));
    heeval.decrypt_result(res, res_ct, n_posns, n_alleles_trials, nctxts_trial);

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Decryption (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

#if (DEBUG)
    cout << "+-------------------------------------------+" << endl;
    cout << "|               Correctness                 |" << endl;
    cout << "+-------------------------------------------+" << endl;

    /* Check the correctness */
    time_start = chrono::high_resolution_clock::now();
    string filename = "data/Aggregation_Matrix_Data/task3_res_plain_" + to_string(encoding) + ".txt";
    ifstream openFile(filename.data());
    vector<vector<int>> plain_geno_res;

    if(openFile.is_open()) {
        read_data(plain_geno_res, filename, n_alleles_trials);

        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        getrusage(RUSAGE_SELF, &usage);
        cout << "Read the precomputed results (seconds) : " << time_diff.count()/(1000000.0) << " s] w/ " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    }
    else{ // read the whole genotypes and precompute their summation, and store them
        vector<vector<int>> prec_geno_res_0;
        string filename0 = "data/Aggregation_Matrix_Data/task2_res_plain_"  + to_string(n_individuals_0) + "_" + to_string(encoding) + "_" + to_string(n_alleles);
        read_data(prec_geno_res_0, filename0, n_alleles);
        cout << prec_geno_res_0[0].size() << endl;

        vector<vector<int>> prec_geno_res_1;
        filename0 = "data/Aggregation_Matrix_Data/task2_res_plain_"  + to_string(encoding) + "_" + to_string(n_individuals_1) + ".txt";
        read_data(prec_geno_res_1, filename0, n_alleles);

        for(int i = 0; i < prec_geno_res_0.size(); ++i){
            vector<int> temp;
            for(int j = 0; j < prec_geno_res_0[0].size(); ++j){
                temp.push_back(prec_geno_res_0[i][j] + prec_geno_res_1[i][j]);
            }
            plain_geno_res.push_back(temp);
        }
        write_data(plain_geno_res, filename);
        cout << "Precompute the results and store them" << endl;
    }

    int nfail = 0;
    for(int allele_i = 0; allele_i < n_alleles_trials; ++allele_i){
        int jend = ((nctxts > nctxts_trial)? nctxts_trial * slot_count: n_posns);
        for(int j = 0; j < jend; ++j){
            if(plain_geno_res[j][allele_i] != res[allele_i][j]){
                cout << j << "-pos: incorrect (" << plain_geno_res[j][allele_i] << ", " << res[allele_i][j] << ")" << endl;
                nfail++;
            }
        }
    }

    if(nfail != 0){
        cout << "> " << nfail << " fails" << endl;
    } else{
        cout << "> All Passed" << endl;
    }
#endif
}



