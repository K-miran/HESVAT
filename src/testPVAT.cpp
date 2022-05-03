/*
 * @file       testPVAT.cpp, cpp file
 * @brief      function for running targeted variant annotation and genotype aggregation
 *
 * @copyright  MIT License
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
#include "testPVAT.h"


using namespace std;

/*
 @Task1: Perform the variant annotation
 @param[in] variant_fp, filename of an annotation vector
 @param[in] impact_fp, filename of a variant loci vector
 @param[in] int n_posns, the number of positions
*/
void run_variant_annotation(const char *variant_fp, const char *impact_fp, int n_posns, long nvariant)
{
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
   
    vector<uint64_t> res;
    
    for(int k = 0; k < n_posns; ++k){    // npo
        res.push_back(variant_vector[k] * impact_vector[k]);
    }
    
    string filename = "data/Annotation_Vector_Data/res_plain_" + to_string(nvariant) + ".txt";
    write_data(res, filename);
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

void run_variant_aggregation(int n_posns, int encoding, int n_individuals, int n_alleles_trials)
{
    long memoryscale = (1 << 20); // 2^20: linux, 2^30: mac
    struct rusage usage;
    
    int n_alleles = 5;

    /* Load data */
    vector<vector<vector<bool>>> geno0;         // bool[n_individuals][n_allels][n_posns]
    vector<vector<vector<vector<bool>>>> geno1; // bool[2][n_individuals][n_allels][n_posns]
    
    if(encoding == 0){
        read_genotype_matrix(geno0, n_posns, encoding, n_individuals, n_alleles);
    } else if(encoding == 1){
        read_genotype_matrix(geno1, n_posns, encoding, n_individuals, n_alleles);
    }
   
    getrusage(RUSAGE_SELF, &usage);
    cout << "Load Data : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
   
    vector<vector<int>> res; // [n_allels][n_posns]
    string filename = "data/Aggregation_Matrix_Data/task2_res_plain_" + to_string(encoding) + "_" + to_string(n_individuals) + ".txt";
    if(encoding == 0){
        colSums(res, geno0, n_individuals, n_alleles_trials, n_posns);// n_posns
    } else if(encoding == 1){
        colSums(res, geno1, n_individuals, n_alleles_trials, n_posns);// n_posns
    }
    write_data(res, filename);
}


/*
 @Task3: Perform the secure variant aggregation with re-encryption

 @param[in] n_posns, the number of positions
 @param[in] encoding, genotype encoding method
 "0" indicates the "variant existence" is stored in the matrix as 0/1 values.  This file contains genotypes encoded with 1 bit per individual.
 "1" indicates "actual genotype" is stored as 0/1/2 values. This file is encoded with 2 bits per individual.
 @param[in] n_alleles_trials, the number of alleles
*/

void run_matrices_aggregation(int n_posns, int encoding, int n_alleles_trials)
{
    long memoryscale = (1 << 20); // 2^20: linux, 2^30: mac
    struct rusage usage;
    
    int n_alleles = 5;
    int n_individuals_0 = 500;
    int n_individuals_1 = 1000;
     
    /* Load data */
    vector<vector<vector<bool>>> geno0_500;             // 0-enc, 500 samples
    vector<vector<vector<bool>>> geno0_1000;            // 0-enc, 1000 samples
    
    vector<vector<vector<vector<bool>>>> geno1_500;     // 1-enc, 500 samples
    vector<vector<vector<vector<bool>>>> geno1_1000;    // 1-enc, 1000 samples
    
    if(encoding == 0){
        read_genotype_matrix(geno0_500, n_posns, encoding, n_individuals_0, n_alleles);
        read_genotype_matrix(geno0_1000, n_posns, encoding, n_individuals_1, n_alleles);
    } else if(encoding == 1){
        read_genotype_matrix(geno1_500, n_posns, encoding, n_individuals_0, n_alleles);
        read_genotype_matrix(geno1_1000, n_posns, encoding, n_individuals_1, n_alleles);
    }
   
    getrusage(RUSAGE_SELF, &usage);
    cout << "Load Data : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    /* Result of the 1st matrix */
    vector<vector<int>> prec_geno_res_0;
    string filename0 = "data/Aggregation_Matrix_Data/task2_res_plain_"  + to_string(encoding)  + "_" + to_string(n_individuals_0) + ".txt";
    
    ifstream openFile0(filename0.data());
    if(openFile0.is_open()) {
        read_data(prec_geno_res_0, filename0, n_alleles);
    } else{
        if(encoding == 0){
            colSums(prec_geno_res_0, geno0_500, n_individuals_0, n_alleles, n_posns);
        } else if(encoding == 1){
            colSums(prec_geno_res_0, geno1_500, n_individuals_0, n_alleles_trials, n_posns);
        }
        write_data(prec_geno_res_0, filename0);
    }
    
    /* Result of the 2nd matrix */
    vector<vector<int>> prec_geno_res_1;
    filename0 = "data/Aggregation_Matrix_Data/task2_res_plain_"  + to_string(encoding)  + "_" + to_string(n_individuals_1) + ".txt";
    
    ifstream openFile1(filename0.data());
    if(openFile1.is_open()) {
        read_data(prec_geno_res_1, filename0, n_alleles);
    } else{
        if(encoding == 0){
            colSums(prec_geno_res_1, geno0_1000, n_individuals_1, n_alleles, n_posns);
        } else if(encoding == 1){
            colSums(prec_geno_res_1, geno1_1000, n_individuals_1, n_alleles_trials, n_posns);
        }
        write_data(prec_geno_res_1, filename0);
    }

    /* Result */
    string filename = "data/Aggregation_Matrix_Data/task3_res_plain_" + to_string(encoding) + ".txt";
    ifstream openFile(filename.data());
    vector<vector<int>> res;
    
    for(int i = 0; i < prec_geno_res_0.size(); ++i){
        vector<int> temp;
        for(int j = 0; j < prec_geno_res_0[0].size(); ++j){
            temp.push_back(prec_geno_res_0[i][j] + prec_geno_res_1[i][j]);
        }
        res.push_back(temp);
    }
    write_data(res, filename);
}
