/*
 * @file       main.cpp, cpp file
 * @brief      test program
 *
 * @date       Apr. 20, 2021
 * @copyright  GNU Pub License
 */

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
#include <sys/resource.h>   // check the memory
#include <omp.h>            // openmp

#include "thread.h"
#include "testSVAT.h"

using namespace std;

// e.g., ./test task1 24 10000 (or 20000 or 50000)
// e.g., ./test task2 24 0 500 1
// e.g., ./test task3 24 0 5


int main(int argc, char *argv[])
{
    string task = argv[1];                  // task1, task2, task3
    size_t NUM_THREADS = atoi(argv[2]);     // number of threads
    
    Thread::initThreadPool(NUM_THREADS);
    
    if(task == "task1"){
        long nvariant = atoi(argv[3]);                      // number of variants
        string dataset_str = "data/Annotation_Vector_Data/";     // file directory
        
        string variant_str = dataset_str + to_string(nvariant) + "_variants_signal_1.bin";   // file for the annotation vector
        string impact_str = dataset_str + "impact_signal_1.bin";                        // file for the variant loci vector
        
        const char *variant_fp = variant_str.c_str();
        const char *impact_fp = impact_str.c_str();
        int n_posns = 3409574;  // number of positions
        
        cout << "+----------------------------------------+" << endl;
        cout << "|            Annotation Vector           |" << endl;
        cout << "|                 (" << nvariant  << ")                |" << endl;
        cout << "+----------------------------------------+" << endl;
        
        run_secure_variant_annotation(variant_fp, impact_fp, n_posns);
    }
    else if(task == "task2"){
        int n_posns = 6601984;                   // number of positions
        int encoding = atoi(argv[3]);           // genotype encoding method: 0 or 1
        int n_individuals = atoi(argv[4]);      // number of individuals: 500 or 1000
        
        if((encoding != 0) && (encoding != 1)){
            throw invalid_argument("genotype encoding is not supportive");
        }
        if((n_individuals != 500) && (n_individuals != 1000)){
            throw invalid_argument("number of individuals is not supportive");
        }
        
        cout << "+------------------------------------+" << endl;
        cout << "|         Aggregation Matrix         |" << endl;
        cout << "|      (enc,nindviduals)=(" << encoding  << "," << n_individuals << ")     |" << endl;
        cout << "+------------------------------------+" << endl;
        
        int n_alleles_trials = atoi(argv[5]);    // number of alleles: 1,2,3,4,5
        
        run_secure_variant_aggregation(n_posns, encoding, n_individuals, n_alleles_trials);
    }
    else if(task == "task3"){
        int n_posns = 6601984;                   // number of positions
        int encoding = atoi(argv[3]);           // genotype encoding method: 0 or 1
        
        if((encoding != 0) && (encoding != 1)){
            throw invalid_argument("genotype encoding is not supportive");
        }
        
        cout << "+-------------------------------------------+" << endl;
        cout << "|   Aggregation Matrix by Proxy Encryption  |" << endl;
        cout << "+-------------------------------------------+" << endl;
        
        int n_alleles_trials = atoi(argv[4]);    // number of alleles: 1,2,3,4,5
    
        run_secure_variant_reencrypt_aggregation(n_posns, encoding, n_alleles_trials);
    }

    return 0;
	
}
