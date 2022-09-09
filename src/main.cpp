/*
 * @file       main.cpp, cpp file
 * @brief      test program
 *
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

#include "utils.h"
#include "thread.h"
#include "testPVAT.h"
#include "testSVAT.h"

using namespace std;

// e.g., ./test task1 24 10000 (./test task1 24 20000 or 50000)
// e.g., ./test task1large 24 DEL
// e.g., ./test task2 24 0 500 1
// e.g., ./test task3 24 0 1

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
        
        vector<uint64_t> res;
        run_secure_variant_annotation(res, variant_fp, impact_fp, n_posns, nvariant);
        
        
        cout << "+----------------------------------------+" << endl;
        cout << "|              Correctness               |" << endl;
        cout << "+----------------------------------------+" << endl;

        string filename = dataset_str + "res_plain_" + to_string(nvariant) + ".txt";
        vector<uint64_t> plain_res;
        ifstream openFile(filename.data());
        if(openFile.is_open()) {
            read_data(plain_res, filename);

        } else{
            run_variant_annotation(variant_fp, impact_fp, n_posns, nvariant);
            read_data(plain_res, filename);
        }

        int nfail = 0;
        for(long k = 0; k < n_posns; ++k){
            if(res[k] != plain_res[k])
            {
                //cout <<  k << ":" << res[k] << " =? " << plain_res[k] << " (" << variant_vector[k] << " * " <<  impact_vector[k] << ")"  << endl;
                nfail++;
            }
        }

        if(nfail != 0){
            cout << "> " << nfail << " fails" << endl;
        } else{
            cout << "> All Passed" << endl;
        }
    }
    else if(task == "task1large"){
        string target = argv[3];                // DEL, INS, SNV

        string dataset_str = "data/REVISION_DATA_VECTORS_14_52_01_07_2022/" + target + "/";
        string length_folder = dataset_str + "chr_lengths.list";     // the folder of chromosome lengths
        vector<vector<string>> chr_info;
        read_data(chr_info, length_folder);
        
        vector<vector<uint64_t>> res;
        run_secure_variant_annotation_large(res, dataset_str, chr_info, target);
        
        cout << "+----------------------------------------+" << endl;
        cout << "|              Correctness               |" << endl;
        cout << "+----------------------------------------+" << endl;
#if 0
        int nchr = res.size(); // res.size(); // number of chromosome int nchr = chr_info.size();
        
        for(int ch = 0; ch < nchr; ch++){
            string chr_id = chr_info[ch][0];                     // chromosome id
            int n_posns = (long) atoi(chr_info[ch][1].c_str());  // number of variants
            
            string variant_str = dataset_str + "VARIANT_SIGNAL/" + chr_id + ".bin"; // file for the annotation vector
            string impact_str = dataset_str + "ANNOTATION_SIGNAL/" + chr_id + ".bin";   // file for the impact vector
            const char *variant_fp = variant_str.c_str();
            const char *impact_fp = impact_str.c_str();

            string filename = dataset_str + "res_plain_" + chr_id + ".txt";
            cout << filename << ": ";

            vector<uint64_t> plain_res;
            ifstream openFile(filename.data());
            if(openFile.is_open()) {
                read_data(plain_res, filename);
                cout << " (read): " ;

            } else{
                if((target == "DEL")||(target == "INS")){
                    run_variant_annotation_large(variant_fp, impact_fp, n_posns, filename);
                } else if (target == "SNV"){
                    run_variant_annotation_large(variant_fp, impact_fp, 4 * n_posns, filename);
                }
                read_data(plain_res, filename);
                cout << " (compute & read): " ;
            }

            int nfail = 0;
            for(long k = 0; k < plain_res.size(); ++k){
                if(res[ch][k] !=  plain_res[k])
                {
                    //cout <<  k << ":" << res[k] << " =? " << plain_res[k] << " (" << variant_vector[ch][k] << " * " <<  impact_vector[ch][k] << ")"  << endl;
                    nfail++;
                }
            }
            
            if(nfail != 0){
                cout << "> " << nfail << " fails" << endl;
            } else{
                cout << "> All Passed" << endl;
            }
        }
#endif
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
