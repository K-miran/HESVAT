/*
 * @file       utils.cpp, cpp file
 * @brief      function for reading and storing data
 *
 * @author     Miran Kim, Arif Harmanci
 * @copyright  GNU Pub License
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream
#include <bitset>
#include <math.h> 

#include "utils.h"
#include "thread.h"

using namespace std;


/*
 Read genotype data consisting of 0 or 1.
 
 @param[in] n_posns, the number of variants
 @param[in] encoding, genotype encoding method
 "0" indicates the "variant existence" is stored in the matrix as 0/1 values.  This file contains genotypes encoded with 1 bit per individual.
 "1" indicates "actual genotype" is stored as 0/1/2 values. This file is encoded with 2 bits per individual.
 @param[in] n_individuals, the number of individuals
 @param[in] n_alleles, the number of allels
 @param[out] res, [n_individuals][n_alleles][n_posns]
*/

void read_genotype_matrix(vector<vector<vector<bool>>>& res, int n_posns, int encoding, int n_individuals, int n_alleles)
{
    if(encoding != 0){
        throw invalid_argument("encoding is not valid");
    }
    
    int n_bits_per_sample = (encoding == 0? 1: 2);
    int n_bytes_per_posn = floor((n_bits_per_sample * n_individuals) / 8) + 1;
    std::cout << "|   n_bits_per_sample = " << n_bits_per_sample<< ", nbytes per positions = " << n_bytes_per_posn << endl;

    // Allocate the genotype matrix
    res.resize(n_individuals, vector<vector<bool>> (n_alleles, vector<bool> (n_posns, 0ULL)));
    cout << "allocate the genotype matrix" << endl;
    
    /* Load variant data */
    string file_str = "data/Aggregation_Matrix_Data/genotype_matrix_1000_" + to_string(encoding) + "_500000_"  + to_string(n_individuals) + ".bin";
    const char *file_fp = file_str.c_str();

    // Allocate the matrix
    // per_posn_per_allele_genotypes = [n_posns][5][n_bytes_per_posn + 1]
    unsigned char ***per_posn_per_allele_genotypes = new unsigned char**[n_posns + 1];
    for(int posn = 0; posn < n_posns + 1; posn++){
        per_posn_per_allele_genotypes[posn] =  new unsigned char*[5];
        for(int allele_i = 0; allele_i < n_alleles; allele_i++){
            per_posn_per_allele_genotypes[posn][allele_i] = new unsigned char[n_bytes_per_posn + 1];
            memset(per_posn_per_allele_genotypes[posn][allele_i], 0, sizeof(unsigned char) * n_bytes_per_posn); // n bytes
        }
    }
    
    FILE *f_geno_matrix;
    f_geno_matrix = fopen(file_fp,"rb");  // r for read, b for binary
  
   
    while(1)
    {
        // first read the first int-type number which is used for position (vectior_pos)
        int vector_pos = 0;
        int n_read = fread(&vector_pos, sizeof(int), 1, f_geno_matrix);   // 4 bytes (32-bits)
        
        if(n_read == 0){
            fprintf(stderr, "Finished reading.\n");
            break;
        }

        // Then read 5 genotype vectors (for each of A,C,G,T,N) into arrays of length n_bytes_per_position
        // unsigned char** per_allele_genotypes = new unsigned char*[5];
        for(int allele_i = 0; allele_i < 5; allele_i++){
            if(fread(per_posn_per_allele_genotypes[vector_pos][allele_i], sizeof(unsigned char), n_bytes_per_posn, f_geno_matrix) != n_bytes_per_posn){
                fprintf(stderr, "Sanity check failed: Could not read the genotypes!\n");
                exit(0);
            }
        }

        // Calcuate the genotype of individuals for each variant and each allele
        //cout << "================= (" << vector_pos << ") ==================" << endl;
        for(int allele_i = 0; allele_i < n_alleles; allele_i++){
            for(int n = 0; n < n_individuals; ++n){
                int byte_i = (n) / 8; // n_bits_per_sample = 1
                int bit_i =  (n) % 8;
                unsigned char buffer = per_posn_per_allele_genotypes[vector_pos][allele_i][byte_i]; // 8 bits
                
                if(buffer != 0x00){
                    res[n][allele_i][vector_pos] = (buffer & (1 << bit_i)) >> bit_i;
    
                    if(res[n][allele_i][vector_pos]  > 1){
                        cout << "error:" << vector_pos <<endl;
                    }
                }
            }
        }
    }
    fclose(f_geno_matrix);
        
    delete [] per_posn_per_allele_genotypes;
}

/*
 Read genotype data consisting of 0,1,2.
 
 @param[in] n_posns, the number of variants
 @param[in] encoding, genotype encoding method
 "0" indicates the "variant existence" is stored in the matrix as 0/1 values.  This file contains genotypes encoded with 1 bit per individual.
 "1" indicates "actual genotype" is stored as 0/1/2 values. This file is encoded with 2 bits per individual.
 @param[in] n_individuals, the number of individuals
 @param[in] n_alleles, the number of allels
 @param[out] res, [n_individuals][n_alleles][n_posns]
*/

void read_genotype_matrix(vector<vector<vector<vector<bool>>>>& res, int n_posns, int encoding, int n_individuals, int n_alleles)
{
    if(encoding != 1){
        throw invalid_argument("encoding is not valide");
    }
    
    int n_bits_per_sample = (encoding == 0? 1: 2);
    int n_bytes_per_posn = floor((n_bits_per_sample * n_individuals) / 8) + 1;
    cout << "n_bits_per_sample = " << n_bits_per_sample << ", nbytes per positions = " << n_bytes_per_posn << endl;

    // Allocate the genotype matrix
    res.resize(2, vector<vector<vector<bool>>> (n_individuals, vector<vector<bool>> (n_alleles, vector<bool> (n_posns, 0ULL))));
    cout << "allocate the genotype matrix" << endl;
    
    /* Load variant data */
    string file_str = "data/Aggregation_Matrix_Data/genotype_matrix_1000_" + to_string(encoding) + "_500000_"  + to_string(n_individuals) + ".bin";
    const char *file_fp = file_str.c_str();

    // Allocate the matrix
    // per_posn_per_allele_genotypes = [n_posns][5][n_bytes_per_posn + 1]
    unsigned char ***per_posn_per_allele_genotypes = new unsigned char**[n_posns + 1];
    for(int posn = 0; posn < n_posns + 1; posn++){
        per_posn_per_allele_genotypes[posn] =  new unsigned char*[5];
        for(int allele_i = 0; allele_i < n_alleles; allele_i++){
            per_posn_per_allele_genotypes[posn][allele_i] = new unsigned char[n_bytes_per_posn + 1];
            memset(per_posn_per_allele_genotypes[posn][allele_i], 0, sizeof(unsigned char) * n_bytes_per_posn); // n bytes
        }
    }
    
    FILE *f_geno_matrix;
    f_geno_matrix = fopen(file_fp,"rb");  // r for read, b for binary
  
    while(1)
    {
        // first read the first int-type number which is used for position (vectior_pos)
        int vector_pos = 0;
        int n_read = fread(&vector_pos, sizeof(int), 1, f_geno_matrix);   // 4 bytes (32-bits)
        
        if(n_read == 0){
            fprintf(stderr, "Finished reading.\n");
            break;
        }

        // Then read 5 genotype vectors (for each of A,C,G,T,N) into arrays of length n_bytes_per_position
        // unsigned char** per_allele_genotypes = new unsigned char*[5];
        for(int allele_i = 0; allele_i < 5; allele_i++){
            if(fread(per_posn_per_allele_genotypes[vector_pos][allele_i], sizeof(unsigned char), n_bytes_per_posn, f_geno_matrix) != n_bytes_per_posn){
                fprintf(stderr, "Sanity check failed: Could not read the genotypes!\n");
                exit(0);
            }
        }
    
        //cout << "================= (" << vector_pos << ") ==================" << endl;
        for(int allele_i = 0; allele_i < n_alleles; allele_i++){
            for(int n = 0; n < n_individuals; ++n){
                int byte_i = (n << 1) / 8; // n * n_bits_per_sample = 2n
                int bit_i =  (n << 1) % 8;
                unsigned char buffer = per_posn_per_allele_genotypes[vector_pos][allele_i][byte_i]; // 8 bits
                
                if(buffer != 0x00){
                    res[0][n][allele_i][vector_pos] = (buffer & (1 << bit_i)) >> bit_i;
                    res[1][n][allele_i][vector_pos] = (buffer & (2 << bit_i)) >> (bit_i + 1);
                  
                    if((res[0][n][allele_i][vector_pos] > 1)||(res[1][n][allele_i][vector_pos] > 1)){
                        cout << "error:" << vector_pos <<endl;
                    }
                }
            }
        }
    }
    fclose(f_geno_matrix);
    
    delete [] per_posn_per_allele_genotypes;
}



/*
 @param[in] data, data matrix
 @param[in] filename, path to be written
*/

void write_data(vector<vector<int>> data, string filename)
{
    ofstream outf_tmp(filename);
    outf_tmp.close();
    
    fstream outf;
    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);   // open the file
    
    size_t nrows = data.size();
    size_t ncols = data[0].size();
    //cout << nrows << "," << ncols << endl;
    
    for(size_t i = 0; i < nrows; ++i){
        for(size_t j = 0; j < ncols; ++j){
            outf << data[i][j] << "\t";
        }
        outf << endl;
    }
    outf.close();
}

void write_data(uint64_t *data, size_t length, string filename)
{
    ofstream outf_tmp(filename);
    outf_tmp.close();
    
    fstream outf;
    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);   // open the file
    
    for(size_t i = 0; i < length; ++i){
        outf << data[i] << endl;
    }
    outf.close();
}

void write_data(vector<uint64_t> data, string filename)
{
    ofstream outf_tmp(filename);
    outf_tmp.close();
    
    fstream outf;
    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);   // open the file
    
    for(size_t i = 0; i < data.size(); ++i){
        outf << data[i] << endl;
    }
    outf.close();
}


void read_data(vector<uint64_t>& res, string filename)
{
    ifstream openFile(filename.data());
    
    string line;
    char split_char = 0x0A; // new line
    if(openFile.is_open()){
        while(getline(openFile, line, split_char)){
            uint64_t num = stol(line);
            res.push_back(num);
        }
    }
}


/*
 @param[in] filename, path of an input file
 @param[out] res, data matrix
 @throws std::invalid_argument if the data is not readable
*/

void read_data(vector<vector<int>>& res, string filename)
{
    char split_char = 0x09; // "\t";
    ifstream openFile(filename.data());

    vector<vector<string>> res_str;
    if(openFile.is_open()) {
        string line;
    
        /*  read each line */
        while(getline(openFile, line)){
            vector<string> vecsline;
            istringstream split(line);
            vector<string> tokens;
            for (string each; getline(split, each, split_char); tokens.push_back(each));
            
            res_str.push_back(tokens);
        }
    } else {
        throw invalid_argument("Error: cannot read file");
    }
    
    for(size_t i = 0; i < res_str.size(); ++i){
        vector<int> temp;
        for(size_t j = 0; j < res_str[i].size(); ++j){
            temp.push_back(atoi(res_str[i][j].c_str()));
        }
        res.push_back(temp);
    }
}


void read_data(vector<vector<string>>& res, string filename)
{
    char split_char = 0x09; // "\t";
    ifstream openFile(filename.data());

    if(openFile.is_open()) {
        string line;
        
        /*  read each line */
        while(getline(openFile, line)){
            vector<string> vecsline;
            istringstream split(line);
            vector<string> tokens;
            for (string each; getline(split, each, split_char); tokens.push_back(each));
            
            res.push_back(tokens);
        }
    } else {
        throw invalid_argument("Error: cannot read file");
    }
}

/*
 @param[in] filename, path of an input file
 @param[in] ncols, the number of the extracted columns
 @param[out] res, data matrix
 @throws std::invalid_argument if the data is not readable
*/

void read_data(vector<vector<int>>& res, string filename, int ncols)
{
    char split_char = 0x09; // "\t";
    ifstream openFile(filename.data());
    
    vector<vector<string>> res_str;
    if(openFile.is_open()) {
        string line;
        
        /*  read each line */
        while(getline(openFile, line)){
            vector<string> vecsline;
            istringstream split(line);
            vector<string> tokens;
            for (string each; getline(split, each, split_char); tokens.push_back(each));
            
            res_str.push_back(tokens);
        }
    } else {
        throw invalid_argument("Error: cannot read file");
    }
    
    for(size_t i = 0; i < res_str.size(); ++i){
        vector<int> temp;
        for(size_t j = 0; j < ncols; ++j){
            temp.push_back(atoi(res_str[i][j].c_str()));
        }
        res.push_back(temp);
    }
}


