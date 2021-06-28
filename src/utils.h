/*
 * @file       utils.h, header file
 * @brief      function for reading and storing data
 *
 * @author     Miran Kim, Arif Harmanci
 * @copyright  GNU Pub License
 */


#ifndef __UTILS__
#define __UTILS__

#include <iostream>
#include <vector>
#include <string>

using namespace std;


void read_genotype_matrix(vector<vector<vector<bool>>>& res, int n_posns, int encoding, int n_individuals, int n_alleles);
void read_genotype_matrix(vector<vector<vector<vector<bool>>>>& res, int n_posns, int encoding, int n_individuals, int n_alleles);


template <typename T>
inline void colSums(vector<vector<int>>& res, vector<vector<vector<T>>> geno, int n_individuals, int n_alleles, int n_posns)
{
    res.resize(n_posns, vector<int> (n_alleles, 0ULL));
    
    for(int j = 0; j < n_posns; ++j){
        int* temp = new int[n_alleles];

        for(int allele_i = 0; allele_i < n_alleles; allele_i++){ // n_alleles
            temp[allele_i] = 0;
            for(int i = 0; i < n_individuals; ++i){
                temp[allele_i] += int(geno[i][allele_i][j]); // summation over n
            }
            if(temp[allele_i] != 0){
                res[j][allele_i] = temp[allele_i];
                cout << "(allele,j)=(" <<  allele_i << "," << j << "): " << temp[allele_i] << endl;
            }
        }
    }
}

inline void colSums(vector<vector<int>>& res, vector<vector<vector<vector<bool>>>> geno, int n_individuals, int n_alleles, int n_posns){
    res.resize(n_posns, vector<int> (n_alleles, 0ULL));
   
    for(int j = 0; j < n_posns; ++j){
        int* temp = new int[n_alleles];

        for(int allele_i = 0; allele_i < n_alleles; allele_i++){ // n_alleles
            temp[allele_i] = 0;
            for(int i = 0; i < n_individuals; ++i){
                temp[allele_i] += ((geno[1][i][allele_i][j] << 1) | (geno[0][i][allele_i][j])); // summation over n
            }
            if(temp[allele_i] != 0){
                res[j][allele_i] = temp[allele_i];
                cout << "(allele,j)=(" <<  allele_i << "," << j << "): " << temp[allele_i] << endl;
            }
        }
    }
}

void write_data(vector<vector<int>> data, string filename);
void read_data(vector<vector<int>>& res, string filename);
void read_data(vector<vector<int>>& res, string filename, int ncols);

#endif // __ALN_UTILS_PLAIN__
