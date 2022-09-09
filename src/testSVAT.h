/*
 * @file       testSVAT.h, header file
 * @brief      function for running secure outsourcing of targeted variant annotation and genotype aggregation
 *
 * @author     Miran Kim
 * @copyright  GNU Pub License
 */

#ifndef __TEST_HESVAT__
#define __TEST_HESVAT__

// Task-1
void run_secure_variant_annotation(vector<uint64_t> &res, const char *variant_fp, const char *impact_fp, int n_posns, long nvariant);
 
void run_secure_variant_annotation_large(vector<vector<uint64_t>> &res,
                                         string dataset_str, vector<vector<string>> chr_info, string target);

// Task-2
void run_secure_variant_aggregation(int n_posns, int encoding, int n_individuals, const int n_alleles_trials = 5);
 
// Task-3
void run_secure_variant_reencrypt_aggregation(int n_posns, int encoding, const int n_alleles_trials = 5);

#endif //
