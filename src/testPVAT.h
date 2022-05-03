/*
 * @file       testVAT.h, header file
 * @brief      function for running targeted variant annotation and genotype aggregation
 *
 * @copyright  MIT License
 */

#ifndef __TEST_VAT__
#define __TEST_VAT__

// Task-1
void run_variant_annotation(const char *variant_fp, const char *impact_fp, int n_posns, long nvariant);
 
// Task-2
void run_variant_aggregation(int n_posns, int encoding, int n_individuals, const int n_alleles_trials = 5);
 
// Task-3
void run_matrices_aggregation(int n_posns, int encoding, const int n_alleles_trials = 5);

#endif 
