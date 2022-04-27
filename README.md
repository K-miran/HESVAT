# HESVAT: Secure Variant Annotation and Aggregation Tool

## Contents

* [Download](#download)
* [Installation](#installation)
    *[Dependencies](#Dependencies)
    *[To Installl HE-SVAT](#To-Install-HE-SVAT)
* [Description of HE-SVAT](#Description-of-HE-SVAT)
    * [Task 1: Secure Annotation](#Task-1-Secure-Annotation)
    * [Task 2: Secure Aggregation](#Task-2-Secure-Aggregation)    
    * [Task 3: Secure Aggregation by Proxy Encryption](#Task-3-Secure-Aggregation-by-proxy-encryption)      
* [Examples](#Examples)
    * [Task 1: Example Run](#Task-1-example-run)
    * [Task 2: Example Run](#Task-2-example-run)   
    * [Task 3: Example Run](#Task-3-example-run)      
    
## Download 

You can download the HESVAT code by clicking on the green `Clone or Download` button and downloading zip or by checking out via git. After that, navigate to the checkout (or download) directory. If you downloaded the zip file of the source, unzip it using "unzip master.zip". 


## Installation 

We recommend to install `HESVAT` into a C++ environment.

### Dependencies 
- homebrew 
- texinfo 
- m4 >=1.4.16
- CMake >= 3.12
- Compiler: g++ version >= 6.0


### To Install HE-SVAT

1. Download the library or clond the repository with the following command: 
    ```
    git clone https://github.com/K-miran/HESVAT.git
    ```
    We will call the directory that you cloned as ${HESVAT_ROOT}.
2. Install the modified SEAL library (SEAL 3.6.4) with the following commands: 
    ```
    cd ${HESVAT_ROOT}/src/external/SEAL-main
    cmake. 
    make 
    ```
3. Install the HESVAT library by running the following commands :
    ```
    cd ${HESVAT_ROOT}
    cmake . 
    make
    ```

## Description of HESVAT

### Task 1: Secure Annotation

**1. Description of dataset**

There are 3 files that contain variant vectors, which are sensitive and need to be encrypted:
- `10000_variants_signal_1.bin`: 10000 variants vectorized into the 3,409,574 positions. 
- `20000_variants_signal_1.bin`: 20000 variants vectorized into the 3,409,574 positions. 
- `50000_variants_signal_1.bin`: 50000 variants vectorized into the 3,409,574 positions. 

The annotation vectors contain the variant existence for 3,409,574 "vectorized" positions. They store 4 bytes for each variant and they are 0/1 values: "1" means there is a variant at the position, otherwise 0.

There is one more file that does *not* need to be encrypted,  `impact_signal_1.bin`. This variant loci vector contains 3,409,574 entries that hold the annotations of all mutations on the target regions. The impact value indicates the impact of the mutation at each position on the vector.

Your directory tree should look like this:
```
data
└── Annotation_Vector_Data
    ├── 10000_variants_signal_1.bin
    ├── 20000_variants_signal_1.bin
    ├── 50000_variants_signal_1.bin
    └── impact_signal_1.bin
```

**2. What do we want to do?**

We would like to perform the constant-ciphertext multiplication with the ciphertext of the annotation vector and the corresponding variant loci vector. 


### Task 2: Secure Aggregation

**1. Description of dataset**

Each genotype matrix file is named as:  `genotype_matrix_[Number of target regions]_[Genotype encoding]_[Number of variants]_[Number of individuals].bin`. The number of individuals is 500 or 1000. The total number of (possible) variant positions that we store in the matrix is n_posns = 6,601,984. 

There are 4 genotype files:
- genotype_matrix_1000_0_500000_500.bin: genotype data for 1000 target regions, with genotype encoding of 0, 500k variants, and 500 individuals. 
- genotype_matrix_1000_1_500000_500.bin: genotype data for 1000 target regions, with genotype encoding of 1, 500k variants, and 500 individuals. 
- genotype_matrix_1000_0_500000_1000.bin: genotype data for 1000 target regions, with genotype encoding of 0, 500k variants, and 1000 individuals. 
- genotype_matrix_1000_1_500000_1000.bin: genotype data for 1000 target regions, with genotype encoding of 1, 500k variants, and 1000 individuals. 

Each file generates 5 genotype matrices (for each of A, C, G, T, N, that is, one per each allele). For example, `genotype_matrix_1000_0_500000_500.bin` file generates 5 genotype matrices, each of which is of size 500*6601984. Here, the row and column dimensions indicate the number of individuals and the number of variant positions, respectively. 

To be specific, we have two types of encoding:
- In the 0-encoding, the "variant existence" is stored in the matrix as 0/1 values.  The value of 1 indicates the individual has the variant and 0 for not. This file contains genotypes encoded with 1 bit per individual. For 500 individuals, we use 500 bits, which can be stored in 63 bytes.
- In the 1-encoding,  the "actual genotype" is stored as 0/1/2 values. This file is encoded with 2 bits per individual. For 500 individuals, we use 1000 bits, but we use 126 bytes.


Your directory tree should look like this:
```
data
└── Aggregation_Matrix_Data
    ├── genotype_matrix_1000_0_500000_500.bin
    ├── genotype_matrix_1000_1_500000_500.bin
    ├── genotype_matrix_1000_0_500000_1000.bin
    └── genotype_matrix_1000_1_500000_1000.bin
```

**2. What do we want to do?**

We encrypt the full genotype matrices by taking the entries in row-major order and encrypting vectors as ciphertexts. Then we compute the frequencies of the mutations by securely aggregating over samples. 


### Task 3: Secure Aggregation by Proxy Encryption

**What do we want to do?**

- Start from two encrypted files (generated by the same encoding) where we assume they come from different sources, so they would be encrypted with different keys. 
- Perform re-encryption operations, so that input ciphertexts of genotype matrices are converted into new ciphertexts decryptable by the same secret key. 
- Do secure aggregation (summation over samples), so we get a vector of length 6,601,984 that can be decryptable by the secret key and each entry is the summation over all individuals.
- Decrypt the resulting vector.


## Examples

### Task 1: Example Run 
The following list of command-line arguments is required after the name of the test program (`foo`):
- Type of Task
- Number of threads
- Number of variants (e.g., 10000, 20000, 50000)

For instance, run the test program with different inputs:
```
./foo task1 4 10000
```

### Task 2: Example Run 
The following list of command-line arguments is required after the name of the test program (`foo`):
- Type of Task
- Number of threads
- Encoding type (e.g., 0, 1)
- Number of individuals (e.g., 500, 1000)
- Number of alleles (e.g., 1, 2, 3, 4, 5)

For instance, run the test program with different inputs:
```
./foo task2 24 0 500 1
```

### Task 3: Example Run 
The following list of command-line arguments is required after the name of the test program (`foo`):
- Type of Task
- Number of threads
- Encoding type (e.g., 0, 1)
- Number of alleles (e.g., 1, 2, 3, 4, 5)

For instance, run the test program with different inputs:
```
./foo task3 24 1 2
```


