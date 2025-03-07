/*
MIT License

Copyright (c) 2024 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Disclaimer: We acknowledge that the TWINE block cipher doesn't adhere to statistical assumptions 
in differential analysis, such as the random sub-key assumption
or Markov cipher assumption. The tool's primary function is to find some bounds
for the security of TWINE against differential and differential-linear cryptanalysis.
*/


#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h> 
#include <math.h>
#include <omp.h>
#include <time.h> 
#include <sys/random.h>
#include "simeck48.h"


#define STEP ((1 << 10) - 1)
#define MAX_STRING_SIZE 64


bool dot_product(const uint32_t *x, const uint32_t *y);
unsigned int init_prng(unsigned int offset);
uint32_t generate_random_24bit();
void binaryToHex(const char *binary, char *hex, int hex_size);
void splitAndConvert(const char *binaryString, uint32_t output[]);
uint64_t dldistinguisher(int R, uint64_t N3, uint32_t *dp, uint32_t *lc);
double run_bunch_of_dldistinguishers(int R, int N1, uint64_t N2, uint64_t N3, uint32_t *dp, uint32_t *lc);

// #######################################################################################################
// #######################################################################################################
// ############################## User must change only the following lines ##############################
const int NUMBER_OF_THREADS = 1; // Number of threads
const int DEG1 = 1;              // Number of bunches per thread: N2 = 2^(DEG1)
const int DEG2 = 24;             // Number of queries per bunch:  N3 = 2^(DEG2)
int NUMBER_OF_EXPERIMENTS = 4;   // Number of independent experiments

int NUMBER_OF_ROUNDS = 11;       // Number of rounds
const char *DP_STR = "000000000000000000010000000000000000000000000000";
const char *LC_STR = "000000000000000000010100000000000000000000001000";
// uint32_t dp[] = {
//     0x20646e,
//     0x726963,
// };
// uint32_t lc[] = {
//     0x000000,
//     0x000000,
// };
// #######################################################################################################
// #######################################################################################################
