/*
Copyright (C) 2023  
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>  // rand(), srand()
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
#include <time.h>    // time()
#include <stdbool.h>
#include <sys/random.h>

#define Nthreads 1
#define STEP ((1 << 9) - 1)

typedef unsigned long long int UINT64;

unsigned int init_prng(unsigned int offset);
int dot_product(const int mask[], const int data[]);
void print_state(int *m);
bool test();
UINT64 bunch_of_diff_lin_tests(int R, UINT64 N3, int* dp, int* lc);
UINT64 parallel_diff_lin_tests(int R, int N1, UINT64 N2, UINT64 N3, int *dp, int *lc);
void convert_hexstr_to_statearray(char hex_str[], int dx[32]);

// #######################################################################################################
// #######################################################################################################
// ############################## User must change only the following lines ##############################
const int DEG1 = 0;
const int DEG2 = 20;
const int NUMBER_OF_EXPERIMENTS = 10;   // Number of independent experiments
const int NUMBER_OF_ROUNDS = 11;       // Number of rounds

char DP_STR[] = "00000000000000a00000000000000000";
char DC_STR[] = "00000000000000020000000000000000";
// #######################################################################################################
// #######################################################################################################




// // #######################################################################################################
// // #######################################################################################################
// // ############################## User must change only the following lines ##############################
// const int DEG1 = 0;
// const int DEG2 = 20;
// const int NUMBER_OF_EXPERIMENTS = 3;   // Number of independent experiments
// const int NUMBER_OF_ROUNDS = 11;       // Number of rounds

// char DP_STR[] = "00000000000000020000000000000000";
// char DC_STR[] = "00000000000000020000000000000000";
// // Expected correlation: 1
// // ##################################################
