#include <stdio.h>
#include <gmp#include <stdio.h>
#include <gmp.h>
#include "Rank_CBC_params.h"

/*
  * the variable nAw denotes n+w
  * the variable nSw denotes n-w
  * the variable nMw denotes n*w
  * the variable nDw denotes n/w
  * the variable n_bc_k denotes binomial coefficient nCk 
*/


int Rank_SDP_instance_list[6][4];
int MAX_DEPTH_log;
mpz_t MAX_DEPTH;
int result_list[7][7];

mpz_t m, n, k, w, r;
mpz_t grover_iteration;
mpz_t grover_processor_num;

void init_Rank_CBC_params(){
    // ROLLO 128 bit
    Rank_SDP_instance_list[0][0] = 67;
    Rank_SDP_instance_list[0][1] = 166;
    Rank_SDP_instance_list[0][2] = 83;
    Rank_SDP_instance_list[0][3] = 7;

    // ROLLO 192 bit
    Rank_SDP_instance_list[1][0] = 79;
    Rank_SDP_instance_list[1][1] = 194;
    Rank_SDP_instance_list[1][2] = 97;
    Rank_SDP_instance_list[1][3] = 8;

    // ROLLO 256 bit
    Rank_SDP_instance_list[2][0] = 97;
    Rank_SDP_instance_list[2][1] = 226;
    Rank_SDP_instance_list[2][2] = 113;
    Rank_SDP_instance_list[2][3] = 9;

    // RQC 128 bit
    Rank_SDP_instance_list[3][0] = 127;
    Rank_SDP_instance_list[3][1] = 113;
    Rank_SDP_instance_list[3][2] = 3;
    Rank_SDP_instance_list[3][3] = 7;

    // ROLLO 192 bit
    Rank_SDP_instance_list[4][0] = 151;
    Rank_SDP_instance_list[4][1] = 149;
    Rank_SDP_instance_list[4][2] = 5;
    Rank_SDP_instance_list[4][3] = 8;

    // ROLLO 256 bit
    Rank_SDP_instance_list[5][0] = 181;
    Rank_SDP_instance_list[5][1] = 179;
    Rank_SDP_instance_list[5][2] = 3;
    Rank_SDP_instance_list[5][3] = 9;

    MAX_DEPTH_log = 95;
    mpz_init_set_ui(MAX_DEPTH, 1);
    for(int index = 1; index <= MAX_DEPTH_log; index++){
        mpz_mul_ui(MAX_DEPTH, MAX_DEPTH, 2);
    }
}