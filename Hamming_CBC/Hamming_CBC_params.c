#include <stdio.h>
#include <gmp.h>
#include "Hamming_CBC_params.h"

/*
  * the variable nAw denotes n+w
  * the variable nSw denotes n-w
  * the variable nMw denotes n*w
  * the variable nDw denotes n/w
  * the variable n_bc_k denotes binomial coefficient nCk 
*/


int SDP_instance_list[9][3];
int MAX_DEPTH_log;
mpz_t MAX_DEPTH;
int Grover_result_list[9][6];
int MMT_BJMM_result_list[9][11];

mpz_t n, k, w, p, ell, ell1, ell2, epsilon;
mpz_t kAell, nSkSell, wSp, kAellD2, pD4Aepsilon, kAellSp, c2epsilon;
mpz_t grover_iteration, QW_MMT_iteration;
mpz_t grover_processor_num, QW_MMT_processor_num;


void init_Hamming_CBC_params(){
    // BIKE 128 bit
    SDP_instance_list[0][0] = 12323 * 2;
    SDP_instance_list[0][1] = 12323;
    SDP_instance_list[0][2] = 134;

    // BIKE 192 bit
    SDP_instance_list[1][0] = 24659 * 2;
    SDP_instance_list[1][1] = 24659;
    SDP_instance_list[1][2] = 199;

    // BIKE 256 bit
    SDP_instance_list[2][0] = 40973 * 2;
    SDP_instance_list[2][1] = 40973;
    SDP_instance_list[2][2] = 264;

    // Classic McEliece 128 bit
    SDP_instance_list[3][0] = 3488;
    SDP_instance_list[3][1] = 2720;
    SDP_instance_list[3][2] = 64;

    // Classic McEliece 192 bit
    SDP_instance_list[4][0] = 4608;
    SDP_instance_list[4][1] = 3360;
    SDP_instance_list[4][2] = 96;

    // Classic McEliece 256 bit
    SDP_instance_list[5][0] = 8192;
    SDP_instance_list[5][1] = 6528;
    SDP_instance_list[5][2] = 128;

    // HQC 128 bit
    SDP_instance_list[6][0] = 17669 * 2;
    SDP_instance_list[6][1] = 17669;
    SDP_instance_list[6][2] = 66 * 2;

    // HQC 192 bit
    SDP_instance_list[7][0] = 35851 * 2;
    SDP_instance_list[7][1] = 35851;
    SDP_instance_list[7][2] = 100 * 2;

    // HQC 256 bit
    SDP_instance_list[8][0] = 57637 * 2;
    SDP_instance_list[8][1] = 57637;
    SDP_instance_list[8][2] = 131 * 2;

    MAX_DEPTH_log = 95;
    mpz_init_set_ui(MAX_DEPTH, 1);
    for(int index = 1; index <= MAX_DEPTH_log; index++){
        mpz_mul_ui(MAX_DEPTH, MAX_DEPTH, 2);
    }

    for(int index = 0; index < 9; index++){
        Grover_result_list[index][0] = SDP_instance_list[index][0];
        Grover_result_list[index][1] = SDP_instance_list[index][1];
        Grover_result_list[index][2] = SDP_instance_list[index][2];
        Grover_result_list[index][3] = 0;
        Grover_result_list[index][4] = MAX_DEPTH_log;
        Grover_result_list[index][5] = 0;

        MMT_BJMM_result_list[index][0] = SDP_instance_list[index][0];
        MMT_BJMM_result_list[index][1] = SDP_instance_list[index][1];
        MMT_BJMM_result_list[index][2] = SDP_instance_list[index][2];
        MMT_BJMM_result_list[index][3] = 0;
        MMT_BJMM_result_list[index][4] = 0;
        MMT_BJMM_result_list[index][5] = 0;
        MMT_BJMM_result_list[index][6] = 0;
        MMT_BJMM_result_list[index][7] = 0;
        MMT_BJMM_result_list[index][8] = 0;
        MMT_BJMM_result_list[index][9] = MAX_DEPTH_log;
        MMT_BJMM_result_list[index][10] = 0;
    }

    mpz_init_set_ui(kAell, 0);
    mpz_init_set_ui(nSkSell, 0);
    mpz_init_set_ui(wSp, 0);
    mpz_init_set_ui(kAellD2, 0);
    mpz_init_set_ui(pD4Aepsilon, 0);
    mpz_init_set_ui(kAellSp, 0);
    mpz_init_set_ui(c2epsilon, 0);
}