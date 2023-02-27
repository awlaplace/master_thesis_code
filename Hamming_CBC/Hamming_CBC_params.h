#ifndef _HAMMING_CBC_PARAMS_H_
#define _HAMMING_CBC_PARAMS_H_

#include <stdio.h>
#include <gmp.h>

// SDP_instance
extern int SDP_instance_list[9][3];

// MAX_DEPTH
extern int MAX_DEPTH_log;
extern mpz_t MAX_DEPTH;

// n, k, w,  G-cost, D-cost, W-cost を格納する配列
extern int Grover_result_list[9][6];
extern int MMT_BJMM_result_list[9][11];

extern mpz_t n, k, w, p, ell, ell1, ell2, epsilon;
extern mpz_t kAell, nSkSell, wSp, kAellD2, pD4Aepsilon, kAellSp, c2epsilon;
extern mpz_t grover_iteration, QW_MMT_iteration;
extern mpz_t grover_processor_num, QW_MMT_processor_num;

void init_Hamming_CBC_params();

#endif