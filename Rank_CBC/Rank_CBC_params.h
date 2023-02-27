#ifndef _RANK_CBC_PARAMS_H_
#define _RANK_CBC_PARAMS_H_

#include <stdio.h>
#include <gmp.h>

// SDP_instance
extern int Rank_SDP_instance_list[6][4];

// MAX_DEPTH
extern int MAX_DEPTH_log;
extern mpz_t MAX_DEPTH;

extern mpz_t m, n, k, w, r;
extern mpz_t grover_iteration;
extern mpz_t grover_processor_num;
extern mpz_t CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Dcost, CGRSGloop_n_large_ancila;
extern mpz_t CGRS_n_large_Gcost, CGRS_n_large_Dcost, CGRS_n_large_Wcost;
extern mpz_t CGRS_n_large_Gcost_r, CGRS_n_large_Dcost_r, CGRS_n_large_Wcost_r;
extern mpz_t CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Dcost, CGRSGloop_m_large_ancila;
extern mpz_t CGRS_m_large_Gcost, CGRS_m_large_Dcost, CGRS_m_large_Wcost;
extern mpz_t CGRS_m_large_Gcost_r, CGRS_m_large_Dcost_r, CGRS_m_large_Wcost_r;

void init_Rank_CBC_params();

#endif