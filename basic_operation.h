#ifndef _BASIC_OPERATION_H_
#define _BASIC_OPERATION_H_

#include <stdio.h>
#include <math.h>
#include <gmp.h>

void binomial_coefficient(mpz_t coefficient, mpz_t x, mpz_t y);

void compute_q_analogue(mpz_t q_analogue, mpz_t m);

void compute_q_binomial(mpz_t q_binomial, mpz_t m, mpz_t r);

void compute_add_cost(mpz_t add_Gcost, mpz_t add_Dcost, mpz_t add_ancila, mpz_t m);

void compute_matmul_cost(mpz_t matmul_Gcost, mpz_t matmul_Dcost, mpz_t matmul_ancila, mpz_t l, mpz_t m, mpz_t n);

void compute_GE_cost(mpz_t GE_Gcost, mpz_t GE_Dcost, mpz_t GE_ancila, mpz_t r, mpz_t n);

void compute_Ham_cost(mpz_t Ham_Gcost, mpz_t Ham_Dcost, mpz_t Ham_ancila, mpz_t m);

void compute_addpow_cost(mpz_t addpow_Gcost, mpz_t addpow_Dcost, mpz_t addpow_ancila, mpz_t m);

void compute_mulpow_cost(mpz_t mulpow_Gcost, mpz_t mulpow_Dcost, mpz_t mulpow_ancila, mpz_t m);

void compute_rankpow_cost(mpz_t rankpow_Gcost, mpz_t rankpow_Dcost, mpz_t rankpow_ancila, mpz_t n, mpz_t m);

void compute_QRA_cost(mpz_t QRA_Gcost, mpz_t QRA_Dcost, mpz_t QRA_ancila, mpz_t n, mpz_t m);

void compute_Dicke_cost(mpz_t Dicke_Gcost, mpz_t Dicke_Dcost, mpz_t Dicke_ancila, mpz_t n, mpz_t r);

void compute_sort_cost(mpz_t sort_Gcost, mpz_t sort_Dcost, mpz_t sort_ancila, mpz_t n, mpz_t k, mpz_t r);

void compute_spDicke_cost(mpz_t spDicke_Gcost, mpz_t spDicke_Dcost, mpz_t spDicke_ancila, mpz_t n, mpz_t k, mpz_t r);

void compute_sp_cost(mpz_t sp_Gcost, mpz_t sp_Dcost, mpz_t sp_ancila, mpz_t Vsize);

void compute_OPFDicke_cost(mpz_t OPFDicke_Gcost, mpz_t OPFDicke_Dcost, mpz_t OPFDicke_ancila, mpz_t x);

void compute_OPF_cost(mpz_t OPF_Gcost, mpz_t OPF_Dcost, mpz_t OPF_ancila, mpz_t Vsize, mpz_t Msize);

void compute_dif_cost(mpz_t dif_Gcost, mpz_t dif_Dcost, mpz_t dif_ancila, mpz_t Vsize);

void compute_matmul_Tcost(mpz_t matmul_Gcost, mpz_t matmul_Dcost, mpz_t matmul_ancila, mpz_t l, mpz_t m, mpz_t n);

void compute_GE_Tcost(mpz_t GE_Gcost, mpz_t GE_Dcost, mpz_t GE_ancila, mpz_t r, mpz_t n);

void compute_Ham_Tcost(mpz_t Ham_Gcost, mpz_t Ham_Dcost, mpz_t Ham_ancila, mpz_t m);

void compute_fac(mpz_t nfac, mpz_t n);

int compute_log(mpz_t cost);

#endif