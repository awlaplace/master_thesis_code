#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "basic_operation.h"


/*
正整数 x, y を引数にとったとき，x と y の2項係数を返す関数
*/
void binomial_coefficient(mpz_t coefficient, mpz_t x, mpz_t y){
    if(mpz_cmp_ui(y, 0) > 0){
        mpz_t z;
        mpz_init(z);
        mpz_sub(z, x, y);
        mpz_add_ui(z, z, 1);

        mpz_mul(coefficient, coefficient, z);

        mpz_t i;
        mpz_init_set_ui(i, 1);

        while(mpz_cmp(i, y) == -1){
            mpz_add_ui(z, z, 1);
            mpz_add_ui(i, i, 1);
            
            mpz_mul(coefficient, coefficient, z);
            mpz_cdiv_q(coefficient, coefficient, i);
        }
    }
}


/*
  q_analogue, m を引数にとったとき，m の q 類似を q_analogue に格納する関数
  * q = 2 のとき，m の q 類似は 2^m - 1
*/
void compute_q_analogue(mpz_t q_analogue, mpz_t m){
    mpz_t loop_index;
    mpz_init_set_ui(loop_index, 0);

    while(mpz_cmp(loop_index, m) <= 0){
        mpz_mul_ui(q_analogue, q_analogue, 2);
        mpz_add_ui(loop_index, loop_index, 1);
    }
    mpz_sub_ui(q_analogue, q_analogue, 1);
}

/*
  q_binomial, m, r を引数にとったとき，
  m と r の q 二項係数を q_binomial に格納する関数
*/
void compute_q_binomial(mpz_t q_binomial, mpz_t m, mpz_t r){
    mpz_t q_binomial_denominator, q_binomial_numerator;
    mpz_init(q_binomial_denominator);
    mpz_init(q_binomial_numerator);

    mpz_t up_counter, down_counter;
    mpz_init_set_ui(up_counter, 1);
    mpz_init_set(down_counter, m);

    mpz_set_ui(q_binomial, 1);

    while(mpz_cmp(up_counter, r) <= 0){
        mpz_set_ui(q_binomial_numerator, 1);
        compute_q_analogue(q_binomial_numerator, down_counter);
        mpz_mul(q_binomial, q_binomial, q_binomial_numerator);

        mpz_set_ui(q_binomial_denominator, 1);
        compute_q_analogue(q_binomial_denominator, up_counter);
        mpz_cdiv_q(q_binomial, q_binomial, q_binomial_denominator);

        mpz_sub_ui(down_counter, down_counter, 1);
        mpz_add_ui(up_counter, up_counter, 1);
    }
}


/*
 add_Gcost, add_Dcost, add_ancila, m を引数にとったとき，
 add_Gcost に 2m, add_Dcost に 2，add_ancila に m を格納する関数
*/
void compute_add_cost(mpz_t add_Gcost, mpz_t add_Dcost, mpz_t add_ancila, mpz_t m){
    mpz_set_ui(add_Gcost, 2);
    mpz_mul(add_Gcost, add_Gcost, m);

    mpz_set_ui(add_Dcost, 2);

    mpz_set(add_ancila, m);
}


/*
 matmul_Gcost, matmul_Dcost, matmul_ancila, l, m, n を引数にとったとき，
 matmul_Gcost に 24 l m n, matmul_Dcost に 16 m，matmul_ancila に l n を格納する関数
*/
void compute_matmul_cost(mpz_t matmul_Gcost, mpz_t matmul_Dcost, mpz_t matmul_ancila, mpz_t l, mpz_t m, mpz_t n){
    mpz_set_ui(matmul_Gcost, 24);
    mpz_mul(matmul_Gcost, matmul_Gcost, l);
    mpz_mul(matmul_Gcost, matmul_Gcost, m);
    mpz_mul(matmul_Gcost, matmul_Gcost, n);

    mpz_set_ui(matmul_Dcost, 16);
    mpz_mul(matmul_Dcost, matmul_Dcost, m);

    mpz_set(matmul_ancila, l);
    mpz_mul(matmul_ancila, matmul_ancila, n);
}


/*
 GE_Gcost, GE_Dcost, GE_ancila, r, n を引数にとったとき，
 GE_Gcost に 4(r - 1)r(3 n - r + 5), GE_Dcost に 16 (r - 1)，GE_ancila に r^2 を格納する関数
*/
void compute_GE_cost(mpz_t GE_Gcost, mpz_t GE_Dcost, mpz_t GE_ancila, mpz_t r, mpz_t n){
    mpz_t rS1, nM3SrA5;
    mpz_init(rS1);
    mpz_sub_ui(rS1, r, 1);
    mpz_init_set_ui(nM3SrA5, 5);
    mpz_addmul_ui(nM3SrA5, n, 3);
    mpz_sub(nM3SrA5, nM3SrA5, r);

    mpz_set_ui(GE_Gcost, 4);
    mpz_mul(GE_Gcost, GE_Gcost, r);
    mpz_mul(GE_Gcost, GE_Gcost, rS1);
    mpz_mul(GE_Gcost, GE_Gcost, nM3SrA5);

    mpz_set_ui(GE_Dcost, 16);
    mpz_mul(GE_Dcost, GE_Dcost, rS1);

    mpz_set(GE_ancila, r);
    mpz_mul(GE_ancila, GE_ancila, r);
}


/*
 Ham_Gcost, Ham_Dcost, Ham_ancila, m を引数にとったとき，
 Ham_Gcost に 51(m - 1), Ham_Dcost に 32，Ham_ancila に 2(m - 1) + log_2 (m) を格納する関数
*/
void compute_Ham_cost(mpz_t Ham_Gcost, mpz_t Ham_Dcost, mpz_t Ham_ancila, mpz_t m){
    mpz_t mS1, cp_m, log_m;
    mpz_init_set(cp_m, m);
    mpz_init_set_ui(mS1, 0);
    mpz_sub_ui(mS1, m, 1);
    mpz_init(log_m);
    int logm = compute_log(cp_m);
    mpz_set_ui(log_m, logm);

    mpz_set_ui(Ham_Gcost, 51);
    mpz_mul(Ham_Gcost, Ham_Gcost, mS1);

    mpz_set_ui(Ham_Dcost, 32);

    mpz_set_ui(Ham_ancila, 2);
    mpz_mul(Ham_ancila, Ham_ancila, mS1);
    mpz_add(Ham_ancila, Ham_ancila, log_m);
}


/*
 addpow_Gcost, addpow_Dcost, addpow_ancila, m を引数にとったとき，
 addpow_Gcost に 26 m - 23, addpow_Dcost に 16，addpow_ancila に 2 m - 1 を格納する関数
*/
void compute_addpow_cost(mpz_t addpow_Gcost, mpz_t addpow_Dcost, mpz_t addpow_ancila, mpz_t m){
    mpz_set_ui(addpow_Gcost, 26);
    mpz_mul(addpow_Gcost, addpow_Gcost, m);
    mpz_sub_ui(addpow_Gcost, addpow_Gcost, 23);

    mpz_set_ui(addpow_Dcost, 16);

    mpz_set_ui(addpow_ancila, 2);
    mpz_mul(addpow_ancila, addpow_ancila, m);
    mpz_sub_ui(addpow_ancila, addpow_ancila, 1);
}


/*
 mulpow_Gcost, mulpow_Dcost, mulpow_ancila, m を引数にとったとき，
 mulpow_Gcost に 12 m (m + 1), mulpow_Dcost に 16 m，mulpow_ancila に m を格納する関数
*/
void compute_mulpow_cost(mpz_t mulpow_Gcost, mpz_t mulpow_Dcost, mpz_t mulpow_ancila, mpz_t m){
    mpz_t mA1;
    mpz_init_set_ui(mA1, 1);
    mpz_add(mA1, mA1, m);

    mpz_set_ui(mulpow_Gcost, 12);
    mpz_mul(mulpow_Gcost, mulpow_Gcost, m);
    mpz_mul(mulpow_Gcost, mulpow_Gcost, mA1);

    mpz_set_ui(mulpow_Dcost, 16);
    mpz_mul(mulpow_Dcost, mulpow_Dcost, m);

    mpz_set(mulpow_ancila, m);
}


/*
 rankpow_Gcost, rankpow_Dcost, rankpow_ancila, n, m を引数にとったとき，
 rankpow_Gcost に 4 (n - 1)n(3 n - m + 5), rankpow_Dcost に 16 (m + 1)，rankpow_ancila に n m を格納する関数
*/
void compute_rankpow_cost(mpz_t rankpow_Gcost, mpz_t rankpow_Dcost, mpz_t rankpow_ancila, mpz_t n, mpz_t m){
    mpz_t nS1, nM3SmA5;
    mpz_init(nS1);
    mpz_sub_ui(nS1, n, 1);
    mpz_init_set_ui(nM3SmA5, 5);
    mpz_addmul_ui(nM3SmA5, n, 3);
    mpz_sub(nM3SmA5, nM3SmA5, m);

    mpz_set_ui(rankpow_Gcost, 4);
    mpz_mul(rankpow_Gcost, rankpow_Gcost, nS1);
    mpz_mul(rankpow_Gcost, rankpow_Gcost, n);
    mpz_mul(rankpow_Gcost, rankpow_Gcost, nM3SmA5);

    mpz_set_ui(rankpow_Dcost, 16);
    mpz_addmul_ui(rankpow_Dcost, m, 16);

    mpz_set(rankpow_ancila, n);
    mpz_mul(rankpow_ancila, m, m);
}


/*
 QRA_Gcost, QRA_Dcost, QRA_ancila, n, m を引数にとったとき，
 QRA_Gcost に n m + n log_2 (n), QRA_Dcost に log_2 (m) + log_2 (n)，QRA_ancila に n m + n log_2 (n) を格納する関数
*/
void compute_QRA_cost(mpz_t QRA_Gcost, mpz_t QRA_Dcost, mpz_t QRA_ancila, mpz_t n, mpz_t m){
    mpz_t cp_m, cp_n, log_m, log_n;
    mpz_init_set(cp_m, m);
    mpz_init(log_m);
    int logm = compute_log(cp_m);
    mpz_set_ui(log_m, logm);
    mpz_init_set(cp_n, n);
    mpz_init(log_n);
    int logn = compute_log(cp_n);
    mpz_set_ui(log_n, logn);

    mpz_set(QRA_Gcost, n);
    mpz_mul(QRA_Gcost, QRA_Gcost, log_n);
    mpz_addmul(QRA_Gcost, n, m);

    mpz_set(QRA_Dcost, log_m);
    mpz_add(QRA_Dcost, QRA_Dcost, log_n);

    mpz_set(QRA_ancila, QRA_Gcost);
}


/*
 Dicke_Gcost, Dicke_Dcost, Dicke_ancila, n, r を引数にとったとき，
 Dicke_Gcost に 1221 n r - 1221 r^2 - 601 n + 4 r + 304, 
 Dicke_Dcost に (27 n r - 12 n - 27 r^2 + 3) / (r - 2), 
 Dicke_ancila に 0 を格納する関数
*/
void compute_Dicke_cost(mpz_t Dicke_Gcost, mpz_t Dicke_Dcost, mpz_t Dicke_ancila, mpz_t n, mpz_t r){
    mpz_t rS2, nSr, c601n, c12n;
    mpz_init(rS2);
    mpz_sub_ui(rS2, r, 2);
    mpz_init(nSr);
    mpz_sub(nSr, n, r);
    mpz_init_set_ui(c601n, 601);
    mpz_mul(c601n, c601n, n);
    mpz_init_set_ui(c12n, 12);
    mpz_mul(c12n, c12n, n);

    mpz_set_ui(Dicke_Gcost, 1221);
    mpz_mul(Dicke_Gcost, Dicke_Gcost, r);
    mpz_mul(Dicke_Gcost, Dicke_Gcost, nSr);
    mpz_sub(Dicke_Gcost, Dicke_Gcost, c601n);
    mpz_addmul_ui(Dicke_Gcost, r, 4);
    mpz_add_ui(Dicke_Gcost, Dicke_Gcost, 304);

    mpz_set_ui(Dicke_Dcost, 27);
    mpz_mul(Dicke_Dcost, Dicke_Dcost, r);
    mpz_mul(Dicke_Dcost, Dicke_Dcost, nSr);
    mpz_sub(Dicke_Dcost, Dicke_Dcost, c12n);
    mpz_cdiv_q(Dicke_Dcost, Dicke_Dcost, rS2);

    mpz_set_ui(Dicke_ancila, 0);
}


/*
 sort_Gcost, sort_Dcost, sort_ancila, n, k, r を引数にとったとき，
 sort_Gcost に 2(n - 1)log_2 (n)(log_2 (n) - 1), 
 sort_Dcost に 2(n - 1)log_2 (n)(log_2 (n) - 1), 
 sort_ancila に (n - 1)log_2 (n)(log_2 (n) - 1)(r + 1) を格納する関数
*/
void compute_sort_cost(mpz_t sort_Gcost, mpz_t sort_Dcost, mpz_t sort_ancila, mpz_t n, mpz_t k, mpz_t r){
    mpz_t nS1, cp_n, log_n, log_nS1, rA1;
    mpz_init(nS1);
    mpz_sub_ui(nS1, n, 1);
    mpz_init_set(cp_n, n);
    mpz_init(log_n);
    int logn = compute_log(cp_n);
    mpz_set_ui(log_n, logn);
    mpz_init(log_nS1);
    mpz_sub_ui(log_nS1, log_n, 1);
    mpz_init_set_ui(rA1, 1);
    mpz_add(rA1, rA1, r);

    mpz_t base_cost;
    mpz_init_set(base_cost, nS1);
    mpz_mul(base_cost, base_cost, log_n);
    mpz_mul(base_cost, base_cost, log_nS1);

    mpz_set(sort_Gcost, base_cost);
    mpz_mul_ui(sort_Gcost, sort_Gcost, 2);

    mpz_set(sort_Dcost, sort_Gcost);

    mpz_set(sort_ancila, rA1);
    mpz_mul(sort_ancila, sort_ancila, base_cost);
}


/*
 spDicke_Gcost, spDicke_Dcost, spDicke_ancila, n, k, r を引数にとったとき，
 spDicke_Gcost に Dicke_Gcost + sort_Gcost, 
 spDicke_Dcost に max{Dicke_Dcost, sort_Dcost}, 
 spDicke_ancila に Dicke_ancila + sort_ancila を格納する関数
*/
void compute_spDicke_cost(mpz_t spDicke_Gcost, mpz_t spDicke_Dcost, mpz_t spDicke_ancila, mpz_t n, mpz_t k, mpz_t r){
    mpz_t Dicke_Gcost, Dicke_Dcost, Dicke_ancila;
    mpz_t sort_Gcost, sort_Dcost, sort_ancila;
    mpz_init_set_ui(Dicke_Gcost, 0);
    mpz_init_set_ui(Dicke_Dcost, 0);
    mpz_init_set_ui(Dicke_ancila, 0);
    mpz_init_set_ui(sort_Gcost, 0);
    mpz_init_set_ui(sort_Dcost, 0);
    mpz_init_set_ui(sort_ancila, 0);

    compute_Dicke_cost(Dicke_Gcost, Dicke_Dcost, Dicke_ancila, n, r);
    compute_sort_cost(sort_Gcost, sort_Dcost, sort_ancila, n, k, r);

    mpz_set(spDicke_Gcost, Dicke_Gcost);
    mpz_add(spDicke_Gcost, spDicke_Gcost, sort_Gcost);

    if(mpz_cmp(Dicke_Dcost, sort_Dcost) > 0){
        mpz_set(spDicke_Dcost, Dicke_Dcost);
    }
    else{
        mpz_set(spDicke_Dcost, sort_Dcost);
    }

    mpz_set(spDicke_ancila, Dicke_ancila);
    mpz_add(spDicke_ancila, spDicke_ancila, sort_ancila);
}


/*
 sp_Gcost, sp_Dcost, sp_ancila, Vsize を引数にとったとき，
 sp_Gcost に Vsize, 
 sp_Dcost に 1, 
 sp_ancila に 0 を格納する関数
*/
void compute_sp_cost(mpz_t sp_Gcost, mpz_t sp_Dcost, mpz_t sp_ancila, mpz_t Vsize){
    mpz_set(sp_Gcost, Vsize);

    mpz_set_ui(sp_Dcost, 1);

    mpz_set_ui(sp_ancila, 0);
}


/*
 OPFDicke_Gcost, OPFDicke_Dcost, OPFDicke_ancila, x を引数にとったとき，
 OPFDicke_Gcost に 48 x - 46, 
 OPFDicke_Dcost に 32, 
 OPFDicke_ancila に x を格納する関数
*/
void compute_OPFDicke_cost(mpz_t OPFDicke_Gcost, mpz_t OPFDicke_Dcost, mpz_t OPFDicke_ancila, mpz_t x){
    mpz_set(OPFDicke_Gcost, x);
    mpz_mul_ui(OPFDicke_Gcost, OPFDicke_Gcost, 48);
    mpz_sub_ui(OPFDicke_Gcost, OPFDicke_Gcost, 46);

    mpz_set_ui(OPFDicke_Dcost, 32);

    mpz_set(OPFDicke_ancila, x);
}


/*
 OPF_Gcost, OPF_Dcost, OPF_ancila, Vsize, Msize を引数にとったとき，
 OPF_Gcost に 4 Msize, 
 OPF_Dcost に 4, 
 OPF_ancila に 0 を格納する関数
*/
void compute_OPF_cost(mpz_t OPF_Gcost, mpz_t OPF_Dcost, mpz_t OPF_ancila, mpz_t Vsize, mpz_t Msize){
    mpz_set_ui(OPF_Gcost, 4);
    mpz_mul(OPF_Gcost, OPF_Gcost, Msize);

    mpz_set_ui(OPF_Dcost, 4);

    mpz_set_ui(OPF_ancila, 0);
}


/*
 dif_Gcost, dif_Dcost, dif_ancila, Vsize を引数にとったとき，
 dif_Gcost に 48 Vsize - 46, 
 dif_Dcost に 32, 
 dif_ancila に Vsize を格納する関数
*/
void compute_dif_cost(mpz_t dif_Gcost, mpz_t dif_Dcost, mpz_t dif_ancila, mpz_t Vsize){
    mpz_set(dif_Gcost, Vsize);
    mpz_mul_ui(dif_Gcost, dif_Gcost, 48);
    mpz_sub_ui(dif_Gcost, dif_Gcost, 46);

    mpz_set_ui(dif_Dcost, 32);

    mpz_set(dif_ancila, Vsize);
}


/*
 matmul_Gcost, matmul_Dcost, matmul_ancila, l, m, n を引数にとったとき，
 matmul_Gcost に 24 l m n, matmul_Dcost に 16 m，matmul_ancila に l n を格納する関数
*/
void compute_matmul_Tcost(mpz_t matmul_Gcost, mpz_t matmul_Dcost, mpz_t matmul_ancila, mpz_t l, mpz_t m, mpz_t n){
    mpz_set_ui(matmul_Gcost, 13);
    mpz_mul(matmul_Gcost, matmul_Gcost, l);
    mpz_mul(matmul_Gcost, matmul_Gcost, m);
    mpz_mul(matmul_Gcost, matmul_Gcost, n);

    mpz_set_ui(matmul_Dcost, 8);
    mpz_mul(matmul_Dcost, matmul_Dcost, m);

    mpz_set(matmul_ancila, l);
    mpz_mul(matmul_ancila, matmul_ancila, n);
}


/*
 GE_Gcost, GE_Dcost, GE_ancila, r, n を引数にとったとき，
 GE_Gcost に 4(r - 1)r(3 n - r + 5), GE_Dcost に 16 (r - 1)，GE_ancila に r^2 を格納する関数
*/
void compute_GE_Tcost(mpz_t GE_Gcost, mpz_t GE_Dcost, mpz_t GE_ancila, mpz_t r, mpz_t n){
    mpz_t rS1, nM3SrA5;
    mpz_init(rS1);
    mpz_sub_ui(rS1, r, 1);
    mpz_init_set_ui(nM3SrA5, 5);
    mpz_addmul_ui(nM3SrA5, n, 3);
    mpz_sub(nM3SrA5, nM3SrA5, r);

    mpz_set_ui(GE_Gcost, 13);
    mpz_mul(GE_Gcost, GE_Gcost, r);
    mpz_mul(GE_Gcost, GE_Gcost, rS1);
    mpz_mul(GE_Gcost, GE_Gcost, nM3SrA5);
    mpz_cdiv_q_ui(GE_Gcost, GE_Gcost, 6);

    mpz_set_ui(GE_Dcost, 8);
    mpz_mul(GE_Dcost, GE_Dcost, rS1);

    mpz_set(GE_ancila, r);
    mpz_mul(GE_ancila, GE_ancila, r);
}


/*
 Ham_Gcost, Ham_Dcost, Ham_ancila, m を引数にとったとき，
 Ham_Gcost に 51(m - 1), Ham_Dcost に 32，Ham_ancila に 2(m - 1) + log_2 (m) を格納する関数
*/
void compute_Ham_Tcost(mpz_t Ham_Gcost, mpz_t Ham_Dcost, mpz_t Ham_ancila, mpz_t m){
    mpz_t mS1, cp_m, log_m;
    mpz_init_set(cp_m, m);
    mpz_init_set_ui(mS1, 0);
    mpz_sub_ui(mS1, m, 1);
    mpz_init(log_m);
    int logm = compute_log(cp_m);
    mpz_set_ui(log_m, logm);

    mpz_set_ui(Ham_Gcost, 26);
    mpz_mul(Ham_Gcost, Ham_Gcost, mS1);

    mpz_set_ui(Ham_Dcost, 16);

    mpz_set_ui(Ham_ancila, 2);
    mpz_mul(Ham_ancila, Ham_ancila, mS1);
    mpz_add(Ham_ancila, Ham_ancila, log_m);
}

/*
 nfac, n を引数にとったとき，n_fac に n! を格納する関数
*/
void compute_fac(mpz_t nfac, mpz_t n){
    mpz_t fac_index;
    mpz_init_set_ui(fac_index, 1);
    mpz_set_ui(nfac, 1);
    while(mpz_cmp(fac_index, n) < 0){
        mpz_mul(nfac, nfac, fac_index);
        mpz_add_ui(fac_index, fac_index, 1);
    }
}


/*
 cost を引数にとったとき，log2 cost を返す関数
*/
int compute_log(mpz_t cost){
    int cost_log = 0;
    while(mpz_cmp_ui(cost, 1) > 0){
        mpz_cdiv_q_ui(cost, cost, 2);
        cost_log++;
    }

    return cost_log;
}