/*
  gcc -c classical_GRS.c -lgmp && gcc -o classical_GRS ../../basic_operation.o ../Rank_CBC_params.o GRS_params.o classical_GRS.o -lgmp && ./classical_GRS
*/

/*
  * the variable nAw denotes n+w
  * the variable nSw denotes n-w
  * the variable nMw denotes n*w
  * the variable nDw denotes n/w
  * the variable n_bc_k denotes binomial coefficient nCk 
*/
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "GRS.h"


int main(){
    init_Rank_CBC_params();

    int rnd_r = 0;

    for(int index = 0; index <= 5; index++){
        mpz_init_set_ui(m, Rank_SDP_instance_list[index][0]);
        mpz_init_set_ui(n, Rank_SDP_instance_list[index][1]);
        mpz_init_set_ui(k, Rank_SDP_instance_list[index][2]);
        mpz_init_set_ui(w, Rank_SDP_instance_list[index][3]);
        mpz_init_set(r, w);

        init_GRS_computational_costs();

        if(index / 3 == 0){
            mpz_t range_end, nSk;
            mpz_init(range_end);
            mpz_init(nSk);
            mpz_sub(nSk, n, k);
            mpz_mul(range_end, m, nSk);
            mpz_cdiv_q(range_end, range_end, n);
            mpz_sub(range_end, m, range_end);
            while(mpz_cmp(r, range_end) < 0){
                init_GRS_computational_costs();

                // compute iteration cost
                mpz_init_set_ui(grover_iteration, 1);
                compute_grover_iteration_n_large(grover_iteration, m, r, w);

                mpz_init_set_ui(grover_processor_num, 1);

                compute_CGRSGloop_n_large(CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Dcost, CGRSGloop_n_large_ancila, m, n, k, r);

                compute_CGRS_n_large_Dcost(CGRS_n_large_Dcost_r, CGRSGloop_n_large_Dcost, grover_iteration, grover_processor_num);
                
                // grover_processor_num の個数を調整
                while(mpz_cmp(CGRS_n_large_Dcost_r, MAX_DEPTH) > 0){
                    mpz_mul_ui(grover_processor_num, grover_processor_num, 2);
                    compute_CGRS_n_large_Dcost(CGRS_n_large_Dcost_r, CGRSGloop_n_large_Dcost, grover_iteration, grover_processor_num);
                }

                compute_CGRS_n_large_Gcost(CGRS_n_large_Gcost_r, CGRSGloop_n_large_Gcost, grover_iteration, grover_processor_num);

                compute_CGRS_n_large_Wcost(CGRS_n_large_Wcost_r, CGRSGloop_n_large_ancila, grover_iteration, grover_processor_num, n, m, k);

                if((mpz_cmp_ui(CGRS_n_large_Gcost, 0) == 0) || (mpz_cmp(CGRS_n_large_Gcost, CGRS_n_large_Gcost_r) > 0)){
                    mpz_set(CGRS_n_large_Gcost, CGRS_n_large_Gcost_r);
                    mpz_set(CGRS_n_large_Dcost, CGRS_n_large_Dcost_r);
                    mpz_set(CGRS_n_large_Wcost, CGRS_n_large_Wcost_r);
                }
                
                mpz_add_ui(r, r, 1);
            }
            int grover_iteration_log = compute_log(grover_iteration);
            int CGRS_n_large_Gcost_log = compute_log(CGRS_n_large_Gcost);
            int CGRS_n_large_Dcost_log = compute_log(CGRS_n_large_Dcost);
            int CGRS_n_large_Wcost_log = compute_log(CGRS_n_large_Wcost);
            gmp_printf("n : %Zd, m : %Zd, k : %Zd, w : %Zd, Gcost_log : %d, Dcost_log : %d, Wcost_log : %d\n", n, m, k, w, CGRS_n_large_Gcost_log, CGRS_n_large_Dcost_log, CGRS_n_large_Wcost_log);
        }
        else{
            mpz_t range_end;
            mpz_init(range_end);
            mpz_sub(range_end, n, k);
            while(mpz_cmp(r, range_end) < 0){
                init_GRS_computational_costs();

                // compute iteration cost
                mpz_init_set_ui(grover_iteration, 1);
                compute_grover_iteration_m_large(grover_iteration, n, r, w);

                mpz_init_set_ui(grover_processor_num, 1);

                compute_CGRSGloop_m_large(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Dcost, CGRSGloop_m_large_ancila, m, n, k, r);

                compute_CGRS_m_large_Dcost(CGRS_m_large_Dcost_r, CGRSGloop_m_large_Dcost, grover_iteration, grover_processor_num);
                
                // grover_processor_num の個数を調整
                while(mpz_cmp(CGRS_m_large_Dcost_r, MAX_DEPTH) > 0){
                    mpz_mul_ui(grover_processor_num, grover_processor_num, 2);
                    compute_CGRS_m_large_Dcost(CGRS_m_large_Dcost_r, CGRSGloop_m_large_Dcost, grover_iteration, grover_processor_num);
                }

                compute_CGRS_m_large_Gcost(CGRS_m_large_Gcost_r, CGRSGloop_m_large_Gcost, grover_iteration, grover_processor_num);

                compute_CGRS_m_large_Wcost(CGRS_m_large_Wcost_r, CGRSGloop_m_large_ancila, grover_iteration, grover_processor_num, n, m, k);
                
                if((mpz_cmp_ui(CGRS_m_large_Gcost, 0) == 0) || (mpz_cmp(CGRS_m_large_Gcost, CGRS_m_large_Gcost_r) > 0)){
                    mpz_set(CGRS_m_large_Gcost, CGRS_m_large_Gcost_r);
                    mpz_set(CGRS_m_large_Dcost, CGRS_m_large_Dcost_r);
                    mpz_set(CGRS_m_large_Wcost, CGRS_m_large_Wcost_r);
                }
                
                mpz_add_ui(r, r, 1);
            }
            int grover_iteration_log = compute_log(grover_iteration);
            int CGRS_m_large_Gcost_log = compute_log(CGRS_m_large_Gcost);
            int CGRS_m_large_Dcost_log = compute_log(CGRS_m_large_Dcost);
            int CGRS_m_large_Wcost_log = compute_log(CGRS_m_large_Wcost);
            gmp_printf("n : %Zd, m : %Zd, k : %Zd, w : %Zd, Gcost_log : %d, Dcost_log : %d, Wcost_log : %d\n", n, m, k, w, CGRS_m_large_Gcost_log, CGRS_m_large_Dcost_log, CGRS_m_large_Wcost_log);
        }
    }
}


void compute_grover_iteration_n_large(mpz_t grover_iteration, mpz_t m, mpz_t r, mpz_t w){
    mpz_t grover_iteration_n_large_denominator, grover_iteration_n_large_numerator;
    mpz_init_set_ui(grover_iteration_n_large_denominator, 1);
    mpz_init_set_ui(grover_iteration_n_large_numerator, 1);

    compute_q_binomial(grover_iteration_n_large_numerator, m, w);
    compute_q_binomial(grover_iteration_n_large_denominator, r, w);

    mpz_set(grover_iteration, grover_iteration_n_large_numerator);
    mpz_cdiv_q(grover_iteration, grover_iteration, grover_iteration_n_large_denominator);
}


void compute_grover_iteration_m_large(mpz_t grover_iteration, mpz_t n, mpz_t r, mpz_t w){
    mpz_t grover_iteration_m_large_denominator, grover_iteration_m_large_numerator;
    mpz_init_set_ui(grover_iteration_m_large_denominator, 1);
    mpz_init_set_ui(grover_iteration_m_large_numerator, 1);

    compute_q_binomial(grover_iteration_m_large_numerator, n, w);
    compute_q_binomial(grover_iteration_m_large_denominator, r, w);

    mpz_set(grover_iteration, grover_iteration_m_large_numerator);
    mpz_cdiv_q(grover_iteration, grover_iteration, grover_iteration_m_large_denominator);
}


void compute_CGRSGloop_n_large(mpz_t CGRSGloop_n_large_Gcost, mpz_t CGRSGloop_n_large_Dcost, mpz_t CGRSGloop_n_large_ancila, mpz_t m, mpz_t n, mpz_t k, mpz_t r){
    mpz_t rankpow_Gcost, rankpow_Dcost, rankpow_ancila;
    mpz_t mulpow_Gcost1, mulpow_Dcost1, mulpow_ancila1;
    mpz_t GE_Gcost, GE_Dcost, GE_ancila;
    mpz_t mulpow_Gcost2, mulpow_Dcost2, mulpow_ancila2;
    mpz_t addpow_Gcost, addpow_Dcost, addpow_ancila;
    mpz_init_set_ui(rankpow_Gcost, 0);
    mpz_init_set_ui(rankpow_Dcost, 0);
    mpz_init_set_ui(rankpow_ancila, 0);
    mpz_init_set_ui(mulpow_Gcost1, 0);
    mpz_init_set_ui(mulpow_Dcost1, 0);
    mpz_init_set_ui(mulpow_ancila1, 0);
    mpz_init_set_ui(GE_Gcost, 0);
    mpz_init_set_ui(GE_Dcost, 0);
    mpz_init_set_ui(GE_ancila, 0);
    mpz_init_set_ui(mulpow_Gcost2, 0);
    mpz_init_set_ui(mulpow_Dcost2, 0);
    mpz_init_set_ui(mulpow_ancila2, 0);
    mpz_init_set_ui(addpow_Gcost, 0);
    mpz_init_set_ui(addpow_Dcost, 0);
    mpz_init_set_ui(addpow_ancila, 0);

    mpz_t nSk, nMr, nSkMnMr, mMnSk, nMrS1;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init(nMr);
    mpz_mul(nMr, n, r);
    mpz_init(nSkMnMr);
    mpz_mul(nSkMnMr, nSk, nMr);
    mpz_init(mMnSk);
    mpz_mul(mMnSk, m, nSk);
    mpz_init(nMrS1);
    mpz_sub(nMrS1, nMr, n);

    compute_rankpow_cost(rankpow_Gcost, rankpow_Dcost, rankpow_ancila, n, m);
    compute_mulpow_cost(mulpow_Gcost1, mulpow_Dcost1, mulpow_ancila1, m);
    compute_GE_cost(GE_Gcost, GE_Dcost, GE_ancila, nMr, mMnSk);
    compute_mulpow_cost(mulpow_Gcost2, mulpow_Dcost2, mulpow_ancila2, m);
    compute_addpow_cost(addpow_Gcost, addpow_Dcost, addpow_ancila, m);

    mpz_mul(mulpow_Gcost1, mulpow_Gcost1, nSkMnMr);
    mpz_mul(mulpow_Dcost1, mulpow_Dcost1, nSkMnMr);
    mpz_mul(mulpow_ancila1, mulpow_ancila1, nSkMnMr);
    mpz_mul(mulpow_Gcost2, mulpow_Gcost2, nMr);
    mpz_mul(mulpow_Dcost2, mulpow_Dcost2, nMr);
    mpz_mul(mulpow_ancila2, mulpow_ancila2, nMr);
    mpz_mul(addpow_Gcost, addpow_Gcost, nMrS1);
    mpz_mul(addpow_Dcost, addpow_Dcost, nMrS1);
    mpz_mul(addpow_ancila, addpow_ancila, nMrS1);

    // CGRSGloop_n_large_Gcost を計算する
    mpz_set(CGRSGloop_n_large_Gcost, rankpow_Gcost);
    mpz_add(CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Gcost, mulpow_Gcost1);
    mpz_add(CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Gcost, GE_Gcost);
    mpz_add(CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Gcost, mulpow_Gcost2);
    mpz_add(CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Gcost, addpow_Gcost);

    // CGRSGloop_n_large_Dcost を計算する
    mpz_set(CGRSGloop_n_large_Dcost, rankpow_Dcost);
    if(mpz_cmp(mulpow_Dcost1, CGRSGloop_n_large_Dcost) > 0){
        mpz_set(CGRSGloop_n_large_Dcost, mulpow_Dcost1);
    }
    if(mpz_cmp(GE_Dcost, CGRSGloop_n_large_Dcost) > 0){
        mpz_set(CGRSGloop_n_large_Dcost, GE_Dcost);
    }
    if(mpz_cmp(mulpow_Dcost2, CGRSGloop_n_large_Dcost) > 0){
        mpz_set(CGRSGloop_n_large_Dcost, mulpow_Dcost2);
    }
    if(mpz_cmp(addpow_Dcost, CGRSGloop_n_large_Dcost) > 0){
        mpz_set(CGRSGloop_n_large_Dcost, addpow_Dcost);
    }

    // CGRSGloop_n_large_ancila を計算する
    mpz_set(CGRSGloop_n_large_ancila, rankpow_ancila);
    mpz_add(CGRSGloop_n_large_ancila, CGRSGloop_n_large_ancila, mulpow_ancila1);
    mpz_add(CGRSGloop_n_large_ancila, CGRSGloop_n_large_ancila, GE_ancila);
    mpz_add(CGRSGloop_n_large_ancila, CGRSGloop_n_large_ancila, mulpow_ancila2);
    mpz_add(CGRSGloop_n_large_ancila, CGRSGloop_n_large_ancila, addpow_ancila);
}


void compute_CGRS_n_large_Gcost(mpz_t CGRS_n_large_Gcost, mpz_t CGRSGloop_n_large_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num){
    mpz_set(CGRS_n_large_Gcost, CGRSGloop_n_large_Gcost);
    mpz_mul(CGRS_n_large_Gcost, CGRS_n_large_Gcost, grover_iteration);
    mpz_mul(CGRS_n_large_Gcost, CGRS_n_large_Gcost, grover_processor_num);
}


void compute_CGRS_n_large_Dcost(mpz_t CGRS_n_large_Dcost, mpz_t CGRSGloop_n_large_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num){
    mpz_set(CGRS_n_large_Dcost, CGRSGloop_n_large_Dcost);
    mpz_mul(CGRS_n_large_Dcost, CGRS_n_large_Dcost, grover_iteration);
    mpz_cdiv_q(CGRS_n_large_Dcost, CGRS_n_large_Dcost, grover_processor_num);
}


void compute_CGRS_n_large_Wcost(mpz_t CGRS_n_large_Wcost, mpz_t CGRSGloop_n_large_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t m, mpz_t k){
    mpz_t CGRS_n_large_input;
    mpz_init_set_ui(CGRS_n_large_input, 0);

    mpz_t nSk, nMm;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init(nMm);
    mpz_mul(nMm, n, m);

    // CGRS_n_large_input を計算する
    mpz_addmul(CGRS_n_large_input, nSk, nMm);
    mpz_addmul(CGRS_n_large_input, nSk, m);
    mpz_add(CGRS_n_large_input, CGRS_n_large_input, nMm);
    
    mpz_set(CGRS_n_large_Wcost, CGRSGloop_n_large_ancila);
    mpz_mul(CGRS_n_large_Wcost, CGRS_n_large_Wcost, grover_iteration);
    mpz_cdiv_q(CGRS_n_large_Wcost, CGRS_n_large_Wcost, grover_processor_num);
    mpz_cdiv_q(CGRS_n_large_Wcost, CGRS_n_large_Wcost, grover_processor_num);
    mpz_add(CGRS_n_large_Wcost, CGRS_n_large_Wcost, CGRS_n_large_input);
}


void compute_CGRSGloop_m_large(mpz_t CGRSGloop_m_large_Gcost, mpz_t CGRSGloop_m_large_Dcost, mpz_t CGRSGloop_m_large_ancila, mpz_t m, mpz_t n, mpz_t k, mpz_t r){
    mpz_t rankpow_Gcost, rankpow_Dcost, rankpow_ancila;
    mpz_t mulpow_Gcost1, mulpow_Dcost1, mulpow_ancila1;
    mpz_t addpow_Gcost1, addpow_Dcost1, addpow_ancila1;
    mpz_t GE_Gcost, GE_Dcost, GE_ancila;
    mpz_t mulpow_Gcost2, mulpow_Dcost2, mulpow_ancila2;
    mpz_t addpow_Gcost2, addpow_Dcost2, addpow_ancila2;
    mpz_t mulpow_Gcost3, mulpow_Dcost3, mulpow_ancila3;
    mpz_t addpow_Gcost3, addpow_Dcost3, addpow_ancila3;
    mpz_init_set_ui(rankpow_Gcost, 0);
    mpz_init_set_ui(rankpow_Dcost, 0);
    mpz_init_set_ui(rankpow_ancila, 0);
    mpz_init_set_ui(mulpow_Gcost1, 0);
    mpz_init_set_ui(mulpow_Dcost1, 0);
    mpz_init_set_ui(mulpow_ancila1, 0);
    mpz_init_set_ui(addpow_Gcost1, 0);
    mpz_init_set_ui(addpow_Dcost1, 0);
    mpz_init_set_ui(addpow_ancila1, 0);
    mpz_init_set_ui(GE_Gcost, 0);
    mpz_init_set_ui(GE_Dcost, 0);
    mpz_init_set_ui(GE_ancila, 0);
    mpz_init_set_ui(mulpow_Gcost2, 0);
    mpz_init_set_ui(mulpow_Dcost2, 0);
    mpz_init_set_ui(mulpow_ancila2, 0);
    mpz_init_set_ui(addpow_Gcost2, 0);
    mpz_init_set_ui(addpow_Dcost2, 0);
    mpz_init_set_ui(addpow_ancila2, 0);
    mpz_init_set_ui(mulpow_Gcost3, 0);
    mpz_init_set_ui(mulpow_Dcost3, 0);
    mpz_init_set_ui(mulpow_ancila3, 0);
    mpz_init_set_ui(addpow_Gcost3, 0);
    mpz_init_set_ui(addpow_Dcost3, 0);
    mpz_init_set_ui(addpow_ancila3, 0);

    mpz_t nSk, nMr, mMnSk, mMnSkMnMr, mMr, mMrS1, mMn;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init(nMr);
    mpz_mul(nMr, n, r);
    mpz_init(mMnSk);
    mpz_mul(mMnSk, m, nSk);
    mpz_init(mMnSkMnMr);
    mpz_mul(mMnSkMnMr, mMnSk, nMr);
    mpz_init(mMr);
    mpz_mul(mMr, m, r);
    mpz_init(mMrS1);
    mpz_sub_ui(mMrS1, mMr, 1);
    mpz_init(mMn);
    mpz_mul(mMn, m, n);

    compute_rankpow_cost(rankpow_Gcost, rankpow_Dcost, rankpow_ancila, n, m);
    compute_mulpow_cost(mulpow_Gcost1, mulpow_Dcost1, mulpow_ancila1, m);
    compute_addpow_cost(addpow_Gcost1, addpow_Dcost1, addpow_ancila1, m);
    compute_GE_cost(GE_Gcost, GE_Dcost, GE_ancila, mMr, mMnSk);
    compute_mulpow_cost(mulpow_Gcost2, mulpow_Dcost2, mulpow_ancila2, m);
    compute_addpow_cost(addpow_Gcost2, addpow_Dcost2, addpow_ancila2, m);
    compute_mulpow_cost(mulpow_Gcost3, mulpow_Dcost3, mulpow_ancila3, m);
    compute_addpow_cost(addpow_Gcost3, addpow_Dcost3, addpow_ancila3, m);

    mpz_mul(mulpow_Gcost1, mulpow_Gcost1, mMnSkMnMr);
    mpz_mul(mulpow_Dcost1, mulpow_Dcost1, mMnSkMnMr);
    mpz_mul(mulpow_ancila1, mulpow_ancila1, mMnSkMnMr);
    mpz_mul(addpow_Gcost1, addpow_Gcost1, mMnSkMnMr);
    mpz_mul(addpow_Dcost1, addpow_Dcost1, mMnSkMnMr);
    mpz_mul(addpow_ancila1, addpow_ancila1, mMnSkMnMr);
    mpz_mul(mulpow_Gcost2, mulpow_Gcost2, mMr);
    mpz_mul(mulpow_Dcost2, mulpow_Dcost2, mMr);
    mpz_mul(mulpow_ancila2, mulpow_ancila2, mMr);
    mpz_mul(addpow_Gcost2, addpow_Gcost2, mMrS1);
    mpz_mul(addpow_Dcost2, addpow_Dcost2, mMrS1);
    mpz_mul(addpow_ancila2, addpow_ancila2, mMrS1);
    mpz_mul(mulpow_Gcost2, mulpow_Gcost2, mMn);
    mpz_mul(mulpow_Dcost2, mulpow_Dcost2, mMn);
    mpz_mul(mulpow_ancila2, mulpow_ancila2, mMn);
    mpz_mul(addpow_Gcost2, addpow_Gcost2, mMn);
    mpz_mul(addpow_Dcost2, addpow_Dcost2, mMn);
    mpz_mul(addpow_ancila2, addpow_ancila2, mMn);

    // CGRSGloop_m_large_Gcost を計算する
    mpz_set(CGRSGloop_m_large_Gcost, rankpow_Gcost);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, mulpow_Gcost1);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, addpow_Gcost1);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, GE_Gcost);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, mulpow_Gcost2);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, addpow_Gcost2);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, mulpow_Gcost3);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, addpow_Gcost3);

    // CGRSGloop_m_large_Dcost を計算する
    mpz_set(CGRSGloop_m_large_Dcost, rankpow_Dcost);
    if(mpz_cmp(mulpow_Dcost1, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, mulpow_Dcost1);
    }
    if(mpz_cmp(addpow_Dcost1, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, addpow_Dcost1);
    }
    if(mpz_cmp(GE_Dcost, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, GE_Dcost);
    }
    if(mpz_cmp(mulpow_Dcost2, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, mulpow_Dcost2);
    }
    if(mpz_cmp(addpow_Dcost2, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, addpow_Dcost2);
    }
    if(mpz_cmp(mulpow_Dcost3, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, mulpow_Dcost3);
    }
    if(mpz_cmp(addpow_Dcost3, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, addpow_Dcost3);
    }

    // CGRSGloop_m_large_ancila を計算する
    mpz_set(CGRSGloop_m_large_ancila, rankpow_ancila);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, mulpow_ancila1);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, addpow_ancila1);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, GE_ancila);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, mulpow_ancila2);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, addpow_ancila2);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, mulpow_ancila3);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, addpow_ancila3);
}


void compute_CGRS_m_large_Gcost(mpz_t CGRS_m_large_Gcost, mpz_t CGRSGloop_m_large_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num){
    mpz_set(CGRS_m_large_Gcost, CGRSGloop_m_large_Gcost);
    mpz_mul(CGRS_m_large_Gcost, CGRS_m_large_Gcost, grover_iteration);
    mpz_mul(CGRS_m_large_Gcost, CGRS_m_large_Gcost, grover_processor_num);
}


void compute_CGRS_m_large_Dcost(mpz_t CGRS_m_large_Dcost, mpz_t CGRSGloop_m_large_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num){
    mpz_set(CGRS_m_large_Dcost, CGRSGloop_m_large_Dcost);
    mpz_mul(CGRS_m_large_Dcost, CGRS_m_large_Dcost, grover_iteration);
    mpz_cdiv_q(CGRS_m_large_Dcost, CGRS_m_large_Dcost, grover_processor_num);
}


void compute_CGRS_m_large_Wcost(mpz_t CGRS_m_large_Wcost, mpz_t CGRSGloop_m_large_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t m, mpz_t k){
    mpz_t CGRS_m_large_input;
    mpz_init_set_ui(CGRS_m_large_input, 0);

    mpz_t nSkMm, nSkMmMn, nMm;
    mpz_init(nSkMm);
    mpz_sub(nSkMm, n, k);
    mpz_mul(nSkMm, nSkMm, m);
    mpz_init(nSkMmMn);
    mpz_mul(nSkMmMn, nSkMm, n);
    mpz_init(nMm);
    mpz_mul(nMm, n, m);

    // CGRS_m_large_input を計算する
    mpz_init_set(CGRS_m_large_input, nSkMmMn);
    mpz_add(CGRS_m_large_input, CGRS_m_large_input, nSkMm);
    mpz_add(CGRS_m_large_input, CGRS_m_large_input, nMm);

    mpz_set(CGRS_m_large_Wcost, CGRSGloop_m_large_ancila);
    mpz_mul(CGRS_m_large_Wcost, CGRS_m_large_Wcost, grover_iteration);
    mpz_mul(CGRS_m_large_Wcost, CGRS_m_large_Wcost, grover_processor_num);
    mpz_mul(CGRS_m_large_Wcost, CGRS_m_large_Wcost, grover_processor_num);
    mpz_add(CGRS_m_large_Wcost, CGRS_m_large_Wcost, CGRS_m_large_input);
}