/*
  gcc -c improved_quantum_Prange.c -lgmp && gcc -o improved_quantum_Prange ../../basic_operation.o ../Hamming_CBC_params.o Prange_params.o improved_quantum_Prange.o -lgmp && ./improved_quantum_Prange
*/

#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "Prange.h"


int main(){
    init_Hamming_CBC_params();

    for(int index = 0; index < 9; index++){
        mpz_init_set_ui(n, SDP_instance_list[index][0]);
        mpz_init_set_ui(k, SDP_instance_list[index][1]);
        mpz_init_set_ui(w, SDP_instance_list[index][2]);

        mpz_init_set_ui(grover_iteration, 1);

        mpz_init_set_ui(grover_processor_num, 1);

        init_Prange_computational_costs();

        // 二項係数が意味を持つ範囲に限定
        if(mpz_cmp(n, w) > 0){
            compute_grover_iteration(grover_iteration, n, k, w);

            compute_CPGloop_cost(CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k);

            compute_QPGloop_cost(QPGloop_Gcost, QPGloop_Dcost, QPGloop_ancila, CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k, w);
            
            // IQPrange_Dcost を計算する
            compute_IQPrange_Dcost(IQPrange_Dcost, QPGloop_Dcost, CPGloop_Dcost, grover_iteration, grover_processor_num, n, k, w);

            // grover_processor_num の個数を調整
            while(mpz_cmp(IQPrange_Dcost, MAX_DEPTH) > 0){
                mpz_mul_ui(grover_processor_num, grover_processor_num, 2);
                compute_IQPrange_Dcost(IQPrange_Dcost, QPGloop_Dcost, CPGloop_Dcost, grover_iteration, grover_processor_num, n, k, w);
            }

            // IQPrange_Gcost を計算する
            compute_IQPrange_Gcost(IQPrange_Gcost, QPGloop_Gcost, CPGloop_Gcost, grover_iteration, grover_processor_num, n, k, w);
            
            // IQPrange_Wcost を計算する
            compute_IQPrange_Wcost(IQPrange_Wcost, QPGloop_ancila, CPGloop_ancila, grover_iteration, grover_processor_num, n, k, w);

            int grover_iteration_log = compute_log(grover_iteration);
            int IQPrange_Gcost_log = compute_log(IQPrange_Gcost);
            int IQPrange_Dcost_log = compute_log(IQPrange_Dcost);
            int IQPrange_Wcost_log = compute_log(IQPrange_Wcost);

            Grover_result_list[index][3] = IQPrange_Gcost_log;
            Grover_result_list[index][4] = IQPrange_Dcost_log;
            Grover_result_list[index][5] = IQPrange_Wcost_log;
        }
        printf("n : %d, k : %d, w : %d, Gcost : %d, Dcost : %d, Wcost : %d\n", Grover_result_list[index][0], Grover_result_list[index][1], Grover_result_list[index][2], Grover_result_list[index][3], Grover_result_list[index][4], Grover_result_list[index][5]);
    }
}


/*
 grover_iteration_cost, n, k, w を引数にとったとき，
 Grover のループ回数を grover_iteration_cost に格納する関数
*/
void compute_grover_iteration(mpz_t grover_iteration, mpz_t n, mpz_t k, mpz_t w){
    // the variable n_bc_w denotes nCw
    // the variable nSk_bc_w denotes n-kCw
    mpz_t nSk, n_bc_w, nSk_bc_w;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init_set_ui(n_bc_w, 1);
    mpz_init_set_ui(nSk_bc_w, 1);

    binomial_coefficient(n_bc_w, n, w);
    binomial_coefficient(nSk_bc_w, nSk, w);

    mpz_cdiv_q(grover_iteration, n_bc_w, nSk_bc_w);
    
    mpz_root(grover_iteration, grover_iteration, 2);
}


/*
 CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k を引数にとったとき，
 古典Prangeアルゴリズムの for 文のループ中の G-cost, D-cost, アンシラビット数を
 それぞれ CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila に格納する関数
*/
void compute_CPGloop_cost(mpz_t CPGloop_Gcost, mpz_t CPGloop_Dcost, mpz_t CPGloop_ancila, mpz_t n, mpz_t k){
    mpz_t matmul_Gcost1, matmul_Dcost1, matmul_ancila1;
    mpz_t GE_Gcost, GE_Dcost, GE_ancila;
    mpz_t matmul_Gcost2, matmul_Dcost2, matmul_ancila2;
    mpz_t Ham_Gcost, Ham_Dcost, Ham_ancila;
    mpz_init_set_ui(matmul_Gcost1, 0);
    mpz_init_set_ui(matmul_Dcost1, 0);
    mpz_init_set_ui(matmul_ancila1, 0);
    mpz_init_set_ui(GE_Gcost, 0);
    mpz_init_set_ui(GE_Dcost, 0);
    mpz_init_set_ui(GE_ancila, 0);
    mpz_init_set_ui(matmul_Gcost2, 0);
    mpz_init_set_ui(matmul_Dcost2, 0);
    mpz_init_set_ui(matmul_ancila2, 0);
    mpz_init_set_ui(Ham_Gcost, 0);
    mpz_init_set_ui(Ham_Dcost, 0);
    mpz_init_set_ui(Ham_ancila, 0);

    mpz_t nSk, one;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init_set_ui(one, 1);

    compute_matmul_cost(matmul_Gcost1, matmul_Dcost1, matmul_ancila1, nSk, n, n);
    compute_GE_cost(GE_Gcost, GE_Dcost, GE_ancila, nSk, n);
    compute_matmul_cost(matmul_Gcost2, matmul_Dcost2, matmul_ancila2, nSk, nSk, one);
    compute_Ham_cost(Ham_Gcost, Ham_Dcost, Ham_ancila, nSk);

    mpz_set(CPGloop_Gcost, matmul_Gcost1);
    mpz_add(CPGloop_Gcost, CPGloop_Gcost, GE_Gcost);
    mpz_add(CPGloop_Gcost, CPGloop_Gcost, matmul_Gcost2);
    mpz_add(CPGloop_Gcost, CPGloop_Gcost, Ham_Gcost);

    mpz_set(CPGloop_Dcost, matmul_Dcost2);
    if(mpz_cmp(GE_Dcost, CPGloop_Dcost) > 0){
        mpz_set(CPGloop_Dcost, GE_Dcost);
    }
    if(mpz_cmp(matmul_Dcost2, CPGloop_Dcost) > 0){
        mpz_set(CPGloop_Dcost, matmul_Dcost2);
    }
    if(mpz_cmp(Ham_Dcost, CPGloop_Dcost) > 0){
        mpz_set(CPGloop_Dcost, Ham_Dcost);
    }

    mpz_set(CPGloop_ancila, matmul_ancila1);
    mpz_add(CPGloop_ancila, CPGloop_ancila, GE_ancila);
    mpz_add(CPGloop_ancila, CPGloop_ancila, matmul_ancila2);
    mpz_add(CPGloop_ancila, CPGloop_ancila, Ham_ancila);
}


/*
 QPGloop_Gcost, QPGloop_Dcost, QPGloop_ancila, CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k, w を引数にとったとき，
 量子Prangeアルゴリズムの for 文のループ中の G-cost, D-cost, アンシラビット数を
 それぞれ QPGloop_Gcost, QPGloop_Dcost, QPGloop_ancila に格納する関数
*/
void compute_QPGloop_cost(mpz_t QPGloop_Gcost, mpz_t QPGloop_Dcost, mpz_t QPGloop_ancila, mpz_t CPGloop_Gcost, mpz_t CPGloop_Dcost, mpz_t CPGloop_ancila, mpz_t n, mpz_t k, mpz_t w){
    mpz_t nMn, nSk;
    mpz_init(nMn);
    mpz_mul(nMn, n, n);
    mpz_init(nSk);
    mpz_sub(nSk, n, k);

    mpz_t n_bc_w, nSk_bc_w;
    mpz_init_set_ui(n_bc_w, 1);
    binomial_coefficient(n_bc_w, n, w);
    mpz_init_set_ui(nSk_bc_w, 1);
    binomial_coefficient(nSk_bc_w, nSk, w);

    mpz_t Vsize, Msize;
    mpz_init(Vsize);
    mpz_cdiv_q(Vsize, n_bc_w, nSk_bc_w);
    mpz_init_set_ui(Msize, 1);

    mpz_t QRA_Gcost, QRA_Dcost, QRA_ancila;
    mpz_t OPF_Gcost, OPF_Dcost, OPF_ancila;
    mpz_t dif_Gcost, dif_Dcost, dif_ancila;
    mpz_init_set_ui(QRA_Gcost, 0);
    mpz_init_set_ui(QRA_Dcost, 0);
    mpz_init_set_ui(QRA_ancila, 0);
    mpz_init_set_ui(OPF_Gcost, 0);
    mpz_init_set_ui(OPF_Dcost, 0);
    mpz_init_set_ui(OPF_ancila, 0);
    mpz_init_set_ui(dif_Gcost, 0);
    mpz_init_set_ui(dif_Dcost, 0);
    mpz_init_set_ui(dif_ancila, 0);

    compute_QRA_cost(QRA_Gcost, QRA_Dcost, QRA_ancila, Vsize, nMn);
    compute_OPF_cost(OPF_Gcost, OPF_Dcost, OPF_ancila, Vsize, Msize);
    compute_dif_cost(dif_Gcost, dif_Dcost, dif_ancila, Vsize);

    mpz_set(QPGloop_Gcost, CPGloop_Gcost);
    mpz_add(QPGloop_Gcost, QPGloop_Gcost, QRA_Gcost);
    mpz_add(QPGloop_Gcost, QPGloop_Gcost, OPF_Gcost);
    mpz_add(QPGloop_Gcost, QPGloop_Gcost, dif_Gcost);

    mpz_set(QPGloop_Dcost, CPGloop_Dcost);
    if(mpz_cmp(QRA_Dcost, QPGloop_Dcost) > 0){
        mpz_set(QPGloop_Dcost, QRA_Dcost);
    }
    if(mpz_cmp(OPF_Dcost, QPGloop_Dcost) > 0){
        mpz_set(QPGloop_Dcost, OPF_Dcost);
    }
    if(mpz_cmp(dif_Dcost, QPGloop_Dcost) > 0){
        mpz_set(QPGloop_Dcost, dif_Dcost);
    }

    mpz_set(QPGloop_ancila, CPGloop_ancila);
    mpz_add(QPGloop_ancila, QPGloop_ancila, QRA_ancila);
    mpz_add(QPGloop_ancila, QPGloop_ancila, OPF_ancila);
    mpz_add(QPGloop_ancila, QPGloop_ancila, dif_ancila);
}


/*
 IQPrange_Gcost, QPGloop_Gcost, CPGloop_Gcost, grover_iteration, grover_processor_num, n, k を引数にとったとき，
 改善版量子Prangeアルゴリズム全体の G-cost を IQPrange_Gcost に格納する関数
*/
void compute_IQPrange_Gcost(mpz_t IQPrange_Gcost, mpz_t QPGloop_Gcost, mpz_t CPGloop_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k, mpz_t w){
    mpz_t QPG_Gcost;
    mpz_init_set(QPG_Gcost, QPGloop_Gcost);
    mpz_mul(QPG_Gcost, QPG_Gcost, grover_iteration);
    mpz_mul(QPG_Gcost, QPG_Gcost, grover_processor_num);

    mpz_t nSk;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);

    mpz_t n_bc_w, nSk_bc_w;
    mpz_init_set_ui(n_bc_w, 1);
    binomial_coefficient(n_bc_w, n, w);
    mpz_init_set_ui(nSk_bc_w, 1);
    binomial_coefficient(nSk_bc_w, nSk, w);

    mpz_t coefficient, log_coefficient;
    mpz_init(coefficient);
    mpz_cdiv_q(coefficient, n_bc_w, nSk_bc_w);
    mpz_init(log_coefficient);
    int logcoefficient = compute_log(coefficient);
    mpz_set_ui(log_coefficient, logcoefficient);

    mpz_add(QPG_Gcost, QPG_Gcost, log_coefficient);

    mpz_t matmul_Gcost, one;
    mpz_init_set_ui(matmul_Gcost, 0);
    mpz_init_set_ui(one, 1);
    compute_matmul_cost(matmul_Gcost, one, one, n, n, one);

    mpz_set(IQPrange_Gcost, QPG_Gcost);
    mpz_add(IQPrange_Gcost, IQPrange_Gcost, QPGloop_Gcost);
    mpz_add(IQPrange_Gcost, IQPrange_Gcost, matmul_Gcost);
}


/*
 IQPrange_Dcost, QPGloop_Dcost, CPGloop_Dcost, grover_iteration, grover_processor_num, n, k を引数にとったとき，
 改善版量子Prangeアルゴリズム全体の D-cost を IQPrange_Dcost に格納する関数
*/
void compute_IQPrange_Dcost(mpz_t IQPrange_Dcost, mpz_t QPGloop_Dcost, mpz_t CPGloop_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k, mpz_t w){
    mpz_t QPG_Dcost;
    mpz_init_set(QPG_Dcost, QPGloop_Dcost);
    mpz_mul(QPG_Dcost, QPG_Dcost, grover_iteration);
    mpz_cdiv_q(QPG_Dcost, QPG_Dcost, grover_processor_num);

    mpz_t matmul_Dcost, one;
    mpz_init_set_ui(matmul_Dcost, 0);
    mpz_init_set_ui(one, 1);
    compute_matmul_cost(one, matmul_Dcost, one, n, n, one);

    mpz_set(IQPrange_Dcost, QPG_Dcost);
    if(mpz_cmp(CPGloop_Dcost, IQPrange_Dcost) > 0){
        mpz_set(IQPrange_Dcost, CPGloop_Dcost);
    }
    if(mpz_cmp(matmul_Dcost, IQPrange_Dcost) > 0){
        mpz_set(IQPrange_Dcost, matmul_Dcost);
    }
}


/*
 IQPrange_Wcost, QPGloop_ancila, CPGloop_ancila, grover_iteration, grover_processor_num, n, k を引数にとったとき，
 改善版量子Prangeアルゴリズム全体の W-cost を IQPrange_Wcost に格納する関数
*/
void compute_IQPrange_Wcost(mpz_t IQPrange_Wcost, mpz_t QPGloop_ancila, mpz_t CPGloop_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k, mpz_t w){
    mpz_t nSk, QPrange_input;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init_set(QPrange_input, nSk);
    mpz_mul(QPrange_input, QPrange_input, n);
    mpz_addmul_ui(QPrange_input, n, 2);
    mpz_sub(QPrange_input, QPrange_input, k);

    mpz_t nMn, log_nMn;
    mpz_init(nMn);
    mpz_mul(nMn, n, n);
    int lognMn = compute_log(nMn);
    mpz_init(log_nMn);
    mpz_set_ui(log_nMn, lognMn);

    mpz_t n_bc_w, nSk_bc_w;
    mpz_init_set_ui(n_bc_w, 1);
    binomial_coefficient(n_bc_w, n, w);
    mpz_init_set_ui(nSk_bc_w, 1);
    binomial_coefficient(nSk_bc_w, nSk, w);

    mpz_t coefficient;
    mpz_init(coefficient);
    mpz_cdiv_q(coefficient, n_bc_w, nSk_bc_w);

    mpz_addmul(QPrange_input, coefficient, log_nMn);

    mpz_t QPG_ancila;
    mpz_init_set(QPG_ancila, QPGloop_ancila);
    mpz_mul(QPG_ancila, QPG_ancila, grover_iteration);
    mpz_mul(QPG_ancila, QPG_ancila, grover_processor_num);
    mpz_mul(QPG_ancila, QPG_ancila, grover_processor_num);

    mpz_t matmul_ancila, one;
    mpz_init_set_ui(matmul_ancila, 0);
    mpz_init_set_ui(one, 1);
    compute_matmul_cost(one, one, matmul_ancila, n, n, one);

    mpz_set(IQPrange_Wcost, QPrange_input);
    mpz_add(IQPrange_Wcost, IQPrange_Wcost, QPG_ancila);
    mpz_add(IQPrange_Wcost, IQPrange_Wcost, CPGloop_ancila);
    mpz_add(IQPrange_Wcost, IQPrange_Wcost, matmul_ancila);
}