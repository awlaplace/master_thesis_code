/*
  gcc -c classical_Prange.c -lgmp && gcc -o classical_Prange ../../basic_operation.o ../Hamming_CBC_params.o Prange_params.o classical_Prange.o -lgmp && ./classical_Prange
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

            // CPrange_Dcost を計算する
            compute_CPrange_Dcost(CPrange_Dcost, CPGloop_Dcost, grover_iteration, grover_processor_num, n, k);

            // grover_processor_num の個数を調整
            while(mpz_cmp(CPrange_Dcost, MAX_DEPTH) > 0){
                mpz_mul_ui(grover_processor_num, grover_processor_num, 2);
                compute_CPrange_Dcost(CPrange_Dcost, CPGloop_Dcost, grover_iteration, grover_processor_num, n, k);
            }

            // CPrange_Gcostを計算する
            compute_CPrange_Gcost(CPrange_Gcost, CPGloop_Gcost, grover_iteration, grover_processor_num, n, k);
            
            // CPrange_Wcostを計算する
            compute_CPrange_Wcost(CPrange_Wcost, CPGloop_ancila, grover_iteration, grover_processor_num, n, k);

            int grover_iteration_log = compute_log(grover_iteration);
            int CPrange_Gcost_log = compute_log(CPrange_Gcost);
            int CPrange_Dcost_log = compute_log(CPrange_Dcost);
            int CPrange_Wcost_log = compute_log(CPrange_Wcost);

            Grover_result_list[index][3] = CPrange_Gcost_log;
            Grover_result_list[index][4] = CPrange_Dcost_log;
            Grover_result_list[index][5] = CPrange_Wcost_log;
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
}


/*
 CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k を引数にとったとき，
 古典Prangeアルゴリズムの for 文のループ回数の G-cost, D-cost, アンシラビット数を
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
 CPrange_Gcost, CPGloop_Gcost, grover_iteration, grover_processor_num, n, k を引数にとったとき，
 古典Prangeアルゴリズムの全体の G-cost を CPrange_Gcost に格納する関数
*/
void compute_CPrange_Gcost(mpz_t CPrange_Gcost, mpz_t CPGloop_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k){
    mpz_t nMn;
    mpz_init(nMn);
    mpz_mul(nMn, n, n);

    mpz_set(CPrange_Gcost, CPGloop_Gcost);
    mpz_mul(CPrange_Gcost, CPrange_Gcost, grover_iteration);
    mpz_mul(CPrange_Gcost, CPrange_Gcost, grover_processor_num);
    mpz_addmul_ui(CPrange_Gcost, nMn, 24);
}


/*
 CPrange_Dcost, CPGloop_Dcost, grover_iteration, grover_processor_num, n, k を引数にとったとき，
 古典Prangeアルゴリズムの全体の GDcost を CPrange_Dcost に格納する関数
*/
void compute_CPrange_Dcost(mpz_t CPrange_Dcost, mpz_t CPGloop_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k){
    mpz_t c16n;
    mpz_init(c16n);
    mpz_mul_ui(c16n, n, 16);
    
    mpz_set(CPrange_Dcost, CPGloop_Dcost);
    mpz_mul(CPrange_Dcost, CPrange_Dcost, grover_iteration);
    mpz_cdiv_q(CPrange_Dcost, CPrange_Dcost, grover_processor_num);
    if(mpz_cmp(c16n, CPrange_Dcost) > 0){
        mpz_set(CPrange_Dcost, c16n);
    }
}


/*
 CPrange_Wcost, CPGloop_ancila, grover_iteration, grover_processor_num, n, k を引数にとったとき，
 古典Prangeアルゴリズムの全体の W-cost を CPrange_Wcost に格納する関数
*/
void compute_CPrange_Wcost(mpz_t CPrange_Wcost, mpz_t CPGloop_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k){
    mpz_t nSk, CPrange_input;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init_set(CPrange_input, nSk);
    mpz_mul(CPrange_input, CPrange_input, n);
    mpz_addmul_ui(CPrange_input, n, 2);
    mpz_sub(CPrange_input, CPrange_input, k);

    mpz_set(CPrange_Wcost, CPGloop_ancila);
    mpz_mul(CPrange_Wcost, CPrange_Wcost, grover_iteration);
    mpz_mul(CPrange_Wcost, CPrange_Wcost, grover_processor_num);
    mpz_mul(CPrange_Wcost, CPrange_Wcost, grover_processor_num);
    mpz_add(CPrange_Wcost, CPrange_Wcost, CPrange_input);
    mpz_add(CPrange_Wcost, CPrange_Wcost, k);
}