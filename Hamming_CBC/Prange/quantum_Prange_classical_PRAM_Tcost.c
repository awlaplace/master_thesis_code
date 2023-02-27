/*
  gcc -c quantum_Prange_classical_PRAM_Tcost.c -lgmp && gcc -o quantum_Prange_classical_PRAM_Tcost ../../basic_operation.o ../Hamming_CBC_params.o Prange_params.o quantum_Prange_classical_PRAM_Tcost.o -lgmp && ./quantum_Prange_classical_PRAM_Tcost
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

            compute_QPGloop_PRAM_cost(QPGloop_Gcost, QPGloop_Dcost, QPGloop_ancila, CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k, w);

            // QPrange_Gcostを計算する
            compute_QPrange_PRAM_Gcost(QPrange_Gcost, QPGloop_Gcost, CPGloop_Gcost, grover_iteration, grover_processor_num, n, k);

            // CPrange_Dcost を計算する
            compute_QPrange_PRAM_Dcost(QPrange_Dcost, QPGloop_Dcost, CPGloop_Dcost, grover_iteration, grover_processor_num, n, k);
            
            // QPrange_Wcostを計算する
            compute_QPrange_PRAM_Wcost(QPrange_Wcost, QPGloop_ancila, CPGloop_ancila, grover_iteration, grover_processor_num, n, k);

            int grover_iteration_log = compute_log(grover_iteration);
            int QPrange_Gcost_log = compute_log(QPrange_Gcost);
            int QPrange_Dcost_log = compute_log(QPrange_Dcost);
            int QPrange_Wcost_log = compute_log(QPrange_Wcost);

            Grover_result_list[index][3] = QPrange_Gcost_log;
            Grover_result_list[index][4] = QPrange_Dcost_log;
            Grover_result_list[index][5] = QPrange_Wcost_log;
        }
        printf("n : %d, k : %d, w : %d, DWcost : %d\n", Grover_result_list[index][0], Grover_result_list[index][1], Grover_result_list[index][2], Grover_result_list[index][4] + Grover_result_list[index][5]);
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
 QPGloop_Gcost, QPGloop_Dcost, QPGloop_ancila, CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k, w を引数にとったとき，
 量子Prangeアルゴリズムの for 文のループ中の古典PRAM演算によるTゲートベースの G-cost, D-cost, アンシラビット数を
 それぞれ QPGloop_Gcost, QPGloop_Dcost, QPGloop_ancila に格納する関数
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

    compute_matmul_Tcost(matmul_Gcost1, matmul_Dcost1, matmul_ancila1, nSk, n, n);
    compute_GE_Tcost(GE_Gcost, GE_Dcost, GE_ancila, nSk, n);
    compute_matmul_Tcost(matmul_Gcost2, matmul_Dcost2, matmul_ancila2, nSk, nSk, one);
    compute_Ham_Tcost(Ham_Gcost, Ham_Dcost, Ham_ancila, nSk);

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
 量子Prangeアルゴリズムの for 文のループ中の古典PRAM演算によるTゲートベースの G-cost, D-cost, アンシラビット数を
 それぞれ QPGloop_Gcost, QPGloop_Dcost, QPGloop_ancila に格納する関数
*/
void compute_QPGloop_PRAM_cost(mpz_t QPGloop_Gcost, mpz_t QPGloop_Dcost, mpz_t QPGloop_ancila, mpz_t CPGloop_Gcost, mpz_t CPGloop_Dcost, mpz_t CPGloop_ancila, mpz_t n, mpz_t k, mpz_t w){
    mpz_set(QPGloop_Gcost, CPGloop_Gcost);

    mpz_set(QPGloop_Dcost, CPGloop_Dcost);

    mpz_set(QPGloop_ancila, CPGloop_ancila);
}


/*
 QPrange_Gcost, QPGloop_Gcost, CPGloop_Gcost, grover_iteration, grover_processor_num, n, k を引数にとったとき，
 量子Prangeアルゴリズムの全体の古典PRAM演算によるTゲートベースの G-cost を QPrange_Gcost に格納する関数
*/
void compute_QPrange_PRAM_Gcost(mpz_t QPrange_Gcost, mpz_t QPGloop_Gcost, mpz_t CPGloop_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k){
    mpz_t QPG_Gcost;
    mpz_init_set(QPG_Gcost, QPGloop_Gcost);
    mpz_mul(QPG_Gcost, QPG_Gcost, grover_iteration);

    mpz_t matmul_Gcost, one;
    mpz_init_set_ui(matmul_Gcost, 0);
    mpz_init_set_ui(one, 1);
    compute_matmul_Tcost(matmul_Gcost, one, one, n, n, one);

    mpz_set(QPrange_Gcost, QPG_Gcost);
    mpz_add(QPrange_Gcost, QPrange_Gcost, QPGloop_Gcost);
    mpz_add(QPrange_Gcost, QPrange_Gcost, matmul_Gcost);
}


/*
 QPrange_Dcost, QPGloop_Dcost, CPGloop_Dcost, grover_iteration, grover_processor_num, n, k を引数にとったとき，
 量子Prangeアルゴリズムの全体の古典PRAM演算によるTゲートベースの D-cost を QPrange_Gcost に格納する関数
*/
void compute_QPrange_PRAM_Dcost(mpz_t QPrange_Dcost, mpz_t QPGloop_Dcost, mpz_t CPGloop_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k){
    mpz_t QPG_Dcost;
    mpz_init_set(QPG_Dcost, QPGloop_Dcost);
    mpz_mul(QPG_Dcost, QPG_Dcost, grover_iteration);

    mpz_t matmul_Dcost, one;
    mpz_init_set_ui(matmul_Dcost, 0);
    mpz_init_set_ui(one, 1);
    compute_matmul_Tcost(one, matmul_Dcost, one, n, n, one);

    mpz_set(QPrange_Dcost, QPG_Dcost);
    if(mpz_cmp(CPGloop_Dcost, QPrange_Dcost) > 0){
        mpz_set(QPrange_Dcost, CPGloop_Dcost);
    }
    if(mpz_cmp(matmul_Dcost, QPrange_Dcost) > 0){
        mpz_set(QPrange_Dcost, matmul_Dcost);
    }
}


/*
 QPrange_Wcost, QPGloop_ancila, CPGloop_ancila, grover_iteration, grover_processor_num, n, k を引数にとったとき，
 量子Prangeアルゴリズムの全体の古典PRAM演算によるTゲートベースの W-cost を QPrange_Gcost に格納する関数
*/
void compute_QPrange_PRAM_Wcost(mpz_t QPrange_Wcost, mpz_t QPGloop_ancila, mpz_t CPGloop_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k){
    mpz_t nSk, QPrange_input;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init_set(QPrange_input, nSk);
    mpz_mul(QPrange_input, QPrange_input, n);
    mpz_addmul_ui(QPrange_input, n, 2);
    mpz_sub(QPrange_input, QPrange_input, k);

    mpz_t QPG_ancila;
    mpz_init_set(QPG_ancila, QPGloop_ancila);

    mpz_t matmul_ancila, one;
    mpz_init_set_ui(matmul_ancila, 0);
    mpz_init_set_ui(one, 1);
    compute_matmul_Tcost(one, one, matmul_ancila, n, n, one);

    mpz_set(QPrange_Wcost, QPrange_input);
    mpz_add(QPrange_Wcost, QPrange_Wcost, QPG_ancila);
    mpz_add(QPrange_Wcost, QPrange_Wcost, CPGloop_ancila);
    mpz_add(QPrange_Wcost, QPrange_Wcost, matmul_ancila);
}