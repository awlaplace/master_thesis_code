/*
  gcc -c quantum_MMT_BJMM.c -lgmp && gcc -o quantum_MMT_BJMM ../../basic_operation.o ../Hamming_CBC_params.o MMT_BJMM_params.o quantum_MMT_BJMM.o -lgmp && ./quantum_MMT_BJMM
*/

#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "MMT_BJMM.h"


int main(){
    init_Hamming_CBC_params();

    for(int index = 0; index < 9; index++){
        mpz_init_set_ui(n, SDP_instance_list[index][0]);
        mpz_init_set_ui(k, SDP_instance_list[index][1]);
        mpz_init_set_ui(w, SDP_instance_list[index][2]);

        mpz_init_set_ui(p, 4);
        mpz_init_set_ui(ell, 7);
        mpz_init_set_ui(ell1, 1);
        mpz_init_set_ui(ell2, 1);

        mpz_init_set_ui(grover_iteration, 1);
        mpz_init_set_ui(QW_MMT_iteration, 1);

        mpz_init_set_ui(grover_processor_num, 1);
        mpz_init_set_ui(QW_MMT_processor_num, 1);

        init_MMT_BJMM_computational_costs();

        int min_QMMTBJMM_Gcost_log = 0;
        int min_QMMTBJMM_Dcost_log = MAX_DEPTH_log + 1;
        int min_QMMTBJMM_Wcost_log = 0;

        for(int int_p = 1; int_p < SDP_instance_list[index][2] / 4; int_p++){
            mpz_set_ui(p, 4 * int_p);
            int int_wSp = SDP_instance_list[index][2] - int_p;
            mpz_set_ui(wSp, int_wSp);
            for(int int_ell = 1; int_ell < SDP_instance_list[index][0] - SDP_instance_list[index][1]; int_ell++){
                mpz_set_ui(ell, int_ell);
                int int_kAell = SDP_instance_list[index][1] + int_ell;
                mpz_set_ui(kAell, int_kAell);
                int int_nSkSell = SDP_instance_list[index][0] - SDP_instance_list[index][1] - int_ell;
                mpz_set_ui(nSkSell, int_nSkSell);
                mpz_cdiv_q_ui(kAellD2, kAell, 2);
                mpz_sub(kAellSp, kAell, p);
                for(int int_epsilon = 0; int_epsilon < int_p; int_epsilon++){
                    mpz_set_ui(epsilon, int_epsilon);
                    mpz_set_ui(pD4Aepsilon, int_p + int_epsilon);
                    mpz_mul_ui(c2epsilon, epsilon, 2);
                    // 二項係数が意味を持つ範囲に限定
                    if((mpz_cmp(n, w) > 0) && (mpz_cmp(kAell, p) > 0) && (mpz_cmp(nSkSell, wSp) > 0) && (mpz_cmp(kAellD2, pD4Aepsilon) > 0) && (mpz_cmp(kAellSp, c2epsilon) > 0)){
                        compute_grover_iteration(grover_iteration, n, k, w, ell, p);

                        compute_QW_MMT_iteration(QW_MMT_iteration, k, ell, p, epsilon);

                        compute_QMMTBJMMQWloop_cost(QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila, k, ell, p, epsilon);

                        compute_QMMTBJMMQWloopout_cost(QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWloopout_ancila, k, ell, p, epsilon);

                        compute_QMMTBJMMGloop_cost(QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Dcost, QMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);

                        // QMMTBJMM_Gcost を計算する
                        compute_QMMTBJMM_Gcost(QMMTBJMM_Gcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                        
                        // QMMTBJMM_Dcost を計算する
                        compute_QMMTBJMM_Dcost(QMMTBJMM_Dcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                        mpz_t cp_QMMTBJMM_Dcost;
                        mpz_init_set(cp_QMMTBJMM_Dcost, QMMTBJMM_Dcost);
                        int cp_QMMTBJMM_Dcost_log = compute_log(cp_QMMTBJMM_Dcost);

                        if(mpz_cmp(QMMTBJMM_Dcost, MAX_DEPTH) >= 0){
                            mpz_cdiv_q(grover_processor_num, QMMTBJMM_Dcost, MAX_DEPTH);
                            compute_QMMTBJMM_Dcost(QMMTBJMM_Dcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                            mpz_t max_grover_processor_num;
                            mpz_init_set(max_grover_processor_num, grover_processor_num);
                            mpz_set_ui(QMMTBJMM_Gcost, 0);
                            compute_QMMTBJMM_Gcost(QMMTBJMM_Gcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                            mpz_t min_QMMTBJMM_Gcost, min_QMMTBJMM_Dcost;
                            mpz_init_set(min_QMMTBJMM_Gcost, QMMTBJMM_Gcost);
                            mpz_init_set(min_QMMTBJMM_Dcost, QMMTBJMM_Dcost);
                            int QW_MMT_processor_cnt = 1;
                            while(QW_MMT_processor_cnt < cp_QMMTBJMM_Dcost_log){
                                mpz_cdiv_q_ui(grover_processor_num, grover_processor_num, 2);
                                mpz_mul_ui(QW_MMT_processor_num, QW_MMT_processor_num, 2);
                                compute_QMMTBJMM_Dcost(QMMTBJMM_Dcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                                mpz_set_ui(QMMTBJMM_Gcost, 0);
                                compute_QMMTBJMM_Gcost(QMMTBJMM_Gcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                                if((mpz_cmp(min_QMMTBJMM_Gcost, QMMTBJMM_Gcost) > 0) && (mpz_cmp(QMMTBJMM_Dcost, MAX_DEPTH) < 0)){
                                    mpz_set(min_QMMTBJMM_Gcost, QMMTBJMM_Gcost);
                                    mpz_set(min_QMMTBJMM_Dcost, QMMTBJMM_Dcost);
                                }
                                QW_MMT_processor_cnt += 1;
                            }
                            mpz_set(QMMTBJMM_Gcost, min_QMMTBJMM_Gcost);
                            mpz_set(QMMTBJMM_Dcost, min_QMMTBJMM_Dcost);
                        }

                        // QMMTBJMM_Wcost を計算する
                        compute_QMMTBJMM_Wcost(QMMTBJMM_Wcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);

                        int QMMTBJMM_Gcost_loop_log = compute_log(QMMTBJMM_Gcost);
                        int QMMTBJMM_Dcost_loop_log = compute_log(QMMTBJMM_Dcost);
                        int QMMTBJMM_Wcost_loop_log = compute_log(QMMTBJMM_Wcost);

                        if((min_QMMTBJMM_Dcost_log == MAX_DEPTH_log + 1) || (QMMTBJMM_Gcost_loop_log < min_QMMTBJMM_Gcost_log)){
                            min_QMMTBJMM_Gcost_log = QMMTBJMM_Gcost_loop_log;
                            min_QMMTBJMM_Dcost_log = QMMTBJMM_Dcost_loop_log;
                            MMT_BJMM_result_list[index][3] = int_p * 4;
                            MMT_BJMM_result_list[index][4] = int_ell;
                            MMT_BJMM_result_list[index][5] = 0;
                            MMT_BJMM_result_list[index][6] = 0;
                            MMT_BJMM_result_list[index][7] = int_epsilon;
                            MMT_BJMM_result_list[index][8] = QMMTBJMM_Gcost_loop_log;
                            MMT_BJMM_result_list[index][9] = QMMTBJMM_Dcost_loop_log;
                            MMT_BJMM_result_list[index][10] = QMMTBJMM_Wcost_loop_log;
                        }
                    }
                }
            }
        }
        printf("n : %d, k : %d, w : %d, p : %d, ell : %d, ell1 : %d, ell2 : %d, epsilon : %d, Gcost : %d, Dcost : %d, Wcost : %d\n", MMT_BJMM_result_list[index][0], MMT_BJMM_result_list[index][1], MMT_BJMM_result_list[index][2], MMT_BJMM_result_list[index][3], MMT_BJMM_result_list[index][4], MMT_BJMM_result_list[index][5], MMT_BJMM_result_list[index][6], MMT_BJMM_result_list[index][7], MMT_BJMM_result_list[index][8], MMT_BJMM_result_list[index][9], MMT_BJMM_result_list[index][10]);
    }
}


/*
 grover_iteration_cost, n, k, w, ell, p を引数にとったとき，
 Grover のループ回数を grover_itration_cost に格納する関数
*/
void compute_grover_iteration(mpz_t grover_iteration, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p){
    // the variable n_bc_w denotes nCw
    // the variable nSk_bc_w denotes n-kCw
    mpz_t nSk, nSkSell, wSp, kAell, n_bc_w, nSkSell_bc_wSp, kAell_bc_p;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init(nSkSell);
    mpz_sub(nSkSell, nSk, ell);
    mpz_init(wSp);
    mpz_sub(wSp, w, p);
    mpz_init(kAell);
    mpz_add(kAell, k, ell);
    mpz_init_set_ui(n_bc_w, 1);
    mpz_init_set_ui(nSkSell_bc_wSp, 1);
    mpz_init_set_ui(kAell_bc_p, 1);

    binomial_coefficient(n_bc_w, n, w);
    binomial_coefficient(nSkSell_bc_wSp, nSkSell, wSp);
    binomial_coefficient(kAell_bc_p, kAell, p);

    mpz_cdiv_q(grover_iteration, n_bc_w, nSkSell_bc_wSp);
    mpz_cdiv_q(grover_iteration, grover_iteration, kAell_bc_p);
    
    mpz_root(grover_iteration, grover_iteration, 2);
}


/*
 QW_MMT_iteration, k, ell, p, epsilon を引数にとったとき，
 MMT/BJMM のループ回数を QW_MMT_iteration に格納する関数
*/
void compute_QW_MMT_iteration(mpz_t QW_MMT_iteration, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t QW_MMT_iteration_denominator, QW_MMT_iteration_numerator, QW_MMT_iteration_numerator1, QW_MMT_iteration_numerator2;
    mpz_init_set_ui(QW_MMT_iteration_denominator, 1);
    mpz_init(QW_MMT_iteration_numerator);
    mpz_init_set_ui(QW_MMT_iteration_numerator1, 1);
    mpz_init_set_ui(QW_MMT_iteration_numerator2, 1);

    mpz_t kAellD2, pD4Aepsilon;
    mpz_init(kAellD2);
    mpz_init(pD4Aepsilon);
    mpz_add(kAellD2, k, ell);
    mpz_cdiv_q_ui(kAellD2, kAellD2, 4);
    mpz_cdiv_q_ui(pD4Aepsilon, p, 4);
    mpz_add(pD4Aepsilon, pD4Aepsilon, epsilon);

    binomial_coefficient(QW_MMT_iteration_denominator, kAellD2, pD4Aepsilon);

    mpz_root(QW_MMT_iteration_denominator, QW_MMT_iteration_denominator, 5);

    mpz_t cp_QW_MMT_iteration_denominator;
    mpz_init_set(cp_QW_MMT_iteration_denominator, QW_MMT_iteration_denominator);
    for(int i = 1; i < 6; i++){
        mpz_mul(QW_MMT_iteration_denominator, QW_MMT_iteration_denominator, cp_QW_MMT_iteration_denominator);
    }

    mpz_t pD2;
    mpz_init_set(pD2, p);
    mpz_cdiv_q_ui(pD2, pD2, 2);

    binomial_coefficient(QW_MMT_iteration_numerator1, p, pD2);

    mpz_t kAellSp, c2epsilon;
    mpz_init(kAellSp);
    mpz_init(c2epsilon);
    mpz_add(kAellSp, k, ell);
    mpz_sub(kAellSp, kAellSp, p);
    mpz_set(c2epsilon, epsilon);
    mpz_mul_ui(c2epsilon, c2epsilon, 2);

    binomial_coefficient(QW_MMT_iteration_numerator2, kAellSp, c2epsilon);

    mpz_mul(QW_MMT_iteration_numerator, QW_MMT_iteration_numerator1, QW_MMT_iteration_numerator2);

    mpz_sqrt(QW_MMT_iteration_numerator, QW_MMT_iteration_numerator);

    mpz_cdiv_q(QW_MMT_iteration, QW_MMT_iteration_denominator, QW_MMT_iteration_numerator);
}


/*
 CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k を引数にとったとき，
 古典Prangeアルゴリズムの while 文のループ回数の G-cost, D-cost, アンシラビット数を
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
 QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila, k, ell, p, epsilon を引数にとったとき，
 量子MMT/BJMMアルゴリズムの量子ウォーク探索の for 文のループ中の G-cost, D-cost, アンシラビット数を
 それぞれ QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila に格納する関数
*/
void compute_QMMTBJMMQWloop_cost(mpz_t QMMTBJMMQWloop_Gcost, mpz_t QMMTBJMMQWloop_Dcost, mpz_t QMMTBJMMQWloop_ancila, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t OPF_Gcost, OPF_Dcost, OPF_ancila;
    mpz_t dif_Gcost, dif_Dcost, dif_ancila;
    mpz_init_set_ui(OPF_Gcost, 0);
    mpz_init_set_ui(OPF_Dcost, 0);
    mpz_init_set_ui(OPF_ancila, 0);
    mpz_init_set_ui(dif_Gcost, 0);
    mpz_init_set_ui(dif_Dcost, 0);
    mpz_init_set_ui(dif_ancila, 0);
    
    mpz_t kAell, kAellD4, pD2, pD4Aepsilon, pM2;
    mpz_init(kAell);
    mpz_init(kAellD4);
    mpz_init(pD2);
    mpz_init(pD4Aepsilon);
    mpz_init(pM2);
    mpz_add(kAell, k, ell);
    mpz_cdiv_q_ui(kAellD4, kAell, 4);
    mpz_cdiv_q_ui(pD2, p, 2);
    mpz_cdiv_q_ui(pD4Aepsilon, p, 4);
    mpz_add(pD4Aepsilon, pD4Aepsilon, epsilon);
    mpz_mul_ui(pM2, p, 2);

    mpz_t Vsize, Msize;
    mpz_init_set_ui(Vsize, 1);
    binomial_coefficient(Vsize, kAellD4, pD4Aepsilon);
    mpz_init_set(Msize, Vsize);
    mpz_root(Msize, Msize, 7);
    mpz_mul(Msize, Msize, Msize);
    mpz_mul(Msize, Msize, Msize);

    mpz_t p_coefficient, cp_p_coefficient;
    mpz_init_set_ui(p_coefficient, 1);
    binomial_coefficient(p_coefficient, p, pD2);
    mpz_root(p_coefficient, p_coefficient, 7);
    mpz_init_set(cp_p_coefficient, p_coefficient);
    mpz_t exp_index;
    mpz_init_set_ui(exp_index, 1);
    while(mpz_cmp(exp_index, pM2) < 0){
        mpz_mul(p_coefficient, p_coefficient, cp_p_coefficient);
        mpz_add_ui(exp_index, exp_index, 1);
    }
    mpz_mul(Msize, Msize, p_coefficient);

    compute_OPF_cost(OPF_Gcost, OPF_Dcost, OPF_ancila, Vsize, Msize);
    compute_dif_cost(dif_Gcost, dif_Dcost, dif_ancila, Vsize);

    mpz_set(QMMTBJMMQWloop_Gcost, OPF_Gcost);
    mpz_add(QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Gcost, dif_Gcost);

    mpz_set(QMMTBJMMQWloop_Dcost, OPF_Dcost);
    if(mpz_cmp(dif_Dcost, QMMTBJMMQWloop_Dcost) > 0){
        mpz_set(QMMTBJMMQWloop_Dcost, dif_Dcost);
    }

    mpz_set(QMMTBJMMQWloop_ancila, OPF_ancila);
    mpz_add(QMMTBJMMQWloop_ancila, QMMTBJMMQWloop_ancila, dif_ancila);
}


/*
 QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWloopput_ancila, k, ell, p, epsilon を引数にとったとき，
 量子MMT/BJMMアルゴリズムの量子ウォーク探索の for 文のループ外の G-cost, D-cost, アンシラビット数を
 それぞれ QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWlooppout_ancila に格納する関数
*/
void compute_QMMTBJMMQWloopout_cost(mpz_t QMMTBJMMQWloopout_Gcost, mpz_t QMMTBJMMQWloopout_Dcost, mpz_t QMMTBJMMQWloopout_ancila, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t sp_Gcost, sp_Dcost, sp_ancila;
    mpz_t QRA_Gcost, QRA_Dcost, QRA_ancila;
    mpz_init_set_ui(sp_Gcost, 0);
    mpz_init_set_ui(sp_Dcost, 0);
    mpz_init_set_ui(sp_ancila, 0);
    mpz_init_set_ui(QRA_Gcost, 0);
    mpz_init_set_ui(QRA_Dcost, 0);
    mpz_init_set_ui(QRA_ancila, 0);

    mpz_t kAell, kAellD4, pD2, pD4Aepsilon, pM2;
    mpz_init(kAell);
    mpz_init(kAellD4);
    mpz_init(pD2);
    mpz_init(pD4Aepsilon);
    mpz_init(pM2);
    mpz_add(kAell, k, ell);
    mpz_cdiv_q_ui(kAellD4, kAell, 4);
    mpz_cdiv_q_ui(pD2, p, 2);
    mpz_cdiv_q_ui(pD4Aepsilon, p, 4);
    mpz_add(pD4Aepsilon, pD4Aepsilon, epsilon);
    mpz_mul_ui(pM2, p, 2);

    mpz_t Vsize, Msize;
    mpz_init_set_ui(Vsize, 1);
    binomial_coefficient(Vsize, kAellD4, pD4Aepsilon);
    mpz_init_set(Msize, Vsize);
    mpz_root(Msize, Msize, 7);
    mpz_mul(Msize, Msize, Msize);
    mpz_mul(Msize, Msize, Msize);

    mpz_t p_coefficient, cp_p_coefficient;
    mpz_init_set_ui(p_coefficient, 1);
    binomial_coefficient(p_coefficient, p, pD2);
    mpz_root(p_coefficient, p_coefficient, 7);
    mpz_init_set(cp_p_coefficient, p_coefficient);
    mpz_t exp_index;
    mpz_init_set_ui(exp_index, 1);
    while(mpz_cmp(exp_index, pM2) < 0){
        mpz_mul(p_coefficient, p_coefficient, cp_p_coefficient);
        mpz_add_ui(exp_index, exp_index, 1);
    }
    mpz_mul(Msize, Msize, p_coefficient);

    compute_sp_cost(sp_Gcost, sp_Dcost, sp_ancila, Vsize);
    compute_QRA_cost(QRA_Gcost, QRA_Dcost, QRA_ancila, Vsize, kAellD4);

    mpz_set(QMMTBJMMQWloopout_Gcost, sp_Gcost);
    mpz_addmul_ui(QMMTBJMMQWloopout_Gcost, QRA_Gcost, 4);

    mpz_set(QMMTBJMMQWloopout_Dcost, sp_Dcost);
    if(mpz_cmp(QRA_Dcost, QMMTBJMMQWloopout_Dcost) > 0){
        mpz_set(QMMTBJMMQWloopout_Dcost, QRA_Dcost);
    }

    mpz_set(QMMTBJMMQWloopout_ancila, sp_ancila);
    mpz_addmul_ui(QMMTBJMMQWloopout_ancila, QRA_ancila, 4);
}


/*
 QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Dcost, QMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon を引数にとったとき，
 量子MMT/BJMMアルゴリズムの Grover の for 文のループ中の G-cost, D-cost, アンシラビット数を
 それぞれ QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Dcost, QMMTBJMMGloop_ancila に格納する関数
*/
void compute_QMMTBJMMGloop_cost(mpz_t QMMTBJMMGloop_Gcost, mpz_t QMMTBJMMGloop_Dcost, mpz_t QMMTBJMMGloop_ancila, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t QRA_Gcost, QRA_Dcost, QRA_ancila;
    mpz_t CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila;
    mpz_t QMMTBJMMQW_Gcost, QMMTBJMMQW_Dcost, QMMTBJMMQW_ancila;
    mpz_t QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila;
    mpz_t QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWloopout_ancila;
    mpz_t OPF_Gcost, OPF_Dcost, OPF_ancila;
    mpz_t dif_Gcost, dif_Dcost, dif_ancila;
    mpz_init_set_ui(QRA_Gcost, 0);
    mpz_init_set_ui(QRA_Dcost, 0);
    mpz_init_set_ui(QRA_ancila, 0);
    mpz_init_set_ui(CPGloop_Gcost, 0);
    mpz_init_set_ui(CPGloop_Dcost, 0);
    mpz_init_set_ui(CPGloop_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQW_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQW_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQW_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_ancila, 0);
    mpz_init_set_ui(OPF_Gcost, 0);
    mpz_init_set_ui(OPF_Dcost, 0);
    mpz_init_set_ui(OPF_ancila, 0);
    mpz_init_set_ui(dif_Gcost, 0);
    mpz_init_set_ui(dif_Dcost, 0);
    mpz_init_set_ui(dif_ancila, 0);

    mpz_t nSw, kAell, nSkSell, wSp;
    mpz_init(nSw);
    mpz_sub(nSw, n, w);
    mpz_init(kAell);
    mpz_add(kAell, k, ell);
    mpz_init(nSkSell);
    mpz_sub(nSkSell, n, k);
    mpz_sub(nSkSell, nSkSell, ell);
    mpz_init(wSp);
    mpz_sub(wSp, w, p);

    mpz_t nfac, nSwfac, wfac, kAell_bc_p, nSkSell_bc_wSp;
    mpz_init_set_ui(nfac, 1);
    compute_fac(nfac, n);
    mpz_init_set_ui(nSwfac, 1);
    compute_fac(nSwfac, nSw);
    mpz_init_set_ui(wfac, 1);
    compute_fac(wfac, w);
    mpz_init_set_ui(kAell_bc_p, 1);
    binomial_coefficient(kAell_bc_p, kAell, p);
    mpz_init_set_ui(nSkSell_bc_wSp, 1);
    binomial_coefficient(nSkSell_bc_wSp, nSkSell, wSp);

    mpz_t Vsize, Msize;
    mpz_init_set(Vsize, nfac);
    mpz_init_set(Msize, nSwfac);
    mpz_mul(Msize, Msize, wfac);
    mpz_mul(Msize, Msize, kAell_bc_p);
    mpz_mul(Msize, Msize, nSkSell_bc_wSp);
    
    compute_QRA_cost(QRA_Gcost, QRA_Dcost, QRA_ancila, Vsize, Msize);
    compute_CPGloop_cost(CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k);
    compute_QMMTBJMMQWloop_cost(QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila, k, ell, p, epsilon);
    compute_QMMTBJMMQWloopout_cost(QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWloopout_ancila, k, ell, p, epsilon);
    compute_OPF_cost(OPF_Gcost, OPF_Dcost, OPF_ancila, Vsize, Msize);
    compute_dif_cost(dif_Gcost, dif_Dcost, dif_ancila, Vsize);

    mpz_set(QMMTBJMMQW_Gcost, QMMTBJMMQWloop_Gcost);
    mpz_mul(QMMTBJMMQW_Gcost, QMMTBJMMQW_Gcost, QW_MMT_iteration);
    mpz_mul(QMMTBJMMQW_Gcost, QMMTBJMMQW_Gcost, QW_MMT_processor_num);
    mpz_add(QMMTBJMMQW_Gcost, QMMTBJMMQW_Gcost, QMMTBJMMQWloopout_Gcost);

    mpz_set(QMMTBJMMQW_Dcost, QMMTBJMMQWloop_Dcost);
    mpz_mul(QMMTBJMMQW_Dcost, QMMTBJMMQW_Dcost, QW_MMT_iteration);
    mpz_cdiv_q(QMMTBJMMQW_Dcost, QMMTBJMMQW_Dcost, QW_MMT_processor_num);
    if(mpz_cmp(QMMTBJMMQWloopout_Dcost, QMMTBJMMQW_Dcost) > 0){
        mpz_set(QMMTBJMMQW_Dcost, QMMTBJMMQWloopout_Dcost);
    }

    mpz_set(QMMTBJMMQW_ancila, QMMTBJMMQWloop_ancila);
    mpz_mul(QMMTBJMMQW_ancila, QMMTBJMMQW_ancila, QW_MMT_iteration);
    mpz_mul(QMMTBJMMQW_ancila, QMMTBJMMQW_ancila, QW_MMT_processor_num);
    mpz_mul(QMMTBJMMQW_ancila, QMMTBJMMQW_ancila, QW_MMT_processor_num);
    mpz_add(QMMTBJMMQW_ancila, QMMTBJMMQW_ancila, QMMTBJMMQWloopout_ancila);

    mpz_set(QMMTBJMMGloop_Gcost, QRA_Gcost);
    mpz_add(QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Gcost, CPGloop_Gcost);
    mpz_add(QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Gcost, QMMTBJMMQW_Gcost);
    mpz_add(QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Gcost, OPF_Gcost);
    mpz_add(QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Gcost, dif_Gcost);

    mpz_set(QMMTBJMMGloop_Dcost, QRA_Dcost);
    if(mpz_cmp(CPGloop_Dcost, QMMTBJMMGloop_Dcost) > 0){
        mpz_set(QMMTBJMMGloop_Dcost, CPGloop_Dcost);
    }
    if(mpz_cmp(QMMTBJMMQW_Dcost, QMMTBJMMGloop_Dcost) > 0){
        mpz_set(QMMTBJMMGloop_Dcost, QMMTBJMMQW_Dcost);
    }
    if(mpz_cmp(OPF_Dcost, QMMTBJMMGloop_Dcost) > 0){
        mpz_set(QMMTBJMMGloop_Dcost, OPF_Dcost);
    }
    if(mpz_cmp(dif_Dcost, QMMTBJMMGloop_Dcost) > 0){
        mpz_set(QMMTBJMMGloop_Dcost, dif_Dcost);
    }

    mpz_set(QMMTBJMMGloop_ancila, QRA_ancila);
    mpz_add(QMMTBJMMGloop_ancila, QMMTBJMMGloop_ancila, CPGloop_ancila);
    mpz_add(QMMTBJMMGloop_ancila, QMMTBJMMGloop_ancila, QMMTBJMMQW_ancila);
    mpz_add(QMMTBJMMGloop_ancila, QMMTBJMMGloop_ancila, OPF_ancila);
    mpz_add(QMMTBJMMGloop_ancila, QMMTBJMMGloop_ancila, dif_ancila);
}


/*
 QMMTBJMM_Gcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon を引数にとったとき，
 量子MMT/BJMMアルゴリズム全体の G-cost を QMMTBJMM_Gcost に格納する関数
*/
void compute_QMMTBJMM_Gcost(mpz_t QMMTBJMM_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t QMMTBJMMG_Gcost, QMMTBJMMG_Dcost, QMMTBJMMG_ancila;
    mpz_t QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Dcost, QMMTBJMMGloop_ancila;
    mpz_t CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila;
    mpz_t matmul_Gcost, matmul_Dcost, matmul_ancila;
    mpz_t QMMTBJMMQW_Gcost, QMMTBJMMQW_Dcost, QMMTBJMMQW_ancila;
    mpz_t QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila;
    mpz_t QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWloopout_ancila;
    mpz_init_set_ui(QMMTBJMMG_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMG_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMG_ancila, 0);
    mpz_init_set_ui(QMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMGloop_ancila, 0);
    mpz_init_set_ui(CPGloop_Gcost, 0);
    mpz_init_set_ui(CPGloop_Dcost, 0);
    mpz_init_set_ui(CPGloop_ancila, 0);
    mpz_init_set_ui(matmul_Gcost, 0);
    mpz_init_set_ui(matmul_Dcost, 0);
    mpz_init_set_ui(matmul_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQW_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQW_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQW_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_ancila, 0);

    mpz_t one;
    mpz_init_set_ui(one, 1);

    compute_QMMTBJMMGloop_cost(QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Dcost, QMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
    compute_CPGloop_cost(CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k);
    compute_matmul_cost(matmul_Gcost, matmul_Dcost, matmul_ancila, n, n, one);
    compute_QMMTBJMMQWloop_cost(QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila, k, ell, p, epsilon);
    compute_QMMTBJMMQWloopout_cost(QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWloopout_ancila, k, ell, p, epsilon);
        
    mpz_t log_nfac, nfac;
    mpz_init_set_ui(nfac, 1);
    compute_fac(nfac, n);
    int lognfac = compute_log(nfac);
    mpz_init_set_ui(log_nfac, lognfac);

    mpz_set(QMMTBJMMG_Gcost, QMMTBJMMGloop_Gcost);
    mpz_mul(QMMTBJMMG_Gcost, QMMTBJMMG_Gcost, grover_iteration);
    mpz_mul(QMMTBJMMG_Gcost, QMMTBJMMG_Gcost, grover_processor_num);
    mpz_add(QMMTBJMMG_Gcost, QMMTBJMMG_Gcost, log_nfac);

    mpz_set(QMMTBJMMQW_Gcost, QMMTBJMMQWloop_Gcost);
    mpz_mul(QMMTBJMMQW_Gcost, QMMTBJMMQW_Gcost, QW_MMT_iteration);
    mpz_mul(QMMTBJMMQW_Gcost, QMMTBJMMQW_Gcost, QW_MMT_processor_num);
    mpz_add(QMMTBJMMQW_Gcost, QMMTBJMMQW_Gcost, QMMTBJMMQWloopout_Gcost);

    mpz_set(QMMTBJMM_Gcost, QMMTBJMMG_Gcost);
    mpz_add(QMMTBJMM_Gcost, QMMTBJMM_Gcost, CPGloop_Gcost);
    mpz_add(QMMTBJMM_Gcost, QMMTBJMM_Gcost, matmul_Gcost);
    mpz_add(QMMTBJMM_Gcost, QMMTBJMM_Gcost, QMMTBJMMQW_Gcost);
}


/*
 QMMTBJMM_Dcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon を引数にとったとき，
 量子MMT/BJMMアルゴリズム全体の D-cost を QMMTBJMM_Dcost に格納する関数
*/
void compute_QMMTBJMM_Dcost(mpz_t QMMTBJMM_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t QMMTBJMMG_Gcost, QMMTBJMMG_Dcost, QMMTBJMMG_ancila;
    mpz_t QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Dcost, QMMTBJMMGloop_ancila;
    mpz_t CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila;
    mpz_t matmul_Gcost, matmul_Dcost, matmul_ancila;
    mpz_t QMMTBJMMQW_Gcost, QMMTBJMMQW_Dcost, QMMTBJMMQW_ancila;
    mpz_t QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila;
    mpz_t QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWloopout_ancila;
    mpz_init_set_ui(QMMTBJMMG_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMG_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMG_ancila, 0);
    mpz_init_set_ui(QMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMGloop_ancila, 0);
    mpz_init_set_ui(CPGloop_Gcost, 0);
    mpz_init_set_ui(CPGloop_Dcost, 0);
    mpz_init_set_ui(CPGloop_ancila, 0);
    mpz_init_set_ui(matmul_Gcost, 0);
    mpz_init_set_ui(matmul_Dcost, 0);
    mpz_init_set_ui(matmul_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQW_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQW_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQW_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_ancila, 0);

    mpz_t one;
    mpz_init_set_ui(one, 1);

    compute_QMMTBJMMGloop_cost(QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Dcost, QMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
    compute_CPGloop_cost(CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k);
    compute_matmul_cost(matmul_Gcost, matmul_Dcost, matmul_ancila, n, n, one);
    compute_QMMTBJMMQWloop_cost(QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila, k, ell, p, epsilon);
    compute_QMMTBJMMQWloopout_cost(QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWloopout_ancila, k, ell, p, epsilon);

    mpz_set(QMMTBJMMG_Dcost, QMMTBJMMGloop_Dcost);
    mpz_mul(QMMTBJMMG_Dcost, QMMTBJMMG_Dcost, grover_iteration);
    mpz_cdiv_q(QMMTBJMMG_Dcost, QMMTBJMMG_Dcost, grover_processor_num);

    mpz_set(QMMTBJMMQW_Dcost, QMMTBJMMQWloop_Dcost);
    mpz_mul(QMMTBJMMQW_Dcost, QMMTBJMMQW_Dcost, QW_MMT_iteration);
    mpz_cdiv_q(QMMTBJMMQW_Dcost, QMMTBJMMQW_Dcost, QW_MMT_processor_num);
    if(mpz_cmp(QMMTBJMMQWloopout_Dcost, QMMTBJMMQW_Dcost) > 0){
        mpz_set(QMMTBJMMQW_Dcost, QMMTBJMMQWloopout_Dcost);
    }

    mpz_set(QMMTBJMM_Dcost, QMMTBJMMG_Dcost);
    if(mpz_cmp(CPGloop_Dcost, QMMTBJMM_Dcost) > 0){
        mpz_set(QMMTBJMM_Dcost, CPGloop_Dcost);
    }
    if(mpz_cmp(matmul_Dcost, QMMTBJMM_Dcost) > 0){
        mpz_set(QMMTBJMM_Dcost, matmul_Dcost);
    }
    if(mpz_cmp(QMMTBJMMQW_Dcost, QMMTBJMM_Dcost) > 0){
        mpz_set(QMMTBJMM_Dcost, QMMTBJMMQW_Dcost);
    }
}


/*
 QMMTBJMM_Wcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon を引数にとったとき，
 量子MMT/BJMMアルゴリズム全体の W-cost を QMMTBJMM_Wcost に格納する関数
*/
void compute_QMMTBJMM_Wcost(mpz_t QMMTBJMM_Wcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t QMMTBJMMG_Gcost, QMMTBJMMG_Dcost, QMMTBJMMG_ancila;
    mpz_t QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Dcost, QMMTBJMMGloop_ancila;
    mpz_t CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila;
    mpz_t matmul_Gcost, matmul_Dcost, matmul_ancila;
    mpz_t QMMTBJMMQW_Gcost, QMMTBJMMQW_Dcost, QMMTBJMMQW_ancila;
    mpz_t QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila;
    mpz_t QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWloopout_ancila;
    mpz_init_set_ui(QMMTBJMMG_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMG_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMG_ancila, 0);
    mpz_init_set_ui(QMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMGloop_ancila, 0);
    mpz_init_set_ui(CPGloop_Gcost, 0);
    mpz_init_set_ui(CPGloop_Dcost, 0);
    mpz_init_set_ui(CPGloop_ancila, 0);
    mpz_init_set_ui(matmul_Gcost, 0);
    mpz_init_set_ui(matmul_Dcost, 0);
    mpz_init_set_ui(matmul_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQW_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQW_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQW_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_ancila, 0);

    mpz_t one;
    mpz_init_set_ui(one, 1);

    compute_QMMTBJMMGloop_cost(QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Dcost, QMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
    compute_CPGloop_cost(CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k);
    compute_matmul_cost(matmul_Gcost, matmul_Dcost, matmul_ancila, n, n, one);
    compute_QMMTBJMMQWloop_cost(QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila, k, ell, p, epsilon);
    compute_QMMTBJMMQWloopout_cost(QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWloopout_ancila, k, ell, p, epsilon);
    
    mpz_t QMMTBJMM_input, coefficient;
    mpz_init_set_ui(QMMTBJMM_input, 0);
    mpz_init_set_ui(coefficient, 1);

    mpz_t kAell, kAellD4, pD4Aepsilon;
    mpz_init(kAell);
    mpz_add(kAell, k, ell);
    mpz_init(kAellD4);
    mpz_cdiv_q_ui(kAellD4, kAell, 4);
    mpz_init(pD4Aepsilon);
    mpz_cdiv_q_ui(pD4Aepsilon, p, 4);
    mpz_add(pD4Aepsilon, pD4Aepsilon, epsilon);

    binomial_coefficient(coefficient, kAellD4, pD4Aepsilon);
    mpz_addmul_ui(QMMTBJMM_input, coefficient, 4);

    mpz_t nSk, log_nfac, nfac;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init_set_ui(nfac, 1);
    compute_fac(nfac, n);
    int lognfac = compute_log(nfac);
    mpz_init_set_ui(log_nfac, lognfac);

    mpz_addmul(QMMTBJMM_input, nSk, n);
    mpz_add(QMMTBJMM_input, QMMTBJMM_input, nSk);
    mpz_add(QMMTBJMM_input, QMMTBJMM_input, log_nfac);

    mpz_add(QMMTBJMM_input, QMMTBJMM_input, n);

    mpz_set(QMMTBJMMG_ancila, QMMTBJMMGloop_ancila);
    mpz_mul(QMMTBJMMG_ancila, QMMTBJMMG_ancila, grover_iteration);
    mpz_mul(QMMTBJMMG_ancila, QMMTBJMMG_ancila, grover_processor_num);
    mpz_mul(QMMTBJMMG_ancila, QMMTBJMMG_ancila, grover_processor_num);

    mpz_set(QMMTBJMMQW_ancila, QMMTBJMMQWloop_ancila);
    mpz_mul(QMMTBJMMQW_ancila, QMMTBJMMQW_ancila, QW_MMT_iteration);
    mpz_mul(QMMTBJMMQW_ancila, QMMTBJMMQW_ancila, QW_MMT_processor_num);
    mpz_mul(QMMTBJMMQW_ancila, QMMTBJMMQW_ancila, QW_MMT_processor_num);
    mpz_add(QMMTBJMMQW_ancila, QMMTBJMMQW_ancila, QMMTBJMMQWloopout_ancila);

    mpz_set(QMMTBJMM_Wcost, QMMTBJMM_input);
    mpz_add(QMMTBJMM_Wcost, QMMTBJMM_Wcost, QMMTBJMMG_ancila);
    mpz_add(QMMTBJMM_Wcost, QMMTBJMM_Wcost, CPGloop_ancila);
    mpz_add(QMMTBJMM_Wcost, QMMTBJMM_Wcost, matmul_ancila);
    mpz_add(QMMTBJMM_Wcost, QMMTBJMM_Wcost, QMMTBJMMQW_ancila);
}