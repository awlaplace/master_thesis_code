/*
  gcc -c improved_quantum_MMT_BJMM.c -lgmp && gcc -o improved_quantum_MMT_BJMM ../../basic_operation.o ../Hamming_CBC_params.o MMT_BJMM_params.o improved_quantum_MMT_BJMM.o -lgmp && ./improved_quantum_MMT_BJMM
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

        int min_IQMMTBJMM_Gcost_log = 0;
        int min_IQMMTBJMM_Dcost_log = MAX_DEPTH_log + 1;
        int min_IQMMTBJMM_Wcost_log = 0;

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

                        compute_IQMMTBJMMQWloop_cost(IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila, k, ell, p, epsilon);

                        compute_IQMMTBJMMQWloopout_cost(IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila, k, ell, p, epsilon);

                        compute_IQMMTBJMMGloop_cost(IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Dcost, IQMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);

                        // IQMMTBJMM_Gcost を計算する
                        compute_IQMMTBJMM_Gcost(IQMMTBJMM_Gcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                        
                        // IQMMTBJMM_Dcost を計算する
                        compute_IQMMTBJMM_Dcost(IQMMTBJMM_Dcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                        mpz_t cp_IQMMTBJMM_Dcost;
                        mpz_init_set(cp_IQMMTBJMM_Dcost, IQMMTBJMM_Dcost);
                        int cp_IQMMTBJMM_Dcost_log = compute_log(cp_IQMMTBJMM_Dcost);

                        // IQMMTBJMM_Wcost を計算する
                        compute_IQMMTBJMM_Wcost(IQMMTBJMM_Wcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);

                        if(mpz_cmp(IQMMTBJMM_Dcost, MAX_DEPTH) >= 0){
                            mpz_cdiv_q(grover_processor_num, IQMMTBJMM_Dcost, MAX_DEPTH);
                            compute_IQMMTBJMM_Dcost(IQMMTBJMM_Dcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                            mpz_t max_grover_processor_num;
                            mpz_init_set(max_grover_processor_num, grover_processor_num);
                            mpz_set_ui(IQMMTBJMM_Gcost, 0);
                            compute_IQMMTBJMM_Gcost(IQMMTBJMM_Gcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                            mpz_t min_IQMMTBJMM_Gcost, min_IQMMTBJMM_Dcost;
                            mpz_init_set(min_IQMMTBJMM_Gcost, IQMMTBJMM_Gcost);
                            mpz_init_set(min_IQMMTBJMM_Dcost, IQMMTBJMM_Dcost);
                            int QW_MMT_processor_cnt = 1;
                            while(QW_MMT_processor_cnt < cp_IQMMTBJMM_Dcost_log){
                                mpz_cdiv_q_ui(grover_processor_num, grover_processor_num, 2);
                                mpz_mul_ui(QW_MMT_processor_num, QW_MMT_processor_num, 2);
                                compute_IQMMTBJMM_Dcost(IQMMTBJMM_Dcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                                mpz_set_ui(IQMMTBJMM_Gcost, 0);
                                compute_IQMMTBJMM_Gcost(IQMMTBJMM_Gcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                                if((mpz_cmp(min_IQMMTBJMM_Gcost, IQMMTBJMM_Gcost) > 0) && (mpz_cmp(IQMMTBJMM_Dcost, MAX_DEPTH) < 0)){
                                    mpz_set(min_IQMMTBJMM_Gcost, IQMMTBJMM_Gcost);
                                    mpz_set(min_IQMMTBJMM_Dcost, IQMMTBJMM_Dcost);
                                }
                                QW_MMT_processor_cnt += 1;
                            }
                            mpz_set(IQMMTBJMM_Gcost, min_IQMMTBJMM_Gcost);
                            mpz_set(IQMMTBJMM_Dcost, min_IQMMTBJMM_Dcost);
                        }

                        int IQMMTBJMM_Gcost_loop_log = compute_log(IQMMTBJMM_Gcost);
                        int IQMMTBJMM_Dcost_loop_log = compute_log(IQMMTBJMM_Dcost);
                        int IQMMTBJMM_Wcost_loop_log = compute_log(IQMMTBJMM_Wcost);

                        if((min_IQMMTBJMM_Dcost_log == MAX_DEPTH_log + 1) || (IQMMTBJMM_Gcost_loop_log < min_IQMMTBJMM_Gcost_log)){
                            min_IQMMTBJMM_Gcost_log = IQMMTBJMM_Gcost_loop_log;
                            min_IQMMTBJMM_Dcost_log = IQMMTBJMM_Dcost_loop_log;
                            MMT_BJMM_result_list[index][3] = int_p * 4;
                            MMT_BJMM_result_list[index][4] = int_ell;
                            MMT_BJMM_result_list[index][5] = 0;
                            MMT_BJMM_result_list[index][6] = 0;
                            MMT_BJMM_result_list[index][7] = int_epsilon;
                            MMT_BJMM_result_list[index][8] = IQMMTBJMM_Gcost_loop_log;
                            MMT_BJMM_result_list[index][9] = IQMMTBJMM_Dcost_loop_log;
                            MMT_BJMM_result_list[index][10] = IQMMTBJMM_Wcost_loop_log;
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
 IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila, k, ell, p, epsilon を引数にとったとき，
 改善版量子MMT/BJMMアルゴリズムの量子ウォーク探索の for 文のループ中の G-cost, D-cost, アンシラビット数を
 それぞれ IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila に格納する関数
*/
void compute_IQMMTBJMMQWloop_cost(mpz_t IQMMTBJMMQWloop_Gcost, mpz_t IQMMTBJMMQWloop_Dcost, mpz_t IQMMTBJMMQWloop_ancila, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon){
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

    mpz_set(IQMMTBJMMQWloop_Gcost, OPF_Gcost);
    mpz_add(IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Gcost, dif_Gcost);

    mpz_set(IQMMTBJMMQWloop_Dcost, OPF_Dcost);
    if(mpz_cmp(dif_Dcost, IQMMTBJMMQWloop_Dcost) > 0){
        mpz_set(IQMMTBJMMQWloop_Dcost, dif_Dcost);
    }

    mpz_set(IQMMTBJMMQWloop_ancila, OPF_ancila);
    mpz_add(IQMMTBJMMQWloop_ancila, IQMMTBJMMQWloop_ancila, dif_ancila);
}


/*
 IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila, k, ell, p, epsilon を引数にとったとき，
 改善版量子MMT/BJMMアルゴリズムの量子ウォーク探索の for 文のループ外の G-cost, D-cost, アンシラビット数を
 それぞれ IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila に格納する関数
*/
void compute_IQMMTBJMMQWloopout_cost(mpz_t IQMMTBJMMQWloopout_Gcost, mpz_t IQMMTBJMMQWloopout_Dcost, mpz_t IQMMTBJMMQWloopout_ancila, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon){
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

    mpz_set(IQMMTBJMMQWloopout_Gcost, sp_Gcost);
    mpz_addmul_ui(IQMMTBJMMQWloopout_Gcost, QRA_Gcost, 4);

    mpz_set(IQMMTBJMMQWloopout_Dcost, sp_Dcost);
    if(mpz_cmp(QRA_Dcost, IQMMTBJMMQWloopout_Dcost) > 0){
        mpz_set(IQMMTBJMMQWloopout_Dcost, QRA_Dcost);
    }

    mpz_set(IQMMTBJMMQWloopout_ancila, sp_ancila);
    mpz_addmul_ui(IQMMTBJMMQWloopout_ancila, QRA_ancila, 4);
}


/*
 IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Dcost, IQMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon を引数にとったとき，
 改善版量子MMT/BJMMアルゴリズムの Grover の for 文のループ中の G-cost, D-cost, アンシラビット数を
 それぞれ IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Dcost, IQMMTBJMMGloop_ancila に格納する関数
*/
void compute_IQMMTBJMMGloop_cost(mpz_t IQMMTBJMMGloop_Gcost, mpz_t IQMMTBJMMGloop_Dcost, mpz_t IQMMTBJMMGloop_ancila, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t QRA_Gcost, QRA_Dcost, QRA_ancila;
    mpz_t CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila;
    mpz_t IQMMTBJMMQW_Gcost, IQMMTBJMMQW_Dcost, IQMMTBJMMQW_ancila;
    mpz_t IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila;
    mpz_t IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila;
    mpz_t OPF_Gcost, OPF_Dcost, OPF_ancila;
    mpz_t dif_Gcost, dif_Dcost, dif_ancila;
    mpz_init_set_ui(QRA_Gcost, 0);
    mpz_init_set_ui(QRA_Dcost, 0);
    mpz_init_set_ui(QRA_ancila, 0);
    mpz_init_set_ui(CPGloop_Gcost, 0);
    mpz_init_set_ui(CPGloop_Dcost, 0);
    mpz_init_set_ui(CPGloop_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQW_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQW_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQW_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_ancila, 0);
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

    mpz_t n_bc_w, kAell_bc_p, nSkSell_bc_wSp;
    mpz_init_set_ui(n_bc_w, 1);
    binomial_coefficient(n_bc_w, n, w);
    mpz_init_set_ui(kAell_bc_p, 1);
    binomial_coefficient(kAell_bc_p, kAell, p);
    mpz_init_set_ui(nSkSell_bc_wSp, 1);
    binomial_coefficient(nSkSell_bc_wSp, nSkSell, wSp);

    mpz_t Vsize, Msize;
    mpz_init_set(Vsize, n_bc_w);
    mpz_cdiv_q(Vsize, Vsize, kAell_bc_p);
    mpz_cdiv_q(Vsize, Vsize, nSkSell_bc_wSp);
    mpz_init_set_ui(Msize, 1);
    
    compute_QRA_cost(QRA_Gcost, QRA_Dcost, QRA_ancila, Vsize, Msize);
    compute_CPGloop_cost(CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k);
    compute_IQMMTBJMMQWloop_cost(IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila, k, ell, p, epsilon);
    compute_IQMMTBJMMQWloopout_cost(IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila, k, ell, p, epsilon);
    compute_OPF_cost(OPF_Gcost, OPF_Dcost, OPF_ancila, Vsize, Msize);
    compute_dif_cost(dif_Gcost, dif_Dcost, dif_ancila, Vsize);

    mpz_set(IQMMTBJMMQW_Gcost, IQMMTBJMMQWloop_Gcost);
    mpz_mul(IQMMTBJMMQW_Gcost, IQMMTBJMMQW_Gcost, QW_MMT_iteration);
    mpz_mul(IQMMTBJMMQW_Gcost, IQMMTBJMMQW_Gcost, QW_MMT_processor_num);
    mpz_add(IQMMTBJMMQW_Gcost, IQMMTBJMMQW_Gcost, IQMMTBJMMQWloopout_Gcost);

    mpz_set(IQMMTBJMMQW_Dcost, IQMMTBJMMQWloop_Dcost);
    mpz_mul(IQMMTBJMMQW_Dcost, IQMMTBJMMQW_Dcost, QW_MMT_iteration);
    mpz_cdiv_q(IQMMTBJMMQW_Dcost, IQMMTBJMMQW_Dcost, QW_MMT_processor_num);
    if(mpz_cmp(IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQW_Dcost) > 0){
        mpz_set(IQMMTBJMMQW_Dcost, IQMMTBJMMQWloopout_Dcost);
    }

    mpz_set(IQMMTBJMMQW_ancila, IQMMTBJMMQWloop_ancila);
    mpz_mul(IQMMTBJMMQW_ancila, IQMMTBJMMQW_ancila, QW_MMT_iteration);
    mpz_mul(IQMMTBJMMQW_ancila, IQMMTBJMMQW_ancila, QW_MMT_processor_num);
    mpz_mul(IQMMTBJMMQW_ancila, IQMMTBJMMQW_ancila, QW_MMT_processor_num);
    mpz_add(IQMMTBJMMQW_ancila, IQMMTBJMMQW_ancila, IQMMTBJMMQWloopout_ancila);

    mpz_set(IQMMTBJMMGloop_Gcost, QRA_Gcost);
    mpz_add(IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Gcost, CPGloop_Gcost);
    mpz_add(IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Gcost, IQMMTBJMMQW_Gcost);
    mpz_add(IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Gcost, OPF_Gcost);
    mpz_add(IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Gcost, dif_Gcost);

    mpz_set(IQMMTBJMMGloop_Dcost, QRA_Dcost);
    if(mpz_cmp(CPGloop_Dcost, IQMMTBJMMGloop_Dcost) > 0){
        mpz_set(IQMMTBJMMGloop_Dcost, CPGloop_Dcost);
    }
    if(mpz_cmp(IQMMTBJMMQW_Dcost, IQMMTBJMMGloop_Dcost) > 0){
        mpz_set(IQMMTBJMMGloop_Dcost, IQMMTBJMMQW_Dcost);
    }
    if(mpz_cmp(OPF_Dcost, IQMMTBJMMGloop_Dcost) > 0){
        mpz_set(IQMMTBJMMGloop_Dcost, OPF_Dcost);
    }
    if(mpz_cmp(dif_Dcost, IQMMTBJMMGloop_Dcost) > 0){
        mpz_set(IQMMTBJMMGloop_Dcost, dif_Dcost);
    }

    mpz_set(IQMMTBJMMGloop_ancila, QRA_ancila);
    mpz_add(IQMMTBJMMGloop_ancila, IQMMTBJMMGloop_ancila, CPGloop_ancila);
    mpz_add(IQMMTBJMMGloop_ancila, IQMMTBJMMGloop_ancila, IQMMTBJMMQW_ancila);
    mpz_add(IQMMTBJMMGloop_ancila, IQMMTBJMMGloop_ancila, OPF_ancila);
    mpz_add(IQMMTBJMMGloop_ancila, IQMMTBJMMGloop_ancila, dif_ancila);
}



/*
 IQMMTBJMM_Gcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon を引数にとったとき，
 改善版量子MMT/BJMMアルゴリズム全体の G-cost を IQMMTBJMM_Gcost に格納する関数
*/
void compute_IQMMTBJMM_Gcost(mpz_t IQMMTBJMM_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t IQMMTBJMMG_Gcost, IQMMTBJMMG_Dcost, IQMMTBJMMG_ancila;
    mpz_t IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Dcost, IQMMTBJMMGloop_ancila;
    mpz_t CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila;
    mpz_t matmul_Gcost, matmul_Dcost, matmul_ancila;
    mpz_t IQMMTBJMMQW_Gcost, IQMMTBJMMQW_Dcost, IQMMTBJMMQW_ancila;
    mpz_t IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila;
    mpz_t IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila;
    mpz_init_set_ui(IQMMTBJMMG_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMG_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMG_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMGloop_ancila, 0);
    mpz_init_set_ui(CPGloop_Gcost, 0);
    mpz_init_set_ui(CPGloop_Dcost, 0);
    mpz_init_set_ui(CPGloop_ancila, 0);
    mpz_init_set_ui(matmul_Gcost, 0);
    mpz_init_set_ui(matmul_Dcost, 0);
    mpz_init_set_ui(matmul_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQW_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQW_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQW_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_ancila, 0);

    mpz_t one;
    mpz_init_set_ui(one, 1);

    compute_IQMMTBJMMGloop_cost(IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Dcost, IQMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
    compute_CPGloop_cost(CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k);
    compute_matmul_cost(matmul_Gcost, matmul_Dcost, matmul_ancila, n, n, one);
    compute_IQMMTBJMMQWloop_cost(IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila, k, ell, p, epsilon);
    compute_IQMMTBJMMQWloopout_cost(IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila, k, ell, p, epsilon);
    
    mpz_t log_nSkSell_bc_wSp, kAell, nSkSell, wSp, n_bc_w, kAell_bc_p, nSkSell_bc_wSp;
    mpz_init(kAell);
    mpz_add(kAell, k, ell);
    mpz_init(nSkSell);
    mpz_sub(nSkSell, n, kAell);
    mpz_init(wSp);
    mpz_sub(wSp, w, p);
    mpz_init_set_ui(n_bc_w, 1);
    binomial_coefficient(n_bc_w, n, w);
    mpz_init_set_ui(kAell_bc_p, 1);
    binomial_coefficient(kAell_bc_p, kAell, p);
    mpz_init_set_ui(nSkSell_bc_wSp, 1);
    binomial_coefficient(nSkSell_bc_wSp, nSkSell, wSp);
    int lognSkSell_bc_wSp = compute_log(nSkSell_bc_wSp);
    mpz_init_set_ui(log_nSkSell_bc_wSp, lognSkSell_bc_wSp);

    mpz_set(IQMMTBJMMG_Gcost, IQMMTBJMMGloop_Gcost);
    mpz_mul(IQMMTBJMMG_Gcost, IQMMTBJMMG_Gcost, grover_iteration);
    mpz_mul(IQMMTBJMMG_Gcost, IQMMTBJMMG_Gcost, grover_processor_num);
    mpz_add(IQMMTBJMMG_Gcost, IQMMTBJMMG_Gcost, log_nSkSell_bc_wSp);

    mpz_set(IQMMTBJMMQW_Gcost, IQMMTBJMMQWloop_Gcost);
    mpz_mul(IQMMTBJMMQW_Gcost, IQMMTBJMMQW_Gcost, QW_MMT_iteration);
    mpz_mul(IQMMTBJMMQW_Gcost, IQMMTBJMMQW_Gcost, QW_MMT_processor_num);
    mpz_add(IQMMTBJMMQW_Gcost, IQMMTBJMMQW_Gcost, IQMMTBJMMQWloopout_Gcost);

    mpz_set(IQMMTBJMM_Gcost, IQMMTBJMMG_Gcost);
    mpz_add(IQMMTBJMM_Gcost, IQMMTBJMM_Gcost, CPGloop_Gcost);
    mpz_add(IQMMTBJMM_Gcost, IQMMTBJMM_Gcost, matmul_Gcost);
    mpz_add(IQMMTBJMM_Gcost, IQMMTBJMM_Gcost, IQMMTBJMMQW_Gcost);
}


/*
 IQMMTBJMM_Dcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon を引数にとったとき，
 改善版量子MMT/BJMMアルゴリズム全体の D-cost を IQMMTBJMM_Dcost に格納する関数
*/
void compute_IQMMTBJMM_Dcost(mpz_t IQMMTBJMM_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t IQMMTBJMMG_Gcost, IQMMTBJMMG_Dcost, IQMMTBJMMG_ancila;
    mpz_t IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Dcost, IQMMTBJMMGloop_ancila;
    mpz_t CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila;
    mpz_t matmul_Gcost, matmul_Dcost, matmul_ancila;
    mpz_t IQMMTBJMMQW_Gcost, IQMMTBJMMQW_Dcost, IQMMTBJMMQW_ancila;
    mpz_t IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila;
    mpz_t IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila;
    mpz_init_set_ui(IQMMTBJMMG_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMG_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMG_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMGloop_ancila, 0);
    mpz_init_set_ui(CPGloop_Gcost, 0);
    mpz_init_set_ui(CPGloop_Dcost, 0);
    mpz_init_set_ui(CPGloop_ancila, 0);
    mpz_init_set_ui(matmul_Gcost, 0);
    mpz_init_set_ui(matmul_Dcost, 0);
    mpz_init_set_ui(matmul_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQW_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQW_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQW_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_ancila, 0);

    mpz_t one;
    mpz_init_set_ui(one, 1);

    compute_IQMMTBJMMGloop_cost(IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Dcost, IQMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
    compute_CPGloop_cost(CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k);
    compute_matmul_cost(matmul_Gcost, matmul_Dcost, matmul_ancila, n, n, one);
    compute_IQMMTBJMMQWloop_cost(IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila, k, ell, p, epsilon);
    compute_IQMMTBJMMQWloopout_cost(IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila, k, ell, p, epsilon);

    mpz_set(IQMMTBJMMG_Dcost, IQMMTBJMMGloop_Dcost);
    mpz_mul(IQMMTBJMMG_Dcost, IQMMTBJMMG_Dcost, grover_iteration);
    mpz_cdiv_q(IQMMTBJMMG_Dcost, IQMMTBJMMG_Dcost, grover_processor_num);

    mpz_set(IQMMTBJMMQW_Dcost, IQMMTBJMMQWloop_Dcost);
    mpz_mul(IQMMTBJMMQW_Dcost, IQMMTBJMMQW_Dcost, QW_MMT_iteration);
    mpz_cdiv_q(IQMMTBJMMQW_Dcost, IQMMTBJMMQW_Dcost, QW_MMT_processor_num);
    if(mpz_cmp(IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQW_Dcost) > 0){
        mpz_set(IQMMTBJMMQW_Dcost, IQMMTBJMMQWloopout_Dcost);
    }

    mpz_set(IQMMTBJMM_Dcost, IQMMTBJMMG_Dcost);
    if(mpz_cmp(CPGloop_Dcost, IQMMTBJMM_Dcost) > 0){
        mpz_set(IQMMTBJMM_Dcost, CPGloop_Dcost);
    }
    if(mpz_cmp(matmul_Dcost, IQMMTBJMM_Dcost) > 0){
        mpz_set(IQMMTBJMM_Dcost, matmul_Dcost);
    }
    if(mpz_cmp(IQMMTBJMMQW_Dcost, IQMMTBJMM_Dcost) > 0){
        mpz_set(IQMMTBJMM_Dcost, IQMMTBJMMQW_Dcost);
    }
}


/*
 IQMMTBJMM_Wcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon を引数にとったとき，
 改善版量子MMT/BJMMアルゴリズム全体の W-cost を IQMMTBJMM_Wcost に格納する関数
*/
void compute_IQMMTBJMM_Wcost(mpz_t IQMMTBJMM_Wcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t IQMMTBJMMG_Gcost, IQMMTBJMMG_Dcost, IQMMTBJMMG_ancila;
    mpz_t IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Dcost, IQMMTBJMMGloop_ancila;
    mpz_t CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila;
    mpz_t matmul_Gcost, matmul_Dcost, matmul_ancila;
    mpz_t IQMMTBJMMQW_Gcost, IQMMTBJMMQW_Dcost, IQMMTBJMMQW_ancila;
    mpz_t IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila;
    mpz_t IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila;
    mpz_init_set_ui(IQMMTBJMMG_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMG_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMG_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMGloop_ancila, 0);
    mpz_init_set_ui(CPGloop_Gcost, 0);
    mpz_init_set_ui(CPGloop_Dcost, 0);
    mpz_init_set_ui(CPGloop_ancila, 0);
    mpz_init_set_ui(matmul_Gcost, 0);
    mpz_init_set_ui(matmul_Dcost, 0);
    mpz_init_set_ui(matmul_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQW_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQW_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQW_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_ancila, 0);

    mpz_t one;
    mpz_init_set_ui(one, 1);

    compute_IQMMTBJMMGloop_cost(IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Dcost, IQMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
    compute_CPGloop_cost(CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k);
    compute_matmul_cost(matmul_Gcost, matmul_Dcost, matmul_ancila, n, n, one);
    compute_IQMMTBJMMQWloop_cost(IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila, k, ell, p, epsilon);
    compute_IQMMTBJMMQWloopout_cost(IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila, k, ell, p, epsilon);
    
    mpz_t IQMMTBJMM_input, coefficient;
    mpz_init_set_ui(IQMMTBJMM_input, 0);
    mpz_init_set_ui(coefficient, 1);

    mpz_t nSk, kAell, kAellD4, pD4Aepsilon;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init(kAell);
    mpz_add(kAell, k, ell);
    mpz_init(kAellD4);
    mpz_cdiv_q_ui(kAellD4, kAell, 4);
    mpz_init(pD4Aepsilon);
    mpz_cdiv_q_ui(pD4Aepsilon, p, 4);
    mpz_add(pD4Aepsilon, pD4Aepsilon, epsilon);

    binomial_coefficient(coefficient, kAellD4, pD4Aepsilon);
    mpz_addmul_ui(IQMMTBJMM_input, coefficient, 4);

    mpz_t nSkSell, wSp, n_bc_w, kAell_bc_p, nSkSell_bc_wSp;
    mpz_init(nSkSell);
    mpz_sub(nSkSell, n, kAell);
    mpz_init(wSp);
    mpz_sub(wSp, w, p);
    mpz_init_set_ui(n_bc_w, 1);
    binomial_coefficient(n_bc_w, n, w);
    mpz_init_set_ui(kAell_bc_p, 1);
    binomial_coefficient(kAell_bc_p, kAell, p);
    mpz_init_set_ui(nSkSell_bc_wSp, 1);
    binomial_coefficient(nSkSell_bc_wSp, nSkSell, wSp);

    mpz_addmul(IQMMTBJMM_input, nSk, n);
    mpz_add(IQMMTBJMM_input, IQMMTBJMM_input, nSk);
    mpz_add(IQMMTBJMM_input, IQMMTBJMM_input, nSkSell_bc_wSp);

    mpz_add(IQMMTBJMM_input, IQMMTBJMM_input, n);

    mpz_set(IQMMTBJMMG_ancila, IQMMTBJMMGloop_ancila);
    mpz_mul(IQMMTBJMMG_ancila, IQMMTBJMMG_ancila, grover_iteration);
    mpz_mul(IQMMTBJMMG_ancila, IQMMTBJMMG_ancila, grover_processor_num);
    mpz_mul(IQMMTBJMMG_ancila, IQMMTBJMMG_ancila, grover_processor_num);

    mpz_set(IQMMTBJMMQW_ancila, IQMMTBJMMQWloop_ancila);
    mpz_mul(IQMMTBJMMQW_ancila, IQMMTBJMMQW_ancila, QW_MMT_iteration);
    mpz_mul(IQMMTBJMMQW_ancila, IQMMTBJMMQW_ancila, QW_MMT_processor_num);
    mpz_mul(IQMMTBJMMQW_ancila, IQMMTBJMMQW_ancila, QW_MMT_processor_num);
    mpz_add(IQMMTBJMMQW_ancila, IQMMTBJMMQW_ancila, IQMMTBJMMQWloopout_ancila);

    mpz_set(IQMMTBJMM_Wcost, IQMMTBJMM_input);
    mpz_add(IQMMTBJMM_Wcost, IQMMTBJMM_Wcost, IQMMTBJMMG_ancila);
    mpz_add(IQMMTBJMM_Wcost, IQMMTBJMM_Wcost, CPGloop_ancila);
    mpz_add(IQMMTBJMM_Wcost, IQMMTBJMM_Wcost, matmul_ancila);
    mpz_add(IQMMTBJMM_Wcost, IQMMTBJMM_Wcost, IQMMTBJMMQW_ancila);
}