/*
  gcc -c Dicke_quantum_MMT_BJMM.c -lgmp && gcc -o Dicke_quantum_MMT_BJMM ../../basic_operation.o ../Hamming_CBC_params.o MMT_BJMM_params.o Dicke_quantum_MMT_BJMM.o -lgmp && ./Dicke_quantum_MMT_BJMM
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

        mpz_init_set_ui(grover_iteration, 1);
        mpz_init_set_ui(QW_MMT_iteration, 1);

        mpz_init_set_ui(grover_processor_num, 1);
        mpz_init_set_ui(QW_MMT_processor_num, 1);

        init_MMT_BJMM_computational_costs();

        int min_DickeQMMTBJMM_Gcost_log = 0;
        int min_DickeQMMTBJMM_Dcost_log = MAX_DEPTH_log + 1;
        int min_DickeQMMTBJMM_Wcost_log = 0;

        for(int int_p = 3; int_p < SDP_instance_list[index][2] / 4; int_p++){
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

                        compute_DickeQMMTBJMMQWloop_cost(DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila, n, k, ell, p, epsilon);

                        compute_DickeQMMTBJMMGloop_cost(DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Dcost, DickeQMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                        
                        // DickeQMMTBJMM_Gcost を計算する
                        compute_DickeQMMTBJMM_Gcost(DickeQMMTBJMM_Gcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                        
                        // DickeQMMTBJMM_Dcost を計算する
                        compute_DickeQMMTBJMM_Dcost(DickeQMMTBJMM_Dcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);

                        mpz_t cp_DickeQMMTBJMM_Dcost;
                        mpz_init_set(cp_DickeQMMTBJMM_Dcost, DickeQMMTBJMM_Dcost);
                        int cp_DickeQMMTBJMM_Dcost_log = compute_log(cp_DickeQMMTBJMM_Dcost);
                        // DickeQMMTBJMM_Wcost を計算する
                        compute_DickeQMMTBJMM_Wcost(DickeQMMTBJMM_Wcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);

                        if(cp_DickeQMMTBJMM_Dcost_log > MAX_DEPTH_log){
                            mpz_cdiv_q(grover_processor_num, DickeQMMTBJMM_Dcost, MAX_DEPTH);
                            compute_DickeQMMTBJMM_Dcost(DickeQMMTBJMM_Dcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                            mpz_t max_grover_processor_num;
                            mpz_init_set(max_grover_processor_num, grover_processor_num);
                            mpz_set_ui(QMMTBJMM_Gcost, 0);
                            compute_DickeQMMTBJMM_Gcost(DickeQMMTBJMM_Gcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                            mpz_t min_DickeQMMTBJMM_Gcost, min_DickeQMMTBJMM_Dcost;
                            mpz_init_set(min_DickeQMMTBJMM_Gcost, DickeQMMTBJMM_Gcost);
                            mpz_init_set(min_DickeQMMTBJMM_Dcost, DickeQMMTBJMM_Dcost);
                            int QW_MMT_processor_cnt = 1;
                            while(QW_MMT_processor_cnt < cp_DickeQMMTBJMM_Dcost_log){
                                mpz_cdiv_q_ui(grover_processor_num, grover_processor_num, 2);
                                mpz_mul_ui(QW_MMT_processor_num, QW_MMT_processor_num, 2);
                                compute_DickeQMMTBJMM_Dcost(DickeQMMTBJMM_Dcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                                mpz_set_ui(DickeQMMTBJMM_Gcost, 0);
                                compute_DickeQMMTBJMM_Gcost(DickeQMMTBJMM_Gcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                                if((mpz_cmp(min_DickeQMMTBJMM_Gcost, DickeQMMTBJMM_Gcost) > 0) && (mpz_cmp(DickeQMMTBJMM_Dcost, MAX_DEPTH) < 0)){
                                    mpz_set(min_DickeQMMTBJMM_Gcost, DickeQMMTBJMM_Gcost);
                                    mpz_set(min_DickeQMMTBJMM_Dcost, DickeQMMTBJMM_Dcost);
                                    mpz_set_ui(DickeQMMTBJMM_Wcost, 0);
                                    compute_DickeQMMTBJMM_Wcost(DickeQMMTBJMM_Wcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
                                }
                                QW_MMT_processor_cnt += 1;
                            }
                            mpz_set(DickeQMMTBJMM_Gcost, min_DickeQMMTBJMM_Gcost);
                            mpz_set(DickeQMMTBJMM_Dcost, min_DickeQMMTBJMM_Dcost);
                        }

                        int DickeQMMTBJMM_Gcost_loop_log = compute_log(DickeQMMTBJMM_Gcost);
                        int DickeQMMTBJMM_Dcost_loop_log = compute_log(DickeQMMTBJMM_Dcost);
                        int DickeQMMTBJMM_Wcost_loop_log = compute_log(DickeQMMTBJMM_Wcost);

                        if((min_DickeQMMTBJMM_Dcost_log == MAX_DEPTH_log + 1) || (DickeQMMTBJMM_Gcost_loop_log < min_DickeQMMTBJMM_Gcost_log)){
                            min_DickeQMMTBJMM_Gcost_log = DickeQMMTBJMM_Gcost_loop_log;
                            min_DickeQMMTBJMM_Dcost_log = DickeQMMTBJMM_Dcost_loop_log;
                            MMT_BJMM_result_list[index][3] = int_p * 4;
                            MMT_BJMM_result_list[index][4] = int_ell;
                            MMT_BJMM_result_list[index][5] = 0;
                            MMT_BJMM_result_list[index][6] = 0;
                            MMT_BJMM_result_list[index][7] = int_epsilon;
                            MMT_BJMM_result_list[index][8] = DickeQMMTBJMM_Gcost_loop_log;
                            MMT_BJMM_result_list[index][9] = DickeQMMTBJMM_Dcost_loop_log;
                            MMT_BJMM_result_list[index][10] = DickeQMMTBJMM_Wcost_loop_log;
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
 G4SP_Gcost, G4SP_Dcost, G4SP_ancila, n, k, ell, p を引数にとったとき，
 G4SP の G-cost, D-cost, アンシラビット数を
 それぞれ G4SP_Gcost, G4SP_Dcost, G4SP_ancila に格納する関数
*/
void compute_G4SP_cost(mpz_t G4SP_Gcost, mpz_t G4SP_Dcost, mpz_t G4SP_ancila, mpz_t n, mpz_t k, mpz_t ell, mpz_t p){
    mpz_t add_Gcost1, add_Dcost1, add_ancila1;
    mpz_t matmul_Gcost, matmul_Dcost, matmul_ancila;
    mpz_t add_Gcost2, add_Dcost2, add_ancila2;
    mpz_init_set_ui(add_Gcost1, 0);
    mpz_init_set_ui(add_Dcost1, 0);
    mpz_init_set_ui(add_ancila1, 0);
    mpz_init_set_ui(matmul_Gcost, 0);
    mpz_init_set_ui(matmul_Dcost, 0);
    mpz_init_set_ui(matmul_ancila, 0);
    mpz_init_set_ui(add_Gcost2, 0);
    mpz_init_set_ui(add_Dcost2, 0);
    mpz_init_set_ui(add_ancila2, 0);

    mpz_t kAell, nSk, one;
    mpz_init(kAell);
    mpz_add(kAell, k, ell);
    mpz_init(nSk);
    mpz_add(nSk, n, k);
    mpz_init_set_ui(one, 1);

    compute_add_cost(add_Gcost1, add_Dcost1, add_ancila1, kAell);
    compute_matmul_cost(matmul_Gcost, matmul_Dcost, matmul_ancila, nSk, kAell, one);
    compute_add_cost(add_Gcost2, add_Dcost2, add_ancila2, nSk);

    mpz_mul_ui(add_Gcost1, add_Gcost1, 6);
    mpz_mul_ui(add_Dcost1, add_Dcost1, 6);
    mpz_mul_ui(add_ancila1, add_ancila1, 6);
    mpz_mul_ui(matmul_Gcost, matmul_Gcost, 2);
    mpz_mul_ui(matmul_Dcost, matmul_Dcost, 2);
    mpz_mul_ui(matmul_ancila, matmul_ancila, 2);
    mpz_addmul_ui(add_Gcost2, add_Gcost2, 2);
    mpz_addmul_ui(add_Dcost2, add_Dcost2, 2);
    mpz_addmul_ui(add_ancila2, add_ancila2, 2);

    mpz_set(G4SP_Gcost, add_Gcost1);
    mpz_add(G4SP_Gcost, G4SP_Gcost, matmul_Gcost);
    mpz_add(G4SP_Gcost, G4SP_Gcost, add_Gcost2);

    mpz_set(G4SP_Dcost, add_Dcost1);
    if(mpz_cmp(matmul_Dcost, G4SP_Dcost) > 0){
        mpz_set(G4SP_Dcost, matmul_Dcost);
    }
    if(mpz_cmp(add_Dcost2, G4SP_Dcost) > 0){
        mpz_set(G4SP_Dcost, add_Dcost2);
    }

    mpz_set(G4SP_ancila, add_ancila1);
    mpz_add(G4SP_ancila, G4SP_ancila, matmul_ancila);
    mpz_add(G4SP_ancila, G4SP_ancila, add_ancila2);
}


/*
 DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila, n, k, ell, p, epsilon を引数にとったとき，
 Dicke状態を用いた量子MMT/BJMMアルゴリズムの量子ウォーク探索の for 文のループ中の G-cost, D-cost, アンシラビット数を
 それぞれ DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila に格納する関数
*/
void compute_DickeQMMTBJMMQWloop_cost(mpz_t DickeQMMTBJMMQWloop_Gcost, mpz_t DickeQMMTBJMMQWloop_Dcost, mpz_t DickeQMMTBJMMQWloop_ancila, mpz_t n, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t OPFDicke_Gcost, OPFDicke_Dcost, OPFDicke_ancila;
    mpz_t dif_Gcost, dif_Dcost, dif_ancila;
    mpz_t G4SP_Gcost, G4SP_Dcost, G4SP_ancila;
    mpz_init_set_ui(OPFDicke_Gcost, 0);
    mpz_init_set_ui(OPFDicke_Dcost, 0);
    mpz_init_set_ui(OPFDicke_ancila, 0);
    mpz_init_set_ui(dif_Gcost, 0);
    mpz_init_set_ui(dif_Dcost, 0);
    mpz_init_set_ui(dif_ancila, 0);
    mpz_init_set_ui(G4SP_Gcost, 0);
    mpz_init_set_ui(G4SP_Dcost, 0);
    mpz_init_set_ui(G4SP_ancila, 0);

    mpz_t kAell, kAellD4, kAellMpD4, pD4Aepsilon;
    mpz_init(kAell);
    mpz_add(kAell, k, ell);
    mpz_init(kAellD4);
    mpz_cdiv_q_ui(kAellD4, kAell, 4);
    mpz_init(kAellMpD4);
    mpz_mul(kAellMpD4, kAell, p);
    mpz_cdiv_q_ui(kAellMpD4, kAellMpD4, 4);
    mpz_init(pD4Aepsilon);
    mpz_cdiv_q_ui(pD4Aepsilon, p, 4);
    mpz_add(pD4Aepsilon, pD4Aepsilon, epsilon);

    compute_OPFDicke_cost(OPFDicke_Gcost, OPFDicke_Dcost, OPFDicke_ancila, kAellMpD4);
    compute_Dicke_cost(dif_Gcost, dif_Dcost, dif_ancila, kAellD4, pD4Aepsilon);
    mpz_mul_ui(dif_Gcost, dif_Gcost, 2);
    mpz_add(dif_Gcost, dif_Gcost, n);
    mpz_mul_ui(dif_Dcost, dif_Dcost, 2);
    mpz_mul_ui(dif_ancila, dif_ancila, 2);
    compute_G4SP_cost(G4SP_Gcost, G4SP_Dcost, G4SP_ancila, n, k, ell, p);

    mpz_set(DickeQMMTBJMMQWloop_Gcost, OPFDicke_Gcost);
    mpz_add(DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Gcost, dif_Gcost);
    mpz_add(DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Gcost, G4SP_Gcost);

    mpz_set(DickeQMMTBJMMQWloop_Dcost, OPFDicke_Dcost);
    if(mpz_cmp(dif_Dcost, DickeQMMTBJMMQWloop_Dcost) > 0){
        mpz_set(DickeQMMTBJMMQWloop_Dcost, dif_Dcost);
    }
    if(mpz_cmp(G4SP_Dcost, DickeQMMTBJMMQWloop_Dcost) > 0){
        mpz_set(DickeQMMTBJMMQWloop_Dcost, G4SP_Dcost);
    }

    mpz_set(DickeQMMTBJMMQWloop_ancila, OPFDicke_ancila);
    mpz_add(DickeQMMTBJMMQWloop_ancila, DickeQMMTBJMMQWloop_ancila, dif_ancila);
    mpz_add(DickeQMMTBJMMQWloop_ancila, DickeQMMTBJMMQWloop_ancila, G4SP_ancila);
}


/*
 DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Dcost, DickeQMMTBJMMGloop_ancila, n, k, ell, p, epsilon を引数にとったとき，
 Dicke状態を用いた量子MMT/BJMMアルゴリズムの Grover の for 文のループ中の G-cost, D-cost, アンシラビット数を
 それぞれ DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Dcost, DickeQMMTBJMMGloop_ancila に格納する関数
*/
void compute_DickeQMMTBJMMGloop_cost(mpz_t DickeQMMTBJMMGloop_Gcost, mpz_t DickeQMMTBJMMGloop_Dcost, mpz_t DickeQMMTBJMMGloop_ancila, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t GE_Gcost, GE_Dcost, GE_ancila;
    mpz_t matmul_Gcost, matmul_Dcost, matmul_ancila;
    mpz_t DickeQMMTBJMMQW_Gcost, DickeQMMTBJMMQW_Dcost, DickeQMMTBJMMQW_ancila;
    mpz_t DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila;
    mpz_t Ham_Gcost, Ham_Dcost, Ham_ancila;
    mpz_t OPFDicke_Gcost, OPFDicke_Dcost, OPFDicke_ancila;
    mpz_t dif_Gcost, dif_Dcost, dif_ancila;
    mpz_init_set_ui(GE_Gcost, 0);
    mpz_init_set_ui(GE_Dcost, 0);
    mpz_init_set_ui(GE_ancila, 0);
    mpz_init_set_ui(matmul_Gcost, 0);
    mpz_init_set_ui(matmul_Dcost, 0);
    mpz_init_set_ui(matmul_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMQW_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQW_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQW_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(Ham_Gcost, 0);
    mpz_init_set_ui(Ham_Dcost, 0);
    mpz_init_set_ui(Ham_ancila, 0);
    mpz_init_set_ui(OPFDicke_Gcost, 0);
    mpz_init_set_ui(OPFDicke_Dcost, 0);
    mpz_init_set_ui(OPFDicke_ancila, 0);
    mpz_init_set_ui(dif_Gcost, 0);
    mpz_init_set_ui(dif_Dcost, 0);
    mpz_init_set_ui(dif_ancila, 0);

    mpz_t nSk, nSkMn, kAellD4, pD4Aepsilon, coefficient, one;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init(nSkMn);
    mpz_mul(nSkMn, nSk, n);
    mpz_init(kAellD4);
    mpz_add(kAellD4, k, ell);
    mpz_cdiv_q_ui(kAellD4, kAellD4, 4);
    mpz_init(pD4Aepsilon);
    mpz_cdiv_q_ui(pD4Aepsilon, p, 4);
    mpz_add(pD4Aepsilon, pD4Aepsilon, epsilon);
    mpz_init_set_ui(coefficient, 1);
    binomial_coefficient(coefficient, kAellD4, pD4Aepsilon);
    mpz_init_set_ui(one, 1);

    compute_GE_cost(GE_Gcost, GE_Dcost, GE_ancila, nSk, n);
    compute_matmul_cost(matmul_Gcost, matmul_Dcost, matmul_ancila, nSk, nSk, one);
    compute_DickeQMMTBJMMQWloop_cost(DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila, n, k, ell, p, epsilon);
    compute_Ham_cost(Ham_Gcost, Ham_Dcost, Ham_ancila, nSk);
    compute_OPFDicke_cost(OPFDicke_Gcost, OPFDicke_Dcost, OPFDicke_ancila, nSkMn);
    compute_Dicke_cost(dif_Gcost, dif_Dcost, dif_ancila, kAellD4, pD4Aepsilon);
    mpz_mul_ui(dif_Gcost, dif_Gcost, 2);
    mpz_add(dif_Gcost, dif_Gcost, n);
    mpz_mul_ui(dif_Dcost, dif_Dcost, 2);
    mpz_mul_ui(dif_ancila, dif_ancila, 2);

    mpz_set(DickeQMMTBJMMQW_Gcost, DickeQMMTBJMMQWloop_Gcost);
    mpz_mul(DickeQMMTBJMMQW_Gcost, DickeQMMTBJMMQW_Gcost, QW_MMT_iteration);
    mpz_mul(DickeQMMTBJMMQW_Gcost, DickeQMMTBJMMQW_Gcost, QW_MMT_processor_num);

    mpz_set(DickeQMMTBJMMQW_Dcost, DickeQMMTBJMMQWloop_Dcost);
    mpz_mul(DickeQMMTBJMMQW_Dcost, DickeQMMTBJMMQW_Dcost, QW_MMT_iteration);
    mpz_cdiv_q(DickeQMMTBJMMQW_Dcost, DickeQMMTBJMMQW_Dcost, QW_MMT_processor_num);

    mpz_set(DickeQMMTBJMMQW_ancila, DickeQMMTBJMMQWloop_ancila);
    mpz_mul(DickeQMMTBJMMQW_ancila, DickeQMMTBJMMQW_ancila, QW_MMT_iteration);
    mpz_mul(DickeQMMTBJMMQW_ancila, DickeQMMTBJMMQW_ancila, QW_MMT_processor_num);
    mpz_mul(DickeQMMTBJMMQW_ancila, DickeQMMTBJMMQW_ancila, QW_MMT_processor_num);

    mpz_set(DickeQMMTBJMMGloop_Gcost, GE_Gcost);
    mpz_add(DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Gcost, matmul_Gcost);
    mpz_add(DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMQW_Gcost);
    mpz_add(DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Gcost, Ham_Gcost);
    mpz_add(DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Gcost, OPFDicke_Gcost);
    mpz_add(DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Gcost, dif_Gcost);

    mpz_set(DickeQMMTBJMMGloop_Dcost, GE_Dcost);
    if(mpz_cmp(matmul_Dcost, DickeQMMTBJMMGloop_Dcost) > 0){
        mpz_set(DickeQMMTBJMMGloop_Dcost, matmul_Dcost);
    }
    if(mpz_cmp(DickeQMMTBJMMQW_Dcost, DickeQMMTBJMMGloop_Dcost) > 0){
        mpz_set(DickeQMMTBJMMGloop_Dcost, DickeQMMTBJMMQW_Dcost);
    }
    if(mpz_cmp(Ham_Dcost, DickeQMMTBJMMGloop_Dcost) > 0){
        mpz_set(DickeQMMTBJMMGloop_Dcost, Ham_Dcost);
    }
    if(mpz_cmp(OPFDicke_Dcost, DickeQMMTBJMMGloop_Dcost) > 0){
        mpz_set(DickeQMMTBJMMGloop_Dcost, OPFDicke_Dcost);
    }
    if(mpz_cmp(dif_Dcost, DickeQMMTBJMMGloop_Dcost) > 0){
        mpz_set(DickeQMMTBJMMGloop_Dcost, dif_Dcost);
    }

    mpz_set(DickeQMMTBJMMGloop_ancila, GE_ancila);
    mpz_add(DickeQMMTBJMMGloop_ancila, DickeQMMTBJMMGloop_ancila, matmul_ancila);
    mpz_add(DickeQMMTBJMMGloop_ancila, DickeQMMTBJMMGloop_ancila, coefficient);
    mpz_add(DickeQMMTBJMMGloop_ancila, DickeQMMTBJMMGloop_ancila, DickeQMMTBJMMQW_ancila);
    mpz_add(DickeQMMTBJMMGloop_ancila, DickeQMMTBJMMGloop_ancila, Ham_ancila);
    mpz_add(DickeQMMTBJMMGloop_ancila, DickeQMMTBJMMGloop_ancila, OPFDicke_ancila);
    mpz_add(DickeQMMTBJMMGloop_ancila, DickeQMMTBJMMGloop_ancila, dif_ancila);
}


/*
 DickeQMMTBJMM_Gcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon を引数にとったとき，
 Dicke状態を用いた量子MMT/BJMMアルゴリズム全体の G-cost を DickeQMMTBJMM_Gcost に格納する関数
*/
void compute_DickeQMMTBJMM_Gcost(mpz_t DickeQMMTBJMM_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t Dicke_Gcost, Dicke_Dcost, Dicke_ancila;
    mpz_t DickeQMMTBJMMG_Gcost, DickeQMMTBJMMG_Dcost, DickeQMMTBJMMG_ancila;
    mpz_t DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Dcost, DickeQMMTBJMMGloop_ancila;
    mpz_t GE_Gcost, GE_Dcost, GE_ancila;
    mpz_t matmul_Gcost1, matmul_Dcost1, matmul_ancila1;
    mpz_t DickeQMMTBJMMQW_Gcost, DickeQMMTBJMMQW_Dcost, DickeQMMTBJMMQW_ancila;
    mpz_t DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila;
    mpz_t Ham_Gcost, Ham_Dcost, Ham_ancila;
    mpz_t matmul_Gcost2, matmul_Dcost2, matmul_ancila2;
    mpz_t spDicke_Gcost, spDicke_Dcost, spDicke_ancila;
    mpz_init_set_ui(Dicke_Gcost, 0);
    mpz_init_set_ui(Dicke_Dcost, 0);
    mpz_init_set_ui(Dicke_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMG_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMG_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMG_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMGloop_ancila, 0);
    mpz_init_set_ui(GE_Gcost, 0);
    mpz_init_set_ui(GE_Dcost, 0);
    mpz_init_set_ui(GE_ancila, 0);
    mpz_init_set_ui(matmul_Gcost1, 0);
    mpz_init_set_ui(matmul_Dcost1, 0);
    mpz_init_set_ui(matmul_ancila1, 0);
    mpz_init_set_ui(DickeQMMTBJMMQW_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQW_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQW_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(Ham_Gcost, 0);
    mpz_init_set_ui(Ham_Dcost, 0);
    mpz_init_set_ui(Ham_ancila, 0);
    mpz_init_set_ui(matmul_Gcost2, 0);
    mpz_init_set_ui(matmul_Dcost2, 0);
    mpz_init_set_ui(matmul_ancila2, 0);
    mpz_init_set_ui(spDicke_Gcost, 0);
    mpz_init_set_ui(spDicke_Dcost, 0);
    mpz_init_set_ui(spDicke_ancila, 0);

    mpz_t nSk, nSkSell, kAellD4, pD4Aepsilon, one;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init_set(nSkSell, nSk);
    mpz_sub(nSkSell, nSkSell, ell);
    mpz_init(kAellD4);
    mpz_add(kAellD4, k, ell);
    mpz_cdiv_q_ui(kAellD4, kAellD4, 4);
    mpz_init(pD4Aepsilon);
    mpz_cdiv_q_ui(pD4Aepsilon, p, 4);
    mpz_add(pD4Aepsilon, pD4Aepsilon, epsilon);
    mpz_init_set_ui(one, 1);

    compute_Dicke_cost(Dicke_Gcost, Dicke_Dcost, Dicke_ancila, kAellD4, pD4Aepsilon);
    compute_DickeQMMTBJMMGloop_cost(DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Dcost, DickeQMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
    compute_GE_cost(GE_Gcost, GE_Dcost, GE_ancila, nSk, n);
    compute_matmul_cost(matmul_Gcost1, matmul_Dcost1, matmul_ancila1, nSk, nSk, one);
    compute_DickeQMMTBJMMQWloop_cost(DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila, n, k, ell, p, epsilon);
    compute_Ham_cost(Ham_Gcost, Ham_Dcost, Ham_ancila, nSk);
    compute_matmul_cost(matmul_Gcost2, matmul_Dcost2, matmul_ancila2, n, n, one);
    compute_spDicke_cost(spDicke_Gcost, spDicke_Dcost, spDicke_ancila, n, k, nSkSell);

    // DickeQMMTBJMMG_Gcost の計算
    mpz_set(DickeQMMTBJMMG_Gcost, DickeQMMTBJMMGloop_Gcost);
    mpz_mul(DickeQMMTBJMMG_Gcost, DickeQMMTBJMMG_Gcost, grover_iteration);
    mpz_mul(DickeQMMTBJMMG_Gcost, DickeQMMTBJMMG_Gcost, grover_processor_num);
    mpz_add(DickeQMMTBJMMG_Gcost, DickeQMMTBJMMG_Gcost, spDicke_Gcost);

    // DickeQMMTBJMMQW_Gcost の計算
    mpz_set(DickeQMMTBJMMQW_Gcost, DickeQMMTBJMMQWloop_Gcost);
    mpz_mul(DickeQMMTBJMMQW_Gcost, DickeQMMTBJMMQW_Gcost, QW_MMT_iteration);
    mpz_mul(DickeQMMTBJMMQW_Gcost, DickeQMMTBJMMQW_Gcost, QW_MMT_processor_num);

    // DickeQMMTBJMM_Gcost の計算
    mpz_set(DickeQMMTBJMM_Gcost, Dicke_Gcost);
    mpz_add(DickeQMMTBJMM_Gcost, DickeQMMTBJMM_Gcost, DickeQMMTBJMMG_Gcost);
    mpz_add(DickeQMMTBJMM_Gcost, DickeQMMTBJMM_Gcost, GE_Gcost);
    mpz_add(DickeQMMTBJMM_Gcost, DickeQMMTBJMM_Gcost, matmul_Gcost1);
    mpz_add(DickeQMMTBJMM_Gcost, DickeQMMTBJMM_Gcost, DickeQMMTBJMMQW_Gcost);
    mpz_add(DickeQMMTBJMM_Gcost, DickeQMMTBJMM_Gcost, Ham_Gcost);
    mpz_add(DickeQMMTBJMM_Gcost, DickeQMMTBJMM_Gcost, matmul_Gcost2);
}


/*
 DickeQMMTBJMM_Dcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon を引数にとったとき，
 Dicke状態を用いた量子MMT/BJMMアルゴリズム全体の D-cost を DickeQMMTBJMM_Dcost に格納する関数
*/
void compute_DickeQMMTBJMM_Dcost(mpz_t DickeQMMTBJMM_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t Dicke_Gcost, Dicke_Dcost, Dicke_ancila;
    mpz_t DickeQMMTBJMMG_Gcost, DickeQMMTBJMMG_Dcost, DickeQMMTBJMMG_ancila;
    mpz_t DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Dcost, DickeQMMTBJMMGloop_ancila;
    mpz_t GE_Gcost, GE_Dcost, GE_ancila;
    mpz_t matmul_Gcost1, matmul_Dcost1, matmul_ancila1;
    mpz_t DickeQMMTBJMMQW_Gcost, DickeQMMTBJMMQW_Dcost, DickeQMMTBJMMQW_ancila;
    mpz_t DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila;
    mpz_t Ham_Gcost, Ham_Dcost, Ham_ancila;
    mpz_t matmul_Gcost2, matmul_Dcost2, matmul_ancila2;
    mpz_t spDicke_Gcost, spDicke_Dcost, spDicke_ancila;
    mpz_init_set_ui(Dicke_Gcost, 0);
    mpz_init_set_ui(Dicke_Dcost, 0);
    mpz_init_set_ui(Dicke_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMG_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMG_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMG_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMGloop_ancila, 0);
    mpz_init_set_ui(GE_Gcost, 0);
    mpz_init_set_ui(GE_Dcost, 0);
    mpz_init_set_ui(GE_ancila, 0);
    mpz_init_set_ui(matmul_Gcost1, 0);
    mpz_init_set_ui(matmul_Dcost1, 0);
    mpz_init_set_ui(matmul_ancila1, 0);
    mpz_init_set_ui(DickeQMMTBJMMQW_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQW_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQW_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(Ham_Gcost, 0);
    mpz_init_set_ui(Ham_Dcost, 0);
    mpz_init_set_ui(Ham_ancila, 0);
    mpz_init_set_ui(matmul_Gcost2, 0);
    mpz_init_set_ui(matmul_Dcost2, 0);
    mpz_init_set_ui(matmul_ancila2, 0);
    mpz_init_set_ui(spDicke_Gcost, 0);
    mpz_init_set_ui(spDicke_Dcost, 0);
    mpz_init_set_ui(spDicke_ancila, 0);

    mpz_t nSk, nSkSell, kAellD4, pD4Aepsilon, one;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init_set(nSkSell, nSk);
    mpz_sub(nSkSell, nSkSell, ell);
    mpz_init(kAellD4);
    mpz_add(kAellD4, k, ell);
    mpz_cdiv_q_ui(kAellD4, kAellD4, 4);
    mpz_init(pD4Aepsilon);
    mpz_cdiv_q_ui(pD4Aepsilon, p, 4);
    mpz_add(pD4Aepsilon, pD4Aepsilon, epsilon);
    mpz_init_set_ui(one, 1);

    compute_Dicke_cost(Dicke_Gcost, Dicke_Dcost, Dicke_ancila, kAellD4, pD4Aepsilon);
    compute_DickeQMMTBJMMGloop_cost(DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Dcost, DickeQMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
    compute_GE_cost(GE_Gcost, GE_Dcost, GE_ancila, nSk, n);
    compute_matmul_cost(matmul_Gcost1, matmul_Dcost1, matmul_ancila1, nSk, nSk, one);
    compute_DickeQMMTBJMMQWloop_cost(DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila, n, k, ell, p, epsilon);
    compute_Ham_cost(Ham_Gcost, Ham_Dcost, Ham_ancila, nSk);
    compute_matmul_cost(matmul_Gcost2, matmul_Dcost2, matmul_ancila2, n, n, one);
    compute_spDicke_cost(spDicke_Gcost, spDicke_Dcost, spDicke_ancila, n, k, nSkSell);

    // DickeQMMTBJMMG_Dcost の計算
    mpz_set(DickeQMMTBJMMG_Dcost, DickeQMMTBJMMGloop_Dcost);
    mpz_mul(DickeQMMTBJMMG_Dcost, DickeQMMTBJMMG_Dcost, grover_iteration);
    mpz_cdiv_q(DickeQMMTBJMMG_Dcost, DickeQMMTBJMMG_Dcost, grover_processor_num);
    if(mpz_cmp(spDicke_Dcost, DickeQMMTBJMMG_Dcost) > 0){
        mpz_set(DickeQMMTBJMMG_Dcost, spDicke_Dcost);
    }

    // DickeQMMTBJMMQW_Dcost の計算
    mpz_set(DickeQMMTBJMMQW_Dcost, DickeQMMTBJMMQWloop_Dcost);
    mpz_mul(DickeQMMTBJMMQW_Dcost, DickeQMMTBJMMQW_Dcost, QW_MMT_iteration);
    mpz_cdiv_q(DickeQMMTBJMMQW_Dcost, DickeQMMTBJMMQW_Dcost, QW_MMT_processor_num);

    // DickeQMMTBJMM_Dcost の計算
    mpz_set(DickeQMMTBJMM_Dcost, Dicke_Dcost);
    if(mpz_cmp(Dicke_Dcost, DickeQMMTBJMM_Dcost) > 0){
        mpz_set(DickeQMMTBJMM_Dcost, Dicke_Dcost);
    }
    if(mpz_cmp(DickeQMMTBJMMG_Dcost, DickeQMMTBJMM_Dcost) > 0){
        mpz_set(DickeQMMTBJMM_Dcost, DickeQMMTBJMMG_Dcost);
    }
    if(mpz_cmp(GE_Dcost, DickeQMMTBJMM_Dcost) > 0){
        mpz_set(DickeQMMTBJMM_Dcost, GE_Dcost);
    }
    if(mpz_cmp(matmul_Dcost1, DickeQMMTBJMM_Dcost) > 0){
        mpz_set(DickeQMMTBJMM_Dcost, matmul_Dcost1);
    }
    if(mpz_cmp(DickeQMMTBJMMQW_Dcost, DickeQMMTBJMM_Dcost) > 0){
        mpz_set(DickeQMMTBJMM_Dcost, DickeQMMTBJMMQW_Dcost);
    }
    if(mpz_cmp(Ham_Dcost, DickeQMMTBJMM_Dcost) > 0){
        mpz_set(DickeQMMTBJMM_Dcost, Ham_Dcost);
    }
    if(mpz_cmp(matmul_Dcost2, DickeQMMTBJMM_Dcost) > 0){
        mpz_set(DickeQMMTBJMM_Dcost, matmul_Dcost2);
    }
}


/*
 DickeQMMTBJMM_Wcost, grover_iteration, grover_processor_num, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon を引数にとったとき，
 Dicke状態を用いた量子MMT/BJMMアルゴリズム全体の W-cost を DickeQMMTBJMM_Wcost に格納する関数
*/
void compute_DickeQMMTBJMM_Wcost(mpz_t DickeQMMTBJMM_Wcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon){
    mpz_t Dicke_Gcost, Dicke_Dcost, Dicke_ancila;
    mpz_t DickeQMMTBJMMG_Gcost, DickeQMMTBJMMG_Dcost, DickeQMMTBJMMG_ancila;
    mpz_t DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Dcost, DickeQMMTBJMMGloop_ancila;
    mpz_t GE_Gcost, GE_Dcost, GE_ancila;
    mpz_t matmul_Gcost1, matmul_Dcost1, matmul_ancila1;
    mpz_t DickeQMMTBJMMQW_Gcost, DickeQMMTBJMMQW_Dcost, DickeQMMTBJMMQW_ancila;
    mpz_t DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila;
    mpz_t Ham_Gcost, Ham_Dcost, Ham_ancila;
    mpz_t matmul_Gcost2, matmul_Dcost2, matmul_ancila2;
    mpz_t spDicke_Gcost, spDicke_Dcost, spDicke_ancila;
    mpz_init_set_ui(Dicke_Gcost, 0);
    mpz_init_set_ui(Dicke_Dcost, 0);
    mpz_init_set_ui(Dicke_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMG_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMG_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMG_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMGloop_ancila, 0);
    mpz_init_set_ui(GE_Gcost, 0);
    mpz_init_set_ui(GE_Dcost, 0);
    mpz_init_set_ui(GE_ancila, 0);
    mpz_init_set_ui(matmul_Gcost1, 0);
    mpz_init_set_ui(matmul_Dcost1, 0);
    mpz_init_set_ui(matmul_ancila1, 0);
    mpz_init_set_ui(DickeQMMTBJMMQW_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQW_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQW_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(Ham_Gcost, 0);
    mpz_init_set_ui(Ham_Dcost, 0);
    mpz_init_set_ui(Ham_ancila, 0);
    mpz_init_set_ui(matmul_Gcost2, 0);
    mpz_init_set_ui(matmul_Dcost2, 0);
    mpz_init_set_ui(matmul_ancila2, 0);
    mpz_init_set_ui(spDicke_Gcost, 0);
    mpz_init_set_ui(spDicke_Dcost, 0);
    mpz_init_set_ui(spDicke_ancila, 0);

    mpz_t nSk, nSkSell, kAellD4, pD4Aepsilon, coefficient, one;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init_set(nSkSell, nSk);
    mpz_sub(nSkSell, nSkSell, ell);
    mpz_init(kAellD4);
    mpz_add(kAellD4, k, ell);
    mpz_cdiv_q_ui(kAellD4, kAellD4, 4);
    mpz_init(pD4Aepsilon);
    mpz_cdiv_q_ui(pD4Aepsilon, p, 4);
    mpz_add(pD4Aepsilon, pD4Aepsilon, epsilon);
    mpz_init_set_ui(coefficient, 1);
    binomial_coefficient(coefficient, kAellD4, pD4Aepsilon);
    mpz_init_set_ui(one, 1);

    compute_Dicke_cost(Dicke_Gcost, Dicke_Dcost, Dicke_ancila, kAellD4, pD4Aepsilon);
    compute_DickeQMMTBJMMGloop_cost(DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Dcost, DickeQMMTBJMMGloop_ancila, QW_MMT_iteration, QW_MMT_processor_num, n, k, w, ell, p, epsilon);
    compute_GE_cost(GE_Gcost, GE_Dcost, GE_ancila, nSk, n);
    compute_matmul_cost(matmul_Gcost1, matmul_Dcost1, matmul_ancila1, nSk, nSk, one);
    compute_DickeQMMTBJMMQWloop_cost(DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila, n, k, ell, p, epsilon);
    compute_Ham_cost(Ham_Gcost, Ham_Dcost, Ham_ancila, nSk);
    compute_matmul_cost(matmul_Gcost2, matmul_Dcost2, matmul_ancila2, n, n, one);
    compute_spDicke_cost(spDicke_Gcost, spDicke_Dcost, spDicke_ancila, n, k, nSkSell);

    // DickeQMMTBJMMG_ancila の計算
    mpz_set(DickeQMMTBJMMG_ancila, DickeQMMTBJMMGloop_ancila);
    mpz_mul(DickeQMMTBJMMG_ancila, DickeQMMTBJMMG_ancila, grover_iteration);
    mpz_mul(DickeQMMTBJMMG_ancila, DickeQMMTBJMMG_ancila, grover_processor_num);
    mpz_mul(DickeQMMTBJMMG_ancila, DickeQMMTBJMMG_ancila, grover_processor_num);
    mpz_add(DickeQMMTBJMMG_ancila, DickeQMMTBJMMG_ancila, spDicke_Gcost);

    // DickeQMMTBJMMQW_ancila の計算
    mpz_set(DickeQMMTBJMMQW_ancila, DickeQMMTBJMMQWloop_ancila);
    mpz_mul(DickeQMMTBJMMQW_ancila, DickeQMMTBJMMQW_ancila, QW_MMT_iteration);
    mpz_mul(DickeQMMTBJMMQW_ancila, DickeQMMTBJMMQW_ancila, QW_MMT_processor_num);
    mpz_mul(DickeQMMTBJMMQW_ancila, DickeQMMTBJMMQW_ancila, QW_MMT_processor_num);

    // DickeQMMTBJMM_input の計算
    mpz_t DickeQMMTBJMM_input;
    mpz_init_set_ui(DickeQMMTBJMM_input, 0);
    mpz_addmul(DickeQMMTBJMM_input, n, nSk);
    mpz_add(DickeQMMTBJMM_input, DickeQMMTBJMM_input, nSk);
    mpz_addmul_ui(DickeQMMTBJMM_input, coefficient, 4);
    mpz_add(DickeQMMTBJMM_input, DickeQMMTBJMM_input, n);

    // DickeQMMTBJMM_Wcost の計算
    mpz_set(DickeQMMTBJMM_Wcost, DickeQMMTBJMM_input);
    mpz_add(DickeQMMTBJMM_Wcost, DickeQMMTBJMM_Wcost, Dicke_ancila);
    mpz_add(DickeQMMTBJMM_Wcost, DickeQMMTBJMM_Wcost, DickeQMMTBJMMG_ancila);
    mpz_add(DickeQMMTBJMM_Wcost, DickeQMMTBJMM_Wcost, GE_ancila);
    mpz_add(DickeQMMTBJMM_Wcost, DickeQMMTBJMM_Wcost, matmul_ancila1);
    mpz_add(DickeQMMTBJMM_Wcost, DickeQMMTBJMM_Wcost, DickeQMMTBJMMQW_ancila);
    mpz_add(DickeQMMTBJMM_Wcost, DickeQMMTBJMM_Wcost, Ham_ancila);
    mpz_add(DickeQMMTBJMM_Wcost, DickeQMMTBJMM_Wcost, matmul_ancila2);
}