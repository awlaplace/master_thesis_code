/*
  gcc -c classical_MMT_BJMM.c -lgmp && gcc -o classical_MMT_BJMM ../../basic_operation.o ../Hamming_CBC_params.o MMT_BJMM_params.o classical_MMT_BJMM.o -lgmp && ./classical_MMT_BJMM
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

        mpz_init_set_ui(grover_processor_num, 1);

        init_MMT_BJMM_computational_costs();

        int min_CMMTBJMM_Gcost_log = 0;
        int min_CMMTBJMM_Dcost_log = MAX_DEPTH_log + 1;
        int min_CMMTBJMM_Wcost_log = 0;

        for(int int_p = 1; int_p < SDP_instance_list[index][2] / 4; int_p++){
            mpz_set_ui(p, 4 * int_p);
            int int_wSp = SDP_instance_list[index][2] - int_p;
            mpz_set_ui(wSp, int_wSp);
            for(int int_ell = 1; int_ell < SDP_instance_list[index][0] - SDP_instance_list[index][1]; int_ell++){
                mpz_set_ui(ell, int_ell);
                for(int int_ell1 = 0; int_ell1 < int_ell; int_ell1++){
                    mpz_set_ui(ell1, int_ell1);
                    int int_ell2 = int_ell - int_ell1;
                    mpz_set_ui(ell2, int_ell2);
                    
                    // 二項係数が意味を持つ範囲に限定
                    if((mpz_cmp(n, w) > 0)){
                        compute_grover_iteration(grover_iteration, n, k, w, ell, p);

                        compute_CMMTBJMMGloop_cost(CMMTBJMMGloop_Gcost, CMMTBJMMGloop_Dcost, CMMTBJMMGloop_ancila, n, k, ell, p, ell1, ell2);
                        
                        // CMMTBJMM_Dcost を計算する
                        compute_CMMTBJMM_Dcost(CMMTBJMM_Dcost, CMMTBJMMGloop_Dcost, grover_iteration, grover_processor_num, n);

                        // grover_processor_num の個数を調整
                        while(mpz_cmp(CMMTBJMM_Dcost, MAX_DEPTH) > 0){
                            mpz_mul_ui(grover_processor_num, grover_processor_num, 2);
                            compute_CMMTBJMM_Dcost(CMMTBJMM_Dcost, CMMTBJMMGloop_Dcost, grover_iteration, grover_processor_num, n);
                        }

                        // CPrange_Gcostを計算する
                        compute_CMMTBJMM_Gcost(CMMTBJMM_Gcost, CMMTBJMMGloop_Gcost, grover_iteration, grover_processor_num, n);
                        
                        // CPrange_Wcostを計算する
                        compute_CMMTBJMM_Wcost(CMMTBJMM_Wcost, CMMTBJMMGloop_ancila, grover_iteration, grover_processor_num, n, k, ell, p);

                        // compute cost_log
                        int CMMTBJMM_Gcost_loop_log = compute_log(CMMTBJMM_Gcost);
                        int CMMTBJMM_Dcost_loop_log = compute_log(CMMTBJMM_Dcost);
                        int CMMTBJMM_Wcost_loop_log = compute_log(CMMTBJMM_Wcost);

                        if((min_CMMTBJMM_Dcost_log == MAX_DEPTH_log + 1) || (CMMTBJMM_Gcost_loop_log < min_CMMTBJMM_Gcost_log)){
                            min_CMMTBJMM_Gcost_log = CMMTBJMM_Gcost_loop_log;
                            min_CMMTBJMM_Dcost_log = CMMTBJMM_Dcost_loop_log;
                            MMT_BJMM_result_list[index][3] = int_p * 4;
                            MMT_BJMM_result_list[index][4] = int_ell;
                            MMT_BJMM_result_list[index][5] = int_ell1;
                            MMT_BJMM_result_list[index][6] = int_ell2;
                            MMT_BJMM_result_list[index][7] = 0;
                            MMT_BJMM_result_list[index][8] = CMMTBJMM_Gcost_loop_log;
                            MMT_BJMM_result_list[index][9] = CMMTBJMM_Dcost_loop_log;
                            MMT_BJMM_result_list[index][10] = CMMTBJMM_Wcost_loop_log;
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
 each_BD_Gcost, each_BD_Dcost, each_BD_ancila, k, ell, p, ellprime, Xsize, Ysize を引数にとったとき，
 Birthday Decoding の G-cost, D-cost, アンシラビット数を
 それぞれ each_BD_Gcost, each_BD_Dcost, each_BD_ancila に格納する関数
*/
void compute_each_BD_cost(mpz_t each_BD_Gcost, mpz_t each_BD_Dcost, mpz_t each_BD_ancila, mpz_t k, mpz_t ell, mpz_t p, mpz_t ellprime, mpz_t Xsize, mpz_t Ysize){
    mpz_t add_Gcost, add_Dcost, add_ancila;
    mpz_t matmul_Gcost, matmul_Dcost, matmul_ancila;
    mpz_init_set_ui(add_Gcost, 0);
    mpz_init_set_ui(add_Dcost, 0);
    mpz_init_set_ui(add_ancila, 0);
    mpz_init_set_ui(matmul_Gcost, 0);
    mpz_init_set_ui(matmul_Dcost, 0);
    mpz_init_set_ui(matmul_ancila, 0);

    mpz_t kAell, one;
    mpz_init(kAell);
    mpz_add(kAell, k, ell);
    mpz_init_set_ui(one, 1);

    compute_add_cost(add_Gcost, add_Dcost, add_ancila, kAell);
    compute_matmul_cost(matmul_Gcost, matmul_Dcost, matmul_ancila, ellprime, kAell, one);

    mpz_t XYsize;
    mpz_init(XYsize);
    mpz_mul(XYsize, Xsize, Ysize);

    mpz_addmul(each_BD_Gcost, matmul_Gcost, Xsize);
    mpz_addmul(each_BD_Gcost, matmul_Gcost, XYsize);
    mpz_addmul(each_BD_Gcost, add_Gcost, Ysize);
    mpz_addmul(each_BD_Gcost, add_Gcost, XYsize);

    mpz_t XsizeMmutmul_Dcost;
    mpz_init(XsizeMmutmul_Dcost);
    mpz_mul(XsizeMmutmul_Dcost, matmul_Dcost, Xsize);

    mpz_set(each_BD_Dcost, add_Dcost);
    mpz_mul(each_BD_Dcost, each_BD_Dcost, Xsize);
    if(mpz_cmp(matmul_Dcost, each_BD_Dcost) > 0){
        mpz_set(each_BD_Dcost, matmul_Dcost);
    }
    mpz_mul(each_BD_Dcost, each_BD_Dcost, Ysize);
    if(mpz_cmp(XsizeMmutmul_Dcost, each_BD_Dcost) > 0){
        mpz_set(each_BD_Dcost, XsizeMmutmul_Dcost);
    }

    mpz_addmul(each_BD_ancila, matmul_ancila, Xsize);
    mpz_addmul(each_BD_ancila, matmul_ancila, XYsize);
    mpz_addmul(each_BD_ancila, add_ancila, Ysize);
    mpz_addmul(each_BD_ancila, add_ancila, XYsize);
}


/*
 BD_Gcost, BD_Dcost, BD_ancila, n, k, ell, p, ell1, ell2 を引数にとったとき，
 G4SP_BD の G-cost, D-cost, アンシラビット数を
 それぞれ BD_Gcost, BD_Dcost, BD_ancila に格納する関数
*/
void compute_BD_cost(mpz_t BD_Gcost, mpz_t BD_Dcost, mpz_t BD_ancila, mpz_t n, mpz_t k, mpz_t ell, mpz_t p, mpz_t ell1, mpz_t ell2){
    mpz_t BD_6_Gcost, BD_6_Dcost, BD_6_ancila;
    mpz_t BD_7_Gcost, BD_7_Dcost, BD_7_ancila;
    mpz_t BD_8_Gcost, BD_8_Dcost, BD_8_ancila;
    mpz_t add_Gcost, add_Dcost, add_ancila;
    mpz_t matmul_Gcost, matmul_Dcost, matmul_ancila;
    mpz_init_set_ui(BD_6_Gcost, 0);
    mpz_init_set_ui(BD_6_Dcost, 0);
    mpz_init_set_ui(BD_6_ancila, 0);
    mpz_init_set_ui(BD_7_Gcost, 0);
    mpz_init_set_ui(BD_7_Dcost, 0);
    mpz_init_set_ui(BD_7_ancila, 0);
    mpz_init_set_ui(BD_8_Gcost, 0);
    mpz_init_set_ui(BD_8_Dcost, 0);
    mpz_init_set_ui(BD_8_ancila, 0);
    mpz_init_set_ui(add_Gcost, 0);
    mpz_init_set_ui(add_Dcost, 0);
    mpz_init_set_ui(add_ancila, 0);
    mpz_init_set_ui(matmul_Gcost, 0);
    mpz_init_set_ui(matmul_Dcost, 0);
    mpz_init_set_ui(matmul_ancila, 0);

    mpz_t kAell, kAellD4, pD4, nSkSell, one;
    mpz_init(kAell);
    mpz_init(kAellD4);
    mpz_init(pD4);
    mpz_add(kAell, k, ell);
    mpz_cdiv_q_ui(kAellD4, kAell, 4);
    mpz_cdiv_q_ui(pD4, p, 4);
    mpz_init_set(nSkSell, n);
    mpz_sub(nSkSell, nSkSell, k);
    mpz_sub(nSkSell, nSkSell, ell);
    mpz_init_set_ui(one, 1);

    mpz_t Vsize;
    mpz_init_set_ui(Vsize, 1);
    binomial_coefficient(Vsize, kAellD4, pD4);

    compute_each_BD_cost(BD_6_Gcost, BD_6_Dcost, BD_6_ancila, k, ell, p, ell1, Vsize, Vsize);
    compute_each_BD_cost(BD_7_Gcost, BD_7_Dcost, BD_7_ancila, k, ell, p, ell1, Vsize, Vsize);
    compute_each_BD_cost(BD_8_Gcost, BD_8_Dcost, BD_8_ancila, k, ell, p, ell2, Vsize, Vsize);
    compute_matmul_cost(matmul_Gcost, matmul_Dcost, matmul_ancila, nSkSell, kAell, one);
    compute_add_cost(add_Gcost, add_Dcost, add_ancila, kAell);

    mpz_set(BD_Gcost, add_Gcost);
    mpz_add(BD_Gcost, BD_Gcost, matmul_Gcost);
    mpz_mul(BD_Gcost, BD_Gcost, Vsize);
    mpz_add(BD_Gcost, BD_Gcost, BD_6_Gcost);
    mpz_add(BD_Gcost, BD_Gcost, BD_7_Gcost);
    mpz_add(BD_Gcost, BD_Gcost, BD_8_Gcost);

    mpz_set(BD_Dcost, kAell);
    mpz_mul(BD_Dcost, BD_Dcost, Vsize);

    mpz_set(BD_ancila, add_ancila);
    mpz_add(BD_ancila, BD_ancila, matmul_ancila);
    mpz_mul(BD_ancila, BD_ancila, Vsize);
    mpz_add(BD_ancila, BD_ancila, BD_6_ancila);
    mpz_add(BD_ancila, BD_ancila, BD_7_ancila);
    mpz_add(BD_ancila, BD_ancila, BD_8_ancila);
}


/*
 CMMTBJMMGloop_Gcost, CMMTBJMMGloop_Dcost, CMMTBJMMGloop_ancila, n, k, ell, p, ell1, ell2 を引数にとったとき，
 古典MMT/BJMMアルゴリズムの for 文のループ回数の G-cost, D-cost, アンシラビット数を
 それぞれ CMMTBJMMGloop_Gcost, CMMTBJMMGloop_Dcost, CMMTBJMMGloop_ancila に格納する関数
*/
void compute_CMMTBJMMGloop_cost(mpz_t CMMTBJMMGloop_Gcost, mpz_t CMMTBJMMGloop_Dcost, mpz_t CMMTBJMMGloop_ancila, mpz_t n, mpz_t k, mpz_t ell, mpz_t p, mpz_t ell1, mpz_t ell2){
    mpz_t CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila;
    mpz_t BD_Gcost, BD_Dcost, BD_ancila;
    mpz_init_set_ui(CPGloop_Gcost, 0);
    mpz_init_set_ui(CPGloop_Dcost, 0);
    mpz_init_set_ui(CPGloop_ancila, 0);
    mpz_init_set_ui(BD_Gcost, 0);
    mpz_init_set_ui(BD_Dcost, 0);
    mpz_init_set_ui(BD_ancila, 0);

    mpz_t nMn;
    mpz_init_set(nMn, n);
    mpz_mul(nMn, nMn, n);
    
    compute_CPGloop_cost(CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila, n, k);
    compute_BD_cost(BD_Gcost, BD_Dcost, BD_ancila, n, k, ell, p, ell1, ell2);

    mpz_set(CMMTBJMMGloop_Gcost, CPGloop_Gcost);
    mpz_add(CMMTBJMMGloop_Gcost, CMMTBJMMGloop_Gcost, BD_Gcost);

    mpz_set(CMMTBJMMGloop_Dcost, BD_Dcost);

    mpz_set(CMMTBJMMGloop_ancila, CPGloop_ancila);
    mpz_add(CMMTBJMMGloop_ancila, CMMTBJMMGloop_ancila, BD_ancila);
    mpz_add(CMMTBJMMGloop_ancila, CMMTBJMMGloop_ancila, nMn);
}


/*
 CMMTBJMM_Gcost, CMMTBJMMGloop_Gcost, grover_iteration, grover_processor_num, n を引数にとったとき，
 古典MMT/BJMMアルゴリズムの全体の G-cost を CMMTBJMM_Gcost に格納する関数
*/
void compute_CMMTBJMM_Gcost(mpz_t CMMTBJMM_Gcost, mpz_t CMMTBJMMGloop_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n){
    mpz_t c24nMn;
    mpz_init_set_ui(c24nMn, 24);
    mpz_mul(c24nMn, c24nMn, n);
    mpz_mul(c24nMn, c24nMn, n);

    mpz_set(CMMTBJMM_Gcost, CMMTBJMMGloop_Gcost);
    mpz_mul(CMMTBJMM_Gcost, CMMTBJMM_Gcost, grover_iteration);
    mpz_mul(CMMTBJMM_Gcost, CMMTBJMM_Gcost, grover_processor_num);
    mpz_add(CMMTBJMM_Gcost, CMMTBJMM_Gcost, c24nMn);
}


/*
 CMMTBJMM_Dcost, CMMTBJMMGloop_Dcost, grover_iteration, grover_processor_num, n を引数にとったとき，
 古典MMT/BJMMアルゴリズムの全体の D-cost を CMMTBJMM_Dcost に格納する関数
*/
void compute_CMMTBJMM_Dcost(mpz_t CMMTBJMM_Dcost, mpz_t CMMTBJMMGloop_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n){
    mpz_t c16n;
    mpz_init_set_ui(c16n, 16);
    mpz_mul(c16n, c16n, n);

    mpz_set(CMMTBJMM_Dcost, CMMTBJMMGloop_Dcost);
    mpz_mul(CMMTBJMM_Dcost, CMMTBJMM_Dcost, grover_iteration);
    mpz_cdiv_q(CMMTBJMM_Dcost, CMMTBJMM_Dcost, grover_processor_num);

    if(mpz_cmp(c16n, CMMTBJMM_Dcost) > 0){
        mpz_set(CMMTBJMM_Dcost, c16n);
    }
}


/*
 CMMTBJMM_Wcost, CMMTBJMMGloop_ancila, grover_iteration, grover_processor_num, n を引数にとったとき，
 古典MMT/BJMMアルゴリズムの全体の W-cost を CMMTBJMM_Wcost に格納する関数
*/
void compute_CMMTBJMM_Wcost(mpz_t CMMTBJMM_Wcost, mpz_t CMMTBJMMGloop_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k, mpz_t ell, mpz_t p){
    mpz_t CMMTBJMM_input;
    
    mpz_t nSk, kAell, kAellD4, pD4;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init(kAell);
    mpz_add(kAell, k, ell);
    mpz_init(kAellD4);
    mpz_cdiv_q_ui(kAellD4, kAell, 4);
    mpz_init(pD4);
    mpz_cdiv_q_ui(pD4, p, 4);

    mpz_init_set_ui(CMMTBJMM_input, 1);
    binomial_coefficient(CMMTBJMM_input, kAellD4, pD4);
    mpz_mul(CMMTBJMM_input, CMMTBJMM_input, kAell);
    mpz_mul_ui(CMMTBJMM_input, CMMTBJMM_input, 4);
    mpz_add(CMMTBJMM_input, CMMTBJMM_input, n);
    mpz_addmul(CMMTBJMM_input, nSk, n);
    mpz_add(CMMTBJMM_input, CMMTBJMM_input, nSk);

    mpz_set(CMMTBJMM_Wcost, CMMTBJMM_input);
    mpz_addmul(CMMTBJMM_Wcost, CMMTBJMMGloop_ancila, grover_iteration);
}