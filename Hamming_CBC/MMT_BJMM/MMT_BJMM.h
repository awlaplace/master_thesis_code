#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "../../basic_operation.h"
#include "../Hamming_CBC_params.h"
#include "MMT_BJMM_params.h"

void compute_grover_iteration(mpz_t grover_iteration, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p);

void compute_QW_MMT_iteration(mpz_t QW_MMT_iteration_cost, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_CPGloop_cost(mpz_t CPGloop_Gcost, mpz_t CPGloop_Dcost, mpz_t CPGloop_ancila, mpz_t n, mpz_t k);

void compute_each_BD_cost(mpz_t each_BD_Gcost, mpz_t each_BD_Dcost, mpz_t each_BD_ancila, mpz_t k, mpz_t ell, mpz_t p, mpz_t ellprime, mpz_t Xsize, mpz_t Ysize);

void compute_BD_cost(mpz_t BD_Gcost, mpz_t BD_Dcost, mpz_t BD_ancila, mpz_t n, mpz_t k, mpz_t ell, mpz_t p, mpz_t ell1, mpz_t ell2);

void compute_CMMTBJMMGloop_cost(mpz_t CMMTBJMMGloop_Gcost, mpz_t CMMTBJMMGloop_Dcost, mpz_t CMMTBJMMGloop_ancila, mpz_t n, mpz_t k, mpz_t ell, mpz_t p, mpz_t ell1, mpz_t ell2);

void compute_CMMTBJMM_Gcost(mpz_t CMMTBJMM_Gcost, mpz_t CMMTBJMMGloop_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n);

void compute_CMMTBJMM_Dcost(mpz_t CMMTBJMM_Dcost, mpz_t CMMTBJMMGloop_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n);

void compute_CMMTBJMM_Wcost(mpz_t CMMTBJMM_Wcost, mpz_t CMMTBJMMGloop_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k, mpz_t ell, mpz_t p);

void compute_QMMTBJMMQWloop_PRAM_cost(mpz_t QMMTBJMMQWloop_Gcost, mpz_t QMMTBJMMQWloop_Dcost, mpz_t QMMTBJMMQWloop_ancila, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_QMMTBJMMGloop_PRAM_cost(mpz_t QMMTBJMMGloop_Gcost, mpz_t QMMTBJMMGloop_Dcost, mpz_t QMMTBJMMGloop_ancila, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_QMMTBJMM_PRAM_Gcost(mpz_t QMMTBJMM_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_QMMTBJMM_PRAM_Dcost(mpz_t QMMTBJMM_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_QMMTBJMM_PRAM_Wcost(mpz_t QMMTBJMM_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_QMMTBJMMQWloop_cost(mpz_t QMMTBJMMQWloop_Gcost, mpz_t QMMTBJMMQWloop_Dcost, mpz_t QMMTBJMMQWloop_ancila, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_QMMTBJMMQWloopout_cost(mpz_t QMMTBJMMQWloopout_Gcost, mpz_t QMMTBJMMQWloopout_Dcost, mpz_t QMMTBJMMQWloopout_ancila, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_QMMTBJMMGloop_cost(mpz_t QMMTBJMMGloop_Gcost, mpz_t QMMTBJMMGloop_Dcost, mpz_t QMMTBJMMGloop_ancila, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_QMMTBJMM_Gcost(mpz_t QMMTBJMM_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_QMMTBJMM_Dcost(mpz_t QMMTBJMM_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_QMMTBJMM_Wcost(mpz_t QMMTBJMM_Wcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_G4SP_cost(mpz_t G4SP_Gcost, mpz_t G4SP_Dcost, mpz_t G4SP_ancila, mpz_t n, mpz_t k, mpz_t ell, mpz_t p);

void compute_DickeQMMTBJMMQWloop_cost(mpz_t DickeQMMTBJMMQWloop_Gcost, mpz_t DickeQMMTBJMMQWloop_Dcost, mpz_t DickeQMMTBJMMQWloop_ancila, mpz_t n, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_DickeQMMTBJMMGloop_cost(mpz_t DickeQMMTBJMMGloop_Gcost, mpz_t DickeQMMTBJMMGloop_Dcost, mpz_t DickeQMMTBJMMGloop_ancila, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_DickeQMMTBJMM_Gcost(mpz_t DickeQMMTBJMM_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_DickeQMMTBJMM_Dcost(mpz_t DickeQMMTBJMM_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_DickeQMMTBJMM_Wcost(mpz_t DickeQMMTBJMM_Wcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_IQMMTBJMMQWloop_cost(mpz_t IQMMTBJMMQWloop_Gcost, mpz_t IQMMTBJMMQWloop_Dcost, mpz_t IQMMTBJMMQWloop_ancila, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_IQMMTBJMMQWloopout_cost(mpz_t IQMMTBJMMQWloopout_Gcost, mpz_t IQMMTBJMMQWloopout_Dcost, mpz_t IQMMTBJMMQWloopout_ancila, mpz_t k, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_IQMMTBJMMGloop_cost(mpz_t QMMTBJMMGloop_Gcost, mpz_t QMMTBJMMGloop_Dcost, mpz_t QMMTBJMMGloop_ancila, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_IQMMTBJMM_Gcost(mpz_t QMMTBJMM_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_IQMMTBJMM_Dcost(mpz_t QMMTBJMM_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);

void compute_IQMMTBJMM_Wcost(mpz_t QMMTBJMM_Wcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t QW_MMT_iteration, mpz_t QW_MMT_processor_num, mpz_t n, mpz_t k, mpz_t w, mpz_t ell, mpz_t p, mpz_t epsilon);