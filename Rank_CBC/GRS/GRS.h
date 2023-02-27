#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "../../basic_operation.h"
#include "../Rank_CBC_params.h"
#include "GRS_params.h"

void compute_grover_iteration_n_large(mpz_t grover_iteration, mpz_t m, mpz_t r, mpz_t w);

void compute_grover_iteration_m_large(mpz_t grover_iteration, mpz_t n, mpz_t r, mpz_t w);

void compute_CGRSGloop_n_large(mpz_t CGRSGloop_n_large_Gcost, mpz_t CGRSGloop_n_large_Dcost, mpz_t CGRSGloop_n_large_ancila, mpz_t m, mpz_t n, mpz_t k, mpz_t r);

void compute_CGRS_n_large_Gcost(mpz_t CGRS_n_large_Gcost, mpz_t CGRSGloop_n_large_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num);

void compute_CGRS_n_large_Dcost(mpz_t CGRS_n_large_Dcost, mpz_t CGRSGloop_n_large_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num);

void compute_CGRS_n_large_Wcost(mpz_t CGRS_n_large_Wcost, mpz_t CGRSGloop_n_large_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t m, mpz_t k);

void compute_CGRS_m_large_Gcost(mpz_t CGRS_m_large_Gcost, mpz_t CGRSGloop_m_large_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num);

void compute_CGRS_m_large_Dcost(mpz_t CGRS_m_large_Dcost, mpz_t CGRSGloop_m_large_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num);

void compute_CGRS_m_large_Wcost(mpz_t CGRS_m_large_Wcost, mpz_t CGRSGloop_m_large_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t m, mpz_t k);

void compute_QGRSGloop_n_large(mpz_t QGRSGloop_n_large_Gcost, mpz_t QGRSGloop_n_large_Dcost, mpz_t QGRSGloop_n_large_ancila, mpz_t m, mpz_t n, mpz_t k,  mpz_t w, mpz_t r);

void compute_QGRS_n_large_Gcost(mpz_t QGRS_n_large_Gcost, mpz_t QGRSGloop_n_large_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t m, mpz_t r);

void compute_QGRS_n_large_Dcost(mpz_t QGRS_n_large_Dcost, mpz_t QGRSGloop_n_large_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t m, mpz_t r);

void compute_QGRS_n_large_Wcost(mpz_t QGRS_n_large_Wcost, mpz_t QGRSGloop_n_large_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t m, mpz_t k, mpz_t r);

void compute_CGRSGloop_m_large(mpz_t CGRSGloop_m_large_Gcost, mpz_t CGRSGloop_m_large_Dcost, mpz_t CGRSGloop_m_large_ancila, mpz_t m, mpz_t n, mpz_t k, mpz_t r);

void compute_QGRSGloop_m_large(mpz_t CGRSGloop_m_large_Gcost, mpz_t CGRSGloop_m_large_Dcost, mpz_t CGRSGloop_m_large_ancila, mpz_t m, mpz_t n, mpz_t k, mpz_t w, mpz_t r);

void compute_QGRS_m_large_Gcost(mpz_t QGRS_m_large_Gcost, mpz_t QGRSGloop_m_large_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t r);

void compute_QGRS_m_large_Dcost(mpz_t QGRS_m_large_Dcost, mpz_t QGRSGloop_m_large_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t r);

void compute_QGRS_m_large_Wcost(mpz_t QGRS_m_large_Wcost, mpz_t QGRSGloop_m_large_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t m, mpz_t k, mpz_t r);

void compute_IQGRSGloop_n_large(mpz_t IQGRSGloop_n_large_Gcost, mpz_t IQGRSGloop_n_large_Dcost, mpz_t IQGRSGloop_n_large_ancila, mpz_t m, mpz_t n, mpz_t k, mpz_t w, mpz_t r);

void compute_IQGRS_n_large_Gcost(mpz_t IQGRS_n_large_Gcost, mpz_t IQGRSGloop_n_large_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t m, mpz_t w, mpz_t r);

void compute_IQGRS_n_large_Dcost(mpz_t IQGRS_n_large_Dcost, mpz_t IQGRSGloop_n_large_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t m, mpz_t w, mpz_t r);

void compute_IQGRS_n_large_Wcost(mpz_t IQGRS_n_large_Wcost, mpz_t IQGRSGloop_n_large_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t m, mpz_t k, mpz_t w, mpz_t r);

void compute_IQGRSGloop_m_large(mpz_t IQGRSGloop_m_large_Gcost, mpz_t IQGRSGloop_m_large_Dcost, mpz_t IQGRSGloop_m_large_ancila, mpz_t m, mpz_t n, mpz_t k, mpz_t w, mpz_t r);

void compute_IQGRS_m_large_Gcost(mpz_t IQGRS_m_large_Gcost, mpz_t IQGRSGloop_m_large_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t w, mpz_t r);

void compute_IQGRS_m_large_Dcost(mpz_t IQGRS_m_large_Dcost, mpz_t IQGRSGloop_m_large_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t w, mpz_t r);

void compute_IQGRS_m_large_Wcost(mpz_t IQGRS_m_large_Wcost, mpz_t IQGRSGloop_m_large_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t m, mpz_t k, mpz_t w, mpz_t r);
