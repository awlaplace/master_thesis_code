#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "../../basic_operation.h"
#include "../Hamming_CBC_params.h"
#include "Prange_params.h"

void compute_grover_iteration(mpz_t grover_iteration, mpz_t n, mpz_t k, mpz_t w);

void compute_CPGloop_cost(mpz_t CPGloop_Gcost, mpz_t CPGloop_Dcost, mpz_t CPGloop_ancila, mpz_t n, mpz_t k);

void compute_QPGloop_PRAM_cost(mpz_t QPGloop_Gcost, mpz_t QPGloop_Dcost, mpz_t QPGloop_ancila, mpz_t CPGloop_Gcost, mpz_t CPGloop_Dcost, mpz_t CPGloop_ancila, mpz_t n, mpz_t Vsize, mpz_t Msize);

void compute_QPGloop_cost(mpz_t QPGloop_Gcost, mpz_t QPGloop_Dcost, mpz_t QPGloop_ancila, mpz_t CPGloop_Gcost, mpz_t CPGloop_Dcost, mpz_t CPGloop_ancila, mpz_t n, mpz_t Vsize, mpz_t Msize);

void compute_CPrange_Gcost(mpz_t CPrange_Gcost, mpz_t CPGloop_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k);

void compute_CPrange_Dcost(mpz_t CPrange_Dcost, mpz_t CPGloop_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k);

void compute_CPrange_Wcost(mpz_t CPrange_Wcost, mpz_t CPGloop_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k);

void compute_QPrange_PRAM_Gcost(mpz_t QPrange_Gcost, mpz_t QPGloop_Gcost, mpz_t CPGloop_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k);

void compute_QPrange_PRAM_Dcost(mpz_t QPrange_Dcost, mpz_t QPGloop_Dcost, mpz_t CPGloop_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k);

void compute_QPrange_PRAM_Wcost(mpz_t QPrange_Wcost, mpz_t QPGloop_ancila, mpz_t CPGloop_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k);

void compute_QPrange_Gcost(mpz_t QPrange_Gcost, mpz_t QPGloop_Gcost, mpz_t CPGloop_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k);

void compute_QPrange_Dcost(mpz_t QPrange_Dcost, mpz_t QPGloop_Dcost, mpz_t CPGloop_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k);

void compute_QPrange_Wcost(mpz_t QPrange_Wcost, mpz_t QPGloop_ancila, mpz_t CPGloop_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k);

void compute_IQPrange_Gcost(mpz_t QPrange_Gcost, mpz_t QPGloop_Gcost, mpz_t CPGloop_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k, mpz_t w);

void compute_IQPrange_Dcost(mpz_t QPrange_Dcost, mpz_t QPGloop_Dcost, mpz_t CPGloop_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k, mpz_t w);

void compute_IQPrange_Wcost(mpz_t QPrange_Wcost, mpz_t QPGloop_ancila, mpz_t CPGloop_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t k, mpz_t w);