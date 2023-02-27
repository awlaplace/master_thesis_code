#ifndef _PRANGE_PARAMS_H_
#define _PRANGE_PARAMS_H_

#include <stdio.h>
#include <gmp.h>

extern mpz_t CPrange_Gcost, CPrange_Dcost, CPrange_Wcost;
extern mpz_t QPrange_Gcost, QPrange_Dcost, QPrange_Wcost;
extern mpz_t IQPrange_Gcost, IQPrange_Dcost, IQPrange_Wcost;
extern mpz_t CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila;
extern mpz_t QPGloop_Gcost, QPGloop_Dcost, QPGloop_ancila;

void init_Prange_computational_costs();

#endif