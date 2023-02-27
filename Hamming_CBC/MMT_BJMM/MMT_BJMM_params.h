#ifndef _MMT_BJMM_PARAMS_H_
#define _MMT_BJMM_PARAMS_H_

#include <stdio.h>
#include <gmp.h>

extern mpz_t CMMTBJMM_Gcost, CMMTBJMM_Dcost, CMMTBJMM_Wcost;
extern mpz_t QMMTBJMM_Gcost, QMMTBJMM_Dcost, QMMTBJMM_Wcost;
extern mpz_t IQMMTBJMM_Gcost, IQMMTBJMM_Dcost, IQMMTBJMM_Wcost;
extern mpz_t DickeQMMTBJMM_Gcost, DickeQMMTBJMM_Dcost, DickeQMMTBJMM_Wcost;
extern mpz_t CMMTBJMMGloop_Gcost, CMMTBJMMGloop_Dcost, CMMTBJMMGloop_ancila;
extern mpz_t QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila;
extern mpz_t QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWloopout_ancila;
extern mpz_t QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Dcost, QMMTBJMMGloop_ancila;
extern mpz_t DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila;
extern mpz_t DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Dcost, DickeQMMTBJMMGloop_ancila;
extern mpz_t IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila;
extern mpz_t IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila;
extern mpz_t IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Dcost, IQMMTBJMMGloop_ancila;

void init_MMT_BJMM_computational_costs();

#endif