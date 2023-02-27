#include <stdio.h>
#include <gmp.h>

mpz_t CMMTBJMM_Gcost, CMMTBJMM_Dcost, CMMTBJMM_Wcost;
mpz_t QMMTBJMM_Gcost, QMMTBJMM_Dcost, QMMTBJMM_Wcost;
mpz_t IQMMTBJMM_Gcost, IQMMTBJMM_Dcost, IQMMTBJMM_Wcost;
mpz_t DickeQMMTBJMM_Gcost, DickeQMMTBJMM_Dcost, DickeQMMTBJMM_Wcost;
mpz_t CMMTBJMMGloop_Gcost, CMMTBJMMGloop_Dcost, CMMTBJMMGloop_ancila;
mpz_t QMMTBJMMQWloop_Gcost, QMMTBJMMQWloop_Dcost, QMMTBJMMQWloop_ancila;
mpz_t QMMTBJMMQWloopout_Gcost, QMMTBJMMQWloopout_Dcost, QMMTBJMMQWloopout_ancila;
mpz_t QMMTBJMMGloop_Gcost, QMMTBJMMGloop_Dcost, QMMTBJMMGloop_ancila;
mpz_t DickeQMMTBJMMQWloop_Gcost, DickeQMMTBJMMQWloop_Dcost, DickeQMMTBJMMQWloop_ancila;
mpz_t DickeQMMTBJMMGloop_Gcost, DickeQMMTBJMMGloop_Dcost, DickeQMMTBJMMGloop_ancila;
mpz_t IQMMTBJMMQWloop_Gcost, IQMMTBJMMQWloop_Dcost, IQMMTBJMMQWloop_ancila;
mpz_t IQMMTBJMMQWloopout_Gcost, IQMMTBJMMQWloopout_Dcost, IQMMTBJMMQWloopout_ancila;
mpz_t IQMMTBJMMGloop_Gcost, IQMMTBJMMGloop_Dcost, IQMMTBJMMGloop_ancila;


void init_MMT_BJMM_computational_costs(){
    mpz_init_set_ui(CMMTBJMM_Gcost, 0);
    mpz_init_set_ui(CMMTBJMM_Dcost, 0);
    mpz_init_set_ui(CMMTBJMM_Wcost, 0);
    mpz_init_set_ui(QMMTBJMM_Gcost, 0);
    mpz_init_set_ui(QMMTBJMM_Dcost, 0);
    mpz_init_set_ui(QMMTBJMM_Wcost, 0);
    mpz_init_set_ui(IQMMTBJMM_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMM_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMM_Wcost, 0);
    mpz_init_set_ui(DickeQMMTBJMM_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMM_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMM_Wcost, 0);
    mpz_init_set_ui(CMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(CMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(CMMTBJMMGloop_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMQWloopout_ancila, 0);
    mpz_init_set_ui(QMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(QMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(QMMTBJMMGloop_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(DickeQMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(DickeQMMTBJMMGloop_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloop_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMQWloopout_ancila, 0);
    mpz_init_set_ui(IQMMTBJMMGloop_Gcost, 0);
    mpz_init_set_ui(IQMMTBJMMGloop_Dcost, 0);
    mpz_init_set_ui(IQMMTBJMMGloop_ancila, 0);
}