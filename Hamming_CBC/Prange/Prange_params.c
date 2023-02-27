#include <stdio.h>
#include <gmp.h>

mpz_t CPrange_Gcost, CPrange_Dcost, CPrange_Wcost;
mpz_t QPrange_Gcost, QPrange_Dcost, QPrange_Wcost;
mpz_t IQPrange_Gcost, IQPrange_Dcost, IQPrange_Wcost;
mpz_t CPGloop_Gcost, CPGloop_Dcost, CPGloop_ancila;
mpz_t QPGloop_Gcost, QPGloop_Dcost, QPGloop_ancila;

void init_Prange_computational_costs(){
    mpz_init_set_ui(CPrange_Gcost, 0);
    mpz_init_set_ui(CPrange_Dcost, 0);
    mpz_init_set_ui(CPrange_Wcost, 0);
    mpz_init_set_ui(QPrange_Gcost, 0);
    mpz_init_set_ui(QPrange_Dcost, 0);
    mpz_init_set_ui(QPrange_Wcost, 0);
    mpz_init_set_ui(IQPrange_Gcost, 0);
    mpz_init_set_ui(IQPrange_Dcost, 0);
    mpz_init_set_ui(IQPrange_Wcost, 0);
    mpz_init_set_ui(CPGloop_Gcost, 0);
    mpz_init_set_ui(CPGloop_Dcost, 0);
    mpz_init_set_ui(CPGloop_ancila, 0);
    mpz_init_set_ui(QPGloop_Gcost, 0);
    mpz_init_set_ui(QPGloop_Dcost, 0);
    mpz_init_set_ui(QPGloop_ancila, 0);
}