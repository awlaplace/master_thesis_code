#include <stdio.h>
#include <gmp.h>

mpz_t CGRS_n_large_Gcost_r, CGRS_n_large_Dcost_r, CGRS_n_large_Wcost_r;
mpz_t CGRS_m_large_Gcost_r, CGRS_m_large_Dcost_r, CGRS_m_large_Wcost_r;
mpz_t CGRS_n_large_Gcost, CGRS_n_large_Dcost, CGRS_n_large_Wcost;
mpz_t CGRS_m_large_Gcost, CGRS_m_large_Dcost, CGRS_m_large_Wcost;
mpz_t CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Dcost, CGRSGloop_n_large_ancila;
mpz_t CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Dcost, CGRSGloop_m_large_ancila;

mpz_t QGRS_n_large_Gcost_r, QGRS_n_large_Dcost_r, QGRS_n_large_Wcost_r;
mpz_t QGRS_m_large_Gcost_r, QGRS_m_large_Dcost_r, QGRS_m_large_Wcost_r;
mpz_t QGRS_n_large_Gcost, QGRS_n_large_Dcost, QGRS_n_large_Wcost;
mpz_t QGRS_m_large_Gcost, QGRS_m_large_Dcost, QGRS_m_large_Wcost;
mpz_t QGRSGloop_n_large_Gcost, QGRSGloop_n_large_Dcost, QGRSGloop_n_large_ancila;
mpz_t QGRSGloop_m_large_Gcost, QGRSGloop_m_large_Dcost, QGRSGloop_m_large_ancila;

mpz_t IQGRS_n_large_Gcost_r, IQGRS_n_large_Dcost_r, IQGRS_n_large_Wcost_r;
mpz_t IQGRS_m_large_Gcost_r, IQGRS_m_large_Dcost_r, IQGRS_m_large_Wcost_r;
mpz_t IQGRS_n_large_Gcost, IQGRS_n_large_Dcost, IQGRS_n_large_Wcost;
mpz_t IQGRS_m_large_Gcost, IQGRS_m_large_Dcost, IQGRS_m_large_Wcost;
mpz_t IQGRSGloop_n_large_Gcost, IQGRSGloop_n_large_Dcost, IQGRSGloop_n_large_ancila;
mpz_t IQGRSGloop_m_large_Gcost, IQGRSGloop_m_large_Dcost, IQGRSGloop_m_large_ancila;

void init_GRS_computational_costs(){
    mpz_init_set_ui(CGRS_n_large_Gcost_r, 0);
    mpz_init_set_ui(CGRS_n_large_Dcost_r, 0);
    mpz_init_set_ui(CGRS_n_large_Wcost_r, 0);
    mpz_init_set_ui(CGRS_m_large_Gcost_r, 0);
    mpz_init_set_ui(CGRS_m_large_Dcost_r, 0);
    mpz_init_set_ui(CGRS_m_large_Wcost_r, 0);
    mpz_init_set_ui(CGRS_n_large_Gcost, 0);
    mpz_init_set_ui(CGRS_n_large_Dcost, 0);
    mpz_init_set_ui(CGRS_n_large_Wcost, 0);
    mpz_init_set_ui(CGRS_m_large_Gcost, 0);
    mpz_init_set_ui(CGRS_m_large_Dcost, 0);
    mpz_init_set_ui(CGRS_m_large_Wcost, 0);
    mpz_init_set_ui(CGRSGloop_n_large_Gcost, 0);
    mpz_init_set_ui(CGRSGloop_n_large_Dcost, 0);
    mpz_init_set_ui(CGRSGloop_n_large_ancila, 0);
    mpz_init_set_ui(CGRSGloop_m_large_Gcost, 0);
    mpz_init_set_ui(CGRSGloop_m_large_Dcost, 0);
    mpz_init_set_ui(CGRSGloop_m_large_ancila, 0);

    mpz_init_set_ui(QGRS_n_large_Gcost_r, 0);
    mpz_init_set_ui(QGRS_n_large_Dcost_r, 0);
    mpz_init_set_ui(QGRS_n_large_Wcost_r, 0);
    mpz_init_set_ui(QGRS_m_large_Gcost_r, 0);
    mpz_init_set_ui(QGRS_m_large_Dcost_r, 0);
    mpz_init_set_ui(QGRS_m_large_Wcost_r, 0);
    mpz_init_set_ui(QGRS_n_large_Gcost, 0);
    mpz_init_set_ui(QGRS_n_large_Dcost, 0);
    mpz_init_set_ui(QGRS_n_large_Wcost, 0);
    mpz_init_set_ui(QGRS_m_large_Gcost, 0);
    mpz_init_set_ui(QGRS_m_large_Dcost, 0);
    mpz_init_set_ui(QGRS_m_large_Wcost, 0);
    mpz_init_set_ui(QGRSGloop_n_large_Gcost, 0);
    mpz_init_set_ui(QGRSGloop_n_large_Dcost, 0);
    mpz_init_set_ui(QGRSGloop_n_large_ancila, 0);
    mpz_init_set_ui(QGRSGloop_m_large_Gcost, 0);
    mpz_init_set_ui(QGRSGloop_m_large_Dcost, 0);
    mpz_init_set_ui(QGRSGloop_m_large_ancila, 0);

    mpz_init_set_ui(IQGRS_n_large_Gcost_r, 0);
    mpz_init_set_ui(IQGRS_n_large_Dcost_r, 0);
    mpz_init_set_ui(IQGRS_n_large_Wcost_r, 0);
    mpz_init_set_ui(IQGRS_m_large_Gcost_r, 0);
    mpz_init_set_ui(IQGRS_m_large_Dcost_r, 0);
    mpz_init_set_ui(IQGRS_m_large_Wcost_r, 0);
    mpz_init_set_ui(IQGRS_n_large_Gcost, 0);
    mpz_init_set_ui(IQGRS_n_large_Dcost, 0);
    mpz_init_set_ui(IQGRS_n_large_Wcost, 0);
    mpz_init_set_ui(IQGRS_m_large_Gcost, 0);
    mpz_init_set_ui(IQGRS_m_large_Dcost, 0);
    mpz_init_set_ui(IQGRS_m_large_Wcost, 0);
    mpz_init_set_ui(IQGRSGloop_n_large_Gcost, 0);
    mpz_init_set_ui(IQGRSGloop_n_large_Dcost, 0);
    mpz_init_set_ui(IQGRSGloop_n_large_ancila, 0);
    mpz_init_set_ui(IQGRSGloop_m_large_Gcost, 0);
    mpz_init_set_ui(IQGRSGloop_m_large_Dcost, 0);
    mpz_init_set_ui(IQGRSGloop_m_large_ancila, 0);
}