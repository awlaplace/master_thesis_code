#ifndef _GRS_PARAMS_H_
#define _GRS_PARAMS_H_

#include <stdio.h>
#include <gmp.h>

extern mpz_t CGRS_n_large_Gcost_r, CGRS_n_large_Dcost_r, CGRS_n_large_Wcost_r;
extern mpz_t CGRS_m_large_Gcost_r, CGRS_m_large_Dcost_r, CGRS_m_large_Wcost_r;
extern mpz_t CGRS_n_large_Gcost, CGRS_n_large_Dcost, CGRS_n_large_Wcost;
extern mpz_t CGRS_m_large_Gcost, CGRS_m_large_Dcost, CGRS_m_large_Wcost;
extern mpz_t CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Dcost, CGRSGloop_n_large_ancila;
extern mpz_t CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Dcost, CGRSGloop_m_large_ancila;

extern mpz_t QGRS_n_large_Gcost_r, QGRS_n_large_Dcost_r, QGRS_n_large_Wcost_r;
extern mpz_t QGRS_m_large_Gcost_r, QGRS_m_large_Dcost_r, QGRS_m_large_Wcost_r;
extern mpz_t QGRS_n_large_Gcost, QGRS_n_large_Dcost, QGRS_n_large_Wcost;
extern mpz_t QGRS_m_large_Gcost, QGRS_m_large_Dcost, QGRS_m_large_Wcost;
extern mpz_t QGRSGloop_n_large_Gcost, QGRSGloop_n_large_Dcost, QGRSGloop_n_large_ancila;
extern mpz_t QGRSGloop_m_large_Gcost, QGRSGloop_m_large_Dcost, QGRSGloop_m_large_ancila;

extern mpz_t IQGRS_n_large_Gcost_r, IQGRS_n_large_Dcost_r, IQGRS_n_large_Wcost_r;
extern mpz_t IQGRS_m_large_Gcost_r, IQGRS_m_large_Dcost_r, IQGRS_m_large_Wcost_r;
extern mpz_t IQGRS_n_large_Gcost, IQGRS_n_large_Dcost, IQGRS_n_large_Wcost;
extern mpz_t IQGRS_m_large_Gcost, IQGRS_m_large_Dcost, IQGRS_m_large_Wcost;
extern mpz_t IQGRSGloop_n_large_Gcost, IQGRSGloop_n_large_Dcost, IQGRSGloop_n_large_ancila;
extern mpz_t IQGRSGloop_m_large_Gcost, IQGRSGloop_m_large_Dcost, IQGRSGloop_m_large_ancila;

void init_GRS_computational_costs();

#endif