/**
 * @file Lagranian1D.h
 * @brief 1D Lagrangian scheme for relativistic Euler equations
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2018-01-17
 */
#ifndef Lagranian1D_H
#define Lagranian1D_H

#include <iostream>
#include <assert.h>
#include <cstdlib>
#include <omp.h>
#include "para.h"
#include "vvector.h"
#include "mvector.h"

#define MAT mvector<mvector<double,3>,3>

typedef mvector<double,3> bU; //(1/D, m/D, E/D), mass D average

class GHOST {
  public:
    bU Con;
    bU Pri;
    double Gamma;
};

class Lagranian1D {
  private:
    typedef bU (*Lfun)(const double t, const double x, double& gamma);
    typedef bU (*NLfun)(bU u, double t, double x);
    typedef vvector<bU> Sol;
    typedef vvector<double> VEC;
    u_int N_x;
    double t_start;
    double t_end;
    double x_start;
    double x_end;
    Lfun initial;
    double CFL;
    Sol Con;//数值解,mass_average conservative variables
    Sol Pri;//数值解,primitive variables
    Sol ReconL_Con, ReconR_Con;//重构mass_average守恒变量,x_{i+\frac12}处左右极限值
    Sol ReconL_Pri, ReconR_Pri;//重构原始变量,x_{i+\frac12}处左右极限值
    VEC Di;//mass on each cell
    VEC mesh;
    VEC Gamma;
    VEC us;
    VEC cs;
    Sol FLUX;
    GHOST ghostl, ghostr;

  public:
    Lagranian1D(u_int Nx, double t_start, double t_end, double x_start, double x_end,
        Lfun initial, double CFL = 0.5)
      : N_x(Nx), t_start(t_start), t_end(t_end), x_start(x_start), x_end(x_end),
      initial(initial), CFL(CFL) {
        Con.resize(N_x);
        Pri.resize(N_x);
        ReconL_Con.resize(N_x+1);//边界加上了
        ReconR_Con.resize(N_x+1);
        ReconL_Pri.resize(N_x+1);//边界加上了
        ReconR_Pri.resize(N_x+1);
        Di.resize(N_x);
        mesh.assign(N_x+1, 0);
        Gamma.assign(N_x, 0);
        us.assign(N_x+1, 0);
        cs.assign(N_x, 0);
        FLUX.resize(N_x+1);
        double h0 = (x_end - x_start) / N_x;
#pragma omp parallel for num_threads(Nthread)
        for(u_int i = 0; i < N_x+1; ++i) {
          mesh[i] = h0*i;
        }
      }

  public:
    void Initial() {
      for(u_int j = 0; j < N_x; ++j) {
        Pri[j] = initial(t_start, 0.5*(mesh[j]+mesh[j+1]), Gamma[j]);
        Con[j] = Pri2Con(Pri[j], Gamma[j]);
        Di[j]  = 1/Con[j][0]*(mesh[j+1] - mesh[j]);
      }
      //std::cout << "Di:\n"  << Di  << std::endl;
      //std::cout << "Con:\n" << Con << std::endl;
      //std::cout << "Pri:\n" << Pri << std::endl;
      //std::cout << "inverse:\n";
      //for(u_int j = 0; j < N_x; ++j) {
        //std::cout << Con2Pri(Con[j], Gamma[j]) << std::endl;
      //}
      //abort();
    }

    void InfiniteBD(Sol& Con) {
      ghostl.Con = Con[0];
      ghostl.Pri = Con2Pri(Con[0], Gamma[0]);
      ghostr.Con = Con[N_x-1];
      ghostr.Pri = Con2Pri(Con[N_x-1], Gamma[N_x-1]);
      ghostl.Gamma = Gamma[0];
      ghostr.Gamma = Gamma[N_x-1];
    }

    double fh(const bU& U, double h, const double Gamma);
    double fhp(const bU& U, double h, const double Gamma);
    /**
     * @brief Con2Pri solve (rho,u,p) from (D, m, E)
     *
     * @param U conservative variables (D, m, E)
     * @param gamma
     *
     * @return primitive variables (rho, u, p)
     */
    bU Con2Pri(const bU& U, const double Gamma);
    bU Pri2Con(const bU& U, const double Gamma);
    double cal_cs(const bU& Con, const bU& Pri, const double Gamma);
    void update_cs(VEC&);
    double cal_max_lambda_Lag(int i);
    double cal_max_lambda_Lag();
    double cal_max_lambda_Eul(int i);
    double cal_max_lambda_Lag(const bU& Con, const bU& Pri, const double Gamma);
    double cal_max_lambda_Eul(const bU& Con, const bU& Pri, const double Gamma);

    void cal_min_max_roe_lam(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double GAMMAL, const double GAMMAR,
        double& lam1, double& lam3);
    double cal_max_roe_lam_Lag(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double GAMMAL, const double GAMMAR);

    double t_step(const double CFL, double& alpha);

    bU F(const bU& CON, const bU& PRI);
    void add_BD_GAMMA(VEC& GAMMAL, VEC& GAMMAR);

    bU LF(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double alpha);
    void cal_flux_LF(Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri, Sol& FLUX, double alpha);

    bU LLF(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double alpha);
    void cal_flux_LLF(Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri, Sol& FLUX);

    bU HLLC(bU& CONL, bU& CONR, bU& PRIL, bU& PRIR, const double Gammal, const double Gammar, double&);
    void cal_flux_HLLC(Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri,
        Sol& FLUX, VEC& us);

    bU TS(bU& CONL, bU& CONR, bU& PRIL, bU& PRIR, const double Gammal, const double Gammar, double& us);
    void cal_flux_TS(Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri,
    Sol& FLUX, vvector<double>& us);

    void cal_us_roeav(Sol& ReconL_Pri, Sol& ReconR_Pri, const VEC&, const VEC&, VEC& us);
    void move_mesh(VEC&, VEC&, double dt, VEC&);

    void Reconstruction(const Sol& sol, const VEC& mesh,
        Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri);

    void ROE_AV_MAT(const bU& PRIL, const bU& PRIR, const double GAMMAL, const double GAMMAR,
        MAT& R, MAT& L);
    bU multiply(const bU& x, MAT& M);
    void Chara_MAT(const bU& PRIL, const double GAMMAL, MAT& R, MAT& L);

    void Euler_forward_LF(double dt, double alpha, VEC& mesh);
    void Euler_forward_LLF(double dt, VEC& mesh);
    void Euler_forward_HLLC(const double dt, VEC& mesh);
    void Euler_forward_TS(const double dt, VEC& mesh);

    void RK2_LF(Sol& Con, Sol& Pri, VEC& mesh, const double dt);
    void RK2_LLF(Sol& Con, Sol& Pri, VEC& mesh, const double dt);
    void RK2_HLLC(Sol& Con, Sol& Pri, VEC& mesh, const double dt);

    void SSP_RK_LF(Sol& Con, Sol& Pri, VEC& mesh, const double dt);
    void SSP_RK_LLF(Sol& Con, Sol& Pri, VEC& mesh, const double dt);
    void SSP_RK_HLLC(Sol& Con, Sol& Pri, VEC& mesh, const double dt);

  public:
    void Solver();

    void update_sol(VEC& mesh, Sol& Con, Sol& Pri, Sol& FLUX, const double dt,
        VEC& mesh1, Sol& Con1, Sol& Pri1);

    void print_con(std::ostream& os);
    void print_pri(std::ostream& os);
    void print_rupe(std::ostream& os);
    friend std::ostream& operator<<(std::ostream&, const Lagranian1D&);

};

#endif

