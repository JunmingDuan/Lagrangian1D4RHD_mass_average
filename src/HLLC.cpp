#include "Lagranian1D.h"

bU Lagranian1D::HLLC(bU& CONL, bU& CONR, bU& PRIL, bU& PRIR, const double Gammal, const double Gammar, double& SM) {
  double SL, SR;//, ui, ci, vl, vr;
  double hl, csl, ul, hr, csr, ur;
  double lam1, lam2, lam3, lam4;
  double roe_lam1, roe_lam3;

  cal_min_max_roe_lam(CONL, CONR, PRIL, PRIR, Gammal, Gammar, roe_lam1, roe_lam3);
  hl = 1+PRIL[2]/PRIL[0]*Gammal/(Gammal-1);
  csl = sqrt(Gammal*PRIL[2]/PRIL[0]/hl);
  ul = PRIL[1];
  lam1 = ul-(1-ul*ul)*csl/(1-ul*csl);//negative speed
  lam2 = ul+(1-ul*ul)*csl/(1+ul*csl);//positive speed
  hr = 1+PRIR[2]/PRIR[0]*Gammar/(Gammar-1);
  csr = sqrt(Gammar*PRIR[2]/PRIR[0]/hr);
  ur = PRIR[1];
  lam3 = ur-(1-ur*ur)*csr/(1-ur*csr);//negative speed
  lam4 = ur+(1-ur*ur)*csr/(1+ur*csr);//positive speed
  SL = std::min(std::min(lam1, lam3), roe_lam1);//left characteristic speed
  SR = std::max(std::max(lam2, lam4), roe_lam3);//right characteristic speed

  bU UL, UR;
  UL[0] = 1./CONL[0];
  UL[1] = CONL[1]/CONL[0];
  UL[2] = CONL[2]/CONL[0];
  UR[0] = 1./CONR[0];
  UR[1] = CONR[1]/CONR[0];
  UR[2] = CONR[2]/CONR[0];
  bU FL, FR;
  FL[0] = UL[0]*PRIL[1];
  FL[1] = UL[1]*PRIL[1] + PRIL[2];
  FL[2] = UL[1];
  FR[0] = UR[0]*PRIR[1];
  FR[1] = UR[1]*PRIR[1] + PRIR[2];
  FR[2] = UR[1];
  bU Uhll = (SR*UR - SL*UL + FL - FR) / (SR - SL);
  Uhll[2] /= Uhll[0];
  Uhll[1] /= Uhll[0];
  Uhll[0] = 1./Uhll[0];
  bU PRI_hll = Con2Pri(Uhll, Gammal);
  SM = PRI_hll[1];
  double PM = PRI_hll[2];
  std::cout << "speed: " << SL << " " << SM << " " << SR << std::endl;
  bU F;
  if(SL-SM > 0) {
    F[0] = 0;
    F[1] = SM;
    F[2] = SM*PM;
    return F;
  }
  else if(SR-SM < 0) {
    F[0] = 0;
    F[1] = SM;
    F[2] = SM*PM;
    return F;
  }
  else {
    F[0] = 0;
    F[1] = SM;
    F[2] = SM*PM;
    return F;
  }
}

/*bU Lagranian1D::HLLC(bU& CONL, bU& CONR, bU& PRIL, bU& PRIR, const double Gammal, const double Gammar, double& SM) {*/
  //double SL, SR;//, ui, ci, vl, vr;
  //double hl, csl, ul, hr, csr, ur;
  //double lam1, lam2, lam3, lam4;
  //double roe_lam1, roe_lam3;

  //cal_min_max_roe_lam(CONL, CONR, PRIL, PRIR, Gammal, Gammar, roe_lam1, roe_lam3);
  //hl = 1+PRIL[2]/PRIL[0]*Gammal/(Gammal-1);
  //csl = sqrt(Gammal*PRIL[2]/PRIL[0]/hl);
  //ul = PRIL[1];
  //lam1 = ul-(1-ul*ul)*csl/(1-ul*csl);//negative speed
  //lam2 = ul+(1-ul*ul)*csl/(1+ul*csl);//positive speed
  //hr = 1+PRIR[2]/PRIR[0]*Gammar/(Gammar-1);
  //csr = sqrt(Gammar*PRIR[2]/PRIR[0]/hr);
  //ur = PRIR[1];
  //lam3 = ur-(1-ur*ur)*csr/(1-ur*csr);//negative speed
  //lam4 = ur+(1-ur*ur)*csr/(1+ur*csr);//positive speed
  ////SL = std::min(lam1, lam3);//left characteristic speed
  ////SR = std::max(lam2, lam4);//right characteristic speed
  //SL = std::min(std::min(lam1, lam3), roe_lam1);//left characteristic speed
  //SR = std::max(std::max(lam2, lam4), roe_lam3);//right characteristic speed
  //{
    ////cal S^* and p^*, S^* is the smaller root of a quadric equation:
    ////F_E^{hll} x^2 - (E^{hll}+F^{hll}_m) x + m^{hll} = 0
    //double gl = 1./sqrt(1-ul*ul);
    //double gr = 1./sqrt(1-ur*ur);

    //double mL = PRIL[0]*hl*gl*gl*ul;
    //double mR = PRIR[0]*hr*gr*gr*ur;
    //double EL = PRIL[0]*hl*gl*gl - PRIL[2];
    //double ER = PRIR[0]*hr*gr*gr - PRIR[2];

    //double coe1 = SR*mL - SL*mR + SL*SR*(ER - EL);
    //double coe2 = ( SR*(mL*PRIL[1]+PRIL[2]) - SL*(mR*PRIR[1]+PRIR[2]) + SL*SR*(mR - mL)
        //+ SR*ER - SL*EL + mL - mR );
    //double coe3 = SR*mR - SL*mL + mL*PRIL[1]+PRIL[2] - mR*PRIR[1]-PRIR[2];
    //double PM;
    //if(fabs(coe1) < 1e-15) SM = coe3/coe2;
    //else SM = (coe2 - sqrt(coe2*coe2 - 4.*coe1*coe3))/2./coe1;
    //PM = (SM*(SL*EL-mL) + PRIL[2] - mL*(SL-PRIL[1])) / (1-SL*SM);
    //bU F;
    //F[0] = 0;
    //F[1] = PM;
    //F[2] = PM*SM;
    //return F;
  //}
//}

void Lagranian1D::cal_flux_HLLC(Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri,
    Sol& FLUX, vvector<double>& us) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    if(i == 0) FLUX[i] = HLLC(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i],
        Gamma[0], Gamma[i], us[i]);
    else if(i == N_x) FLUX[i] = HLLC(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i],
        Gamma[i-1], Gamma[N_x-1], us[i]);
    else FLUX[i] = HLLC(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i],
        Gamma[i-1], Gamma[i], us[i]);
  }
}

