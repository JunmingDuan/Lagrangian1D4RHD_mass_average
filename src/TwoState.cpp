#include "Lagranian1D.h"

bU Lagranian1D::TS(bU& CONL, bU& CONR, bU& PRIL, bU& PRIR, const double Gammal, const double Gammar, double& us) {
  double pl = PRIL[2];
  double pr = PRIR[2];
  double ul = PRIL[1];
  double ur = PRIR[1];
  double El = CONL[2];
  double Er = CONR[2];
  double al = cal_cs(CONL, PRIL, Gammal);
  double ar = cal_cs(CONR, PRIR, Gammar);
  double zl = 1./CONL[0] * al;
  double zr = 1./CONR[0] * ar;
  double ps;

  double coe1 = Er*ul*zr - El*ur*zl - pl*ul*ur + pr*ul*ur - (Er*zl*zr)*CONL[0] + (El*zl*zr)*CONR[0]
    - (pr*ur*zl)*CONL[0] + (pl*ul*zr)*CONR[0];
  double coe2 = El*zl - Er*zr + pl*ul + pl*ur - pr*ul - pr*ur + CONL[1]*ur*zl - CONR[1]*ul*zr
    + pr*zl*CONL[0] - pl*zr*CONR[0] + CONR[1]*zl*zr*CONL[0] - CONL[1]*zl*zr*CONR[0];
  double coe3 = pr - pl - CONL[1]*zl + CONR[1]*zr;
  double us1(10), us2(10);
  if(fabs(coe1) < 1e-15) us = coe3/coe2;
  else {
    us1 = (-coe2 - sqrt(coe2*coe2 - 4.*coe1*coe3))/2./coe1;
    us2 = (-coe2 + sqrt(coe2*coe2 - 4.*coe1*coe3))/2./coe1;
  }
  if(fabs(us1) <= 1) us = us1;
  else if(fabs(us2) <= 1) us = us2;
  else us = 100;
  ps = ((El*zl + pl*ul)*us - CONL[1]*zl - pl) / ((ul - zl*CONL[0])*us - 1);

  std::cout << ul << " " << ur << std::endl;
  std::cout << zl << " " << zr << std::endl;
  std::cout << pl << " " << pr << std::endl;
  std::cout << CONL[0] << " " << CONR[0] << std::endl;
  std::cout << "coe123: " << coe1 << " " << coe2 << " " << coe3 << std::endl;
  std::cout << "us1, us2: " << us1 << " " << us2 << std::endl;
  std::cout << "us, ps: " << us << " " << ps << std::endl;

  bU F;
  F[0] = -us;
  F[1] = ps;
  F[2] = ps*us;
  return F;
}

void Lagranian1D::cal_flux_TS(Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri,
    Sol& FLUX, vvector<double>& us) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    if(i == 0) FLUX[i] = TS(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i],
        Gamma[0], Gamma[i], us[i]);
    else if(i == N_x) FLUX[i] = TS(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i],
        Gamma[i-1], Gamma[N_x-1], us[i]);
    else FLUX[i] = TS(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i],
        Gamma[i-1], Gamma[i], us[i]);
  }
}

