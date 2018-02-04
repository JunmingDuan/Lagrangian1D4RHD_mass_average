#include "Lagranian1D.h"

std::ostream& operator<<(std::ostream& os, const Lagranian1D& H) {
  os.setf(std::ios::scientific);
  os.precision(16);
  for(u_int i=0; i < H.N_x; i++) {
    os << 0.5*(H.mesh[i]+H.mesh[i+1]) << " " << H.Pri[i] << "\n";
  }
  return os;
}

void Lagranian1D::print_con(std::ostream& os) {
  os.setf(std::ios::scientific);
  os.precision(16);
  for(u_int i=0; i < N_x; i++) {
    os << 0.5*(mesh[i]+mesh[i+1]) << " " << Con[i] << "\n";
  }
}

void Lagranian1D::print_pri(std::ostream& os) {
  os.setf(std::ios::scientific);
  os.precision(16);
  for(u_int i=0; i < N_x; i++) {
    os << 0.5*(mesh[i]+mesh[i+1]) << " " << Pri[i] << "\n";
  }
}

void Lagranian1D::print_rupe(std::ostream& os) {
  os.setf(std::ios::scientific);
  os.precision(16);
  for(u_int i=0; i < N_x; i++) {
    os << 0.5*(mesh[i]+mesh[i+1]) << " "
      << Pri[i][0] << " " << Pri[i][1] << " " << Pri[i][2] << " "
      << Pri[i][2]/Pri[i][0]/(Gamma[i]-1) << "\n";
  }
}

double Lagranian1D::fh(const bU& U, double h, const double Gamma) {
  double h2 = h*h;
  double s2h2 = U[1]*U[1] + h2;
  return h2 + (Gamma-1)*h - Gamma*U[2]*sqrt(s2h2) + Gamma*U[1]*U[1];
}

double Lagranian1D::fhp(const bU& U, double h, const double Gamma) {
  double h2 = h*h;
  double s2h2 = U[1]*U[1] + h2;
  return 2*h + (Gamma-1) - Gamma*U[2]*h/sqrt(s2h2);
}

//solve a nonlinear equation by Newton method to obtain h
bU Lagranian1D::Con2Pri(const bU& U, const double Gamma) {
  bU prim;
  u_int ite(0), MAXITE(10);
  double eps = 1e-15;
  double a = 1, b = Gamma*U[2];
  double h(0.5*(1.+b)), h1(h), y(fh(U,h, Gamma));
  while(fabs(y) > eps && ite < MAXITE) {
    h1 = h - y/fhp(U, h, Gamma);
    if(fabs(h1-h) < eps) { h = h1; break; }
    if(h1 < a) h1 = a;
    if(h1 > b) h1 = b;
    y = fh(U, h1, Gamma);
    h = h1;
    ite++;
  }
  if(h < 1 || h != h) {
    std::cout << U[0] << " " << U[1] << " " << U[2] << std::endl;
    std::cout << "ite: " << ite << " ; h: " << h << " ; y: " << y << std::endl;
    std::cout << "derivative: " << fhp(U, 0.5*b, Gamma) << std::endl;
    //std::cout << "p1: " << 0.5*b - fp(U,0.5*b,Gamma)/fpp(U, 0.5*b, Gamma) << std::endl;
    std::cout << "f(p1): " << fh(U,h, Gamma) << std::endl;
    //std::cout << U[1]/(U[2]+h) << std::endl;
    abort();
  }

  prim[1] = U[1]/sqrt(h*h+U[1]*U[1]);
  double gamma2 = 1./(1-prim[1]*prim[1]);
  prim[0] = 1/U[0]/sqrt(gamma2);
  prim[2] = (h-1)*(Gamma-1)/Gamma*prim[0];
  //std::cout << "ite: " << ite << " ; final h: " << h << std::endl;
  //std::cout << "rho: " << prim[0] << std::endl;
  //abort();

  return prim;
}

bU Lagranian1D::Pri2Con(const bU& U, const double Gamma) {
  bU Con;
  double gamma = 1./sqrt(1-U[1]*U[1]);
  double h = 1 + U[2]/U[0]*Gamma/(Gamma-1);
  Con[0] = 1./(gamma*U[0]);
  Con[1] = h*gamma*U[1];
  Con[2] = h*gamma - U[2]*Con[0];

  return Con;
}

void Lagranian1D::update_cs(VEC& cs) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    cs[i] = cal_cs(Con[i], Pri[i], Gamma[i]);
  }
}

double Lagranian1D::cal_cs
(const bU& Con, const bU& Pri, const double Gamma) {
  double h = 1+Pri[2]/Pri[0]*Gamma/(Gamma-1);
  return sqrt(Gamma*Pri[2]/Pri[0]/h);
}

double Lagranian1D::cal_max_lambda_Lag(int i) {
  double u = Pri[i][1];
  return fabs(1-u*u)*cs[i]/(1+fabs(u)*cs[i]);
}

double Lagranian1D::cal_max_lambda_Eul(int i) {
  double u = Pri[i][1];
  return fabs(1-u*u)*cs[i]/(1+fabs(u)*cs[i])+fabs(u);
}

double Lagranian1D::cal_max_lambda_Lag
(const bU& Con, const bU& Pri, const double Gamma) {
  double h = 1+Pri[2]/Pri[0]*Gamma/(Gamma-1);
  double cs = sqrt(Gamma*Pri[2]/Pri[0]/h);
  double u = Pri[1];
  return fabs(1-u*u)*cs/(1+fabs(u)*cs);
}

double Lagranian1D::cal_max_lambda_Eul
(const bU& Con, const bU& Pri, const double Gamma) {
  double h = 1+Pri[2]/Pri[0]*Gamma/(Gamma-1);
  double cs = sqrt(Gamma*Pri[2]/Pri[0]/h);
  double u = Pri[1];
  return fabs(1-u*u)*cs/(1+fabs(u)*cs)+fabs(cs);
}

double Lagranian1D::cal_max_lambda_Lag() {
  double a(0);
  for(u_int i = 0; i < N_x; ++i) {
    a = std::max(a, cal_max_lambda_Lag(i));
  }
  return a;
}

bU Lagranian1D::F(const bU& CON, const bU& PRI) {
  bU tmp;
  tmp[0] = -PRI[1];
  tmp[1] = PRI[2];
  tmp[2] = PRI[1]*PRI[2];
  return tmp;
}

void Lagranian1D::update_sol(VEC& mesh, Sol& Con, Sol& Pri, Sol& FLUX, const double dt,
    VEC& mesh1, Sol& Con1, Sol& Pri1) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con1[i] = Con[i] - dt*(FLUX[i+1] - FLUX[i]) / Di[i];
    Pri1[i] = Con2Pri(Con1[i], Gamma[i]);
  }
}

double Lagranian1D::t_step(const double CFL, double& alpha) {
  double a(1), tmp_lam, hi, tmp_t;
  double roe_lam1, roe_lam3;
  alpha = 0;
  for(u_int i = 0; i < N_x; ++i) {
    hi = mesh[i+1] - mesh[i];
    tmp_lam = cal_max_lambda_Eul(i);
    alpha = std::max(alpha, tmp_lam);
    tmp_t = hi/alpha;
    a = std::min(a, tmp_t);
    //roe_average_characteristic_speed
    if(i != N_x-1) {
      cal_min_max_roe_lam(Con[i], Con[i+1], Pri[i], Pri[i+1], Gamma[i], Gamma[i+1], roe_lam1, roe_lam3);
      tmp_lam = std::max(tmp_lam, std::max(fabs(roe_lam1), fabs(roe_lam3)));
      tmp_t = hi/alpha;
      a = std::min(a, tmp_t);
    }
  }
  return CFL*a;
}

void Lagranian1D::cal_min_max_roe_lam(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double GAMMAL, const double GAMMAR,
    double& lam1, double& lam3) {
    double hl, hr, kl, kr, gl, gr, ul, ur;
    double v0, v1, v2;
    hl = 1 + GAMMAL/(GAMMAL-1)*PRIL[2]/PRIL[0];
    hr = 1 + GAMMAR/(GAMMAR-1)*PRIR[2]/PRIR[0];
    kl = sqrt(PRIL[0]*hl);
    kr = sqrt(PRIR[0]*hr);
    ul = PRIL[1];
    ur = PRIR[1];
    gl = 1./sqrt(1-ul*ul);
    gr = 1./sqrt(1-ur*ur);
    v0 = (kl*gl + kr*gr)/(kl+kr);
    v1 = (kl*gl*ul + kr*gr*ur)/(kl+kr);
    v2 = (kl*PRIL[2]/PRIL[0]/hl + kr*PRIR[2]/PRIR[0]/hr) / (kl+kr);
    double va = - v0*v0 + v1*v1;
    double s2 = 0.5*GAMMAL*v2*(1-va) - 0.5*(GAMMAL-1)*(1+va);
    double s = sqrt(s2);
    double e = -va;
    double y = sqrt((1-GAMMAL*v2)*e + s2);

    lam1 = ((1-GAMMAL*v2)*v0*v1 - s*y) / ((1-GAMMAL*v2)*v0*v0 + s2);
    lam3 = ((1-GAMMAR*v2)*v0*v1 + s*y) / ((1-GAMMAL*v2)*v0*v0 + s2);
}

double Lagranian1D::cal_max_roe_lam_Lag(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double GAMMAL, const double GAMMAR) {
    double hl, hr, kl, kr, gl, gr, ul, ur;
    double v0, v1, v2;
    hl = 1 + GAMMAL/(GAMMAL-1)*PRIL[2]/PRIL[0];
    hr = 1 + GAMMAR/(GAMMAR-1)*PRIR[2]/PRIR[0];
    kl = sqrt(PRIL[0]*hl);
    kr = sqrt(PRIR[0]*hr);
    ul = PRIL[1];
    ur = PRIR[1];
    gl = 1./sqrt(1-ul*ul);
    gr = 1./sqrt(1-ur*ur);
    v0 = (kl*gl + kr*gr)/(kl+kr);
    v1 = (kl*gl*ul + kr*gr*ur)/(kl+kr);
    v2 = (kl*PRIL[2]/PRIL[0]/hl + kr*PRIR[2]/PRIR[0]/hr) / (kl+kr);
    double va = - v0*v0 + v1*v1;
    double s2 = 0.5*GAMMAL*v2*(1-va) - 0.5*(GAMMAL-1)*(1+va);
    double s = sqrt(s2);
    double e = -va;
    double y = sqrt((1-GAMMAL*v2)*e + s2);

    double lam1 = ((1-GAMMAL*v2)*v0*v1 - s*y) / ((1-GAMMAL*v2)*v0*v0 + s2);
    double lam2 = v1/v0;
    double lam3 = ((1-GAMMAR*v2)*v0*v1 + s*y) / ((1-GAMMAL*v2)*v0*v0 + s2);
    return std::max(fabs(lam1-lam2), fabs(lam3-lam2));

}

