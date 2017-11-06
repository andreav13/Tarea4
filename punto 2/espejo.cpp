#include <iostream>
#include <fstream>
#include <cmath> 
#include "Random64.h"

using namespace std;

const int Lx=600;
const int Ly=200;
const int Q=5;
const double W0=1/3.;
const double tau=0.5;
const double Utau=1./tau;
const double UmUtau=1-Utau;


class LatticeBoltzmann{
private:
  double w[Q];
  int V[2][Q]; //V[alpha][i] alpha=1 es x, alpha=0 es y
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; //f[ix][iy][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, bool UseNew, double sigma);
  double Jx(int ix, int iy);
  double Jy(int ix, int iy);
  double Ccelda(int ix, int iy);
  double feq(int i, int ix, int iy, double rho0, double Jx0, double Jy0);
  void Inicie(double rho0, double Jx0, double Jy0);
  double GetSigma(int ix, int iy, int t);
  void Colisione(int t);
  void Adveccione(void);
  void Imprimase(char const * NombreArchivo, int t);
  void ImprimaUnaLinea(char const * NombreArchivo,int t);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  w[0]=W0;
  w[1]=w[2]=w[3]=w[4]=1/6.;

  V[0][0]=0;
  V[1][0]=0;

  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1;  V[0][4]=0;
  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;   V[1][4]=-1;
  
}

double LatticeBoltzmann::rho(int ix, int iy, bool UseNew, double sigma){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
    suma+=fnew[ix][iy][i];
    else
      suma+=f[ix][iy][i];
  return suma+1/.2*sigma;
}

double LatticeBoltzmann::Jx(int ix, int iy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    suma+=V[0][i]*f[ix][iy][i];
  return suma;
}

double LatticeBoltzmann::Jy(int ix, int iy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    suma+=V[1][i]*f[ix][iy][i];
  return suma;
}

double LatticeBoltzmann::Ccelda(int ix, int iy){
  return 0.5;
}

double LatticeBoltzmann::feq(int i, int ix, int iy, double rho0, double Jx0, double Jy0){
  if(i==0)
    return (1-3*Ccelda(ix,iy)*Ccelda(ix,iy)*(1-W0))*rho0;
  else
    return w[i]*(3*Ccelda(ix,iy)*Ccelda(ix,iy)*rho0+3*(V[0][i]*Jx0+V[1][i]*Jy0));
  
}

void LatticeBoltzmann::Inicie(double rho0, double Jx0, double Jy0){
int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){
	f[ix][iy][i]=feq(i,ix,iy,rho0,Jx0,Jy0);
      }
}

double LatticeBoltzmann::GetSigma(int ix, int iy, int t){
  double A=10, lambda=10, omega=2*M_PI*Ccelda(ix,iy)/lambda;
  if(ix==0)
    return -A/omega*cos(omega*t);
  else
    return 0;
}

void LatticeBoltzmann::Colisione(int t){ //de f a fnew
  int ix,iy,i; double rho0,Jx0,Jy0,sigma;
  double A=10, lambda=10, omega=2*M_PI*Ccelda(ix,iy)/lambda;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){ //para cada celda

      if(ix>=100 and ix<=200 and (ix-50)*(ix-50)+(iy-100)*(iy-100)>100*100 )
	//condicion de espejo
      else{
	sigma=GetSigma(ix,iy,t);
	rho0=rho(ix,iy,false,sigma); Jx0=Jx(ix,iy); Jy0=Jy(ix,iy); //calculo campos
	//ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
	for(i=0;i<Q;i++) //para cada direccion
	  fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(i,ix,iy,rho0,Jx0,Jy0); //evoluciono
      }
    }
}
   
void LatticeBoltzmann::Adveccione(void){ //de fnew a f
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}

void LatticeBoltzmann::Imprimase(char const * NombreArchivo, int t){
  double rho0,Jx0,Jy0;
  ofstream MiArchivo(NombreArchivo);
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      double sigma = GetSigma(ix,iy,t);
      rho0=rho(ix,iy,true,sigma); Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy);
      //ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
      MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
    }MiArchivo<<endl;
  }
  MiArchivo.close();
}

void LatticeBoltzmann::ImprimaUnaLinea(char const * NombreArchivo,int t){
  ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0;
  int ix=Lx/2;
  for(int iy=0;iy<Ly;iy++){
    double sigma = GetSigma(ix,iy,t);
    rho0=rho(ix,iy,true,sigma);   Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy);
    //ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
    MiArchivo<<iy<<" "<<rho0<<endl;
  }
  MiArchivo.close();
}


//-----------------------Funciones Globales-----------------------


int main(void){

  LatticeBoltzmann Ondas;
  int t,tmax=400;

  double rho0=0,Jx0=0,Jy0=0;
  
  //Inicie
  Ondas.Inicie(rho0,Jx0,Jy0);

  //Corra
  for(t=0;t<=tmax;t++){
    Ondas.Colisione(t);
    Ondas.Adveccione();
  }
  
  Ondas.Imprimase("espejo.dat", t);
  
  
  return 0;
}
