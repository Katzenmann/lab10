#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
//---------------------------------------
using namespace std;
//---------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N, double);
void initialize(double* const u1, double* const u0, const double dx,const double dt, const double xmin,
                const int N);
void step(double* const u2, const double* const u1,const double* const u0,
          const double dt, const double dx, const int N);

//---------------------------------------
int main(){

  const double tEnd = 0.15 ;

  //const double M_PI=3.141;
  const int N  = 64;
  const double xmin = 0;
  const double xmax = 1;
  const double dx = (xmax-xmin)/N ;
  double dt = 0.001;
  double t = 0;
  const int Na = 10;
  const int Nk = int(tEnd/Na/dt);

  double* u0 = new double[N];
  double* u1 = new double[N];
  double* u2 = new double[N];
  double* h;
  stringstream strm;

  initialize(u1,u0,dx,dt, xmin,N);

  writeToFile(u0, "u_0", dx, xmin, N, t);

  cout << "Nk = " << Nk << endl;

  
  for(int i=1; i<=Na; i++)
  {
   for(int j=0; j<Nk; j++){

	  // step + swap here
     step(u2, u1, u0, dt, dx, N);
  double* temp=u0;
  u0=u1; 
  u1=u2;
  u2=temp;
  
      t +=dt;
   }
   strm.str("");
   strm << "u_" << i;
   writeToFile(u0, strm.str(), dx, xmin, N, t);
  }

  cout << "t = " << t << endl;

  delete[] u0;
  delete[] u1;
  delete[] u2;
  return 0;
}
//-----------------------------------------------
void step(double* const u2, const double* const u1,const double* const u0,
          const double dt, const double dx, const int N)
{
//zeitlich
  u2[0]=u0[0]-u1[0]*dt/dx*(u1[1]-u1[N-1]);
for(int i=1; i<N-1; i++){
 
  u2[i]=u0[i]-u1[i]*dt/dx*(u1[i+1]-u1[i-1]);
  
}
  u2[N-1]=u0[N-1]-u1[N-1]*dt/dx*(u1[0]-u1[N-2]);
}
//-----------------------------------------------
void initialize(double* const u1, double* const u0, const double dx,
                const double dt, const double xmin,  const int N)
{
  //räumlich
   double u,ux, uxx;
   for(int i=0; i<N; i++)
   {
     
     double x = xmin + i*dx;
     u1[i]=sin(2*M_PI*x);
     u0[i]=u1[i]+2*M_PI*dt*u1[i]*cos(2*M_PI*x)+1/2*dt*dt*(2*u1[1]*4*M_PI*M_PI*cos(2*M_PI*x)*cos(2*M_PI*x)-u1[i]*u1[i]*4*M_PI*M_PI*u1[i]);
   }
}
//-----------------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N, double t)
{
   ofstream out(s.c_str());
   for(int i=0; i<N; i++){
     double x = xmin + i * dx;
     out << x << "\t" << u[i] <<"\t"<< x+sin(2*M_PI*x)*t<< "\t"<<sin(2*M_PI*x)<< endl;
   }
   out.close();
}
