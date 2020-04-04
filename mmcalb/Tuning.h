#ifndef Tuning_h
#define Tuning_h 1

#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"

//=================================================//
//=============     Tuning     ====================//
//=================================================//

class Tuning{

 public:
  Tuning();
  ~Tuning();
  double calcf2t_zt(double * P,
		    double xf, double xpf,
		    double yf, double ypf);

  double calcf2t_ang(double* P,
			    double xf, double xpf, 
			    double yf, double fpf,
			    double zt);  

  double calcRasterCor(double a, double b, double c);

  //--- calcf2t_zt --------//
  const int nMatTz=nnz; 
  const int nXf=nn;
  const int nXpf=nn;
  const int nYf=nn;
  const int nYpf=nn;
  double Y,x; 
  int npar;
  int a,b,c,d,e;
  
  //---- calcf2t_ang -----//
  const int nMatT=nn;
  const int nZt=nn;

  //--- calcRasterCor ----//
  
};

Tuning::Tuning(){};
Tuning::~Tuning(){};



double Tuning::calcf2t_zt(double* P, double xf, double xpf, 
                 double yf, double ypf){


  
  Y=0.;
  x=1.; 
  npar=0;
  a=0,b=0,c=0,d=0;
  
  for (int n=0;n<nMatTz+1;n++){
  	for (d=0;d<n+1;d++){
	  for (c=0;c<n+1;c++){
	    for (b=0;b<n+1;b++){
	      for (a=0;a<n+1;a++){
		
		if (a+b+c+d==n){
		  if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf){
		    x = pow(xf,double(a))*pow(xpf,double(b))*
		      pow(yf,double(c))*pow(ypf,double(d));
		  }
		  else{
		    x = 0.;
		  }
		  Y += x*P[npar];
		  npar++;
		}
		
	      }
	    }
	  }
  	}
  }
  
  return Y;
}

//============ calcf2t_ang ============//

 double Tuning::calcf2t_ang(double* P, double xf, double xpf, 
			      double yf, double ypf, double zt)

{
  
  // ------------------------------------------------ //
  // ----- 4rd order using xf, xpf, yf, ypf, zt ----- //
  // ------------------------------------------------ //

  
  Y=0.;
  x=1.; 
  npar=0;
  a=0,b=0,c=0,d=0,e=0;
  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));

		}
		else{
		  x = 0.;
		}
		Y += x*P[npar]; 
	      npar++;
	      }
	      
	    }
	  }
	}
      }    
    }
  }

  return Y; 
  
};

double Tuning::calcRasterCor(double a, double b, double c){
  return a * b + c;
}


#endif
