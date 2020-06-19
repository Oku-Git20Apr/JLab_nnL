#include <iostream>
#include <fstream>
using namespace std;
#include "TApplication.h"
#include "eff_z2d.h"
#include "Param.h"
#include "Tree.h"
#include "TMath.h"

double s2f1_off(int i,char const* ARM,char const* MODE, int KINE);
double Calc_ras(double a,double b,double c){return  a *b + c;};  
double calcf2t_ang(double* P,double xf, double xpf, double yf, double fpf,double z);
double calcf2t_zt(double* P, double xf, double xpf, double yf, double ypf);
double calcf2t_mom(double* P, double xf, double xpf, double yf, double ypf, double zt);
double Num_Al(double a);
double calcf_pathl(double* P, double xf, double xpf, double yf, double ypf, double zt);
double calcf2t_3rd(double* P, double xf, double xpf, double yf, double ypf, double zt);

double F_Voigt( double *x, double *par )
  {
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
    // par[3] : lorentz fwhm
    double val = par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4);
    return val;
  }
double F_mmnoAC( double *x, double *par )// Npar=10
  {
    // par[0] : scale 
    // par[1] : location
    // par[2] : gaussian sigma
    // par[3] : lorentz fwhm
    double val = par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4) + par[4]*TMath::Gaus(x[0],par[5],par[6],1)+par[7]*TMath::Gaus(x[0],par[8],par[9],1);
    return val;
  }


// #################################################
double s2f1_off(int i,char const* ARM,char const* MODE, int KINE){


  double RS2_offset[16],LS2_offset[16];
  if(*MODE=='H' && KINE==2){
 
 double  RS2_off_H2[16]={-16911.4,-16864.3,-16900,-16897,-16873.8,-16868.4,-16901.1,-16876.8,-16895.4,-16860.9,-16893.1,-16884.4,-16847.3,-16842.7,-16836.9,-16882.6};
 double  LS2_off_H2[16]={-25336.9,-25386.6,-25367.5,-25392.3,-25391.1,-25386.2,-25422,-25428.9,-25417.3,-25426.8,-25438.7,-25383.4,-25396,-25418.5,-25436.4,-26082.1};
 
  LS2_offset[i]=LS2_off_H2[i];
  RS2_offset[i]=RS2_off_H2[i];
  }


  if(*MODE=='H' && KINE==1){
    
    //double  RS2_off_H1[16]={-16911.4,-16864.9,-16900,-16897.6,-16874.8,-16869.3,-16901.1,-16876.8,-16895.6,-16860.3,-16892.6,-16885,-16847.3,-16843.3,-16838.4,-16882.6};
    //double  LS2_off_H1[16]={-25336.9,-25385.7,-25367,-25392.2,-25391,-25386.3,-25422,-25428.9,-25415.2,-25425,-25438,-25381,-25394.4,-25417.5,-25432.8,-26082.1};

double  RS2_off_H1[16]={-16828.7,-16863,-16894,-16893.3,-16870.9,-16867.2,-16900.3,-16876.8,-16895.6,-16861.6,-16895,-16890.7,-16854.6,-16852.9,-16850.5,-16861.9};
double  LS2_off_H1[16]={-25335,-25385.6,-25367,-25392.1,-25391.7,-25386.4,-25422.1,-25428.9,-25414.9,-25424.7,-25436.9, -25381.2,-25390,-25413.4,-25428.7,-26640.8};
  LS2_offset[i]=LS2_off_H1[i];
  RS2_offset[i]=RS2_off_H1[i];
  }

 double s2f1_offset; 
 if(*ARM=='R')s2f1_offset=RS2_offset[i];
 else  if(*ARM=='L')s2f1_offset=LS2_offset[i];
 else {s2f1_offset=0.;cout<<"false read out !!"<<endl;}

  return s2f1_offset;

}



// ###################################################
double calcf2t_mom(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt){
  // ----------------------------------------------------------------- //
  // ------ 4rd order using xf, xpf, yf, ypf, zt, xt, xpt, yt, ytp --- //
  // ----------------------------------------------------------------- //

  const int nMatT=nnp;  
  const int nXf=nnp;
  const int nXpf=nnp;
  const int nYf=nnp;
  const int nYpf=nnp;
  const int nZt=nnp;
  //  const int nXt=nnp;
  //  const int nXpt=nnp;
  //  const int nYt=nnp;
  //  const int nYpt=nnp;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
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
  
}



// ####################################################
double calcf2t_ang(double* P, double xf, double xpf, 
		     double yf, double ypf, double zt){
// ####################################################

  // ------------------------------------------------ //
  // ----- 4rd order using xf, xpf, yf, ypf, zt ----- //
  // ------------------------------------------------ //
  
  const int nMatT=nn;  
  const int nXf=nn;
  const int nXpf=nn;
  const int nYf=nn;
  const int nYpf=nn;
  const int nZt=nn;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
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
  
}
      	


// #################################################
double calcf2t_zt(double* P, double xf, double xpf, 
                 double yf, double ypf){
// ###############################################

  int nnz=3;  

  const int nMatT=nnz; 
  const int nXf=nnz;
  const int nXpf=nnz;
  const int nYf=nnz;
  const int nYpf=nnz;

  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0;
  
  for (int n=0;n<nMatT+1;n++){
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



// #################################################
void tuning::GetACParam(){

  cout<<"==================================="<<endl;
  cout<<"========== GetACParam ============="<<endl;
  cout<<"==================================="<<endl;
  // taken by /ac/param/offset_ac.dat 
  string pname="../ac/param/offset_ac.dat";
  ifstream ifp(pname.c_str(),ios::in);
  if (ifp.fail()){ cerr << "failed open files" <<pname.c_str()<<endl; exit(1);}
  cout<<" Param file : "<<pname.c_str()<<endl;
  
  string buf;
  int AC,Seg;
  double off,pe;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >> AC >> Seg >> off >> pe;

    if(AC==1){
      ac1_off[Seg]=off;
      ac1_1pe[Seg]=pe;
    }else if(AC==2){
      ac2_off[Seg]=off;
      ac2_1pe[Seg]=pe;
    }else{
      cout<<"Error :"<<endl; exit(1);
    }

  }

  }



// #################################################
double tuning::AC_npe(int nac, int seg, double adc){



  double npe,ac_off,ac_1pe;

  if(nac==1){
    ac_off=ac1_off[seg];
    ac_1pe=ac1_1pe[seg];
  }else if(nac==2){
    ac_off=ac2_off[seg];
    ac_1pe=ac2_1pe[seg];
  }else {
    cout<<"Error : falid Get AC parameters "<<endl; exit(1);}

  //  npe=(adc)/(ac_1pe - ac_off); // Just correct gain
  if(nac==1)npe=(adc)/(ac_1pe - ac_off)*2.0;     // Gogami AC DB was changed gain 400 -> 200
  else if(nac==2)  npe=(adc)/(ac_1pe - ac_off);
    // in this case, we need scale gain parameter 2 times
  return npe;  
}



// #################################################
double calcf_pathl(double* P, double xf, double xpf, double yf, double ypf, double zt){
// #################################################



  const int nMatT=nnc; 
  const int nXf=nnc;
  const int nXpf=nnc;
  const int nYf=nnc;
  const int nYpf=nnc;
  const int nZt=nnc;

  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  
  for (int n=0;n<nMatT+1;n++){
    for (e=0;e<n+1;e++){
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
}


//////////////////////////////////////////////////
double calcf2t_3rd(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt){
  // ------------------------------------------------ //
  // ----- 3rd order using xf, xpf, yf, ypf, zt ----- //
  // ------------------------------------------------ //
  const int nMatT=3;  
  const int nXf=3;
  const int nXpf=3;
  const int nYf=3;
  const int nYpf=3;
  const int nZt=3;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  
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
  
}

///////////////////////////////////////////////////////////////////////////
void tuning::matrix(string mtparam){

  cout<<endl;
  cout<<"==============================="<<endl;
  cout<<"=== Input Matrix Parameters ==="<<endl;
  cout<<"==============================="<<endl;

  string buf;
  int s=0;
  ifstream ifp(Form("%s",mtparam.c_str()),ios::in);
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >>param_mt[s];
    cout<<param_mt[s]<<endl;
    s++;
  }

  for(int i=0;i<12;i++)MT_p[i]=false;

  //======= Tuning selection flag =====================//
  //--------- RHRS ------------------------//
  MT_p[0] = true;  // RHRS z correction
  MT_p[1] = true;  // RHRS raster correction
  MT_p[2] = true;  // RHRS theta correction
  MT_p[3] = true;  // RHRS phi correction
  //--------- LHRS -----------------------//
  MT_p[4] = true;  // LHRS z correction
  MT_p[5] = true;  // LHRS raster correction
  MT_p[6] = true;  // LHRS theta correction
  MT_p[7] = true;  // LHRS phi correction
  //-------- momentum calibration ---------//
  MT_p[8] = true; // RHRS momentum correction  
  MT_p[9] = true; // LHRS momentum correction  
  ploss = true;  // Energy Loss
  //================================================//

  //  MT_p[10] = true; // RHRS path length correction  
  //  MT_p[11] = true; // LHRS path length correction
  MT_p[11] = false; // LHRS path length correction
  MT_p[10] = false; // RHRS path length correction  
  cout<<endl;
  MTParam_R();cout<<" Input RHRS Matrix parameter "<<endl;
  MTParam_L();cout<<" Input LHRS Matrix parameter "<<endl;
  MTParam_G();cout<<"Input Gogami parameter "<<endl;
  MTP_mom();cout<<"Input Mom parameter "<<endl;


  cout<<endl;
  
  cout<<"======== Correction Parameters ========="<<endl;
  if(MT_p[0])cout<<" RHRS z      correction "<<endl;
  else     cout<<" RHRS z                    no correction "<<endl;
  if(MT_p[1])cout<<" RHRS raster correction "<<endl;
  else     cout<<" RHRS raster               no correction "<<endl;
  if(MT_p[2])cout<<" RHRS theta  correction "<<endl;
  else     cout<<" RHRS theta                no correction "<<endl;
  if(MT_p[3])cout<<" RHRS phi    correction "<<endl;
  else     cout<<" RHRS phi                  no correction "<<endl;
  if(MT_p[4])cout<<" LHRS z      correction "<<endl;
  else     cout<<" LHRS z                    no correction "<<endl;
  if(MT_p[5])cout<<" LHRS raster correction "<<endl;
  else     cout<<" LHRS raster               no correction "<<endl;
  if(MT_p[6])cout<<" LHRS theta  correction "<<endl;
  else     cout<<" LHRS theta                no correction "<<endl;
  if(MT_p[7])cout<<" LHRS phi    correction "<<endl;
  else     cout<<" LHRS phi                  no correction "<<endl;
  if(MT_p[8])cout<<" RHRS mom    correction "<<endl;
  else     cout<<" RHRS mom                  no correction "<<endl;
  if(MT_p[9])cout<<" LHRS mom    correction "<<endl;
  else     cout<<" LHRS mom                  no correction "<<endl;
  if(ploss)cout<<" Energy Los  correction "<<endl;
  else     cout<<" Energy Los                no correction "<<endl;
 if(MT_p[10])cout<<" RHRS PathL  correction "<<endl;
  else     cout<<" RHRS PathL                no correction "<<endl;
 if(MT_p[11])cout<<" LHRS PathL  correction "<<endl;
  else     cout<<" LHRS PathL                no correction "<<endl;

  cout<<endl;
  
}

///////////////////////////////////////////////////////////////////////////

void tuning::MTParam_R(){

  //=================//
  //==== RHRS =======//
  //=================//


  //====== RHRS z parameters ======//

    char name_Mzt[500];
    sprintf(name_Mzt, param_mt[0].c_str()); // optimized
    ifstream Mzt(name_Mzt);
   if (Mzt.fail()){ cerr << "failed open files" <<name_Mzt<<endl; exit(1);}
   for(int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt >> par >> p >> p >> p >> p;
    Pzt[i]=par;
    //    cout<<"R Mzt : "<<Pzt[i]<<endl;
   }
  Mzt.close();

  
  //====== RHRS raster paramters =======//
    char name_Mras[500];
    sprintf(name_Mras, param_mt[1].c_str()); // optimized
    //    cout<<"RHRS Raster parameters file: "<<name_Mras<<endl;
  ifstream Mras(name_Mras);
   if (Mras.fail()){ cerr << "failed open files " <<name_Mras<<endl; exit(1);}
  for (int i=0;i<nParamT_ras;i++){

    Mras >> Pras[i];
    // ---- raster x is tuned by this macro ---- //
    if(i==1 || i==3) Pras[i] = 0.0;
  }


  Mras.close();    

  
  //===== RHRS theta parameters ======// 
    char name_Mxpt[500];
    sprintf(name_Mxpt, param_mt[2].c_str()); // optimized
  ifstream Mxpt(name_Mxpt);
   if (Mxpt.fail()){ cerr << "failed open files " <<name_Mxpt<<endl; exit(1);}
    for(int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt >> par >> p >> p >> p >> p >> p;
    Pxpt[i]  = par;
  }
  Mxpt.close();  
 //===== RHRS phi parameters ======//
  char name_Mypt[500];
    sprintf(name_Mypt, param_mt[3].c_str()); // optimized  
    ifstream Mypt(name_Mypt);
    if(Mypt.fail()){ cerr << "failed open files " <<name_Mypt<<endl; exit(1);}
    for (int i=0;i<nParamT;i++){
      double par=0.;
      int p=0;
      Mypt >> par >> p >> p >> p >> p >> p;
      Pypt[i]  = par;
    }
  Mypt.close();    

 //===== RHRS Path Length  parameters ======//

  if(MT_p[10]){
  char name_Mpl[500];
    sprintf(name_Mpl, param_mt[10].c_str()); // optimized  
    ifstream Mpl(name_Mpl);
    if(Mpl.fail()){ cerr << "failed open files " <<name_Mpl<<endl; exit(1);}
    for (int i=0;i<nParamTc;i++){
      double par=0.;
      int p=0;
      Mpl >> par >> p >> p >> p >> p >> p;
      Ppl[i]  = par;
    }
  Mpl.close();    
  }
  
};
//////////////////////////////////////////////////////////////

void tuning::MTParam_L(){

  //=================//
  //===== LHRS ======//
  //=================//

  
  //====== LHRS z parameters ======//  
  char name_Mzt_L[500];
  sprintf(name_Mzt_L,param_mt[4].c_str()); // optimized
  ifstream Mzt_L(name_Mzt_L);
  if (Mzt_L.fail()){ cerr << "failed open files " <<name_Mzt_L<<endl; exit(1);}
  for (int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt_L >> par >> p >> p >> p >> p;
    Pzt_L[i]=par;
  }
 Mzt_L.close();

  //====== LHRS raster paramters =======//
    char name_Mras_L[500];
    sprintf(name_Mras_L, param_mt[5].c_str()); // optimized
    //    cout<<"LHRS Raster parameters file: "<<name_Mras_L<<endl;
  ifstream Mras_L(name_Mras_L);
  if (Mras_L.fail()){ cerr << "failed open files " <<name_Mras_L<<endl; exit(1);}
  for (int i=0;i<nParamT_ras;i++){

    Mras_L >> Pras_L[i];
    // ---- raster x is tuned by this macro ---- //
    if(i==1 || i==3) Pras_L[i] = 0.0;
  }
  
  Mras_L.close();    

 
 //===== LHRS theta parameters ======// 
  char name_Mxpt_L[500];
    sprintf(name_Mxpt_L, param_mt[6].c_str()); // optimized
  ifstream Mxpt_L(name_Mxpt_L);
  if (Mxpt_L.fail()){ cerr << "failed open files " <<name_Mxpt_L<<endl; exit(1);}
  //  cout<<"LHRS theta parameters file: "<<name_Mxpt_L<<endl;  
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt_L >> par >> p >> p >> p >> p >> p;
    //   cout<<"LHRS theta : "<<par<<endl;
    Pxpt_L[i]  = par;
  }
  Mxpt_L.close();

  
 //===== LHRS phi parameters ===x==//
  char name_Mypt_L[500];
    sprintf(name_Mypt_L, param_mt[7].c_str()); // optimized
  ifstream Mypt_L(name_Mypt_L);
  if (Mypt_L.fail()){ cerr << "failed open files " <<name_Mypt_L<<endl; exit(1);}
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mypt_L >> par >> p >> p >> p >> p >> p;
    Pypt_L[i]  = par;    
  }
  Mypt_L.close();    


 //===== LHRS Path Length  parameters ======//

  if(MT_p[11]){
  char name_Mpl_L[500];
    sprintf(name_Mpl_L, param_mt[11].c_str()); // optimized  
    ifstream Mpl_L(name_Mpl_L);
    if(Mpl_L.fail()){ cerr << "failed open files " <<name_Mpl_L<<endl; exit(1);}
    for (int i=0;i<nParamTc;i++){
      double par=0.;
      int p=0;
      Mpl_L >> par >> p >> p >> p >> p >> p;
      Ppl_L[i]  = par;
    }
  Mpl_L.close();      

  }

}

//========================================================//


void tuning::MTP_mom(){

  //====== RHRS Momentum parameters ========//
    char name_Mpt[500];
    sprintf(name_Mpt, param_mt[8].c_str()); // optimized
    ifstream Mpt(name_Mpt);
  if (Mpt.fail()){ cerr << "failed open files " <<name_Mpt<<endl; exit(1);}
   for(int i=0;i<nParamTp;i++){
    double par=0.;
    int p=0;
    Mpt >> par >> p >> p >> p >> p >> p;
    Prp[i]=par;
    Opt_par_R[i]=par;
    Opt_par[i]=par;
   }
  Mpt.close();

  
  //====== LHRS Momentum parameters ========//

    char name_Mpt_L[500];
    sprintf(name_Mpt_L, param_mt[9].c_str()); // optimized
    ifstream Mpt_L(name_Mpt_L);
  if (Mpt_L.fail()){ cerr << "failed open files " <<name_Mpt_L<<endl; exit(1);}
   for(int i=0;i<nParamTp;i++){
    double par=0.;
    int p=0;
    Mpt_L >> par >> p >> p >> p >> p >> p;
    Plp[i]=par;
    Opt_par_L[i]=par;
    Opt_par[i+nParamTp]=par;  // Both momentum paramters
   }
  Mpt_L.close();

  
}

///////////////////////////////////////////////////////////////////////////

void tuning::MTParam_G(){


  cout<<"================================"<<endl;
  cout<<"=======Gogami Param ============"<<endl;
  cout<<"================================"<<endl;

//  int nParamT_3=3;
  char name_MctimeL[100];
  char name_MctimeR[100];
  sprintf(name_MctimeL,"../goga_mac/Rootfiles/matrices/ctimeL.dat"); 
  sprintf(name_MctimeR,"../goga_mac/Rootfiles/matrices/ctimeR.dat"); 
  ifstream MctimeL(name_MctimeL);
  ifstream MctimeR(name_MctimeR);
  //  double PctimeL[nParamT_3];
  //  double PctimeR[nParamT_3];
  for (int i=0;i<nParamTc;i++){
    double par = 0.0;
    int p = 0;
    MctimeL >> par >> p >> p >> p >> p >> p; 
    PctimeL[i]=par;
    
    par = 0.0;
    p   = 0;
    MctimeR >> par >> p >> p >> p >> p >> p; 
    PctimeR[i]=par;
  }
  MctimeR.close();


}



///////////////////////////////////////////////////////////////////////////

void tuning::Calib(int rt, int lt ){



  //======= Nomalization ==================//
  R_tr_x[rt]    = (R_tr_x[rt]-XFPm)/XFPr;
  R_tr_th[rt]   = (R_tr_th[rt]-XpFPm)/XpFPr;
  R_tr_y[rt]    = (R_tr_y[rt]-YFPm)/YFPr;
  R_tr_ph[rt]   = (R_tr_ph[rt]-YpFPm)/YpFPr;
  R_tr_vz[rt]   = (R_tr_vz[rt]-Ztm)/Ztr;
  R_tr_tg_th[rt]= (R_tr_tg_th[rt] - Xptm)/Xptr;
  R_tr_tg_ph[rt]= (R_tr_tg_ph[rt] - Yptm)/Yptr;

  R_p = (R_p - PRm)/PRr;

  L_tr_x[lt]    = (L_tr_x[lt]-XFPm)/XFPr; 
  L_tr_th[lt]   = (L_tr_th[lt]-XpFPm)/XpFPr;
  L_tr_y[lt]    = (L_tr_y[lt]-YFPm)/YFPr;
  L_tr_vz[lt]   = (L_tr_vz[lt]-Ztm)/Ztr;
  L_tr_ph[lt]   = (L_tr_ph[lt]-YpFPm)/YpFPr;
  L_tr_tg_th[lt]= (L_tr_tg_th[lt] - Xptm)/Xptr;
  L_tr_tg_ph[lt]= (L_tr_tg_ph[lt] - Yptm)/Yptr;  

  L_p = (L_p - PLm)/PLr;

  //========================================//
  
  if(MT_p[0]) R_tr_vz[rt]   = calcf2t_zt(Pzt, R_tr_x[rt], R_tr_th[rt], R_tr_y[rt], R_tr_ph[rt]); // nomalized
  if(MT_p[4]) L_tr_vz[lt]   = calcf2t_zt(Pzt_L, L_tr_x[lt], L_tr_th[lt], L_tr_y[lt], L_tr_ph[lt]); //nomalized


    //======== Raster Correction ==========================//    

    RasterCor = Calc_ras(R_Ras_x, Pras[2], Pras[0]);
    RasterCor = RasterCor/tan(hrs_ang);
    
    R_tr_vz[rt]  = R_tr_vz[rt]*Ztr +Ztm; // scaled     
    if(MT_p[1])    R_tr_vz[rt]  = R_tr_vz[rt] + RasterCor; // correction
    R_tr_vz[rt]  = (R_tr_vz[rt]-Ztm)/Ztr;    // nomalization     
    RasterCor_L  = Calc_ras(L_Ras_x, Pras_L[2], Pras_L[0]);
    RasterCor_L  = RasterCor_L/tan(hrs_ang);
    L_tr_vz[lt]  = L_tr_vz[lt]*Ztr +Ztm;     // scaled
    if(MT_p[5])    L_tr_vz[lt]  = L_tr_vz[lt] + RasterCor_L;
    L_tr_vz[lt]  =  (L_tr_vz[lt]  -  Ztm)/Ztr;    // nomalization

    //====================================================//

    


    if(MT_p[2])    R_tr_tg_th[rt]  = calcf2t_ang(Pxpt,   R_tr_x[rt], R_tr_th[rt], R_tr_y[rt], R_tr_ph[rt],R_tr_vz[rt]); // nomalized
    if(MT_p[3])    R_tr_tg_ph[rt]  = calcf2t_ang(Pypt,   R_tr_x[rt], R_tr_th[rt], R_tr_y[rt], R_tr_ph[rt],R_tr_vz[rt]); // nomalized
    if(MT_p[6])    L_tr_tg_th[lt]  = calcf2t_ang(Pxpt_L, L_tr_x[lt], L_tr_th[lt], L_tr_y[lt], L_tr_ph[lt], L_tr_vz[lt]); // nomalized
    if(MT_p[7])    L_tr_tg_ph[lt]  = calcf2t_ang(Pypt_L, L_tr_x[lt], L_tr_th[lt], L_tr_y[lt], L_tr_ph[lt], L_tr_vz[lt]); // nomalized   

    if(MT_p[8])    R_p = calcf2t_mom(Opt_par_R, R_tr_x[rt], R_tr_th[rt], R_tr_y[rt], R_tr_ph[rt],R_tr_vz[rt]);
    if(MT_p[9])    L_p = calcf2t_mom(Opt_par_L, L_tr_x[lt], L_tr_th[lt], L_tr_y[lt], L_tr_ph[lt],L_tr_vz[lt]);

    
    //========== Scaled at FP ==================//
    R_tr_x[rt]  = R_tr_x[rt]  * XFPr + XFPm;
    R_tr_th[rt] = R_tr_th[rt] * XpFPr + XpFPm;
    R_tr_y[rt]  = R_tr_y[rt]  * YFPr + YFPm;
    R_tr_ph[rt] = R_tr_ph[rt] * YpFPr + YpFPm;

    L_tr_x[lt]  = L_tr_x[lt]  * XFPr + XFPm;
    L_tr_th[lt] = L_tr_th[lt] * XpFPr + XpFPm;
    L_tr_y[lt]  = L_tr_y[lt]  * YFPr + YFPm;
    L_tr_ph[lt] = L_tr_ph[lt] * YpFPr + YpFPm;    

    //=========== Scaled at Taget =============//

    R_tr_vz[rt]     = R_tr_vz[rt] * Ztr + Ztm; // scaled
    R_tr_tg_th[rt]  = R_tr_tg_th[rt] * Xptr + Xptm; // scaled
    R_tr_tg_ph[rt]  = R_tr_tg_ph[rt] * Yptr + Yptm; // scaled
    R_p             = R_p * PRr + PRm; // scaled
    L_tr_vz[lt]     = L_tr_vz[lt] * Ztr + Ztm; // scaled
    L_tr_tg_th[lt]  = L_tr_tg_th[lt] * Xptr + Xptm;  // scaled    
    L_tr_tg_ph[lt]  = L_tr_tg_ph[lt] * Yptr + Yptm;  // scaled    
    L_p             = L_p * PLr + PLm; // scaled    
    

    // Lp = 2.2 GeV mode //
    if(Lp_scale)L_p=2.21807/2.1*L_p;
    //L_p=2.2/2.1*L_p;


    //=========== Energy Loss ===================//
    B_p     = B_p + Eloss(0.0,0,"B");
    R_p     = R_p + Eloss(R_tr_tg_ph[rt],R_tr_vz[rt],"R");
    L_p     = L_p + Eloss(L_tr_tg_ph[lt],L_tr_vz[lt],"L");


    
}
////////////////////////////////////////////////////////////////

/* class "tuning" */
tuning::tuning(){
  set= new Setting();
  set->Initialize();
}

tuning::~tuning(){}

void tuning::SetRoot(string ifname){
  tnew = new TChain("tree");
	cout << "SetRoot" <<endl;
//  fnew = new TFile(Form("%s",ifname.c_str()),"recreate");
//	cout << "SetRoot" <<endl;
//  tnew =new TTree("tnew",ifname.c_str());
//	cout << "SetRoot" <<endl;
//  tnew = tree->CloneTree(0);

	cout << "SetRoot" <<endl;
  add_tree(ifname);
	cout << "SetRoot" <<endl;
  pack_tree();
	cout << "SetRoot" <<endl;
  readtreeHRSR();
	cout << "SetRoot" <<endl;
  readtreeHRSL();

}
/////////////////////////////
void tuning::SetRunList(string ifname){

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    add_tree(runname);
    //    cout<<buf<<endl;
  }
  ENum=tree->GetEntries();
  cout<<"Events: "<<ENum<<endl; 

  pack_tree();
  readtreeHRSR();
  readtreeHRSL();
}

void tuning::ReadParam(string name){

  param = new ParamMan(name.c_str());
  cout<<"param name : "<<name<<endl;
  if(param -> SetVal())cout<<"F1TDC parameter setted"<<endl; 
  tdc_time=param->F1Res();
  coin_offset=param->GetF1CoinOffset();
  cout<<"coin off : "<<coin_offset<<endl;
  coin_shift=param->GetF1ShiftOffset();
  coin_shift = coin_shift*tdc_time;
  cout<<"coin shift : "<<coin_shift<<endl;
}
////////////////////////////////////////////////////////////////////////////
void tuning::CoinCalc(int RS2_seg, int LS2_seg, int rhit, int lhit){

  
  convertF1TDCR(param);
  convertF1TDCL(param);
  PathCalib(rhit,lhit);
  
  double Beta_R=R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+MK*MK);
  double Beta_L=L_tr_p[lhit]/sqrt(L_tr_p[lhit]*L_tr_p[lhit]+Me*Me);
  
 

  double tof_r  = RS2_F1time[RS2_seg] - R_pathl/(Beta_R*LightVelocity);
  double tof_l  = LS2_F1time[LS2_seg] - L_pathl/(Beta_L*LightVelocity);
  double tof_rc = RS2_F1time_c[RS2_seg] - R_pathl/(Beta_R*LightVelocity);
  double tof_lc = LS2_F1time_c[LS2_seg] - L_pathl/(Beta_L*LightVelocity);

//  double tof_lg = LS2_F1time_g[LS2_seg] - L_pathl/(Beta_L*LightVelocity);

  
    tr.RS2T_ref=RF1Ref[0];
    tr.RS2B_ref=RF1Ref[1];
    tr.LS2T_ref=LF1Ref[0];
    tr.LS2B_ref=LF1Ref[1];
    tr.RS2T_F1[RS2_seg]=RS2T_F1[RS2_seg];
    tr.RS2B_F1[RS2_seg]=RS2B_F1[RS2_seg];
    tr.LS2T_F1[LS2_seg]=LS2T_F1[LS2_seg];
    tr.LS2B_F1[LS2_seg]=LS2B_F1[LS2_seg];
    tr.RS2T_F1_c[RS2_seg]=RS2T_F1_c[RS2_seg];
    tr.RS2B_F1_c[RS2_seg]=RS2B_F1_c[RS2_seg];
    tr.LS2T_F1_c[LS2_seg]=LS2T_F1_c[LS2_seg];
    tr.LS2B_F1_c[LS2_seg]=LS2B_F1_c[LS2_seg];
    tr.RS2T_F1_b[RS2_seg]=RS2T_F1_b[RS2_seg];
    tr.RS2B_F1_b[RS2_seg]=RS2B_F1_b[RS2_seg];
    tr.LS2T_F1_b[LS2_seg]=LS2T_F1_b[LS2_seg];
    tr.LS2B_F1_b[LS2_seg]=LS2B_F1_b[LS2_seg];         
    tr.Rtof[RS2_seg]=tof_r;
    tr.Ltof[LS2_seg]=tof_l;
    

 
  if(RS2_F1time[RS2_seg]!=-9999. && LS2_F1time[LS2_seg]!=-9999.){
    ct       = - tof_rc + tof_lc - coin_offset;
    tr.ct_b  = - tof_r + tof_l - coin_offset;
    tr.ct_c  = - tof_rc + tof_lc - coin_offset;
    //    tr.ct_g    = - tof_rc + tof_lg - coin_offset;
  }else if(RS2_F1time[RS2_seg]!=-9999. && LS2_F1time[LS2_seg]==-9999.){
    tr.ct_b  = - tof_r + tof_l - coin_offset;
    tr.ct_c  = - tof_rc + tof_lc - coin_offset - coin_shift;
    ct       = - tof_rc + tof_lc - coin_offset - coin_shift;
    //    tr.ct_g    = - tof_rc + tof_lg - coin_offset - coin_shift;
  }else{
    ct=-1000;
    tr.ct_c =-1000;
    tr.ct_b =-1000;
    //    tr.ct_g =-1000;
  }


  if(RS2_seg<0 || LS2_seg<0){
    tr.ct_c =-1000;
    tr.ct_b =-1000;
    //    tr.ct_g =-1000;
  }


  tr.ct_g=-1000.;
  tr.ct_gb=-1000.;
  //tr.ct_g=CoinCalc_gogami(RS2_seg,LS2_seg,rhit,lhit);
  ct=CoinCalc_gogami(RS2_seg,LS2_seg,rhit,lhit);

//Changed
  
}


////////////////////////////////////////////////////////////////////////////
double tuning::CoinCalc_gogami(int RS2_seg, int LS2_seg,int rhit, int lhit){


//beta
  double beta_R  = R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+Mpi*Mpi);
  double beta_L  = L_tr_p[lhit]/sqrt(L_tr_p[lhit]*L_tr_p[lhit]+Me*Me);  

//length
  double LenL  = L_tr_pathl[lhit];// - L_s2_trpath[lhit];
  double LenR  = R_tr_pathl[rhit];// - R_s2_trpath[rhit];
  
//correction timing
  double cor_L = (LenL-3.18)/(beta_L*LightVelocity);
  double cor_R = (LenR- R_s2_trpath[rhit])/(beta_R*LightVelocity);



//  double tref_L  = LTDC_F1FirstHit[40]       * tdc_time;
  double tref_R  = RTDC_F1FirstHit[9]        * tdc_time;
  //  double rf      = LTDC_F1FirstHit[47]       * tdc_time;
  //  double rf_R    = RTDC_F1FirstHit[15]        * tdc_time;


 double timeL_R = RTDC_F1FirstHit[RS2_seg+16] * tdc_time;
 double timeR_R = RTDC_F1FirstHit[RS2_seg+48] * tdc_time; 
// double timeL_L = LTDC_F1FirstHit[LS2_seg] * tdc_time;
// double timeR_L = LTDC_F1FirstHit[LS2_seg+48] * tdc_time; 


 double toffset_R = -364.6-150.; // for H2_1
// double toffset_L = 1762.0;

// double meantime_L = tref_L - (timeL_L+timeR_L)/2.0 + toffset_L + cor_L;
 double meantime_R = tref_R - (timeL_R+timeR_R)/2.0 + toffset_R + cor_R;

 // meantime_R=100;RS2_seg=7;LS2_seg=7;


 double s2_tzero_R[16]={0,-1.02037,0.0046854,0.0534834,-0.534372,-0.60597,0.343139,-0.293262,0.267898,-0.666823,0.272364,0.0969059,-0.893806,-1.01129,-1.13495,-0.784991};
 double s2_tzero_L[16]={0,2.005,0.654413,1.34976,0.0290891,0.187557,0.00499421,-0.914343,-1.24058,-0.535878,-0.77564,2.22918,0.804909,0.607826,-0.635764,0};

 meantime_R= meantime_R - s2_tzero_R[RS2_seg] -s2_tzero_L[LS2_seg];

 double yfp_cor_R =R_tr_y[rhit]*-0.182869 +R_tr_ph[rhit]*-0.0211276;
 double yfp_cor_L = -10.0* L_tr_y[lhit] -28.0* L_tr_ph[lhit];

 // cout<<"yfp_cor_R "<<yfp_cor_R<<" yfp_corL "<<yfp_cor_L<<endl;
 
 meantime_R= meantime_R + yfp_cor_R + yfp_cor_L;
 meantime_R= meantime_R - cor_L +75.4;
 
 tr.yp_cor=0.0;
 tr.yp_cor= + yfp_cor_R + yfp_cor_L;




  //======= Nomalization ==================//
  R_tr_x[rhit]    = (R_tr_x[rhit]-XFPm)/XFPr;
  R_tr_th[rhit]   = (R_tr_th[rhit]-XpFPm)/XpFPr;
  R_tr_y[rhit]    = (R_tr_y[rhit]-YFPm)/YFPr;
  R_tr_ph[rhit]   = (R_tr_ph[rhit]-YpFPm)/YpFPr;
  R_tr_vz[rhit]   = (R_tr_vz[rhit]-Ztm)/Ztr;


  L_tr_x[lhit]    = (L_tr_x[lhit]-XFPm)/XFPr; 
  L_tr_th[lhit]   = (L_tr_th[lhit]-XpFPm)/XpFPr;
  L_tr_y[lhit]    = (L_tr_y[lhit]-YFPm)/YFPr;
  L_tr_ph[lhit]   = (L_tr_ph[lhit]-YpFPm)/YpFPr;
  L_tr_vz[lhit]   = (L_tr_vz[lhit]-Ztm)/Ztr;




 double ctimecorR = calcf2t_3rd(PctimeR, R_tr_x[rhit],R_tr_th[rhit],R_tr_y[rhit],R_tr_ph[rhit],R_tr_vz[rhit]);
 double ctimecorL = calcf2t_3rd(PctimeL, L_tr_x[lhit],L_tr_th[lhit],L_tr_y[lhit],L_tr_ph[lhit],L_tr_vz[lhit]);


    //========== Scaled at FP ==================//
    R_tr_x[rhit]  = R_tr_x[rhit]  * XFPr + XFPm;
    R_tr_th[rhit] = R_tr_th[rhit] * XpFPr + XpFPm;
    R_tr_y[rhit]  = R_tr_y[rhit]  * YFPr + YFPm;
    R_tr_ph[rhit] = R_tr_ph[rhit] * YpFPr + YpFPm;

    L_tr_x[lhit]  = L_tr_x[lhit]  * XFPr + XFPm;
    L_tr_th[lhit] = L_tr_th[lhit] * XpFPr + XpFPm;
    L_tr_y[lhit]  = L_tr_y[lhit]  * YFPr + YFPm;
    L_tr_ph[lhit] = L_tr_ph[lhit] * YpFPr + YpFPm;    



 double kcenter = 3.122;
 kcenter=0.0;
 double ctime = - meantime_R - kcenter + ctimecorL + ctimecorR;

 tr.ctimecorR=ctimecorR;
 tr.ctimecorL=ctimecorL;
 // double time_rf = rf - meantime;
 // double time_rf_R = rf - meantime_R -tref_R;
 // double ctime = - meantime_R + mean_time - kcenter;





//  cout<<"======== check coin time ====== "<<endl;
//  cout<<"s2R off "<<s2_tzero_R[RS2_seg]<<" s2L off "<<s2_tzero_L[LS2_seg]<<" meantime_R "<<meantime_R<<" ctime "<<ctime<<endl;
//  cout<<" meantime_L "<<meantime_L<<" meantime_R "<<meantime_R<<" corL "<<cor_L<<" corR "<<cor_R<<" ctime "<<ctime<<endl;
//  cout<<" ctimecorL "<<ctimecorL<<" ctimecorR "<<ctimecorR<<" yfp_cor_R "<<yfp_cor_R<<" yfp_cor_L "<<yfp_cor_L<<endl;

 ctime=ctime -1.4;


 tr.ct_gb=-ctime;

 if(3600.0<ctime && ctime<3665){
   ctime = ctime - 3637.88 - 12.76;
   ctime = ctime - 12.0-3.1;
 }

 ctime=-ctime;

 return ctime;


}

///////////////////////////////////////////////////////////////////////////
//double tuning::CoinCalc_c(int RS2_seg, int LS2_seg, int rhit, int lhit){
//
//  double cointime=0.0;
//
//  
//  convertF1TDCR(param);
//  convertF1TDCL(param);
//
//  //  double Rpathl=R_tr_pathl[rhit]+R_s2_trpath[rhit];
//  //  double Lpathl=L_tr_pathl[lhit]+L_s2_trpath[lhit];
//
//  double Beta_R=R_p/sqrt(R_p*R_p+MK*MK);
//  double Beta_L=L_p/sqrt(L_p*L_p+Me*Me);
//  double tof_r=RS2_F1time[RS2_seg] - R_pathl/(Beta_R*LightVelocity);
//  double tof_l=LS2_F1time[LS2_seg] - L_pathl/(Beta_L*LightVelocity);
//  double tof_rc=RS2_F1time_c[RS2_seg] - R_pathl/(Beta_R*LightVelocity);
//  double tof_lc=LS2_F1time_c[LS2_seg] - L_pathl/(Beta_L*LightVelocity);
//
//    
//  if(RS2_F1time[RS2_seg]!=-9999. &&LS2_F1time[LS2_seg]!=-9999.){
//    cointime= - tof_r + tof_l - coin_offset;
//    //    tr.ct_c= - (rtof[RS2_seg] - R_pathl/(Beta_R*LightVelocity))
//    //      + (ltof[LS2_seg] - L_pathl/(Beta_L*LightVelocity))        - coin_offset;
//    tr.ct_c= - tof_rc + tof_lc -coin_offset;
//  }
//  else{
//    cointime=-1000;
//    tr.ct_c =-1000;
//  }
//
//  
//  if(tof_r!=tof_rc)
//    cout<<" ct "<<cointime<<" ct_c "<<tr.ct_c<<endl;
//  
//  return cointime;
//  
//}

///////////////////////////////////////////////////////////////////////////

void tuning::PathCalib(int rhit, int lhit){


  R_pathl=0.0;
  L_pathl=0.0;
    
  R_pathl= R_tr_pathl[rhit] ;// + R_s2_trpath[rhit];
  L_pathl= L_tr_pathl[lhit] ;// + L_s2_trpath[lhit];



  
  //  tr.Rpathl=R_pathl;
  //  tr.Lpathl=L_pathl;

  R_pathl       = (R_pathl - PaRm )/PaRr;
  L_pathl       = (L_pathl - PaLm )/PaLr;  
  
  R_tr_x[rhit]  = (R_tr_x[rhit]-XFPm)/XFPr;
  R_tr_th[rhit] = (R_tr_th[rhit]-XpFPm)/XpFPr;
  R_tr_y[rhit]  = (R_tr_y[rhit]-YFPm)/YFPr;
  R_tr_ph[rhit] = (R_tr_ph[rhit]-YpFPm)/YpFPr;
  R_tr_vz[rhit] = (R_tr_vz[rhit] - Ztm)/Ztr;
  
  L_tr_x[lhit]  = (L_tr_x[lhit]-XFPm)/XFPr;
  L_tr_th[lhit] = (L_tr_th[lhit]-XpFPm)/XpFPr;
  L_tr_y[lhit]  = (L_tr_y[lhit]-YFPm)/YFPr;
  L_tr_ph[lhit] = (L_tr_ph[lhit]-YpFPm)/YpFPr;
  L_tr_vz[lhit] = (L_tr_vz[lhit] - Ztm)/Ztr;      
  

  
  //==== Calc Path Length =========//

  if(MT_p[10])  R_pathl = calcf_pathl(Ppl  ,R_tr_x[rhit],R_tr_th[rhit],R_tr_y[rhit],R_tr_ph[rhit],R_tr_vz[rhit]); // ns
  if(MT_p[11])  L_pathl = calcf_pathl(Ppl_L,L_tr_x[lhit],L_tr_th[lhit],L_tr_y[lhit],L_tr_ph[lhit],L_tr_vz[lhit]); // ns

  
  R_tr_x[rhit]  = R_tr_x[rhit] * XFPr + XFPm;
  R_tr_th[rhit] = R_tr_th[rhit] * XpFPr + XpFPm;
  R_tr_y[rhit]  = R_tr_y[rhit] * YFPr + YFPm;
  R_tr_ph[rhit] = R_tr_ph[rhit] * YpFPr   + YpFPm;
  R_tr_vz[rhit] = R_tr_vz[rhit] * Ztr + Ztm;
  
  L_tr_x[lhit]  = L_tr_x[lhit] * XFPr + XFPm;
  L_tr_th[lhit] = L_tr_th[lhit] * XpFPr + XpFPm;
  L_tr_y[lhit]  = L_tr_y[lhit] * YFPr + YFPm;
  L_tr_ph[lhit] = L_tr_ph[lhit] * YpFPr   + YpFPm;
  L_tr_vz[lhit] = L_tr_vz[lhit] * Ztr + Ztm;
  
  R_pathl = R_pathl * PaRr + PaRm + R_s2_trpath[rhit] ;
  L_pathl = L_pathl * PaLr + PaLm + L_s2_trpath[lhit] ;

  
  tr.Rpathl = R_pathl;
  tr.Lpathl = L_pathl;



}


//////////////////////////////////////////////////////////////////

double tuning::Eloss(double yp,double z,char const* arm){

  double hrs_ang=13.2*3.14159/180.;  
  double x;
  

  //----- Original coordinate  -------//
  // Definition by K.N. Suzuki  (fixed Oct. 23rd, 2019)//
  // R-HRS : right hand coordinate (Unticlockwise rotation)//
  // L-HRS : left  hand coordinate (    Clockwise rotation)//
  
  if(*arm=='R')        x = - hrs_ang - yp; //yp : phi [rad] RHRS
  else if(*arm=='L')   x = - hrs_ang + yp; //yp : phi [rad] LHRS
  else x=0.0;
  double ph[3],pl[2];
  double dEloss=0.0;
  bool high;
  double dEloss_h = 0.0;
  double dEloss_l = 0.0;
  if(z>0.08)high=false;
  else high=true;
  
    //==== thickness 0.400 mm ========//

  if(*arm=='R'){
    ph[0] = -1.3175;
    ph[1] = -4.6151;
    ph[2] = 2.0369;
    pl[0] = 3.158e-2;
    pl[1] = 4.058e-1;
  }else if(*arm=='L'){
    ph[0] = -1.3576;
    ph[1] = -4.5957;
    ph[2] = 2.0909;
    pl[0] = 6.2341e-3;
    pl[1] = 4.0336e-1;
  }else ph[0]=0.;ph[1]=0.;ph[2]=0.;pl[0]=0.;pl[1]=0.;//Okuyama

 
  if(high){
    dEloss_h = ph[0]*sin(ph[1]*x)+ph[2];    
    dEloss = dEloss_h;
  }else{
    dEloss_l = pl[0]*x +pl[1];    
    dEloss = dEloss_l;}
  //==== thickness 0.4 mm in beam energy loss ======//
  if(*arm=='B')dEloss=0.184; //[MeV/c]
  dEloss=dEloss/1000.; // [GeV/c]
  return dEloss;

  
}

#if 0
void tuning::SetBranch(){

	T->SetBranchStatus("*",0);

//------ Right Arm -------------//

 T->SetBranchStatus("RTDC.F1FirstHit",1);
 T->SetBranchAddress("RTDC.F1FirstHit",RF1); 
 T->SetBranchStatus("R.s2.t_pads",1);
 T->SetBranchAddress("R.s2.t_pads",Rs2tpads);
 T->SetBranchStatus("R.s2.trpad",1);
 T->SetBranchAddress("R.s2.trpad",Rs2trpad);
 T->SetBranchStatus("R.a1.a_c",1);
 T->SetBranchAddress("R.a1.a_c",Ra1a_c);
 T->SetBranchStatus("R.a2.a_c",1);
 T->SetBranchAddress("R.a2.a_c",Ra2a_c); 
// T->SetBranchStatus("R.a1.asum_c",1);
// T->SetBranchAddress("R.a1.asum_c",&R_a1_asum_p);
// T->SetBranchStatus("R.a2.asum_c",1);
// T->SetBranchAddress("R.a2.asum_c",&R_a2_asum_p);
 T->SetBranchStatus("R.a1.asum_p",1);
 T->SetBranchAddress("R.a1.asum_p",&R_a1_asum_p);
 T->SetBranchStatus("R.a2.asum_p",1);
 T->SetBranchAddress("R.a2.asum_p",&R_a2_asum_p);
//change 
 T->SetBranchStatus("R.vdc.u1.time",1);
 T->SetBranchAddress("R.vdc.u1.time",Ru1_time);
 T->SetBranchStatus("Ndata.R.vdc.u1.time",1);
 T->SetBranchAddress("Ndata.R.vdc.u1.time",&NRu1_time);
 // path length//
 T->SetBranchStatus("R.s2.trpath",1); 
 T->SetBranchAddress("R.s2.trpath",rs2pathl); 
 T->SetBranchStatus("R.tr.pathl",1);  
 T->SetBranchAddress("R.tr.pathl",rtrpathl);
 // Target positon information //
 T->SetBranchStatus("R.tr.p",1);
 T->SetBranchAddress("R.tr.p",Rp);
 T->SetBranchStatus("R.tr.vz",1);    
 T->SetBranchAddress("R.tr.vz",Rvz); 

 //------ Left Arm a---------------//
 T->SetBranchStatus("LTDC.F1FirstHit",1);
 T->SetBranchAddress("LTDC.F1FirstHit",LF1); 
 T->SetBranchStatus("L.s2.t_pads",1);
 T->SetBranchAddress("L.s2.t_pads",Ls2tpads);
 T->SetBranchStatus("L.s2.trpad",1);
 T->SetBranchAddress("L.s2.trpad",Ls2trpad);
  // path length//
 T->SetBranchStatus("L.s2.trpath",1); 
 T->SetBranchAddress("L.s2.trpath",ls2pathl); 
 T->SetBranchStatus("L.tr.pathl",1);   
 T->SetBranchAddress("L.tr.pathl",ltrpathl);
 T->SetBranchStatus("L.tr.p",1);
 T->SetBranchAddress("L.tr.p",Lp);  
 T->SetBranchStatus("L.tr.vz",1);    
 T->SetBranchAddress("L.tr.vz",Lvz);
}
#endif

void tuning::SetParam(){

	mt = Mp;//target mass
	mh = ML;//hypernuclei
cout << "mt = " << mt << endl;

	min_Lp = 1.8;
	max_Lp = 2.8;

   //min_coin=-10;
   min_coin=-20.0;
   max_coin=20.0;
   //min_coin_c=-10;
   min_coin_c=-20.0;
   max_coin_c=20.0;
//change
//   min_coin_c=-100;
//   max_coin_c=1000.0;
 min_ac1=0.0;
 max_ac1=5000.;
 min_ac2=0.0;
 max_ac2=20000.;
 min_adc=-500.0;
 max_adc=20000.;

	//=== AC Threshold variable ===//
	 //th1_max=2000.;
	 //ac1_adc[0]=400.;
	 //ac1_adc[1]=420.;
	 //ac1_adc[2]=440.;
	 //ac1_adc[3]=460.;
	 //ac1_adc[4]=480.;
	 //ac1_adc[5]=500.;
	 //ac1_adc[6]=520.;
	 //ac1_adc[7]=540.;
	 //ac1_adc[8]=560.;
	 //ac1_adc[9]=580.;

     //th2_max=6000.;
     //ac2l_adc[0]=600.;
     //ac2l_adc[1]=650.;
     //ac2l_adc[2]=700.;
     //ac2l_adc[3]=750.;
     //ac2l_adc[4]=800.;
     //ac2l_adc[5]=850.;
     //ac2l_adc[6]=900.;
     //ac2l_adc[7]=950.;
     //ac2l_adc[8]=1000.;
     //ac2l_adc[9]=1050.;

     //ac2u_adc[0]=1500.;
     //ac2u_adc[1]=1600.;
     //ac2u_adc[2]=1700.;
     //ac2u_adc[3]=1800.;
     //ac2u_adc[4]=1900.;
     //ac2u_adc[5]=2000.;
     //ac2u_adc[6]=2500.;
     //ac2u_adc[7]=3000.;
     //ac2u_adc[8]=3500.;
     //ac2u_adc[9]=4000.;
     //----------------NPE---------------//
	 th1_max=2.0;
	 zver[0]=0.;
	 for(int i=1; i<nth; i++) zver[i]=zver[i-1]+5.;
	 zver_diff[0]=0.;
	 for(int i=1; i<nth; i++) zver_diff[i]=zver_diff[i-1]+1.;

	 ac1_adc[0]=0.0;
	 for(int i=1; i<nth; i++) ac1_adc[i]=ac1_adc[i-1]+0.25;
	 //for(int i=1; i<nth; i++) ac1_adc[i]=ac1_adc[i-1]+2.;
	 //ac1_adc[1]=0.4;
	 //ac1_adc[2]=0.6;
	 //ac1_adc[3]=0.8;
	 //ac1_adc[4]=1.0;
	 //ac1_adc[5]=1.2;
	 //ac1_adc[6]=1.4;
	 //ac1_adc[7]=1.6;
	 //ac1_adc[8]=1.8;
	 //ac1_adc[9]=2.0;

     th2_max=20.0;
     th2_min=0.0;
     ac2l_adc[0]=0.0;//best = 3.0? ac2l[30]
	 for(int i=1; i<nth; i++) ac2l_adc[i]=ac2l_adc[i-1]+0.2;
	 //for(int i=1; i<nth; i++) ac2l_adc[i]=ac2l_adc[i-1]+1;
//     ac2l_adc[0]=0.4;
     //ac2l_adc[1]=0.8;
     //ac2l_adc[2]=1.2;
     //ac2l_adc[3]=1.6;
     //ac2l_adc[4]=2.0;
     //ac2l_adc[5]=2.4;
     //ac2l_adc[6]=2.8;
     //ac2l_adc[7]=3.2;
     //ac2l_adc[8]=3.6;
     //ac2l_adc[9]=4.0;

     ac2u_adc[0]=0.0;//best = 18? ac2u[16]
	 for(int i=1; i<nth; i++) ac2u_adc[i]=ac2u_adc[i-1]+0.5;
	 //for(int i=1; i<nth; i++) ac2u_adc[i]=ac2u_adc[i-1]+3;
//     ac2u_adc[0]=14.0;
     //ac2u_adc[1]=15.0;
     //ac2u_adc[2]=16.0;
     //ac2u_adc[3]=17.0;
     //ac2u_adc[4]=18.0;
     //ac2u_adc[5]=19.0;
     //ac2u_adc[6]=20.0;
     //ac2u_adc[7]=21.0;
     //ac2u_adc[8]=22.0;
     //ac2u_adc[9]=23.0;
	// th1_max=2.0;

	// ac1_adc[0]=0.;
	// ac1_adc[1]=0.5;
	// ac1_adc[2]=1.0;
	// ac1_adc[3]=1.5;
	// ac1_adc[4]=2.0;
	// ac1_adc[5]=2.5;
	// ac1_adc[6]=3.0;
	// ac1_adc[7]=3.5;
	// ac1_adc[8]=4.0;
	// ac1_adc[9]=4.5;

    // th2_max=6.0;
    // th2_min=2.0;
    // ac2l_adc[0]=0.;
    // ac2l_adc[1]=0.4;
    // ac2l_adc[2]=0.8;
    // ac2l_adc[3]=1.2;
    // ac2l_adc[4]=1.6;
    // ac2l_adc[5]=2.0;
    // ac2l_adc[6]=2.4;
    // ac2l_adc[7]=2.8;
    // ac2l_adc[8]=3.2;
    // ac2l_adc[9]=3.6;

    // ac2u_adc[0]=4.0;
    // ac2u_adc[1]=5.0;
    // ac2u_adc[2]=6.0;
    // ac2u_adc[3]=7.0;
    // ac2u_adc[4]=8.0;
    // ac2u_adc[5]=9.0;
    // ac2u_adc[6]=10.0;
    // ac2u_adc[7]=11.0;
    // ac2u_adc[8]=12.0;
    // ac2u_adc[9]=13.0;
     //-----------tuning_ac_oku_10.pdf--------------//
	 //ac1_adc[0]=0.;
	 //ac1_adc[1]=0.1;
	 //ac1_adc[2]=0.2;
	 //ac1_adc[3]=0.3;
	 //ac1_adc[4]=0.4;
	 //ac1_adc[5]=0.5;
	 //ac1_adc[6]=0.6;
	 //ac1_adc[7]=0.7;
	 //ac1_adc[8]=0.8;
	 //ac1_adc[9]=0.9;

     //th2_max=6.0;
     //th2_min=2.0;
     //ac2l_adc[0]=0.;
     //ac2l_adc[1]=0.4;
     //ac2l_adc[2]=0.8;
     //ac2l_adc[3]=1.2;
     //ac2l_adc[4]=0.6;
     //ac2l_adc[5]=2.0;
     //ac2l_adc[6]=2.4;
     //ac2l_adc[7]=2.8;
     //ac2l_adc[8]=3.2;
     //ac2l_adc[9]=3.6;

     //ac2u_adc[0]=6.0;
     //ac2u_adc[1]=6.2;
     //ac2u_adc[2]=6.4;
     //ac2u_adc[3]=6.6;
     //ac2u_adc[4]=6.8;
     //ac2u_adc[5]=7.0;
     //ac2u_adc[6]=7.2;
     //ac2u_adc[7]=7.4;
     //ac2u_adc[8]=7.6;
     //ac2u_adc[9]=7.8;
	//---Kaon Cut ----//
 	ac1_kcut=100.;
 	ac2_kcut_min=1000.;
 	ac2_kcut_max=5000.;


}
////////////////////////////////////////////////////////////

void tuning::MakeHist(){
  cout<<"Make Hist "<<endl;
  tree_out = new TTree("T","T");
  //`tree_out ->Branch("branch name",variable ,"branch name/type");

  tree_out ->Branch("pid_cut"        ,&tr.pid_cut      ,"pid_cut/I"     );
  tree_out ->Branch("ct_cut"        ,&tr.ct_cut      ,"ct_cut/I"     );
  tree_out ->Branch("z_cut"        ,&tr.z_cut      ,"z_cut/I"     );
  tree_out ->Branch("nrun"        ,&tr.nrun      ,"nrun/I"     );
  tree_out ->Branch("nev"        ,&tr.nev      ,"nev/I"     );
  tree_out ->Branch("ntr_r",&tr.ntrack_r ,"ntr_r/I");
  tree_out ->Branch("ntr_l",&tr.ntrack_l ,"ntr_l/I");
  tree_out ->Branch("mm",&tr.missing_mass ,"missing_mass/D");
  tree_out ->Branch("mm_b",&tr.missing_mass_b ,"missing_mass_b/D");
  tree_out ->Branch("mm_L",&tr.missing_mass_L ,"missing_mass_L/D");
  tree_out ->Branch("mm_nnL",&tr.missing_mass_nnL ,"missing_mass_nnL/D");
  tree_out ->Branch("mm_H3L",&tr.missing_mass_H3L ,"missing_mass_H3L/D");
  tree_out ->Branch("mm_cut",&tr.missing_mass_cut ,"missing_mass_cut/D");
  tree_out ->Branch("mm_MgL",&tr.missing_mass_MgL ,"missing_mass_MgL/D");
  tree_out ->Branch("mm_MgL_acc",&tr.missing_mass_MgL_acc ,"missing_mass_MgL_acc/D");
  tree_out ->Branch("mm_acc",&tr.missing_mass_acc ,"missing_mass_acc/D");
  tree_out ->Branch("runnum",&runnum ,"runnum/I");
  tree_out ->Branch("ct_b",&tr.ct_b ,"ct_b/D");
  tree_out ->Branch("ct_c",&tr.ct_c ,"ct_c/D");
  tree_out ->Branch("ct_g",&tr.ct_g ,"ct_g/D");
  tree_out ->Branch("yp_cor",&tr.yp_cor ,"yp_cor/D");
  tree_out ->Branch("ctimecorR",&tr.ctimecorR ,"ctimecorR/D");
  tree_out ->Branch("ctimecorL",&tr.ctimecorL ,"ctimecorL/D");
  tree_out ->Branch("ct_gb",&tr.ct_gb ,"ct_gb/D");
  tree_out ->Branch("rtof"        ,tr.Rtof      ,"rtof[16]/D"     );
  tree_out ->Branch("ltof"        ,tr.Ltof      ,"ltof[16]/D"     );  
  tree_out ->Branch("RS2T"        ,tr.RS2T_F1      ,"RS2T_F1[16]/D"     );
  tree_out ->Branch("RS2B"        ,tr.RS2B_F1      ,"RS2B_F1[16]/D"     );
  tree_out ->Branch("LS2T"        ,tr.LS2T_F1      ,"LS2T_F1[16]/D"     );
  tree_out ->Branch("LS2B"        ,tr.LS2B_F1      ,"LS2B_F1[16]/D"     );
  tree_out ->Branch("RS2T_c"        ,tr.RS2T_F1_c      ,"RS2T_F1_c[16]/D"     );
  tree_out ->Branch("RS2B_c"        ,tr.RS2B_F1_c      ,"RS2B_F1_c[16]/D"     );
  tree_out ->Branch("LS2T_c"        ,tr.LS2T_F1_c      ,"LS2T_F1_c[16]/D"     );
  tree_out ->Branch("LS2B_c"        ,tr.LS2B_F1_c      ,"LS2B_F1_c[16]/D"     );
  tree_out ->Branch("RS2T_b"        ,tr.RS2T_F1_b      ,"RS2T_F1_b[16]/D"     );
  tree_out ->Branch("RS2B_b"        ,tr.RS2B_F1_b      ,"RS2B_F1_b[16]/D"     );
  tree_out ->Branch("LS2T_b"        ,tr.LS2T_F1_b      ,"LS2T_F1_b[16]/D"     );
  tree_out ->Branch("LS2B_b"        ,tr.LS2B_F1_b      ,"LS2B_F1_b[16]/D"     );      
  tree_out ->Branch("RS2T_ref"        ,&tr.RS2T_ref   ,"RS2T_ref/D"     );
  tree_out ->Branch("RS2B_ref"        ,&tr.RS2B_ref   ,"RS2B_ref/D"     );
  tree_out ->Branch("LS2T_ref"        ,&tr.LS2T_ref   ,"LS2T_ref/D"     );
  tree_out ->Branch("LS2B_ref"        ,&tr.LS2B_ref   ,"LS2B_ref/D"     );  
  tree_out ->Branch("ct"   ,&tr.coin_time ,"coin_time/D");
  tree_out ->Branch("Rp"        ,&tr.momR      ,"momR/D"     );
  tree_out ->Branch("Lp"        ,&tr.momL      ,"momL/D"     );
  tree_out ->Branch("Rs2_pad",tr.Rs2_pad,"Rs2_pad[100]/I");
  tree_out ->Branch("Ls2_pad",tr.Ls2_pad,"Ls2_pad[100]/I");

  tree_out ->Branch("Rth_fp"          ,&tr.RXpFP        ,"RXpFP/D"       );  
  tree_out ->Branch("Lth_fp"          ,&tr.LXpFP        ,"LXpFP/D"       );  
  tree_out ->Branch("Rph_fp"          ,&tr.RYpFP        ,"RYpFP/D"       );  
  tree_out ->Branch("Lph_fp"          ,&tr.LYpFP        ,"LYpFP/D"       );
  tree_out ->Branch("Rx_fp"          ,&tr.RXFP        ,"RXFP/D"       );
  tree_out ->Branch("Lx_fp"          ,&tr.LXFP        ,"LXFP/D"       );  
  tree_out ->Branch("Ry_fp"          ,&tr.RYFP        ,"RYFP/D"       );  
  tree_out ->Branch("Ly_fp"          ,&tr.LYFP        ,"LYFP/D"       );
    
  tree_out ->Branch("Rth"          ,&tr.RXpt        ,"RXpt/D"       );  
  tree_out ->Branch("Lth"          ,&tr.LXpt        ,"LXpt/D"       );  
  tree_out ->Branch("Rph"          ,&tr.RYpt        ,"RYpt/D"       );  
  tree_out ->Branch("Lph"          ,&tr.LYpt        ,"LYpt/D"       );
  tree_out ->Branch("Rx"          ,&tr.RXt        ,"RXt/D"       );
  tree_out ->Branch("Lx"          ,&tr.LXt        ,"LXt/D"       );  
  tree_out ->Branch("Ry"          ,&tr.RYt        ,"RYt/D"       );
  tree_out ->Branch("Ly"          ,&tr.LYt        ,"LYt/D"       );    
  tree_out ->Branch("Rz"          ,&tr.zR        ,"zR/D"       );
  tree_out ->Branch("Lz"          ,&tr.zL        ,"zL/D"       );

  tree_out ->Branch("ac1_sum"     ,&tr.AC1_sum   ,"AC1_sum/D"  );
  tree_out ->Branch("ac2_sum"     ,&tr.AC2_sum   ,"AC2_sum/D"  );
  tree_out ->Branch("ac1_npe_sum"     ,&tr.AC1_npe_sum   ,"AC1_npe_sum/D"  );
  tree_out ->Branch("ac2_npe_sum"     ,&tr.AC2_npe_sum   ,"AC2_npe_sum/D"  );
  tree_out ->Branch("ac1_npe"     ,tr.AC1_npe   ,"AC1_npe[24]/D"  );
  tree_out ->Branch("ac2_npe"     ,tr.AC2_npe   ,"AC2_npe[26]/D"  );    
  tree_out ->Branch("ct_acc"     ,&tr.ct_acc   ,"ct_acc/D"  );
  tree_out ->Branch("Rs0ra_p"     ,&tr.Rs0ra_p   ,"Rs0ra_p/D"  );
  tree_out ->Branch("Rs0la_p"     ,&tr.Rs0la_p   ,"Rs0la_p/D"  );
  tree_out ->Branch("Rs2ra_p"     ,tr.Rs2ra_p   ,"Rs2ra_p[16]/D"  );
  tree_out ->Branch("Rs2la_p"     ,tr.Rs2la_p   ,"Rs2la_p[16]/D"  );
  tree_out ->Branch("Ls2ra_p"     ,tr.Ls2ra_p   ,"Ls2ra_p[16]/D"  );
  tree_out ->Branch("Ls2la_p"     ,tr.Ls2la_p   ,"Ls2la_p[16]/D"  );  
  tree_out->Branch("Bp"     ,&tr.Bp   ,"Bp/D"  );
  //  tree_out->Branch("Lp"     ,tr.Lp   ,"Lp[100]/D"  );
  //  tree_out->Branch("Rp"     ,tr.Rp   ,"Rp[100]/D"  );
  tree_out->Branch("Bp_c"     ,&tr.Bp_c   ,"Bp_c/D"  );
  tree_out->Branch("Lp_c"     ,tr.Lp_c   ,"Lp_c[100]/D"  );
  tree_out->Branch("Rp_c"     ,tr.Rp_c   ,"Rp_c[100]/D"  );
  tree_out ->Branch("trig"     ,&tr.trig   ,"trig/D"  );
  tree_out->Branch("dpe"     ,&tr.dpe   ,"dpe/D"  );
  tree_out->Branch("dpe_"     ,tr.dpe_   ,"dpe_[10]/D"  );
  tree_out->Branch("dpk"     ,tr.dpk   ,"dpk[10]/D"  );

//////////////
//Parameters//
//////////////

  min_mm=-0.1;//GeV/c^2
  max_mm=0.2;//GeV/c^2
  bin_mm=(max_mm-min_mm)/0.002; //Counts/2 MeV
  bin_mm=(int)bin_mm;
	iter_ac1=30;//iteration number

/////////////
//// AC  ////
/////////////

	  npe_sum_a1 = new TH1F("npe_sum_a1","NPE SUM A1",2000,0.,40.);
//	  npe_sum_a2 = new TH1F("npe_sum_a2","NPE SUM A2",4000,0.,80.);
	  npe_sum_a2 = new TH1F("npe_sum_a2","NPE SUM A2 (Kaon)",300,0.,30.);
////--------------------------------------------//
	  h_pisr1 = new TH1F("h_pisr1","Pion Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_ksr1 = new TH1F("h_ksr1","Kaon Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_psr1 = new TH1F("h_psr1","Proton Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_Lsr1 = new TH1F("h_Lsr1","Lambda Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_Ssr1 = new TH1F("h_Ssr1","Signa Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_pitot1 = new TH1F("h_pitot1","Pion Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_ktot1 = new TH1F("h_ktot1","Kaon Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_ptot1 = new TH1F("h_ptot1","Proton Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_Ltot1 = new TH1F("h_Ltot1","Lambda Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_Stot1 = new TH1F("h_Stot1","Sigma Survival Ratio",nth-1,zver[0],zver[nth-1]);
//--------------------------------------------//
	  h_pisr2l = new TH1F("h_pisr2l","Pion Survival Ratio",nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_ksr2l = new TH1F("h_ksr2l","Kaon Survival Ratio",nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_psr2l = new TH1F("h_psr2l","Proton Survival Ratio",nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_Lsr2l = new TH1F("h_Lsr2l","Lambda Survival Ratio",nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_Ssr2l = new TH1F("h_Ssr2l","Sigma Survival Ratio",nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_pitot2l = new TH1F("h_pitot2l","Pion Survival Ratio",nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_ktot2l = new TH1F("h_ktot2l","Kaon Survival Ratio",nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_ptot2l = new TH1F("h_ptot2l","Proton Survival Ratio",nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_Ltot2l = new TH1F("h_Ltot2l","Lambda Survival Ratio",nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_Stot2l = new TH1F("h_Stot2l","Sigma Survival Ratio",nth-1,zver_diff[0],zver_diff[nth-1]);
//--------------------------------------------//
	  h_ksr11 = new TH1F("h_ksr11","Kaon Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_ksr111 = new TH1F("h_ksr111","Kaon Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_ksr1111 = new TH1F("h_ksr1111","Kaon Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_ktot11 = new TH1F("h_ktot11","Kaon Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_ktot111 = new TH1F("h_ktot111","Kaon Survival Ratio",nth-1,zver[0],zver[nth-1]);
	  h_ktot1111 = new TH1F("h_ktot1111","Kaon Survival Ratio",nth-1,zver[0],zver[nth-1]);
//2D Eff.
	  h_pisr2d = new TH2F("h_pisr2d","Pion Survival Ratio",nth-1,zver[0],zver[nth-1],nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_ksr2d = new TH2F("h_ksr2d","Kaon Survival Ratio",nth-1,zver[0],zver[nth-1],nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_psr2d = new TH2F("h_psr2d","Proton Survival Ratio",nth-1,zver[0],zver[nth-1],nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_Lsr2d = new TH2F("h_Lsr2d","Lambda Survival Ratio",nth-1,zver[0],zver[nth-1],nth-1,zver_diff[0],zver_diff[nth-1]);
	  h_Ssr2d = new TH2F("h_Ssr2d","Sigma Survival Ratio",nth-1,zver[0],zver[nth-1],nth-1,zver_diff[0],zver_diff[nth-1]);
//--------------------------------------------//
	  h_pisr1->Sumw2();
	  h_pisr1->SetLineColor(kRed);
	  h_ksr1->Sumw2();
	  h_ksr1->SetLineColor(kRed);
	  h_psr1->Sumw2();
	  h_psr1->SetLineColor(kRed);
	  h_Lsr1->Sumw2();
	  h_Lsr1->SetLineColor(kAzure);
	  h_Ssr1->Sumw2();
	  h_Ssr1->SetLineColor(kCyan);
	  h_pitot1->Sumw2();
	  h_pitot1->SetLineColor(kBlack);
	  h_ktot1->Sumw2();
	  h_ktot1->SetLineColor(kBlack);
	  h_ptot1->Sumw2();
	  h_ptot1->SetLineColor(kBlack);
	  h_Ltot1->Sumw2();
	  h_Ltot1->SetLineColor(kAzure);
	  h_Stot1->Sumw2();
	  h_Stot1->SetLineColor(kCyan);
	  h_pisr2l->Sumw2();
	  h_pisr2l->SetLineColor(kRed);
	  h_ksr2l->Sumw2();
	  h_ksr2l->SetLineColor(kRed);
	  h_psr2l->Sumw2();
	  h_psr2l->SetLineColor(kRed);
	  h_Lsr2l->Sumw2();
	  h_Lsr2l->SetLineColor(kAzure);
	  h_Ssr2l->Sumw2();
	  h_Ssr2l->SetLineColor(kCyan);
	  h_pitot2l->Sumw2();
	  h_pitot2l->SetLineColor(kBlack);
	  h_ktot2l->Sumw2();
	  h_ktot2l->SetLineColor(kBlack);
	  h_ptot2l->Sumw2();
	  h_ptot2l->SetLineColor(kBlack);
	  h_Ltot2l->Sumw2();
	  h_Ltot2l->SetLineColor(kAzure);
	  h_Stot2l->Sumw2();
	  h_Stot2l->SetLineColor(kCyan);
  
	  h_ksr11->Sumw2();
	  h_ksr11->SetLineColor(kRed);
	  h_ksr111->Sumw2();
	  h_ksr111->SetLineColor(kRed);
	  h_ksr1111->Sumw2();
	  h_ksr1111->SetLineColor(kRed);
	  h_ktot11->Sumw2();
	  h_ktot11->SetLineColor(kBlack);
	  h_ktot111->Sumw2();
	  h_ktot111->SetLineColor(kBlack);
	  h_ktot1111->Sumw2();
	  h_ktot1111->SetLineColor(kBlack);
/////////////
//// BPM ////
/////////////
  h_rbay_rbax = new TH2D("h_rbay_rbax","h_rbay_rbax",200,-4,4,200,-3,7);
  h_rbby_rbbx = new TH2D("h_rbby_rbbx","h_rbby_rbbx",200,-4,4,200,-3,7);
  h_rby_rbx   = new TH2D("h_rby_rbx"  ,"h_rby_rbx"  ,200,-6,4,200,-6,4);
  set->SetTH2(h_rbay_rbax,"BPM A"         ,"X","Y");
  set->SetTH2(h_rbby_rbbx,"BPM B"         ,"X","Y");
  set->SetTH2(h_rby_rbx  ,"Raster Pattern","X","Y");

//////////////
//// LHRS ////
//////////////
  h_L_trig = new TH1D("h_L_trig","h_L_trig",10,0,10);
  set->SetTH1(h_L_trig,"Trigger Flag","Trig No.","Counts");

  h_L_tr_n      = new TH1D("h_L_tr_n"     ,"h_L_tr_n"     ,15 ,    0,  15);
  h_L_tr_ch2    = new TH1D("h_L_tr_ch2"   ,"h_L_tr_ch2"   ,400,    0,0.03);
  h_L_p         = new TH1D("h_L_p"        ,"h_L_p"        ,400,  1.9, 2.3);
  h_L_pathl     = new TH1D("h_L_pathl"    ,"h_L_pathl"    ,400, 25.2,26.3);
  h_L_px        = new TH1D("h_L_px"       ,"h_L_px"       ,400, 0.35, 0.6);
  h_L_py        = new TH1D("h_L_py"       ,"h_L_py"       ,400, -0.2, 0.2);
  h_L_pz        = new TH1D("h_L_pz"       ,"h_L_pz"       ,400, 1.85,2.25);
  h_L_tgy       = new TH1D("h_L_tgy"      ,"h_L_tgy"      ,400,-0.06,0.06);
  h_L_tgth      = new TH1D("h_L_tgth"     ,"h_L_tgth"     ,400, -0.1, 0.1);
  h_L_tgph      = new TH1D("h_L_tgph"     ,"h_L_tgph"     ,400,-0.06,0.06);
  h_L_vx        = new TH1D("h_L_vx"       ,"h_L_vx"       ,400,-0.005,0.002);
  h_L_vy        = new TH1D("h_L_vy"       ,"h_L_vy"       ,400,-0.004,0.003);
  h_L_vz        = new TH1D("h_L_vz"       ,"h_L_vz"       ,400,-0.25,0.25);
  h_L_y_x       = new TH2D("h_L_y_x"      ,"h_L_y_x"      ,200,   -1,  1 ,200,-0.1,0.1);
  h_L_th_x      = new TH2D("h_L_th_x"     ,"h_L_th_x"     ,200,   -1,  1 ,200,-0.2,0.2);
  h_L_ph_y      = new TH2D("h_L_ph_y"     ,"h_L_ph_y"     ,200, -0.1, 0.1,200,-0.1,0.1);
  h_L_tgph_tgth = new TH2D("h_L_tgph_tgth","h_L_tgph_tgth",200, -0.1, 0.1,200,-0.06,0.06);
  set->SetTH1(h_L_tr_n     ,"No. of Tracks"           ,"No. of Tracks"   ,"Counts");
  set->SetTH1(h_L_tr_ch2   ,"Tracking #chi^{2}"       ,"#chi^{2}"        ,"Counts");
  set->SetTH1(h_L_p        ,"Track Momentum"          ,"p (GeV/#it{c})"  ,"Counts");
  set->SetTH1(h_L_pathl    ,"Track Path Length"       ,"l (m)"           ,"Counts");
  set->SetTH1(h_L_px       ,"Momentum X"              ,"px (GeV/#it{c})" ,"Counts");
  set->SetTH1(h_L_py       ,"Momentum Y"              ,"py (GeV/#it{c})" ,"Counts");
  set->SetTH1(h_L_pz       ,"Momentum Z"              ,"pz (GeV/#it{c})" ,"Counts");
  set->SetTH1(h_L_tgy      ,"Target Plane Y"          ,"y_{t} (m)"       ,"Counts");
  set->SetTH1(h_L_tgth     ,"Target Plane #theta"     ,"#theta_{t} (rad)","Counts");
  set->SetTH1(h_L_tgph     ,"Target Plane #phi"       ,"#phi_{t} (rad)"  ,"Counts");
  set->SetTH1(h_L_vx       ,"Vertex X"                ,"x_{v} (m)"       ,"Counts");
  set->SetTH1(h_L_vy       ,"Vertex Y"                ,"y_{v} (m)"       ,"Counts");
  set->SetTH1(h_L_vz       ,"Vertex Z"                ,"z_{v} (m)"       ,"Counts");
  set->SetTH2(h_L_y_x      ,"Focal Plane Y v.s X"     ,"X (m)"           ,"Y (m)");
  set->SetTH2(h_L_th_x     ,"Focal Plane #theta v.s X","X (m)"           ,"#theta (rad)");
  set->SetTH2(h_L_ph_y     ,"Focal Plane #phi v.s Y"  ,"Y (m)"           ,"#phi (rad)");
  set->SetTH2(h_L_tgph_tgth,"Target #phi v.s #theta"  ,"#theta_{t} (rad)","#phi_{t} (rad)");

  h_L_beta        = new TH1D("h_L_beta"       ,"h_L_beta"       ,400,   0,  2); 
  h_L_m2          = new TH1D("h_L_m2"         ,"h_L_m2"         ,400,-0.5,2.5); 
  h_L_beta_p      = new TH2D("h_L_beta_p"     ,"h_L_beta_p"     ,200, 1.9,2.3,200,   0,  2); 
  h_L_beta_m2     = new TH2D("h_L_beta_m2"    ,"h_L_beta_m2"    ,200,-0.5,  2,200,   0,  2); 
  h_L_dedx_p      = new TH2D("h_L_dedx_p"     ,"h_L_dedx_p"     ,200, 1.9,2.3,200,   0, 10); 
  h_L_dedx_m2     = new TH2D("h_L_dedx_m2"    ,"h_L_dedx_m2"    ,200,-0.5,  2,200,   0, 10); 
  h_L_s0_dedx     = new TH1D("h_L_s0_dedx"    ,"h_L_s0_dedx"    ,400,   0, 10); 
  h_L_s0_beta_x   = new TH2D("h_L_s0_beta_x"  ,"h_L_s0_beta_x"  ,200,  -1,  1,200,   0,  2); 
  h_L_s0_dedx_x   = new TH2D("h_L_s0_dedx_x"  ,"h_L_s0_dedx_x"  ,200,  -1,  1,200,   0, 10); 
  h_L_s2_pad      = new TH1D("h_L_s2_pad"     ,"h_L_s2_pad"     , 18,  -1, 17); 
  h_L_s2_dedx     = new TH1D("h_L_s2_dedx"    ,"h_L_s2_dedx"    ,400,   0, 10); 
  h_L_s2_beta_x   = new TH2D("h_L_s2_beta_x"  ,"h_L_s2_beta_x"  ,200,  -1,  1,200,   0,  2); 
  h_L_s2_dedx_x   = new TH2D("h_L_s2_dedx_x"  ,"h_L_s2_dedx_x"  ,200,  -1,  1,200,   0, 10); 
  h_L_s2_beta_pad = new TH2D("h_L_s2_beta_pad","h_L_s2_beta_pad", 16,  0, 16,200,    0,  2); 
  h_L_s2_dedx_pad = new TH2D("h_L_s2_dedx_pad","h_L_s2_dedx_pad", 16,  0, 16,200,    0, 10); 
  set->SetTH1(h_L_beta       ,"Track beta"                    ,"#beta"                ,"Counts");
  set->SetTH1(h_L_m2         ,"Mass Square"                   ,"M^{2} (GeV^{2}/c^{4})","Counts");
  set->SetTH2(h_L_beta_p     ,"#beta v.s Momentum"            ,"p (GeV/c)"            ,"#beta",0.0);
  set->SetTH2(h_L_beta_m2    ,"#beta v.s Mass Square"         ,"M^{2} (GeV^{2}/c^{4})","#beta");
  set->SetTH2(h_L_dedx_p     ,"Energy Deposit v.s Momentum"   ,"p (GeV/c)"            ,"dE/dx ()");
  set->SetTH2(h_L_dedx_m2    ,"Energy Deposit v.s Mass Square","M^{2} (GeV^{2}/c^{4})","dE/dx ()");
  set->SetTH1(h_L_s0_dedx    ,"Energy Deposit (S0)"           ,"dE/dx"                ,"Counts");
  set->SetTH2(h_L_s0_beta_x  ,"#beta v.s X-pos (S0)"          ,"X (m)"                ,"#beta");
  set->SetTH2(h_L_s0_dedx_x  ,"Energy Deposit (S0) v.s X-pos" ,"X (m)"                ,"dE/dx ()");
  set->SetTH1(h_L_s2_pad     ,"Hit Paddle (S2)"               ,"Paddle No."           ,"Counts");
  set->SetTH1(h_L_s2_dedx    ,"Energy Deposit (S2)"           ,"dE/dx"                ,"Counts");
  set->SetTH2(h_L_s2_beta_x  ,"#beta v.s X-pos (S2)"          ,"X (m)"                ,"#beta");
  set->SetTH2(h_L_s2_dedx_x  ,"Energy Deposit (S2) v.s X-pos" ,"X (m)"                ,"dE/dx ()");
  set->SetTH2(h_L_s2_beta_pad,"#beta v.s Paddle (S2)"         ,"Paddle No."           ,"#beta");
  set->SetTH2(h_L_s2_dedx_pad,"Energy Deposit (S2) v.s Paddle","Paddle No."           ,"dE/dx ()");

  h_L_tgt       = new TH1D("h_L_tgt"      ,"h_L_tgt"       ,40000,-2000,2000);
  h_L_s2pad_tgt = new TH2D("h_L_s2pad_tgt","h_L_s2pad_tgt" ,200,-1000,1000, 16,    0,  16);
  h_L_p_tgt     = new TH2D("h_L_p_tgt"    ,"h_L_p_tgt"     ,200,-1000,1000,200,  1.9, 2.3);
  h_L_pathl_tgt = new TH2D("h_L_pathl_tgt","h_L_pathl_tgt" ,200,-1000,1000,200, 25.2,26.3);
  h_L_tgy_tgt   = new TH2D("h_L_tgy_tgt"  ,"h_L_tgy_tgt"   ,200,-1000,1000,200,-0.06,0.06);
  h_L_tgth_tgt  = new TH2D("h_L_tgth_tgt" ,"h_L_tgth_tgt"  ,200,-1000,1000,200, -0.1, 0.1);
  h_L_tgph_tgt  = new TH2D("h_L_tgph_tgt" ,"h_L_tgph_tgt"  ,200,-1000,1000,200,-0.06,0.06);
  h_L_x_tgt     = new TH2D("h_L_x_tgt"    ,"h_L_x_tgt"     ,200,-1000,1000,200,   -1,   1);
  h_L_y_tgt     = new TH2D("h_L_y_tgt"    ,"h_L_y_tgt"     ,200,-1000,1000,200, -0.1, 0.1);
  set->SetTH1(h_L_tgt      ,"Time at Target (S2-RF)","Time (ns)"         ,"Counts");
  set->SetTH2(h_L_s2pad_tgt,"S2 Paddle v.s Time at Target (S2-RF)"       ,"Time (ns)","Paddle No.");
  set->SetTH2(h_L_p_tgt    ,"Momentum v.s Time at Target (S2-RF)"        ,"Time (ns)","Momentum (GeV/c)");
  set->SetTH2(h_L_pathl_tgt,"Path Length v.s Time at Target (S2-RF)"     ,"Time (ns)","L (m)");
  set->SetTH2(h_L_tgy_tgt  ,"Y at Target v.s Time at Target (S2-RF)"     ,"Time (ns)","Y_{t}");
  set->SetTH2(h_L_tgth_tgt ,"#theta at Target v.s Time at Target (S2-RF)","Time (ns)","#theta_{t}");
  set->SetTH2(h_L_tgph_tgt ,"#phi at Target v.s Time at Target (S2-RF)"  ,"Time (ns)","#phi_{t}");
  set->SetTH2(h_L_x_tgt    ,"X at FP v.s Time at Target (S2-RF)"         ,"Time (ns)","X");
  set->SetTH2(h_L_y_tgt    ,"Y at FP v.s Time at Target (S2-RF)"         ,"Time (ns)","Y");



//////////////
//// RHRS ////
//////////////
  h_R_trig = new TH1D("h_R_trig","h_R_trig",10,0,10);
  set->SetTH1(h_R_trig,"Trigger Flag","Trig No.","Counts");

  h_R_tr_n      = new TH1D("h_R_tr_n"     ,"h_R_tr_n"     ,15 ,    0,  15);
  h_R_tr_ch2    = new TH1D("h_R_tr_ch2"   ,"h_R_tr_ch2"   ,400,    0,0.03);
  h_R_p         = new TH1D("h_R_p"        ,"h_R_p"        ,400,  1.7,1.95);
  h_R_pathl     = new TH1D("h_R_pathl"    ,"h_R_pathl"    ,400, 25.2,26.3);
  h_R_px        = new TH1D("h_R_px"       ,"h_R_px"       ,400, -0.5,-0.3);
  h_R_py        = new TH1D("h_R_py"       ,"h_R_py"       ,400, -0.2, 0.2);
  h_R_pz        = new TH1D("h_R_pz"       ,"h_R_pz"       ,400,  1.6,1.95);
  h_R_tgy       = new TH1D("h_R_tgy"      ,"h_R_tgy"      ,400,-0.06,0.06);
  h_R_tgth      = new TH1D("h_R_tgth"     ,"h_R_tgth"     ,400, -0.1, 0.1);
  h_R_tgph      = new TH1D("h_R_tgph"     ,"h_R_tgph"     ,400,-0.06,0.06);
  h_R_vx        = new TH1D("h_R_vx"       ,"h_R_vx"       ,400,-0.005,0.002);
  h_R_vy        = new TH1D("h_R_vy"       ,"h_R_vy"       ,400,-0.004,0.003);
  h_R_vz        = new TH1D("h_R_vz"       ,"h_R_vz"       ,400,-0.25,0.25);
  h_R_y_x       = new TH2D("h_R_y_x"      ,"h_R_y_x"      ,200,   -1,  1 ,200,-0.1,0.1);
  h_R_th_x      = new TH2D("h_R_th_x"     ,"h_R_th_x"     ,200,   -1,  1 ,200,-0.2,0.2);
  h_R_ph_y      = new TH2D("h_R_ph_y"     ,"h_R_ph_y"     ,200, -0.1, 0.1,200,-0.1,0.1);
  h_R_tgph_tgth = new TH2D("h_R_tgph_tgth","h_R_tgph_tgth",200, -0.1, 0.1,200,-0.06,0.06);
  set->SetTH1(h_R_tr_n  ,"No. of Tracks"           ,"No. of Tracks"   ,"Counts");
  set->SetTH1(h_R_tr_ch2,"Tracking #chi^{2}"       ,"#chi^{2}"        ,"Counts");
  set->SetTH1(h_R_p     ,"Track Momentum"          ,"p (GeV/#it{c})"  ,"Counts");
  set->SetTH1(h_R_pathl ,"Track Path Length"       ,"l (m)"           ,"Counts");
  set->SetTH1(h_R_px    ,"Momentum X"              ,"px (GeV/#it{c})" ,"Counts");
  set->SetTH1(h_R_py    ,"Momentum Y"              ,"py (GeV/#it{c})" ,"Counts");
  set->SetTH1(h_R_pz    ,"Momentum Z"              ,"pz (GeV/#it{c})" ,"Counts");
  set->SetTH1(h_R_tgy   ,"Target Plane Y"          ,"y_{t} (m)"       ,"Counts");
  set->SetTH1(h_R_tgth  ,"Target Plane #theta"     ,"#theta_{t} (rad)","Counts");
  set->SetTH1(h_R_tgph  ,"Target Plane #phi"       ,"#phi_{t} (rad)"  ,"Counts");
  set->SetTH1(h_R_vx    ,"Vertex X"                ,"x_{v} (m)"       ,"Counts");
  set->SetTH1(h_R_vy    ,"Vertex Y"                ,"y_{v} (m)"       ,"Counts");
  set->SetTH1(h_R_vz    ,"Vertex Z"                ,"z_{v} (m)"       ,"Counts");
  set->SetTH2(h_R_y_x   ,"Focal Plane Y v.s X"     ,"X (m)"           ,"Y (m)");
  set->SetTH2(h_R_th_x  ,"Focal Plane #theta v.s X","X (m)"           ,"#theta (rad)");
  set->SetTH2(h_R_ph_y  ,"Focal Plane #phi v.s Y"  ,"Y (m)"           ,"#phi (rad)");
  set->SetTH2(h_R_tgph_tgth,"Target #phi v.s #theta"  ,"#theta_{t} (rad)","#phi_{t} (rad)");

  h_R_beta      = new TH1D("h_R_beta"     ,"h_R_beta"     ,400,   0,  2); 
  h_R_m2        = new TH1D("h_R_m2"       ,"h_R_m2"       ,400,-0.5,2.5); 
  h_R_beta_p    = new TH2D("h_R_beta_p"   ,"h_R_beta_p"   ,200, 1.7,1.95,200,   0,   2); 
  h_R_beta_m2   = new TH2D("h_R_beta_m2"  ,"h_R_beta_m2"  ,200,-0.5,   2,200,   0,   2); 
  h_R_dedx_p    = new TH2D("h_R_dedx_p"   ,"h_R_dedx_p"   ,200, 1.7,1.95,200,   0,  10); 
  h_R_dedx_m2   = new TH2D("h_R_dedx_m2"  ,"h_R_dedx_m2"  ,200,-0.5,   2,200,   0,  10); 
  h_R_s0_dedx   = new TH1D("h_R_s0_dedx"  ,"h_R_s0_dedx"  ,400,   0,  10); 
  h_R_s0_beta_x = new TH2D("h_R_s0_beta_x","h_R_s0_beta_x",200,  -1,   1,200,   0,  10); 
  h_R_s0_dedx_x = new TH2D("h_R_s0_dedx_x","h_R_s0_dedx_x",200,  -1,   1,200,   0,   2); 
  h_R_s2_pad      = new TH1D("h_R_s2_pad"     ,"h_R_s2_pad"     , 18,  -1, 18); 
  h_R_s2_dedx     = new TH1D("h_R_s2_dedx"    ,"h_R_s2_dedx"    ,400,   0, 10); 
  h_R_s2_beta_x   = new TH2D("h_R_s2_beta_x"  ,"h_R_s2_beta_x"  ,200,  -1,  1,200,   0,  2); 
  h_R_s2_dedx_x   = new TH2D("h_R_s2_dedx_x"  ,"h_R_s2_dedx_x"  ,200,  -1,  1,200,   0, 10); 
  h_R_s2_beta_pad = new TH2D("h_R_s2_beta_pad","h_R_s2_beta_pad", 16,  0, 16,200,   0,  2); 
  h_R_s2_dedx_pad = new TH2D("h_R_s2_dedx_pad","h_R_s2_dedx_pad", 16,  0, 16,200,   0, 10); 
  h_R_a1_sum    = new TH1D("h_R_a1_sum"   ,"h_R_a1_sum"   ,400,   0,6000);
  h_R_a1_sum_x  = new TH2D("h_R_a1_sum_x" ,"h_R_a1_sum_x" ,200,  -1,   1,200,   0,6000); 
  h_R_a1_sum_p  = new TH2D("h_R_a1_sum_p" ,"h_R_a1_sum_p" ,200, 1.7,1.95,200,   0,6000); 
  h_R_a1_sum_m2 = new TH2D("h_R_a1_sum_m2","h_R_a1_sum_m2",200,-0.5, 2.5,200,   0,6000); 
  h_R_a2_sum    = new TH1D("h_R_a2_sum"   ,"h_R_a2_sum"   ,400,   0,30000);
  h_R_a2_sum_x  = new TH2D("h_R_a2_sum_x" ,"h_R_a2_sum_x" ,200,  -1,   1,200,   0,30000); 
  h_R_a2_sum_p  = new TH2D("h_R_a2_sum_p" ,"h_R_a2_sum_p" ,200, 1.7,1.95,200,   0,30000); 
  h_R_a2_sum_m2 = new TH2D("h_R_a2_sum_m2","h_R_a2_sum_m2",200,-0.5, 2.5,200,   0,30000); 
  set->SetTH1(h_R_beta       ,"Track beta"                        ,"#beta"               ,"Counts");
  set->SetTH1(h_R_m2         ,"Mass Square"                       ,"M^{2} (GeV^{2}/c^{4}","Counts");
  set->SetTH2(h_R_beta_p     ,"#beta v.s Momentum"                ,"p (GeV/c)"           ,"#beta");
  set->SetTH2(h_R_beta_m2    ,"#beta v.s Mass Square"             ,"M^{2} (GeV^{2}/c^{4}","#beta");
  set->SetTH2(h_R_dedx_p     ,"Energy Deposit v.s Momentum"       ,"p (GeV/c)"           ,"dE/dx ()");
  set->SetTH2(h_R_dedx_m2    ,"Energy Deposit v.s Mass Square"    ,"M^{2} (GeV^{2}/c^{4}","dE/dx ()");
  set->SetTH1(h_R_s0_dedx    ,"Energy Deposit (S0)"               ,"dE/dx"               ,"Counts");
  set->SetTH2(h_R_s0_beta_x  ,"#beta v.s X-pos (S0)"              ,"X (m)"               ,"#beta");
  set->SetTH2(h_R_s0_dedx_x  ,"Energy Deposit (S0) v.s X-pos"     ,"X (m)"               ,"dE/dx ()");
  set->SetTH1(h_R_s2_pad     ,"Hit Paddle (S2)"                   ,"Paddle No."          ,"Counts");
  set->SetTH1(h_R_s2_dedx    ,"Energy Deposit (S2)"               ,"dE/dx"               ,"Counts");
  set->SetTH2(h_R_s2_beta_x  ,"#beta v.s X-pos (S2)"              ,"X (m)"               ,"#beta");
  set->SetTH2(h_R_s2_dedx_x  ,"Energy Deposit (S2) v.s X-pos"     ,"X (m)"               ,"dE/dx ()");
  set->SetTH2(h_R_s2_beta_pad,"#beta v.s Paddle (S2)"             ,"Paddle No."          ,"#beta");
  set->SetTH2(h_R_s2_dedx_pad,"Energy Deposit (S2) v.s Paddle"    ,"Paddle No."          ,"dE/dx ()");
  set->SetTH1(h_R_a1_sum     ,"Cherenkov SUM (A1)"                ,""                    ,"Counts");
  set->SetTH2(h_R_a1_sum_x   ,"Cherenkov SUM v.s X-pos (A1)"      ,"X (m)"               ,"");
  set->SetTH2(h_R_a1_sum_p   ,"Cherenkov SUM v.s Momentum (A1)"   ,"p (GeV/c)"           ,"");
  set->SetTH2(h_R_a1_sum_m2  ,"Cherenkov SUM v.s Mass Square (A1)","M^{2} (GeV^{2}/c^{4}","");
  set->SetTH1(h_R_a2_sum     ,"Cherenkov SUM (A2)"                ,""                    ,"Counts");
  set->SetTH2(h_R_a2_sum_x   ,"Cherenkov SUM v.s X-pos (A2)"      ,"X (m)"               ,"");
  set->SetTH2(h_R_a2_sum_p   ,"Cherenkov SUM v.s Momentum (A2)"   ,"p (GeV/c)"           ,"");
  set->SetTH2(h_R_a2_sum_m2  ,"Cherenkov SUM v.s Mass Square (A2)","M^{2} (GeV^{2}/c^{4}","");

  h_R_tgt       = new TH1D("h_R_tgt"      ,"h_R_tgt"       ,40000,-2000,2000);
  h_R_s2pad_tgt = new TH2D("h_R_s2pad_tgt","h_R_s2pad_tgt" ,200,-1000,1000, 16,    0,  16);
  h_R_p_tgt     = new TH2D("h_R_p_tgt"    ,"h_R_p_tgt"     ,200,-1000,1000,200,  1.7,1.95);
  h_R_pathl_tgt = new TH2D("h_R_pathl_tgt","h_R_pathl_tgt" ,200,-1000,1000,200, 25.2,26.3);
  h_R_tgy_tgt   = new TH2D("h_R_tgy_tgt"  ,"h_R_tgy_tgt"   ,200,-1000,1000,200,-0.06,0.06);
  h_R_tgth_tgt  = new TH2D("h_R_tgth_tgt" ,"h_R_tgth_tgt"  ,200,-1000,1000,200,-0.1,0.1);
  h_R_tgph_tgt  = new TH2D("h_R_tgph_tgt" ,"h_R_tgph_tgt"  ,200,-1000,1000,200,-0.06,0.06);
  h_R_x_tgt     = new TH2D("h_R_x_tgt"    ,"h_R_x_tgt"     ,200,-1000,1000,200,  -1,  1);
  h_R_y_tgt     = new TH2D("h_R_y_tgt"    ,"h_R_y_tgt"     ,200,-1000,1000,200,-0.1,0.1);
  set->SetTH1(h_R_tgt      ,"Time at Target (S2-RF)","Time (ns)"         ,"Counts");
  set->SetTH2(h_R_s2pad_tgt,"S2 Paddle v.s Time at Target (S2-RF)"       ,"Time (ns)","Paddle No.");
  set->SetTH2(h_R_p_tgt    ,"Momentum v.s Time at Target (S2-RF)"        ,"Time (ns)","Momentum (GeV/c)");
  set->SetTH2(h_R_pathl_tgt,"Path Length v.s Time at Target (S2-RF)"     ,"Time (ns)","L (m)");
  set->SetTH2(h_R_tgy_tgt  ,"Y at Target v.s Time at Target (S2-RF)"     ,"Time (ns)","Y_{t}");
  set->SetTH2(h_R_tgth_tgt ,"#theta at Target v.s Time at Target (S2-RF)","Time (ns)","#theta_{t}");
  set->SetTH2(h_R_tgph_tgt ,"#phi at Target v.s Time at Target (S2-RF)"  ,"Time (ns)","#phi_{t}");
  set->SetTH2(h_R_x_tgt    ,"X at FP v.s Time at Target (S2-RF)"         ,"Time (ns)","X");
  set->SetTH2(h_R_y_tgt    ,"Y at FP v.s Time at Target (S2-RF)"         ,"Time (ns)","Y");

/////////////////////
//// Coincidence ////
/////////////////////
  h_ct       = new TH1D("h_ct"      ,"h_ct"      ,1000, -20, 20);//to adjust offset
  h_Rs2      = new TH1D("h_Rs2"      ,"h_Rs2"      ,4000, -100, 100);
  h_Ls2      = new TH1D("h_Ls2"      ,"h_Ls2"      ,4000, -100, 100);
  h_ct_acc       = new TH1D("h_ct_acc"  ,"h_ct_acc"   ,4000, -80, 80);//to adjust offset 
  h_ct_wK    = new TH1D("h_ct_wK"   ,"h_ct_wK"   ,1000, -20, 20); 
  h_ct_wK_z  = new TH1D("h_ct_wK_z" ,"h_ct_wK_z" ,4000, -80, 80); 
  h_ct_wK_z_all  = new TH1D("h_ct_wK_z_all" ,"h_ct_wK_z_all",4000, -80, 80); 
  h_ct_wK_acc    = new TH1D("h_ct_wK_acc"   ,"h_ct_wK_acc"   ,4000, -80, 80); 
  h_ct_wK_z_acc  = new TH1D("h_ct_wK_z_acc" ,"h_ct_wK_z_acc" ,4000, -80, 80); 
  h_Rs2x_ct  = new TH2D("h_Rs2x_ct" ,"h_Rs2x_ct" , 200, -20, 20,200,   -1,  1); 
  h_ct_Rp  = new TH2D("h_ct_Rp" ,"h_ct_Rp" ,400,  1.7,1.95,4000, -80, 80);//to adjust offset 
  h_Ls2x_ct  = new TH2D("h_Ls2x_ct" ,"h_Ls2x_ct" , 200, -20, 20,200,   -1,  1); 
  h_a1sum_ct = new TH2D("h_a1sum_ct","h_a1sum_ct", 200, -20, 20,200,    0,6000); 
  h_a2sum_ct = new TH2D("h_a2sum_ct","h_a2sum_ct", 200, -20, 20,200,    0,30000); 

  //--------- Missing mass range ---=====--------------------//



  
  h_mm       = new TH1D("h_mm"      ,"h_mm"      , bin_mm,min_mm,max_mm);  //range bin=2 MeV
  h_mm_acc   = new TH1D("h_mm_acc"  ,"h_mm_acc"  , bin_mm,min_mm,max_mm);  //range bin=2 MeV
  h_peak_mm       = new TH1D("h_peak_mm"      ,"h_mm_peak"      , bin_mm,min_mm,max_mm); //bin=2 MeV
  h_mm_pi       = new TH1D("h_mm_pi"      ,"h_mm_pi"      , bin_mm,min_mm,max_mm); //Lambda Pion mass range bin=2 MeV
  h_mm_Al      = new TH1D("h_mm_Al","h_mm_Al",bin_mm,min_mm,max_mm); //Alminium mass bin=4 MeV
  h_mm_Al_acc      = new TH1D("h_mm_Al_acc","h_mm_Al_acc",bin_mm,min_mm,max_mm); //Alminium mass bin=4 MeV
  h_mm_Al_bg      = new TH1D("h_mm_Al_bg","h_mm_Al_bg",bin_mm,min_mm,max_mm); //Alminium mass bin=4 MeV
  h_peak_Al      = new TH1D("h_peak_Al","h_peak_Al",bin_mm,min_mm,max_mm); //Alminium mass bin=4 MeV

  h_mm_MgL      = new TH1D("h_mm_MgL","h_mm_MgL",bin_mm,min_mm,max_mm); //Alminium mass bin=4 MeV
  h_mm_MgL_acc      = new TH1D("h_mm_MgL_acc","h_mm_MgL_acc",bin_mm,min_mm,max_mm); //Mg mass bin=4 MeV
  h_peak_MgL      = new TH1D("h_peak_MgL","h_peak_MgL",bin_mm,min_mm,max_mm); //MgL mass bin=4 MeV


  h_mm_L       = new TH1D("h_mm_L"      ,"h_mm_L"      , bin_mm,min_mm,max_mm ); //Lambda mass range bin=2 MeV
  h_mm_L_ec       = new TH1D("h_mm_L_ec"      ,"h_mm_L_ec"      , bin_mm,min_mm,max_mm); //Lambda mass range bin=2 MeV  
  h_mm_nnL       = new TH1D("h_mm_nnL"      ,"h_mm_nnL"      , bin_mm,min_mm,max_mm); //nnL mass range bin=2 MeV
  h_mm_H3L       = new TH1D("h_mm_H3L"      ,"h_mm_H3L"      , bin_mm,min_mm,max_mm); //H3L mass range bin=2 MeV  
  h_acc_L       = new TH1D("h_acc_L"      ,"h_acc_L"      , bin_mm,min_mm,max_mm); //Lambda mass ACC  bin=2 MeV
  h_acc_nnL       = new TH1D("h_acc_nnL"      ,"h_acc_nnL"      , bin_mm,min_mm,max_mm); //nnL mass ACC bin=2 MeV
  h_acc_H3L       = new TH1D("h_acc_H3L"      ,"h_acc_H3L"      , bin_mm,min_mm,max_mm); //H3L mass ACC bin=2 MeV  
  h_peak_L       = new TH1D("h_peak_L"      ,"h_peak_L"      , bin_mm,min_mm,max_mm); //Lambda mass range bin=2 MeV
  h_peak_nnL       = new TH1D("h_peak_nnL"      ,"h_peak_nnL"      , bin_mm,min_mm,max_mm); //nnL mass range bin=2 MeV
  h_peak_H3L       = new TH1D("h_peak_H3L"      ,"h_peak_H3L"      , bin_mm,min_mm,max_mm); //H3L mass range bin=2 MeV  


  h_Rz     = new TH1D("h_Rz", "h_Rz",1000,-0.2,0.2);
  h_Rz_c   = new TH1D("h_Rz_c", "h_Rz_c",1000,-0.2,0.2);
  h_Rz_cut   = new TH1D("h_Rz_cut", "h_Rz_cut",1000,-0.2,0.2);
  h_Lz     = new TH1D("h_Lz", "h_Lz",1000,-0.2,0.2);
  h_Lz_c   = new TH1D("h_Lz_c", "h_Lz_c",1000,-0.2,0.2);

  h_Rth     = new TH1D("h_Rth", "h_Rth",1000,-0.1,0.1);
  h_Rth_c   = new TH1D("h_Rth_c", "h_Rth_c",1000,-0.1,0.1);
  h_Lth     = new TH1D("h_Lth", "h_Lth",1000,-0.1,0.1);
  h_Lth_c   = new TH1D("h_Lth_c", "h_Lth_c",1000,-0.1,0.1);

  h_Rph     = new TH1D("h_Rph", "h_Rph",1000,-0.1,0.1);
  h_Rph_c   = new TH1D("h_Rph_c", "h_Rph_c",1000,-0.1,0.1);
  h_Lph     = new TH1D("h_Lph", "h_Lph",1000,-0.1,0.1);
  h_Lph_c   = new TH1D("h_Lph_c", "h_Lph_c",1000,-0.1,0.1);

  h_Rp     = new TH1D("h_Rp", "h_Rp",1000,1.5,2.5);
  h_Rp_c   = new TH1D("h_Rp_c", "h_Rp_c",1000,1.5,2.5);
  h_Lp     = new TH1D("h_Lp", "h_Lp",1000,1.8,2.8);
  h_Lp_c   = new TH1D("h_Lp_c", "h_Lp_c",1000,1.8,2.8);
  

  h_Lp_mm    = new TH2D("h_Lp_mm"   ,"h_Lp_mm"   , bin_mm,min_mm,max_mm,bin_Lp,min_Lp,max_Lp); 
  //h_mmall    = new TH1D("h_mmall"   ,"h_mmall"   , 100,-1,1); 
  h_mmall    = new TH1D("h_mmall"   ,"h_mmall"   , bin_mm,min_mm,max_mm); 
  h_mmfoil   = new TH1D("h_mmfoil"  ,"h_mmfoil"  , bin_mm,min_mm,max_mm); 
  h_mmbg     = new TH1D("h_mmbg"    ,"h_mmbg"    , bin_mm,min_mm,max_mm); 
  h_mmallbg  = new TH1D("h_mmallbg" ,"h_mmallbg" , bin_mm,min_mm,max_mm); 
  h_mmfoilbg = new TH1D("h_mmfoilbg","h_mmfoilbg", bin_mm,min_mm,max_mm); 

  h_Ll_mm    = new TH2D("h_Ll_mm"   ,"h_Ll_mm"   , bin_mm,min_mm,max_mm,200,  25.2, 26.3); 
  h_Ltgy_mm  = new TH2D("h_Ltgy_mm" ,"h_Ltgy_mm" , bin_mm,min_mm,max_mm,200, -0.06, 0.06); 
  h_Ltgth_mm = new TH2D("h_Ltgth_mm","h_Ltgth_mm", bin_mm,min_mm,max_mm,200,  -0.1,  0.1); 
  h_Ltgph_mm = new TH2D("h_Ltgph_mm","h_Ltgph_mm", bin_mm,min_mm,max_mm,200, -0.06, 0.06); 
  h_Lvx_mm   = new TH2D("h_Lvx_mm"  ,"h_Lvx_mm"  , bin_mm,min_mm,max_mm,200,-0.005,0.002); 
  h_Lvy_mm   = new TH2D("h_Lvy_mm"  ,"h_Lvy_mm"  , bin_mm,min_mm,max_mm,200,-0.004,0.003); 
  h_Lvz_mm   = new TH2D("h_Lvz_mm"  ,"h_Lvz_mm"  , bin_mm,min_mm,max_mm,200, -0.25, 0.25); 
  h_Lx_mm    = new TH2D("h_Lx_mm"   ,"h_Lx_mm"   , bin_mm,min_mm,max_mm,200,    -1,   1); 
  h_Ly_mm    = new TH2D("h_Ly_mm"   ,"h_Ly_mm"   , bin_mm,min_mm,max_mm,200,  -0.1, 0.1); 
  h_Lth_mm   = new TH2D("h_Lth_mm"  ,"h_Lth_mm"  , bin_mm,min_mm,max_mm,200,  -0.2, 0.2); 
  h_Lph_mm   = new TH2D("h_Lph_mm"  ,"h_Lph_mm"  , bin_mm,min_mm,max_mm,200,  -0.1, 0.1); 
  h_Rp_mm    = new TH2D("h_Rp_mm"   ,"h_Rp_mm"   , bin_mm,min_mm,max_mm,200,   1.7, 1.95); 
  h_Rl_mm    = new TH2D("h_Rl_mm"   ,"h_Rl_mm"   , bin_mm,min_mm,max_mm,200,  25.2, 26.3); 
  h_Rtgy_mm  = new TH2D("h_Rtgy_mm" ,"h_Rtgy_mm" , bin_mm,min_mm,max_mm,200, -0.06, 0.06); 
  h_Rtgth_mm = new TH2D("h_Rtgth_mm","h_Rtgth_mm", bin_mm,min_mm,max_mm,200,  -0.1,  0.1); 
  h_Rtgph_mm = new TH2D("h_Rtgph_mm","h_Rtgph_mm", bin_mm,min_mm,max_mm,200, -0.06, 0.06); 
  h_Rvx_mm   = new TH2D("h_Rvx_mm"  ,"h_Rvx_mm"  , bin_mm,min_mm,max_mm,200,-0.005,0.002); 
  h_Rvy_mm   = new TH2D("h_Rvy_mm"  ,"h_Rvy_mm"  , bin_mm,min_mm,max_mm,200,-0.004,0.003); 
  h_Rvz_mm   = new TH2D("h_Rvz_mm"  ,"h_Rvz_mm"  , bin_mm,min_mm,max_mm,200, -0.25, 0.25); 
  h_Rx_mm    = new TH2D("h_Rx_mm"   ,"h_Rx_mm"   , bin_mm,min_mm,max_mm,200,    -1,    1); 
  h_Ry_mm    = new TH2D("h_Ry_mm"   ,"h_Ry_mm"   , bin_mm,min_mm,max_mm,200,  -0.1,  0.1); 
  h_Rth_mm   = new TH2D("h_Rth_mm"  ,"h_Rth_mm"  , bin_mm,min_mm,max_mm,200,  -0.2,  0.2); 
  h_Rph_mm   = new TH2D("h_Rph_mm"  ,"h_Rph_mm"  , bin_mm,min_mm,max_mm,200,  -0.1,  0.1); 
  h_Rp_Lp    = new TH2D("h_Rp_Lp"   ,"h_Rp_Lp"   , bin_Lp,min_Lp,max_Lp,200,   1.7, 1.95); 
  h_m2_mm    = new TH2D("h_m2_mm"   ,"h_m2_mm"   , 100,  -0.4, 1.4,bin_mm,min_mm,max_mm); 
 // h_m2_mm    = new TH2D("h_m2_mm"   ,"h_m2_mm"   ,bin_mm,min_mm,max_mm,400,0.,20.); 
  h_m2_ac    = new TH2D("h_m2_ac"   ,"h_m2_ac"   , 100,-0.4,1.4,400,  0., 20.); 


  set->SetTH1(h_ct      ,"Coincidence Time"                      ,"Cointime (ns)"           ,"Counts");
  h_ct->SetMinimum(0.8);
  set->SetTH1(h_ct_wK   ,"Coincidence Time (w/ K cut)"           ,"Cointime (ns)"           ,"Counts",1,3001,3);
  set->SetTH1(h_ct_wK_z ,"Coincidence Time (w/ K cut & Gas)"     ,"Cointime (ns)"           ,"Counts",1,3001,5);
  set->SetTH1(h_ct_wK_z_all ,"Coincidence Time (w/ K cut & Gas)"     ,"Cointime (ns)"           ,"Counts",1,3001,5);
  set->SetTH1(h_ct_wK_acc   ,"Coincidence Time ACC (w/ K cut)"           ,"Cointime (ns)"           ,"Counts",1,3001,3);
  set->SetTH1(h_ct_wK_z_acc ,"Coincidence Time ACC(w/ K cut & Gas)"     ,"Cointime (ns)"           ,"Counts",1,3001,5);
  set->SetTH2(h_Rs2x_ct ,"RHRS S2 X-pos v.s Cointime"            ,"Cointime (ns)"           ,"X (m)");
  set->SetTH2(h_Ls2x_ct ,"LHRS S2 X-pos v.s Cointime"            ,"Cointime (ns)"           ,"X (m)");
  set->SetTH2(h_a1sum_ct,"RHRS A1 SUM v.s Cointime"              ,"Cointime (ns)"           ,"");
  set->SetTH2(h_a2sum_ct,"RHRS A2 SUM v.s Cointime"              ,"Cointime (ns)"           ,"");
  set->SetTH1(h_mm      ,"Lambda Binding Energy w/o AC cut"          ,"Missing mass [GeV/c^2]","Counts/2 MeV");
  set->SetTH1(h_mm_acc   ,"Lambda Binding Energy ACC"          ,"Missing mass [GeV/c^2]","Counts/2 MeV");
  set->SetTH1(h_mm_Al_acc  ,"#Alminium Missing Mass(ACC)"          ,"Missing mass [GeV/c^2]","Counts/2 MeV");
  set->SetTH1(h_mm_Al_bg  ,"#Missing Mass( Al Back Ground)"          ,"Missing mass [GeV/c^2]","Counts/2 MeV");
    set->SetTH1(h_mm_MgL_acc  ,"#Mg27L Missing Mass(ACC)"          ,"Missing mass [GeV/c^2]","Counts/2 MeV");
  set->SetTH1(h_mm_pi     ,"Lambda(Pi mass) Binding Energy w/o AC cut"          ,"Missing mass [GeV/c^2]","Counts/2 MeV");
  set->SetTH1(h_mm_pi,"Lambda (Pi mass) Binding Energy w/o AC cut","Missing mass [GeV/c^2]","Counts/2 MeV");
  set->SetTH2(h_Ly_mm   ,"LHRS FP Y v.s B_{Lambda}"             ,"-B_{Lambda} (GeV/c^{2})","Y_{FP} (m)");
  set->SetTH2(h_Lth_mm  ,"LHRS FP theta v.s B_{Lambda}"        ,"-B_{Lambda} (GeV/c^{2})","theta_{FP} (rad)");
  set->SetTH2(h_Lph_mm  ,"LHRS FP phi v.s B_{Lambda}"          ,"-B_{Lambda} (GeV/c^{2})","phi_{FP} (rad)");
  set->SetTH2(h_Rp_Lp   ,"RHRS momentum v.s LHRS momentum"       ,"Lp (GeV/c)"              ,"Rp (GeV/c)");

//  TF1* fAl_R=new TF1("fAl_R","gausn(0)",-0.135,-0.115);


// min_vdc=-0.2e-6;
// max_vdc= 1.2e-6;
// bin_vdc=(max_vdc-min_vdc)/tdc_time*1.0e6;
// bin_vdc=(int)bin_vdc;
// min_s0=-10;
// max_s0=10000;
// bin_s0=int(max_s0-min_s0);
//
//
 min_s2=-10;
 max_s2=5000;
 bin_s2=max_s2-min_s2;
        bin_coin=(int)(max_coin-min_coin)/tdc_time;
        bin_coin_c=(int)((max_coin_c-min_coin_c)/tdc_time);
cout<<"tdc"<<tdc_time<<endl;
cout<<"max coin"<<max_coin_c<<endl;
cout<<"min coin"<<min_coin_c<<endl;
cout<<"bin coin"<<bin_coin_c<<endl;
//        bin_coin_c=(int)(max_coin_c-min_coin_c)/tdc_time;
//        bin_coin_c=(int)750;
//////////////////////////////////////////////////////////////should be changed
//        bin_beta=6000;
//	bin_adc=(int)max_adc-min_adc;
//	bin_ac1=(int)(max_ac1-min_ac1)*3; 
//	bin_ac2=(int)(max_ac2-min_ac2)*3; 

 hcoin_tc=new TH1F("hcoin_tc","Coincidence time w/ Path Length Correction  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
for (int i=0;i<nth;i++){
 hcoin_t2[i]=new TH1F(Form("hcoin_t2[%d]",i), Form("Coincidence %lf<AC2<%lf  cut",ac1_adc[i],th2_max),bin_coin_c,min_coin_c,max_coin_c);
   set->SetTH1(hcoin_t2[i],"Coincidence time ","","");
 hcoin_k_ac1[i]=new TH1F(Form("hcoin_k_ac1[%d]",i), Form("Cointime (Kaon) AC1<%lf  cut",ac1_adc[i]),bin_coin_c,min_coin_c,max_coin_c);
   set->SetTH1(hcoin_k_ac1[i],Form("Cointime (Kaon) AC1<%lf  cut",ac1_adc[i]),"","");
 hcoin_k_ac2[i]=new TH1F(Form("hcoin_k_ac2[%d]",i), Form("Cointime (Kaon) 1000<AC2<%lf  cut",ac2l_adc[i]),bin_coin_c,min_coin_c,max_coin_c);
   set->SetTH1(hcoin_k_ac2[i],Form("Cointime (Kaon) %lf<AC2<4000  cut",ac2l_adc[i]),"","");
 }


    //-------No AC1 cut------------//
	hcoin_k_fom_noAC1=new TH1F("hcoin_k_fom_noAC1", "No AC1 cut",bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_k_fom_noAC1,"No AC1 cut","","");
	hcoin_bg_fom_noAC1=new TH1F("hcoin_bg_fom_noAC1", "No AC1 cut", bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_bg_fom_noAC1,"No AC1 cut","","");
	hcoin_wo_bg_fom_noAC1=new TH1F("hcoin_wo_bg_fom_noAC1", "No AC1 cut", bin_coin_c,min_coin_c,max_coin_c);

	hmm_L_fom_noAC1=new TH1F("hmm_L_fom_noAC1", "No AC1 cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_L_fom_noAC1,"No AC1 cut","","");
	hmm_bg_fom_noAC1=new TH1F("hmm_bg_fom_noAC1", "No AC1 cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_bg_fom_noAC1,"No AC1 cut","","");
	hmm_wo_bg_fom_noAC1=new TH1F("hmm_wo_bg_fom_noAC1", "No AC1 cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_wo_bg_fom_noAC1,"No AC1 cut","","");
    //-------No AC2 cut------------//
	hcoin_k_fom_noAC2=new TH1F("hcoin_k_fom_noAC2", "No AC2 cut",bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_k_fom_noAC2,"No AC2 cut","","");
	hcoin_bg_fom_noAC2=new TH1F("hcoin_bg_fom_noAC2", "No AC2 cut", bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_bg_fom_noAC2,"No AC2 cut","","");
	hcoin_wo_bg_fom_noAC2=new TH1F("hcoin_wo_bg_fom_noAC2", "No AC2 cut", bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_wo_bg_fom_noAC2,"No AC2 cut","","");

	hmm_L_fom_noAC2=new TH1F("hmm_L_fom_noAC2", "No AC2 cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_L_fom_noAC2,"No AC2 cut","","");
	hmm_bg_fom_noAC2=new TH1F("hmm_bg_fom_noAC2", "No AC2 cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_bg_fom_noAC2,"No AC2 cut","","");
	hmm_wo_bg_fom_noAC2=new TH1F("hmm_wo_bg_fom_noAC2", "No AC2 cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_wo_bg_fom_noAC2,"No AC2 cut","","");
    //-------No AC cut------------//
	hcoin_k_fom_noAC=new TH1F("hcoin_k_fom_noAC", "No AC cut",bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_k_fom_noAC,"No AC cut","","");
	hcoin_bg_fom_noAC=new TH1F("hcoin_bg_fom_noAC", "No AC cut", bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_bg_fom_noAC,"No AC cut","","");
	hcoin_wo_bg_fom_noAC=new TH1F("hcoin_wo_bg_fom_noAC", "No AC cut", bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_wo_bg_fom_noAC,"No AC cut","","");

	hmm_L_fom_noAC=new TH1F("hmm_L_fom_noAC", "No AC cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_L_fom_noAC,"No AC cut","","");
	hmm_bg_fom_noAC=new TH1F("hmm_bg_fom_noAC", "No AC cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_bg_fom_noAC,"No AC cut","","");
	hmm_wo_bg_fom_noAC=new TH1F("hmm_wo_bg_fom_noAC", "No AC cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_wo_bg_fom_noAC,"No AC cut","","");
	hmm_pi_fom_noAC=new TH1F("hmm_pi_fom_noAC", "No AC cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_pi_fom_noAC,"No AC cut (Pion Selected)","","");
	hmm_pibg_fom_noAC=new TH1F("hmm_pibg_fom_noAC", "No AC cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_pibg_fom_noAC,"No AC cut","","");
	hmm_pi_wobg_fom_noAC=new TH1F("hmm_pi_wobg_fom_noAC", "No AC cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_pi_wobg_fom_noAC,"No AC cut","","");
    //-------No Z cut------------//
for(int i=0;i<nth;i++){
	hcoin_k_fom_noZ[i]=new TH1F(Form("hcoin_k_fom_noZ[%d]",i), "No Z cut",bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_k_fom_noZ[i],"No Z cut","","");
	hcoin_bg_fom_noZ[i]=new TH1F(Form("hcoin_bg_fom_noZ[%d]",i), "No Z cut", bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_bg_fom_noZ[i],"No Z cut","","");
	hcoin_wo_bg_fom_noZ[i]=new TH1F(Form("hcoin_wo_bg_fom_noZ[%d]",i), "No Z cut", bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_wo_bg_fom_noZ[i],"No Z cut","","");

	hmm_L_fom_noZ[i]=new TH1F(Form("hmm_L_fom_noZ[%d]",i), "No Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_L_fom_noZ[i],"No Z cut","","");
	hmm_bg_fom_noZ[i]=new TH1F(Form("hmm_bg_fom_noZ[%d]",i), "No Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_bg_fom_noZ[i],"No Z cut","","");
	hmm_wo_bg_fom_noZ[i]=new TH1F(Form("hmm_wo_bg_fom_noZ[%d]",i), "No Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_wo_bg_fom_noZ[i],"No Z cut","","");
	hmm_pi_fom_noZ[i]=new TH1F(Form("hmm_pi_fom_noZ[%d]",i), "No Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_pi_fom_noZ[i],"No Z cut (Pion Selected)","","");
	hmm_pibg_fom_noZ[i]=new TH1F(Form("hmm_pibg_fom_noZ[%d]",i), "No Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_pibg_fom_noZ[i],"No Z cut","","");
	hmm_pi_wobg_fom_noZ[i]=new TH1F(Form("hmm_pi_wobg_fom_noZ[%d]",i), "No Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_pi_wobg_fom_noZ[i],"No Z cut","","");
}//for i

	for (int i=0;i<nth;i++){
		for (int j=0;j<nth;j++){
//			for (int l=0;l<nth;l++){
				int l=0;
				hcoin_k_fom[i][j][l]=new TH1F(Form("hcoin_k_fom[%d][%d][%d]",i,j,l), Form("Cointime (Kaon) z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),bin_coin_c,min_coin_c,max_coin_c);
				set->SetTH1(hcoin_k_fom[i][j][l],Form("Cointime (Kaon) z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),"","");
				hcoin_bg_fom[i][j][l]=new TH1F(Form("hcoin_bg_fom[%d][%d][%d]",i,j,l), Form("Cointime (Kaon) z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),bin_coin_c,min_coin_c,max_coin_c);
				set->SetTH1(hcoin_bg_fom[i][j][l],Form("Cointime (Kaon) z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),"","");
				hcoin_wo_bg_fom[i][j][l]=new TH1F(Form("hcoin_wo_bg_fom[%d][%d][%d]",i,j,l), Form("Cointime (Kaon) z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),bin_coin_c,min_coin_c,max_coin_c);
				set->SetTH1(hcoin_wo_bg_fom[i][j][l],Form("Cointime (Kaon) z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),"","");


				hmm_L_fom[i][j][l]=new TH1F(Form("hmm_L_fom[%d][%d][%d]",i,j,l), Form("Missing Mass z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),bin_mm,min_mm,max_mm);
				set->SetTH1(hmm_L_fom[i][j][l],Form("Missing Mass z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),"","");
				hmm_bg_fom[i][j][l]=new TH1F(Form("hmm_bg_fom[%d][%d][%d]",i,j,l), Form("Missing Mass z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),bin_mm,min_mm,max_mm);
				set->SetTH1(hmm_bg_fom[i][j][l],Form("Missing Mass z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),"","");
				hmm_wo_bg_fom[i][j][l]=new TH1F(Form("hmm_wo_bg_fom[%d][%d][%d]",i,j,l), Form("Missing Mass z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),bin_mm,min_mm,max_mm);
				set->SetTH1(hmm_wo_bg_fom[i][j][l],Form("Missing Mass z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),"","");
				hmm_pi_fom[i][j][l]=new TH1F(Form("hmm_pi_fom[%d][%d][%d]",i,j,l), Form("Missing Mass (Pion Selected) z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),bin_mm,min_mm,max_mm);
				set->SetTH1(hmm_pi_fom[i][j][l],Form("Missing Mass (Pion Selected) z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),"","");
				hmm_pibg_fom[i][j][l]=new TH1F(Form("hmm_pibg_fom[%d][%d][%d]",i,j,l), Form("Missing Mass z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),bin_mm,min_mm,max_mm);
				set->SetTH1(hmm_pibg_fom[i][j][l],Form("Missing Mass z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),"","");
				hmm_pi_wobg_fom[i][j][l]=new TH1F(Form("hmm_pi_wobg_fom[%d][%d][%d]",i,j,l), Form("Missing Mass z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),bin_mm,min_mm,max_mm);
				set->SetTH1(hmm_pi_wobg_fom[i][j][l],Form("Missing Mass z_sum<%lf, z_diff<%lf  cut",zver[i],zver_diff[j]),"","");
//			}//for l
		}//for j
	}//for i




 hmm=new TH1F("hmm","hcoin",bin_mm,min_mm,max_mm);
 set->SetTH1(hmm,"Mass w/o AC tuning","Mass [GeV]","Counts/2 MeV"); 
      hmm->GetXaxis()->SetRangeUser(1.0,1.3);
      hmm->GetYaxis()->SetRangeUser(0.0,450.0);

 hcoin_k=new TH1F("hcoin_k","Coincidence time w/ Correction Kaon Cut  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 hcoin_pi=new TH1F("hcoin_pi","Coincidence time w/ Correction Pion  Cut S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 hcoin_p=new TH1F("hcoin_p","Coincidence time w/ Correction Proton  Cut S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 




 facc_kc=new TF1("facc_kc","[0]",min_coin_c,max_coin_c);
 facc_kc->SetNpx(2000);
 fk_kc=new TF1("fk_kc","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fk_kc->SetNpx(2000);
 fpi_pic=new TF1("fpi_pic","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fpi_pic->SetNpx(2000);
 fp_pc=new TF1("fp_pc","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fp_pc->SetNpx(2000);
 for(int i=0;i<nth;i++){
 fac[i]=new TF1(Form("fac[%d]",i),"[0]",min_coin_c,max_coin_c);
 fac[i]->SetNpx(2000);
 fkk[i]=new TF1(Form("fkk[%d]",i),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fkk[i]->SetNpx(2000);
}
} // MakeHist()

////////////////////////////////////////////////////////////
void tuning::Filling(){

//--------------------------//
	time_t start, end;
	start = time(NULL);
	time(&start);
//--------------------------//

//////////////////////////////
/////mode = H, kine = 1 ////// in hrs_tuning.h
//////////////////////////////
// pathl_off=0.0;
// pathl_off=-498.+30.-3.0+0.5;
// s2_offset=-500.0+25.;
// coin_offset=-41.35+498.;

	int ev = 0;
	for(int k=0;k<ENum;k++){
		//T->GetEntry(k);
		tree->GetEntry(k);

	if(k==ev*100000){
cout << "Event (Fill) : " << k << "/" << ENum << endl;
	ev += 1;
	}

  bool L_Tr = false; // LHRS Tracking Chi2 cut
  bool L_FP = false; // LHRS FP plane cut
  bool R_Tr = false; // RHRS Tracking Chi2 cut
  bool R_FP = false; // RHRS FP plane cut
  bool Kaon = false; // Kaon cut
  bool zcut = false; // z-vertex cut
  bool LHRS = true;  // do LHRS analysis 
  bool RHRS = true;  // do RHRS analysis
/////////////////////
//// Coincidence ////
/////////////////////


	   


    if(LHRS && RHRS && R_evtype==5){
      int NLtr = (int)L_tr_n;  if(NLtr>MAX) NLtr = MAX;
      int NRtr = (int)R_tr_n;  if(NRtr>MAX) NRtr = MAX;
      
      for(int lt=0;lt<NLtr;lt++){
        L_Tr = L_FP = false;
        if( L_tr_chi2[lt]<0.01 ) L_Tr = true;
        if( L_tr_th[lt]<0.17*L_tr_x[lt]+0.025
         && L_tr_th[lt]>0.17*L_tr_x[lt]-0.035
         && L_tr_th[lt]<0.40*L_tr_x[lt]+0.130 ) L_FP = true;
	
        for(int rt=0;rt<NRtr;rt++){
        R_Tr = R_FP = false;
        // FP and chi2 cuts
        if( R_tr_chi2[rt]<0.01 ) R_Tr = true;
        if( R_tr_th[rt]<0.17*R_tr_x[rt]+0.025
         && R_tr_th[rt]>0.17*R_tr_x[rt]-0.035
         && R_tr_th[rt]<0.40*R_tr_x[rt]+0.130 ) R_FP = true;

//#ifdef F1TDC
      convertF1TDCR(param);
      R_s0_t = RS0_F1time[0];
      for(int i=0;i<16;i++){
        if(RS2_F1time[i]>-9999.)R_s2_t[i] =  RS2_F1time[i];
        else R_s2_t[i] = -99.;
      }
//#endif
	  Kaon = false; // Kaon cut 
	  zcut = false; // z-vertex cut
	  
  
	    //---- Initialization ----//
	    tr.Lp[lt] =-100.;
	    tr.Lp[rt] =-100.;
	    tr.Bp     =-100.;
	    tr.dpe     = -100.;
	    tr.dpk[rt] = -100.;
	    tr.dpe_[lt]= -100.;
	    
	    tr.Lp[lt] = L_p;
	    tr.Rp[rt] = R_p;
	    tr.Bp     = B_p;
	    tr.ct_c=-1000.;
	    tr.ct_g=-1000.;
	    tr.pid_cut = 0;
	    tr.ct_cut  = 0;
	    tr.z_cut   = 0;
	    tr.Lp_c[lt] = -100.;
	    tr.Rp_c[rt] = -100.;
	    tr.Bp_c     = -100.;
	    tr.missing_mass=-100000.;
	    tr.coin_time=-1000000.;
	    tr.missing_mass_acc =-100000.;
	    tr.missing_mass_L   =-100000.;
	    tr.missing_mass_nnL =-100000.;
	    tr.missing_mass_H3L =-100000.;
	    tr.missing_mass_cut =-100000.;
	    tr.missing_mass_Al  =-100000.;
	    tr.missing_mass_Lb  =-100000.;
	    tr.missing_mass_nnLb=-100000.;
	    tr.missing_mass_b   =-100000.;
	    tr.missing_mass_Al=-100000.;
	    tr.missing_mass_MgL=-100000.;
	    tr.missing_mass_MgL_acc =-100000.;
	    tr.missing_mass_Al_bg=-100000.;
	    tr.Rpathl=-100.; tr.Lpathl=-100.;
	    tr.Rpathl_c=-100.; tr.Lpathl_c=-100.;
	    ct=-1000.0;

	    tr.AC1_npe_sum=0.0;
	    tr.AC2_npe_sum=0.0;
	    for(int seg=0;seg<24;seg++)
	      tr.AC1_npe[seg]=0.0;
	    for(int seg=0;seg<26;seg++)
	      tr.AC2_npe[seg]=0.0;	  

	    
	  //==== AC ADC convert ch to npe =======//
	  //	  tr.AC1_npe_sum=R_a1_asum_p/400.;
	  //	  tr.AC2_npe_sum=R_a2_asum_p/400.;

	  for(int seg=0;seg<24;seg++){
	    tr.AC1_npe[seg]=AC_npe(1,seg,R_a1_a_p[seg]);
	    tr.AC1_npe_sum+=tr.AC1_npe[seg];
	  }
	  for(int seg=0;seg<26;seg++){
	    tr.AC2_npe[seg]=AC_npe(2,seg,R_a2_a_p[seg]);
	    tr.AC2_npe_sum+=tr.AC2_npe[seg];
	  }    
	  
	  npe_sum_a1->Fill(tr.AC1_npe_sum);
	  //npe_sum_a2->Fill(tr.AC2_npe_sum);




          //Kaon = true; // Without AC CUT
//2020 Okuyama	  	  if( R_a1_asum_p<200 && R_a2_asum_p>1000 && R_a2_asum_p<4000) Kaon = true;
	  // if( R_a1_asum_p<a1_th && R_a2_asum_p>a2_th) Kaon = true;
	  //	  if( R_a1_asum_p<1.0 && R_a2_asum_p>3.0 && R_a2_asum_p<7.0) Kaon = true;	  
	  if( tr.AC1_npe_sum < th1_max && tr.AC2_npe_sum > th2_min) Kaon = true;
	  //	  if(fabs(R_tr_vz[rt])<0.1
	  //         && fabs(L_tr_vz[lt])<0.1 && fabs(R_tr_vz[rt] - L_tr_vz[lt])<0.03)zcut=true;




	    
	  if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1)zcut=true;

	  


	  if( L_Tr && L_FP && R_Tr && R_FP ){



	    B_p     = HALLA_p/1000.0;// [GeV/c]	    
	    L_p     = L_tr_p[lt];
	    R_p     = R_tr_p[rt];
	    
	    //==== Energy Loss calibration ======//

//	    double B_pc,R_pc,L_pc;

	    tr.dpe     = Eloss(0.0,R_tr_vz[0],"B");
	    tr.dpk[rt] = Eloss(R_tr_tg_ph[rt],R_tr_vz[rt],"R");
	    tr.dpe_[lt]= Eloss(L_tr_tg_ph[lt],L_tr_vz[lt],"L");
	  //  
	  //  R_pc = R_p + tr.dpk[rt];
	  //  L_pc = L_p + tr.dpe_[lt];
	  //  B_pc = B_p - tr.dpe;

	    //===================================//	    
//	    double B_E     = sqrt( Me*Me + B_p*B_p );
            int L_s2pad = (int)L_s2_trpad[lt];
            double L_E     = sqrt( Me*Me + L_p*L_p );
//            double L_betae = L_p / sqrt(Me*Me + L_p*L_p);
            int R_s2pad    = (int)R_s2_trpad[rt];
            double R_E     = sqrt( MK*MK + R_p*R_p );
//	    double R_Epi   = sqrt( Mpi*Mpi + R_p*R_p );
            double R_betaK = R_p / sqrt(MK*MK + R_p*R_p);
//	    double R_betaPi =R_p/ sqrt(Mpi*Mpi + R_p*R_p);


	    CoinCalc(R_s2pad,L_s2pad,rt,lt);
	    //	     double test =CoinCalc_c(R_s2pad,L_s2pad,rt,lt);

	     //	     cout<<"ct "<<ct<<" ct_c "<<test<<endl;
	    //	    double ct =CoinCalc(R_s2pad,L_s2pad,rt,lt);



	    
	    //================= ===================== ======================================//
            double L_tgt = L_s2_t[L_s2pad] - (L_tr_pathl[lt] + L_s2_trpath[lt])/c;
            double R_tgt = R_s2_t[R_s2pad] - (R_tr_pathl[rt] + R_s2_trpath[rt])/R_betaK/c;
	   
	    //            double R_tgt_pi = R_s2_t[R_s2pad] - (R_tr_pathl[rt] + R_s2_trpath[rt])/R_betaPi/c;
	    //	    double ct = L_tgt - R_tgt;
	    //ct = L_tgt - R_tgt -1.6; // nnL_small4 
	    //================= ===================== ======================================//


	    if(Kaon)tr.pid_cut=1;
	    if(fabs(ct)<1.0)tr.ct_cut=1;
	    if(zcut)tr.z_cut=1;

            h_ct   ->Fill( ct );
	    h_Rs2  ->Fill(R_tgt);
	    h_Ls2  ->Fill(L_tgt);
	    
            if( Kaon ) h_ct_wK->Fill( ct );
            h_Ls2x_ct ->Fill( ct, L_s2_trx[lt] );
            h_Rs2x_ct ->Fill( ct, R_s2_trx[rt] );
            h_a1sum_ct->Fill( ct, R_a1_asum_p );
            h_a2sum_ct->Fill( ct, R_a2_asum_p );

	    h_Rz->Fill(R_tr_vz[rt]);
	    
	    h_Rth->Fill(R_tr_tg_th[rt]);
	    h_Rph->Fill(R_tr_tg_ph[rt]);
	    h_Rp->Fill(R_p);
	    h_Lz->Fill(L_tr_vz[lt]);
	    h_Lth->Fill(L_tr_tg_th[lt]);	    	    
	    h_Lph->Fill(L_tr_tg_ph[lt]);	    	    
	    h_Lp->Fill(L_p);	    

	     
	    //======== w/o momentum correction ============//

 	    double Ee_b = sqrt( Me*Me + B_p*B_p );
	    double L_Eb = sqrt( Me*Me + L_p*L_p );
	    double R_Eb = sqrt( MK*MK + R_p*R_p );
	    
	    //==== Right Hand Coordinate ========//

	    double R_pz_b=R_p/sqrt(1.0*1.0 + pow(R_tr_tg_th[rt], 2.0) + pow( R_tr_tg_ph[rt],2.0));
	    double R_px_b=R_pz_b * R_tr_tg_th[rt];
	    double R_py_b=R_pz_b * R_tr_tg_ph[rt];
	    double L_pz_b=L_p/sqrt(1.0*1.0 + pow(L_tr_tg_th[lt], 2.0) + pow( L_tr_tg_ph[lt],2.0));
	    double L_px_b=L_pz_b * L_tr_tg_th[lt];
	    double L_py_b=L_pz_b * L_tr_tg_ph[lt];

            TVector3 L_vb, R_vb, B_vb; // Energy loss correction

	    B_vb.SetXYZ(0.0,0.0,B_p);
	    L_vb.SetXYZ(L_px_b, L_py_b, L_pz_b);
	    R_vb.SetXYZ(R_px_b, R_py_b, R_pz_b);
	    R_vb.RotateX(  13.2/180.*PI );
	    L_vb.RotateX( -13.2/180.*PI );




	    
	    
	    double mass_b, mm_b;
//		double mm_Lb;
            mass_b = sqrt( (Ee_b + mt - L_Eb - R_Eb)*(Ee_b + mt - L_Eb - R_Eb)
			 - (B_vb - L_vb - R_vb)*(B_vb - L_vb - R_vb) );

	    mm_b=mass_b - mh;
	    //mm_b=mm_b*1000.; // GeV -> MeV


	    //============================//
	    //=====  calibration =========//
	    //===========================//



	    Calib(rt, lt);

	    //	    tr.ct_c=CoinCalc_c(R_s2pad,L_s2pad,rt,lt);

	    h_Rz_c->Fill(R_tr_vz[rt]);
	    h_Rth_c->Fill(R_tr_tg_th[rt]);
	    h_Rph_c->Fill(R_tr_tg_ph[rt]);
	    h_Rp_c->Fill(R_p);
	    h_Lz_c->Fill(L_tr_vz[lt]);
	    h_Lth_c->Fill(L_tr_tg_th[lt]);	    	    
	    h_Lph_c->Fill(L_tr_tg_ph[lt]);	    	    
	    h_Lp_c->Fill(L_p);
	    tr.Lp_c[lt] = L_p;
	    tr.Rp_c[rt] = R_p;
	    tr.Bp_c     = B_p;
	    



	    //======= W/ Matrix calibraiton ==========================//

            double Ee;

	    Ee =sqrt(B_p*B_p + Me*Me);
	    R_E =sqrt(R_p*R_p + MK*MK);
	    L_E =sqrt(L_p*L_p + Me*Me);


	    //===== Right Hand Coordinate ====//
	    

	    double R_pz = R_p/sqrt(1.0*1.0 + pow(tan(R_tr_tg_th[rt]), 2.0) + pow(tan( R_tr_tg_ph[rt]),2.0) );
	    double R_px = R_pz * tan (R_tr_tg_th[rt] );
	    double R_py = R_pz * tan( R_tr_tg_ph[rt] );

	    double L_pz = L_p/sqrt(1.0*1.0 + pow(tan( L_tr_tg_th[lt] ), 2.0) + pow(tan( L_tr_tg_ph[lt]),2.0));
	    double L_px = L_pz * tan( L_tr_tg_th[lt] );
	    double L_py = L_pz * tan( L_tr_tg_ph[lt] );




            TVector3 L_v, R_v, B_v;
	    B_v.SetXYZ(0.0,0.0,B_p);
	    L_v.SetXYZ(L_px, L_py, L_pz);
	    R_v.SetXYZ(R_px, R_py, R_pz);	    
	    R_v.RotateX(  13.2/180.*PI );
	    L_v.RotateX( -13.2/180.*PI );

	    //======= W/ Matrix & Energy Loss calibraiton ============//

            TVector3 L_vc, R_vc, B_vc;
	    B_vc.SetXYZ(0.0,0.0,B_p);
	    L_vc.SetXYZ(L_px, L_py, L_pz);
	    R_vc.SetXYZ(R_px, R_py, R_pz);
	    R_vc.RotateX(  13.2/180.*PI );
	    L_vc.RotateX( -13.2/180.*PI );
	    double Eec =sqrt(B_p*B_p + Me*Me);
	    double R_Ec =sqrt(R_p*R_p + MK*MK);
	    double L_Ec =sqrt(L_p*L_p + Me*Me);


	   	    
        double mass,mm,mass_L,mass_nnL,mm_L,mm_nnL,mm_Al,mass_Al,mass_MgL;
//		double mass_c, mm_c, mm2, mass2;
	    double mass_pc, mass_H3L,mm_H3L,mm_MgL;

	    
            mass = sqrt( (Ee + mt - L_E - R_E)*(Ee + mt - L_E - R_E)
                              - (B_v - L_v - R_v)*(B_v - L_v - R_v) );

//Changed
        int s2pad = (int)R_s2_trpad[rt];
        double beta = -99.;
		double m2=-99.;
//        if( R_s2_t[s2pad]>0 && R_s0_t>0 && s2pad>=0 ){
	    double p    = R_tr_p[rt];//GeV
        double path = R_s2_trpath[rt] - R_s0_trpath[rt];//m
        beta = path / ( R_s2_t[s2pad] - R_s0_t ) / c/2.;
		//cout<<"p="<<p<<endl;
		//cout<<"c="<<c<<endl;
		//cout<<"path="<<path<<endl;
//		cout<<"S2time="<<R_s2_t[s2pad]<<endl;
//		cout<<"S0time="<<R_s0_t<<endl;
//		cout<<"beta="<<beta<<endl;
        //double beta=R_tr_p[rt]/sqrt(R_tr_p[rt]*R_tr_p[rt]+MK*MK);
        m2 = ( 1./beta/beta - 1. ) * p * p;
		if(k%100000==0){cout<<"m2="<<m2<<endl;}
//		cout<<"m2="<<m2<<endl;
//			double beta2 = R_v*R_v/(R_E*R_E);//beta^2
//			m2 = R_v*R_v*(1./beta2-1.);
//		}

   
//            mass2= sqrt( (Ee + mt - L_E - R_Epi)*(Ee + mt - L_E - R_Epi)
//                              - (B_v - L_v - R_v)*(B_v - L_v - R_v) );

	    
            mass_pc = sqrt( (Eec + mt - L_Ec - R_Ec)*(Eec + mt - L_Ec - R_Ec)
                              - (B_vc - L_vc - R_vc)*(B_vc - L_vc - R_vc) );


	    
	    mm=mass - mh;
            //mm2=mass2 - mh;

	   // mm = mm*1000.; // GeV -> MeV
	    //mm2 = mm2*100.; // GeV ->MeV 

	    //=== w/ matrix tuning ======//
	    
	    // Lambda Mass //
           mass_L = sqrt( (Ee + Mp - L_E - R_E)*(Ee + Mp - L_E - R_E)
                              - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
	   mm_L=mass_L - ML;
	   mm_L = mm_L*1000.;
	    // nnL Mass //
           mass_nnL = sqrt( (Ee + MT - L_E - R_E)*(Ee + MT - L_E - R_E)
                              - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
	   mm_nnL=mass_nnL - MnnL;
	   mm_nnL = mm_nnL*1000.;
	    // H3L Mass //
           mass_H3L = sqrt( (Ee + MHe3 - L_E - R_E)*(Ee + MHe3 - L_E - R_E)
                              - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
	   mm_H3L=mass_H3L - MH3L;	   
	   mm_H3L = mm_H3L*1000.;
	   
	    // Alminium Mass //
           mass_Al = sqrt( (Ee + MAl - L_E - R_E)*(Ee + MAl - L_E - R_E)
                              - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
	   mm_Al=mass_Al - MAl;
	   mm_Al = mm_Al*1000.;
	   
	   // Mg27L Mass //
           mass_MgL = sqrt( (Ee + MAl - L_E - R_E)*(Ee + MAl - L_E - R_E)
                              - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
	   mm_MgL=mass_MgL - MMgL;	   
	   mm_MgL = mm_MgL*1000.;
	   
	    
	    if( Kaon && (fabs(ct-30.)<10. || fabs(ct+30.)<10.) ){
              h_mmallbg->Fill( mm );
              if( fabs( L_tr_vz[lt] + 0.125 ) < 0.015 || fabs( L_tr_vz[lt] - 0.125 ) < 0.015 ){ 
                h_mmfoilbg->Fill( mm );
              }
	      //              if( fabs( L_tr_vz[lt] ) < 0.1 ){
	      if(zcut ){ 
                h_mmbg->Fill( mm );
              }
		  //}
	    }




//-------------------------------------------//
//-------------------------------------------//






////-------------------------------------------//
////-------------AC tuning---------------------//
////-------------------------------------------//



	if(zcut){//no AC cut
	hcoin_k_fom_noAC->Fill(ct);
	//h_m2_mm->Fill(m2,mm);
	//h_m2_mm->Fill(tr.AC2_npe_sum,mm);
	//h_m2_ac->Fill(m2,tr.AC2_npe_sum);
				if((-100.<ct && ct <-20.) || (20.<ct && ct<100.)){
					double ct_ = ct;
				        while(1){
					  if(-20.<ct && ct<20.){
						 hcoin_bg_fom_noAC->Fill(ct);
						 hmm_bg_fom_noAC->Fill(mm);
						 hmm_pibg_fom_noAC->Fill(mm);
						 break;}
					       else if(ct<-20.){ct=ct+40.;}
					       else if(20.<ct){ct=ct-40.;}
					 }
					ct = ct_;
					}//cointime
					if(fabs(ct)<1.){
						hmm_L_fom_noAC->Fill(mm);
				}
					if(fabs(ct-3.05)<0.7){
						hmm_pi_fom_noAC->Fill(mm);//MM if pion
					}
	}


//-------------------------------------------//
//-------------No AC1 cut--------------------//
//-------------------------------------------//
	cut_ac2=false;
   	//if(ac2l_adc[30]<tr.AC2_npe_sum && tr.AC2_npe_sum < ac2u_adc[16])cut_ac2=true;
   	if(3.<tr.AC2_npe_sum && tr.AC2_npe_sum < 20.)cut_ac2=true;
	if(zcut && cut_ac2){//no AC1 cut
	hcoin_k_fom_noAC1->Fill(ct);
				if((-100.<ct && ct <-20.) || (20.<ct && ct<100.)){
					double ct_ = ct;
				        while(1){
					  if(-20.<ct && ct<20.){
						 hcoin_bg_fom_noAC1->Fill(ct);
						 hmm_bg_fom_noAC1->Fill(mm);
						 break;}
					       else if(ct<-20.){ct=ct+40.;}
					       else if(20.<ct){ct=ct-40.;}
					 }
					ct = ct_;
					}//cointime
					if(fabs(ct)<1.){
						hmm_L_fom_noAC1->Fill(mm);
				}
	}
//-------------------------------------------//
//-------------No AC2 cut--------------------//
//-------------------------------------------//

				cut_ac1=false;
		    	//if(tr.AC1_npe_sum < ac1_adc[16])cut_ac1=true;
		    	if(tr.AC1_npe_sum < 3.75 && tr.AC2_npe_sum > 0.)cut_ac1=true;
	if(zcut&&cut_ac1){//no AC2 cut
	hcoin_k_fom_noAC2->Fill(ct);
				if((-100.<ct && ct <-20.) || (20.<ct&&ct<100.)){
					double ct_ = ct;
				        while(1){
					  if(-20.<ct && ct<20.){
						 hcoin_bg_fom_noAC2->Fill(ct);
						 hmm_bg_fom_noAC2->Fill(mm);
						 break;}
					       else if(ct<-20.){ct=ct+40.;}
					       else if(20.<ct){ct=ct-40.;}
					 }
					ct = ct_;
					}//cointime

					if(fabs(ct)<1.){
						hmm_L_fom_noAC2->Fill(mm);
				}
	}

//--------------Z vertex---------------//
	h_L_vz->Fill(R_tr_vz[rt]+L_tr_vz[lt]);
	h_R_vz->Fill(R_tr_vz[rt]-L_tr_vz[lt]);



//-------------------------------------------//
//-------------No Z cut        --------------//
//-------------------------------------------//
for(int i=0;i<nth;i++){
	cut_ac1=false;
	cut_ac2=false;
	zcut=false;
	if(tr.AC1_npe_sum<3.75)cut_ac1=true;
	if(tr.AC2_npe_sum >3. && tr.AC2_npe_sum < 20.)cut_ac2=true;
	if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<(zver_diff[i]/1000))zcut=true;
	if(zcut && cut_ac1 && cut_ac2){//no AC cut
	hcoin_k_fom_noZ[i]->Fill(ct);
				if((-100.<ct && ct <-20.) || (20.<ct && ct<100.)){
					double ct_ = ct;
				        while(1){
					  if(-20.<ct && ct<20.){
						 hcoin_bg_fom_noZ[i]->Fill(ct);
						 hmm_bg_fom_noZ[i]->Fill(mm);
						 hmm_pibg_fom_noZ[i]->Fill(mm);
						 break;}
					       else if(ct<-20.){ct=ct+40.;}
					       else if(20.<ct){ct=ct-40.;}
					 }
					ct = ct_;
					}//cointime
					if(fabs(ct)<1.){
						hmm_L_fom_noZ[i]->Fill(mm);
				}
					if(fabs(ct-3.05)<0.7){
						hmm_pi_fom_noZ[i]->Fill(mm);//MM if pion
					}
	}
}

//-------------------------------------------------------------------//
//-------------------------------------------------------------------//
	for (int i=0;i<nth;i++){
		for (int j=0;j<nth;j++){
//			for (int l=0;l<nth;l++){
				int l=0;
//cout<<"j="<<j<<endl;

//				ac1_adc[i]=3.75;//to be consistent with No AC2 cut
//				ac2l_adc[j]=0.;//to be consistent with No AC1 cut
//Slice
				ac2u_adc[l]=20.;//to be consistent with No AC1 cut


				cut_ac1=false;
				cut_ac2=false;
				zcut=false;
//Okuyama				//if(R_a1_asum_p<ac1_adc[i])cut_ac1=true;
				if(tr.AC1_npe_sum<3.75)cut_ac1=true;
				//if(ac2l_adc[j]<R_a2_asum_p && R_a2_asum_p < 4000)cut_ac2=true;
//Okuyama		    	//if(ac2l_adc[j]<R_a2_asum_p && R_a2_asum_p < ac2u_adc[l])cut_ac2=true;
		    	if(tr.AC2_npe_sum > 3. && tr.AC2_npe_sum < 20.)cut_ac2=true;
				if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<(zver_diff[i]/1000) && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<(zver[j]/1000))zcut=true;
//Changed			
//				if(ac2l_adc[j]<tr.AC2_npe_sum && tr.AC2_npe_sum < ac2l_adc[j]+10.)cut_ac2=true;
				//if( zcut && cut_ac1 && cut_ac2){
				//}
				if( zcut && cut_ac1 && cut_ac2){
				hcoin_k_fom[i][j][l]->Fill(ct);
			    //cout<<"hcoin_k_fom is filled" << endl;


		//-------------------------------------------//
		//---------Accidental B.G.-------------------//
		//-------------------------------------------//
				if((-100.<ct && ct <-20.) || (20.<ct&&ct<100.)){
					double ct_ = ct;
				        while(1){
					  if(-20.<ct && ct<20.){
						 hcoin_bg_fom[i][j][l]->Fill(ct);
						 hmm_bg_fom[i][j][l]->Fill(mm);
						 hmm_pibg_fom[i][j][l]->Fill(mm);//MM if pion
						 break;}
					       else if(ct<-20.){ct=ct+40.;}
					       else if(20.<ct){ct=ct-40.;}
					 }
						  ct = ct_;
				}
		//-------------------------------------------//



					if(fabs(ct)<1.){
					//if(fabs(ct-mean_k[i][j][l])<sig_k[i][j][l]){
						hmm_L_fom[i][j][l]->Fill(mm);
					}//cointime
					if(fabs(ct-3.05)<0.7){
						// def_sig_pi=0.443; def_mean_pi=3.0;
						hmm_pi_fom[i][j][l]->Fill(mm);//MM if pion
					}
					//}
				}//if cut condition
//			}//for l
		}//for j
	}//for i	
//-------------------------------------------//
//-------------------------------------------//


	    if( Kaon && fabs(ct)<1. ){	      
              h_mmall ->Fill( mm );
              if( fabs( L_tr_vz[lt] + 0.125 ) < 0.015 || fabs( L_tr_vz[lt] - 0.125 ) < 0.015 ){ 
		tr.missing_mass_Al=mm_Al;
		tr.missing_mass_MgL=mm_MgL;
		
		h_mmfoil->Fill( mm );
		
              }
              if( fabs( L_tr_vz[lt]  ) < 0.1 ){ 
                h_Lp_mm   ->Fill( mm, L_tr_p[lt] );
                h_Ll_mm   ->Fill( mm, L_tr_pathl[lt] );
                h_Ltgy_mm ->Fill( mm, L_tr_tg_y[lt] );
                h_Ltgth_mm->Fill( mm, L_tr_tg_th[lt] );
                h_Ltgph_mm->Fill( mm, L_tr_tg_ph[lt] );
                h_Lvx_mm  ->Fill( mm, L_tr_vx[lt] );
                h_Lvy_mm  ->Fill( mm, L_tr_vy[lt] );
                h_Lvz_mm  ->Fill( mm, L_tr_vz[lt] );
                h_Lx_mm   ->Fill( mm, L_tr_x[lt] );
                h_Ly_mm   ->Fill( mm, L_tr_y[lt] );
                h_Lth_mm  ->Fill( mm, L_tr_th[lt] );
                h_Lph_mm  ->Fill( mm, L_tr_ph[lt] );
              }
              if( fabs( R_tr_vz[rt] ) < 0.1 ){ 
                h_Rp_mm   ->Fill( mm, R_tr_p[rt] );
                h_Rl_mm   ->Fill( mm, R_tr_pathl[rt] );
                h_Rtgy_mm ->Fill( mm, R_tr_tg_y[rt] );
                h_Rtgth_mm->Fill( mm, R_tr_tg_th[rt] );
                h_Rtgph_mm->Fill( mm, R_tr_tg_ph[rt] );
                h_Rvx_mm  ->Fill( mm, R_tr_vx[rt] );
                h_Rvy_mm  ->Fill( mm, R_tr_vy[rt] );
                h_Rvz_mm  ->Fill( mm, R_tr_vz[rt] );
                h_Rx_mm   ->Fill( mm, R_tr_x[rt] );
                h_Ry_mm   ->Fill( mm, R_tr_y[rt] );
                h_Rth_mm  ->Fill( mm, R_tr_th[rt] );
                h_Rph_mm  ->Fill( mm, R_tr_ph[rt] );
                h_Rp_Lp   ->Fill( L_tr_p[lt], R_tr_p[rt] );
                h_ct_Rp->Fill(R_tr_p[rt],ct);
              }

	      //	      if(fabs(R_tr_vz[rt]-0.125)<0.01 || fabs(R_tr_vz[rt] +0.125 )<0.01){
	      if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && (fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0 >0.125 ) ){
		tr.missing_mass_Al_bg=mm;
		h_mm_Al_bg->Fill(mm);
		h_Rz_cut->Fill(R_tr_vz[rt]);
	      }
	      //}

	      	      
	      //              if( fabs( L_tr_vz[lt] ) < 0.1 && fabs( R_tr_vz[rt] ) < 0.1 ){
	      if(zcut){

		//                h_mm      ->Fill( mm );

	       tr.missing_mass_cut = mm;
	       tr.missing_mass_L = mm_L;
	       tr.missing_mass_nnL = mm_nnL;
	       tr.missing_mass_H3L = mm_H3L;
	       tr.missing_mass_b=mm_b;


	       
                h_mm_L    ->Fill( mm_L );
                h_mm_L_ec    ->Fill( mass_pc);		
                h_mm_nnL  ->Fill( mm_nnL );
		h_mm_H3L  ->Fill( mm_H3L );
                h_ct_wK_z->Fill( ct );                
        
	      }
	      //}
	    
		    


	    
	    } // if Kaon

				    

              if((Kaon && fabs(ct)<1.0) && ((-0.15<(L_tr_vz[lt]) && (L_tr_vz[lt])<-0.1) || (( 0.1<(L_tr_vz[lt])) && (L_tr_vz[lt])<0.15)) &&  (fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025) && ((-0.15<(R_tr_vz[rt]) && (R_tr_vz[rt])<-0.1) ||( 0.1<(R_tr_vz[rt]) && (R_tr_vz[rt])<0.15)))h_mm_MgL->Fill(mm_MgL);//h_mm_Al->Fill(mm_Al);

	      if((Kaon) && (((-35<ct) && (ct<-15.0)) || ((15.0<ct) && (ct<35))) && (((-0.15<(L_tr_vz[lt])) && ((L_tr_vz[lt])<-0.1)) || ( (0.1<(L_tr_vz[lt])) && ((L_tr_vz[lt])<0.15))) && (fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025) && (((-0.15<(R_tr_vz[rt]-0.01)) && ((R_tr_vz[rt])<-0.1)) ||( (0.1<(R_tr_vz[rt])) && ((R_tr_vz[rt])<0.15)))){
		tr.missing_mass_MgL_acc=mm_MgL;
		
		h_mm_MgL_acc->Fill(mm_MgL);
	      }

	      
	      if( Kaon && ((-35<ct && ct<-15.0) || (15.0<ct && ct<35)) && zcut){
		 //		 fabs( L_tr_vz[lt] ) < 0.1 && fabs( R_tr_vz[rt] ) < 0.1 &&fabs( L_tr_vz[lt] ) < 0.1){
                h_acc_nnL     ->Fill(mm_nnL);
		h_acc_H3L     ->Fill(mm_H3L);
                h_acc_L       ->Fill(mm_L);
                h_ct_wK_z_acc ->Fill( ct );
	     //}
	     }

	 
              double ctime=-1000.;
	     //--------------------------------------------------------------------------------//
              if( Kaon && zcut){
		  //		  fabs( L_tr_vz[lt] ) < 0.1 && fabs( R_tr_vz[rt] ) < 0.1 &&fabs( L_tr_vz[lt] ) < 0.1){
               h_ct_wK_z_all->Fill(ct);
            

              if((-35<ct && ct <-15) || (15<ct && ct<53)){
	     
	       ctime=ct;
	       
              while(1){
	       if(-1.0<ctime && ctime<1.0){
		 h_ct_acc->Fill(ctime);
                 h_ct_acc->Fill(ctime-36);
		 break;}
	       else if(ctime<-1.0){ctime=ctime+2;}
	       else if(1.0<ctime){ctime=ctime-2;}
	      }
	      }
	      //}
	      }
	
	
             tr.missing_mass = mm          ; tr.coin_time =ct         ;
	     tr.momR         = R_tr_p[0]  ; tr.momL      =L_tr_p[0] ;
	     tr.zR           = R_tr_vz[0] ; tr.zL        =L_tr_vz[0];
	     //	     tr.AC1_sum      = R_a1_asum_p/400. ; tr.AC2_sum   =R_a2_asum_p/400.;
	     tr.AC1_sum      = R_a1_asum_p ; tr.AC2_sum   =R_a2_asum_p;
	     tr.ct_acc=ctime;
	    	    // tree_out->Fill();
	  
    	      //--------------------------------------------------------------------------------------//


	     //	     if( fabs( L_tr_vz[lt]  ) < 0.1 && fabs( R_tr_vz[rt]  ) < 0.1 &&fabs(ct)<1.0)
	     if( zcut && fabs(ct)<1.0)
	       h_mm->Fill( mm ); //No Kaon Cut
	     //	     if( fabs( L_tr_vz[lt]  ) < 0.1 && fabs( R_tr_vz[rt]  ) < 0.1 && 2.0<ct && ct<4.0)
	     if( zcut && 2.0<ct && ct<4.0)
	       h_mm_pi->Fill( mm ); //No Kaon Cut
	     //	     if( fabs( L_tr_vz[lt]  ) < 0.1 && fabs( R_tr_vz[rt]  ) < 0.1
	     if(  zcut && ((-35<ct && ct<-15.0) || (15.0<ct && ct<35))){
	       h_mm_acc->Fill( mm ); //No Kaon Cut
	       tr.missing_mass_acc=mm;
	     }




////MISSING MASS w/o cut////////////
	if(ct<1.0 && -1.0<ct) hmm->Fill(mm);




          } // if L_Tr && L_FP && R_Tr && R_FP
        } // for NRtr
      } // for NLtr
	  //tree_out->Fill();
    } // if LHRS && RHRS
    if(k%100000==0){
      end = time(NULL);
      time(&end);
      double diff = difftime(end,start);
      double esttime = diff * ENum / (k+1) - diff;
      cout<<k<<" / "<<ENum<<" : "<<Form("%.0lf sec passed,  %.0lf sec left",diff,esttime)<<endl;
    }
  } // for ENum


}



////////////////////////////////////////////////////////////


void tuning::Fitting(){
//OLD
 //def_sig_p=0.852; def_mean_p=0.0;
 //def_sig_pi=0.443; def_mean_pi=11;
 //def_sig_k=0.644; def_mean_k=8.;
 //def_acc=27.7;

 def_sig_p=0.852; def_mean_p=-8.0;
 def_sig_pi=0.443; def_mean_pi=3.0;
 def_sig_k=0.644; def_mean_k=0.0;
 def_acc=27.7;

// fpi_pic->SetParameter(1,def_mean_pi);
// fpi_pic->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
// fpi_pic->SetParameter(2,def_sig_pi);
// fpi_pic->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
// fpi_pic->SetParameter(3,def_acc);
//
// fp_pc->SetParameter(1,def_mean_p);
// fp_pc->SetParLimits(1,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
// fp_pc->SetParameter(2,def_sig_p);
// fp_pc->SetParLimits(2,0.8*def_sig_p,1.2*def_sig_p);
// fp_pc->SetParameter(3,def_acc);
//
//// hcoin_k->Fit("facc_kc","Rq","",min_coin_c,min_coin_c+3.0);
// h_ct_wK->Fit("facc_kc","Rq","",min_coin_c,min_coin_c+3.0);
// def_acc_k=facc_kc->GetParameter(0);
//
// fk_kc->SetParameter(1,def_mean_k);
// fk_kc->SetParLimits(1,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
// fk_kc->SetParameter(2,def_sig_k);
// fk_kc->SetParLimits(2,0.8*def_sig_k,1.2*def_sig_k);
// fk_kc->FixParameter(3,def_acc_k);
//
//// hcoin_k->Fit("fk_kc","Rq","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
// h_ct_wK->Fit("fk_kc","Rq","",-1,1);
// h_ct_wK->Fit("fpi_pic","Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
// h_ct_wK->Fit("fp_pc","Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
//
// def_num_k=fk_kc->GetParameter(0);
// def_mean_k=fk_kc->GetParameter(1);
// def_sig_k=fk_kc->GetParameter(2);
//
// def_num_p=fp_pc->GetParameter(0);
// def_mean_p=fp_pc->GetParameter(1);
// def_sig_p=fp_pc->GetParameter(2);
//
// def_num_pi=fpi_pic->GetParameter(0);
// def_mean_pi=fpi_pic->GetParameter(1);
// def_sig_pi=fpi_pic->GetParameter(2);
 
//for(int i=0;i<nth;i++){
// noise[i]=signal[i]=0.;
//cout << "hcoin_k_ac2[" <<i<<"] fitting start" << endl;
//cout << min_coin_c << endl; 
// hcoin_k_ac2[i]->Fit(Form("fac[%d]",i),"Rq","",min_coin_c,min_coin_c+3.0);
//cout << "hcoin_k_ac2[" <<i<<"] fitting end" << endl;
// noise[i]=fac[i]->GetParameter(0);
//
// fkk[i]->SetParameter(1,def_mean_k);
// fkk[i]->SetParLimits(1,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
// fkk[i]->SetParameter(2,def_sig_k);
// fkk[i]->SetParLimits(2,0.8*def_sig_k,1.2*def_sig_k);
// fkk[i]->FixParameter(3,noise[i]);
//
// hcoin_k_ac2[i]->Fit(Form("fkk[%d]",i),"Rq","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
//
// signal[i]=fkk[i]->GetParameter(0);
//cout << "Cut[" << i << "] " << "S = " << signal[i] << "/ N = " << noise[i] << "... S*S/N = " << signal[i]*signal[i]/noise[i] << endl;
//}//for i
}//Fitting
//////////////////////////////////////////////////////////////////

void tuning::ACtune(){ 
  cout<<"==========================="<<endl;
  cout<<"========= Tuning =========="<<endl;
  cout<<"==========================="<<endl;

 def_sig_L=0.003; def_mean_L=0.0;
 def_sig_S=0.004; def_mean_S=MS0-ML;



//-----Main FOM histgram------//
 int h3_fom_bin_ac1 = 10;
 int h3_fom_bin_ac2l = 10;
 int h3_fom_bin_ac2u = 10;
 //h3_fom = new TH3D("h3_fom","",10,400.,600.,10,600.,1050.,25,1500.,4000.);
 h3_fom = new TH3D("h3_fom","",h3_fom_bin_ac1,ac1_adc[0],ac1_adc[9],h3_fom_bin_ac2l,ac2l_adc[0],ac2l_adc[9],h3_fom_bin_ac2u,ac2u_adc[0],ac2u_adc[9]);
 //h3_fom = new TH3D("h3_fom","",10,400.,580.,10,600.,1050.,10,1500.,4000.);
 

 //-------------------------------------------------//
 //----No AC1 cut---//
 //Efficiency 100%//
 
//cout<<"BG subtraction cointime"<<endl;
 hcoin_bg_fom_noAC1->Scale(40./160.);
 hcoin_wo_bg_fom_noAC1->Add(hcoin_k_fom_noAC1,hcoin_bg_fom_noAC1,1.0,-1.0);
 fp_noAC1=new TF1("fp_noAC1","gausn(0)",min_coin_c,max_coin_c);
 fp_noAC1->SetNpx(2000);
 fpi_noAC1 =new TF1("fpi_noAC1","gausn(0)",min_coin_c,max_coin_c);
 fpi_noAC1->SetNpx(2000);
 fk_noAC1=new TF1("fk_noAC1","gausn(0)",min_coin_c,max_coin_c);
 fk_noAC1->SetNpx(2000);
//cout<<"fp fit start"<<endl;
 hcoin_wo_bg_fom_noAC1->Fit("fp_noAC1","Rq","0",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
// n_p_noAC1=fp_noAC1->GetParameter(0);
 mean_p_noAC1=fp_noAC1->GetParameter(1);
 sig_p_noAC1=fp_noAC1->GetParameter(2);
 center_p=mean_p_noAC1;
 range_p=1.435;
 n_p_noAC1=hcoin_wo_bg_fom_noAC1->Integral(hcoin_wo_bg_fom_noAC1->FindBin(center_p-range_p),hcoin_wo_bg_fom_noAC1->FindBin(center_p+range_p));
//cout<<"fpi fit start"<<endl;
 hcoin_wo_bg_fom_noAC1->Fit("fpi_noAC1","Rq","0",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
// n_pi_noAC1=fpi_noAC1->GetParameter(0);
 mean_pi_noAC1=fpi_noAC1->GetParameter(1);
 sig_pi_noAC1=fpi_noAC1->GetParameter(2);
 center_pi=mean_pi_noAC1;
 range_pi=0.7;
 n_pi_noAC1=hcoin_wo_bg_fom_noAC1->Integral(hcoin_wo_bg_fom_noAC1->FindBin(center_pi-range_pi),hcoin_wo_bg_fom_noAC1->FindBin(center_pi+range_pi));
//cout<<"fk fit start"<<endl;
 //hcoin_wo_bg_fom_noAC1->Fit("fk_noAC1","Rq","0",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_wo_bg_fom_noAC1->Fit("fk_noAC1","Rq","0",-1,1);
 //n_k_noAC1=fk_noAC1->GetParameter(0);
 mean_k_noAC1=fk_noAC1->GetParameter(1);
// mean_k_noAC1=0;
 sig_k_noAC1=fk_noAC1->GetParameter(2);

 center_k=mean_k_noAC1;
 range_k=1.;
 n_k_noAC1=hcoin_wo_bg_fom_noAC1->Integral(hcoin_wo_bg_fom_noAC1->FindBin(center_k-range_k),hcoin_wo_bg_fom_noAC1->FindBin(center_k+range_k));

 //-----Fitting as a whole function---------//

 fcoin_noAC1 =new TF1("fcoin_noAC1","gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin_noAC1->SetNpx(2000);
 fcoin_noAC1->SetTitle("Cointime w/o AC1 cut;Cointime [ns];Counts [1/56 ns]");
 fcoin_noAC1->SetParameters(n_pi_noAC1,mean_pi_noAC1,sig_pi_noAC1,n_k_noAC1,mean_k_noAC1,sig_k_noAC1,n_p_noAC1,mean_p_noAC1,sig_p_noAC1);
 //hcoin_wo_bg_fom_noAC1->Fit("fcoin_woAC1","Rq","0",min_coin_c,max_coin_c);
 //n_pi_noAC1=fcoin_noAC1->GetParameter(0);//Npi_nocut
 //mean_pi_noAC1=fcoin_noAC1->GetParameter(1);
 //sig_pi_noAC1=fcoin_noAC1->GetParameter(2);
 //n_k_noAC1=fcoin_noAC1->GetParameter(3);//Nk_nocut
 //mean_k_noAC1=fcoin_noAC1->GetParameter(4);
 //sig_k_noAC1=fcoin_noAC1->GetParameter(5);
 //n_p_noAC1=fcoin_noAC1->GetParameter(6);//Np_nocut
 //mean_p_noAC1=fcoin_noAC1->GetParameter(7);
 //sig_p_noAC1=fcoin_noAC1->GetParameter(8);
cout<<"n_pi_noAC1="<<n_pi_noAC1<<"n_k_noAC1="<<n_k_noAC1<<"n_p_noAC1="<<n_p_noAC1
<<"mean_pi_noAC1="<<mean_pi_noAC1<<"sig_pi_noAC1="<<sig_pi_noAC1<<"mean_k_noAC1="<<mean_k_noAC1<<"sig_k_noAC1="<<sig_k_noAC1<<"mean_p_noAC1="<<mean_p_noAC1<<"sig_p_noAC1="<<sig_p_noAC1<<endl;

 //------- Get Error Paramters ---//
 n_pi_err_noAC1=fpi_noAC1->GetParError(0); 
 n_k_err_noAC1=fk_noAC1->GetParError(3); 
 n_p_err_noAC1=fp_noAC1->GetParError(6); 
 

 hmm_bg_fom_noAC1->Scale(1./80.);
 hmm_wo_bg_fom_noAC1->Add(hmm_L_fom_noAC1,hmm_bg_fom_noAC1,1.0,-1.0);


// fmmbg_noAC1=new TF1("fmmbg_noAC1","pol3(0)",min_mm,max_mm);
// fmmbg_noAC1->SetNpx(2000);
 fL_noAC1=new TF1("fL_noAC1","gausn(0)",min_mm,max_mm);
 fL_noAC1->SetNpx(2000);
 fL_noAC1->SetParLimits(2,0.,0.01);
 fS_noAC1=new TF1("fS_noAC1","gausn(0)",min_mm,max_mm);
 fS_noAC1->SetNpx(2000);
 fS_noAC1->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
 
 fmm_noAC1=new TF1("fmm_noAC1","gausn(0)+gausn(3)",min_mm,max_mm);
 fmm_noAC1->SetNpx(2000);
 fmm_noAC1->SetTitle("Missing Mass w/o AC cut;Coin time [ns];Counts [1/56 ns]");



 //------- Fitting ----------//

//cout<<"fmmbg fit start"<<endl;
// hmm_bg_fom_noAC1->Fit("fmmbg_noAC1","Rq","",min_mm,max_mm);
// double d = fmmbg_noAC1->GetParameter(0);
// double c = fmmbg_noAC1->GetParameter(1);
// double b = fmmbg_noAC1->GetParameter(2);
// double a = fmmbg_noAC1->GetParameter(3);

//cout<<"fL fit start"<<endl;
 hmm_wo_bg_fom_noAC1->Fit("fL_noAC1","Rq","",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
// n_L_noAC1=fL_noAC1->GetParameter(0);
 mean_L_noAC1=fL_noAC1->GetParameter(1);
 sig_L_noAC1=fL_noAC1->GetParameter(2);
 mean_L_noAC1=def_mean_L;
 sig_L_noAC1=def_sig_L;
 center_L=mean_L_noAC1;
 range_L=2*sig_L_noAC1;
 n_L_noAC1=hmm_wo_bg_fom_noAC1->Integral(hmm_wo_bg_fom_noAC1->FindBin(center_L-range_L),hmm_wo_bg_fom_noAC1->FindBin(center_L+range_L));
 cout<<"n_L"<<n_L_noAC1<<endl;

//cout<<"fS fit start"<<endl;
 hmm_wo_bg_fom_noAC1->Fit("fS_noAC1","Rq","",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
// n_S_noAC1=fS_noAC1->GetParameter(0);
 mean_S_noAC1=fS_noAC1->GetParameter(1);
 mean_S_noAC1=def_mean_S;
 sig_S_noAC1=def_sig_S;
 sig_S_noAC1=fS_noAC1->GetParameter(2);
 center_S=mean_S_noAC1;
 range_S=2*sig_S_noAC1;
 n_S_noAC1=hmm_wo_bg_fom_noAC1->Integral(hmm_wo_bg_fom_noAC1->FindBin(center_S-range_S),hmm_wo_bg_fom_noAC1->FindBin(center_S+range_S));
 cout<<"n_S"<<n_S_noAC1<<endl;
//cout<<"mean_S"<<mean_S_noAC1<<endl;
//cout<<"sig_S"<<sig_S_noAC1<<endl;


 //--Fitting again as a total function--//
// fmm_noAC1->SetParameters(n_L_noAC1,mean_L_noAC1,sig_L_noAC1,n_S_noAC1,mean_S_noAC1,sig_S_noAC1,d,c,b,a);
//// hmm_L_fom_noAC1->Fit(Form("fmm[%d][%d][%d]",i,j,l),"Rq","",min_mm,max_mm);
// n_L_noAC1=fmm_noAC1->GetParameter(0);
// mean_L_noAC1=fmm_noAC1->GetParameter(1);
// sig_L_noAC1=fmm_noAC1->GetParameter(2);
// n_S_noAC1=fmm_noAC1->GetParameter(3);
// mean_S_noAC1=fmm_noAC1->GetParameter(4);
// sig_S_noAC1=fmm_noAC1->GetParameter(5);
// //------- Get Error Paramters ---//
// n_L_err_noAC1=fmm_noAC1->GetParError(0); 
// n_S_err_noAC1=fmm_noAC1->GetParError(3); 
//
//
//
// signal_noAC1=noise_noAC1=0.;
// d = fmm_noAC1->GetParError(6); 
// c = fmm_noAC1->GetParError(7); 
// b = fmm_noAC1->GetParError(8); 
// a = fmm_noAC1->GetParError(9); 
 //-------------------------------------------------//
 //-------------------------------------------------//
 //----No AC2 cut---//
 
//cout<<"BG subtraction cointime"<<endl;
 hcoin_bg_fom_noAC2->Scale(40./160.);
 hcoin_wo_bg_fom_noAC2->Add(hcoin_k_fom_noAC2,hcoin_bg_fom_noAC2,1.0,-1.0);
 fp_noAC2=new TF1("fp_noAC2","gausn(0)",min_coin_c,max_coin_c);
 fp_noAC2->SetNpx(2000);
 fpi_noAC2 =new TF1("fpi_noAC2","gausn(0)",min_coin_c,max_coin_c);
 fpi_noAC2->SetNpx(2000);
 fk_noAC2=new TF1("fk_noAC2","gausn(0)",min_coin_c,max_coin_c);
 fk_noAC2->SetNpx(2000);
//cout<<"fp fit start"<<endl;
 hcoin_wo_bg_fom_noAC2->Fit("fp_noAC2","Rq","0",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
// n_p_noAC2=fp_noAC2->GetParameter(0);
 mean_p_noAC2=fp_noAC2->GetParameter(1);
 sig_p_noAC2=fp_noAC2->GetParameter(2);
 n_p_noAC2=hcoin_wo_bg_fom_noAC2->Integral(hcoin_wo_bg_fom_noAC2->FindBin(center_p-range_p),hcoin_wo_bg_fom_noAC2->FindBin(center_p+range_p));
//cout<<"fpi fit start"<<endl;
 hcoin_wo_bg_fom_noAC2->Fit("fpi_noAC2","Rq","0",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
// n_pi_noAC2=fpi_noAC2->GetParameter(0);
 mean_pi_noAC2=fpi_noAC2->GetParameter(1);
 sig_pi_noAC2=fpi_noAC2->GetParameter(2);
 n_pi_noAC2=hcoin_wo_bg_fom_noAC2->Integral(hcoin_wo_bg_fom_noAC2->FindBin(center_pi-range_pi),hcoin_wo_bg_fom_noAC2->FindBin(center_pi+range_pi));
//cout<<"fk fit start"<<endl;
 //hcoin_wo_bg_fom_noAC2->Fit("fk_noAC2","Rq","0",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_wo_bg_fom_noAC2->Fit("fk_noAC2","Rq","0",-1,1);
 //n_k_noAC2=fk_noAC2->GetParameter(0);
 mean_k_noAC2=fk_noAC2->GetParameter(1);
 sig_k_noAC2=fk_noAC2->GetParameter(2);

 n_k_noAC2=hcoin_wo_bg_fom_noAC2->Integral(hcoin_wo_bg_fom_noAC2->FindBin(center_k-range_k),hcoin_wo_bg_fom_noAC2->FindBin(center_k+range_k));

 //----- Fitting as a whole function---------//

 fcoin_noAC2 =new TF1("fcoin_noAC2","gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin_noAC2->SetNpx(2000);
 fcoin_noAC2->SetTitle("Cointime w/o AC2 cut;Cointime [ns];Counts [1/56 ns]");
 fcoin_noAC2->SetParameters(n_pi_noAC2,mean_pi_noAC2,sig_pi_noAC2,n_k_noAC2,mean_k_noAC2,sig_k_noAC2,n_p_noAC2,mean_p_noAC2,sig_p_noAC2);
 //hcoin_wo_bg_fom_noAC2->Fit("fcoin_woAC2","Rq","0",min_coin_c,max_coin_c);
 //n_pi_noAC2=fcoin_noAC2->GetParameter(0);//Npi_nocut
 //mean_pi_noAC2=fcoin_noAC2->GetParameter(1);
 //sig_pi_noAC2=fcoin_noAC2->GetParameter(2);
 //n_k_noAC2=fcoin_noAC2->GetParameter(3);//Nk_nocut
 //mean_k_noAC2=fcoin_noAC2->GetParameter(4);
 //sig_k_noAC2=fcoin_noAC2->GetParameter(5);
 //n_p_noAC2=fcoin_noAC2->GetParameter(6);//Np_nocut
 //mean_p_noAC2=fcoin_noAC2->GetParameter(7);
 //sig_p_noAC2=fcoin_noAC2->GetParameter(8);
cout<<"n_pi_noAC2="<<n_pi_noAC2<<"n_k_noAC2="<<n_k_noAC2<<"n_p_noAC2="<<n_p_noAC2
<<"mean_pi_noAC2="<<mean_pi_noAC2<<"sig_pi_noAC2="<<sig_pi_noAC2<<"mean_k_noAC2="<<mean_k_noAC2<<"sig_k_noAC2="<<sig_k_noAC2<<"mean_p_noAC2="<<mean_p_noAC2<<"sig_p_noAC2="<<sig_p_noAC2<<endl;

 //------- Get Error Paramters ---//
 n_pi_err_noAC2=fpi_noAC2->GetParError(0); 
 n_k_err_noAC2=fk_noAC2->GetParError(3); 
 n_p_err_noAC2=fp_noAC2->GetParError(6); 
cout << "after nth thing..." << endl;

 hmm_bg_fom_noAC2->Scale(1./80.);
 hmm_wo_bg_fom_noAC2->Add(hmm_L_fom_noAC2,hmm_bg_fom_noAC2,1.0,-1.0);


// fmmbg_noAC2=new TF1("fmmbg_noAC2","pol3(0)",min_mm,max_mm);
// fmmbg_noAC2->SetNpx(2000);
 fL_noAC2=new TF1("fL_noAC2","gausn(0)",min_mm,max_mm);
 fL_noAC2->SetNpx(2000);
 fL_noAC2->SetParLimits(2,0.,0.01);
 fS_noAC2=new TF1("fS_noAC2","gausn(0)",min_mm,max_mm);
 fS_noAC2->SetNpx(2000);
 fS_noAC2->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
 
 fmm_noAC2=new TF1("fmm_noAC2","gausn(0)+gausn(3)",min_mm,max_mm);
 fmm_noAC2->SetNpx(2000);
 fmm_noAC2->SetTitle("Missing Mass w/o AC cut;Coin time [ns];Counts [1/56 ns]");



 //------- Fitting ----------//

//cout<<"fmmbg fit start"<<endl;
// hmm_bg_fom_noAC2->Fit("fmmbg_noAC2","Rq","",min_mm,max_mm);
// double d = fmmbg_noAC2->GetParameter(0);
// double c = fmmbg_noAC2->GetParameter(1);
// double b = fmmbg_noAC2->GetParameter(2);
// double a = fmmbg_noAC2->GetParameter(3);

//cout<<"fL fit start"<<endl;
 hmm_wo_bg_fom_noAC2->Fit("fL_noAC2","Rq","",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
// n_L_noAC2=fL_noAC2->GetParameter(0);
 mean_L_noAC2=fL_noAC2->GetParameter(1);
 sig_L_noAC2=fL_noAC2->GetParameter(2);
 mean_L_noAC2=def_mean_L;
 sig_L_noAC2=def_sig_L;
 n_L_noAC2=hmm_wo_bg_fom_noAC2->Integral(hmm_wo_bg_fom_noAC2->FindBin(center_L-range_L),hmm_wo_bg_fom_noAC2->FindBin(center_L+range_L));
 cout<<"n_L"<<n_L_noAC2<<endl;

//cout<<"fS fit start"<<endl;
 hmm_wo_bg_fom_noAC2->Fit("fS_noAC2","Rq","",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
// n_S_noAC2=fS_noAC2->GetParameter(0);
 mean_S_noAC2=fS_noAC2->GetParameter(1);
 sig_S_noAC2=fS_noAC2->GetParameter(2);
 mean_S_noAC2=def_mean_S;
 sig_S_noAC2=def_sig_S;
 n_S_noAC2=hmm_wo_bg_fom_noAC2->Integral(hmm_wo_bg_fom_noAC2->FindBin(center_S-range_S),hmm_wo_bg_fom_noAC2->FindBin(center_S+range_S));
 cout<<"n_S"<<n_S_noAC2<<endl;
//cout<<"mean_S"<<mean_S_noAC2<<endl;
//cout<<"sig_S"<<sig_S_noAC2<<endl;


 //--Fitting again as a total function--//
// fmm_noAC2->SetParameters(n_L_noAC2,mean_L_noAC2,sig_L_noAC2,n_S_noAC2,mean_S_noAC2,sig_S_noAC2,d,c,b,a);
//// hmm_L_fom_noAC2->Fit(Form("fmm[%d][%d][%d]",i,j,l),"Rq","",min_mm,max_mm);
// n_L_noAC2=fmm_noAC2->GetParameter(0);
// mean_L_noAC2=fmm_noAC2->GetParameter(1);
// sig_L_noAC2=fmm_noAC2->GetParameter(2);
// n_S_noAC2=fmm_noAC2->GetParameter(3);
// mean_S_noAC2=fmm_noAC2->GetParameter(4);
// sig_S_noAC2=fmm_noAC2->GetParameter(5);
// //------- Get Error Paramters ---//
// n_L_err_noAC2=fmm_noAC2->GetParError(0); 
// n_S_err_noAC2=fmm_noAC2->GetParError(3); 
//
//
//
// signal_noAC2=noise_noAC2=0.;
// d = fmm_noAC2->GetParError(6); 
// c = fmm_noAC2->GetParError(7); 
// b = fmm_noAC2->GetParError(8); 
// a = fmm_noAC2->GetParError(9); 




//-----No AC Cut-----//
 hcoin_bg_fom_noAC->Scale(40./160.);
 hcoin_wo_bg_fom_noAC->Add(hcoin_k_fom_noAC,hcoin_bg_fom_noAC,1.0,-1.0);
 fp_noAC=new TF1("fp_noAC","gausn(0)",min_coin_c,max_coin_c);
 fp_noAC->SetNpx(2000);
 fpi_noAC =new TF1("fpi_noAC","gausn(0)",min_coin_c,max_coin_c);
 fpi_noAC->SetNpx(2000);
 fk_noAC=new TF1("fk_noAC","gausn(0)",min_coin_c,max_coin_c);
 fk_noAC->SetNpx(2000);
//cout<<"fp fit start"<<endl;
 hcoin_wo_bg_fom_noAC->Fit("fp_noAC","Rq","0",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
// n_p_noAC=fp_noAC->GetParameter(0);
 mean_p_noAC=fp_noAC->GetParameter(1);
 sig_p_noAC=fp_noAC->GetParameter(2);
 n_p_noAC=hcoin_wo_bg_fom_noAC->Integral(hcoin_wo_bg_fom_noAC->FindBin(center_p-range_p),hcoin_wo_bg_fom_noAC->FindBin(center_p+range_p));
//cout<<"fpi fit start"<<endl;
 hcoin_wo_bg_fom_noAC->Fit("fpi_noAC","Rq","0",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
// n_pi_noAC=fpi_noAC->GetParameter(0);
 mean_pi_noAC=fpi_noAC->GetParameter(1);
 sig_pi_noAC=fpi_noAC->GetParameter(2);
 n_pi_noAC=hcoin_wo_bg_fom_noAC->Integral(hcoin_wo_bg_fom_noAC->FindBin(center_pi-range_pi),hcoin_wo_bg_fom_noAC->FindBin(center_pi+range_pi));
//cout<<"fk fit start"<<endl;
 //hcoin_wo_bg_fom_noAC->Fit("fk_noAC","Rq","0",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_wo_bg_fom_noAC->Fit("fk_noAC","Rq","0",-1,1);
 //n_k_noAC=fk_noAC->GetParameter(0);
 mean_k_noAC=fk_noAC->GetParameter(1);
// mean_k_noAC=0;
 sig_k_noAC=fk_noAC->GetParameter(2);

 n_k_noAC=hcoin_wo_bg_fom_noAC->Integral(hcoin_wo_bg_fom_noAC->FindBin(center_k-range_k),hcoin_wo_bg_fom_noAC->FindBin(center_k+range_k));

 //-----Fitting as a whole function---------//

 fcoin_noAC =new TF1("fcoin_noAC","gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin_noAC->SetNpx(2000);
 fcoin_noAC->SetTitle("Cointime w/o AC cut;Cointime [ns];Counts [1/56 ns]");
 fcoin_noAC->SetParameters(n_pi_noAC,mean_pi_noAC,sig_pi_noAC,n_k_noAC,mean_k_noAC,sig_k_noAC,n_p_noAC,mean_p_noAC,sig_p_noAC);
 //hcoin_wo_bg_fom_noAC->Fit("fcoin_woAC","Rq","0",min_coin_c,max_coin_c);
 //n_pi_noAC=fcoin_noAC->GetParameter(0);//Npi_nocut
 //mean_pi_noAC=fcoin_noAC->GetParameter(1);
 //sig_pi_noAC=fcoin_noAC->GetParameter(2);
 //n_k_noAC=fcoin_noAC->GetParameter(3);//Nk_nocut
 //mean_k_noAC=fcoin_noAC->GetParameter(4);
 //sig_k_noAC=fcoin_noAC->GetParameter(5);
 //n_p_noAC=fcoin_noAC->GetParameter(6);//Np_nocut
 //mean_p_noAC=fcoin_noAC->GetParameter(7);
 //sig_p_noAC=fcoin_noAC->GetParameter(8);
cout<<"n_pi_noAC="<<n_pi_noAC<<"n_k_noAC="<<n_k_noAC<<"n_p_noAC="<<n_p_noAC
<<"mean_pi_noAC="<<mean_pi_noAC<<"sig_pi_noAC="<<sig_pi_noAC<<"mean_k_noAC="<<mean_k_noAC<<"sig_k_noAC="<<sig_k_noAC<<"mean_p_noAC="<<mean_p_noAC<<"sig_p_noAC="<<sig_p_noAC<<endl;

 

//----------------------------------------------//
//--	Missing Mass  Start     ----------------//
//----------------------------------------------//
 hmm_bg_fom_noAC->Scale(1./80.);
 hmm_pibg_fom_noAC->Scale(0.7/80.);
 hmm_wo_bg_fom_noAC->Add(hmm_L_fom_noAC,hmm_bg_fom_noAC,1.0,-1.0);
 hmm_pi_wobg_fom_noAC->Add(hmm_pi_fom_noAC,hmm_pibg_fom_noAC,1.0,-1.0);


// fmmbg_noAC=new TF1("fmmbg_noAC","gausn(0)+gausn(3)",min_mm,max_mm);
 //fmmbg_noAC=new TF1("fmmbg_noAC",F_Voigt,min_mm,max_mm,4);
 fmmbg_noAC=new TF1("fmmbg_noAC","pol4",-0.05,0.15);
// fmmbg_noAC->SetParameters(5,0.05,0.05,0.01);
// fmmbg_noAC->SetParLimits(0,0.,100000.);//positive
// fmmbg_noAC->SetParLimits(3,0.,100.);//positive
 fmmbg_noAC->SetNpx(2000);
 //fmmbg_noAC->SetParameters(100,0.05,0.03,10,0.05,0.03);//test.list
// fmmbg_noAC->SetParameters(300,0.05,1.2,30,0.1,0.02);//small.list
// fmmbg_noAC->SetParameter(1,0.05);
// fmmbg_noAC->SetParameter(2,0.03);
 fL_noAC=new TF1("fL_noAC","gausn(0)",min_mm,max_mm);
 fL_noAC->SetNpx(2000);
 fL_noAC->SetParLimits(2,0.,0.01);
 fS_noAC=new TF1("fS_noAC","gausn(0)",min_mm,max_mm);
 fS_noAC->SetNpx(2000);
 fS_noAC->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
 
 fmm_noAC=new TF1("fmm_noAC","gausn(0)+gausn(3)+pol4(6)",-0.05,0.15);
// fmm_noAC=new TF1("fmm_noAC",F_mmnoAC,-0.1,0.15,10);
 fmm_noAC->SetNpx(2000);
 fmm_noAC->SetTitle("Missing Mass w/o AC cut;Coin time [ns];Counts [1/56 ns]");
 fmm_noAC->SetParLimits(0,0.,1000000.);//positive
 fmm_noAC->SetParLimits(3,0.,1000000.);//positive
// fmm_noAC->SetParLimits(3,0.,100.);//positive//Voigt
// fmm_noAC->SetParLimits(4,0.,100000.);//positive//Voigt
// fmm_noAC->SetParLimits(7,0.,100000.);//positive//Voigt
// fmm_noAC->SetParameter(1,def_mean_L);
// fmm_noAC->SetParameter(4,def_mean_S);



 //------- Fitting ----------//

cout<<"fmmbg fit start"<<endl;
// hmm_bg_fom_noAC->Fit("fmmbg_noAC","Rq","",min_mm,max_mm);
// hmm_pi_wobg_fom_noAC->Fit("fmmbg_noAC","Rq","",-0.05,0.15);
//double fmmbga = fmmbg_noAC->GetParameter(0);
//double fmmbgb = fmmbg_noAC->GetParameter(1);
//double fmmbgc = fmmbg_noAC->GetParameter(2);
//double fmmbgd = fmmbg_noAC->GetParameter(3);
//cout<<"0:1:2:3:4:5="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<endl;//"::"<<e<<"::"<<f<<endl;
// fmm_noAC->FixParameter(1,b);//mean
// fmm_noAC->FixParameter(2,c);//sigma
// fmm_noAC->FixParameter(3,d);//lg
//double e = fmmbg_noAC->GetParameter(4);
//double f = fmmbg_noAC->GetParameter(5);
// fmm_noAC->SetParameter(0,a);
// fmm_noAC->SetParameter(1,b);
// fmm_noAC->SetParameter(2,c);
// fmm_noAC->SetParameter(3,d);
// fmm_noAC->SetParameter(4,e);
// fmm_noAC->SetParameter(5,f);
// fmm_noAC->SetParameter(6,500);
 fmm_noAC->SetParameter(1,def_mean_L);
 fmm_noAC->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
 fmm_noAC->SetParameter(2,def_sig_L);
 fmm_noAC->SetParLimits(2,0.,2*def_sig_L);
// fmm_noAC->SetParameters(9,100);
 fmm_noAC->SetParameter(4,def_mean_S);
 fmm_noAC->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
 fmm_noAC->SetParameter(5,def_sig_S);
 fmm_noAC->SetParLimits(5,0.,2*def_sig_S);
 hmm_wo_bg_fom_noAC->Fit("fmm_noAC","","",-0.05,0.15);
double fmm_noACpar0 = fmm_noAC->GetParameter(0);cout<<"fmm_noAC[0]="<<fmm_noACpar0<<endl;//area(L)
double fmm_noACpar1 = fmm_noAC->GetParameter(1);cout<<"fmm_noAC[1]="<<fmm_noACpar1<<endl;//mean(L)
double fmm_noACpar2 = fmm_noAC->GetParameter(2);cout<<"fmm_noAC[2]="<<fmm_noACpar2<<endl;//sigma(L)
double fmm_noACpar3 = fmm_noAC->GetParameter(3);cout<<"fmm_noAC[3]="<<fmm_noACpar3<<endl;//area(S)
double fmm_noACpar4 = fmm_noAC->GetParameter(4);cout<<"fmm_noAC[4]="<<fmm_noACpar4<<endl;//mean(S)
double fmm_noACpar5 = fmm_noAC->GetParameter(5);cout<<"fmm_noAC[5]="<<fmm_noACpar5<<endl;//sigma(S)
double fmm_noACpar6 = fmm_noAC->GetParameter(6);cout<<"fmm_noAC[6]="<<fmm_noACpar6<<endl;//poly_const
double fmm_noACpar7 = fmm_noAC->GetParameter(7);cout<<"fmm_noAC[7]="<<fmm_noACpar7<<endl;//poly_x
double fmm_noACpar8 = fmm_noAC->GetParameter(8);cout<<"fmm_noAC[8]="<<fmm_noACpar8<<endl;//poly_x^2
double fmm_noACpar9 = fmm_noAC->GetParameter(9);cout<<"fmm_noAC[9]="<<fmm_noACpar9<<endl;//poly_x^3
double fmm_noACpar10 = fmm_noAC->GetParameter(10);cout<<"fmm_noAC[10]="<<fmm_noACpar10<<endl;//poly_x^4
 fmmbg_noAC->SetParameters(fmm_noACpar6,fmm_noACpar7,fmm_noACpar8,fmm_noACpar9,fmm_noACpar10);
//cout<<"0:1:2:3:4:5(as a total func)="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<"::"<<e<<"::"<<f<<endl;
//
//cout<<"fL fit start"<<endl;
// hmm_wo_bg_fom_noAC->Fit("fL_noAC","Rq","",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
// n_L_noAC=fL_noAC->GetParameter(0);
// mean_L_noAC=fL_noAC->GetParameter(1);
// sig_L_noAC=fL_noAC->GetParameter(2);
 mean_L_noAC=def_mean_L;
 mean_S_noAC=def_mean_S;
 sig_L_noAC=def_sig_L;
 sig_S_noAC=def_sig_S;
 n_L_noAC=hmm_wo_bg_fom_noAC->Integral(hmm_wo_bg_fom_noAC->FindBin(center_L-range_L),hmm_wo_bg_fom_noAC->FindBin(center_L+range_L));
cout<<"before(L):: "<<n_L_noAC<<endl;
 double integralL=fmmbg_noAC->Integral(center_L-range_L,center_L+range_L);
 if(integralL>0)n_L_noAC-=integralL;
cout<<"after(L):: "<<n_L_noAC<<endl;
//n_L_noAC-=(pow(mean_L_noAC+2*sig_L_noAC,5)-pow(mean_L_noAC-2*sig_L_noAC,5))*a/5;
//n_L_noAC-=(pow(mean_L_noAC+2*sig_L_noAC,4)-pow(mean_L_noAC-2*sig_L_noAC,4))*b/4;
//n_L_noAC-=(pow(mean_L_noAC+2*sig_L_noAC,3)-pow(mean_L_noAC-2*sig_L_noAC,3))*c/3;
//n_L_noAC-=(pow(mean_L_noAC+2*sig_L_noAC,2)-pow(mean_L_noAC-2*sig_L_noAC,2))*d/2;
//n_L_noAC-=(pow(mean_L_noAC+2*sig_L_noAC,1)-pow(mean_L_noAC-2*sig_L_noAC,1))*e;
//
 n_S_noAC=hmm_wo_bg_fom_noAC->Integral(hmm_wo_bg_fom_noAC->FindBin(center_S-range_S),hmm_wo_bg_fom_noAC->FindBin(center_S+range_S));
cout<<"before(S):: "<<n_S_noAC<<endl;
 double integralS=fmmbg_noAC->Integral(center_S-range_S,center_S+range_S);
 if(integralS>0)n_S_noAC-=integralS;
cout<<"after(S):: "<<n_S_noAC<<endl;
//n_S_noAC-=(pow(mean_S_noAC+2*sig_S_noAC,5)-pow(mean_S_noAC-2*sig_S_noAC,5))*a/5;
//n_S_noAC-=(pow(mean_S_noAC+2*sig_S_noAC,4)-pow(mean_S_noAC-2*sig_S_noAC,4))*b/4;
//n_S_noAC-=(pow(mean_S_noAC+2*sig_S_noAC,3)-pow(mean_S_noAC-2*sig_S_noAC,3))*c/3;
//n_S_noAC-=(pow(mean_S_noAC+2*sig_S_noAC,2)-pow(mean_S_noAC-2*sig_S_noAC,2))*d/2;
//n_S_noAC-=(pow(mean_S_noAC+2*sig_S_noAC,1)-pow(mean_S_noAC-2*sig_S_noAC,1))*e;

 cout<<"n_L"<<n_L_noAC<<endl;
cout<<"mean_L"<<mean_L_noAC<<endl;
cout<<"sig_L"<<sig_L_noAC<<endl;

////cout<<"fS fit start"<<endl;
// hmm_wo_bg_fom_noAC->Fit("fS_noAC","Rq","",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
//// n_S_noAC=fS_noAC->GetParameter(0);
// mean_S_noAC=fS_noAC->GetParameter(1);
// mean_S_noAC=def_mean_S;
// sig_S_noAC=def_sig_S;
// sig_S_noAC=fS_noAC->GetParameter(2);
// n_S_noAC=hmm_wo_bg_fom_noAC->Integral(hmm_wo_bg_fom_noAC->FindBin(mean_S_noAC-2*sig_S_noAC),hmm_wo_bg_fom_noAC->FindBin(mean_S_noAC+2*sig_S_noAC));
 cout<<"n_S"<<n_S_noAC<<endl;
cout<<"mean_S"<<mean_S_noAC<<endl;
cout<<"sig_S"<<sig_S_noAC<<endl;


//-----No AC Cut-----//
for(int i=0;i<nth;i++){
 hcoin_bg_fom_noZ[i]->Scale(40./160.);
 hcoin_wo_bg_fom_noZ[i]->Add(hcoin_k_fom_noZ[i],hcoin_bg_fom_noZ[i],1.0,-1.0);
 fp_noZ[i]=new TF1("fp_noZ[i]","gausn(0)",min_coin_c,max_coin_c);
 fp_noZ[i]->SetNpx(2000);
 fpi_noZ[i] =new TF1("fpi_noZ[i]","gausn(0)",min_coin_c,max_coin_c);
 fpi_noZ[i]->SetNpx(2000);
 fk_noZ[i]=new TF1("fk_noZ[i]","gausn(0)",min_coin_c,max_coin_c);
 fk_noZ[i]->SetNpx(2000);
//cout<<"fp fit start"<<endl;
 hcoin_wo_bg_fom_noZ[i]->Fit("fp_noZ[i]","Rq","0",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
// n_p_noZ[i]=fp_noZ[i]->GetParameter(0);
 mean_p_noZ[i]=fp_noZ[i]->GetParameter(1);
 sig_p_noZ[i]=fp_noZ[i]->GetParameter(2);
 n_p_noZ[i]=hcoin_wo_bg_fom_noZ[i]->Integral(hcoin_wo_bg_fom_noZ[i]->FindBin(center_p-range_p),hcoin_wo_bg_fom_noZ[i]->FindBin(center_p+range_p));
//cout<<"fpi fit start"<<endl;
 hcoin_wo_bg_fom_noZ[i]->Fit("fpi_noZ[i]","Rq","0",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
// n_pi_noZ[i]=fpi_noZ[i]->GetParameter(0);
 mean_pi_noZ[i]=fpi_noZ[i]->GetParameter(1);
 sig_pi_noZ[i]=fpi_noZ[i]->GetParameter(2);
 n_pi_noZ[i]=hcoin_wo_bg_fom_noZ[i]->Integral(hcoin_wo_bg_fom_noZ[i]->FindBin(center_pi-range_pi),hcoin_wo_bg_fom_noZ[i]->FindBin(center_pi+range_pi));
//cout<<"fk fit start"<<endl;
 //hcoin_wo_bg_fom_noZ[i]->Fit("fk_noZ[i]","Rq","0",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_wo_bg_fom_noZ[i]->Fit("fk_noZ[i]","Rq","0",-1,1);
 //n_k_noZ[i]=fk_noZ[i]->GetParameter(0);
 mean_k_noZ[i]=fk_noZ[i]->GetParameter(1);
// mean_k_noZ[i]=0;
 sig_k_noZ[i]=fk_noZ[i]->GetParameter(2);

 n_k_noZ[i]=hcoin_wo_bg_fom_noZ[i]->Integral(hcoin_wo_bg_fom_noZ[i]->FindBin(center_k-range_k),hcoin_wo_bg_fom_noZ[i]->FindBin(center_k+range_k));

 //-----Fitting as a whole function---------//

 fcoin_noZ[i] =new TF1("fcoin_noZ[i]","gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin_noZ[i]->SetNpx(2000);
 fcoin_noZ[i]->SetTitle("Cointime w/o AC cut;Cointime [ns];Counts [1/56 ns]");
 fcoin_noZ[i]->SetParameters(n_pi_noZ[i],mean_pi_noZ[i],sig_pi_noZ[i],n_k_noZ[i],mean_k_noZ[i],sig_k_noZ[i],n_p_noZ[i],mean_p_noZ[i],sig_p_noZ[i]);
 //hcoin_wo_bg_fom_noZ[i]->Fit("fcoin_woAC","Rq","0",min_coin_c,max_coin_c);
 //n_pi_noZ[i]=fcoin_noZ[i]->GetParameter(0);//Npi_nocut
 //mean_pi_noZ[i]=fcoin_noZ[i]->GetParameter(1);
 //sig_pi_noZ[i]=fcoin_noZ[i]->GetParameter(2);
 //n_k_noZ[i]=fcoin_noZ[i]->GetParameter(3);//Nk_nocut
 //mean_k_noZ[i]=fcoin_noZ[i]->GetParameter(4);
 //sig_k_noZ[i]=fcoin_noZ[i]->GetParameter(5);
 //n_p_noZ[i]=fcoin_noZ[i]->GetParameter(6);//Np_nocut
 //mean_p_noZ[i]=fcoin_noZ[i]->GetParameter(7);
 //sig_p_noZ[i]=fcoin_noZ[i]->GetParameter(8);
cout<<"n_pi_noZ[i]="<<n_pi_noZ[i]<<"n_k_noZ[i]="<<n_k_noZ[i]<<"n_p_noZ[i]="<<n_p_noZ[i]
<<"mean_pi_noZ[i]="<<mean_pi_noZ[i]<<"sig_pi_noZ[i]="<<sig_pi_noZ[i]<<"mean_k_noZ[i]="<<mean_k_noZ[i]<<"sig_k_noZ[i]="<<sig_k_noZ[i]<<"mean_p_noZ[i]="<<mean_p_noZ[i]<<"sig_p_noZ[i]="<<sig_p_noZ[i]<<endl;

 

//----------------------------------------------//
//--	Missing Mass  Start     ----------------//
//----------------------------------------------//
 hmm_bg_fom_noZ[i]->Scale(1./80.);
 hmm_pibg_fom_noZ[i]->Scale(0.7/80.);
 hmm_wo_bg_fom_noZ[i]->Add(hmm_L_fom_noZ[i],hmm_bg_fom_noZ[i],1.0,-1.0);
 hmm_pi_wobg_fom_noZ[i]->Add(hmm_pi_fom_noZ[i],hmm_pibg_fom_noZ[i],1.0,-1.0);


// fmmbg_noZ[i]=new TF1("fmmbg_noZ[i]","gausn(0)+gausn(3)",min_mm,max_mm);
 //fmmbg_noZ[i]=new TF1("fmmbg_noZ[i]",F_Voigt,min_mm,max_mm,4);
 fmmbg_noZ[i]=new TF1("fmmbg_noZ[i]","pol4",-0.05,0.15);
// fmmbg_noZ[i]->SetParameters(5,0.05,0.05,0.01);
// fmmbg_noZ[i]->SetParLimits(0,0.,100000.);//positive
// fmmbg_noZ[i]->SetParLimits(3,0.,100.);//positive
 fmmbg_noZ[i]->SetNpx(2000);
 //fmmbg_noZ[i]->SetParameters(100,0.05,0.03,10,0.05,0.03);//test.list
// fmmbg_noZ[i]->SetParameters(300,0.05,1.2,30,0.1,0.02);//small.list
// fmmbg_noZ[i]->SetParameter(1,0.05);
// fmmbg_noZ[i]->SetParameter(2,0.03);
 fL_noZ[i]=new TF1("fL_noZ[i]","gausn(0)",min_mm,max_mm);
 fL_noZ[i]->SetNpx(2000);
 fL_noZ[i]->SetParLimits(2,0.,0.01);
 fS_noZ[i]=new TF1("fS_noZ[i]","gausn(0)",min_mm,max_mm);
 fS_noZ[i]->SetNpx(2000);
 fS_noZ[i]->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
 
 fmm_noZ[i]=new TF1("fmm_noZ[i]","gausn(0)+gausn(3)+pol4(6)",-0.05,0.15);
// fmm_noZ[i]=new TF1("fmm_noZ[i]",F_mmnoZ[i],-0.1,0.15,10);
 fmm_noZ[i]->SetNpx(2000);
 fmm_noZ[i]->SetTitle("Missing Mass w/o AC cut;Coin time [ns];Counts [1/56 ns]");
 fmm_noZ[i]->SetParLimits(0,0.,1000000.);//positive
 fmm_noZ[i]->SetParLimits(3,0.,1000000.);//positive
// fmm_noZ[i]->SetParLimits(3,0.,100.);//positive//Voigt
// fmm_noZ[i]->SetParLimits(4,0.,100000.);//positive//Voigt
// fmm_noZ[i]->SetParLimits(7,0.,100000.);//positive//Voigt
// fmm_noZ[i]->SetParameter(1,def_mean_L);
// fmm_noZ[i]->SetParameter(4,def_mean_S);



 //------- Fitting ----------//

cout<<"fmmbg fit start"<<endl;
// hmm_bg_fom_noZ[i]->Fit("fmmbg_noZ[i]","Rq","",min_mm,max_mm);
// hmm_pi_wobg_fom_noZ[i]->Fit("fmmbg_noZ[i]","Rq","",-0.05,0.15);
//double fmmbga = fmmbg_noZ[i]->GetParameter(0);
//double fmmbgb = fmmbg_noZ[i]->GetParameter(1);
//double fmmbgc = fmmbg_noZ[i]->GetParameter(2);
//double fmmbgd = fmmbg_noZ[i]->GetParameter(3);
//cout<<"0:1:2:3:4:5="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<endl;//"::"<<e<<"::"<<f<<endl;
// fmm_noZ[i]->FixParameter(1,b);//mean
// fmm_noZ[i]->FixParameter(2,c);//sigma
// fmm_noZ[i]->FixParameter(3,d);//lg
//double e = fmmbg_noZ[i]->GetParameter(4);
//double f = fmmbg_noZ[i]->GetParameter(5);
// fmm_noZ[i]->SetParameter(0,a);
// fmm_noZ[i]->SetParameter(1,b);
// fmm_noZ[i]->SetParameter(2,c);
// fmm_noZ[i]->SetParameter(3,d);
// fmm_noZ[i]->SetParameter(4,e);
// fmm_noZ[i]->SetParameter(5,f);
// fmm_noZ[i]->SetParameter(6,500);
 fmm_noZ[i]->SetParameter(1,def_mean_L);
 fmm_noZ[i]->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
 fmm_noZ[i]->SetParameter(2,def_sig_L);
 fmm_noZ[i]->SetParLimits(2,0.,2*def_sig_L);
// fmm_noZ[i]->SetParameters(9,100);
 fmm_noZ[i]->SetParameter(4,def_mean_S);
 fmm_noZ[i]->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
 fmm_noZ[i]->SetParameter(5,def_sig_S);
 fmm_noZ[i]->SetParLimits(5,0.,2*def_sig_S);
 hmm_wo_bg_fom_noZ[i]->Fit("fmm_noZ[i]","","",-0.05,0.15);
double fmm_noZpar0 = fmm_noZ[i]->GetParameter(0);cout<<"fmm_noZ[0]="<<fmm_noZpar0<<endl;//area(L)
double fmm_noZpar1 = fmm_noZ[i]->GetParameter(1);cout<<"fmm_noZ[1]="<<fmm_noZpar1<<endl;//mean(L)
double fmm_noZpar2 = fmm_noZ[i]->GetParameter(2);cout<<"fmm_noZ[2]="<<fmm_noZpar2<<endl;//sigma(L)
double fmm_noZpar3 = fmm_noZ[i]->GetParameter(3);cout<<"fmm_noZ[3]="<<fmm_noZpar3<<endl;//area(S)
double fmm_noZpar4 = fmm_noZ[i]->GetParameter(4);cout<<"fmm_noZ[4]="<<fmm_noZpar4<<endl;//mean(S)
double fmm_noZpar5 = fmm_noZ[i]->GetParameter(5);cout<<"fmm_noZ[5]="<<fmm_noZpar5<<endl;//sigma(S)
double fmm_noZpar6 = fmm_noZ[i]->GetParameter(6);cout<<"fmm_noZ[6]="<<fmm_noZpar6<<endl;//poly_const
double fmm_noZpar7 = fmm_noZ[i]->GetParameter(7);cout<<"fmm_noZ[7]="<<fmm_noZpar7<<endl;//poly_x
double fmm_noZpar8 = fmm_noZ[i]->GetParameter(8);cout<<"fmm_noZ[8]="<<fmm_noZpar8<<endl;//poly_x^2
double fmm_noZpar9 = fmm_noZ[i]->GetParameter(9);cout<<"fmm_noZ[9]="<<fmm_noZpar9<<endl;//poly_x^3
double fmm_noZpar10 = fmm_noZ[i]->GetParameter(10);cout<<"fmm_noZ[10]="<<fmm_noZpar10<<endl;//poly_x^4
 fmmbg_noZ[i]->SetParameters(fmm_noZpar6,fmm_noZpar7,fmm_noZpar8,fmm_noZpar9,fmm_noZpar10);
//cout<<"0:1:2:3:4:5(as a total func)="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<"::"<<e<<"::"<<f<<endl;
//
//cout<<"fL fit start"<<endl;
// hmm_wo_bg_fom_noZ[i]->Fit("fL_noZ[i]","Rq","",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
// n_L_noZ[i]=fL_noZ[i]->GetParameter(0);
// mean_L_noZ[i]=fL_noZ[i]->GetParameter(1);
// sig_L_noZ[i]=fL_noZ[i]->GetParameter(2);
 mean_L_noZ[i]=def_mean_L;
 mean_S_noZ[i]=def_mean_S;
 sig_L_noZ[i]=def_sig_L;
 sig_S_noZ[i]=def_sig_S;
 n_L_noZ[i]=hmm_wo_bg_fom_noZ[i]->Integral(hmm_wo_bg_fom_noZ[i]->FindBin(center_L-range_L),hmm_wo_bg_fom_noZ[i]->FindBin(center_L+range_L));
cout<<"before(L):: "<<n_L_noZ[i]<<endl;
 double integralL=fmmbg_noZ[i]->Integral(center_L-range_L,center_L+range_L);
 if(integralL>0)n_L_noZ[i]-=integralL;
cout<<"after(L):: "<<n_L_noZ[i]<<endl;
//n_L_noZ[i]-=(pow(mean_L_noZ[i]+2*sig_L_noZ[i],5)-pow(mean_L_noZ[i]-2*sig_L_noZ[i],5))*a/5;
//n_L_noZ[i]-=(pow(mean_L_noZ[i]+2*sig_L_noZ[i],4)-pow(mean_L_noZ[i]-2*sig_L_noZ[i],4))*b/4;
//n_L_noZ[i]-=(pow(mean_L_noZ[i]+2*sig_L_noZ[i],3)-pow(mean_L_noZ[i]-2*sig_L_noZ[i],3))*c/3;
//n_L_noZ[i]-=(pow(mean_L_noZ[i]+2*sig_L_noZ[i],2)-pow(mean_L_noZ[i]-2*sig_L_noZ[i],2))*d/2;
//n_L_noZ[i]-=(pow(mean_L_noZ[i]+2*sig_L_noZ[i],1)-pow(mean_L_noZ[i]-2*sig_L_noZ[i],1))*e;
//
 n_S_noZ[i]=hmm_wo_bg_fom_noZ[i]->Integral(hmm_wo_bg_fom_noZ[i]->FindBin(center_S-range_S),hmm_wo_bg_fom_noZ[i]->FindBin(center_S+range_S));
cout<<"before(S):: "<<n_S_noZ[i]<<endl;
 double integralS=fmmbg_noZ[i]->Integral(center_S-range_S,center_S+range_S);
 if(integralS>0)n_S_noZ[i]-=integralS;
cout<<"after(S):: "<<n_S_noZ[i]<<endl;
//n_S_noZ[i]-=(pow(mean_S_noZ[i]+2*sig_S_noZ[i],5)-pow(mean_S_noZ[i]-2*sig_S_noZ[i],5))*a/5;
//n_S_noZ[i]-=(pow(mean_S_noZ[i]+2*sig_S_noZ[i],4)-pow(mean_S_noZ[i]-2*sig_S_noZ[i],4))*b/4;
//n_S_noZ[i]-=(pow(mean_S_noZ[i]+2*sig_S_noZ[i],3)-pow(mean_S_noZ[i]-2*sig_S_noZ[i],3))*c/3;
//n_S_noZ[i]-=(pow(mean_S_noZ[i]+2*sig_S_noZ[i],2)-pow(mean_S_noZ[i]-2*sig_S_noZ[i],2))*d/2;
//n_S_noZ[i]-=(pow(mean_S_noZ[i]+2*sig_S_noZ[i],1)-pow(mean_S_noZ[i]-2*sig_S_noZ[i],1))*e;

 cout<<"n_L"<<n_L_noZ[i]<<endl;
cout<<"mean_L"<<mean_L_noZ[i]<<endl;
cout<<"sig_L"<<sig_L_noZ[i]<<endl;

////cout<<"fS fit start"<<endl;
// hmm_wo_bg_fom_noZ[i]->Fit("fS_noZ[i]","Rq","",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
//// n_S_noZ[i]=fS_noZ[i]->GetParameter(0);
// mean_S_noZ[i]=fS_noZ[i]->GetParameter(1);
// mean_S_noZ[i]=def_mean_S;
// sig_S_noZ[i]=def_sig_S;
// sig_S_noZ[i]=fS_noZ[i]->GetParameter(2);
// n_S_noZ[i]=hmm_wo_bg_fom_noZ[i]->Integral(hmm_wo_bg_fom_noZ[i]->FindBin(mean_S_noZ[i]-2*sig_S_noZ[i]),hmm_wo_bg_fom_noZ[i]->FindBin(mean_S_noZ[i]+2*sig_S_noZ[i]));
 cout<<"n_S"<<n_S_noZ[i]<<endl;
cout<<"mean_S"<<mean_S_noZ[i]<<endl;
cout<<"sig_S"<<sig_S_noZ[i]<<endl;
}//for i
//----------------------------------------------//
//--	Missing Mass  End     ------------------//
//----------------------------------------------//
 //----------------------------------------------------------------//
 //----------------------------------------------------------------//
 //----------------------------------------------------------------//

	for(int i=0;i<nth;i++){
		for(int j=0;j<nth;j++){
//			for(int l=0;l<nth;l++){

			int l=0;
			
//-----Background subtraction-----//
//------------cointime------------//
//cout<<"BG subtraction cointime"<<endl;
 hcoin_bg_fom[i][j][l]->Scale(40./160.);
 hcoin_wo_bg_fom[i][j][l]->Add(hcoin_k_fom[i][j][l],hcoin_bg_fom[i][j][l],1.0,-1.0);

 fp[i][j][l]=new TF1(Form("fp[%d][%d][%d]",i,j,l),"gausn(0)",min_coin_c,max_coin_c);
 fp[i][j][l]->SetNpx(2000);
 fpi[i][j][l] =new TF1(Form("fpi[%d][%d][%d]",i,j,l),"gausn(0)",min_coin_c,max_coin_c);
 fpi[i][j][l]->SetNpx(2000);
 fk[i][j][l]=new TF1(Form("fk[%d][%d][%d]",i,j,l),"gausn(0)",min_coin_c,max_coin_c);
 fk[i][j][l]->SetNpx(2000);
//cout<<"fp fit start"<<endl;
 hcoin_wo_bg_fom[i][j][l]->Fit(Form("fp[%d][%d][%d]",i,j,l),"Rq","0",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 //n_p[i][j][l]=fp[i][j][l]->GetParameter(0);
 mean_p[i][j][l]=fp[i][j][l]->GetParameter(1);
 sig_p[i][j][l]=fp[i][j][l]->GetParameter(2);
 n_p[i][j][l]=hcoin_wo_bg_fom[i][j][l]->Integral(hcoin_wo_bg_fom[i][j][l]->FindBin(center_p-range_p),hcoin_wo_bg_fom[i][j][l]->FindBin(center_p+range_p));
//cout<<"fpi fit start"<<endl;
 hcoin_wo_bg_fom[i][j][l]->Fit(Form("fpi[%d][%d][%d]",i,j,l),"Rq","0",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 //n_pi[i][j][l]=fpi[i][j][l]->GetParameter(0);
 mean_pi[i][j][l]=fpi[i][j][l]->GetParameter(1);
cout << mean_pi[i][j][l] << endl;
 sig_pi[i][j][l]=fpi[i][j][l]->GetParameter(2);
n_pi[i][j][l]=hcoin_wo_bg_fom[i][j][l]->Integral(hcoin_wo_bg_fom[i][j][l]->FindBin(center_pi-range_pi),hcoin_wo_bg_fom[i][j][l]->FindBin(center_pi+range_pi));
//cout<<"fk fit start"<<endl;
// hcoin_wo_bg_fom[i][j][l]->Fit(Form("fk[%d][%d][%d]",i,j,l),"Rq","0",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_wo_bg_fom[i][j][l]->Fit(Form("fk[%d][%d][%d]",i,j,l),"Rq","0",-1,1);
// n_k[i][j][l]=fk[i][j][l]->GetParameter(0);
 mean_k[i][j][l]=fk[i][j][l]->GetParameter(1);
 sig_k[i][j][l]=fk[i][j][l]->GetParameter(2);
 n_k[i][j][l]=hcoin_wo_bg_fom[i][j][l]->Integral(hcoin_wo_bg_fom[i][j][l]->FindBin(center_k-range_k),hcoin_wo_bg_fom[i][j][l]->FindBin(center_k+range_k));

 // n_k[i][th1][1]=hcoin_k_fom[i][th1]->Integral(hcoin_k_fom[i][th1]->FindBin(-3*sig_k[i][th1][1]+mean_k[i][th1][1])
 //		     ,hcoin_k_fom[i][th1]->FindBin(+3*sig_k[i][th1][1]+mean_k[i][th1][1]));
 fcoin[i][j][l] =new TF1(Form("fcoin[%d][%d][%d]",i,j,l),"gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin[i][j][l]->SetNpx(2000);
 fcoin[i][j][l]->SetTitle(Form("Coin-Time w Z cut ( z_sum <%lf ch && z_diff < %lf ch);Coin time [ns];Counts [1/56 ns]",zver[i],zver_diff[j]));
 fcoin[i][j][l]->SetParameters(n_pi[i][j][l],mean_pi[i][j][l],sig_pi[i][j][l],n_k[i][j][l],mean_k[i][j][l],sig_k[i][j][l],n_p[i][j][l],mean_p[i][j][l],sig_p[i][j][l]);
 //hcoin_wo_bg_fom[i][j][l]->Fit(Form("fcoin[%d][%d][%d]",i,j,l),"Rq","0",min_coin_c,max_coin_c);
 //n_pi[i][j][l]=fcoin[i][j][l]->GetParameter(0);//Npi_nocut
 //mean_pi[i][j][l]=fcoin[i][j][l]->GetParameter(1);
 //sig_pi[i][j][l]=fcoin[i][j][l]->GetParameter(2);
 //n_k[i][j][l]=fcoin[i][j][l]->GetParameter(3);//Nk_nocut
 //mean_k[i][j][l]=fcoin[i][j][l]->GetParameter(4);
 //sig_k[i][j][l]=fcoin[i][j][l]->GetParameter(5);
 //n_p[i][j][l]=fcoin[i][j][l]->GetParameter(6);//Np_nocut
 //mean_p[i][j][l]=fcoin[i][j][l]->GetParameter(7);
 //sig_p[i][j][l]=fcoin[i][j][l]->GetParameter(8);

cout<<"n_pi["<<i<<"]["<<j<<"]["<<l<<"]="<<n_pi[i][j][l]<<endl;
cout<<"n_k["<<i<<"]["<<j<<"]["<<l<<"]="<<n_k[i][j][l]<<endl;
cout<<"n_p["<<i<<"]["<<j<<"]["<<l<<"]="<<n_p[i][j][l]<<endl;

//----------------------------------------------//
//--	Missing Mass  Start     ----------------//
//----------------------------------------------//
 hmm_bg_fom[i][j][l]->Scale(1./80.);
 hmm_pibg_fom[i][j][l]->Scale(0.7/80.);
 hmm_wo_bg_fom[i][j][l]->Add(hmm_L_fom[i][j][l],hmm_bg_fom[i][j][l],1.0,-1.0);
 hmm_pi_wobg_fom[i][j][l]->Add(hmm_pi_fom[i][j][l],hmm_pibg_fom[i][j][l],1.0,-1.0);


// fmmbg[i][j][l]=new TF1("fmmbg[i][j][l]","gausn(0)+gausn(3)",min_mm,max_mm);
 //fmmbg[i][j][l]=new TF1("fmmbg[i][j][l]",F_Voigt,min_mm,max_mm,4);
 fmmbg[i][j][l]=new TF1("fmmbg[i][j][l]","pol4",-0.05,0.15);
// fmmbg[i][j][l]->SetParameters(5,0.05,0.05,0.01);
// fmmbg[i][j][l]->SetParLimits(0,0.,100000.);//positive
// fmmbg[i][j][l]->SetParLimits(3,0.,100.);//positive
 fmmbg[i][j][l]->SetNpx(2000);
 //fmmbg[i][j][l]->SetParameters(100,0.05,0.03,10,0.05,0.03);//test.list
// fmmbg[i][j][l]->SetParameters(300,0.05,1.2,30,0.1,0.02);//small.list
// fmmbg[i][j][l]->SetParameter(1,0.05);
// fmmbg[i][j][l]->SetParameter(2,0.03);
 fL[i][j][l]=new TF1("fL[i][j][l]","gausn(0)",min_mm,max_mm);
 fL[i][j][l]->SetNpx(2000);
 fL[i][j][l]->SetParLimits(2,0.,0.01);
 fS[i][j][l]=new TF1("fS[i][j][l]","gausn(0)",min_mm,max_mm);
 fS[i][j][l]->SetNpx(2000);
 fS[i][j][l]->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
 
 fmm[i][j][l]=new TF1("fmm[i][j][l]","gausn(0)+gausn(3)+pol4(6)",-0.05,0.15);
// fmm[i][j][l]=new TF1("fmm[i][j][l]",F_mmnoAC,-0.1,0.15,10);
 fmm[i][j][l]->SetNpx(2000);
 fmm[i][j][l]->SetTitle("Missing Mass w/o AC cut;Coin time [ns];Counts [1/56 ns]");
 fmm[i][j][l]->SetParLimits(0,0.,1000000.);//positive
 fmm[i][j][l]->SetParLimits(3,0.,1000000.);//positive
// fmm[i][j][l]->SetParLimits(3,0.,100.);//positive//Voigt
// fmm[i][j][l]->SetParLimits(4,0.,100000.);//positive//Voigt
// fmm[i][j][l]->SetParLimits(7,0.,100000.);//positive//Voigt
// fmm[i][j][l]->SetParameter(1,def_mean_L);
// fmm[i][j][l]->SetParameter(4,def_mean_S);



 //------- Fitting ----------//

//cout<<"fmmbg fit start"<<endl;
// hmm_bg_fom[i][j][l]->Fit("fmmbg[i][j][l]","Rq","",min_mm,max_mm);
// hmm_pi_wobg_fom[i][j][l]->Fit("fmmbg[i][j][l]","Rq","",-0.05,0.15);
//double fmmbga = fmmbg[i][j][l]->GetParameter(0);
//double fmmbgb = fmmbg[i][j][l]->GetParameter(1);
//double fmmbgc = fmmbg[i][j][l]->GetParameter(2);
//double fmmbgd = fmmbg[i][j][l]->GetParameter(3);
//cout<<"0:1:2:3:4:5="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<endl;//"::"<<e<<"::"<<f<<endl;
// fmm[i][j][l]->FixParameter(1,b);//mean
// fmm[i][j][l]->FixParameter(2,c);//sigma
// fmm[i][j][l]->FixParameter(3,d);//lg
//double e = fmmbg[i][j][l]->GetParameter(4);
//double f = fmmbg[i][j][l]->GetParameter(5);
// fmm[i][j][l]->SetParameter(0,a);
// fmm[i][j][l]->SetParameter(1,b);
// fmm[i][j][l]->SetParameter(2,c);
// fmm[i][j][l]->SetParameter(3,d);
// fmm[i][j][l]->SetParameter(4,e);
// fmm[i][j][l]->SetParameter(5,f);
// fmm[i][j][l]->SetParameter(6,500);
 fmm[i][j][l]->SetParameter(1,def_mean_L);
 fmm[i][j][l]->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
 fmm[i][j][l]->SetParameter(2,def_sig_L);
 fmm[i][j][l]->SetParLimits(2,0.,2*def_sig_L);
// fmm[i][j][l]->SetParameters(9,100);
 fmm[i][j][l]->SetParameter(4,def_mean_S);
 fmm[i][j][l]->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
 fmm[i][j][l]->SetParameter(5,def_sig_S);
 fmm[i][j][l]->SetParLimits(5,0.,2*def_sig_S);
 hmm_wo_bg_fom[i][j][l]->Fit("fmm[i][j][l]","","",-0.05,0.15);
//double fmmpar0 = fmm[i][j][l]->GetParameter(0);//cout<<"fmm[0]="<<fmmpar0<<endl;//area(L)
//double fmmpar1 = fmm[i][j][l]->GetParameter(1);//cout<<"fmm[1]="<<fmmpar1<<endl;//mean(L)
//double fmmpar2 = fmm[i][j][l]->GetParameter(2);//cout<<"fmm[2]="<<fmmpar2<<endl;//sigma(L)
//double fmmpar3 = fmm[i][j][l]->GetParameter(3);//cout<<"fmm[3]="<<fmmpar3<<endl;//area(S)
//double fmmpar4 = fmm[i][j][l]->GetParameter(4);//cout<<"fmm[4]="<<fmmpar4<<endl;//mean(S)
//double fmmpar5 = fmm[i][j][l]->GetParameter(5);//cout<<"fmm[5]="<<fmmpar5<<endl;//sigma(S)
double fmmpar6 = fmm[i][j][l]->GetParameter(6);//cout<<"fmm[6]="<<fmmpar6<<endl;//poly_const
double fmmpar7 = fmm[i][j][l]->GetParameter(7);//cout<<"fmm[7]="<<fmmpar7<<endl;//poly_x
double fmmpar8 = fmm[i][j][l]->GetParameter(8);//cout<<"fmm[8]="<<fmmpar8<<endl;//poly_x^2
double fmmpar9 = fmm[i][j][l]->GetParameter(9);//cout<<"fmm[9]="<<fmmpar9<<endl;//poly_x^3
double fmmpar10 = fmm[i][j][l]->GetParameter(10);//cout<<"fmm[10]="<<fmmpar10<<endl;//poly_x^4
 fmmbg[i][j][l]->SetParameters(fmmpar6,fmmpar7,fmmpar8,fmmpar9,fmmpar10);
//cout<<"0:1:2:3:4:5(as a total func)="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<"::"<<e<<"::"<<f<<endl;
//
//cout<<"fL fit start"<<endl;
// hmm_wo_bg_fom[i][j][l]->Fit("fL[i][j][l]","Rq","",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
// n_L[i][j][l]=fL[i][j][l]->GetParameter(0);
// mean_L[i][j][l]=fL[i][j][l]->GetParameter(1);
// sig_L[i][j][l]=fL[i][j][l]->GetParameter(2);
 mean_L[i][j][l]=def_mean_L;
 mean_S[i][j][l]=def_mean_S;
 sig_L[i][j][l]=def_sig_L;
 sig_S[i][j][l]=def_sig_S;
 n_L[i][j][l]=hmm_wo_bg_fom[i][j][l]->Integral(hmm_wo_bg_fom[i][j][l]->FindBin(center_L-range_L),hmm_wo_bg_fom[i][j][l]->FindBin(center_L+range_L));
cout<<"before(L):: "<<n_L[i][j][l]<<endl;
 double integralL=fmmbg[i][j][l]->Integral(center_L-range_L,center_L+range_L);
 if(integralL>0)n_L[i][j][l]-=integralL;
else cout<<"negative BG: ("<<i<<","<<j<<","<<l<<")"<<endl;
 n_L[i][j][l]-=fmmbg[i][j][l]->Integral(center_L-range_L,center_L+range_L);
cout<<"after(L):: "<<n_L[i][j][l]<<endl;
//n_L[i][j][l]-=(pow(mean_L[i][j][l]+2*sig_L[i][j][l],5)-pow(mean_L[i][j][l]-2*sig_L[i][j][l],5))*a/5;
//n_L[i][j][l]-=(pow(mean_L[i][j][l]+2*sig_L[i][j][l],4)-pow(mean_L[i][j][l]-2*sig_L[i][j][l],4))*b/4;
//n_L[i][j][l]-=(pow(mean_L[i][j][l]+2*sig_L[i][j][l],3)-pow(mean_L[i][j][l]-2*sig_L[i][j][l],3))*c/3;
//n_L[i][j][l]-=(pow(mean_L[i][j][l]+2*sig_L[i][j][l],2)-pow(mean_L[i][j][l]-2*sig_L[i][j][l],2))*d/2;
//n_L[i][j][l]-=(pow(mean_L[i][j][l]+2*sig_L[i][j][l],1)-pow(mean_L[i][j][l]-2*sig_L[i][j][l],1))*e;
//
 n_S[i][j][l]=hmm_wo_bg_fom[i][j][l]->Integral(hmm_wo_bg_fom[i][j][l]->FindBin(center_S-range_S),hmm_wo_bg_fom[i][j][l]->FindBin(center_S+range_S));
cout<<"before(S):: "<<n_S[i][j][l]<<endl;
 double integralS=fmmbg[i][j][l]->Integral(center_S-range_S,center_S+range_S);
 if(integralS>0)n_S[i][j][l]-=integralS;
else cout<<"negative BG: ("<<i<<","<<j<<","<<l<<")"<<endl;
cout<<"after(S):: "<<n_S[i][j][l]<<endl;
//n_S[i][j][l]-=(pow(mean_S[i][j][l]+2*sig_S[i][j][l],5)-pow(mean_S[i][j][l]-2*sig_S[i][j][l],5))*a/5;
//n_S[i][j][l]-=(pow(mean_S[i][j][l]+2*sig_S[i][j][l],4)-pow(mean_S[i][j][l]-2*sig_S[i][j][l],4))*b/4;
//n_S[i][j][l]-=(pow(mean_S[i][j][l]+2*sig_S[i][j][l],3)-pow(mean_S[i][j][l]-2*sig_S[i][j][l],3))*c/3;
//n_S[i][j][l]-=(pow(mean_S[i][j][l]+2*sig_S[i][j][l],2)-pow(mean_S[i][j][l]-2*sig_S[i][j][l],2))*d/2;
//n_S[i][j][l]-=(pow(mean_S[i][j][l]+2*sig_S[i][j][l],1)-pow(mean_S[i][j][l]-2*sig_S[i][j][l],1))*e;

 cout<<"n_L"<<n_L[i][j][l]<<endl;
cout<<"mean_L"<<mean_L[i][j][l]<<endl;
cout<<"sig_L"<<sig_L[i][j][l]<<endl;

////cout<<"fS fit start"<<endl;
// hmm_wo_bg_fom[i][j][l]->Fit("fS[i][j][l]","Rq","",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
//// n_S[i][j][l]=fS[i][j][l]->GetParameter(0);
// mean_S[i][j][l]=fS[i][j][l]->GetParameter(1);
// mean_S[i][j][l]=def_mean_S;
// sig_S[i][j][l]=def_sig_S;
// sig_S[i][j][l]=fS[i][j][l]->GetParameter(2);
// n_S[i][j][l]=hmm_wo_bg_fom[i][j][l]->Integral(hmm_wo_bg_fom[i][j][l]->FindBin(mean_S[i][j][l]-2*sig_S[i][j][l]),hmm_wo_bg_fom[i][j][l]->FindBin(mean_S[i][j][l]+2*sig_S[i][j][l]));
 cout<<"n_S"<<n_S[i][j][l]<<endl;
cout<<"mean_S"<<mean_S[i][j][l]<<endl;
cout<<"sig_S"<<sig_S[i][j][l]<<endl;
//----------------------------------------------//
//--	Missing Mass  End     ------------------//
//----------------------------------------------//

	
	//----------------------------------------------//
	//-----				DEBUG                  -----//
	//----------------------------------------------//
//	n_pi[i][j][l]=0.;n_pi_noAC1=0.;n_pi_noAC2=0.;//Pion ignoring
//	n_k[i][j][l]=0.;n_k_noAC1=0.;n_k_noAC2=0.;//Kaon ignoring
//	n_p[i][j][l]=0.;n_p_noAC1=0.;n_p_noAC2=0.;//Proton ingoring
	//----------------------------------------------//
	//-----				DEBUG                  -----//
	//----------------------------------------------//
	
//Changed
//	ac2l_adc[j]+=5.;	
//	if(ac2l_adc[j]>20)break;
	if(n_pi[i][j][l]>0.){}else{n_pi[i][j][l]=1.;}
	if(n_k[i][j][l]>0.){}else{n_k[i][j][l]=1.;}
	if(n_p[i][j][l]>0.){}else{n_p[i][j][l]=1.;}
	if(n_L[i][j][l]>0.){}else{n_L[i][j][l]=1.;}
	if(n_S[i][j][l]>0.){}else{n_S[i][j][l]=1.;}

	if(n_pi_noZ[i]>0.){}else{n_pi_noZ[i]=1.;}
	if(n_k_noZ[i]>0.){}else{n_k_noZ[i]=1.;}
	if(n_p_noZ[i]>0.){}else{n_p_noZ[i]=1.;}
	if(n_L_noZ[i]>0.){}else{n_L_noZ[i]=1.;}
	if(n_S_noZ[i]>0.){}else{n_S_noZ[i]=1.;}
//With Weight
	//h_pisr12l->Fill(ac1_adc[i],ac2l_adc[j],n_pi[i][j][l]/n_pi_noAC);
	//h_ksr12l->Fill(ac1_adc[i],ac2l_adc[j],n_k[i][j][l]/n_k_noAC);
	//h_psr12l->Fill(ac1_adc[i],ac2l_adc[j],n_p[i][j][l]/n_p_noAC);
	//h_Lsr12l->Fill(ac1_adc[i],ac2l_adc[j],n_L[i][j][l]/n_L_noAC);
	//h_Ssr12l->Fill(ac1_adc[i],ac2l_adc[j],n_S[i][j][l]/n_S_noAC);
	for(int f=0;f<1000*n_pi[i][j][l]/n_pi_noZ[i];f++){h_pisr2d->Fill(zver[j],zver_diff[i]);}
	for(int f=0;f<1000*n_k[i][j][l]/n_k_noZ[i];f++){h_ksr2d->Fill(zver[j],zver_diff[i]);}
	for(int f=0;f<1000*n_p[i][j][l]/n_p_noZ[i];f++){h_psr2d->Fill(zver[j],zver_diff[i]);}
	for(int f=0;f<1000*n_L[i][j][l]/n_L_noZ[i];f++){h_Lsr2d->Fill(zver[j],zver_diff[i]);}
	for(int f=0;f<1000*n_S[i][j][l]/n_S_noZ[i];f++){h_Ssr2d->Fill(zver[j],zver_diff[i]);}
//
//	if(n_pi[i][j][l]>n_pi_noAC){
//		if(j==0 && l==0){for(int fill=0;fill<n_pi_noAC;fill++){h_pisr1->Fill(ac1_adc[i]);}}
//	}
//	else {if(j==0 && l==0){for(int fill=0;fill<n_pi[i][j][l];fill++){h_pisr1->Fill(ac1_adc[i]);}}}
//	if(n_pi[i][j][l]>n_pi_noAC){
//		//if(i==0 && l==0){for(int fill=0;fill<n_pi_noAC;fill++){h_pisr2l->Fill(ac2l_adc[j]);}}
//		if(l==0){for(int fill=0;fill<n_pi_noAC;fill++){h_pisr12l->Fill(ac1_adc[i],ac2l_adc[j]);h_pisr2l->Fill(ac2l_adc[j])}}
//		if(i==0 && j==0){for(int fill=0;fill<n_pi_noAC;fill++){h_pisr2u->Fill(ac2u_adc[l]);}}
//	}
//	else {
//		//if(i==0 && l==0){for(int fill=0;fill<n_pi[i][j][l];fill++){h_pisr2l->Fill(ac2l_adc[j]);}}
//		if(l==0){for(int fill=0;fill<n_pi[i][j][l];fill++){h_pisr12l->Fill(ac1_adc[i],ac2l_adc[j]);h_pisr2l->Fill(ac2l_adc[j]);}}
//		if(i==0 && j==0){for(int fill=0;fill<n_pi[i][j][l];fill++){h_pisr2u->Fill(ac2u_adc[l]);}}
//	}
//
//	if(n_k[i][j][l]>n_k_noAC){
//		if(j==0 && l==0){for(int fill=0;fill<n_k_noAC;fill++){h_ksr1->Fill(ac1_adc[i]);}}
//	}
//	else {if(j==0 && l==0){for(int fill=0;fill<n_k[i][j][l];fill++){h_ksr1->Fill(ac1_adc[i]);}}}
//	if(n_k[i][j][l]>n_k_noAC){
//		//if(i==0 && l==0){for(int fill=0;fill<n_k_noAC;fill++){h_ksr2l->Fill(ac2l_adc[j]);}}
//		if(l==0){for(int fill=0;fill<n_k_noAC;fill++){h_ksr12l->Fill(ac1_adc[i],ac2l_adc[j]);h_ksr2l->Fill(ac2l_adc[j]);}}
//		if(i==0 && j==0){for(int fill=0;fill<n_k_noAC;fill++){h_ksr2u->Fill(ac2u_adc[l]);}}
//	}
//	else {
//		//if(i==0 && l==0){for(int fill=0;fill<n_k[i][j][l];fill++){h_ksr2l->Fill(ac2l_adc[j]);}}
//		if(l==0){for(int fill=0;fill<n_k[i][j][l];fill++){h_ksr12l->Fill(ac1_adc[i],ac2l_adc[j]);h_ksr2l->Fill(ac2l_adc[j]);}}
//		if(i==0 && j==0){for(int fill=0;fill<n_k[i][j][l];fill++){h_ksr2u->Fill(ac2u_adc[l]);}}
//	}
//
//	if(n_p[i][j][l]>n_p_noAC){
//		if(j==0 && l==0){for(int fill=0;fill<n_p_noAC;fill++){h_psr1->Fill(ac1_adc[i]);}}
//	}
//	else {if(j==0 && l==0){for(int fill=0;fill<n_p[i][j][l];fill++){h_psr1->Fill(ac1_adc[i]);}}}
//	if(n_p[i][j][l]>n_p_noAC){
//		//if(i==0 && l==0){for(int fill=0;fill<n_p_noAC;fill++){h_psr2l->Fill(ac2l_adc[j]);}}
//		if(l==0){for(int fill=0;fill<n_p_noAC;fill++){h_psr12l->Fill(ac1_adc[i],ac2l_adc[j]);h_psr2l->Fill(ac2l_adc[j]);}}
//		if(i==0 && j==0){for(int fill=0;fill<n_p_noAC;fill++){h_psr2u->Fill(ac2u_adc[l]);}}
//	}
//	else {
//		//if(i==0 && l==0){for(int fill=0;fill<n_p[i][j][l];fill++){h_psr2l->Fill(ac2l_adc[j]);}}
//		if(l==0){for(int fill=0;fill<n_p[i][j][l];fill++){h_psr12l->Fill(ac1_adc[i],ac2l_adc[j]);h_psr2l->Fill(ac2l_adc[j]);}}
//		if(i==0 && j==0){for(int fill=0;fill<n_p[i][j][l];fill++){h_psr2u->Fill(ac2u_adc[l]);}}
//	}
//
//	if(n_L[i][j][l]>n_L_noAC){
//		if(j==0 && l==0){for(int fill=0;fill<n_L_noAC;fill++){h_Lsr1->Fill(ac1_adc[i]);}}
//	}
//	else {if(j==0 && l==0){for(int fill=0;fill<n_L[i][j][l];fill++){h_Lsr1->Fill(ac1_adc[i]);}}}
//	if(n_L[i][j][l]>n_L_noAC){
//		//if(i==0 && l==0){for(int fill=0;fill<n_L_noAC;fill++){h_Lsr2l->Fill(ac2l_adc[j]);}}
//		if(l==0){for(int fill=0;fill<n_L_noAC;fill++){h_Lsr12l->Fill(ac1_adc[i],ac2l_adc[j]);h_Lsr2l->Fill(ac2l_adc[j]);}}
//		if(i==0 && j==0){for(int fill=0;fill<n_L_noAC;fill++){h_Lsr2u->Fill(ac2u_adc[l]);}}
//	}
//	else {
//		//if(i==0 && l==0){for(int fill=0;fill<n_L[i][j][l];fill++){h_Lsr2l->Fill(ac2l_adc[j]);}}
//		if(l==0){for(int fill=0;fill<n_L[i][j][l];fill++){h_Lsr12l->Fill(ac1_adc[i],ac2l_adc[j]);h_Lsr2l->Fill(ac2l_adc[j]);}}
//		if(i==0 && j==0){for(int fill=0;fill<n_L[i][j][l];fill++){h_Lsr2u->Fill(ac2u_adc[l]);}}
//	}
//
//	if(n_S[i][j][l]>n_S_noAC){
//		if(j==0 && l==0){for(int fill=0;fill<n_S_noAC;fill++){h_Ssr1->Fill(ac1_adc[i]);}}
//	}
//	else {if(j==0 && l==0){for(int fill=0;fill<n_S[i][j][l];fill++){h_Ssr1->Fill(ac1_adc[i]);}}}
//	if(n_S[i][j][l]>n_S_noAC){
//		//if(i==0 && l==0){for(int fill=0;fill<n_S_noAC;fill++){h_Ssr2l->Fill(ac2l_adc[j]);}}
//		if(l==0){for(int fill=0;fill<n_S_noAC;fill++){h_Ssr12l->Fill(ac1_adc[i],ac2l_adc[j]);h_Ssr2l->Fill(ac2l_adc[j]);}}
//		if(i==0 && j==0){for(int fill=0;fill<n_S_noAC;fill++){h_Ssr2u->Fill(ac2u_adc[l]);}}
//	}
//	else {
//		//if(i==0 && l==0){for(int fill=0;fill<n_S[i][j][l];fill++){h_Ssr2l->Fill(ac2l_adc[j]);}}
//		if(l==0){for(int fill=0;fill<n_S[i][j][l];fill++){h_Ssr12l->Fill(ac1_adc[i],ac2l_adc[j]);h_Ssr2l->Fill(ac2l_adc[j]);}}
//		if(i==0 && j==0){for(int fill=0;fill<n_S[i][j][l];fill++){h_Ssr2u->Fill(ac2u_adc[l]);}}
//	}



//------------06.06, 1D Effciency---------------//
//	j==20 (zver=0.1), i==25 (zver_diff=0.025) 
	if(n_pi[i][j][l]>n_pi_noZ[i]){
		if(i==25 && l==0){for(int fill=0;fill<n_pi_noZ[i];fill++){h_pisr1->Fill(zver[j]);}}
	}
	else {if(i==25 && l==0){for(int fill=0;fill<n_pi[i][j][l];fill++){h_pisr1->Fill(zver[j]);}}}
	if(n_pi[i][j][l]>n_pi_noZ[i]){
		if(j==20 && l==0){for(int fill=0;fill<n_pi_noZ[i];fill++){h_pisr2l->Fill(zver_diff[i]);}}
	}
	else {
		if(j==20 && l==0){for(int fill=0;fill<n_pi[i][j][l];fill++){h_pisr2l->Fill(zver_diff[i]);}}
	}

	if(n_k[i][j][l]>n_k_noZ[i]){
		if(i==25 && l==0){for(int fill=0;fill<n_k_noZ[i];fill++){h_ksr1->Fill(zver[j]);}}
	}
	else {if(i==25 && l==0){for(int fill=0;fill<n_k[i][j][l];fill++){h_ksr1->Fill(zver[j]);}}}
	if(n_k[i][j][l]>n_k_noZ[i]){
		if(j==20 && l==0){for(int fill=0;fill<n_k_noZ[i];fill++){h_ksr2l->Fill(zver_diff[i]);}}
	}
	else {
		if(j==20 && l==0){for(int fill=0;fill<n_k[i][j][l];fill++){h_ksr2l->Fill(zver_diff[i]);}}
	}

	if(n_p[i][j][l]>n_p_noZ[i]){
		if(i==25 && l==0){for(int fill=0;fill<n_p_noZ[i];fill++){h_psr1->Fill(zver[j]);}}
	}
	else {if(i==25 && l==0){for(int fill=0;fill<n_p[i][j][l];fill++){h_psr1->Fill(zver[j]);}}}
	if(n_p[i][j][l]>n_p_noZ[i]){
		if(j==20 && l==0){for(int fill=0;fill<n_p_noZ[i];fill++){h_psr2l->Fill(zver_diff[i]);}}
	}
	else {
		if(j==20 && l==0){for(int fill=0;fill<n_p[i][j][l];fill++){h_psr2l->Fill(zver_diff[i]);}}
	}

	if(n_L[i][j][l]>n_L_noZ[i]){
		if(i==25 && l==0){for(int fill=0;fill<n_L_noZ[i];fill++){h_Lsr1->Fill(zver[j]);}}
	}
	else {if(i==25 && l==0){for(int fill=0;fill<n_L[i][j][l];fill++){h_Lsr1->Fill(zver[j]);}}}
	if(n_L[i][j][l]>n_L_noZ[i]){
		if(j==20 && l==0){for(int fill=0;fill<n_L_noZ[i];fill++){h_Lsr2l->Fill(zver_diff[i]);}}
	}
	else {
		if(j==20 && l==0){for(int fill=0;fill<n_L[i][j][l];fill++){h_Lsr2l->Fill(zver_diff[i]);}}
	}

	if(n_S[i][j][l]>n_S_noZ[i]){
		if(i==25 && l==0){for(int fill=0;fill<n_S_noZ[i];fill++){h_Ssr1->Fill(zver[j]);}}
	}
	else {if(i==25 && l==0){for(int fill=0;fill<n_S[i][j][l];fill++){h_Ssr1->Fill(zver[j]);}}}
	if(n_S[i][j][l]>n_S_noZ[i]){
		if(j==20 && l==0){for(int fill=0;fill<n_S_noZ[i];fill++){h_Ssr2l->Fill(zver_diff[i]);}}
	}
	else {
		if(j==20 && l==0){for(int fill=0;fill<n_S[i][j][l];fill++){h_Ssr2l->Fill(zver_diff[i]);}}
	}
//------------06.06, 1D Effciency---------------//
	if(n_k[i][j][l]>n_k_noZ[i]){
		if(i==10 && l==0){for(int fill=0;fill<n_k_noZ[i];fill++){h_ksr11->Fill(zver[j]);}}
	}
	else {if(i==10 && l==0){for(int fill=0;fill<n_k[i][j][l];fill++){h_ksr11->Fill(zver[j]);}}}
	if(n_k[i][j][l]>n_k_noZ[i]){
		if(i==50 && l==0){for(int fill=0;fill<n_k_noZ[i];fill++){h_ksr111->Fill(zver[j]);}}
	}
	else {if(i==50 && l==0){for(int fill=0;fill<n_k[i][j][l];fill++){h_ksr111->Fill(zver[j]);}}}
	if(n_k[i][j][l]>n_k_noZ[i]){
		if(i==99 && l==0){for(int fill=0;fill<n_k_noZ[i];fill++){h_ksr1111->Fill(zver[j]);}}
	}
	else {if(i==99 && l==0){for(int fill=0;fill<n_k[i][j][l];fill++){h_ksr1111->Fill(zver[j]);}}}


	
	if(i==25 && l==0){
	for(int fill=0;fill<n_pi_noZ[i];fill++) h_pitot1->Fill(zver[j]);
	for(int fill=0;fill<n_k_noZ[i];fill++) h_ktot1->Fill(zver[j]);
	for(int fill=0;fill<n_p_noZ[i];fill++) h_ptot1->Fill(zver[j]);
	for(int fill=0;fill<n_L_noZ[i];fill++) h_Ltot1->Fill(zver[j]);
	for(int fill=0;fill<n_S_noZ[i];fill++) h_Stot1->Fill(zver[j]);
			}//denominator
	if(i==10 && l==0)for(int fill=0;fill<n_k_noZ[i];fill++) h_ktot11->Fill(zver[j]);
	if(i==50 && l==0)for(int fill=0;fill<n_k_noZ[i];fill++) h_ktot111->Fill(zver[j]);
	if(i==99 && l==0)for(int fill=0;fill<n_k_noZ[i];fill++) h_ktot1111->Fill(zver[j]);
	if(j==20 && l==0){
	for(int fill=0;fill<n_pi_noZ[i];fill++) h_pitot2l->Fill(zver_diff[i]);
	for(int fill=0;fill<n_k_noZ[i];fill++) h_ktot2l->Fill(zver_diff[i]);
	for(int fill=0;fill<n_p_noZ[i];fill++) h_ptot2l->Fill(zver_diff[i]);
	for(int fill=0;fill<n_L_noZ[i];fill++) h_Ltot2l->Fill(zver_diff[i]);
	for(int fill=0;fill<n_S_noZ[i];fill++) h_Stot2l->Fill(zver_diff[i]);
			}//denominator

//			}//for l
		}//for j
	}//for i

ofstream fout("SR_z2d.dat");
cout<<"Z Efficiency is filled."<<endl;
//fout<<n_pi_noAC<<" "<<n_k_noAC<<" "<<n_p_noAC<<" "<<n_L_noAC<<" "<<n_S_noAC<<endl;
		double th1,th2;
for(int i=0;i<nth;i++){
	for(int j=0;j<nth;j++){
		int l=0;
		th1=zver[j];
		th2=zver_diff[i];
		fout<<n_pi[i][j][l]/n_pi_noZ[i]<<" "<<n_k[i][j][l]/n_k_noZ[i]<<" "<<n_p[i][j][l]/n_p_noZ[i]<<" "<<n_L[i][j][l]/n_L_noZ[i]<<" "<<n_S[i][j][l]/n_S_noZ[i]<<" z_sum<"<<th1<<", z_diff<"<<th2<<endl;
	}
}
		


cout << "TEfficiency!" << endl;
//-----------Z_sum-------------------//
cout<<"pEff1:"<<endl;
if(TEfficiency::CheckConsistency(*h_pisr1,*h_pitot1,"w")){
pEff1 = new TEfficiency(*h_pisr1,*h_pitot1);
//pEff1->Write();
}
cout<<"pEff2:"<<endl;
if(TEfficiency::CheckConsistency(*h_ksr1,*h_ktot1,"w")){
pEff2 = new TEfficiency(*h_ksr1,*h_ktot1);
//pEff2->Write();
}
cout<<"pEff3:"<<endl;
if(TEfficiency::CheckConsistency(*h_psr1,*h_ptot1,"w")){
pEff3 = new TEfficiency(*h_psr1,*h_ptot1);
//pEff3->Write();
}
cout<<"pEff4:"<<endl;
if(TEfficiency::CheckConsistency(*h_Lsr1,*h_Ltot1,"w")){
pEff4 = new TEfficiency(*h_Lsr1,*h_Ltot1);
//pEff4->Write();
}
cout<<"pEff5:"<<endl;
if(TEfficiency::CheckConsistency(*h_Ssr1,*h_Stot1,"w")){
pEff5 = new TEfficiency(*h_Ssr1,*h_Stot1);
//pEff5->Write();
}
//-----------Z_diff-------------------//
cout<<"pEff6:"<<endl;
if(TEfficiency::CheckConsistency(*h_pisr2l,*h_pitot2l,"w")){
pEff6 = new TEfficiency(*h_pisr2l,*h_pitot2l);
//pEff6->Write();
}
cout<<"pEff7:"<<endl;
if(TEfficiency::CheckConsistency(*h_ksr2l,*h_ktot2l,"w")){
pEff7 = new TEfficiency(*h_ksr2l,*h_ktot2l);
//pEff7->Write();
}
cout<<"pEff8:"<<endl;
if(TEfficiency::CheckConsistency(*h_psr2l,*h_ptot2l,"w")){
pEff8 = new TEfficiency(*h_psr2l,*h_ptot2l);
//pEff8->Write();
}
cout<<"pEff9:"<<endl;
if(TEfficiency::CheckConsistency(*h_Lsr2l,*h_Ltot2l,"w")){
pEff9 = new TEfficiency(*h_Lsr2l,*h_Ltot2l);
//pEff9->Write();
}
cout<<"pEff10:"<<endl;
if(TEfficiency::CheckConsistency(*h_Ssr2l,*h_Stot2l,"w")){
pEff10 = new TEfficiency(*h_Ssr2l,*h_Stot2l);
//pEff10->Write();
}

cout<<"pEff11:"<<endl;
if(TEfficiency::CheckConsistency(*h_ksr11,*h_ktot11,"w")){
pEff11 = new TEfficiency(*h_ksr11,*h_ktot11);
//pEff11->Write();
}
cout<<"pEff111:"<<endl;
if(TEfficiency::CheckConsistency(*h_ksr111,*h_ktot111,"w")){
pEff111 = new TEfficiency(*h_ksr111,*h_ktot111);
//pEff111->Write();
}
cout<<"pEff1111:"<<endl;
if(TEfficiency::CheckConsistency(*h_ksr1111,*h_ktot1111,"w")){
pEff1111 = new TEfficiency(*h_ksr1111,*h_ktot1111);
//pEff1111->Write();
}

cout << "TGraphErrors!" << endl;
//gcoin_pi_sr[0] = new TGraphErrors(99, ac1_adc, fom_pi1, 0, err_fom_pi1);
gcoin_pi_sr[0] = new TGraphErrors(nth-1, ac1_adc, fom_pi1, 0, 0);
gcoin_pi_sr[0]->SetMarkerColor(kOrange);
gcoin_pi_sr[0]->SetMarkerStyle(22);
gcoin_pi_sr[0]->SetMarkerSize(1);
//gcoin_k_sr[0] = new TGraphErrors(99, ac1_adc, fom_k1, 0, err_fom_k1);
gcoin_k_sr[0] = new TGraphErrors(nth-1, ac1_adc, fom_k1, 0, 0);
gcoin_k_sr[0]->SetMarkerColor(kGreen);
gcoin_k_sr[0]->SetMarkerStyle(22);
gcoin_k_sr[0]->SetMarkerSize(1);
//gcoin_p_sr[0] = new TGraphErrors(99, ac1_adc, fom_p1, 0, err_fom_p1);
gcoin_p_sr[0] = new TGraphErrors(nth-1, ac1_adc, fom_p1, 0, 0);
gcoin_p_sr[0]->SetMarkerColor(kRed);
gcoin_p_sr[0]->SetMarkerStyle(22);
gcoin_p_sr[0]->SetMarkerSize(1);
//gcoin_pi_sr[1] = new TGraphErrors(100, ac2l_adc, fom_pi2l, 0, 0);
//gcoin_pi_sr[2] = new TGraphErrors(100, ac2u_adc, fom_pi2u, 0, 0);
//gcoin_pi_sr[3] = new TGraphErrors(10, ac2u_adc, fom_pi1, 0, 0);
//gcoin_pi_sr[4] = new TGraphErrors(10, ac2u_adc, fom_pi1, 0, 0);
//gcoin_pi_sr[5] = new TGraphErrors(10, ac2u_adc, fom_pi1, 0, 0);
//gcoin_pi_sr[6] = new TGraphErrors(10, ac2u_adc, fom_pi1, 0, 0);
//gcoin_pi_sr[7] = new TGraphErrors(10, ac2u_adc, fom_pi1, 0, 0);
//gcoin_pi_sr[8] = new TGraphErrors(10, ac2u_adc, fom_pi1, 0, 0);
//gcoin_pi_sr[9] = new TGraphErrors(10, ac2u_adc, fom_pi1, 0, 0);




cout << "After fill" << endl;

	
}

///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

void tuning::Draw(){

cout<< "Draw Start" << endl;

//Test
c1 = new TCanvas("c1","c1",800.,800.);
c2 = new TCanvas("c2","c2",800.,800.);
c3 = new TCanvas("c3","c3",800.,800.);
c4 = new TCanvas("c4","c4",800.,800.);
c5 = new TCanvas("c5","c5",800.,800.);
c6 = new TCanvas("c6","c6",800.,800.);
c7 = new TCanvas("c7","c7",800.,800.);
c8 = new TCanvas("c8","c8",800.,800.);
c9 = new TCanvas("c9","c9",800.,800.);
c10 = new TCanvas("c10","c10",800.,800.);
c11 = new TCanvas("c11","c11",800.,800.);
c12 = new TCanvas("c12","c12",800.,800.);
c13 = new TCanvas("c13","c13",800.,800.);
c14 = new TCanvas("c14","c14",800.,800.);
c15 = new TCanvas("c15","c15",800.,800.);
c16 = new TCanvas("c16","c16",800.,800.);
c17 = new TCanvas("c17","c17",800.,800.);
c18 = new TCanvas("c18","c18",800.,800.);
c19 = new TCanvas("c19","c19",800.,800.);

c1->Divide(3,2);
//cout << "start hcoin_tc" << endl;
c1->cd(1);hcoin_wo_bg_fom_noAC1->Draw("");
c1->cd(2);hcoin_wo_bg_fom_noAC2->Draw("");
c1->cd(3);hcoin_wo_bg_fom_noAC->Draw("");
double ymax = (hcoin_wo_bg_fom_noAC->GetBinContent(hcoin_wo_bg_fom_noAC->GetMaximumBin()));
TLine *tl1, *tl2, *tl3, *tl4, *tl5, *tl6;
	tl1 = new TLine(center_pi-range_pi,-20,center_pi-range_pi,0.5*ymax);
	tl1->SetLineWidth(1);
	tl1->SetLineColor(kOrange);
	tl1->Draw("same");
	tl2 = new TLine(center_pi+range_pi,-20,center_pi+range_pi,0.5*ymax);
	tl2->SetLineWidth(1);
	tl2->SetLineColor(kOrange);
	tl2->Draw("same");
	tl3 = new TLine(center_p-range_p,-20,center_p-range_p,0.5*ymax);
	tl3->SetLineWidth(1);
	tl3->SetLineColor(kRed);
	tl3->Draw("same");
	tl4 = new TLine(center_p+range_p,-20,center_p+range_p,0.5*ymax);
	tl4->SetLineWidth(1);
	tl4->SetLineColor(kRed);
	tl4->Draw("same");
	tl5 = new TLine(center_k-range_k,-20,center_k-range_k,0.5*ymax);
	tl5->SetLineWidth(1);
	tl5->SetLineColor(kGreen);
	tl5->Draw("same");
	tl6 = new TLine(center_k+range_k,-20,center_k+range_k,0.5*ymax);
	tl6->SetLineWidth(1);
	tl6->SetLineColor(kGreen);
	tl6->Draw("same");
c1->cd(4);hmm_L_fom_noAC1->Draw("");
c1->cd(5);hmm_L_fom_noAC2->Draw("");
c1->cd(6);hmm_L_fom_noAC->Draw("");
    ymax = (hmm_L_fom_noAC->GetBinContent(hmm_L_fom_noAC->GetMaximumBin()));
TLine *tl7, *tl8, *tl9, *tl10;
	tl7 = new TLine(center_L-range_L,-20,center_L-range_L,0.5*ymax);
	tl7->SetLineWidth(1);
	tl7->SetLineColor(kAzure);
	tl7->Draw("same");
	tl8 = new TLine(center_L+range_L,-20,center_L+range_L,0.5*ymax);
	tl8->SetLineWidth(1);
	tl8->SetLineColor(kAzure);
	tl8->Draw("same");
	tl9 = new TLine(center_S-range_S,-20,center_S-range_S,0.5*ymax);
	tl9->SetLineWidth(1);
	tl9->SetLineColor(kCyan);
	tl9->Draw("same");
	tl10 = new TLine(center_S+range_S,-20,center_S+range_S,0.5*ymax);
	tl10->SetLineWidth(1);
	tl10->SetLineColor(kCyan);
	tl10->Draw("same");
//c1->cd(1);h_ct->Draw("");
//c1->cd(2);hcoin_wo_bg_fom[nth][nth][nth]->Draw("");fcoin[nth][nth][nth]->Draw("same");
//c1->cd(3);h_ct_wK->Draw("");fk_kc->SetLineColor(kRed);fk_kc->Draw("same");
//c1->cd(4);h_mmall->Draw("");
////c1->cd(5);h_mmfoil->Draw("");
////c1->cd(6);h_mm->Draw("");
////c1->cd(4);hmm_bg_fom[0][0][0]->Draw("");
//c1->cd(4);hcoin_wo_bg_fom[8][8][8]->Draw("");fcoin[8][8][8]->SetLineColor(kRed);fcoin[8][8][8]->Draw("same");
////c1->cd(5);hmm_wo_bg_fom[0][0][0]->Draw("");
//c1->cd(5);hmm_L_fom[8][8][8]->Draw("");fmm[8][8][8]->SetLineColor(kRed);fmm[8][8][8]->Draw("same");
//c1->cd(6);h3_fom->Draw("box2 z");
////h3_fom->GetXaxis()->SetRange(ac1_adc[0],ac1_adc[9]);
//TH2D *hProjectionyz = (TH2D*) h3_fom->Project3D("yz");
////h3_fom->GetYaxis()->SetRange(ac2l_adc[0],ac2l_adc[9]);
//TH2D *hProjectionxz = (TH2D*) h3_fom->Project3D("xz");
////h3_fom->GetZaxis()->SetRange(ac2u_adc[0],ac2u_adc[9]);
//TH2D *hProjectionxy = (TH2D*) h3_fom->Project3D("xy");
//c1->cd(7);hProjectionyz->Draw("colz");
//c1->cd(8);hProjectionxz->Draw("colz");
//c1->cd(9);hProjectionxy->Draw("colz");
//c1->cd(1);hcoin_pi->Draw("");fpi_pic->SetLineColor(kRed);fpi_pic->Draw("same");
//c1->cd(2);hcoin_k->Draw("");fk_kc->SetLineColor(kRed);fk_kc->Draw("same");
//c1->cd(3);hcoin_p->Draw("");fp_pc->SetLineColor(kRed);fp_pc->Draw("same");
//for(int i=1;i<=3;i++){
//c1->cd(i+3);hcoin_k_ac2[i-1]->Draw("");fkk[i-1]->SetLineColor(kRed);fkk[i-1]->Draw("same");
//}

c2->Divide(3,2);
cout<<"c2 start"<<endl;
c2->cd(1);
hmm_wo_bg_fom_noAC1->Draw("");
c2->cd(2);
hmm_wo_bg_fom_noAC2->Draw("");
c2->cd(3);
hmm_wo_bg_fom_noAC->Draw("");
fmm_noAC->SetLineColor(kRed);
fmmbg_noAC->SetLineColor(kGreen);
fmm_noAC->Draw("same");
fmmbg_noAC->Draw("same");
//hmm_L_fom_noAC->Draw("");
//hmm_pi_wobg_fom_noAC->SetLineColor(kRed);hmm_pi_wobg_fom_noAC->Draw("same");
c2->cd(4);
hmm_L_fom_noAC1->Draw("");hmm_bg_fom_noAC1->SetLineColor(kRed);hmm_bg_fom_noAC1->Draw("same");
c2->cd(5);
hmm_L_fom_noAC2->Draw("");hmm_bg_fom_noAC2->SetLineColor(kRed);hmm_bg_fom_noAC2->Draw("same");
c2->cd(6);
hmm_L_fom_noAC->Draw("");hmm_bg_fom_noAC->SetLineColor(kRed);hmm_bg_fom_noAC->Draw("same");
c3->cd()->DrawFrame(0.,0.,500,1.2);//K,pi,p vs SUM 
pEff1->SetLineColor(kOrange);pEff1->Draw("same");
pEff2->SetLineColor(kGreen);pEff2->Draw("same");
pEff3->SetLineColor(kRed);pEff3->Draw("same");
c4->cd()->DrawFrame(0.,0.,500,1.2);
pEff4->SetLineColor(kAzure);pEff4->Draw("same");
pEff5->SetLineColor(kCyan);pEff5->Draw("same");
c5->cd()->DrawFrame(0.,0.,500,1.2);
pEff2->SetLineColor(kGreen);pEff2->Draw("same");
pEff4->SetLineColor(kAzure);pEff4->Draw("same");
pEff5->SetLineColor(kCyan);pEff5->Draw("same");
c6->cd()->DrawFrame(0.,0.,100,1.2);//K,pi,p vs DIFF 
pEff6->SetLineColor(kOrange);pEff6->Draw("same");
pEff7->SetLineColor(kGreen);pEff7->Draw("same");
pEff8->SetLineColor(kRed);pEff8->Draw("same");
c7->cd()->DrawFrame(0.,0.,100,1.2);
pEff9->SetLineColor(kAzure);pEff9->Draw("same");
pEff10->SetLineColor(kCyan);pEff10->Draw("same");
c8->cd()->DrawFrame(0.,0.,100,1.2);
pEff7->SetLineColor(kGreen);pEff7->Draw("same");
pEff9->SetLineColor(kAzure);pEff9->Draw("same");
pEff10->SetLineColor(kCyan);pEff10->Draw("same");
c9->cd()->DrawFrame(0.,0.,500,1.2);
pEff2->SetLineColor(kGreen);pEff2->Draw("same");
pEff11->SetLineColor(kGreen+1);pEff11->Draw("same");
pEff111->SetLineColor(kGreen+2);pEff111->Draw("same");
pEff1111->SetLineColor(kGreen+3);pEff1111->Draw("same");
c16->Divide(2,2);
c16->cd(1);h_pisr2d->Draw("colz");
c16->cd(2);h_ksr2d->Draw("colz");
c16->cd(3);h_Lsr2d->Draw("colz");
c16->cd(4);h_Ssr2d->Draw("colz");
c17->Divide(2,2);
c17->cd(1);h_pisr2d->Draw("surf3z");
c17->cd(2);h_ksr2d->Draw("surf3z");
c17->cd(3);h_Lsr2d->Draw("surf3z");
c17->cd(4);h_Ssr2d->Draw("surf3z");
c18->cd();
h_Lsr2d->Draw("colz");
//c19->cd();
//h_Ssr12l->Draw("colz");
//h_ptot1->Draw("sh");h_psr1->Draw("samesh");
//npe_sum_a1->Draw("");
//hcoin_wo_bg_fom[80][0][0]->Draw("");fcoin[80][0][0]->SetLineColor(kRed);fcoin[80][0][0]->Draw("same");
//c13->cd()->SetLogy(1);
//npe_sum_a2->Draw("");
//hcoin_wo_bg_fom[nth][nth][nth]->Draw("");fcoin[nth][nth][nth]->SetLineColor(kRed);fcoin[nth][nth][nth]->Draw("same");
//c2->Divide(3,3);
//c3->Divide(3,3);
//c4->Divide(3,3);
//for(int i=0;i<3;i++){
//	for(int j=0;j<3;j++){
//		for(int l=0;l<3;l++){
//			if(i==0)c2->cd(j*3+l+1);
//			if(i==1)c3->cd(j*3+l+1);
//			if(i==2)c4->cd(j*3+l+1);
//			hcoin_wo_bg_fom[i][j][l]->Draw("");fcoin[i][j][l]->SetLineColor(kRed);fcoin[i][j][l]->Draw("same");
//		}
//	}
//}
//c5->Divide(9,3);
//c6->Divide(9,3);
//c7->Divide(9,3);
//c8->Divide(9,3);
//c9->Divide(9,3);
//c10->Divide(9,3);
//c11->Divide(9,3);
//c12->Divide(9,3);
//c13->Divide(9,3);
//for(int i=0;i<9;i++){
//	for(int j=0;j<3;j++){
//		for(int l=0;l<9;l++){
//			//c3->cd(i*9+j*3+l+1);
//			if(i==0)c5->cd(j*9+l+1);
//			if(i==1)c6->cd(j*9+l+1);
//			if(i==2)c7->cd(j*9+l+1);
//			if(i==3)c8->cd(j*9+l+1);
//			if(i==4)c9->cd(j*9+l+1);
//			if(i==5)c10->cd(j*9+l+1);
//			if(i==6)c11->cd(j*9+l+1);
//			if(i==7)c12->cd(j*9+l+1);
//			if(i==8)c13->cd(j*9+l+1);
//			hmm_L_fom[i][j][l]->Draw("");fmm[i][j][l]->SetLineColor(kRed);fmm[i][j][l]->Draw("same");
//		}
//	}
//}

}
//////////////////////////////////////////////////////////




void tuning::Print(string ofname){



  cout<<"Print is starting"<<endl;
  cout<<"pdf name: "<<ofname<<endl;
 c1->Print(Form("%s[",ofname.c_str()));
 c1->Print(Form("%s",ofname.c_str()));
 c2->Print(Form("%s",ofname.c_str()));
 c3->Print(Form("%s",ofname.c_str()));
 c4->Print(Form("%s",ofname.c_str()));
 c5->Print(Form("%s",ofname.c_str()));
 c6->Print(Form("%s",ofname.c_str()));
 c7->Print(Form("%s",ofname.c_str()));
 c8->Print(Form("%s",ofname.c_str()));
 c9->Print(Form("%s",ofname.c_str()));
 c10->Print(Form("%s",ofname.c_str()));
 c11->Print(Form("%s",ofname.c_str()));
 c12->Print(Form("%s",ofname.c_str()));
 c13->Print(Form("%s",ofname.c_str()));
 c14->Print(Form("%s",ofname.c_str()));
 c15->Print(Form("%s",ofname.c_str()));
 c16->Print(Form("%s",ofname.c_str()));
 c17->Print(Form("%s",ofname.c_str()));
 c18->Print(Form("%s",ofname.c_str())); 
 c19->Print(Form("%s",ofname.c_str())); 
 c19->Print(Form("%s]",ofname.c_str()));
//      cout<<"c18 is done"<<endl;       
// c21->Print(Form("%s",ofname.c_str()));  
// c22->Print(Form("%s",ofname.c_str())); 
// c23->Print(Form("%s",ofname.c_str()));
// c24->Print(Form("%s",ofname.c_str()));
// c25->Print(Form("%s",ofname.c_str())); 
// c30->Print(Form("%s",ofname.c_str()));
// c30->Print(Form("%s]",ofname.c_str()));
//// c32->Print(Form("%s",ofname.c_str()));
//// c32->Print(Form("%s]",ofname.c_str()));
 
    
 cout<<"Print is done "<<endl;
   


}


///////////////////////////////////////////////////////////////

void tuning::Write(){


// gL_ac1[fom_max_th2]->SetName(Form("gL_ac1_%d",fom_max_th2));
// gL_ac1[fom_max_th2]->Write();
// gS_ac1[fom_max_th2]->SetName(Form("gS_ac1_%d",fom_max_th2)); 
// gS_ac1[fom_max_th2]->Write();
// gL_FOM_ac1[fom_max_th2]->SetName(Form("gL_FOM_ac1_%d",fom_max_th2));
// gL_FOM_ac1[fom_max_th2]->Write();
// 
// gL_ac2[fom_max_th1]->SetName(Form("gL_ac2_%d",fom_max_th1));  
// gL_ac2[fom_max_th1]->Write();
// gS_ac2[fom_max_th1]->SetName(Form("gS_ac2_%d",fom_max_th1)); 
// gS_ac2[fom_max_th1]->Write();  
// gL_FOM_ac2[fom_max_th1]->SetName(Form("gL_FOM_ac2_%d",fom_max_th1));
// gL_FOM_ac2[fom_max_th1]->Write();  
//
 //for(int i=0;i<3;i++){
// gSN_k_ac1[i][i]->SetFillColor(i+1);
// gSN_k_ac1[i][i]->SetMarkerColor(i+1);
// gSN_k_ac1[i][i]->SetFillStyle(3005);
// gSN_k_ac2[i][i]->SetFillColor(i+1);
// gSN_k_ac2[i][i]->SetMarkerColor(i+1);
// gSN_k_ac2[i][i]->SetFillStyle(3005);
// // TGraphErrors* gsum_k_ac1[100][100];
// gsum_k_ac1[i][i]->SetFillColor(i+1);
// gsum_k_ac1[i][i]->SetMarkerColor(i+1);
// gsum_k_ac1[i][i]->SetFillStyle(3005);
// gsum_k_ac2[i][i]->SetFillColor(i+1);
// gsum_k_ac2[i][i]->SetMarkerColor(i+1);
// gsum_k_ac2[i][i]->SetFillStyle(3005); 
// grate_k_ac1[i][i]->SetFillColor(i+1);
// grate_k_ac1[i][i]->SetMarkerColor(i+1);
// grate_k_ac1[i][i]->SetFillStyle(3005);
// grate_k_ac2[i][i]->SetFillColor(i+1);
// grate_k_ac2[i][i]->SetMarkerColor(i+1);
// grate_k_ac2[i][i]->SetFillStyle(3005);
//  //----- Pion -----------//
// grate_pi_ac1[i][i]->SetFillColor(i+1);
// grate_pi_ac1[i][i]->SetMarkerColor(i+1);
// grate_pi_ac1[i][i]->SetFillStyle(3005);
// grate_pi_ac2[i][i]->SetFillColor(i+1);
// grate_pi_ac2[i][i]->SetMarkerColor(i+1);
// grate_pi_ac2[i][i]->SetFillStyle(3005);  
//
//  //----- Proton -----------//
// grate_p_ac1[i][i]->SetFillColor(i+1);
// grate_p_ac1[i][i]->SetMarkerColor(i+1);
// grate_p_ac1[i][i]->SetFillStyle(3005);
// grate_p_ac2[i][i]->SetFillColor(i+1);
// grate_p_ac2[i][i]->SetMarkerColor(i+1);
// grate_p_ac2[i][i]->SetFillStyle(3005);  
// 
// 
// gSN_k_ac1[i][i]->SetName(Form("gSN_k_ac1_%d",i));
// gSN_k_ac1[i][i]->Write(); 
// gSN_k_ac2[i][i]->SetName(Form("gSN_k_ac2_%d",i)); 
// gSN_k_ac2[i][i]->Write(); 
// grate_k_ac1[i][i]->SetName(Form("grate_k_ac1_%d",i));
// grate_k_ac1[i][i]->Write(); 
// grate_k_ac2[i][i]->SetName(Form("grate_k_ac2_%d",i)); 
// grate_k_ac2[i][i]->Write();
// gsum_k_ac1[i][i]->SetName(Form("gsum_k_ac1_%d",i));
// gsum_k_ac1[i][i]->Write(); 
// gsum_k_ac2[i][i]->SetName(Form("gsum_k_ac2_%d",i)); 
// gsum_k_ac2[i][i]->Write();
//
// //---- proton ----//
// grate_p_ac1[i][i]->SetName(Form("grate_p_ac1_%d",i));
// grate_p_ac1[i][i]->Write(); 
// grate_p_ac2[i][i]->SetName(Form("grate_p_ac2_%d",i)); 
// grate_p_ac2[i][i]->Write(); 
//
// //---- pion ----//
// grate_pi_ac1[i][i]->SetName(Form("grate_pi_ac1_%d",i));
// grate_pi_ac1[i][i]->Write(); 
// grate_pi_ac2[i][i]->SetName(Form("grate_pi_ac2_%d",i)); 
// grate_pi_ac2[i][i]->Write(); 
// 
// 
// hcoin_t1[i]->Write();
// hcoin_t2[i]->Write();
// 
// }
//
//
// gSN_k_ac1[fom_max_th2][fom_max_th2]->SetName(Form("gSN_k_ac1_%d",fom_max_th2));
// gSN_k_ac1[fom_max_th2][fom_max_th2]->Write(); 
// gSN_k_ac2[fom_max_th1][fom_max_th1]->SetName(Form("gSN_k_ac2_%d",fom_max_th1)); 
// gSN_k_ac2[fom_max_th1][fom_max_th1]->Write(); 
// grate_k_ac1[fom_max_th2][fom_max_th2]->SetName(Form("grate_k_ac1_%d",fom_max_th2));
// grate_k_ac1[fom_max_th2][fom_max_th2]->Write(); 
// grate_k_ac2[fom_max_th1][fom_max_th1]->SetName(Form("grate_k_ac2_%d",fom_max_th1)); 
// grate_k_ac2[fom_max_th1][fom_max_th1]->Write(); 
// 
// facc_kc->Write();
// fk_kc->Write();
// fpi_pic->Write();
// fp_pc->Write();
// set->SetTH1(hmm_ac1_all_p[fom_th1][fom_max_th2],"hmm_ac1_all_p","Mass [GeV]","Counts/2 MeV");
// hmm_ac1_all_p[fom_th1][fom_max_th2]->Write();
// set->SetTH1(hmm_ac2_all_p[fom_th2][fom_max_th1],"hmm_ac2_all_p","Mass [GeV]","Counts/2 MeV"); 
// hmm_ac2_all_p[fom_th2][fom_max_th1]->Write();
// hmm->Write();
// hmm_acc->Write(); 
// hmm_p->Write();
// fL_all->Write();
// fL_p->Write();
// fS_all->Write();
// fS_p->Write(); 
// 
// hmm_fom->Write();
// hmm_fom_acc->Write();
// hmm_fom_p->Write();
// fL_fom->Write();
// fL_fom_p->Write();
// fS_fom->Write();
// fS_fom_p->Write(); 
//
// hmm_ac1_all_p[fom_th1][fom_max_th2]->Write();
// hmm_ac2_all_p[fom_th2][fom_max_th1]->Write();
//
// fLam[fom_th1][fom_max_th2][0]->Write();
// fSig[fom_th1][fom_max_th2][0]->Write();
// fLam_p->SetParameters(L0[fom_th1][fom_max_th2][0],L1[fom_th1][fom_max_th2][0],L2[fom_th1][fom_max_th2][0]);
// fSig_p->SetParameters(S0[fom_th1][fom_max_th2][0],S1[fom_th1][fom_max_th2][0],S2[fom_th1][fom_max_th2][0]);
// fLam_p->Write();
// fSig_p->Write(); 
// fk_fom->Write();
// 
 tnew->Write();
pEff1->Write();
pEff2->Write();
pEff3->Write();
pEff4->Write();
pEff5->Write();
pEff6->Write();
pEff7->Write();
pEff8->Write();
pEff9->Write();
h_pisr1->Write();
h_ksr1->Write();
h_psr1->Write();
h_pisr2l->Write();
h_ksr2l->Write();
h_psr2l->Write();
h_pisr2u->Write();
h_ksr2u->Write();
h_psr2u->Write();
h_pitot1->Write();
h_ktot1->Write();
h_ptot1->Write();
h_pitot2l->Write();
h_ktot2l->Write();
h_ptot2l->Write();
h_pitot2u->Write();
h_ktot2u->Write();
h_ptot2u->Write();
npe_sum_a1->Write();
npe_sum_a2->Write();
hcoin_k_fom_noAC1->Write();
hcoin_k_fom_noAC2->Write();
hcoin_bg_fom_noAC1->Write();
hcoin_bg_fom_noAC2->Write();
hcoin_wo_bg_fom_noAC1->Write();
hcoin_wo_bg_fom_noAC2->Write();
for(int i=0;i<nth;i++){
	for(int j=0;j<nth;j++){
		for(int l=0;l<nth;l++){
			hcoin_k_fom[i][j][l]->Write();
			hcoin_bg_fom[i][j][l]->Write();
			hcoin_wo_bg_fom[i][j][l]->Write();
		}
	}
 }



// hAC->Write();
// hAC2->Write(); 
// hcoin_tc->Write();
// hcoin_fom->Write();
// hRu1_time_s->Write(); 
// hRu1_time_c->Write();
//
// 
 //fnew->Close();
}



// #################################################
//=======================================================//
//================     Main       =======================//
//=======================================================//


int main(int argc, char** argv){

//  gStyle->SetOptFit(111111111);
//  int ch;
  //string ifname = "/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/scripts/ita_scripts/run_list/Lambda_test.list";
  //string ofname = "/pdf/hydro1_AC_eff_test.pdf";
 string ifname = "../small.list";//Run111157~111220
// string ifname = "../small2.list";//Run111157~111220 & Run111480~111542
//  string ifname = "../test.list";//for debug
//  string runlistname;
  string pname = "./Lambda_H1.param";
  string mtparam = "../matrix/matrix_new.list";
//  string root_init;
//  string root_end;
//  string pdf_init;
//  string pdf_end;
//  string ofname = "./test.pdf";
  string root_name = "./test.root";
  string print_name = "./test_print.pdf";
 // bool output_flag = false;
 // bool output_tree_flag = false;
 // bool draw_flag  = true;
 // bool coin_flag  = false;
  bool print_flag = true;
  bool root_flag  = false;
 // //  bool ac2_min=true;
 // string itel;  
 // string pngname;

  
//  extern char *optarg;
//  while((ch=getopt(argc,argv,"h:f:w:s:n:r:i:o:bcop:GHT12"))!=-1){
//    switch(ch){
//    case 'f':
//      ifname = optarg;
//      cout<<"input filename : "<<ifname<<endl;
//      break;
//    case 's':
//      ifname = optarg;
//      cout<<"input filename : "<<ifname<<endl;
//      root_flag = true;
//      draw_flag = false;            
//      break;
//
//      
//    case 'w':
//      print_flag = true;
//       draw_flag = false;
//      print_name = optarg;
//      cout<<"output PDF filename : "<<print_name<<endl;
//      break;
//
//
//    case 'r':
//      root_flag = true;
//      draw_flag = false;      
//      root_name = optarg;
//      cout<<"output root filename : "<<root_name<<endl;      
//      break;
//
//    case 'o':
//      root_flag=true;
//      draw_flag=false;
//      print_flag=true;
//      ofname = optarg;
//      root_name="./../rootfiles/ACtuning/" + ofname + ".root";
//      print_name="./../pdf/ACtuning/" +ofname + ".pdf";
//      break;
//      
//    case 'i':
//      itel= optarg;
//      nth= atoi(itel.c_str());
//      break;
//
//      
//    case 'b':
//      draw_flag = false;
//      cout<<"BACH MODE!"<<endl;
//      break;
//  
//    case 'G':
//    mode="G";
//      break;
//  
//    case 'H':
//    mode="H";
//      break;
//
//    case 'T':
//      mode="T";    
//	break;
//
//    case 'U':
//      ac2_min=false;    
//	break;
//
//  case '1':
//    tdc_time=56.23e-3;//[ns]
//    kine=1;
//      break;
//
//  case '2':
//    tdc_time=58e-3;//[ns]
//    kine=2;
//      break;
//
//
//    case 'h':
//      cout<<"-f : input root filename"<<endl;
//      cout<<"-w : output pdf filename"<<endl;
//      cout<<"-r : output root filename"<<endl;      
//      cout<<"-o : output pdf & root  filename"<<endl;
//      cout<<"-1 : F1TDC resolution 56 ns"<<endl;
//      cout<<"-2 : F1TDC resolution 58 ns"<<endl;      
//      cout<<"-T or H or G : Mode of root"<<endl;      
//      return 0;
//      break;
//    case '?':
//      cout<<"unknown option...."<<endl;
//      return 0;
//      break;
//    default:
//      cout<<"type -h to see help!!"<<endl;
//      return 0;
//    }
//  }

//mode = "H";
//kine = 1;
    //tdc_time=56.23e-3;//[ns]

  
  tuning* AC=new tuning();
//  TApplication theApp("App", &argc, argv);
cout << "Start SetRunList" << endl;
  AC->matrix(mtparam);
  if(root_flag)AC->SetRoot(ifname);//root_name?
//cout << "Start SetBranch" << endl;
//  AC->SetBranch();
  AC->GetACParam();
cout << "Start SetParam" << endl;
  AC->SetParam();
cout << "Start MakeHist" << endl;
  AC->MakeHist();
  AC->ReadParam(pname);
cout << "Start SetParam" << endl;
  AC->SetRunList(ifname);
cout << "Start Fill" << endl;
  AC->Filling();
	cout << "AC->Fill() is done" << endl;
	cout << "fabs(-1)= " << fabs(-1) << endl;
  AC->Fitting();
	cout << "AC->Fitting() is done" << endl;
  AC->ACtune();
	cout << "AC->ACtune() is done" << endl;
  AC->Draw();
  if(print_flag)AC->Print(print_name);
  if(root_flag)AC->Write();
//  AC->Comment();
  
  cout<<"=========== Output files ============="<<endl;
//  cout<<"output rootfile: "<<root_name<<endl;
  
  
	
//  TApplication *theApp =new TApplication("App",&argc,argv);
// if(draw_flag==0)gROOT->SetBatch(1);
// if(draw_flag==0)gSystem->Exit(1);
// theApp->Run();
 return 0;


}
