#include <iostream>
#include <fstream>
using namespace std;
#include "TApplication.h"
#include "h2all.h"
#include "Param.h"
#include "Tree.h"
#include "TMath.h"
#include "TROOT.h"//FitSlice, gROOT->FindObject

#define F1TDC
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

}//s2f1_off()



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
  
}//calcf2t_mom()



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
  
}//calcf2t_ang()
      	


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
}//calcf2t_zt



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

  }//GetACParam()



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
}//AC_npe



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
}//calcf_pathl


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
  
}//calcf2t_3rd()

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
  
}//matrix()

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
  
}//MTParam_R()
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

}//MTParam_L()

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
  
}//MTP_mom()

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
//  sprintf(name_MctimeL,"../ctimeL.dat"); 
//  sprintf(name_MctimeR,"../ctimeR.dat"); 
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


}//MTParam_G()

///////////////////////////////////////////////////////////////////////////

void tuning::Calib(int rt, int lt ){

  // ==== Initialization ======//
  R_p = 0.0;
  L_p = 0.0;



  //======= Nomalization ==================//
  R_tr_x[rt]    = (R_tr_x[rt]-XFPm)/XFPr;
  R_tr_th[rt]   = (R_tr_th[rt]-XpFPm)/XpFPr;
  R_tr_y[rt]    = (R_tr_y[rt]-YFPm)/YFPr;
  R_tr_ph[rt]   = (R_tr_ph[rt]-YpFPm)/YpFPr;
  R_tr_vz[rt]   = (R_tr_vz[rt]-Ztm)/Ztr;
  R_tr_tg_th[rt]= (R_tr_tg_th[rt] - Xptm)/Xptr;
  R_tr_tg_ph[rt]= (R_tr_tg_ph[rt] - Yptm)/Yptr;

 //  R_p = (R_p - PRm)/PRr;
  R_p = (R_tr_p[rt] - PRm)/PRr;

  L_tr_x[lt]    = (L_tr_x[lt]-XFPm)/XFPr; 
  L_tr_th[lt]   = (L_tr_th[lt]-XpFPm)/XpFPr;
  L_tr_y[lt]    = (L_tr_y[lt]-YFPm)/YFPr;
  L_tr_vz[lt]   = (L_tr_vz[lt]-Ztm)/Ztr;
  L_tr_ph[lt]   = (L_tr_ph[lt]-YpFPm)/YpFPr;
  L_tr_tg_th[lt]= (L_tr_tg_th[lt] - Xptm)/Xptr;
  L_tr_tg_ph[lt]= (L_tr_tg_ph[lt] - Yptm)/Yptr;  

  //  L_p = (L_p - PLm)/PLr;
  L_p = (L_tr_p[lt] - PLm)/PLr;

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
    B_p     = B_p - Eloss(0.0,0,"B");
    R_p     = R_p + Eloss(R_tr_tg_ph[rt],R_tr_vz[rt],"R");
    L_p     = L_p + Eloss(L_tr_tg_ph[lt],L_tr_vz[lt],"L");


    
}//Calib()
////////////////////////////////////////////////////////////////

//Constructer
tuning::tuning(){
  set= new Setting();
  set->Initialize();
}//tuning()

//Destructer
tuning::~tuning(){}

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
}//SetRunList()

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
}//ReadParam()

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
  ct_test=ct;
  tr.ct_itabashi=ct;
  //tr.ct_g=CoinCalc_gogami(RS2_seg,LS2_seg,rhit,lhit);
  ct=CoinCalc_gogami(RS2_seg,LS2_seg,rhit,lhit);
 // ct_test=CoinCalc_c(RS2_seg,LS2_seg,rhit,lhit);

//Changed
  
}//CoinCalc()


////////////////////////////////////////////////////////////////////////////
double tuning::CoinCalc_gogami(int RS2_seg, int LS2_seg,int rhit, int lhit){

//okuyama (2020/11/19)
  int dataflag = 1;//H_1
  //int dataflag = 2;//H_2
  double kcenter = 3.122;

//beta
  double beta_R  = R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+Mpi*Mpi);
  //double beta_R  = R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+MK*MK);
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



 //double ctime = - meantime_R + ctimecorL + ctimecorR;
 double ctime = - meantime_R;//okuyama (2020/11/19)
 double ctime_before = - meantime_R;


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
 ctime_before=ctime_before -1.4;


 tr.ct_gb=-ctime;
double pion_pos = 3.1;//ns
//K. Okuyama (2020.11.19)
//Cointime collection from the out of range was performed
//by adjusting the pion position.

if(dataflag==1){//H2-1
	if(fabs(ctime - 3662.43)<25.){
 	   ctime = ctime - 3662.43 - pion_pos;
 	}
}else{//H2-2
	if(fabs(ctime - 3662.39)<25.){
 	   ctime = ctime - 3662.39 - pion_pos;
 	}
	if(fabs(ctime - 2984.54)<25.){
 	   ctime = ctime - 2984.54 - pion_pos;
 	}
	if(fabs(ctime - 2520.03)<25.){
 	   ctime = ctime - 2520.03 - pion_pos;
 	}
	if(fabs(ctime - 2181.05)<25.){
 	   ctime = ctime - 2181.05 - pion_pos;
 	}
	if(fabs(ctime - 2017.86)<25.){
 	   ctime = ctime - 2017.86 - pion_pos;
 	}
	if(fabs(ctime - 1151.75)<25.){
 	   ctime = ctime - 1151.75 - pion_pos;
 	}
	if(fabs(ctime - 687.25)<25.){
 	   ctime = ctime - 687.25 - pion_pos;
 	}
	if(fabs(ctime - 348.28)<25.){
 	   ctime = ctime - 348.28 - pion_pos;
 	}
	if(fabs(ctime - 185.09)<25.){
 	   ctime = ctime - 185.09 - pion_pos;
 	}
	//if(fabs(ctime + 681.04)<25.){
 	//   ctime = ctime + 681.04 - pion_pos;
 	//}
	//if(fabs(ctime + 1145.54)<25.){
 	//   ctime = ctime + 1145.54 - pion_pos;
 	//}
}

/////////////////////////	
// OLD COINTIME OFFSET //
/////////////////////////	
// if(3600.0<ctime && ctime<3665){
//   ctime = ctime - 3637.88 - 12.76;
//   ctime = ctime - 12.0-3.1;
// }
////else ctime=-9999.;
// if(3600.0<ctime_before && ctime_before<3665){
//   ctime_before = ctime_before - 3637.88 - 12.76;
//   ctime_before = ctime_before - 12.0-3.1;
// }
//
//
//	    // ----- ctime offset ------ 
//	    if(dataflag==1){
//	      if(fabs(ctime-3650.67)<25.0){
//	    	ctime = ctime - 3650.67;
//	      }
//	    }
//	    else{
//	      if(fabs(ctime+1469.0)<25.0){
//	    	ctime = ctime +1469.0;
//	      }
//	      else if (fabs(ctime-348.9)<25.0){
//	    	ctime = ctime -348.9;
//	      }
//	      else if (fabs(ctime-3638.00)<25.0){
//	    	ctime = ctime -3638.00;
//	      }
//	    }
//
 ctime=-ctime;
 ctime_before=-ctime_before;
tr.ct_orig=ctime;

 return ctime;


}//CoinCalc_gogami

///////////////////////////////////////////////////////////////////////////
double tuning::CoinCalc_c(int RS2_seg, int LS2_seg, int rhit, int lhit){

  double cointime=0.0;

  
  convertF1TDCR(param);
  convertF1TDCL(param);

  //  double Rpathl=R_tr_pathl[rhit]+R_s2_trpath[rhit];
  //  double Lpathl=L_tr_pathl[lhit]+L_s2_trpath[lhit];

  double Beta_R=R_p/sqrt(R_p*R_p+MK*MK);
  double Beta_L=L_p/sqrt(L_p*L_p+Me*Me);
  double tof_r=RS2_F1time[RS2_seg] - R_pathl/(Beta_R*LightVelocity);
  double tof_l=LS2_F1time[LS2_seg] - L_pathl/(Beta_L*LightVelocity);
  double tof_rc=RS2_F1time_c[RS2_seg] - R_pathl/(Beta_R*LightVelocity);
  double tof_lc=LS2_F1time_c[LS2_seg] - L_pathl/(Beta_L*LightVelocity);

    
  if(RS2_F1time[RS2_seg]!=-9999. &&LS2_F1time[LS2_seg]!=-9999.){
    cointime= - tof_r + tof_l - coin_offset;
    //    tr.ct_c= - (rtof[RS2_seg] - R_pathl/(Beta_R*LightVelocity))
    //      + (ltof[LS2_seg] - L_pathl/(Beta_L*LightVelocity))        - coin_offset;
    tr.ct_c= - tof_rc + tof_lc -coin_offset;
  }
  else{
    cointime=-1000;
    tr.ct_c =-1000;
  }

  
  if(tof_r!=tof_rc)
    cout<<" ct "<<cointime<<" ct_c "<<tr.ct_c<<endl;
  
  return cointime;
  
}//CoinCalc_c

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



}//PathCalib()


//////////////////////////////////////////////////////////////////

double tuning::Eloss(double yp,double z,char const* arm){

  double hrs_ang=13.2*3.14159/180.;  
  double x;
  

  //----- Original coordinate  -------//
  // Definition by K.N. Suzuki  (fixed Oct. 23rd, 2019)//
  // R-HRS : right hand coordinate (Unticlockwise rotation)//
  // L-HRS : left  hand coordinate (    Clockwise rotation)//
  
  if(*arm=='R')        x = - hrs_ang + yp; //yp : phi [rad] RHRS
  else if(*arm=='L')   x = - hrs_ang - yp; //yp : phi [rad] LHRS
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
  }else{ph[0]=0.;ph[1]=0.;ph[2]=0.;pl[0]=0.;pl[1]=0.;}//Okuyama

 
  if(high){
    dEloss_h = ph[0]*sin(ph[1]*x)+ph[2];    
    dEloss = dEloss_h;
  }else{
    dEloss_l = pl[0]*x +pl[1];    
    dEloss = dEloss_l;}
  //==== thickness 0.4 mm in beam energy loss ======//
  if(*arm=='B')dEloss=0.184; //[MeV/c]

if(z<-0.1){
    if(*arm=='B')dEloss  = 0.1345;
    if(*arm=='R')dEloss += 0.803;
    if(*arm=='L')dEloss += 0.809;
  }
  else if(z>0.1){
    if(*arm=='B')dEloss  = 0.264;
    if(*arm=='R')dEloss += 0.313;
    if(*arm=='L')dEloss += 0.313;
    
  }

  dEloss=dEloss/1000.; // [GeV/c]
  return dEloss;

  
}//Eloss()

void tuning::SetParam(){

	mt = Mp;//target mass
	mh = ML;//hypernuclei
cout << "mt (target) = " << mt << endl;
cout << "mh (hyper) = " << mh << endl;
     //----------------Histogram parameters---------------//
	 min_Lp = 1.8;
	 max_Lp = 2.8;
     min_coin_c=-20.0;
     max_coin_c=20.0;
     min_mm=-0.1;//GeV/c^2
     max_mm=0.2;//GeV/c^2
     bin_mm=(max_mm-min_mm)/0.001; //Counts/2 MeV
     bin_mm=(int)bin_mm;
	 iter_ac1=30;//iteration number
     min_s2=-10;
     max_s2=5000;
     bin_s2=max_s2-min_s2;
     bin_coin_c=(int)((max_coin_c-min_coin_c)/0.056);
cout<<"tdc"<<tdc_time<<endl;
cout<<"max coin"<<max_coin_c<<endl;
cout<<"min coin"<<min_coin_c<<endl;
cout<<"bin coin"<<bin_coin_c<<endl;

     //----------------ana_Lambda's Notation---------------//
	 th1_max=2.0;
     th2_max=20.0;
     th2_min=0.0;

     //----------------Cut parameters---------------//
	 zver[0]=0.;
	 for(int i=1; i<nth; i++) zver[i]=zver[i-1]+0.005;//SUM
	 ac1_adc[0]=0.0;
	 for(int i=1; i<nth; i++) ac1_adc[i]=ac1_adc[i-1]+0.25;
     ac2l_adc[0]=0.0;//best = 3.0? ac2l[30]
	 for(int i=1; i<nth; i++) ac2l_adc[i]=ac2l_adc[i-1]+0.2;
     ac2u_adc[0]=0.0;//best = 18? ac2u[16]
	 for(int i=1; i<nth; i++) ac2u_adc[i]=ac2u_adc[i-1]+0.5;
	 //---Kaon Cut ----//
 	 ac1_kcut=100.;
 	 ac2_kcut_min=1000.;
 	 ac2_kcut_max=5000.;


}//SetParam()
////////////////////////////////////////////////////////////

void tuning::MakeHist(){
  cout<<"Make Hist "<<endl;
	file_out = new TFile("h2all_1_temp.root","recreate");
	//file_out = new TFile("h2all_1_2020Nov.root","recreate");
	tree_out = new TTree("tree_out","tree_out");
	//`tree_out ->Branch("branch name",variable ,"branch name/type");
	
	tree_out ->Branch("pid_cut"        ,&tr.pid_cut      ,"pid_cut/I"     );
	tree_out ->Branch("ct_cut"        ,&tr.ct_cut      ,"ct_cut/I"     );
	tree_out ->Branch("z_cut"        ,&tr.z_cut      ,"z_cut/I"     );
    tree_out ->Branch("nrun"        ,&tr.nrun      ,"nrun/I"     );
    tree_out ->Branch("nev"        ,&tr.nev      ,"nev/I"     );
    tree_out ->Branch("runnum",&runnum ,"runnum/I");
	//  tree_out ->Branch("ntr_r",&tr.ntrack_r ,"ntr_r/I");
	//  tree_out ->Branch("ntr_l",&tr.ntrack_l ,"ntr_l/I");
	tree_out ->Branch("mm",&tr.missing_mass ,"missing_mass/D");
	tree_out ->Branch("mm_b",&tr.missing_mass_b ,"missing_mass_b/D");
	tree_out ->Branch("mm_L",&tr.missing_mass_L ,"missing_mass_L/D");
	tree_out ->Branch("mm_nnL",&tr.missing_mass_nnL ,"missing_mass_nnL/D");
	tree_out ->Branch("mm_H3L",&tr.missing_mass_H3L ,"missing_mass_H3L/D");
	tree_out ->Branch("mm_cut",&tr.missing_mass_cut ,"missing_mass_cut/D");
	tree_out ->Branch("mm_MgL",&tr.missing_mass_MgL ,"missing_mass_MgL/D");
	tree_out ->Branch("mm_MgL_acc",&tr.missing_mass_MgL_acc ,"missing_mass_MgL_acc/D");
	tree_out ->Branch("mm_acc",&tr.missing_mass_acc ,"missing_mass_acc/D");
	tree_out ->Branch("ct_b",&tr.ct_b ,"ct_b/D");
	tree_out ->Branch("ct_c",&tr.ct_c ,"ct_c/D");
	tree_out ->Branch("ct_g",&tr.ct_g ,"ct_g/D");
	tree_out ->Branch("ct_gb",&tr.ct_gb ,"ct_gb/D");
	tree_out ->Branch("yp_cor",&tr.yp_cor ,"yp_cor/D");
	tree_out ->Branch("ctimecorR",&tr.ctimecorR ,"ctimecorR/D");
	tree_out ->Branch("ctimecorL",&tr.ctimecorL ,"ctimecorL/D");
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
	tree_out ->Branch("momR"        ,&tr.momR      ,"momR/D"     );
	tree_out ->Branch("momL"        ,&tr.momL      ,"momL/D"     );
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
	tree_out->Branch("Bp_c"     ,&tr.Bp_c   ,"Bp_c/D"  );
	tree_out->Branch("Lp_c"     ,&tr.Lp_c   ,"Lp_c/D"  );
	tree_out->Branch("Rp_c"     ,&tr.Rp_c   ,"Rp_c/D"  );
	tree_out->Branch("Bp_before"     ,&tr.Bp_before   ,"Bp_before/D"  );
	tree_out->Branch("Lp_before"     ,&tr.Lp_before   ,"Lp_before/D"  );
	tree_out->Branch("Rp_before"     ,&tr.Rp_before   ,"Rp_before/D"  );
	tree_out ->Branch("trig"     ,&tr.trig   ,"trig/D"  );
	tree_out->Branch("dpe"     ,&tr.dpe   ,"dpe/D"  );
	tree_out->Branch("dpe_"     ,tr.dpe_   ,"dpe_[10]/D"  );
	tree_out->Branch("dpk"     ,tr.dpk   ,"dpk[10]/D"  );
	tree_out->Branch("ct_orig"     ,&tr.ct_orig   ,"ct_orig/D"  );
	tree_out->Branch("ct_itabashi"     ,&tr.ct_itabashi   ,"ct_itabashi/D"  );
	tree_out->Branch("tr.ntrack_l"  ,&tr.ntrack_l   ,"tr.ntrack_l/I"  );
	tree_out->Branch("tr.ntrack_r"  ,&tr.ntrack_r   ,"tr.ntrack_r/I"  );

//ADD 2020/8/12
    tree_out->Branch("L.tr.chi2"  ,&tr.chi2_l   ,"L.tr.chi2/D"  );
    tree_out->Branch("L.tr.x"  ,&tr.x_l   ,"L.tr.x/D"  );
    tree_out->Branch("L.tr.y"  ,&tr.y_l   ,"L.tr.y/D"  );
    tree_out->Branch("L.tr.th"  ,&tr.th_l   ,"L.tr.th/D"  );
    tree_out->Branch("L.tr.ph"  ,&tr.ph_l   ,"L.tr.ph/D"  );
    tree_out->Branch("L.tr.p"  ,&tr.mom_l   ,"L.tr.p/D"  );
    tree_out->Branch("L.tr.tg_th"  ,&tr.tg_th_l   ,"L.tr.tg_th/D"  );
    tree_out->Branch("L.tr.tg_ph"  ,&tr.tg_ph_l   ,"L.tr.tg_ph/D"  );
    tree_out->Branch("L.tr.vz"  ,&tr.vz_l   ,"L.tr.vz/D"  );

    tree_out->Branch("R.tr.chi2"  ,&tr.chi2_r   ,"R.tr.chi2/D"  );
    tree_out->Branch("R.tr.x"  ,&tr.x_r   ,"R.tr.x/D"  );
    tree_out->Branch("R.tr.y"  ,&tr.y_r   ,"R.tr.y/D"  );
    tree_out->Branch("R.tr.th"  ,&tr.th_r   ,"R.tr.th/D"  );
    tree_out->Branch("R.tr.ph"  ,&tr.ph_r   ,"R.tr.ph/D"  );
    tree_out->Branch("R.tr.p"  ,&tr.mom_r   ,"R.tr.p/D"  );
    tree_out->Branch("R.tr.tg_th"  ,&tr.tg_th_r   ,"R.tr.tg_th/D"  );
    tree_out->Branch("R.tr.tg_ph"  ,&tr.tg_ph_r   ,"R.tr.tg_ph/D"  );
    tree_out->Branch("R.tr.vz"  ,&tr.vz_r   ,"R.tr.vz/D"  );
//ADD 2020/11/21
    tree_out->Branch("R.tr.pathl"  ,&tr.pathl_r   ,"R.tr.pathl/D"  );
    tree_out->Branch("L.tr.pathl"  ,&tr.pathl_l   ,"L.tr.pathl/D"  );
    tree_out->Branch("R.s2.trpath"  ,&tr.trpath_r   ,"R.s2.trpath/D"  );
    tree_out->Branch("L.s2.trpath"  ,&tr.trpath_l   ,"L.s2.trpath/D"  );
    tree_out->Branch("theta_gk_cm"  ,&tr.theta_gk_cm   ,"theta_gk_cm/D"  );
    tree_out->Branch("Qsq"  ,&tr.Qsq   ,"Qsq/D"  );
    tree_out->Branch("W"  ,&tr.W   ,"W/D"  );
    tree_out->Branch("eps"  ,&tr.eps   ,"eps/D"  );


/////////////
//// Z  /////
/////////////

	npe_sum_a1 = new TH1F("npe_sum_a1","NPE SUM A1",2000,0.,40.);
//	  npe_sum_a2 = new TH1F("npe_sum_a2","NPE SUM A2",4000,0.,80.);
	npe_sum_a2 = new TH1F("npe_sum_a2","NPE SUM A2 (Kaon)",300,0.,30.);
////--------------------------------------------//
  
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
cout<<"LHRS Set TH"<<endl;
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
    h_L_vz2       = new TH1D("h_L_vz2"      ,"h_L_vz2"      ,400,-0.25,0.25);
    h_L_y_x       = new TH2D("h_L_y_x"      ,"h_L_y_x"      ,200,   -1,  1 ,200,-0.1,0.1);
    h_L_th_x      = new TH2D("h_L_th_x"     ,"h_L_th_x"     ,200,   -1,  1 ,200,-0.2,0.2);
    h_L_ph_y      = new TH2D("h_L_ph_y"     ,"h_L_ph_y"     ,200, -0.1, 0.1,200,-0.1,0.1);
    h_L_tgph_tgth = new TH2D("h_L_tgph_tgth","L_th : L_ph (w/ theta_ee Cut)",200, -0.1, 0.1,200,-0.06,0.06);
//    h_L_tgph_tgth2 = new TH2D("h_L_tgph_tgth2","L_th : L_ph (w/ Z, AC, Kaon(ct) Cut)",200, -0.1, 0.1,200,-0.06,0.06);
    h_L_tgph_tgth2 = new TH2D("h_L_tgph_tgth2","L_th : L_ph (w/ VP Flux Cut)",200, -0.1, 0.1,200,-0.06,0.06);
    h_L_tgph_tgth3 = new TH2D("h_L_tgph_tgth3","L_th : L_ph (w/ Z Cut)",200, -0.1, 0.1,200,-0.06,0.06);
    h_L_tgph_tgth4 = new TH2D("h_L_tgph_tgth4","L_th : L_ph (w/ AC Cut)",200, -0.1, 0.1,200,-0.06,0.06);
    h_L_tgph_tgth5 = new TH2D("h_L_tgph_tgth5","L_th : L_ph (w/ Z, AC Cut)",200, -0.1, 0.1,200,-0.06,0.06);
    h_L_tgph_tgth6 = new TH2D("h_L_tgph_tgth6","L_th : L_ph (w/ ZZ_zone1 Cut)",200, -0.1, 0.1,200,-0.06,0.06);
    h_L_tgph_tgth7 = new TH2D("h_L_tgph_tgth7","L_th : L_ph (w/ ZZ_zone2 Cut)",200, -0.1, 0.1,200,-0.06,0.06);
    h_L_tgph_tgth8 = new TH2D("h_L_tgph_tgth8","L_th : L_ph (w/ ZZ_zone3 Cut)",200, -0.1, 0.1,200,-0.06,0.06);
    h_L_tgph_tgth9 = new TH2D("h_L_tgph_tgth9","L_th : L_ph (w/ ZZ_zone4 Cut)",200, -0.1, 0.1,200,-0.06,0.06);
    h_L_tgph_tgth10 = new TH2D("h_L_tgph_tgth10","L_th : L_ph (Al front)",200, -0.1, 0.1,200,-0.06,0.06);
    h_L_tgph_tgth11 = new TH2D("h_L_tgph_tgth11","L_th : L_ph (Al back)",200, -0.1, 0.1,200,-0.06,0.06);
    h_L_tgph_tgth12 = new TH2D("h_L_tgph_tgth12","L_th : L_ph (No Cut)",200, -0.1, 0.1,200,-0.06,0.06);
    h_LHRS_phth = new TH2D("h_LHRS_phth","L_th : L_ph (LHRS frame), No Cut",200, 0., 0.1,200,0.,2*PI);
    h_LHRS_phth2 = new TH2D("h_LHRS_phth2","L_th : L_ph (LHRS frame), w/ Z, AC",200, 0., 0.1,200,0.,2*PI);
    h_original_phth = new TH2D("h_original_phth","L_th : L_ph (original frame), No Cut",200, 0.1, 0.35,200,1.,2.);
    h_original_phcosth = new TH2D("h_original_phcosth","L_th : L_ph (original frame), 2.1<pL<2.15",200, 0.95, 1.,200,1.,2.);
    h_original_phcosth2 = new TH2D("h_original_phcosth2","L_th : L_ph (original frame), pL>2.16",200, 0.95, 1.,200,1.,2.);
    h_original_phth2 = new TH2D("h_original_phth2","L_th : L_ph (original frame), w/ Z, AC Cut",200, 0.1, 0.35,200,1.,2.);
    h_original_phth3 = new TH2D("h_original_phth3","L_th : L_ph (original frame), top quality",200, 0.1, 0.35,200,1.,2.);
    h_original_phth4 = new TH2D("h_original_phth4","L_th : L_ph (original frame), w/ Z Cut",200, 0.1, 0.35,200,1.,2.);
    h_original_phth5 = new TH2D("h_original_phth5","L_th : L_ph (original frame), w/ Z Cut, p>2.18",200, 0.1, 0.35,200,1.,2.);
    h_original_phth6 = new TH2D("h_original_phth6","L_th : L_ph (original frame), w/ Z Cut, 2.1<p<2.16",200, 0.1, 0.35,200,1.,2.);
    h_original_phth7 = new TH2D("h_original_phth7","L_th : L_ph (original frame), w/ Z Cut, 2.10<p<2.12",200, 0.1, 0.35,200,1.,2.);
    h_original_phth8 = new TH2D("h_original_phth8","L_th : L_ph (original frame), w/ Z Cut, 2.12<p<2.14",200, 0.1, 0.35,200,1.,2.);
    h_original_phth9 = new TH2D("h_original_phth9","L_th : L_ph (original frame), w/ Z Cut, 2.18<p<2.2",200, 0.1, 0.35,200,1.,2.);
    h_original_phth10 = new TH2D("h_original_phth10","L_th : L_ph (original frame), w/ Z Cut, 2.2<p<2.22",200, 0.1, 0.35,200,1.,2.);
    h_L_z_x       = new TH2D("h_L_z_x"      ,"h_L_z_x"      ,200,   -1.,1.,  200,-0.25,0.25);
    h_L_z_y       = new TH2D("h_L_z_y"      ,"h_L_z_y"      ,200,   -0.1,0.1,  200,-0.25,0.25);
    h_L_th_y      = new TH2D("h_L_th_y"     ,"h_L_th_y"     ,200,   -0.1,  0.1 ,200,-0.1,0.1);
    h_L_ph_x      = new TH2D("h_L_ph_x"     ,"h_L_ph_x"     ,200, -1., 1.,200,-0.1,0.1);
    h_L_th_z      = new TH2D("h_L_th_z"      ,"h_L_th_z"      ,200,   -0.25,  0.25 ,200,-0.1,0.1);
    h_L_ph_z      = new TH2D("h_L_ph_z"      ,"h_L_ph_z"      ,200,   -0.25,  0.25 ,200,-0.1,0.1);
    h_L_tgth_x      = new TH2D("h_L_tgth_x"     ,"h_L_tgth_x"     ,200,   -1,  1 ,200,-0.1,0.1);
    h_L_tgph_y      = new TH2D("h_L_tgph_y"     ,"h_L_tgph_y"     ,200, -0.1, 0.1,200,-0.1,0.1);
    h_L_tgth_y      = new TH2D("h_L_tgth_y"     ,"h_L_tgth_y"     ,200,   -0.1,  0.1 ,200,-0.1,0.1);
    h_L_tgph_x      = new TH2D("h_L_tgph_x"     ,"h_L_tgph_x"     ,200, -1., 1.,200,-0.1,0.1);
    h_L_tgth_z      = new TH2D("h_L_tgth_z"      ,"h_L_tgth_z"      ,200,   -0.25,  0.25 ,200,-0.1,0.1);
    h_L_tgph_z      = new TH2D("h_L_tgph_z"      ,"h_L_tgph_z"      ,200,   -0.25,  0.25 ,200,-0.1,0.1);
    h_L_corR_tgth   = new TH2D("h_L_corR_tgth",	"h_L_corR_tgth",	200, -0.1, 0.1, 200, -5., 5.);	
    h_L_corR_tgph   = new TH2D("h_L_corR_tgph",	"h_L_corR_tgph",	200, -0.1, 0.1, 200, -5., 5.);	
    h_L_corR_th   = new TH2D("h_L_corR_th",	"h_L_corR_th",	200, -0.1, 0.1, 200, -5., 5.);	
    h_L_corR_ph   = new TH2D("h_L_corR_ph",	"h_L_corR_ph",	200, -0.1, 0.1, 200, -5., 5.);	
    h_L_corR_x   = new TH2D("h_L_corR_x",	"h_L_corR_x",	200, -1., 1., 200, -5., 5.);	
    h_L_corR_y   = new TH2D("h_L_corR_y",	"h_L_corR_y",	200, -0.1, 0.1, 200, -5., 5.);	
    h_L_corR_z   = new TH2D("h_L_corR_z",	"h_L_corR_z",	200, -0.25, 0.25, 200, -5., 5.);	
    h_L_corR_zdiff   = new TH2D("h_L_corR_zdiff",	"h_L_corR_zdiff",	200, -0.25, 0.25, 200, -5., 5.);	
    h_L_corL_tgth   = new TH2D("h_L_corL_tgth",	"h_L_corL_tgth",	200, -0.1, 0.1, 200, -5., 5.);	
    h_L_corL_tgph   = new TH2D("h_L_corL_tgph",	"h_L_corL_tgph",	200, -0.1, 0.1, 200, -5., 5.);	
    h_L_corL_th   = new TH2D("h_L_corL_th",	"h_L_corL_th",	200, -0.1, 0.1, 200, -5., 5.);	
    h_L_corL_ph   = new TH2D("h_L_corL_ph",	"h_L_corL_ph",	200, -0.1, 0.1, 200, -5., 5.);	
    h_L_corL_x   = new TH2D("h_L_corL_x",	"h_L_corL_x",	200, -1., 1., 200, -5., 5.);	
    h_L_corL_y   = new TH2D("h_L_corL_y",	"h_L_corL_y",	200, -0.1, 0.1, 200, -5., 5.);	
    h_L_corL_z   = new TH2D("h_L_corL_z",	"h_L_corL_z",	200, -0.25, 0.25, 200, -5., 5.);	
    h_L_corL_zdiff   = new TH2D("h_L_corL_zdiff",	"h_L_corL_zdiff",	200, -0.25, 0.25, 200, -5., 5.);	
    h_L_ct_zdiff   = new TH2D("h_L_ct_zdiff",	"h_L_ct_zdiff",	200, -0.25, 0.25, bin_coin_c, min_coin_c, max_coin_c);	
    h_L_corsum_zdiff   = new TH2D("h_L_corsum_zdiff",	"h_L_corsum_zdiff",	200, -0.25, 0.25, 200, -5., 5.);	
    h_L_cordiff_zdiff   = new TH2D("h_L_cordiff_zdiff",	"h_L_cordiff_zdiff",	200, -0.25, 0.25, 200, -5., 5.);	
    h_L_corL_corR   = new TH2D("h_L_corL_corR",	"h_L_corL_corR",	200, -5., 5., 200, -5., 5.);	

    
 
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
    set->SetTH2(h_L_z_x      ,"Focal Plane Z v.s X"     ,"X (m)"           ,"Z (m)");
    set->SetTH2(h_L_z_y      ,"Focal Plane Z v.s Y"     ,"Y (m)"           ,"Z (m)");
    set->SetTH2(h_L_th_x     ,"Focal Plane #theta v.s X","X (m)"           ,"#theta (rad)");
    set->SetTH2(h_L_ph_y     ,"Focal Plane #phi v.s Y"  ,"Y (m)"           ,"#phi (rad)");
    set->SetTH2(h_L_th_y     ,"Focal Plane #theta v.s Y","Y (m)"           ,"#theta (rad)");
    set->SetTH2(h_L_ph_x     ,"Focal Plane #phi v.s X"  ,"X (m)"           ,"#phi (rad)");
    set->SetTH2(h_L_th_z     ,"Focal Plane #theta v.s Z","Z (m)"           ,"#theta (rad)");
    set->SetTH2(h_L_ph_z     ,"Focal Plane #phi v.s Z"  ,"Z (m)"           ,"#phi (rad)");
    set->SetTH2(h_L_tgph_tgth,"Target #phi v.s #theta (w/ theta_ee Cut)"  ,"#theta_{t} (rad)","#phi_{t} (rad)");

    h_L_beta        = new TH1D("h_L_beta"       ,"h_L_beta"       ,400,   0,  2); 
    //h_L_m2          = new TH1D("h_L_m2"         ,"h_L_m2"         ,400,-0.5,2.5); 
    h_L_m2          = new TH1D("h_L_m2"         ,"h_L_m2"         ,400,-0.5,1.5); 
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
cout<<"RHRS Set TH"<<endl;
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
    h_R_vz2       = new TH1D("h_R_vz2"      ,"h_R_vz2"      ,400,-0.25,0.25);
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
cout<<"Coin Set TH"<<endl;
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
    h_Rz2     = new TH1D("h_Rz2", "h_Rz2",1000,-0.2,0.2);
    h_Rz_c   = new TH1D("h_Rz_c", "h_Rz_c",1000,-0.2,0.2);
    h_Rz_cut   = new TH1D("h_Rz_cut", "h_Rz_cut",1000,-0.2,0.2);
    h_Lz     = new TH1D("h_Lz", "h_Lz",1000,-0.2,0.2);
    h_Lz2     = new TH1D("h_Lz2", "h_Lz2",1000,-0.2,0.2);
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
    h_Lp_top   = new TH1D("h_Lp_top", "mom_L (w/ Lambda Cut)",1000,1.9,2.3);
    
    h_theta_ee = new TH1D("h_theta_ee", "theta_ee",1000,0.1,0.35);
    h_theta_ee2 = new TH1D("h_theta_ee2", "theta_ee (w/ Z_Diff Cut)",1000,0.1,0.35);
    h_theta_ee3 = new TH1D("h_theta_ee3", "theta_ee (w/ Z Cut)",1000,0.1,0.35);
    h_theta_ee4 = new TH1D("h_theta_ee4", "theta_ee (w/ Z, AC Cut)",1000,0.1,0.35);
    h_phi_ee = new TH1D("h_phi_ee", "phi_ee",1000,0.,2*PI);
    h_phi_ee2 = new TH1D("h_phi_ee2", "phi_ee (w/ Z_Diff Cut)",1000,0.,2*PI);
    h_phi_ee3 = new TH1D("h_phi_ee3", "phi_ee (w/ Z Cut)",1000,0.,2*PI);
    h_phi_ee4 = new TH1D("h_phi_ee4", "phi_ee (w/ Z, AC Cut)",1000,0.,2*PI);
    h_theta_ek = new TH1D("h_theta_ek", "theta_ek",1000,0.1,0.35);
    h_phi_ek = new TH1D("h_phi_ek", "phi_ek",1000,3*PI/2-1.,3*PI/2+1.);
    h_theta_g = new TH1D("h_theta_g", "theta_g",1000,0.1,0.35);
    h_phi_g = new TH1D("h_phi_g", "phi_g",1000,3*PI/2-1.,3*PI/2+1.);
    h_theta_gk_lab = new TH1D("h_theta_gk_lab", "theta_gk_lab",1000,0.,0.2);
    h_theta_gk_cm = new TH1D("h_theta_gk_cm", "theta_gk_cm",1000,0.,0.3);
    h_cos_gk_lab = new TH1D("h_cos_gk_lab", "cos_gk_lab",1000,0.97,1.0);
    h_cos_gk_cm = new TH1D("h_cos_gk_cm", "cos_gk_cm",1000,0.8,1.0);
    h_mom_g = new TH1D("h_mom_g", "mom_g",1000,1.8,2.5);
    h_qsq = new TH1D("h_qsq", "Q^2",1000,0.,0.8);
    h_w = new TH1D("h_w", "W",1000,0.,0.8);
    h_pL_pR = new TH2D("h_pL_pR", "pR:pL",30,1.73,1.93,30,1.95,2.25);
    h_pL_pR2 = new TH2D("h_pL_pR2", "pR:pL (bestcut)",30,1.73,1.93,30,1.95,2.25);
    h_pL_pR3 = new TH2D("h_pL_pR3", "pR:pL (top-quality)",30,1.73,1.93,30,1.95,2.25);
    h_thph_ek = new TH2D("h_thph_ek", "theta_ek:phi_ek" ,1000,0.1,0.35,1000,3*PI/2-1.,3*PI/2+1.);
    h_thph_g = new TH2D("h_thph_g", "theta_g:phi_g" ,1000,0.1,0.35,1000,3*PI/2-1.,3*PI/2+1.);
    h_theta_ee_p = new TH2D("h_theta_ee_p", "theta_ee:mom (w/ theta_ee Cut)",1000,0.1,0.35,1000,2.0,2.25);
    h_theta_ee_p2 = new TH2D("h_theta_ee_p2", "theta_ee:mom (w/ VP Flux Cut)",1000,0.1,0.35,1000,2.0,2.25);
    //h_theta_ee_p2 = new TH2D("h_theta_ee_p2", "theta_ee:mom (w/ Z_Diff Cut)",1000,0.1,0.35,1000,2.0,2.25);
    h_theta_ee_p3 = new TH2D("h_theta_ee_p3", "theta_ee:mom (w/ Z Cut)",1000,0.1,0.35,1000,2.0,2.25);
    h_theta_ee_p4 = new TH2D("h_theta_ee_p4", "theta_ee:mom (w/ Z, AC Cut)",1000,0.1,0.35,1000,2.0,2.25);
    h_phi_ee_p = new TH2D("h_phi_ee_p", "phi_ee:mom (w/ Z Cut)",1000,1.0,2.0,1000,2.0,2.25);
    h_vpflux = new TH1D("h_vpflux", "VP Flux [/GeV/sr] (top quality)",1000,0.001,0.006);
    h_vpflux2 = new TH1D("h_vpflux2", "VP Flux [/GeV/sr] (acceptance)",1000,0.001,0.006);
    h_vpflux3 = new TH1D("h_vpflux3", "VP Flux [/GeV/sr] (w/ Z Cut)",1000,0.001,0.006);
    h_vpflux4 = new TH1D("h_vpflux4", "VP Flux [/GeV/sr] (w/ Z, AC Cut)",1000,0.001,0.006);

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
    h_zz    = new TH2D("h_zz"   ,"h_zz"   , 400,-0.25,0.25,400,  -0.25, 0.25); 


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


cout<<"Others Set TH"<<endl;
    hcoin_tc=new TH1F("hcoin_tc","Coincidence time w/ Path Length Correction  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
		for (int i=0;i<nth;i++){
			hcoin_t2[i]=new TH1F(Form("hcoin_t2[%d]",i), Form("Coincidence %lf<AC2<%lf  cut",ac1_adc[i],th2_max),bin_coin_c,min_coin_c,max_coin_c);
    		set->SetTH1(hcoin_t2[i],"Coincidence time ","","");
    		hcoin_k_ac1[i]=new TH1F(Form("hcoin_k_ac1[%d]",i), Form("Cointime (Kaon) AC1<%lf  cut",ac1_adc[i]),bin_coin_c,min_coin_c,max_coin_c);
    		set->SetTH1(hcoin_k_ac1[i],Form("Cointime (Kaon) AC1<%lf  cut",ac1_adc[i]),"","");
    		hcoin_k_ac2[i]=new TH1F(Form("hcoin_k_ac2[%d]",i), Form("Cointime (Kaon) 1000<AC2<%lf  cut",ac2l_adc[i]),bin_coin_c,min_coin_c,max_coin_c);
    		set->SetTH1(hcoin_k_ac2[i],Form("Cointime (Kaon) %lf<AC2<4000  cut",ac2l_adc[i]),"","");
		}


    //-------No Z cut------------//
	hcoin_k_fom_noZ=new TH1F("hcoin_k_fom_noZ", "No Z cut",bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_k_fom_noZ,"No Z cut","","");
	hcoin_bg_fom_noZ=new TH1F("hcoin_bg_fom_noZ", "No Z cut", bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_bg_fom_noZ,"No Z cut","","");
	hcoin_wo_bg_fom_noZ=new TH1F("hcoin_wo_bg_fom_noZ", "No Z cut", bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_wo_bg_fom_noZ,"No Z cut","","");
	hcoin_wo_bg_fom_noZ2=new TH1F("hcoin_wo_bg_fom_noZ2", "No Z cut", bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_wo_bg_fom_noZ2,"No Z cut","","");

	hmm_L_fom_noZ=new TH1F("hmm_L_fom_noZ", "No Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_L_fom_noZ,"No Z cut","","");
	hmm_bg_fom_noZ=new TH1F("hmm_bg_fom_noZ", "No Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_bg_fom_noZ,"No Z cut","","");
	hmm_wo_bg_fom_noZ=new TH1F("hmm_wo_bg_fom_noZ", "No Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_wo_bg_fom_noZ,"No Z cut","","");
	hmm_pi_fom_noZ=new TH1F("hmm_pi_fom_noZ", "No Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_pi_fom_noZ,"No Z cut (Pion Selected)","","");
	hmm_pibg_fom_noZ=new TH1F("hmm_pibg_fom_noZ", "No Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_pibg_fom_noZ,"No Z cut","","");
	hmm_pi_wobg_fom_noZ=new TH1F("hmm_pi_wobg_fom_noZ", "No Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_pi_wobg_fom_noZ,"No Z cut","","");
	hmm_pi_fom_best=new TH1F("hmm_pi_fom_best", "w/ Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_pi_fom_best,"w/ Z cut (Pion Selected)","","");
	hmm_Al_fom_best=new TH1F("hmm_Al_fom_best", "w/ Z cut",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm_Al_fom_best,"w/ Z cut (Al Selected)","","");

    //-------Z vertex test------------//
    h_zz1    = new TH2D("h_zz1"   ,"h_zz1"   ,bin_coin_c,min_coin_c,max_coin_c,bin_mm,min_mm,max_mm); 
    h_zz2    = new TH2D("h_zz2"   ,"h_zz2"   ,bin_coin_c,min_coin_c,max_coin_c,bin_mm,min_mm,max_mm); 
    h_zz3    = new TH2D("h_zz3"   ,"h_zz3"   ,bin_coin_c,min_coin_c,max_coin_c,bin_mm,min_mm,max_mm); 
    h_zz4    = new TH2D("h_zz4"   ,"h_zz4"   ,bin_coin_c,min_coin_c,max_coin_c,bin_mm,min_mm,max_mm); 
    h_z1        = new TH1D("h_z1"       ,"h_z1"       ,bin_coin_c,min_coin_c,max_coin_c);
    h_z2        = new TH1D("h_z2"       ,"h_z2"       ,bin_coin_c,min_coin_c,max_coin_c);
    h_z3        = new TH1D("h_z3"       ,"h_z3"       ,bin_coin_c,min_coin_c,max_coin_c);
    h_z4        = new TH1D("h_z4"       ,"h_z4"       ,bin_coin_c,min_coin_c,max_coin_c);
    h_z11        = new TH1D("h_z11"       ,"h_z11"       ,bin_mm,min_mm,max_mm);
    h_z22        = new TH1D("h_z22"       ,"h_z22"       ,bin_mm,min_mm,max_mm);
    h_z33        = new TH1D("h_z33"       ,"h_z33"       ,bin_mm,min_mm,max_mm);
    h_z44        = new TH1D("h_z44"       ,"h_z44"       ,bin_mm,min_mm,max_mm);
    hcoin_k_fom_Zdiff=new TH1F("hcoin_k_fom_Zdiff", "Zdiff cut",bin_coin_c,min_coin_c,max_coin_c);
    set->SetTH1(hcoin_k_fom_Zdiff,"No Z cut","","");
    hmm_L_fom_Zdiff=new TH1F("hmm_L_fom_Zdiff", "Zdiff cut",bin_mm,min_mm,max_mm);
    set->SetTH1(hmm_L_fom_Zdiff,"Zdiff cut","","");
    hcoin_k_fom_Zsum=new TH1F("hcoin_k_fom_Zsum", "Zsum cut",bin_coin_c,min_coin_c,max_coin_c);
    set->SetTH1(hcoin_k_fom_Zsum,"No Z cut","","");
    hmm_L_fom_Zsum=new TH1F("hmm_L_fom_Zsum", "Zsum cut",bin_mm,min_mm,max_mm);
    set->SetTH1(hmm_L_fom_Zsum,"Zsum cut","","");
    hct_test=new TH1F("hct_test", "Cointime_before",bin_coin_c,min_coin_c,max_coin_c);
    set->SetTH1(hct_test,"Cointime_before","","");
    hct_test2=new TH1F("hct_test2", "Cointime_after",bin_coin_c,min_coin_c,max_coin_c);
    set->SetTH1(hct_test2,"Cointime_after","","");
    hct_test3=new TH1F("hct_test3", "Itabashi_Cointime",bin_coin_c,min_coin_c,max_coin_c);
    set->SetTH1(hct_test3,"Itabashi_Cointime","","");
    h_ctct    = new TH2F("h_ctct"   ,"ct(w/o corr.):ct(w/ corr.)"  ,7145 ,-200.,200.,7145,-200.,200.); 
    h_ct_x    = new TH2F("h_ct_x"   ,"Cointime vs X"  ,bin_coin_c,min_coin_c,max_coin_c,200,-1.,1.); 
    h_ct_y    = new TH2F("h_ct_y"   ,"Cointime vs Y"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.1,0.1); 
    h_ct_th    = new TH2F("h_ct_th"   ,"Cointime vs theta"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.1,0.1); 
    h_ct_ph    = new TH2F("h_ct_ph"   ,"Cointime vs phi"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.1,0.1); 
    h_ct_xy    = new TH2F("h_ct_xy"   ,"Cointime vs XY"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.05,0.05); 
    h_ct_xth    = new TH2F("h_ct_xth"   ,"Cointime vs X*theta"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.01,0.05); 
    h_ct_yth    = new TH2F("h_ct_yth"   ,"Cointime vs Y*theta"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.002,0.002); 
    h_ct_xph    = new TH2F("h_ct_xph"   ,"Cointime vs X*phi"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.02,0.02); 
    h_ct_yph    = new TH2F("h_ct_yph"   ,"Cointime vs Y*phi"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.002,0.002); 
    h_ct_thph    = new TH2F("h_ct_thph"   ,"Cointime vs theta*phi"  ,bin_coin_c,min_coin_c,max_coin_c,20,-0.002,0.002); 
    h_ct_tgth    = new TH2F("h_ct_tgth"   ,"Cointime vs theta_t"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.1,0.1); 
    h_ct_tgph    = new TH2F("h_ct_tgph"   ,"Cointime vs phi_t"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.1,0.1); 
    h_ct_vz    = new TH2F("h_ct_vz"   ,"Cointime vs vertex_z"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.15,0.15); 
    h_ct_tgthtgph    = new TH2F("h_ct_tgthtgph"   ,"Cointime vs theta_t*phi_t"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.005,0.005); 
    h_ct_tgthz    = new TH2F("h_ct_tgthz"   ,"Cointime vs theta_t*vertex_z"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.015,0.015); 
    h_ct_tgphz    = new TH2F("h_ct_tgphz"   ,"Cointime vs phi_t*vertex_z"  ,bin_coin_c,min_coin_c,max_coin_c,200,-0.015,0.005); 

	h_gbetaR = new TH1D("h_gbetaR", "h_gbetaR", 1000, 0.995,1.0);
	h_gbetaL = new TH1D("h_gbetaL", "h_gbetaL", 1000, 0.998,1.003);
	h_gLenR = new TH1D("h_gLenR", "LenR [m]", 1000, 24., 28.);
	h_gLenL = new TH1D("h_gLenL", "LenL [m]", 1000, 24., 28.);
	h_gpR = new TH1D("h_gpR", "R_tr_p[rt] [GeV/c]", 1000, 1.6, 2.0);
	h_gpL = new TH1D("h_gpL", "L_tr_p[lt] [GeV/c]", 1000, 1.9, 2.3);
	h_gcorR = new TH1D("h_gcorR", "cor_R", 1000, 70., 80.);
	h_gcorL = new TH1D("h_gcorL", "cor_L", 1000, 70., 80.);
	h_gcorLR = new TH2D("h_gcorLR", "cor_L:cor_R", 1000, 73., 77., 1000, 73., 77.);
	h_gtref_R = new TH1D("hgtref_R", "tref_R", 1000, -1000., 4000.);
	h_gtimeR_R = new TH1D("h_gtimeR_R", "timeR_R", 1000, -1000., 4000.);
	h_gtimeL_R = new TH1D("h_gtimeL_R", "timeL_R", 1000, -1000., 4000.);
	h_gmeantime = new TH1D("h_gmeantime", "meantime_R", 1000, -200., 200.);
	h_gctcorR = new TH1D("h_gctcorR", "ctimecorR", 1000, -10., 10.);
	h_gctcorL = new TH1D("h_gctcorL", "ctimecorL", 1000, -10., 10.);


	for (int i=0;i<nth;i++){
//		for (int j=0;j<nth;j++){
//			for (int l=0;l<nth;l++){
				int j=0; int l=0;
				hcoin_k_fom[i][j][l]=new TH1F(Form("hcoin_k_fom[%d][%d][%d]",i,j,l), Form("Cointime (Kaon) (zR+zL)/2<%lf  cut",zver[i]),bin_coin_c,min_coin_c,max_coin_c);
				set->SetTH1(hcoin_k_fom[i][j][l],Form("Cointime (Kaon) (zR+zL)/2<%lf cut",zver[i]),"","");
				hcoin_bg_fom[i][j][l]=new TH1F(Form("hcoin_bg_fom[%d][%d][%d]",i,j,l), Form("Cointime (Kaon) AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),bin_coin_c,min_coin_c,max_coin_c);
				set->SetTH1(hcoin_bg_fom[i][j][l],Form("Cointime (Kaon) AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),"","");
				hcoin_wo_bg_fom[i][j][l]=new TH1F(Form("hcoin_wo_bg_fom[%d][%d][%d]",i,j,l), Form("Cointime (Kaon) AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),bin_coin_c,min_coin_c,max_coin_c);
				set->SetTH1(hcoin_wo_bg_fom[i][j][l],Form("Cointime (Kaon) AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),"","");


				hmm_L_fom[i][j][l]=new TH1F(Form("hmm_L_fom[%d][%d][%d]",i,j,l), Form("Missing Mass AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),bin_mm,min_mm,max_mm);
				set->SetTH1(hmm_L_fom[i][j][l],Form("Missing Mass AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),"","");
				hmm_bg_fom[i][j][l]=new TH1F(Form("hmm_bg_fom[%d][%d][%d]",i,j,l), Form("Missing Mass AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),bin_mm,min_mm,max_mm);
				set->SetTH1(hmm_bg_fom[i][j][l],Form("Missing Mass AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),"","");
				hmm_wo_bg_fom[i][j][l]=new TH1F(Form("hmm_wo_bg_fom[%d][%d][%d]",i,j,l), Form("Missing Mass AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),bin_mm,min_mm,max_mm);
				set->SetTH1(hmm_wo_bg_fom[i][j][l],Form("Missing Mass AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),"","");
				hmm_pi_fom[i][j][l]=new TH1F(Form("hmm_pi_fom[%d][%d][%d]",i,j,l), Form("Missing Mass (Pion Selected) AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),bin_mm,min_mm,max_mm);
				set->SetTH1(hmm_pi_fom[i][j][l],Form("Missing Mass (Pion Selected) AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),"","");
				hmm_pibg_fom[i][j][l]=new TH1F(Form("hmm_pibg_fom[%d][%d][%d]",i,j,l), Form("Missing Mass AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),bin_mm,min_mm,max_mm);
				set->SetTH1(hmm_pibg_fom[i][j][l],Form("Missing Mass AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),"","");
				hmm_pi_wobg_fom[i][j][l]=new TH1F(Form("hmm_pi_wobg_fom[%d][%d][%d]",i,j,l), Form("Missing Mass AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),bin_mm,min_mm,max_mm);
				set->SetTH1(hmm_pi_wobg_fom[i][j][l],Form("Missing Mass AC1<%lf, %lf<AC2<%lf  cut",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]),"","");

				hcoin_pi[i][j][l]=new TH1F(Form("hcoin_pi[%d][%d][%d]",i,j,l),"Coincidence time (Pion Regenerated) ",bin_coin_c,min_coin_c,max_coin_c);
//			}//for l
//		}//for j
	}//for i




    hmm=new TH1F("hmm","hcoin",bin_mm,min_mm,max_mm);
	set->SetTH1(hmm,"Mass w/o AC tuning","Mass [GeV]","Counts/2 MeV"); 
    hmm->GetXaxis()->SetRangeUser(1.0,1.3);
    hmm->GetYaxis()->SetRangeUser(0.0,450.0);

    hcoin_k=new TH1F("hcoin_k","Coincidence time w/ Correction Kaon Cut  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
    //hcoin_pi=new TH1F("hcoin_pi","Coincidence time w/ Correction Pion  Cut S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
    hcoin_p=new TH1F("hcoin_p","Coincidence time w/ Correction Proton  Cut S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
    hcoin_pi_noZ=new TH1F("hcoin_pi_noZ","Coincidence time (Pion Regenerated) ",bin_coin_c,min_coin_c,max_coin_c);
 
//    facc_kc=new TF1("facc_kc","[0]",min_coin_c,max_coin_c);
//    facc_kc->SetNpx(2000);
//    fk_kc=new TF1("fk_kc","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
//    fk_kc->SetNpx(2000);
//    fpi_pic=new TF1("fpi_pic","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
//    fpi_pic->SetNpx(2000);
//    fp_pc=new TF1("fp_pc","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
//    fp_pc->SetNpx(2000);
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

	int ev = 0;
	for(int k=0;k<ENum;k++){
		tree->GetEntry(k);

        tr.trig=R_evtype;//trigget 5 (Coin.)
        tr.nrun=(int)runnum;//runnum --> tritium???????.root
        tr.nev=k;//Entry

	if(k==ev*100000){
cout << "Event (Fill) : " << k << "/" << ENum << endl;
	ev += 1;
	}
 
///Cointime///
////////////////////from ana_Lambda.cc line 978 
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
      tr.ntrack_l=NLtr;
      tr.ntrack_r=NRtr;
	  h_L_tr_n->Fill(NLtr);
	  h_R_tr_n->Fill(NRtr);
	if(NLtr>1){//Multi-track
		NLtr = 0;
		NRtr = 0;
	}
      
      for(int lt=0;lt<NLtr;lt++){
        L_Tr = L_FP = false;
        // FP and chi2 cuts
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

#ifdef F1TDC
      convertF1TDCR(param);
      R_s0_t = RS0_F1time[0];
      for(int i=0;i<16;i++){
        if(RS2_F1time[i]>-9999.)R_s2_t[i] =  RS2_F1time[i];
        else R_s2_t[i] = -99.;
      }
#endif
	  
  
	    //---- Initialization ----//
	    tr.pid_cut = 0;
	    tr.ct_cut  = 0;
	    tr.z_cut   = 0;
	    tr.missing_mass=-100000.;
	    tr.missing_mass_b   =-100000.;
	    tr.missing_mass_L   =-100000.;
	    tr.missing_mass_nnL =-100000.;
	    tr.missing_mass_H3L =-100000.;
	    tr.missing_mass_cut =-100000.;
	    tr.missing_mass_MgL=-100000.;
	    tr.missing_mass_MgL_acc =-100000.;
	    tr.missing_mass_acc =-100000.;
		tr.ct_b=-1000;
	    tr.ct_c=-1000.;
	    tr.ct_g=-1000.;
	    tr.ct_gb=-1000.;
		tr.yp_cor=-1000.;
		tr.ctimecorR=-1000.;
 		tr.ctimecorL=-1000.;
		tr.Rtof[rt]=-1000.;
		tr.Ltof[lt]=-1000.;
		tr.RS2T_ref=-100000.;
    	tr.RS2B_ref=-100000.;
    	tr.LS2T_ref=-100000.;
    	tr.LS2B_ref=-100000.;
    	tr.RS2T_F1[rt]=-100000.;
    	tr.RS2B_F1[rt]=-100000.;
    	tr.LS2T_F1[lt]=-100000.;
    	tr.LS2B_F1[lt]=-100000.;
    	tr.RS2T_F1_c[rt]=-100000.;
    	tr.RS2B_F1_c[rt]=-100000.;
    	tr.LS2T_F1_c[lt]=-100000.;
    	tr.LS2B_F1_c[lt]=-100000.;
    	tr.RS2T_F1_b[rt]=-100000.;
    	tr.RS2B_F1_b[rt]=-100000.;
    	tr.LS2T_F1_b[lt]=-100000.;
    	tr.LS2B_F1_b[lt]=-100000.;
	    tr.momR=-1000.;
		tr.momL=-1000.;

	    ct=-1000.0;


	    tr.dpe     = -100.;
	    tr.dpk[rt] = -100.;
	    tr.dpe_[lt]= -100.;
	    
	    tr.missing_mass_Al  =-100000.;
	    tr.missing_mass_Lb  =-100000.;
	    tr.missing_mass_nnLb=-100000.;
	    tr.missing_mass_Al=-100000.;
	    tr.missing_mass_Al_bg=-100000.;
	    tr.Rpathl=-100.; tr.Lpathl=-100.;
	    tr.Rpathl_c=-100.; tr.Lpathl_c=-100.;
	  //==== tree_out for mixed event analysis =======//
	    tr.Lp_c = -100.;
	    tr.Rp_c = -100.;
	    tr.Bp_c = -100.;
	    tr.Lp_before = -100.;
	    tr.Rp_before = -100.;
	    tr.Bp_before = -100.;
	  //==== relevant to tree =======//
		L_p = 0.;
		R_p = 0.;

		tr.chi2_l = -1000.;
		tr.x_l    = -1000.;
		tr.y_l    = -1000.;
		tr.th_l   = -1000.;
		tr.ph_l   = -1000.;
		tr.mom_l  = -1000.;
		tr.tg_th_l= -1000.;
		tr.tg_ph_l= -1000.;
		tr.vz_l   = -1000.;
		tr.chi2_r = -1000.;
		tr.x_r    = -1000.;
		tr.y_r    = -1000.;
		tr.th_r   = -1000.;
		tr.ph_r   = -1000.;
		tr.mom_r  = -1000.;
		tr.tg_th_r= -1000.;
		tr.tg_ph_r= -1000.;
		tr.vz_r   = -1000.;
		tr.pathl_r= -1000.; 
    	tr.pathl_l= -1000.; 
    	tr.trpath_r= -1000.;
    	tr.trpath_l= -1000.;
    	tr.theta_gk_cm= -1000.;
    	tr.Qsq= -1000.;
    	tr.W= -1000.;
    	tr.eps= -1000.;
		tr.ct_orig= -1000.;
		tr.ct_itabashi = -1000.;

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
	  npe_sum_a2->Fill(tr.AC2_npe_sum);

	  Kaon = false; // Kaon cut 
	  zcut=false; // Vertex Z cut
	  cut_ac1=false;  //AC1 cut
	  cut_ac2=false;  //AC2 cut
	  bestcut=false;  //Z cut && AC cut
	  if( tr.AC1_npe_sum < th1_max && tr.AC2_npe_sum > th2_min) Kaon = true;
	  if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1)zcut=true;
	  if(tr.AC1_npe_sum<3.75)cut_ac1=true;
	  if(3.<tr.AC2_npe_sum && tr.AC2_npe_sum < 20.)cut_ac2=true;
	  if(tr.AC1_npe_sum<3.75 && 3.<tr.AC2_npe_sum && tr.AC2_npe_sum < 20. && fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1)bestcut=true;

	  if( L_Tr && L_FP && R_Tr && R_FP ){



	    B_p     = HALLA_p/1000.0;// [GeV/c]	    
	    L_p     = L_tr_p[lt];
	    R_p     = R_tr_p[rt];
	    
	    //==== Energy Loss calibration ======//

//	    double B_pc,R_pc,L_pc;

	    tr.dpe     = Eloss(0.0,R_tr_vz[0],"B");
//if(tr.dpe>0.)cout<<"ELoss(B) = "<<tr.dpe<<endl;
	    tr.dpk[rt] = Eloss(R_tr_tg_ph[rt],R_tr_vz[rt],"R");
//if(tr.dpk[rt]>0.)cout<<"ELoss(R) = "<<tr.dpk[rt]<<endl;
	    tr.dpe_[lt]= Eloss(L_tr_tg_ph[lt],L_tr_vz[lt],"L");
//if(tr.dpe_[lt]>0.)cout<<"ELoss(L) = "<<tr.dpe_[lt]<<endl;
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

		/* -----When CoinCalc_gagami is used----------*/
		h_gpR->Fill(R_tr_p[rt]);
		h_gpL->Fill(L_tr_p[lt]);
//kazuki


	    
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

	    
	    h_Rth->Fill(R_tr_tg_th[rt]);
	    h_Rph->Fill(R_tr_tg_ph[rt]);
	    h_Rp->Fill(R_p);
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

	    tr.Lp_before = L_p;
	    tr.Rp_before = R_p;
	    tr.Bp_before = B_p;

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
	    tr.Lp_c = L_p;
	    tr.Rp_c = R_p;
	    tr.Bp_c = B_p;
	    
	  //==== tree_out for mixed event analysis =======//
	  //Fill again after event selection
		tr.chi2_l = L_tr_chi2[lt];
		tr.x_l    = L_tr_x[lt];
		tr.y_l    = L_tr_y[lt];
		tr.th_l   = L_tr_th[lt];
		tr.ph_l   = L_tr_ph[lt];
		tr.mom_l  = L_tr_p[lt];
		tr.tg_th_l= L_tr_tg_th[lt];
		tr.tg_ph_l= L_tr_tg_ph[lt];
		tr.vz_l   = L_tr_vz[lt];

		tr.chi2_r = R_tr_chi2[rt];
		tr.x_r    = R_tr_x[rt];
		tr.y_r    = R_tr_y[rt];
		tr.th_r   = R_tr_th[rt];
		tr.ph_r   = R_tr_ph[rt];
		tr.mom_r  = R_tr_p[rt];
		tr.tg_th_r= R_tr_tg_th[rt];
		tr.tg_ph_r= R_tr_tg_ph[rt];
		tr.vz_r   = R_tr_vz[rt];
		tr.pathl_r= R_tr_pathl[rt]; 
    	tr.pathl_l= L_tr_pathl[lt]; 
    	tr.trpath_r= R_s2_trpath[rt];
    	tr.trpath_l= L_s2_trpath[lt];



	    //======= W/ Matrix calibraiton ==========================//

            double Ee;

	    Ee =sqrt(B_p*B_p + Me*Me);
	    R_E =sqrt(R_p*R_p + MK*MK);
	    L_E =sqrt(L_p*L_p + Me*Me);


	    //===== Right Hand Coordinate ====//
	    //th and phi are originally meant tan(theta) and tan(phi),
	    //so, they should not be treated like tan(R_tr_tr_th) //2020.6.30 Okuyama

	    double R_pz = R_p/sqrt(1.0*1.0 + pow((R_tr_tg_th[rt]), 2.0) + pow(( R_tr_tg_ph[rt]),2.0) );
	    double R_px = R_pz * (R_tr_tg_th[rt] );
	    double R_py = R_pz * ( R_tr_tg_ph[rt] );

	    double L_pz = L_p/sqrt(1.0*1.0 + pow(( L_tr_tg_th[lt] ), 2.0) + pow(( L_tr_tg_ph[lt]),2.0));
	    double L_px = L_pz * ( L_tr_tg_th[lt] );
	    double L_py = L_pz * ( L_tr_tg_ph[lt] );




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

	    
		//Missing mass calc.
            mass = sqrt( (Ee + mt - L_E - R_E)*(Ee + mt - L_E - R_E)
                              - (B_v - L_v - R_v)*(B_v - L_v - R_v) );

//Changed
	    h_Rz->Fill(R_tr_vz[rt]);
	    if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025)h_Lz->Fill((L_tr_vz[lt]+R_tr_vz[rt])/2.);

//================================================================//
//======== KINEMATICS START ======================================//
//========   --> smart way is perfomed in kinematics.C ===========//
//================================================================//
//================================================================//

        int s2pad = (int)R_s2_trpad[rt];
        double beta = -99.;
		double m2=-99.;
//        if( R_s2_t[s2pad]>0 && R_s0_t>0 && s2pad>=0 )
	    double p    = R_tr_p[rt];//GeV
        double path = R_s2_trpath[rt] - R_s0_trpath[rt];//m
        beta = path / ( R_s2_t[s2pad] - R_s0_t ) / c;
		//cout<<"p="<<p<<endl;
		//cout<<"c="<<c<<endl;
		//cout<<"path="<<path<<endl;
//		cout<<"S2time="<<R_s2_t[s2pad]<<endl;
//		cout<<"S0time="<<R_s0_t<<endl;
//		cout<<"beta="<<beta<<endl;
        //double beta=R_tr_p[rt]/sqrt(R_tr_p[rt]*R_tr_p[rt]+MK*MK);
	    double pL    = L_tr_p[lt];//GeV
	    double pR    = R_tr_p[rt];//GeV
		double theta = L_tr_tg_th[lt];
		double theta_R = R_tr_tg_th[rt];
		double phi = L_tr_tg_ph[lt];
		double phi_R = R_tr_tg_ph[rt];
		double phi0=13.2*PI/180;//rad
		//double theta_L = acos((-tan(phi)*sin(phi0)+cos(phi0))/(sqrt(1+tan(theta)*tan(theta)+tan(phi)*tan(phi))));//LHRS frame
		double theta_L = acos(1./(sqrt(1+theta*theta+phi*phi)));//LHRS frame
		double theta_ee = acos((-phi*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
		double theta_ek = acos((phi_R*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
		//double phi_L = atan((phi*cos(phi0)+sin(phi0))/theta);//LHRS frame
		double phi_L = 0.;//LHRS frame
		double phi_RHRS = 0.;//RHRS frame
		double phi_ee = 0.;//original frame
		double phi_ek = 0.;//original frame

		double Escat = sqrt(pL*pL+Me*Me); 
		double Einc = 4.3; 
		double omega = Einc - Escat;
		if(theta>0.&&phi>0.){phi_L = atan(phi/theta);}//LHRS frame
		if(theta<0.&&phi>0.){phi_L = atan(phi/theta)+PI;}
		if(theta<0.&&phi<0.){phi_L = atan(phi/theta)+PI;}
		if(theta>0.&&phi<0.){phi_L = atan(phi/theta)+2*PI;}
		if((theta/(-phi*sin(phi0)+cos(phi0)))>0.&&((phi*cos(phi0)+sin(phi0))/(-phi*sin(phi0)+cos(phi0)))>0.){phi_ee = atan((phi*cos(phi0)+sin(phi0))/theta);}//original frame
		if((theta/(-phi*sin(phi0)+cos(phi0)))<0.&&((phi*cos(phi0)+sin(phi0))/(-phi*sin(phi0)+cos(phi0)))>0.){phi_ee = atan((phi*cos(phi0)+sin(phi0))/theta)+PI;}
		if((theta/(-phi*sin(phi0)+cos(phi0)))<0.&&((phi*cos(phi0)+sin(phi0))/(-phi*sin(phi0)+cos(phi0)))<0.){phi_ee = atan((phi*cos(phi0)+sin(phi0))/theta)+PI;}
		if((theta/(-phi*sin(phi0)+cos(phi0)))>0.&&((phi*cos(phi0)+sin(phi0))/(-phi*sin(phi0)+cos(phi0)))<0.){phi_ee = atan((phi*cos(phi0)+sin(phi0))/theta)+2*PI;}
		//if((theta/(-phi*sin(phi0)+cos(phi0)))>0.&&(theta/(phi*cos(phi0)+sin(phi0)))>0.){phi_ee = atan((phi*cos(phi0)+sin(phi0))/theta);}//original frame
		//if((theta/(-phi*sin(phi0)+cos(phi0)))<0.&&(theta/(phi*cos(phi0)+sin(phi0)))>0.){phi_ee = atan((phi*cos(phi0)+sin(phi0))/theta)+PI;}
		//if((theta/(-phi*sin(phi0)+cos(phi0)))<0.&&(theta/(phi*cos(phi0)+sin(phi0)))<0.){phi_ee = atan((phi*cos(phi0)+sin(phi0))/theta)+PI;}
		//if((theta/(-phi*sin(phi0)+cos(phi0)))>0.&&(theta/(phi*cos(phi0)+sin(phi0)))<0.){phi_ee = atan((phi*cos(phi0)+sin(phi0))/theta)+2*PI;}
		
		if(theta_R>0.&&phi_R>0.){phi_RHRS = atan(phi_R/theta_R);}//RHRS frame
		if(theta_R<0.&&phi_R>0.){phi_RHRS = atan(phi_R/theta_R)+PI;}
		if(theta_R<0.&&phi_R<0.){phi_RHRS = atan(phi_R/theta_R)+PI;}
		if(theta_R>0.&&phi_R<0.){phi_RHRS = atan(phi_R/theta_R)+2*PI;}
		if((theta_R/(phi_R*sin(phi0)+cos(phi0)))>0.&&(theta_R/(phi_R*cos(phi0)-sin(phi0)))>0.){phi_ek = atan((phi_R*cos(phi0)-sin(phi0))/theta_R);}//original frame
		if((theta_R/(phi_R*sin(phi0)+cos(phi0)))<0.&&(theta_R/(phi_R*cos(phi0)-sin(phi0)))>0.){phi_ek = atan((phi_R*cos(phi0)-sin(phi0))/theta_R)+PI;}
		if((theta_R/(phi_R*sin(phi0)+cos(phi0)))<0.&&(theta_R/(phi_R*cos(phi0)-sin(phi0)))<0.){phi_ek = atan((phi_R*cos(phi0)-sin(phi0))/theta_R)+PI;}
		if((theta_R/(phi_R*sin(phi0)+cos(phi0)))>0.&&(theta_R/(phi_R*cos(phi0)-sin(phi0)))<0.){phi_ek = atan((phi_R*cos(phi0)-sin(phi0))/theta_R)+2*PI;}

		double mom_g=sqrt(pL*pL*sin(theta_ee)*sin(theta_ee)+(4.3-pL*cos(theta_ee))*(4.3-pL*cos(theta_ee)));
		double Qsq=mom_g*mom_g-omega*omega;
		double phi_g = 0.; 
		if((theta/(-phi*sin(phi0)+cos(phi0)))<0.&&((phi*cos(phi0)+sin(phi0))/(-phi*sin(phi0)+cos(phi0)))<0.){phi_g=atan(phi*cos(phi0)+sin(phi0)/theta);}//original frame
		if((theta/(-phi*sin(phi0)+cos(phi0)))>0.&&((phi*cos(phi0)+sin(phi0))/(-phi*sin(phi0)+cos(phi0)))<0.){phi_g=atan(phi*cos(phi0)+sin(phi0)/theta)+PI;}
		if((theta/(-phi*sin(phi0)+cos(phi0)))>0.&&((phi*cos(phi0)+sin(phi0))/(-phi*sin(phi0)+cos(phi0)))>0.){phi_g=atan(phi*cos(phi0)+sin(phi0)/theta)+PI;}
		if((theta/(-phi*sin(phi0)+cos(phi0)))<0.&&((phi*cos(phi0)+sin(phi0))/(-phi*sin(phi0)+cos(phi0)))>0.){phi_g=atan(phi*cos(phi0)+sin(phi0)/theta)+2*PI;}
		//if((theta/(-phi*sin(phi0)+cos(phi0)))<0.&&(theta/(phi*cos(phi0)+sin(phi0)))<0.){phi_g=atan(pL*(phi*cos(phi0)+sin(phi0)/theta)/mom_g);}//original frame
		//if((theta/(-phi*sin(phi0)+cos(phi0)))>0.&&(theta/(phi*cos(phi0)+sin(phi0)))<0.){phi_g=atan(pL*(phi*cos(phi0)+sin(phi0)/theta)/mom_g)+PI;}
		//if((theta/(-phi*sin(phi0)+cos(phi0)))>0.&&(theta/(phi*cos(phi0)+sin(phi0)))>0.){phi_g=atan(pL*(phi*cos(phi0)+sin(phi0)/theta)/mom_g)+PI;}
		//if((theta/(-phi*sin(phi0)+cos(phi0)))<0.&&(theta/(phi*cos(phi0)+sin(phi0)))>0.){phi_g=atan(pL*(phi*cos(phi0)+sin(phi0)/theta)/mom_g)+2*PI;}
		//double theta_g = asin(pL*sin(theta_ee)/mom_g);
		double theta_g = acos((4.3-pL*cos(theta_ee))/mom_g);
		double pgpR=mom_g*sin(theta_g)*cos(phi_g)*pR*sin(theta_ek)*cos(phi_ek)+mom_g*sin(theta_g)*sin(phi_g)*pR*sin(theta_ek)*sin(phi_ek)+mom_g*cos(theta_g)*pR*cos(theta_ek);
		double theta_gk_lab=acos(pgpR/mom_g/pR);
	
		//--Rotation & Lorentz--//
		//double ekl=sqrt(pR*pR+MK*MK);
		//double pkx=pR*sin(theta_ek)*cos(phi_ek);
		//double pky=pR*sin(theta_ek)*sin(phi_ek)*cos(theta_g)+pR*cos(theta_ek)*sin(theta_g);
		//double pkz=-pR*sin(theta_ek)*sin(phi_ek)*sin(theta_g)+pR*cos(theta_ek)*cos(theta_g);
		double beta_cm=mom_g/(omega+Mp);
		double gamma_cm=1./(sqrt(1-beta_cm*beta_cm));
		//double ekcm=gamma_cm*ekl-gamma_cm*beta_cm*pkz;
		//double pkxcm=pkx;
		//double pkycm=pky;
		//double pkzcm=-gamma_cm*beta_cm*ekl+gamma_cm*pkz;
		//double theta_gk_cm=acos(pkzcm/(sqrt(pkxcm*pkxcm+pkycm*pkycm+pkzcm*pkzcm)));
		//--Lorentz Transformation in another frame--//
		double theta_gk_cm=atan(pR*sin(theta_gk_lab)/(-gamma_cm*beta_cm*sqrt(pR*pR+MK*MK)+gamma_cm*pR*cos(theta_gk_lab)));

//2020/11/21
		double W = sqrt((omega+Mp)*(omega+Mp)-mom_g*mom_g);
		double q2=Qsq+omega*omega;
		double eps=1/(1+2*(q2/Qsq)*tan(theta_ee/2)*tan(theta_ee/2));
    	tr.theta_gk_cm= theta_gk_cm;
    	tr.Qsq= Qsq;
    	tr.W= W;
    	tr.eps= eps;


		double A=Me*Me*omega*omega/(4*Einc*Einc*Escat*Escat);
		double sinterm=sin(theta_ee/2)*sin(theta_ee/2);
		double a1=((Einc*Einc+Escat*Escat)/(2*Einc*Einc))/(A+sinterm);
		double a2=(Escat/Einc)*A/((A+sinterm)*(A+sinterm));
		double a3=((Einc+Escat)*(Einc+Escat)/(4*Einc*Einc))/(omega*omega/(4*Einc*Escat)+sinterm);
		double vpflux=(a1-a2-a3)/(137*4*PI*PI*omega);
		Ng_det+=vpflux;
//cout test
		//if(k%100000==0){
		//cout<<"Me="<<Me<<endl;
		//cout<<"Einc="<<Einc<<endl;
		//cout<<"Escat="<<Escat<<endl;
		//cout<<"omega="<<omega<<endl;
		//cout<<"theta="<<theta_ee*180/PI<<endl;
		//cout<<"A="<<A<<endl;
		//cout<<"sinterm="<<sinterm<<endl;
		//cout<<"a1="<<a1<<endl;
		//cout<<"a2="<<a2<<endl;
		//cout<<"a3="<<a3<<endl;
		//cout<<"vpflux="<<vpflux<<endl;
		//}
		h_L_tgph_tgth12->Fill(theta,phi);//No Cut
		h_original_phth->Fill(theta_ee,phi_ee);//No Cut
		h_LHRS_phth->Fill(theta_L,phi_L);//No Cut
		if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1){
		h_theta_ee ->Fill(theta_ee);
		h_phi_ee ->Fill(phi_ee);
		h_theta_ek ->Fill(theta_ek);
		h_phi_ek ->Fill(phi_ek);
		h_theta_g ->Fill(theta_g);
		h_phi_g ->Fill(phi_g);
		h_thph_ek->Fill(theta_ek,phi_ek);
		h_thph_g->Fill(theta_g,phi_g);
		h_mom_g->Fill(mom_g);
		h_qsq->Fill(Qsq);
		if(fabs(ct)<1.)h_pL_pR->Fill(pR,pL);
		if(bestcut&&fabs(ct)<1.)h_pL_pR2->Fill(pR,pL);
		h_theta_gk_lab->Fill(theta_gk_lab);
		h_theta_gk_cm->Fill(theta_gk_cm);
		h_cos_gk_lab->Fill(cos(theta_gk_lab));
		h_cos_gk_cm->Fill(cos(theta_gk_cm));
		}
		//h_theta_ee_p ->Fill(theta_ee,pL);
		//h_theta_ee_p ->SetBinContent(theta_ee,pL,vpflux*100000);
		//h_vpflux ->Fill(vpflux);
		if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025){//Z_Diff
		//h_theta_ee2 ->Fill(theta_ee);
		h_phi_ee2 ->Fill(phi_ee);
		//h_theta_ee_p2 ->Fill(theta_ee,pL);
		//h_theta_ee_p2 ->SetBinContent(theta_ee,pL,vpflux*100000);
		//h_vpflux2 ->Fill(vpflux);
		}
		if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1){//Z Cut
		h_theta_ee3 ->Fill(theta_ee);
		h_phi_ee3 ->Fill(phi_ee);
		h_theta_ee_p3 ->Fill(theta_ee,pL);
		h_phi_ee_p ->Fill(phi_ee,pL);
		//h_theta_ee_p3 ->SetBinContent(theta_ee,pL,vpflux*100000);
		h_vpflux3 ->Fill(vpflux);
		}
		if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && tr.AC1_npe_sum<3.75 && 3.<tr.AC2_npe_sum && tr.AC2_npe_sum < 20.){//Z, AC Cut
		h_theta_ee4 ->Fill(theta_ee);
		h_phi_ee4 ->Fill(phi_ee);
		//h_theta_ee_p4 ->Fill(theta_ee,pL);
		h_theta_ee_p4 ->SetBinContent(theta_ee,pL,vpflux*100000);
		h_vpflux4 ->Fill(vpflux);
		}
		if(fabs(theta_ee-0.2475)<0.0075&&bestcut){//theta_ee Cut
		h_L_tgph_tgth->Fill(theta,phi);
		h_theta_ee_p ->Fill(theta_ee,pL);
		}
		if(fabs(vpflux-0.0028)<0.0002&&bestcut){//VP Flux Cut
		h_L_tgph_tgth2->Fill(theta,phi);
		h_theta_ee_p2 ->Fill(theta_ee,pL);
		}
		bool top_quality = false;
		bool acceptance = false;//6msr & 0.2452GeV/c
		if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(theta_ee-0.225)<0.0125&&fabs(phi_ee-1.6)<0.125&&fabs(pL-2.125)<0.025)top_quality=true;//Âç¶ã¾ú
		if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(theta_ee-0.225)<0.025&&fabs(phi_ee-1.6)<0.25&&fabs(pL-2.1)<0.1226)acceptance=true;//
		if(fabs(theta_ee-0.225)<0.0125&&fabs(phi_ee-1.6)<0.125&&fabs(pL-2.1)<0.05&&fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025){h_Lz2->Fill((L_tr_vz[lt]+R_tr_vz[rt])/2.);h_Rz2->Fill(R_tr_vz[rt]);}
		if(top_quality&&fabs(ct)<1.)h_pL_pR3->Fill(pR,pL);
		if(acceptance){
		Ng_det_acc+=vpflux;
		h_vpflux2 ->Fill(vpflux);
		}
		if(top_quality){
		h_vpflux ->Fill(vpflux);
		h_original_phth3->Fill(theta_ee,phi_ee);
		Ng_det_top+=vpflux;
		}
		if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1){
		h_original_phth4->Fill(theta_ee,phi_ee);
		if(L_tr_p[lt]>2.18)h_original_phth5->Fill(theta_ee,phi_ee);
		if(L_tr_p[lt]>2.1&&L_tr_p[lt]<2.15)h_original_phcosth->Fill(cos(theta_ee),phi_ee);
		if(L_tr_p[lt]>2.16)h_original_phcosth2->Fill(cos(theta_ee),phi_ee);
		if(L_tr_p[lt]>2.1&&L_tr_p[lt]<2.16)h_original_phth6->Fill(theta_ee,phi_ee);
		if(L_tr_p[lt]>2.1&&L_tr_p[lt]<2.12)h_original_phth7->Fill(theta_ee,phi_ee);
		if(L_tr_p[lt]>2.12&&L_tr_p[lt]<2.14)h_original_phth8->Fill(theta_ee,phi_ee);
		if(L_tr_p[lt]>2.18&&L_tr_p[lt]<2.2)h_original_phth9->Fill(theta_ee,phi_ee);
		if(L_tr_p[lt]>2.2&&L_tr_p[lt]<2.22)h_original_phth10->Fill(theta_ee,phi_ee);

		}

		//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && tr.AC1_npe_sum<3.75 && 3.<tr.AC2_npe_sum && tr.AC2_npe_sum < 20.&& fabs(ct)<1.)h_L_tgph_tgth2->Fill(theta,phi);
		if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1)h_L_tgph_tgth3->Fill(theta,phi);
		if(tr.AC1_npe_sum<3.75 && 3.<tr.AC2_npe_sum && tr.AC2_npe_sum < 20.)h_L_tgph_tgth4->Fill(theta,phi);
		if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && tr.AC1_npe_sum<3.75 && 3.<tr.AC2_npe_sum && tr.AC2_npe_sum < 20.){h_L_tgph_tgth5->Fill(theta,phi);h_original_phth2->Fill(theta_ee,phi_ee);h_LHRS_phth2->Fill(theta_L,phi_L);}//bestcut

		if(fabs(R_tr_vz[rt]-0.12)<0.02 && fabs(L_tr_vz[lt]+0.12)<0.02){h_L_tgph_tgth6->Fill(theta,phi);}
		if(fabs(R_tr_vz[rt]+0.12)<0.02 && fabs(L_tr_vz[lt]-0.12)<0.02){h_L_tgph_tgth7->Fill(theta,phi);}
		if(fabs(R_tr_vz[rt])<0.1 && fabs(L_tr_vz[lt]-0.12)<0.02){h_L_tgph_tgth8->Fill(theta,phi);}
		if(fabs(L_tr_vz[lt])<0.1 && fabs(R_tr_vz[rt]-0.12)<0.02){h_L_tgph_tgth9->Fill(theta,phi);}
		if(fabs(L_tr_vz[lt]+0.12)<0.02 && fabs(R_tr_vz[rt]+0.12)<0.02){h_L_tgph_tgth10->Fill(theta,phi);}
		if(fabs(L_tr_vz[lt]-0.12)<0.02 && fabs(R_tr_vz[rt]-0.12)<0.02){h_L_tgph_tgth11->Fill(theta,phi);}


	    h_Rth_c->Fill(R_tr_tg_th[rt]);
	    h_Rph_c->Fill(R_tr_tg_ph[rt]);
	    h_Rp_c->Fill(R_p);
	    h_Lth_c->Fill(L_tr_tg_th[lt]);	    	    
	    h_Lph_c->Fill(L_tr_tg_ph[lt]);	    	    
	    h_Lp_c->Fill(L_p);
        if( R_s2_t[s2pad]>0 && R_s0_t>0 && s2pad>=0 ){
        m2 = ( 1./beta/beta - 1. ) * p * p;
		h_R_m2->Fill(m2); 
//		cout<<"m2="<<m2<<endl;
//			double beta2 = R_v*R_v/(R_E*R_E);//beta^2
//			m2 = R_v*R_v*(1./beta2-1.);
		}
//================================================================//
//kinematics calc. and kinematical cut were performed.
//================================================================//


	    
	    mm=mass - mh;
            //mm2=mass2 - mh;

	   // mm = mm*1000.; // GeV -> MeV
	    //mm2 = mm2*100.; // GeV ->MeV 


//-------------------------------------------//
//-----Cointime correlation------------------//
//-------------------------------------------//
		//if(bestcut){
		h_ct_x->Fill(ct,L_tr_x[lt]);	
		h_ct_y->Fill(ct,L_tr_y[lt]);	
		h_ct_th->Fill(ct,L_tr_th[lt]);	
		h_ct_ph->Fill(ct,L_tr_ph[lt]);	
		h_ct_xy->Fill(ct,L_tr_x[lt]*L_tr_y[lt]);	
		h_ct_thph->Fill(ct,L_tr_th[lt]*L_tr_ph[lt]);	
		h_ct_xth->Fill(ct,L_tr_x[lt]*L_tr_th[lt]);	
		h_ct_yth->Fill(ct,L_tr_y[lt]*L_tr_th[lt]);	
		h_ct_xph->Fill(ct,L_tr_x[lt]*L_tr_ph[lt]);	
		h_ct_yph->Fill(ct,L_tr_y[lt]*L_tr_ph[lt]);	
		h_ct_tgth->Fill(ct,L_tr_tg_th[lt]);	
		h_ct_tgph->Fill(ct,L_tr_tg_ph[lt]);	
		h_ct_vz->Fill(ct,L_tr_vz[lt]);	
		h_ct_tgthtgph->Fill(ct,L_tr_tg_th[lt]*L_tr_tg_ph[lt]);	
		h_ct_tgthz->Fill(ct,L_tr_tg_th[lt]*L_tr_vz[lt]);	
		h_ct_tgphz->Fill(ct,L_tr_tg_ph[lt]*L_tr_vz[lt]);	
		//}
//-------------------------------------------//
//-------------------------------------------//


 double ctimecorR = calcf2t_3rd(PctimeR, R_tr_x[rt],R_tr_th[rt],R_tr_y[rt],R_tr_ph[rt],R_tr_vz[rt]);
 double ctimecorL = calcf2t_3rd(PctimeL, L_tr_x[lt],L_tr_th[lt],L_tr_y[lt],L_tr_ph[lt],L_tr_vz[lt]);
		//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1)//gas region
		//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + 0.12)<0.02 && fabs(L_tr_vz[lt] + 0.12)<0.02)//Al front
		//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] - 0.12)<0.02 && fabs(L_tr_vz[lt] - 0.12)<0.02)//Al back
		//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + 0.12)<0.02 && fabs(L_tr_vz[lt] + 0.12)<0.02 && fabs(ct-3.1)<0.7)//Al front && pion
//-------------------------------------------//
//-----FP Correlation------------------------//
//-------------------------------------------//
		if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025){
		h_L_y_x->Fill(L_tr_x[lt],L_tr_y[lt]);
		h_L_z_x->Fill(L_tr_x[lt],L_tr_vz[lt]);
		h_L_z_y->Fill(L_tr_y[lt],L_tr_vz[lt]);
		h_L_th_x->Fill(L_tr_x[lt],L_tr_th[lt]);
		h_L_ph_y->Fill(L_tr_y[lt],L_tr_ph[lt]);
		h_L_th_y->Fill(L_tr_y[lt],L_tr_th[lt]);
		h_L_ph_x->Fill(L_tr_x[lt],L_tr_ph[lt]);
		h_L_th_z->Fill(L_tr_vz[lt],L_tr_th[lt]);
		h_L_ph_z->Fill(L_tr_vz[lt],L_tr_ph[lt]);
		h_L_tgth_y->Fill(L_tr_y[lt],L_tr_tg_th[lt]);
		h_L_tgph_x->Fill(L_tr_x[lt],L_tr_tg_ph[lt]);
   	    h_L_tgth_z->Fill(L_tr_vz[lt],L_tr_tg_th[lt]);
   	    h_L_tgph_z->Fill(L_tr_vz[lt],L_tr_tg_ph[lt]);
   	    h_L_tgph_y->Fill(L_tr_y[lt],L_tr_tg_ph[lt]);
   	    h_L_tgth_x->Fill(L_tr_x[lt],L_tr_tg_th[lt]);

   	    h_L_corR_tgth->Fill(L_tr_tg_th[lt],ctimecorR);
   	    h_L_corR_tgph->Fill(L_tr_tg_ph[lt],ctimecorR);
   	    h_L_corR_th->Fill(L_tr_th[lt],ctimecorR);
   	    h_L_corR_ph->Fill(L_tr_ph[lt],ctimecorR);
   	    h_L_corR_x->Fill(L_tr_x[lt],ctimecorR);
   	    h_L_corR_y->Fill(L_tr_y[lt],ctimecorR);
   	    h_L_corR_z->Fill(L_tr_vz[lt],ctimecorR);
   	    h_L_corR_zdiff->Fill(R_tr_vz[rt]-L_tr_vz[lt],ctimecorR);
   	    h_L_corL_tgth->Fill(L_tr_tg_th[lt],ctimecorL);
   	    h_L_corL_tgph->Fill(L_tr_tg_ph[lt],ctimecorL);
   	    h_L_corL_th->Fill(L_tr_th[lt],ctimecorL);
   	    h_L_corL_ph->Fill(L_tr_ph[lt],ctimecorL);
   	    h_L_corL_x->Fill(L_tr_x[lt],ctimecorL);
   	    h_L_corL_y->Fill(L_tr_y[lt],ctimecorL);
   	    h_L_corL_z->Fill(L_tr_vz[lt],ctimecorL);
   	    h_L_corL_zdiff->Fill(R_tr_vz[rt]-L_tr_vz[lt],ctimecorL);

   	    h_L_ct_zdiff->Fill(R_tr_vz[rt]-L_tr_vz[lt],ct);
   	    h_L_corsum_zdiff->Fill(R_tr_vz[rt]-L_tr_vz[lt],ctimecorL+ctimecorR);
   	    h_L_cordiff_zdiff->Fill(R_tr_vz[rt]-L_tr_vz[lt],ctimecorR-ctimecorL);
		h_L_corL_corR->Fill(ctimecorR,ctimecorL);
		}
//-------------------------------------------//
//-----FP Correlation------------------------//
//-------------------------------------------//



		if(fabs(mm)<0.008){
		h_Lp_top->Fill(pL);
		}

//-------------------------------------------//
//-------------No Z cut at all--------------//
//-------------------------------------------//
	if(cut_ac1 && cut_ac2){//no Z cut
//	if(tr.AC1_npe_sum<3.75 && tr.AC2_npe_sum>3. && tr.AC2_npe_sum<20. && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1)//no Z cut,w/ diff cut

	hcoin_k_fom_noZ->Fill(ct);
	tr.ct_den=ct;
	//hct_test->Fill(ct_test);
	//h_ctct->Fill(ct,ct_test);
	h_m2_mm->Fill(m2,mm);
	//h_m2_mm->Fill(tr.AC2_npe_sum,mm);
	h_m2_ac->Fill(m2,tr.AC2_npe_sum);
//--------------Z vertex---------------//
	h_R_vz->Fill(R_tr_vz[rt]);
	h_L_vz->Fill(L_tr_vz[lt]);
	if(fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1)h_R_vz2->Fill(R_tr_vz[rt]-L_tr_vz[lt]);
	if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025)h_L_vz2->Fill((R_tr_vz[rt]+L_tr_vz[lt])/2);
	//h_zz->Fill(R_tr_vz[rt],L_tr_vz[lt]);
	//if((L_tr_vz[lt]<(R_tr_vz[rt]-0.025) && -0.08<R_tr_vz[rt] && R_tr_vz[rt] <0.1 && -0.1<L_tr_vz[lt] && L_tr_vz[lt]<0.08) || (L_tr_vz[lt]>(R_tr_vz[rt]+0.025) && -0.1<R_tr_vz[rt] && R_tr_vz[rt]<0.08 && -0.08<L_tr_vz[lt] && L_tr_vz[lt]<0.1)){h_zz->Fill(R_tr_vz[rt],L_tr_vz[lt]);h_z4->Fill(ct);if(fabs(ct-10.)<1){h_z44->Fill(mm);}}
	if(fabs(R_tr_vz[rt]-0.12)<0.02 && fabs(L_tr_vz[lt]+0.12)<0.02){h_zz1->Fill(ct_test,mm);h_z1->Fill(ct_test);if(fabs(ct_test)<1){h_z11->Fill(mm);}}
	if(fabs(R_tr_vz[rt]+0.12)<0.02 && fabs(L_tr_vz[lt]-0.12)<0.02){h_zz2->Fill(ct_test,mm);h_z2->Fill(ct_test);if(fabs(ct_test)<1){h_z22->Fill(mm);}}
	if(fabs(R_tr_vz[rt])<0.1 && fabs(L_tr_vz[lt]-0.12)<0.02){h_zz3->Fill(ct_test,mm);h_z3->Fill(ct_test);if(fabs(ct_test)<1){h_z33->Fill(mm);}}
	if(fabs(L_tr_vz[lt])<0.1 && fabs(R_tr_vz[rt]-0.12)<0.02){h_zz4->Fill(ct_test,mm);h_z4->Fill(ct_test);if(fabs(ct_test)<1){h_z44->Fill(mm);}}
	if(fabs(R_tr_vz[rt]-L_tr_vz[lt])>0.2)hcoin_k_fom_Zdiff->Fill(ct);
	if(fabs(R_tr_vz[rt]+L_tr_vz[lt])/2.>0.5)hcoin_k_fom_Zsum->Fill(ct);

					//if((-100.<ct && ct <-20.) || (20.<ct && ct<100.)){
				//}
				if(20.<ct && ct<100.){
					double ct_ = ct;
				        while(1){
					  if(-20.<ct && ct<20.){
						 hcoin_bg_fom_noZ->Fill(ct);
						 hmm_bg_fom_noZ->Fill(mm);
						 hmm_pibg_fom_noZ->Fill(mm);
						 break;}
					       else if(ct<-20.){ct=ct+40.;}
					       else if(20.<ct){ct=ct-40.;}
					 }
					ct = ct_;
					}//cointime
					if(fabs(ct)<1.){
						hmm_L_fom_noZ->Fill(mm);
						if(fabs(R_tr_vz[rt]-L_tr_vz[lt])>0.2)hmm_L_fom_Zdiff->Fill(mm);
						if(fabs(R_tr_vz[rt]+L_tr_vz[lt])/2.>0.5)hmm_L_fom_Zsum->Fill(mm);
						if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025&&(fabs(fabs(R_tr_vz[rt]+L_tr_vz[lt])/2.-0.12)<0.01||fabs(fabs(R_tr_vz[rt]+L_tr_vz[lt])/2.+0.12)<0.01))hmm_Al_fom_best->Fill(mm);
				}
					if(fabs(ct-3.0)<1.0){//3.05//0.7
						hmm_pi_fom_noZ->Fill(mm);//MM if pion
						if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1)hmm_pi_fom_best->Fill(mm);
					}
	}//No Z Cut, w/ AC Cut







//-------------------------------------------------------------------//
//--------------(i,j,l) running--------------------------------------//
//-------------------------------------------------------------------//
	for (int i=0;i<nth;i++){
//		for (int j=0;j<nth;j++){
//			for (int l=0;l<nth;l++){
				int j=0; int l=0;
//cout<<"j="<<j<<endl;



				zcut=false;
				tr.cut_par[i]=zver[i];
			//	if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<zver[i]/1000 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1)zcut=true;
				if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<zver[i])zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && pL<2.12  &&i==10)zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(theta_ee-0.225)<0.025&&fabs(phi_ee-1.6)<0.25&&fabs(pL-2.125)<0.025&&i==20)zcut=true;//Âç¶ã¾ú
				////if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(theta_ee-0.225)<0.0125&&fabs(phi_ee-1.6)<0.125&&fabs(pL-2.125)<0.025&&i==30)zcut=true;//LambdaÂç¶ã¾ú
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(theta_ee-0.225)<0.025&&fabs(phi_ee-1.6)<0.25&&fabs(pL-2.075)<0.025&&i==30)zcut=true;//SigmaÂç¶ã¾ú
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && pL>2.1 && pL<2.15  &&i==31)zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && pL>2.15 && pL<2.25  &&i==32)zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && pL>2.1 && pL<2.125  &&i==33)zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && pL>2.125 && pL<2.15  &&i==34)zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && pL>2.15 && pL<2.175  &&i==35)zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && pL>2.175 && pL<2.2  &&i==36)zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && pL>2.1 && pL<2.12  &&i==37)zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && pL>2.12 && pL<2.14  &&i==38)zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && pL>2.14 && pL<2.16  &&i==39)zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && pL>2.16 && pL<2.18  &&i==40)zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && pL>2.18 && pL<2.2  &&i==41)zcut=true;
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(theta_ee-0.225)<0.0125&&fabs(phi_ee-1.6)<0.125&&fabs(pL-2.125)<0.025&&i==51)zcut=true;//Âç¶ã¾ú
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(theta_ee-0.25)<0.0125&&fabs(phi_ee-1.6)<0.125&&fabs(pL-2.125)<0.025&&i==52)zcut=true;//Âç¶ã¾ú
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(theta_ee-0.25)<0.0125&&fabs(phi_ee-1.85)<0.125&&fabs(pL-2.125)<0.025&&i==53)zcut=true;//Âç¶ã¾ú
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(theta_ee-0.225)<0.0125&&fabs(phi_ee-1.85)<0.125&&fabs(pL-2.125)<0.025&&i==54)zcut=true;//Âç¶ã¾ú
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(cos(theta_ee)-0.9748)<0.003&&fabs(phi_ee-1.6)<0.125&&fabs(pL-2.125)<0.025&&i==56)zcut=true;//Âç¶ã¾ú
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(cos(theta_ee)-0.9808)<0.003&&fabs(phi_ee-1.6)<0.125&&fabs(pL-2.125)<0.025&&i==57)zcut=true;//Âç¶ã¾ú
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(cos(theta_ee)-0.9808)<0.003&&fabs(phi_ee-1.85)<0.125&&fabs(pL-2.125)<0.025&&i==58)zcut=true;//Âç¶ã¾ú
				//if(fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025 && fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1&&fabs(cos(theta_ee)-0.9748)<0.003&&fabs(phi_ee-1.85)<0.125&&fabs(pL-2.125)<0.025&&i==59)zcut=true;//Âç¶ã¾ú
				if( zcut && cut_ac1 && cut_ac2){
				hcoin_k_fom[i][j][l]->Fill(ct);
				tr.ct_eff[i]=ct;
			    //cout<<"hcoin_k_fom is filled" << endl;


		//-------------------------------------------//
		//---------Accidental B.G.-------------------//
		//-------------------------------------------//
				//if((-100.<ct && ct <-20.) || (20.<ct&&ct<100.)){
				//}
				if(20.<ct&&ct<100.){
					double ct_ = ct;
//kazuki
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
					if(fabs(ct-3.05)<1.0){
						// def_sig_pi=0.443; def_mean_pi=3.0;
						hmm_pi_fom[i][j][l]->Fill(mm);//MM if pion
						//if(i==20)hmm_pi_fom_best->Fill(mm);
					}
					//}
				}//if cut condition
//			}//for l
//		}//for j
	}//for i	
//-------------------------------------------//
//-------------------------------------------//

				    

    	//--------------------------From ana_Lambda.cc----------------------------------------//
    	//------------------------------------------------------------------------------------//
    	//------------------------------------------------------------------------------------//
    	//-----------------------Okuyama doesn't use these.-----------------------------------//

   
//            mass2= sqrt( (Ee + mt - L_E - R_Epi)*(Ee + mt - L_E - R_Epi)
//                              - (B_v - L_v - R_v)*(B_v - L_v - R_v) );

	    
            mass_pc = sqrt( (Eec + mt - L_Ec - R_Ec)*(Eec + mt - L_Ec - R_Ec)
                              - (B_vc - L_vc - R_vc)*(B_vc - L_vc - R_vc) );
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
		 if(ct<1.0 && -1.0<ct) hmm->Fill(mm);
    	//------------------------------------------------------------------------------------//
    	//--------------------------From ana_Lambda.cc----------------------------------------//


	  //==== tree_out for mixed event analysis =======//
	 	 tr.chi2_l = L_tr_chi2[lt];
	 	 tr.x_l    = L_tr_x[lt];
	 	 tr.y_l    = L_tr_y[lt];
	 	 tr.th_l   = L_tr_th[lt];
	 	 tr.ph_l   = L_tr_ph[lt];
	 	 tr.mom_l  = L_tr_p[lt];
	 	 tr.tg_th_l= L_tr_tg_th[lt];
	 	 tr.tg_ph_l= L_tr_tg_ph[lt];
	 	 tr.vz_l   = L_tr_vz[lt];

	 	 tr.chi2_r = R_tr_chi2[rt];
	 	 tr.x_r    = R_tr_x[rt];
	 	 tr.y_r    = R_tr_y[rt];
	 	 tr.th_r   = R_tr_th[rt];
	 	 tr.ph_r   = R_tr_ph[rt];
	 	 tr.mom_r  = R_tr_p[rt];
	 	 tr.tg_th_r= R_tr_tg_th[rt];
	 	 tr.tg_ph_r= R_tr_tg_ph[rt];
	 	 tr.vz_r   = R_tr_vz[rt];
	
         tr.missing_mass = mm;
	     tr.momR         = R_tr_p[0];
		 tr.momL         = L_tr_p[0];
	     tr.zR           = R_tr_vz[0] ; 
		 tr.zL           = L_tr_vz[0];
	     tr.AC1_sum      = R_a1_asum_p ;
         tr.AC2_sum      = R_a2_asum_p;
	     tr.ct_acc		 = ctime;


	     tree_out->Fill();//Entry <--> Coin. Tracking

          } // if L_Tr && L_FP && R_Tr && R_FP
        } // for NRtr
      } // for NLtr
	  //tree_out->Fill();//Entry <--> Single Tracking
    } // if LHRS && RHRS
    if(k%100000==0){
      end = time(NULL);
      time(&end);
      double diff = difftime(end,start);
      double esttime = diff * ENum / (k+1) - diff;
      cout<<k<<" / "<<ENum<<" : "<<Form("%.0lf sec passed,  %.0lf sec left",diff,esttime)<<endl;
    }
  } // for ENum


}//Filling()



////////////////////////////////////////////////////////////


void tuning::Fitting(){
}//Fitting
//////////////////////////////////////////////////////////////////

void tuning::ACtune(){ 
  cout<<"==========================="<<endl;
  cout<<"========= Tuning =========="<<endl;
  cout<<"==========================="<<endl;

 def_sig_L=0.003; def_mean_L=0.0;
 def_sig_S=0.004; def_mean_S=MS0-ML;
 def_sig_p=0.852; def_mean_p=-8.0;
 def_sig_pi=0.443; def_mean_pi=3.0;
 def_sig_k=0.644; def_mean_k=0.0;
 def_acc=27.7;


//-----No AC Cut-----//
 hcoin_bg_fom_noZ->Scale(40./80.);
 hcoin_wo_bg_fom_noZ->Add(hcoin_k_fom_noZ,hcoin_bg_fom_noZ,1.0,-1.0);
 fp_noZ=new TF1("fp_noZ","gausn(0)",min_coin_c,max_coin_c);
 fp_noZ->SetNpx(2000);
 fpi_noZ =new TF1("fpi_noZ","gausn(0)+gausn(3)",def_mean_pi-3*def_sig_pi,def_mean_pi+6*def_sig_pi);
 fpi_noZ->SetNpx(2000);
 fk_noZ=new TF1("fk_noZ","gausn(0)",min_coin_c,max_coin_c);
 fk_noZ->SetNpx(2000);
//cout<<"fp fit start"<<endl;
 hcoin_wo_bg_fom_noZ->Fit("fp_noZ","Rq0","0",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
// n_p_noZ=fp_noZ->GetParameter(0);
 mean_p_noZ=fp_noZ->GetParameter(1);
 sig_p_noZ=fp_noZ->GetParameter(2);
 //center_p=mean_p_noZ;
 //range_p=2*sig_p_noZ;
 //range_p=1.435;
 center_p=def_mean_p;
 range_p=2*def_sig_p;
 n_p_noZ=hcoin_wo_bg_fom_noZ->Integral(hcoin_wo_bg_fom_noZ->FindBin(center_p-range_p),hcoin_wo_bg_fom_noZ->FindBin(center_p+range_p));
//cout<<"fpi fit start"<<endl;
 fpi_noZ->SetParameters(10000.,def_mean_pi,def_sig_pi,2000.,def_mean_pi,def_sig_pi);
 //fpi_noZ->SetParameters(1000.,def_mean_pi,def_sig_pi);
 hcoin_wo_bg_fom_noZ->Fit("fpi_noZ","Rq0","0",def_mean_pi-3*def_sig_pi,def_mean_pi+5*def_sig_pi);
 n_pi_noZ=fpi_noZ->Integral(def_mean_pi-5*def_sig_pi,def_mean_pi+5*def_sig_pi);
 n_pi_noZ=n_pi_noZ/((max_coin_c-min_coin_c)/bin_coin_c);
 mean_pi_noZ=fpi_noZ->GetParameter(1);
 sig_pi_noZ=fpi_noZ->GetParameter(2);
 hcoin_pi_noZ->FillRandom("fpi_noZ",n_pi_noZ*1000.);
 hcoin_pi_noZ->Scale(1./1000.);
 hcoin_wo_bg_fom_noZ2->Add(hcoin_wo_bg_fom_noZ,hcoin_pi_noZ,1.,-1.);
 //center_pi=mean_pi_noZ;
 //range_pi=2*sig_pi_noZ;
 //range_pi=0.7;
 center_pi=def_mean_pi;
 range_pi=2*def_sig_pi;
 n_pi_noZ=hcoin_wo_bg_fom_noZ->Integral(hcoin_wo_bg_fom_noZ->FindBin(center_pi-range_pi),hcoin_wo_bg_fom_noZ->FindBin(center_pi+range_pi));
//cout<<"fk fit start"<<endl;
 //hcoin_wo_bg_fom_noZ->Fit("fk_noZ","Rq0","0",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_wo_bg_fom_noZ->Fit("fk_noZ","Rq0","0",-1,1);
 //n_k_noZ=fk_noZ->GetParameter(0);
 mean_k_noZ=fk_noZ->GetParameter(1);
// mean_k_noZ=0;
 sig_k_noZ=fk_noZ->GetParameter(2);
 //center_k=mean_k_noZ;
 center_k=0.;
 //range_k=2*sig_k_noZ;
 range_k=1.;

 n_k_noZ=hcoin_wo_bg_fom_noZ->Integral(hcoin_wo_bg_fom_noZ->FindBin(center_k-range_k),hcoin_wo_bg_fom_noZ->FindBin(center_k+range_k));

 //-----Fitting as a whole function---------//

 fcoin_noZ =new TF1("fcoin_noZ","gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin_noZ->SetNpx(2000);
 fcoin_noZ->SetTitle("Cointime w/o AC cut;Cointime [ns];Counts [1/56 ns]");
 fcoin_noZ->SetParameters(n_pi_noZ,mean_pi_noZ,sig_pi_noZ,n_k_noZ,mean_k_noZ,sig_k_noZ,n_p_noZ,mean_p_noZ,sig_p_noZ);
 //hcoin_wo_bg_fom_noZ->Fit("fcoin_woAC","Rq0","0",min_coin_c,max_coin_c);
 //n_pi_noZ=fcoin_noZ->GetParameter(0);//Npi_nocut
 //mean_pi_noZ=fcoin_noZ->GetParameter(1);
 //sig_pi_noZ=fcoin_noZ->GetParameter(2);
 //n_k_noZ=fcoin_noZ->GetParameter(3);//Nk_nocut
 //mean_k_noZ=fcoin_noZ->GetParameter(4);
 //sig_k_noZ=fcoin_noZ->GetParameter(5);
 //n_p_noZ=fcoin_noZ->GetParameter(6);//Np_nocut
 //mean_p_noZ=fcoin_noZ->GetParameter(7);
 //sig_p_noZ=fcoin_noZ->GetParameter(8);
cout<<"n_pi_noZ="<<n_pi_noZ<<"n_k_noZ="<<n_k_noZ<<"n_p_noZ="<<n_p_noZ
<<"mean_pi_noZ="<<mean_pi_noZ<<"sig_pi_noZ="<<sig_pi_noZ<<"mean_k_noZ="<<mean_k_noZ<<"sig_k_noZ="<<sig_k_noZ<<"mean_p_noZ="<<mean_p_noZ<<"sig_p_noZ="<<sig_p_noZ<<endl;

 //------- Get Error Paramters ---//
 n_pi_err_noZ=fpi_noZ->GetParError(0); 
 n_k_err_noZ=fk_noZ->GetParError(3); 
 n_p_err_noZ=fp_noZ->GetParError(6); 
 

//----------------------------------------------//
//--	Missing Mass  Start     ----------------//
//----------------------------------------------//
 hmm_bg_fom_noZ->Scale(2./80.);
 hmm_pibg_fom_noZ->Scale(2.0/80.);
 hmm_wo_bg_fom_noZ->Add(hmm_L_fom_noZ,hmm_bg_fom_noZ,1.0,-1.0);
 hmm_pi_wobg_fom_noZ->Add(hmm_pi_fom_noZ,hmm_pibg_fom_noZ,1.0,-1.0);


// fmmbg_noZ=new TF1("fmmbg_noZ","gausn(0)+gausn(3)",min_mm,max_mm);
 //fmmbg_noZ=new TF1("fmmbg_noZ",F_Voigt,min_mm,max_mm,4);
 fmmbg_noZ=new TF1("fmmbg_noZ","pol4",-0.05,0.15);
// fmmbg_noZ->SetParameters(5,0.05,0.05,0.01);
// fmmbg_noZ->SetParLimits(0,0.,100000.);//positive
// fmmbg_noZ->SetParLimits(3,0.,100.);//positive
 fmmbg_noZ->SetNpx(2000);
 //fmmbg_noZ->SetParameters(100,0.05,0.03,10,0.05,0.03);//test.list
// fmmbg_noZ->SetParameters(300,0.05,1.2,30,0.1,0.02);//small.list
// fmmbg_noZ->SetParameter(1,0.05);
// fmmbg_noZ->SetParameter(2,0.03);
 fL_noZ=new TF1("fL_noZ","gausn(0)",min_mm,max_mm);
 fL_noZ->SetNpx(2000);
 fL_noZ->SetParLimits(2,0.,0.01);
 fS_noZ=new TF1("fS_noZ","gausn(0)",min_mm,max_mm);
 fS_noZ->SetNpx(2000);
 fS_noZ->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
 
 fmm_noZ=new TF1("fmm_noZ","gausn(0)+gausn(3)+pol4(6)",-0.05,0.15);
 fmm_noZ->SetNpx(2000);
 fmm_noZ->SetTitle("Missing Mass w/o AC cut;Coin time [ns];Counts [1/56 ns]");
 fmm_noZ->SetParLimits(0,0.,1000000.);//positive
 fmm_noZ->SetParLimits(3,0.,1000000.);//positive
// fmm_noZ->SetParLimits(3,0.,100.);//positive//Voigt
// fmm_noZ->SetParLimits(4,0.,100000.);//positive//Voigt
// fmm_noZ->SetParLimits(7,0.,100000.);//positive//Voigt
// fmm_noZ->SetParameter(1,def_mean_L);
// fmm_noZ->SetParameter(4,def_mean_S);

 hmm_wo_bg_fom_noZ->Fit("fL_noZ","Rq0","0",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
 mean_L_noZ=fL_noZ->GetParameter(1);
 sig_L_noZ=fL_noZ->GetParameter(2);
 center_L=def_mean_L;
 range_L=2*def_sig_L;

 hmm_wo_bg_fom_noZ->Fit("fS_noZ","Rq0","0",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
 mean_S_noZ=fS_noZ->GetParameter(1);
 sig_S_noZ=fS_noZ->GetParameter(2);
 center_S=def_mean_S;
 range_S=2*def_sig_S;


 //------- Fitting ----------//

cout<<"fmmbg fit start"<<endl;
// hmm_bg_fom_noZ->Fit("fmmbg_noZ","Rq0","",min_mm,max_mm);
// hmm_pi_wobg_fom_noZ->Fit("fmmbg_noZ","Rq0","",-0.05,0.15);
//double fmmbga = fmmbg_noZ->GetParameter(0);
//double fmmbgb = fmmbg_noZ->GetParameter(1);
//double fmmbgc = fmmbg_noZ->GetParameter(2);
//double fmmbgd = fmmbg_noZ->GetParameter(3);
//cout<<"0:1:2:3:4:5="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<endl;//"::"<<e<<"::"<<f<<endl;
// fmm_noZ->FixParameter(1,b);//mean
// fmm_noZ->FixParameter(2,c);//sigma
// fmm_noZ->FixParameter(3,d);//lg
//double e = fmmbg_noZ->GetParameter(4);
//double f = fmmbg_noZ->GetParameter(5);
// fmm_noZ->SetParameter(0,a);
// fmm_noZ->SetParameter(1,b);
// fmm_noZ->SetParameter(2,c);
// fmm_noZ->SetParameter(3,d);
// fmm_noZ->SetParameter(4,e);
// fmm_noZ->SetParameter(5,f);
// fmm_noZ->SetParameter(6,500);
 fmm_noZ->SetParameter(1,def_mean_L);
 fmm_noZ->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
 fmm_noZ->SetParameter(2,def_sig_L);
 fmm_noZ->SetParLimits(2,0.,2*def_sig_L);
// fmm_noZ->SetParameters(9,100);
 fmm_noZ->SetParameter(4,def_mean_S);
 fmm_noZ->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
 fmm_noZ->SetParameter(5,def_sig_S);
 fmm_noZ->SetParLimits(5,0.,2*def_sig_S);
 hmm_wo_bg_fom_noZ->Fit("fmm_noZ","Rq0","0",-0.05,0.15);
double fmm_noZpar0 = fmm_noZ->GetParameter(0);cout<<"fmm_noZ[0]="<<fmm_noZpar0<<endl;//area(L)
double fmm_noZpar1 = fmm_noZ->GetParameter(1);cout<<"fmm_noZ[1]="<<fmm_noZpar1<<endl;//mean(L)
double fmm_noZpar2 = fmm_noZ->GetParameter(2);cout<<"fmm_noZ[2]="<<fmm_noZpar2<<endl;//sigma(L)
double fmm_noZpar3 = fmm_noZ->GetParameter(3);cout<<"fmm_noZ[3]="<<fmm_noZpar3<<endl;//area(S)
double fmm_noZpar4 = fmm_noZ->GetParameter(4);cout<<"fmm_noZ[4]="<<fmm_noZpar4<<endl;//mean(S)
double fmm_noZpar5 = fmm_noZ->GetParameter(5);cout<<"fmm_noZ[5]="<<fmm_noZpar5<<endl;//sigma(S)
double fmm_noZpar6 = fmm_noZ->GetParameter(6);cout<<"fmm_noZ[6]="<<fmm_noZpar6<<endl;//poly_const
double fmm_noZpar7 = fmm_noZ->GetParameter(7);cout<<"fmm_noZ[7]="<<fmm_noZpar7<<endl;//poly_x
double fmm_noZpar8 = fmm_noZ->GetParameter(8);cout<<"fmm_noZ[8]="<<fmm_noZpar8<<endl;//poly_x^2
double fmm_noZpar9 = fmm_noZ->GetParameter(9);cout<<"fmm_noZ[9]="<<fmm_noZpar9<<endl;//poly_x^3
double fmm_noZpar10 = fmm_noZ->GetParameter(10);cout<<"fmm_noZ[10]="<<fmm_noZpar10<<endl;//poly_x^4
 fmmbg_noZ->SetParameters(fmm_noZpar6,fmm_noZpar7,fmm_noZpar8,fmm_noZpar9,fmm_noZpar10);
//cout<<"0:1:2:3:4:5(as a total func)="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<"::"<<e<<"::"<<f<<endl;
//
//cout<<"fL fit start"<<endl;
// hmm_wo_bg_fom_noZ->Fit("fL_noZ","Rq0","",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
// n_L_noZ=fL_noZ->GetParameter(0);
// mean_L_noZ=fL_noZ->GetParameter(1);
// sig_L_noZ=fL_noZ->GetParameter(2);
 mean_L_noZ=def_mean_L;
 mean_S_noZ=def_mean_S;
 sig_L_noZ=def_sig_L;
 sig_S_noZ=def_sig_S;
 n_L_noZ=hmm_wo_bg_fom_noZ->Integral(hmm_wo_bg_fom_noZ->FindBin(center_L-range_L),hmm_wo_bg_fom_noZ->FindBin(center_L+range_L));
cout<<"before(L):: "<<n_L_noZ<<endl;
 double integralL=fmmbg_noZ->Integral(center_L-range_L,center_L+range_L);
 integralL=integralL/(2*range_L/(hmm_wo_bg_fom_noZ->FindBin(center_L+range_L)-hmm_wo_bg_fom_noZ->FindBin(center_L-range_L)));
cout<<"integralL="<<integralL<<endl;
 if(integralL>0)n_L_noZ=n_L_noZ-integralL;
cout<<"after(L):: "<<n_L_noZ<<endl;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,5)-pow(mean_L_noZ-2*sig_L_noZ,5))*a/5;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,4)-pow(mean_L_noZ-2*sig_L_noZ,4))*b/4;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,3)-pow(mean_L_noZ-2*sig_L_noZ,3))*c/3;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,2)-pow(mean_L_noZ-2*sig_L_noZ,2))*d/2;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,1)-pow(mean_L_noZ-2*sig_L_noZ,1))*e;
//
 n_S_noZ=hmm_wo_bg_fom_noZ->Integral(hmm_wo_bg_fom_noZ->FindBin(center_S-range_S),hmm_wo_bg_fom_noZ->FindBin(center_S+range_S));
cout<<"before(S):: "<<n_S_noZ<<endl;
 double integralS=fmmbg_noZ->Integral(center_S-range_S,center_S+range_S);
 integralS=integralS/(2*range_S/(hmm_wo_bg_fom_noZ->FindBin(center_S+range_S)-hmm_wo_bg_fom_noZ->FindBin(center_S-range_S)));
cout<<"integralS="<<integralS<<endl;
 if(integralS>0)n_S_noZ-=integralS;
cout<<"after(S):: "<<n_S_noZ<<endl;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,5)-pow(mean_S_noZ-2*sig_S_noZ,5))*a/5;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,4)-pow(mean_S_noZ-2*sig_S_noZ,4))*b/4;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,3)-pow(mean_S_noZ-2*sig_S_noZ,3))*c/3;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,2)-pow(mean_S_noZ-2*sig_S_noZ,2))*d/2;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,1)-pow(mean_S_noZ-2*sig_S_noZ,1))*e;

 cout<<"n_L"<<n_L_noZ<<endl;
cout<<"mean_L"<<mean_L_noZ<<endl;
cout<<"sig_L"<<sig_L_noZ<<endl;

////cout<<"fS fit start"<<endl;
// hmm_wo_bg_fom_noZ->Fit("fS_noZ","Rq0","",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
//// n_S_noZ=fS_noZ->GetParameter(0);
// mean_S_noZ=fS_noZ->GetParameter(1);
// mean_S_noZ=def_mean_S;
// sig_S_noZ=def_sig_S;
// sig_S_noZ=fS_noZ->GetParameter(2);
// n_S_noZ=hmm_wo_bg_fom_noZ->Integral(hmm_wo_bg_fom_noZ->FindBin(mean_S_noZ-2*sig_S_noZ),hmm_wo_bg_fom_noZ->FindBin(mean_S_noZ+2*sig_S_noZ));
 cout<<"n_S"<<n_S_noZ<<endl;
cout<<"mean_S"<<mean_S_noZ<<endl;
cout<<"sig_S"<<sig_S_noZ<<endl;
//----------------------------------------------//
//--	Missing Mass  End     ------------------//
//----------------------------------------------//
 //----------------------------------------------------------------//
 //----------------------------------------------------------------//
 //----------------------------------------------------------------//

	for(int i=0;i<nth;i++){
//		for(int j=0;j<nth;j++){
//			for(int l=0;l<nth;l++){

			int j=0; int l=0;
			
//-----Background subtraction-----//
//------------cointime------------//
//cout<<"BG subtraction cointime"<<endl;
 hcoin_bg_fom[i][j][l]->Scale(40./80.);
 hcoin_wo_bg_fom[i][j][l]->Add(hcoin_k_fom[i][j][l],hcoin_bg_fom[i][j][l],1.0,-1.0);

 fp[i][j][l]=new TF1(Form("fp[%d][%d][%d]",i,j,l),"gausn(0)",min_coin_c,max_coin_c);
 fp[i][j][l]->SetNpx(2000);
 fpi[i][j][l] =new TF1(Form("fpi[%d][%d][%d]",i,j,l),"gausn(0)+gausn(3)",min_coin_c,max_coin_c);
 fpi[i][j][l]->SetNpx(2000);
 fk[i][j][l]=new TF1(Form("fk[%d][%d][%d]",i,j,l),"gausn(0)",min_coin_c,max_coin_c);
 fk[i][j][l]->SetNpx(2000);
//cout<<"fp fit start"<<endl;
 hcoin_wo_bg_fom[i][j][l]->Fit(Form("fp[%d][%d][%d]",i,j,l),"Rq0","0",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 //n_p[i][j][l]=fp[i][j][l]->GetParameter(0);
 mean_p[i][j][l]=fp[i][j][l]->GetParameter(1);
 sig_p[i][j][l]=fp[i][j][l]->GetParameter(2);
 n_p[i][j][l]=hcoin_wo_bg_fom[i][j][l]->Integral(hcoin_wo_bg_fom[i][j][l]->FindBin(center_p-range_p),hcoin_wo_bg_fom[i][j][l]->FindBin(center_p+range_p));
//cout<<"fpi fit start"<<endl;
 fpi[i][j][l]->SetParameters(10000.,def_mean_pi,def_sig_pi,2000.,def_mean_pi,def_sig_pi);
 hcoin_wo_bg_fom[i][j][l]->Fit(Form("fpi[%d][%d][%d]",i,j,l),"Rq0","0",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 n_pi[i][j][l]=fpi[i][j][l]->Integral(min_coin_c,max_coin_c);
 n_pi[i][j][l]=n_pi[i][j][l]/((max_coin_c-min_coin_c)/bin_coin_c);
 //n_pi[i][j][l]=n_pi[i][j][l]/(2*range_pi/(hcoin_wo_bg_fom[i][j][l]->FindBin(center_pi+range_pi)-hcoin_wo_bg_fom[i][j][l]->FindBin(center_pi-range_pi)));
 mean_pi[i][j][l]=fpi[i][j][l]->GetParameter(1);
cout << mean_pi[i][j][l] << endl;
 sig_pi[i][j][l]=fpi[i][j][l]->GetParameter(2);
 hcoin_pi[i][j][l]->FillRandom(Form("fpi[%d][%d][%d]",i,j,l),n_pi[i][j][l]);
 //hcoin_wo_bg_fom[i][j][l]->Add(hcoin_wo_bg_fom[i][j][l],hcoin_pi[i][j][l],1.,-1.);
n_pi[i][j][l]=hcoin_wo_bg_fom[i][j][l]->Integral(hcoin_wo_bg_fom[i][j][l]->FindBin(center_pi-range_pi),hcoin_wo_bg_fom[i][j][l]->FindBin(center_pi+range_pi));
//cout<<"fk fit start"<<endl;
// hcoin_wo_bg_fom[i][j][l]->Fit(Form("fk[%d][%d][%d]",i,j,l),"Rq0","0",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_wo_bg_fom[i][j][l]->Fit(Form("fk[%d][%d][%d]",i,j,l),"Rq0","0",-1,1);
// n_k[i][j][l]=fk[i][j][l]->GetParameter(0);
 mean_k[i][j][l]=fk[i][j][l]->GetParameter(1);
 sig_k[i][j][l]=fk[i][j][l]->GetParameter(2);
 n_k[i][j][l]=hcoin_wo_bg_fom[i][j][l]->Integral(hcoin_wo_bg_fom[i][j][l]->FindBin(center_k-range_k),hcoin_wo_bg_fom[i][j][l]->FindBin(center_k+range_k));

 // n_k[i][th1][1]=hcoin_k_fom[i][th1]->Integral(hcoin_k_fom[i][th1]->FindBin(-3*sig_k[i][th1][1]+mean_k[i][th1][1])
 //		     ,hcoin_k_fom[i][th1]->FindBin(+3*sig_k[i][th1][1]+mean_k[i][th1][1]));
 fcoin[i][j][l] =new TF1(Form("fcoin[%d][%d][%d]",i,j,l),"gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin[i][j][l]->SetNpx(2000);
 fcoin[i][j][l]->SetTitle(Form("Coin-Time w AC cut (AC1 Cut <%lf ch && %lf ch < AC2 Cut < %lf ch);Coin time [ns];Counts [1/56 ns]",ac1_adc[i],ac2l_adc[j],ac2u_adc[l]));
 fcoin[i][j][l]->SetParameters(n_pi[i][j][l],mean_pi[i][j][l],sig_pi[i][j][l],n_k[i][j][l],mean_k[i][j][l],sig_k[i][j][l],n_p[i][j][l],mean_p[i][j][l],sig_p[i][j][l]);
 //hcoin_wo_bg_fom[i][j][l]->Fit(Form("fcoin[%d][%d][%d]",i,j,l),"Rq0","0",min_coin_c,max_coin_c);
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
 hmm_bg_fom[i][j][l]->Scale(2./80.);
 hmm_pibg_fom[i][j][l]->Scale(2.0/80.);
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
// hmm_bg_fom[i][j][l]->Fit("fmmbg[i][j][l]","Rq0","",min_mm,max_mm);
// hmm_pi_wobg_fom[i][j][l]->Fit("fmmbg[i][j][l]","Rq0","",-0.05,0.15);
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
 hmm_wo_bg_fom[i][j][l]->Fit("fmm[i][j][l]","Rq0","0",-0.05,0.15);
//double fmmpar0 = fmm[i][j][l]->GetParameter(0);//cout<<"fmm[0]="<<fmmpar0<<endl;//area(L)
//double fmmpar1 = fmm[i][j][l]->GetParameter(1);//cout<<"fmm[1]="<<fmmpar1<<endl;//mean(L)
//double fmmpar2 = fmm[i][j][l]->GetParameter(2);//cout<<"fmm[2]="<<fmmpar2<<endl;//sigma(L)
//double fmmpar3 = fmm[i][j][l]->GetParameter(3);//cout<<"fmm[3]="<<fmmpar3<<endl;//area(S)
//double fmmpar4 = fmm[i][j][l]->GetParameter(4);//cout<<"fmm[4]="<<fmmpar4<<endl;//mean(S)
//double fmmpar5 = fmm[i][j][l]->GetParameter(5);//cout<<"fmm[5]="<<fmmpar5<<endl;//sigma(S)
mean_L[i][j][l]=fmm[i][j][l]->GetParameter(1);
sig_L[i][j][l]=fmm[i][j][l]->GetParameter(2);
mean_S[i][j][l]=fmm[i][j][l]->GetParameter(4);
sig_S[i][j][l]=fmm[i][j][l]->GetParameter(5);
double fmmpar6 = fmm[i][j][l]->GetParameter(6);//cout<<"fmm[6]="<<fmmpar6<<endl;//poly_const
double fmmpar7 = fmm[i][j][l]->GetParameter(7);//cout<<"fmm[7]="<<fmmpar7<<endl;//poly_x
double fmmpar8 = fmm[i][j][l]->GetParameter(8);//cout<<"fmm[8]="<<fmmpar8<<endl;//poly_x^2
double fmmpar9 = fmm[i][j][l]->GetParameter(9);//cout<<"fmm[9]="<<fmmpar9<<endl;//poly_x^3
double fmmpar10 = fmm[i][j][l]->GetParameter(10);//cout<<"fmm[10]="<<fmmpar10<<endl;//poly_x^4
 fmmbg[i][j][l]->SetParameters(fmmpar6,fmmpar7,fmmpar8,fmmpar9,fmmpar10);
//cout<<"0:1:2:3:4:5(as a total func)="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<"::"<<e<<"::"<<f<<endl;
//
//cout<<"fL fit start"<<endl;
// hmm_wo_bg_fom[i][j][l]->Fit("fL[i][j][l]","Rq0","",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
// n_L[i][j][l]=fL[i][j][l]->GetParameter(0);
// mean_L[i][j][l]=fL[i][j][l]->GetParameter(1);
// sig_L[i][j][l]=fL[i][j][l]->GetParameter(2);
 n_L[i][j][l]=hmm_wo_bg_fom[i][j][l]->Integral(hmm_wo_bg_fom[i][j][l]->FindBin(center_L-range_L),hmm_wo_bg_fom[i][j][l]->FindBin(center_L+range_L));
cout<<"before(L):: "<<n_L[i][j][l]<<endl;
 double integralL=fmmbg[i][j][l]->Integral(center_L-range_L,center_L+range_L);
 integralL=integralL/(2*range_L/(hmm_wo_bg_fom[i][j][l]->FindBin(center_L+range_L)-hmm_wo_bg_fom[i][j][l]->FindBin(center_L-range_L)));
cout<<"integralL="<<integralL<<endl;
 if(integralL>0)n_L[i][j][l]=n_L[i][j][l]-integralL;
else cout<<"negative BG: ("<<i<<","<<j<<","<<l<<")"<<endl;
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
 integralS=integralS/(2*range_S/(hmm_wo_bg_fom[i][j][l]->FindBin(center_S+range_S)-hmm_wo_bg_fom[i][j][l]->FindBin(center_S-range_S)));
cout<<"integralS="<<integralS<<endl;
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
// hmm_wo_bg_fom[i][j][l]->Fit("fS[i][j][l]","Rq0","",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
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

	
	
	if(n_pi[i][j][l]>0.){}else{n_pi[i][j][l]=1.;}
	if(n_k[i][j][l]>0.){}else{n_k[i][j][l]=1.;}
	if(n_p[i][j][l]>0.){}else{n_p[i][j][l]=1.;}
	if(n_L[i][j][l]>0.){}else{n_L[i][j][l]=1.;}
	if(n_S[i][j][l]>0.){}else{n_S[i][j][l]=1.;}









	

//			}//for l
//		}//for j
	}//for i

//ofstream fout(Form("SR_z_%d.dat",(int)zver[0]));
ofstream fout("SR_z.dat");
cout<<"AC Efficiency is filled."<<endl;
fout<<"i=10: pL<2.12, Full Range (Sigma)"<<endl;
fout<<"i=20: 2.1<pL<2.15, 6msr cut, Top(Lambda)"<<endl;
fout<<"i=30: 2.05<pL<2.1, 6msr cut, Top(Sigma)"<<endl;
fout<<"i=31,32: pL is divide into 2, (w/o theta_ee,phi_ee cut)"<<endl;
fout<<"i=33,34,35,36: pL is divide into 4, (w/o theta_ee,phi_ee cut)"<<endl;
fout<<"i=37,38,39,40,41: pL is divide into 5, (w/o theta_ee,phi_ee cut)"<<endl;
fout<<"i=51,52,53,54: (theta,phi) center shift, (2.1<pL<2.5)"<<endl;

fout<<n_pi_noZ<<" "<<n_k_noZ<<" "<<n_p_noZ<<" "<<n_L_noZ<<" "<<n_S_noZ<<endl;
for(int i=0;i<nth;i++){
//	for(int j=0;j<nth;j++){
		int j=0;	int l=0;
		//fout<<n_pi[i][j][l]/n_pi_noZ<<" "<<n_k[i][j][l]/n_k_noZ<<" "<<n_p[i][j][l]/n_p_noZ<<" "<<n_L[i][j][l]/n_L_noZ<<" "<<n_S[i][j][l]/n_S_noZ<<"|Ave(z)|<"<<zver[i] <<endl;
		fout<<n_pi[i][j][l]<<" "<<n_k[i][j][l]<<" "<<n_p[i][j][l]<<" "<<n_L[i][j][l]<<" "<<n_S[i][j][l]<<"|Ave(z)|<"<<zver[i] <<endl;
//	}
}
		




}//ACtune()

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
c20 = new TCanvas("c20","c20",800.,800.);
c21 = new TCanvas("c21","c21",800.,800.);
c22 = new TCanvas("c22","c22",800.,800.);
c23 = new TCanvas("c23","c23",800.,800.);
c24 = new TCanvas("c24","c24",800.,800.);
c25 = new TCanvas("c25","c25",800.,800.);
c26 = new TCanvas("c26","c26",800.,800.);
c27 = new TCanvas("c27","c27",800.,800.);
c28 = new TCanvas("c28","c28",800.,800.);
c29 = new TCanvas("c29","c29",800.,800.);
c30 = new TCanvas("c30","c30",800.,800.);
c31 = new TCanvas("c31","c31",800.,800.);
c32 = new TCanvas("c32","c32",800.,800.);
c33 = new TCanvas("c33","c33",800.,800.);
c34 = new TCanvas("c34","c34",800.,800.);
c35 = new TCanvas("c35","c35",800.,800.);
c36 = new TCanvas("c36","c36",800.,800.);
c37 = new TCanvas("c37","c37",800.,800.);
c38 = new TCanvas("c38","c38",800.,800.);
c39 = new TCanvas("c39","c39",800.,800.);

c1->Divide(3,2);
c1->cd(1);
hcoin_wo_bg_fom_noZ->Draw("");
c1->cd(2);
hcoin_wo_bg_fom_noZ2->Draw("");
//c1->cd(2);hcoin_wo_bg_fom[20][0][0]->Draw("");
//hcoin_pi[20][0][0]->SetLineColor(kOrange);
//hcoin_pi[20][0][0]->Draw("same");

//cout << "start hcoin_tc" << endl;
c1->cd(3);hcoin_wo_bg_fom_noZ->Draw("");
double ymax = (hcoin_wo_bg_fom_noZ->GetBinContent(hcoin_wo_bg_fom_noZ->GetMaximumBin()));
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
c1->cd(4)->DrawFrame(1.,-100.,5.,9000.);
hcoin_pi_noZ->Draw("same");
c1->cd(5)->DrawFrame(1.,-100.,5.,9000.);
fpi_noZ->SetLineColor(kAzure);
fpi_noZ->Draw("same");

c1->cd(6);hmm_L_fom_noZ->Draw("");
    ymax = (hmm_L_fom_noZ->GetBinContent(hmm_L_fom_noZ->GetMaximumBin()));
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

c2->Divide(2,2);
cout<<"c2 start"<<endl;
c2->cd(1);
hmm_wo_bg_fom_noZ->Draw("");
//hmm_wo_bg_fom[20][0][0]->SetLineColor(kRed);
//hmm_wo_bg_fom[20][0][0]->Draw("same");
c2->cd(2);
hcoin_wo_bg_fom_noZ->Draw("");
//hcoin_wo_bg_fom[20][0][0]->SetLineColor(kRed);
//hcoin_wo_bg_fom[20][0][0]->Draw("same");
c2->cd(3);
hmm_wo_bg_fom_noZ->Draw("");
//hmm_wo_bg_fom[20][0][0]->SetLineColor(kRed);
//hmm_wo_bg_fom[20][0][0]->Draw("same");
hmm_wo_bg_fom[99][0][0]->SetLineColor(kGreen);
hmm_wo_bg_fom[99][0][0]->Draw("same");
c2->cd(4);
hcoin_wo_bg_fom_noZ->Draw("");
hcoin_wo_bg_fom[20][0][0]->SetLineColor(kRed);
hcoin_wo_bg_fom[20][0][0]->Draw("same");
hcoin_wo_bg_fom[99][0][0]->SetLineColor(kGreen);
hcoin_wo_bg_fom[99][0][0]->Draw("same");
c3->Divide(2,2);
c3->cd(1);
hcoin_k_fom_noZ->Draw("");
c3->cd(2);
hmm_L_fom_noZ->Draw("");
c3->cd(3);
hcoin_k_fom_noZ->Draw("");
hcoin_k_fom_Zdiff->SetLineColor(kAzure);hcoin_k_fom_Zdiff->Draw("same");
hcoin_k_fom_Zsum->SetLineColor(kRed);hcoin_k_fom_Zsum->Draw("same");
c3->cd(4);
hmm_L_fom_noZ->Draw("");
hmm_L_fom_Zdiff->SetLineColor(kAzure);hmm_L_fom_Zdiff->Draw("same");
hmm_L_fom_Zsum->SetLineColor(kRed);hmm_L_fom_Zsum->Draw("same");
c4->Divide(2,2);
cout<<"c4 start"<<endl;
c4->cd(3)->SetLogz(1);
h_zz->Draw("colz");
c4->cd(4);
h_m2_ac->Draw("colz");
c5->Divide(2,2);
c5->cd(1);
hct_test2->Draw("");
c5->cd(2);
hct_test->Draw("");
c5->cd(3);
h_ctct->Draw("colz");
c5->cd(4)->SetLogz(1);
h_ctct->Draw("colz");

c8->Divide(2,2);
c8->cd(1);
h_zz1->Draw("colz");
c8->cd(2);
h_zz2->Draw("colz");
c8->cd(3);
h_zz3->Draw("colz");
c8->cd(4);
h_zz4->Draw("colz");
c9->Divide(2,2);
c9->cd(1);
h_z1->Draw("");
c9->cd(2);
h_z2->Draw("");
c9->cd(3);
h_z3->Draw("");
c9->cd(4);
h_z4->Draw("");
c10->Divide(2,2);
c10->cd(1);
h_z11->Draw("");
c10->cd(2);
h_z22->Draw("");
c10->cd(3);
h_z33->Draw("");
c10->cd(4);
h_z44->Draw("");
c11->cd()->DrawFrame(0.,0.,0.5,1.2);//K,L,S vs th1//Diff
cout<<"c11 start"<<endl;
pEff4->SetLineColor(kAzure);pEff4->Draw("same");
pEff5->SetLineColor(kCyan);pEff5->Draw("same");
c12->cd()->DrawFrame(0.,0.,0.5,1.2);//K,L,S vs th1//Diff
cout<<"c12 start"<<endl;
pEff2->SetLineColor(kGreen);pEff2->Draw("same");
pEff4->SetLineColor(kAzure);pEff4->Draw("same");
pEff5->SetLineColor(kCyan);pEff5->Draw("same");

c13->Divide(2,2);
c13->cd(1);
//h_theta_ee->Draw("");
h_L_tgph_tgth->Draw("colz");
c13->cd(2);
//h_theta_ee2->Draw("");
h_L_tgph_tgth2->Draw("colz");
c13->cd(3);
h_theta_ee_p->Draw("colz");
c13->cd(4);
h_theta_ee_p2->Draw("colz");

//c14->Divide(2,2);
//c14->cd(1);
//h_theta_ee3->Draw("");
//c14->cd(2);
//h_theta_ee4->Draw("");
//c14->cd(3);
//h_theta_ee_p3->Draw("colz");
//c14->cd(4);
//h_theta_ee_p4->Draw("colz");

c14->cd();
h_qsq->Draw("");

c15->Divide(2,2);
c15->cd(1);
h_phi_ee->Draw("");
c15->cd(2);
h_phi_ee2->Draw("");
c15->cd(3);
h_phi_ee3->Draw("");
c15->cd(4);
h_phi_ee4->Draw("");

c16->Divide(3,2);
c16->cd(1);
h_theta_ek->Draw("");
c16->cd(2);
h_phi_ek->Draw("");
c16->cd(3);
h_theta_g->Draw("");
c16->cd(4);
h_phi_g->Draw("");
c16->cd(5);
h_thph_ek->Draw("colz");
c16->cd(6);
h_thph_g->Draw("colz");

c17->Divide(2,2);
c17->cd(1);
h_theta_gk_lab->Draw("");
c17->cd(2);
h_theta_gk_cm->Draw("");
c17->cd(3);
h_cos_gk_lab->Draw("");
c17->cd(4);
h_cos_gk_cm->Draw("");

c18->Divide(2,2);
c18->cd(1);
h_mom_g->Draw("");
c18->cd(2);
h_pL_pR->Draw("colz");
c18->cd(3);
h_pL_pR2->Draw("colz");
c18->cd(4);
h_pL_pR3->Draw("colz");


c19->Divide(2,2);
c19->cd(1);
h_Lz2->Draw("");
c19->cd(2);
h_Rz2->Draw("");
c19->cd(3)->SetLogy(1);
h_Lz->Scale(1./10.);
h_Lz->Draw("");
h_Lz2->SetLineColor(kAzure);
h_Lz2->Draw("same");
c19->cd(4)->SetLogy(1);
h_Rz->Scale(1./10.);
h_Rz->Draw("");
h_Rz2->SetLineColor(kAzure);
h_Rz2->Draw("same");


c20->Divide(2,2);
c20->cd(1);
hmm_wo_bg_fom_noZ->Draw("");
fmm_noZ->SetLineColor(kRed);
fmm_noZ->Draw("same");
c20->cd(2);
hmm_wo_bg_fom[20][0][0]->Draw("");
fmm[20][0][0]->SetLineColor(kRed);
fmm[20][0][0]->Draw("same");
c20->cd(3);
hmm_wo_bg_fom[10][0][0]->Draw("");
fmm[10][0][0]->SetLineColor(kRed);
fmm[10][0][0]->Draw("same");
c20->cd(4);
hmm_wo_bg_fom[30][0][0]->Draw("");
fmm[30][0][0]->SetLineColor(kRed);
fmm[30][0][0]->Draw("same");
//c16->Divide(2,2);
//c16->cd(1);
//h_vpflux->SetStats(1);
//h_vpflux->Draw("");
//cout<<"h_vpflux mean"<<endl;
//cout<<"Mean(x)="<<h_vpflux->GetMean(1)<<endl;
////cout<<"Mean(y)="<<h_vpflux->GetMean(2)<<endl;
//c16->cd(2);
//cout<<"h_vpflux2 mean"<<endl;
//cout<<"Mean(x)="<<h_vpflux2->GetMean(1)<<endl;
////cout<<"Mean(y)="<<h_vpflux2->GetMean(2)<<endl;
//h_vpflux2->SetStats(1);
//h_vpflux2->Draw("");
//c16->cd(3);
//cout<<"h_vpflux3 mean"<<endl;
//cout<<"Mean(x)="<<h_vpflux3->GetMean(1)<<endl;
////cout<<"Mean(y)="<<h_vpflux3->GetMean(2)<<endl;
//h_vpflux3->SetStats(1);
//h_vpflux3->Draw("");
//c16->cd(4);
//cout<<"h_vpflux4 mean"<<endl;
//cout<<"Mean(x)="<<h_vpflux4->GetMean(1)<<endl;
////cout<<"Mean(y)="<<h_vpflux4->GetMean(2)<<endl;
//h_vpflux4->SetStats(1);
//h_vpflux4->Draw("");

//c17->Divide(2,2);
//c17->cd(1)->SetLogz(1);
//TProfile *pfx = h_ct_x->ProfileX();
////h_ct_x->Draw("colz");
//pfx->Draw("");
//c17->cd(2)->SetLogz(1);
//TProfile *pfy = h_ct_y->ProfileX();
//pfy->Draw("");
////h_ct_y->Draw("colz");
//c17->cd(3)->SetLogz(1);
//TProfile *pfth = h_ct_th->ProfileX();
//pfth->Draw("");
////h_ct_th->Draw("colz");
//c17->cd(4)->SetLogz(1);
//TProfile *pfph = h_ct_ph->ProfileX();
////h_ct_ph->Draw("colz");
//pfph->Draw("");
//
//c18->Divide(2,3);
//c18->cd(1);
////h_ct_xy->Draw("colz");
//TProfile *pfxy = h_ct_xy->ProfileY();
//pfxy->Draw("");
//c18->cd(2);
//cout<<"4845 line here"<<endl;
////h_ct_thph->Draw("colz");
//h_ct_thph->FitSlicesX();
////h_ct_thph->Draw("");
//c18->cd(3);
////h_ct_xth->Draw("colz");
//h_ct_xth->FitSlicesX();
////TH1F *h_ct_xth_1=(TH1F*)gDirectory->Get("h_ct_xth_1");
//h_ct_xth->Draw("");
////h_ct_xth_1->Draw("same");
//c18->cd(4);
////h_ct_yph->Draw("colz");
//TProfile *pfyph = h_ct_yph->ProfileY();
//pfyph->Draw("");
//c18->cd(5);
////h_ct_xph->Draw("colz");
//TProfile *pfxph = h_ct_xph->ProfileY();
//pfxph->Draw("");
//c18->cd(6);
////h_ct_yth->Draw("colz");
//TProfile *pfyth = h_ct_yth->ProfileY();
//pfyth->Draw("");
//
//c19->Divide(2,3);
//c19->cd(1)->SetLogz(1);
//h_ct_xy->Draw("colz");
//c19->cd(2)->SetLogz(1);
//h_ct_thph->Draw("colz");
//c19->cd(3)->SetLogz(1);
//h_ct_xth->Draw("colz");
//c19->cd(4)->SetLogz(1);
//h_ct_yph->Draw("colz");
//c19->cd(5)->SetLogz(1);
//h_ct_xph->Draw("colz");
//c19->cd(6)->SetLogz(1);
//h_ct_yth->Draw("colz");
//
//c20->Divide(2,3);
//c20->cd(1)->SetLogz(1);
//h_ct_tgth->Draw("colz");
//c20->cd(2)->SetLogz(1);
//h_ct_tgph->Draw("colz");
//c20->cd(3)->SetLogz(1);
//h_ct_vz->Draw("colz");
//c20->cd(4)->SetLogz(1);
//h_ct_tgthtgph->Draw("colz");
//c20->cd(5)->SetLogz(1);
//h_ct_tgthz->Draw("colz");
//c20->cd(6)->SetLogz(1);
//h_ct_tgphz->Draw("colz");
//
//
c21->Divide(2,2);
c21->cd(1);
h_ct_thph->Draw("colz");
c21->cd(2);
fpi_slice = new TF1("fpi_slice","gausn(0)+pol0(3)",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
fpi_slice->SetNpx(2000);
fpi_slice->SetParameters(1000.,def_mean_pi,def_sig_pi,20.);
//h_ct_thph->FitSlicesY(0,h_ct_thph->FindBin(def_mean_pi-3*def_sig_pi),h_ct_thph->FindBin(def_mean_pi+3*def_sig_pi),20);//more than 20 bins
h_ct_thph->FitSlicesX(fpi_slice);
TH1F *h_ct_thph_0=(TH1F*)gDirectory->Get("h_ct_thph_0");
TH1F *h_ct_thph_1=(TH1F*)gDirectory->Get("h_ct_thph_1");
TH1F *h_ct_thph_2=(TH1F*)gDirectory->Get("h_ct_thph_2");
//h_ct_thph_1->GetYaxis()->SetRangeUser(-0.5,0.5);
h_ct_thph_0->Draw("");
c21->cd(3);
h_ct_thph_1->Draw("");
c21->cd(4);
h_ct_thph_2->Draw("");



c22->Divide(2,2);
c22->cd(1);
h_theta_ee_p3->Draw("colz");
c22->cd(2);
h_phi_ee_p->Draw("colz");
c22->cd(3);
h_original_phth5->Draw("colz");
c22->cd(4);
h_original_phth6->Draw("colz");

c23->Divide(2,2);
c23->cd(1);
h_original_phth7->Draw("colz");
c23->cd(2);
h_original_phth8->Draw("colz");
c23->cd(3);
h_original_phth9->Draw("colz");
c23->cd(4);
h_original_phth10->Draw("colz");

c24->Divide(2,2);
c24->cd(1);
h_original_phcosth->Draw("colz");
c24->cd(2);
h_original_phcosth2->Draw("colz");
//
//c22->cd();
//h_ct_thph->Draw("");
//h_ct_thph_1->Draw("same");
//c14->cd();
//h_R_m2->Draw("");
//c15->cd();
//h_L_tgph_tgth->Draw("colz");
//c17->cd()->DrawFrame(-2.,0.,2.,500.);
//hcoin_wo_bg_fom_noZ->Draw("same");
//hcoin_wo_bg_fom[20][0][0]->SetLineColor(kRed);
//hcoin_wo_bg_fom[20][0][0]->Draw("same");
//hcoin_wo_bg_fom[99][0][0]->SetLineColor(kGreen);
//hcoin_wo_bg_fom[99][0][0]->Draw("same");
//
//c18->cd()->DrawFrame(-0.03,0.,0.03,500.);
//hmm_wo_bg_fom_noZ->Draw("same");
//hmm_wo_bg_fom[20][0][0]->SetLineColor(kRed);
//hmm_wo_bg_fom[20][0][0]->Draw("same");
//hmm_wo_bg_fom[99][0][0]->SetLineColor(kGreen);
//hmm_wo_bg_fom[99][0][0]->Draw("same");
//
//c19->Divide(2,2);
//c19->cd(1);
//hmm_wo_bg_fom[0][0][0]->Draw("");
//fmmbg[0][0][0]->SetLineColor(kGreen);
//fmmbg[0][0][0]->Draw("same");
//c19->cd(2);
//hmm_wo_bg_fom[20][0][0]->Draw("");
//fmmbg[20][0][0]->SetLineColor(kGreen);
//fmmbg[20][0][0]->Draw("same");
//c19->cd(3);
//hmm_wo_bg_fom[50][0][0]->Draw("");
//fmmbg[50][0][0]->SetLineColor(kGreen);
//fmmbg[50][0][0]->Draw("same");
//c19->cd(4);
//hmm_wo_bg_fom[99][0][0]->Draw("");
//fmmbg[99][0][0]->SetLineColor(kGreen);
//fmmbg[99][0][0]->Draw("same");

//c20->Divide(2,3);
//c20->cd(1);
//h_gbetaL->Draw("");
//c20->cd(2);
//h_gbetaR->Draw("");
//c20->cd(3);
//h_gLenL->Draw("");
//c20->cd(4);
//h_gLenR->Draw("");
//c20->cd(5);
//h_gpL->Draw("");
//c20->cd(6);
//h_gpR->Draw("");
//
//c21->Divide(2,2);
//c21->cd(1);
//h_gcorL->Draw("");
//c21->cd(2);
//h_gcorR->Draw("");
//c21->cd(3);
//h_gcorLR->Draw("colz");
//
//c22->Divide(2,3);
//c22->cd(1);
//h_gtref_R->Draw("");
//c22->cd(2);
//h_gmeantime->Draw("");
//c22->cd(3);
//h_gtimeL_R->Draw("");
//c22->cd(4);
//h_gtimeR_R->Draw("");
//c22->cd(5);
//h_gctcorL->Draw("");
//c22->cd(6);
//h_gctcorR->Draw("");
//
//
//c23->Divide(2,3);
//c23->cd(1)->SetLogy(1);
//h_gtref_R->Draw("");
//c23->cd(2)->SetLogy(1);
//h_gmeantime->Draw("");
//c23->cd(3)->SetLogy(1);
//h_gtimeL_R->Draw("");
//c23->cd(4)->SetLogy(1);
//h_gtimeR_R->Draw("");
//c23->cd(5)->SetLogy(1);
//h_gctcorL->Draw("");
//c23->cd(6)->SetLogy(1);
//h_gctcorR->Draw("");
//
//c24->Divide(2,2);
//c24->cd(1);
//h_L_tgph_tgth12->Draw("colz");//No Cut
//c24->cd(2);
//h_L_tgph_tgth5->Draw("colz");//bestcut
//c24->cd(3);
//h_LHRS_phth->Draw("colz");
//c24->cd(4);
//h_LHRS_phth2->Draw("colz");
//
//c25->Divide(2,2);
//c25->cd(1);
//h_original_phth->Draw("colz");
//c25->cd(2);
//h_original_phth2->Draw("colz");
//c25->cd(3);
//h_LHRS_phth->Draw("colz");
//c25->cd(4);
//h_LHRS_phth2->Draw("colz");



c26->cd();
hmm_wo_bg_fom[20][0][0]->Draw("");//top_quality
cout<<"N(Lambda)_top"<<n_L[20][0][0]<<endl;
cout<<"N(Sigma)_top"<<n_S[20][0][0]<<endl;

c27->cd();
hmm_wo_bg_fom_noZ->Draw("");//acceptance
cout<<"N(Lambda)"<<n_L_noZ<<endl;
cout<<"N(Sigma)"<<n_S_noZ<<endl;

c28->Divide(2,2);
c28->cd(1);
h_vpflux->Draw("");//top
c28->cd(2);
h_original_phth3->Draw("colz");//top	
c28->cd(3);
h_Lp_top->Draw("");//Z
c28->cd(4);
h_original_phth4->Draw("colz");

c29->cd();
hmm_wo_bg_fom[20][0][0]->Draw("");//top_quality
fmm[20][0][0]->SetLineColor(kRed);//top_quality
fmm[20][0][0]->Draw("same");//top_quality

c30->cd();
//hmm_wo_bg_fom_noZ->Draw("");//acceptance
//fmm_noZ->SetLineColor(kRed);
//fmm_noZ->Draw("same");
hmm_L_fom_noZ->Draw("");//acceptance
hmm_bg_fom_noZ->SetLineColor(kGreen);
hmm_bg_fom_noZ->Draw("same");//acceptance

c31->Divide(2,2);
c31->cd(1);
h_vpflux->Draw("");//top
c31->cd(2);
h_vpflux2->Draw("");//acc
c31->cd(3);
h_vpflux3->Draw("");//all (Z Only)
cout<<"Ne(det,top)(Entry)"<<h_vpflux->GetEntries()<<endl;
cout<<"Ne(det,top)(Integral)"<<h_vpflux->Integral()<<endl;
cout<<"Ne(det,acc)(Entry)"<<h_vpflux2->GetEntries()<<endl;
cout<<"Ne(det,acc)(Integral)"<<h_vpflux2->Integral()<<endl;
cout<<"Ne(det,all)(Entry)"<<h_vpflux3->GetEntries()<<endl;
cout<<"Ne(det,all)(Integral)"<<h_vpflux3->Integral()<<endl;
cout<<"Ng(det,top)"<<Ng_det_top<<endl;
cout<<"Ng(det,acc)"<<Ng_det_acc<<endl;
cout<<"Ng(det,all)"<<Ng_det<<endl;

//c25->Divide(2,2);
//c25->cd(1);
//h_L_tgph_tgth6->Draw("colz");
//c25->cd(2);
//h_L_tgph_tgth7->Draw("colz");
//c25->cd(3);
//h_L_tgph_tgth8->Draw("colz");
//c25->cd(4);
//h_L_tgph_tgth9->Draw("colz");
//
//c26->Divide(2,2);
//c26->cd(1);
//h_L_tgph_tgth10->Draw("colz");
//c26->cd(2);
//h_L_tgph_tgth11->Draw("colz");
//c26->cd(3);
//h_L_y_x->Draw("colz");
//
//
//c27->Divide(2,2);
//c27->cd(1);
//h_L_z_x->Draw("colz");
//c27->cd(2);
//h_L_z_y->Draw("colz");
//c27->cd(3);
//h_L_th_y->Draw("colz");
//c27->cd(4);
//h_L_ph_y->Draw("colz");
//
//
//c28->Divide(2,2);
//c28->cd(1);
//h_L_th_x->Draw("colz");
//c28->cd(2);
//h_L_ph_x->Draw("colz");
//c28->cd(3);
//h_L_th_z->Draw("colz");
//c28->cd(4);
//h_L_ph_z->Draw("colz");

//cout<<"n_pi_noZ="<<n_pi_noZ<<endl;
//cout<<"n_pi[20][0][0]="<<n_pi[20][0][0]<<endl;
//cout<<"n_pi[99][0][0]="<<n_pi[99][0][0]<<endl;
//cout<<"n_k_noZ="<<n_k_noZ<<endl;
//cout<<"n_k[20][0][0]="<<n_k[20][0][0]<<endl;
//cout<<"n_k[99][0][0]="<<n_k[99][0][0]<<endl;
//cout<<"n_L_noZ="<<n_L_noZ<<endl;
//cout<<"n_L[20][0][0]="<<n_L[20][0][0]<<endl;
//cout<<"n_L[99][0][0]="<<n_L[99][0][0]<<endl;
//cout<<"n_S_noZ="<<n_S_noZ<<endl;
//cout<<"n_S[20][0][0]="<<n_S[20][0][0]<<endl;
//cout<<"n_S[99][0][0]="<<n_S[99][0][0]<<endl;
//
//cout<<"Results in the case of Best Cut"<<endl;
//cout<<"n_pi:mean_pi:sig_pi"<<n_pi[20][0][0]<<" : "<<mean_pi[20][0][0]<<" : "<<sig_pi[20][0][0]<<endl;
//cout<<"n_k:mean_k:sig_k"<<n_k[20][0][0]<<" : "<<mean_k[20][0][0]<<" : "<<sig_k[20][0][0]<<endl;
//cout<<"n_p:mean_p:sig_p"<<n_p[20][0][0]<<" : "<<mean_p[20][0][0]<<" : "<<sig_p[20][0][0]<<endl;
//cout<<"n_L:mean_L:sig_L"<<n_L[20][0][0]<<" : "<<mean_L[20][0][0]<<" : "<<sig_L[20][0][0]<<endl;
//cout<<"n_S:mean_S:sig_S"<<n_S[20][0][0]<<" : "<<mean_S[20][0][0]<<" : "<<sig_S[20][0][0]<<endl;
//
//cout<<"Al front evNum (L) = "<<h_L_vz->Integral(h_L_vz->FindBin(-0.14),h_L_vz->FindBin(-0.1))<<endl;
//cout<<"Al front evNum (R) = "<<h_R_vz->Integral(h_R_vz->FindBin(-0.14),h_R_vz->FindBin(-0.1))<<endl;
//cout<<"Al back evNum (L) = "<<h_L_vz->Integral(h_L_vz->FindBin(0.1),h_L_vz->FindBin(0.14))<<endl;
//cout<<"Al back evNum (R) = "<<h_R_vz->Integral(h_R_vz->FindBin(0.1),h_R_vz->FindBin(0.14))<<endl;
//cout<<"H2 evNum (L) = "<<h_L_vz->Integral(h_L_vz->FindBin(-0.1),h_L_vz->FindBin(0.1))<<endl;
//cout<<"H2 evNum (R) = "<<h_R_vz->Integral(h_R_vz->FindBin(-0.1),h_R_vz->FindBin(0.1))<<endl;
//
//c29->Divide(2,2);
//c29->cd(1);
//fpi_z1 = new TF1("fpi_z1","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
//fpi_z2 = new TF1("fpi_z2","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
//fpi_z3 = new TF1("fpi_z3","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
//fpi_z4 = new TF1("fpi_z4","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
//fpi_z1->SetNpx(2000);
//fpi_z2->SetNpx(2000);
//fpi_z3->SetNpx(2000);
//fpi_z4->SetNpx(2000);
//fpi_z1->SetParameters(1000.,6,def_sig_pi,20.);
//fpi_z2->SetParameters(1000.,0,def_sig_pi,20.);
//fpi_z3->SetParameters(1000.,def_mean_pi,def_sig_pi,20.);
//fpi_z4->SetParameters(1000.,def_mean_pi,def_sig_pi,20.);
//h_z1->Fit("fpi_z1","","",min_coin_c,max_coin_c);
//h_z2->Fit("fpi_z2","","",min_coin_c,max_coin_c);
//h_z3->Fit("fpi_z3","","",min_coin_c,max_coin_c);
//h_z4->Fit("fpi_z4","","",min_coin_c,max_coin_c);
//double fitpar,fitpar2,fitpar3;
//fitpar=fpi_z1->GetParameter(0);//N
//fitpar2=fpi_z1->GetParameter(2);//sig
//fitpar3=fpi_z1->GetParameter(3);//const
//cout<<"S_1="<<fitpar/0.056<<", N_1="<<fitpar3*6*fitpar2/0.056<<endl;
//fitpar=fpi_z2->GetParameter(0);//N
//fitpar2=fpi_z2->GetParameter(2);//sig
//fitpar3=fpi_z2->GetParameter(3);//const
//cout<<"S_2="<<fitpar/0.056<<", N_2="<<fitpar3*6*fitpar2/0.056<<endl;
//fitpar=fpi_z3->GetParameter(0);//N
//fitpar2=fpi_z3->GetParameter(2);//sig
//fitpar3=fpi_z3->GetParameter(3);//const
//cout<<"S_3="<<fitpar/0.056<<", N_3="<<fitpar3*6*fitpar2/0.056<<endl;
//fitpar=fpi_z4->GetParameter(0);//N
//fitpar2=fpi_z4->GetParameter(2);//sig
//fitpar3=fpi_z4->GetParameter(3);//const
//cout<<"S_4="<<fitpar/0.056<<", N_4="<<fitpar3*6*fitpar2/0.056<<endl;
//c29->cd(1);
//h_z1->Draw("");
//fpi_z1->SetLineColor(kRed);
//fpi_z1->Draw("same");
//c29->cd(2);
//h_z2->Draw("");
//fpi_z2->SetLineColor(kRed);
//fpi_z2->Draw("same");
//c29->cd(3);
//h_z3->Draw("");
//fpi_z3->SetLineColor(kRed);
//fpi_z3->Draw("same");
//c29->cd(4);
//h_z4->Draw("");
//fpi_z4->SetLineColor(kRed);
//fpi_z4->Draw("same");

//c30->Divide(2,2);
//c30->cd(1);
//h_L_corL_tgth->Draw("colz");
//c30->cd(2);
//h_L_corR_tgth->Draw("colz");
//c30->cd(3);
//h_L_corL_tgph->Draw("colz");
//c30->cd(4);
//h_L_corR_tgph->Draw("colz");
//c31->Divide(2,2);
//c31->cd(1);
//h_L_corL_th->Draw("colz");
//c31->cd(2);
//h_L_corR_th->Draw("colz");
//c31->cd(3);
//h_L_corL_ph->Draw("colz");
//c31->cd(4);
//h_L_corR_ph->Draw("colz");
//c32->Divide(2,3);
//c32->cd(1);
//h_L_corL_x->Draw("colz");
//c32->cd(2);
//h_L_corR_x->Draw("colz");
//c32->cd(3);
//h_L_corL_y->Draw("colz");
//c32->cd(4);
//h_L_corR_y->Draw("colz");
//c32->cd(5);
//h_L_corL_z->Draw("colz");
//c32->cd(6);
//h_L_corR_z->Draw("colz");
//
//
//c34->Divide(2,2);
//c34->cd(1);
//h_L_corL_zdiff->Draw("colz");
//c34->cd(2);
//h_L_corR_zdiff->Draw("colz");
//c34->cd(3);
//h_L_tgth_x->Draw("colz");
//c34->cd(4);
//h_L_tgph_y->Draw("colz");
//
//
//c35->Divide(2,2);
//c35->cd(1);
//h_L_tgth_y->Draw("colz");
//c35->cd(2);
//h_L_tgph_x->Draw("colz");
//c35->cd(3);
//h_L_tgth_z->Draw("colz");
//c35->cd(4);
//h_L_tgph_z->Draw("colz");
//
//c36->Divide(2,2);
//c36->cd(1);
//h_L_ct_zdiff->Draw("colz");
//c36->cd(2);
//h_L_corsum_zdiff->Draw("colz");
//c36->cd(3);
//h_L_cordiff_zdiff->Draw("colz");
//c36->cd(4);
//h_L_corL_corR->Draw("colz");
//
//c37->Divide(2,2);
//c37->cd(1)->SetLogz(1);
//h_L_ct_zdiff->Draw("colz");
//c37->cd(2)->SetLogz(1);
//h_L_corsum_zdiff->Draw("colz");
//c37->cd(3)->SetLogz(1);
//h_L_cordiff_zdiff->Draw("colz");
//c37->cd(4)->SetLogz(1);
//h_L_corL_corR->Draw("colz");

c38->Divide(2,2);
c38->cd(1);
 def_sig_p=0.852; def_mean_p=-8.0;
 def_sig_pi=0.443; def_mean_pi=3.0;
 def_sig_k=0.644; def_mean_k=0.0;
 def_acc=27.7;
hct_test->Draw("");
 facc_kc=new TF1("facc_kc","[0]",min_coin_c,max_coin_c);
 facc_kc->SetNpx(2000);
 fk_kc=new TF1("fk_kc","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fk_kc->SetNpx(2000);
 fpi_pic=new TF1("fpi_pic","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fpi_pic->SetNpx(2000);
 fp_pc=new TF1("fp_pc","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fp_pc->SetNpx(2000);
//hct_test->Fit("facc_kc","","",min_coin_c,min_coin_c+8.);
double cosnt = facc_kc->GetParameter(0);
fk_kc->SetParameters(10,def_mean_k,def_sig_k);
fk_kc->SetParLimits(0,0.,1000.);
fp_pc->SetParLimits(0,0.,1000.);
fpi_pic->SetParLimits(0,0.,1000.);
fpi_pic->SetParameters(100,def_mean_pi,def_sig_pi);
fp_pc->SetParameters(100,def_mean_p,def_sig_p);
fk_kc->SetParameter(3,cosnt);
fp_pc->SetParameter(3,cosnt);
fpi_pic->SetParameter(3,cosnt);
//hct_test->Fit("fk_kc","","",-1.,1.);
//hct_test->Fit("fp_pc","","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
//hct_test->Fit("fpi_pic","","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
fk_kc->Draw("same");
fp_pc->Draw("same");
fpi_pic->Draw("same");
c38->cd(2);
 facc_kc2=new TF1("facc_kc2","[0]",min_coin_c,max_coin_c);
 facc_kc2->SetNpx(2000);
 fk_kc2=new TF1("fk_kc2","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fk_kc2->SetNpx(2000);
 fpi_pic2=new TF1("fpi_pic2","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fpi_pic2->SetNpx(2000);
 fp_pc2=new TF1("fp_pc2","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fp_pc2->SetNpx(2000);
hct_test2->SetLineColor(kRed);hct_test2->Draw("");
//hct_test2->Fit("facc_kc2","","",min_coin_c,min_coin_c+8.);
cosnt = facc_kc2->GetParameter(0);
fk_kc2->SetParameters(10,def_mean_k,def_sig_k);
fpi_pic2->SetParameters(100,def_mean_pi,def_sig_pi);
fp_pc2->SetParameters(100,def_mean_p,def_sig_p);
fk_kc2->SetParameter(3,cosnt);
fp_pc2->SetParameter(3,cosnt);
fpi_pic2->SetParameter(3,cosnt);
fk_kc2->SetParLimits(0,0.,1000.);
fp_pc2->SetParLimits(0,0.,1000.);
fpi_pic2->SetParLimits(0,0.,1000.);
//hct_test2->Fit("fk_kc2","","",-1.,1.);
//hct_test2->Fit("fp_pc2","","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
//hct_test2->Fit("fpi_pic2","","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
fk_kc2->Draw("same");
fp_pc2->Draw("same");
fpi_pic2->Draw("same");
c38->cd(3);
 facc_kc3=new TF1("facc_kc3","[0]",min_coin_c,max_coin_c);
 facc_kc3->SetNpx(2000);
 fk_kc3=new TF1("fk_kc3","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fk_kc3->SetNpx(2000);
 fpi_pic3=new TF1("fpi_pic3","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fpi_pic3->SetNpx(2000);
 fp_pc3=new TF1("fp_pc3","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fp_pc3->SetNpx(2000);
hct_test3->SetLineColor(kBlue);hct_test3->Draw("");
//hct_test3->Fit("facc_kc3","","",min_coin_c,min_coin_c+8.);
cosnt = facc_kc3->GetParameter(0);
fk_kc3->SetParameters(10,def_mean_k,def_sig_k);
fpi_pic3->SetParameters(100,def_mean_pi,def_sig_pi);
fp_pc3->SetParameters(100,def_mean_p,def_sig_p);
fk_kc3->SetParameter(3,cosnt);
fp_pc3->SetParameter(3,cosnt);
fpi_pic3->SetParameter(3,cosnt);
fk_kc3->SetParLimits(0,0.,1000.);
fp_pc3->SetParLimits(0,0.,1000.);
fpi_pic3->SetParLimits(0,0.,1000.);
//hct_test3->Fit("fk_kc3","","",-1.,1.);
//hct_test3->Fit("fp_pc3","","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
//hct_test3->Fit("fpi_pic3","","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
fk_kc3->Draw("same");
fp_pc3->Draw("same");
fpi_pic3->Draw("same");
cout<<"N(before)="<<hct_test->Integral(hct_test->FindBin(min_coin_c),hct_test->FindBin(max_coin_c))<<endl;
cout<<"N(after)="<<hct_test2->Integral(hct_test2->FindBin(min_coin_c),hct_test2->FindBin(max_coin_c))<<endl;;
c39->Divide(2,2);
c39->cd(1);
hct_test->Draw("");
c39->cd(2);
hct_test2->Draw("");
c39->cd(3);
hct_test3->Draw("");
c39->cd(4);
hct_test->SetLineColor(kBlack);
hct_test2->SetLineColor(kRed);
hct_test3->SetLineColor(kAzure);
hct_test2->Draw("");
hct_test3->Draw("same");
hct_test->Draw("same");

}//Draw()
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
 c20->Print(Form("%s",ofname.c_str()));
 c21->Print(Form("%s",ofname.c_str()));
 c22->Print(Form("%s",ofname.c_str()));
 c23->Print(Form("%s",ofname.c_str()));
 c24->Print(Form("%s",ofname.c_str()));
 c25->Print(Form("%s",ofname.c_str()));
 c26->Print(Form("%s",ofname.c_str()));
 c27->Print(Form("%s",ofname.c_str()));
 c28->Print(Form("%s",ofname.c_str()));
 c29->Print(Form("%s",ofname.c_str()));
 c30->Print(Form("%s",ofname.c_str()));
 c31->Print(Form("%s",ofname.c_str()));
 c32->Print(Form("%s",ofname.c_str()));
 c33->Print(Form("%s",ofname.c_str()));
 c34->Print(Form("%s",ofname.c_str()));
 c35->Print(Form("%s",ofname.c_str()));
 c36->Print(Form("%s",ofname.c_str()));
 c37->Print(Form("%s",ofname.c_str()));
 c38->Print(Form("%s",ofname.c_str()));
 c39->Print(Form("%s",ofname.c_str()));
 c39->Print(Form("%s]",ofname.c_str()));
    
 cout<<"Print is done "<<endl;
}//Print()


///////////////////////////////////////////////////////////////

void tuning::Write(){
cout<<"before tree_out"<<endl;
 tree_out->Write();
cout<<"before file_out"<<endl;
 file_out->Write();
cout<<"between file_out"<<endl;
 file_out->Close();
cout<<"after file_out"<<endl;
//fnew->Close();
}//Write()



// #################################################
//=======================================================//
//================     Main       =======================//
//=======================================================//


int main(int argc, char** argv){

	cout<<"nth = "<<nth<<endl;
  string ifname = "../small_1.list";//Run111157~111220
  //string ifname = "../small2.list";//Run111157~111220 & Run111480~111542
//    string ifname = "../test.list";//for debug
    string pname = "./Lambda_H1.param";//Run111157~111220
    //string pname = "./Lambda_H2.param";//Run111480~111542
    string mtparam = "../matrix/matrix_new.list";
    string print_name = "./test_print.pdf";
    bool print_flag = false;//.pdf
    bool root_flag  = true;//.root
    bool draw_flag  = false;//Draw canvas

  
    tuning* AC=new tuning();
    TApplication theApp("App", &argc, argv);
cout << "Start SetRunList" << endl;
    AC->matrix(mtparam);
cout << "Start SetBranch" << endl;
//    AC->SetBranch();
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
    AC->Fitting();
cout << "AC->Fitting() is done" << endl;
    AC->ACtune();
cout << "AC->ACtune() is done" << endl;
    if(draw_flag)AC->Draw();
    if(print_flag)AC->Print(print_name);
    if(root_flag)AC->Write();
  
cout<<"=========== Output files ============="<<endl;
cout<<"h2all_1_temp.root was created!"<<endl;
cout<<"If you want to confirm, please rename it to h2all_1_Nov.root"<<endl;
  
  
	
   if(draw_flag==0)gROOT->SetBatch(1);
   if(draw_flag==0)gSystem->Exit(1);
   if(root_flag)theApp.Run();
   return 0;
}//main
