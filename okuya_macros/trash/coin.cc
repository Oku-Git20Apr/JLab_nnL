void tuning::MTParam_G(){


  cout<<"================================"<<endl;
  cout<<"=======Gogami Param ============"<<endl;
  cout<<"================================"<<endl;

//  int nParamT_3=3;
  char name_MctimeL[100];
  char name_MctimeR[100];
//  sprintf(name_MctimeL,"../goga_mac/Rootfiles/matrices/ctimeL.dat"); 
//  sprintf(name_MctimeR,"../goga_mac/Rootfiles/matrices/ctimeR.dat"); 
  sprintf(name_MctimeL,"../ctimeL.dat"); 
  sprintf(name_MctimeR,"../ctimeR.dat"); 
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



 double ctime = - meantime_R + ctimecorL + ctimecorR;

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

 if(3600.0<ctime && ctime<3665){
   ctime = ctime - 3637.88 - 12.76;
   ctime = ctime - 12.0-3.1;
 }

 ctime=-ctime;
 return ctime;


}


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
