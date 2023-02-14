#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <TChain.h>
#include <TMinuit.h>
#include <iostream>
#include <fstream>
//#include "Param.h"

const int nfoil = 10;
double fcent[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
		       0.00, 0.025, 0.05, 0.10, 0.125};
double fcent_real[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
			    0.000, 0.025, 0.050, 0.100, 0.125};
double selection_width = 0.01; // event selection width for z
const double l0 = 1.03; // [m]
//const double l0 = 103; // [cm]
double dth[nfoil];
const double step = 0.492 * 2.54;
const int nrow = 11; // the number of row in SSX pattern
const int ncol = 8;  // the number of column in SSY pattern
const int nsshole = nrow*ncol; // the number of holes to consider
const double hrs_ang = 13.2/180.*3.14;
double l[nfoil];
//double refx[nsshole];
//double refy[nsshole];
double ssx_cent_real[nsshole];
double ssy_cent_real[nsshole];
double conv_inch_to_cm =2.54;

bool rasy =true;


void ang_acc_L(){


  //  string main = "ang_LHRS_4th_0914_0";
  //  string main = "ang_LHRS_sieve_init";
  //    string main = "ang_LHRS_5th_0601";
  //  string main = "test";
  //  string main = "test_LHRS_wRas_test";
  
  string main = "test_LHRS_woRas_test";
  if(rasy)  main = "ang_LHRS_wRasy"; // raster y 
  string end_root = ".root";
  string end_dat = ".dat";
  string ifname="../../rootfiles/angcalib/" + main + ".root";
  //  string ofname="./param/" + main + ".dat";

  TFile* f = new TFile(ifname.c_str() );
  TTree* t = (TTree*)f->Get("T");

  double ss_x,ss_y,xp[100],yp[100], vz[100],gs;
  int zfoil;
  t -> SetBranchAddress("L.tr.tg_th_tuned",xp);
  t -> SetBranchAddress("L.tr.tg_ph_tuned",yp);
  t -> SetBranchAddress("ss_x",&ss_x);
  t -> SetBranchAddress("ss_y",&ss_y);
  t -> SetBranchAddress("L.cer.asum_c",&gs);
  t -> SetBranchAddress("L.tr.vz_tuned",vz);
  t -> SetBranchAddress("zfoil",&zfoil);

  TMarker* mark[nsshole];
  TMarker* mark_real[nfoil][nsshole];  
  double w[nfoil][nsshole];
  //  double ssx_off[nfoil][nsshole];
  //  double ssy_off[nfoil][nsshole];
  //  bool flag[nfoil][nsshole];
  double refx[nsshole],refy[nsshole],refx_real[nfoil][nsshole],refy_real[nfoil][nsshole];
  double ssy_cent_real[nsshole],ssx_cent_real[nsshole];
  double ref_xp[nsshole][nfoil],ref_yp[nsshole][nfoil];
  

  TMarker* Mss[nsshole][nfoil];
  TMarker* Mang[nsshole][nfoil];

  TMarker* Mss_real[nsshole][nfoil];
  TMarker* Mang_real[nsshole][nfoil];  
  
    for(int i=0;i<nfoil;i++){
      l[i] = 0;
      l[i]=(l0-fcent_real[i]/cos(hrs_ang))*100.;
      dth[i] = atan(l0*sin(hrs_ang)/(l0*cos(hrs_ang) -fcent_real[i]));
    }
    
		  
    bool RHRS =false;
    TLine* xp_ref_line[nrow][nfoil];
    TLine* yp_ref_line[ncol][nfoil];
    TLine* yp_ref_line2[ncol][nfoil];
    TLine* ssx_line[nrow][nfoil];
    TLine* ssy_line[ncol][nfoil];
    int nhole=0;
    for(int k=0;k<nfoil;k++){
      nhole=0;
      for(int i=0; i<ncol ; i++){
	for(int j=0; j<nrow; j++){
	  ssy_cent_real[i] = -3.0*step + step*i;
	  if(j%2==0)ssy_cent_real[i] = ssy_cent_real[i] - step/2.0;
	  ssx_cent_real[j] = 5.0*step - step*j;
	  refx[nhole] = ssx_cent_real[j];
	  refy[nhole] = ssy_cent_real[i];
	  
	  if(!RHRS ) ref_yp[nhole][k] = - cos(dth[k])/(l[k]/refy[nhole] - sin(dth[k]));
	  else       ref_yp[nhole][k] = - cos(dth[k])/(l[k]/refy[nhole] + sin(dth[k]));
	  
	  double lx;
	  if(refy[nhole]>0)lx=sqrt(pow(l[k],2.0) + pow(refy[nhole],2.0) + 2.0*l[k]*refy[nhole]*sin(dth[k]));
	  else     lx=sqrt(pow(l[k],2.0) + pow(refy[nhole],2.0) - 2.0*l[k]*refy[nhole]*sin(dth[k]));
	  ref_xp[nhole][k] = - refx[nhole]/lx;

	  Mss[nhole][k] = new TMarker(refy[nhole],refx[nhole],28);
	  Mss[nhole][k] -> SetMarkerColor(2);
	  Mang[nhole][k] = new TMarker(ref_yp[nhole][k],ref_xp[nhole][k],28);
	  Mang[nhole][k] -> SetMarkerColor(2);	  
	  if(i==0) xp_ref_line[j][k] = new TLine(ref_xp[nhole][k],0,ref_xp[nhole][k],300.);
	  if(j==0) yp_ref_line[i][k] = new TLine(ref_yp[nhole][k],0,ref_yp[nhole][k],500.);
	  else if(j==1)
	    yp_ref_line2[i][k] = new TLine(ref_yp[nhole][k],0,ref_yp[nhole][k],500.);
	  
	  nhole++;
	}
      }
    }


   //======== Loop ============//
   
   TH1D* hssx[nsshole][nfoil];
   TH1D* hssy[nsshole][nfoil];
   TH1D* hxp[nsshole][nfoil];
   TH1D* hyp[nsshole][nfoil];
   TH1D* hssx_c[nsshole][nfoil];
   TH1D* hssy_c[nsshole][nfoil];
   TH1D* hxp_c[nsshole][nfoil];
   TH1D* hyp_c[nsshole][nfoil];
   
   
   int    bin_xp =  1000;
   double min_xp = -0.1;
   double max_xp =  0.1;
   int    bin_yp =  1000;
   double min_yp = -0.1;
   double max_yp =  0.1;
   int    bin_ssx =  1000;
   double min_ssx = -10.;
   double max_ssx =  10.;
   int    bin_ssy =  1000;
   double min_ssy = -10.;
   double max_ssy =  10.;      

   TH1D* hssx_sum_c[nfoil];
   TH1D* hssy_sum_c[nfoil];
   TH1D* hxp_sum_c[nfoil];
   TH1D* hyp_sum_c[nfoil];

   TH1D* hxp_diff[nfoil];
   TH1D* hyp_diff[nfoil];
   TH1D* hssx_diff[nfoil];
   TH1D* hssy_diff[nfoil];
   
   TH1D* hxp_[nfoil];
   TH1D* hyp_[nfoil];
   TH1D* hssx_[nfoil];
   TH1D* hssy_[nfoil];
   
   //const int nrow = 11; // the number of row in SSX pattern
   //const int ncol = 8;  // the number of column in SSY pattern

   
   TH2D* hss= new TH2D("hss","SS with cut; ssy [cm]  ; ssx [cm]",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);   
   TH2D* hss_cut= new TH2D("hss_cut","SS with Event Select ; ssy [cm]  ; ssx [cm]",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);

   TH2D* hang = new TH2D("hang","Xp vs Yp; Yp [rad]; Xp [rad] ",bin_yp,min_yp,max_yp,bin_xp,min_xp,max_xp);
   TH2D* hang_cut = new TH2D("hang","Xp vs Yp with cut; Yp [rad]; Xp [rad] ",bin_yp,min_yp,max_yp,bin_xp,min_xp,max_xp);   



   
   int ihole =0;
   
   for(int j=0 ; j<nfoil ; j++){
     ihole =0;

     hssx_sum_c[j] = new TH1D(Form("hssx_sum_c_%d",j),"SSX Diff sum; Diff [cm]; Counts",bin_ssx,min_ssx,max_ssx);
     hssy_sum_c[j] = new TH1D(Form("hssy_sum_c_%d",j),"SSY Diff sum; Diff [cm]; Counts",bin_ssy,min_ssy,max_ssy);
     hxp_sum_c[j]  = new TH1D(Form("hxp_sum_c_%d",j),"Xp Diff sum; Diff [rad]; Counts",bin_xp,min_xp,max_xp);
     hyp_sum_c[j]  = new TH1D(Form("hyp_sum_c_%d",j),"Yp Diff sum ; Diff [rad]; Counts",bin_yp,min_yp,max_yp);
     hxp_[j] = new TH1D(Form("hxp_%d",j),Form("Xp Hist foil %d",j),bin_xp,min_xp,max_xp);
     hyp_[j] = new TH1D(Form("hyp_%d",j),Form("Yp Hist foil %d",j),bin_yp,min_yp,max_yp);

     hssx_[j] = new TH1D(Form("hssx_%d",j),Form("SSX Hist foil %d; ssx [cm] ;Counts",j),bin_ssx,min_ssx,max_ssx);
     hssy_[j] = new TH1D(Form("hssy_%d",j),Form("SSY Hist foil %d; ssy [cm] ;Counts",j),bin_ssy,min_ssy,max_ssy);
   
     for(int col=0 ; col<ncol ; col++){
       for(int row=0 ; row<nrow ; row++){
	 hssx[ihole][j] = new TH1D(Form("hssx_%d_%d",ihole,j),"", bin_ssx,min_ssx,max_ssx);
	 hssy[ihole][j] = new TH1D(Form("hssy_%d_%d",ihole,j),"", bin_ssy,min_ssy,max_ssy);
	 hxp[ihole][j] = new TH1D(Form("hxp_%d_%d",ihole,j),"", bin_xp,min_xp,max_xp);
	 hyp[ihole][j] = new TH1D(Form("hyp_%d_%d",ihole,j),"", bin_yp,min_yp,max_yp);
	 hssx_c[ihole][j] = new TH1D(Form("hssx_c_%d_%d",ihole,j),"", bin_ssx,min_ssx,max_ssx);
	 hssy_c[ihole][j] = new TH1D(Form("hssy_c_%d_%d",ihole,j),"", bin_ssy,min_ssy,max_ssy);
	 hxp_c[ihole][j] = new TH1D(Form("hxp_c_%d_%d",ihole,j),"", bin_xp,min_xp,max_xp);
	 hyp_c[ihole][j] = new TH1D(Form("hyp_c_%d_%d",ihole,j),"", bin_yp,min_yp,max_yp);	 
	 ihole++;
       }
     }
   }
   
  
   int ENum = t ->GetEntries();
   cout<<"ENum : "<<ENum<<endl;
   double selec_widthx = 0.75; // cm
   double selec_widthy = 0.75; // cm
   for(int k=0;k<ENum;k++){

     t->GetEntry(k);


     // cut condition //

     //     if(gs>2000.){
     if(gs>2000.){
       

     hss->Fill(ss_y,ss_x);
     hang->Fill(yp[0],xp[0]);
     for(int j=0 ; j<nfoil ; j++){
       if(fcent[j]-selection_width<vz[0] && vz[0]<fcent[j]+selection_width){
	 ihole =0;
	 for(int col=0 ; col<ncol ; col++){
	   for(int row=0 ; row<nrow ; row++){


	     // ==== Select Hole Position =======//
	     if(pow(ss_x-(refx[ihole]),2.0)/pow(selec_widthx,2.0)
		+ pow(ss_y-(refy[ihole]),2.0)/pow(selec_widthy,2.0)<0.64){


	       if(zfoil==5)  hang_cut->Fill(yp[0],xp[0]);
	       if(zfoil==5)  hss_cut->Fill(ss_y,ss_x);
	       
	       hxp_[j] ->Fill(xp[0]);
	       hyp_[j] ->Fill(yp[0]);
	       hssx_[j] ->Fill(ss_x);
	       hssy_[j] ->Fill(ss_y);
	       hxp[ihole][j]   ->Fill(xp[0]);
	       hyp[ihole][j]   ->Fill(yp[0]);
	       hssx[ihole][j]  ->Fill(ss_x);
	       hssy[ihole][j]  ->Fill(ss_y);

	       hxp_c[ihole][j] ->Fill(xp[0] - ref_xp[ihole][j]);
	       hyp_c[ihole][j] ->Fill(yp[0] - ref_yp[ihole][j]);
	       hssx_c[ihole][j]->Fill(ss_x -  refx[ihole]     );
	       hssy_c[ihole][j]->Fill(ss_y -  refy[ihole]     );	       

	       hssx_sum_c[j] -> Fill((ss_x -  refx[ihole])     ); // cm
	       hssy_sum_c[j] -> Fill((ss_y -  refy[ihole])     ); // cm
	       hxp_sum_c[j]  -> Fill((xp[0] - ref_xp[ihole][j])     );  // rad
	       hyp_sum_c[j]  -> Fill((yp[0] - ref_yp[ihole][j])     );  // rad
	       
	       
	     } // END if
	     
	     ihole++;
	   }
	 }
	    
	 
       } // end if
     } // END j
     }// END cut
   } // END loop
   

   //==== Fitting zfoil ==============//
   TGraphErrors* gXp_mean =new TGraphErrors();
   gXp_mean ->SetName("gXp_mean");
   gXp_mean ->SetTitle("Xp Sum Diff ; #foil ; diff [mrad]");
   gXp_mean ->GetYaxis()->SetRangeUser(-1,1);
   gXp_mean ->SetMarkerStyle(7);
   gXp_mean ->SetMarkerColor(4);
   TGraphErrors* gYp_mean =new TGraphErrors();
   gYp_mean ->SetName("gYp_mean");
   gYp_mean ->SetTitle("Yp Sum Diff ; #foil ; diff [mrad]");
   gYp_mean ->GetYaxis()->SetRangeUser(-1,1);
   gYp_mean ->SetMarkerStyle(7);
   gYp_mean ->SetMarkerColor(2);   
   TGraphErrors* gSSX_mean =new TGraphErrors();
   gSSX_mean ->SetName("gSSX_mean");
   gSSX_mean ->SetTitle("SSX Sum Diff ; #foil ; diff [mm]");
   gSSX_mean ->GetYaxis()->SetRangeUser(-10,10);
   gSSX_mean ->SetMarkerStyle(7);
   gSSX_mean ->SetMarkerColor(4);
   
   TGraphErrors* gSSY_mean =new TGraphErrors();
   gSSY_mean ->SetName("gSSY_mean");
   gSSY_mean ->SetTitle("SSY Sum Diff ; #foil ; diff [mm]");
   gSSY_mean ->GetYaxis()->SetRangeUser(-10,10);
   gSSY_mean ->SetMarkerStyle(7);
   gSSY_mean ->SetMarkerColor(2);

   
   TGraphErrors* gXp_sig =new TGraphErrors();
   gXp_sig ->SetName("gXp_sig");
   gXp_sig ->SetTitle("Xp Sum sigma ; #foil ; resolution [mrad]");
   gXp_sig ->GetYaxis()->SetRangeUser(1,3);
   gXp_sig ->SetMarkerStyle(7);
   gXp_sig ->SetMarkerColor(4);
   TGraphErrors* gYp_sig =new TGraphErrors();
   gYp_sig ->SetName("gYp_sig");
   gYp_sig ->SetTitle("Yp Sum sigma ; #foil ; resolution [mrad]");
   gYp_sig ->GetYaxis()->SetRangeUser(1,3);
   gYp_sig ->SetMarkerStyle(7);
   gYp_sig ->SetMarkerColor(2);   
   TGraphErrors* gSSX_sig =new TGraphErrors();
   gSSX_sig ->SetName("gSSX_sig");
   gSSX_sig ->SetTitle("SSX sigma ; #foil ; resolution [mm]");
   gSSX_sig ->GetYaxis()->SetRangeUser(-10,10);
   gSSX_sig ->SetMarkerStyle(7);
   gSSX_sig ->SetMarkerColor(4);
   TGraphErrors* gSSY_sig =new TGraphErrors();
   gSSY_sig ->SetName("gSSY_sig");      
   gSSY_sig ->SetTitle("SSY sigma ; #foil ; resolution [mm]");
   gSSY_sig ->GetYaxis()->SetRangeUser(-10,10);
   gSSY_sig ->SetMarkerStyle(7);
   gSSY_sig ->SetMarkerColor(4);

   
   TF1* fxp = new TF1("fxp","gaus(0)",-1.0,1.0);
   TF1* fyp = new TF1("fyp","gaus(0)",-1.0,1.0);
   TF1* fssx = new TF1("fssx","gaus(0)",-1,1);
   TF1* fssy = new TF1("fssy","gaus(0)",-1,1);
   
   double xp_sum_mean,yp_sum_mean, xp_sum_sig, yp_sum_sig;
   double xp_sum_mean_err,yp_sum_mean_err, xp_sum_sig_err, yp_sum_sig_err;
   double ssx_sum_mean,ssy_sum_mean, ssx_sum_sig, ssy_sum_sig;
   double ssx_sum_mean_err,ssy_sum_mean_err, ssx_sum_sig_err, ssy_sum_sig_err;   

   for(int z=0;z<nfoil;z++){
     hxp_sum_c[z] ->Fit("fxp","QR","QR",-1.0,1.0);
     hyp_sum_c[z] ->Fit("fyp","QR","QR",-1.0,1.0);
     hssx_sum_c[z] ->Fit("fssx","QR","QR",-1,1);
     hssy_sum_c[z] ->Fit("fssy","QR","QR",-1,1);
     
     xp_sum_mean = fxp->GetParameter(1)*1000.;
     yp_sum_mean = fyp->GetParameter(1)*1000.;
     xp_sum_sig  = fxp->GetParameter(2)*1000.;
     yp_sum_sig  = fyp->GetParameter(2)*1000.;

     xp_sum_mean_err = fxp->GetParError(1)*1000.;
     yp_sum_mean_err = fyp->GetParError(1)*1000.;
     xp_sum_sig_err  = fxp->GetParError(2)*1000.;
     yp_sum_sig_err  = fyp->GetParError(2)*1000.;


     ssx_sum_mean = fssx->GetParameter(1)*10.;
     ssy_sum_mean = fssy->GetParameter(1)*10.;
     ssx_sum_sig  = fssx->GetParameter(2)*10.;
     ssy_sum_sig  = fssy->GetParameter(2)*10.;

     ssx_sum_mean_err = fssx->GetParError(1)*10.;
     ssy_sum_mean_err = fssy->GetParError(1)*10.;
     ssx_sum_sig_err  = fssx->GetParError(2)*10.;
     ssy_sum_sig_err  = fssy->GetParError(2)*10.;


     gXp_mean ->SetPoint(z,z,xp_sum_mean);
     gXp_sig  ->SetPoint(z,z,xp_sum_sig );
     gXp_mean ->SetPointError(z,0.0,xp_sum_mean_err);
     gXp_sig  ->SetPointError(z,0.0,xp_sum_sig_err );     

     gYp_mean ->SetPoint(z,z,yp_sum_mean);
     gYp_sig  ->SetPoint(z,z,yp_sum_sig );
     gYp_mean ->SetPointError(z,0.0,yp_sum_mean_err);
     gYp_sig  ->SetPointError(z,0.0,yp_sum_sig_err );     


     gSSX_mean ->SetPoint(z,z,ssx_sum_mean);
     gSSX_sig  ->SetPoint(z,z,ssx_sum_sig );
     gSSX_mean ->SetPointError(z,0.0,ssx_sum_mean_err);
     gSSX_sig  ->SetPointError(z,0.0,ssx_sum_sig_err );     

     gSSY_mean ->SetPoint(z,z,ssy_sum_mean);
     gSSY_sig  ->SetPoint(z,z,ssy_sum_sig );
     gSSY_mean ->SetPointError(z,0.0,ssy_sum_mean_err);
     gSSY_sig  ->SetPointError(z,0.0,ssy_sum_sig_err );          

     
   }
   

   //==== Fitting Hole Position =====//

   TGraphErrors* gXp_diff[nfoil];
   TGraphErrors* gYp_diff[nfoil];
   TGraphErrors* gSSX_diff[nfoil];
   TGraphErrors* gSSY_diff[nfoil];
   TGraphErrors* gTheta_diff[nfoil];
   
   double ssx_mean[nsshole][nfoil];
   double ssy_mean[nsshole][nfoil];
   double ssx_mean_err[nsshole][nfoil];
   double ssy_mean_err[nsshole][nfoil];
   double ssx_peak[nsshole][nfoil];
   double ssy_peak[nsshole][nfoil];

   double xp_mean[nsshole][nfoil];
   double yp_mean[nsshole][nfoil];
   double xp_mean_err[nsshole][nfoil];
   double yp_mean_err[nsshole][nfoil];
   double xp_peak[nsshole][nfoil];
   double yp_peak[nsshole][nfoil];
   
   TF1* fXp[nsshole][nfoil];
   TF1* fYp[nsshole][nfoil];
   TF1* fSSX[nsshole][nfoil];
   TF1* fSSY[nsshole][nfoil];

   double Theta_diff[nsshole][nfoil];
   TH2D* h_xp_diff[nfoil];
   TH2D* h_yp_diff[nfoil];
   
   TH1D* hxp_diff_ = new TH1D("hxp_diff","Xp Sum Diff; diff [mrad] ; Hole Counts",100,-10,10);
   TH1D* hyp_diff_ = new TH1D("hyp_diff","Yp Sum Diff; diff [mrad] ; Hole Counts",100,-10,10);


   
   
   for(int z=0;z<nfoil;z++){
     gXp_diff[z] =  new TGraphErrors();
     gXp_diff[z]->SetName(Form("gXp_diff_%d",z));
     gXp_diff[z]->SetMarkerStyle(7);
     gXp_diff[z]->SetMarkerColor(z+1);
     gYp_diff[z] =  new TGraphErrors();
     gYp_diff[z]->SetName(Form("gYp_diff_%d",z));
     gYp_diff[z]->SetMarkerStyle(7);
     gYp_diff[z]->SetMarkerColor(z+1);

     gSSX_diff[z] =  new TGraphErrors();
     gSSX_diff[z]->SetName(Form("gSSX_diff_%d",z));
     gSSX_diff[z]->SetMarkerStyle(7);
     gSSX_diff[z]->SetMarkerColor(z+1);
     gSSY_diff[z] =  new TGraphErrors();
     gSSY_diff[z]->SetName(Form("gSSY_diff_%d",z));
     gSSY_diff[z]->SetMarkerStyle(7);
     gSSY_diff[z]->SetMarkerColor(z+1);

     gSSX_diff[z]->SetTitle("SSX Mean Diff; # hole ; diff [mm]");
     gSSY_diff[z]->SetTitle("SSY Mean Diff; # hole ; diff [mm]");
     gXp_diff[z]->SetTitle("Xp Mean Diff; # hole ; diff [mrad]");
     gYp_diff[z]->SetTitle("Yp Mean Diff; # hole ; diff [mrad]");

     gTheta_diff[z] = new TGraphErrors();
     gTheta_diff[z] -> SetName(Form("gTheta_diff_%d",z));
     gTheta_diff[z]->SetMarkerStyle(7);
     gTheta_diff[z]->SetMarkerColor(z+1);

     
     hxp_diff[z] = new TH1D(Form("hxp_diff_%d",z),Form("Xp Sum Diff (Foil %d); diff [mrad] ; Hole Counts",z),100,-10,10);
     hyp_diff[z] = new TH1D(Form("hyp_diff_%d",z),Form("Yp Sum Diff (Foil %d); diff [mrad] ; Hole Counts",z),100,-10,10);

     hssx_diff[z] = new TH1D(Form("hssx_diff_%d",z),Form("SSX Sum Diff (Foil %d); diff [mm]; Hole Counts",z),100,-10,10);
     hssy_diff[z] = new TH1D(Form("hssy_diff_%d",z),Form("SSY Sum Diff (Foil %d); diff [mm]; Hole Counts",z),100,-10,10);

     
     h_xp_diff[z]= new TH2D(Form("h_xp_diff_%d",z),"",2*ncol,-0.05,0.05,2*nrow,-0.08,0.08);
     h_yp_diff[z]= new TH2D(Form("h_yp_diff_%d",z),"",2*ncol,-0.05,0.05,2*nrow,-0.08,0.08);
     
     if(z==9){
     gSSX_diff[z]->SetMarkerColor(z+2);
     gSSY_diff[z]->SetMarkerColor(z+2);
     gXp_diff[z]->SetMarkerColor(z+2);
     gYp_diff[z]->SetMarkerColor(z+2);          
     }
     
     for(int ihole =0;ihole<nsshole;ihole++){
       
       fXp[ihole][z]  = new TF1(Form("fXp_%d_%d",ihole,z),"gaus(0)",-3,3);
       fYp[ihole][z]  = new TF1(Form("fYp_%d_%d",ihole,z),"gaus(0)",-3,3);
       fSSX[ihole][z] = new TF1(Form("fSSX_%d_%d",ihole,z),"gaus(0)",-3,3);
       fSSY[ihole][z] = new TF1(Form("fSSY_%d_%d",ihole,z),"gaus(0)",-3,3);

    
       ssx_mean[ihole][z] = hssx_c[ihole][z]->GetBinCenter(hssx_c[ihole][z]->GetMaximumBin());
       ssx_peak[ihole][z] = hssx_c[ihole][z]->GetBinContent(hssx_c[ihole][z]->GetMaximumBin());
       ssy_mean[ihole][z] = hssy_c[ihole][z]->GetBinCenter(hssy_c[ihole][z]->GetMaximumBin());
       ssy_peak[ihole][z] = hssy_c[ihole][z]->GetBinContent(hssy_c[ihole][z]->GetMaximumBin());       

       fSSX[ihole][z]->SetParameter(0,ssx_peak[ihole][z]);
       fSSX[ihole][z]->SetParameter(1,ssx_mean[ihole][z]);
       fSSY[ihole][z]->SetParameter(0,ssy_peak[ihole][z]);
       fSSY[ihole][z]->SetParameter(1,ssy_mean[ihole][z]);       


       if(ssx_peak[ihole][z]>0){
	 hssx_c[ihole][z] ->Fit(Form("fSSX_%d_%d",ihole,z),"QR","QR",-1+ssx_mean[ihole][z], ssx_mean[ihole][z] +1);
	 ssx_mean[ihole][z]     = fSSX[ihole][z]->GetParameter(1);
	 ssx_mean_err[ihole][z] = fSSX[ihole][z]->GetParError(1);
       }
       
       if(ssy_peak[ihole][z]>0){
	 hssy_c[ihole][z] ->Fit(Form("fSSY_%d_%d",ihole,z),"QR","QR",-1+ssy_mean[ihole][z], ssy_mean[ihole][z] +1);
       ssy_mean[ihole][z]     = fSSY[ihole][z]->GetParameter(1);
       ssy_mean_err[ihole][z] = fSSY[ihole][z]->GetParError(1);
       }

       if(ssx_peak[ihole][z]>0){
	 gSSX_diff[z]->SetPoint(ihole,ihole,ssx_mean[ihole][z]*10.); // mm
	 gSSX_diff[z]->SetPointError(ihole,0.0,ssx_mean_err[ihole][z]*10.); //mm
       }else gSSX_diff[z]->SetPoint(ihole,ihole,-10.);
       if(ssy_peak[ihole][z]>0){
	 gSSY_diff[z]->SetPoint(ihole,ihole,ssy_mean[ihole][z]*10.); // mm
	 gSSY_diff[z]->SetPointError(ihole,0.0,ssy_mean_err[ihole][z]*10.); // mm
       }else gSSY_diff[z]->SetPoint(ihole,ihole,-10.);


       Mss_real[ihole][z] = new TMarker(ssy_mean[ihole][z]+ refy[ihole],ssx_mean[ihole][z] + refx[ihole],28);
       Mss_real[ihole][z]->SetMarkerColor(1);

       hssx_diff[z] -> Fill(ssx_mean[ihole][z]);
       hssy_diff[z] -> Fill(ssy_mean[ihole][z]);

//       hxp_diff[z]->Fill(xp_mean[ihole][z]*1000.);
//       hyp_diff[z]->Fill(yp_mean[ihole][z]*1000.);       

       
       //refx[ihole]
       
       //======== Xp Yp Fitting =======//

       xp_mean[ihole][z] = hxp_c[ihole][z]->GetBinCenter(hxp_c[ihole][z]->GetMaximumBin());
       xp_peak[ihole][z] = hxp_c[ihole][z]->GetBinContent(hxp_c[ihole][z]->GetMaximumBin());
       yp_mean[ihole][z] = hyp_c[ihole][z]->GetBinCenter(hyp_c[ihole][z]->GetMaximumBin());
       yp_peak[ihole][z] = hyp_c[ihole][z]->GetBinContent(hyp_c[ihole][z]->GetMaximumBin());       

       fXp[ihole][z]->SetParameter(0,xp_peak[ihole][z]);
       fXp[ihole][z]->SetParameter(1,xp_mean[ihole][z]);
       fYp[ihole][z]->SetParameter(0,yp_peak[ihole][z]);
       fYp[ihole][z]->SetParameter(1,yp_mean[ihole][z]);       


       
       if(xp_peak[ihole][z]>0){
	 hxp_c[ihole][z] ->Fit(Form("fXp_%d_%d",ihole,z),"QR","QR",-1+xp_mean[ihole][z], xp_mean[ihole][z] +1);
	 xp_mean[ihole][z] = fXp[ihole][z]->GetParameter(1);
	 xp_mean_err[ihole][z] = fXp[ihole][z]->GetParError(1);
       }
       
       if(yp_peak[ihole][z]>0){
	 hyp_c[ihole][z] ->Fit(Form("fYp_%d_%d",ihole,z),"QR","QR",-1+yp_mean[ihole][z], yp_mean[ihole][z] +1);
       yp_mean[ihole][z] = fYp[ihole][z]->GetParameter(1);
       yp_mean_err[ihole][z] = fYp[ihole][z]->GetParError(1);
       }

       
       if(xp_peak[ihole][z]>0){
	 gXp_diff[z]->SetPoint(ihole,ihole,xp_mean[ihole][z]*1000.); // mrad
	 //	 gXp_diff[z]->SetPointError(ihole,0.0,xp_mean_err[ihole][z]);
       }else gXp_diff[z]->SetPoint(ihole,ihole,-10.);
       if(yp_peak[ihole][z]>0){
	 gYp_diff[z]->SetPoint(ihole,ihole,yp_mean[ihole][z]*1000.); // mrad
	 //	 gYp_diff[z]->SetPointError(ihole,0.0,yp_mean_err[ihole][z]);
       }else gYp_diff[z]->SetPoint(ihole,ihole,-10.);


       //==== dTheta Calculation ======//

       double Xp  = xp_mean[ihole][z] + ref_xp[ihole][z];
       double dXp = xp_mean[ihole][z];
       double Yp  = yp_mean[ihole][z] + ref_yp[ihole][z];
       double dYp = yp_mean[ihole][z];
       
       double hrs_ang = 13.2*3.14/180.;
       double theta     =  acos(( -Yp*sin(hrs_ang)+ cos(hrs_ang) )/sqrt(1.+Xp*Xp + Yp*Yp));
       double cos_theta = ( -Yp*sin(hrs_ang)+ cos(hrs_ang) )/sqrt(1.+Xp*Xp + Yp*Yp);
       double d_theta = -1./sqrt(1.- cos_theta*cos_theta)*
	 ( cos_theta*Xp/(1.+Xp*Xp + Yp*Yp) * dXp // dX'
	   + (-sin(hrs_ang)/sqrt(1. + Xp*Xp + Yp*Yp)  +  cos_theta/(1. + Xp*Xp + Yp*Yp)*Yp )*dYp );    // dY'
       

       if(xp_peak[ihole][z]<=0 || yp_peak[ihole][z]<=0 )
	 d_theta =-10.;

       gTheta_diff[z]->SetPoint(ihole,ihole,d_theta*1000.); // mrad
       Theta_diff[ihole][z] = d_theta * 1000.;


	 
       hxp_diff[z]->Fill(xp_mean[ihole][z]*1000.);
       hyp_diff[z]->Fill(yp_mean[ihole][z]*1000.);       
       hxp_diff_->Fill(xp_mean[ihole][z]*1000.);
       hyp_diff_->Fill(yp_mean[ihole][z]*1000.);
       
       Mang_real[ihole][z] = new TMarker(yp_mean[ihole][z]+ ref_yp[ihole][z],xp_mean[ihole][z] + ref_xp[ihole][z],28);
       Mang_real[ihole][z]->SetMarkerColor(1);

       if(fabs(xp_mean[ihole][z]*1000.)<5.)
	 h_xp_diff[z] -> Fill(ref_xp[ihole][z], ref_yp[ihole][z],xp_mean[ihole][z]*1000.);
       if(fabs(yp_mean[ihole][z]*1000.)<5.)
       h_yp_diff[z] -> Fill(ref_xp[ihole][z], ref_yp[ihole][z],yp_mean[ihole][z]*1000.);
              
     }
   }



   TH1D* hTheta_diff[nfoil];

   for(int z=0;z<nfoil;z++){
     hTheta_diff[z] = new TH1D(Form("hTheta_diff_%d",z),"Theta Diff ; Diff [mrad] ; Hole Counts ",200,-5.,5.);
     for(int ihole =0;ihole<nsshole;ihole++)
       hTheta_diff[z] ->Fill(Theta_diff[ihole][z]);
   }
   
   

   TCanvas* c0 =new TCanvas("c0","c0");
   c0->cd();
   hang_cut->Draw("colz");
   for(int ihole =0;ihole<nsshole;ihole++){
     Mang[ihole][5]->Draw("same");
     Mang_real[ihole][5]->Draw("same");
     }


   TCanvas* c1 =new TCanvas("c1","c1");
   c1->cd();
   hss_cut->Draw("colz");
   for(int ihole =0;ihole<nsshole;ihole++){
     Mss[ihole][5]->Draw("same");
     Mss_real[ihole][5]->Draw("same");
     
     }


   TCanvas* c2 = new TCanvas("c2","c2");
   c2->Divide(1,2);
   c2->cd(1);
   hxp_[5] ->Draw();
   for(int irow=0;irow<nrow;irow++){
     xp_ref_line[irow][5]->SetLineColor(2);
     xp_ref_line[irow][5]->Draw("same");
   }

   c2->cd(2);
   hyp_[5] ->Draw();
   for(int icol=0;icol<ncol;icol++){
     yp_ref_line[icol][5]->SetLineColor(2);
     yp_ref_line[icol][5]->Draw("same");
     yp_ref_line2[icol][5]->SetLineColor(2);
     yp_ref_line2[icol][5]->Draw("same");     
   }

   
   string ofrname = "./hist/ang_acc_L.root";
   if(rasy) ofrname = "./hist/ang_acc_L_wRasy.root";
   TFile * fr =new TFile(ofrname.c_str(),"recreate");


   for(int z=0;z<nfoil;z++){     
     for(int ihole =0;ihole<nsshole;ihole++){
       hxp[ihole][z]->Write();
       hyp[ihole][z]->Write();
       hxp_c[ihole][z]->Write();
       hyp_c[ihole][z]->Write();
       hssx[ihole][z]->Write();
       hssx_c[ihole][z]->Write();
       hssy[ihole][z]->Write();
       hssy_c[ihole][z]->Write();
       
     }
     gXp_diff[z]->Write();
     gYp_diff[z]->Write();
     gSSX_diff[z]->Write();
     gSSY_diff[z]->Write();
     hssx_sum_c[z]->Write();
     hssy_sum_c[z]->Write();
     hxp_sum_c[z]->Write();
     hyp_sum_c[z]->Write();
     hxp_diff[z]->Write();
     hyp_diff[z]->Write();
     hxp_[z] ->Write();
     hyp_[z] ->Write();
     hssx_[z] ->Write();
     hssy_[z] ->Write();
     h_xp_diff[z]->Write();
     h_yp_diff[z]->Write();
     gTheta_diff[z]->Write();
     hTheta_diff[z] ->Write();
   }

   gXp_mean->Write();
   gYp_mean->Write();
   gXp_sig->Write();
   gYp_sig->Write();
   gSSX_mean->Write();
   gSSY_mean->Write();
   gSSX_sig->Write();
   gSSY_sig->Write();   
   hss->Write();
   hss_cut->Write();
   hang->Write();
   hang_cut->Write();
   hxp_diff_->Write();
   hyp_diff_->Write();
}
  

