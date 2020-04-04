using namespace std;
#include "Setting.cc"


void makePDF(){

  bool Lambda=true;
  bool nnL=false;
  
  
  gROOT->SetBatch(1);
  Setting* set=new Setting();
  set->Initialize();
  //  string name="Lambda_small_H_Ltuned_mm";
  string name="Lambda_small_H_0th";  
    //"Lambda_small1_mm_Ltuned";


  string root=".root";
  string pdf=".pdf";
  string ifname="../rootfiles/mmass/ana_Lambda/";
  ifname=ifname+name+root;
  TFile* fin=new TFile(ifname.c_str(),"read");
 
  //=============================//
  //======== Get Hist ===========//
  //============================//

  
  //====================//
  //==== MMass Hist ====//
  //====================//
  
  
  TH1D* h_mm=(TH1D*)fin->Get("h_mm");
  set->SetTH1(h_mm,"Lambda Mass w/o AC cut","Mass [GeV]","Counts/2 MeV");  
  TH1D* h_mm_acc=(TH1D*)fin->Get("h_mm_acc");
  TH1D* h_mm_L=(TH1D*)fin->Get("h_mm_L");
  set->SetTH1(h_mm,"Lambda Mass w/ AC cut","Mass [GeV]","Counts/2 MeV");    
  TH1D* h_acc_L=(TH1D*)fin->Get("h_acc_L");
  TH1D* h_mm_nnL=(TH1D*)fin->Get("h_mm_nnL");
  set->SetTH1(h_mm,"nnL Mass w/ AC cut","Mass [GeV]","Counts/2 MeV");    
  TH1D* h_acc_nnL=(TH1D*)fin->Get("h_acc_nnL");
  TH1D* h_peak_L=(TH1D*)fin->Get("h_peak_L");
  set->SetTH1(h_peak_L,"Lambda Mass Peak hist w/ AC cut","mass [GeV]","Counts/2 MeV");    
  TH1D* h_peak_nnL=(TH1D*)fin->Get("h_peak_nnL");
  set->SetTH1(h_peak_nnL,"nnL Mass Peak hist w/ AC cut","mass [GeV]","Counts/2 MeV");


  //====================//
  //==== LHRS Hist =====//
  //====================//

  //------- mass vs L-FP ------------//
  TH2D* h_Lx_mm=(TH2D*)fin->Get("h_Lx_mm");
    set->SetTH2(h_Lx_mm,"Mass vs L-FPx","mass [GeV]","FPx [m]");
  TH2D* h_Ly_mm=(TH2D*)fin->Get("h_Ly_mm");  
    set->SetTH2(h_Ly_mm,"Mass vs L-FPy","mass [GeV]","FPy [m]");    
  TH2D* h_Lth_mm=(TH2D*)fin->Get("h_Lth_mm");  
    set->SetTH2(h_Lth_mm,"Mass vs L-FPth","mass [GeV]","FPth [rad]");        
  TH2D* h_Lph_mm=(TH2D*)fin->Get("h_Lph_mm");
     set->SetTH2(h_Lph_mm,"Mass vs L-FPph","mass [GeV]","FPph [rad]");
     //--- mass vs L-Target ---------//
  TH2D* h_Lvz_mm=(TH2D*)fin->Get("h_Lvz_mm");
    set->SetTH2(h_Lvz_mm,"Mass vs L-Target z","mass [GeV]","z [m]");
  TH2D* h_Ltgth_mm=(TH2D*)fin->Get("h_Ltgth_mm");  
    set->SetTH2(h_Ltgth_mm,"Mass vs L-Target theta","mass [GeV]","theta [rad]");
  TH2D* h_Ltgph_mm=(TH2D*)fin->Get("h_Ltgph_mm");  
    set->SetTH2(h_Ltgph_mm,"Mass vs L-Target phi","mass [GeV]","phi [rad]");        
  TH2D* h_Lp_mm=(TH2D*)fin->Get("h_Lp_mm");
     set->SetTH2(h_Lp_mm,"Mass vs Momentum","mass [GeV]","p [GeV]");

  //====================//
  //==== RHRS Hist =====//
  //====================//

     
  //------- mass vs R-FP ------------//     
  TH2D* h_Rx_mm=(TH2D*)fin->Get("h_Rx_mm");
    set->SetTH2(h_Rx_mm,"Mass vs R-FPx","mass [GeV]","FPx [m]");
  TH2D* h_Ry_mm=(TH2D*)fin->Get("h_Ry_mm");  
    set->SetTH2(h_Ry_mm,"Mass vs R-FPy","mass [GeV]","FPy [m]");    
  TH2D* h_Rth_mm=(TH2D*)fin->Get("h_Rth_mm");  
    set->SetTH2(h_Rth_mm,"Mass vs R-FPth","mass [GeV]","FPth [rad]");        
  TH2D* h_Rph_mm=(TH2D*)fin->Get("h_Rph_mm");
     set->SetTH2(h_Rph_mm,"Mass vs R-FPph","mass [GeV]","FPph [rad]");


     //--- mass vs R-Target ---------//
  TH2D* h_Rvz_mm=(TH2D*)fin->Get("h_Rvz_mm");
    set->SetTH2(h_Rvz_mm,"Mass vs R-Target z","mass [GeV]","z [m]");
  TH2D* h_Rtgth_mm=(TH2D*)fin->Get("h_Rtgth_mm");  
    set->SetTH2(h_Rtgth_mm,"Mass vs R-Target theta","mass [GeV]","theta [rad]");
  TH2D* h_Rtgph_mm=(TH2D*)fin->Get("h_Rtgph_mm");  
    set->SetTH2(h_Rtgph_mm,"Mass vs R-Target phi","mass [GeV]","phi [rad]");        
  TH2D* h_Rp_mm=(TH2D*)fin->Get("h_Rp_mm");
     set->SetTH2(h_Rp_mm,"Mass vs Momentum","mass [GeV]","p [GeV]");



 
     TF1* fL=new TF1("fL","gausn(0)",1.13,1.16);
     fL->SetLineColor(2);
     fL->SetParNames("const","mean","sigma");
     h_peak_L->Fit("fL","Rq","",1.13,1.16);
     gStyle->SetOptStat(000000011);
     //     gStyle->SetStatW(0.2);
     //     gStyle->SetOptFit();
          gStyle->SetOptFit(0010);
     TF1* fnnL=new TF1("fnnL","gausn(0)",3.02,3.07);
     fnnL->SetLineColor(2);
     fnnL->SetParNames("const","mean","sigma");     
     h_peak_nnL->Fit("fnnL","Rq","",3.02,3.07);
          gStyle->SetOptFit(0010);

     

     //=============================//
     //======= Draw Hist ===========//
     //=============================//
     double min_mm,max_mm;
     if(Lambda){min_mm=1.05; max_mm=1.35;}
     if(nnL){min_mm=2.95; max_mm=1.25;}

  
  TCanvas *c30 = new TCanvas("c30","Mssing Mass 1",1000,800);
  c30->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(0);
  h_mm->SetLineColor(2);
  h_mm->SetFillStyle(3002);
  h_mm->SetFillColor(2);
  h_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
  h_mm      ->Draw();
  h_mm_acc->SetLineColor(4);
  h_mm_acc->SetFillStyle(3002);
  h_mm_acc->SetFillColor(4);  
  h_mm_acc ->Draw("same"); 

  
  TCanvas *c31 = new TCanvas("c31","Mssing Mass w/ AC cut",1000,800);
  c31->Divide(2,1);
  c31->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(0);
  h_mm_L->SetLineColor(2);
  h_mm_L->SetFillStyle(3002);
  h_mm_L->SetFillColor(2);
  h_mm_L->GetXaxis()->SetRangeUser(1.05,1.35);
  h_mm_L      ->Draw();
  h_acc_L->SetLineColor(4);
  h_acc_L->SetFillStyle(3002);
  h_acc_L->SetFillColor(4);
  h_acc_L ->Draw("same");

  c31->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(0);
  h_mm_nnL->SetLineColor(2);
  h_mm_nnL->SetFillStyle(3002);
  h_mm_nnL->SetFillColor(2);
  h_mm_nnL->GetXaxis()->SetRangeUser(2.95,3.25);  
  h_mm_nnL      ->Draw();
  h_acc_nnL->SetLineColor(4);
  h_acc_nnL->SetFillStyle(3002);
  h_acc_nnL->SetFillColor(4);  
  h_acc_nnL ->Draw("same");

  
   TCanvas *c32 = new TCanvas("c32","Mssing Mass Peak w/ AC cut",1000,800);
   c32->Divide(2,1);
   h_peak_L->GetXaxis()->SetRangeUser(1.05,1.35);
   h_peak_nnL->GetXaxis()->SetRangeUser(2.95,3.25);
   c32->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(0);// gStyle->SetOptFit(0111);
   h_peak_L      ->Draw();  fL->Draw("same");
   c32->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(0);// gStyle->SetOptFit(0111);
   h_peak_nnL    ->Draw();  fnnL->Draw("same");



   
   TCanvas *c33 = new TCanvas("c33","Mssing Mass vs L-FP",1000,800);
    c33->Divide(2,2);
    h_Lx_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Ly_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Lth_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Lph_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Lx_mm->GetYaxis()->SetRangeUser(-0.6,0.6);
    h_Ly_mm->GetYaxis()->SetRangeUser(-0.06,0.06);
    h_Lth_mm->GetYaxis()->SetRangeUser(-0.15,0.15);
    h_Lph_mm->GetYaxis()->SetRangeUser(-0.06,0.06);
    c33->cd(1) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Lx_mm   ->Draw("colz");    gPad->SetGridx(1);
    c33->cd(2) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Ly_mm   ->Draw("colz");    gPad->SetGridx(1);
    c33->cd(3) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Lth_mm  ->Draw("colz");    gPad->SetGridx(1);
    c33->cd(4) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Lph_mm  ->Draw("colz");    gPad->SetGridx(1);


  TCanvas *c34 = new TCanvas("c34","Mssing Mass vs R-FP",1000,800);
    c34->Divide(2,2);
    h_Rx_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Ry_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Rth_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Rph_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Rx_mm->GetYaxis()->SetRangeUser(-0.6,0.6);
    h_Ry_mm->GetYaxis()->SetRangeUser(-0.06,0.06);
    h_Rth_mm->GetYaxis()->SetRangeUser(-0.15,0.15);
    h_Rph_mm->GetYaxis()->SetRangeUser(-0.06,0.06);    
    c34->cd(1) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rx_mm   ->Draw("colz");    gPad->SetGridx(1);
    c34->cd(2) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Ry_mm   ->Draw("colz");    gPad->SetGridx(1);
    c34->cd(3) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rth_mm  ->Draw("colz");    gPad->SetGridx(1);
    c34->cd(4) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rph_mm  ->Draw("colz");    gPad->SetGridx(1);


   TCanvas *c35 = new TCanvas("c35","Mssing Mass vs L-Target",1000,800);
    c35->Divide(2,2);
    h_Lvz_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Ltgth_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Ltgph_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Lp_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Lvz_mm->GetYaxis()->SetRangeUser(-0.15,0.15);
    h_Ltgth_mm->GetYaxis()->SetRangeUser(-0.10,0.10);
    h_Ltgph_mm->GetYaxis()->SetRangeUser(-0.06,0.06);
    h_Lp_mm->GetYaxis()->SetRangeUser(1.95,2.25);    
    c35->cd(1) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Lvz_mm    ->Draw("colz");    gPad->SetGridx(1);
    c35->cd(2) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Ltgth_mm  ->Draw("colz");    gPad->SetGridx(1);
    c35->cd(3) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Ltgph_mm  ->Draw("colz");    gPad->SetGridx(1);
    c35->cd(4) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Lp_mm     ->Draw("colz");    gPad->SetGridx(1);


  TCanvas *c36 = new TCanvas("c36","Mssing Mass vs R-Target",1000,800);

    h_Rvz_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Rtgth_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Rtgph_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Rp_mm->GetXaxis()->SetRangeUser(min_mm,max_mm);
    h_Rvz_mm->GetYaxis()->SetRangeUser(-0.15,0.15);
    h_Rtgth_mm->GetYaxis()->SetRangeUser(-0.10,0.10);
    h_Rtgph_mm->GetYaxis()->SetRangeUser(-0.06,0.06);
    h_Rp_mm->GetYaxis()->SetRangeUser(1.7,1.95);
    c36->Divide(2,2);
    c36->cd(1) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rvz_mm    ->Draw("colz");    gPad->SetGridx(1);
    c36->cd(2) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rtgth_mm  ->Draw("colz");    gPad->SetGridx(1);
    c36->cd(3) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rtgph_mm  ->Draw("colz");    gPad->SetGridx(1);
    c36->cd(4) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rp_mm     ->Draw("colz");    gPad->SetGridx(1);



    string ofname="../pdf/mmass/";
    ofname=ofname+name+pdf;
    
    c30->Print(Form("%s[",ofname.c_str()));
    c30->Print(Form("%s",ofname.c_str()));
    c31->Print(Form("%s",ofname.c_str()));
    c32->Print(Form("%s",ofname.c_str()));
    c33->Print(Form("%s",ofname.c_str()));
    c34->Print(Form("%s",ofname.c_str()));
    c35->Print(Form("%s",ofname.c_str()));
    c36->Print(Form("%s",ofname.c_str()));     
    c30->Print(Form("%s]",ofname.c_str()));

    cout<<"input rootfile: "<<ifname<<endl;
    cout<<"output pdffile: "<<ofname<<endl;

        gSystem->Exit(1);
    
}

