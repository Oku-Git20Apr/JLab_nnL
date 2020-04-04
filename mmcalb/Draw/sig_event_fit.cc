

void sig_event_fit(){

  TFile* f1=new TFile(Form("../../rootfiles/mmass/ana_Lambda/Lambda_small_H_0th.root"));
  TFile* f2=new TFile(Form("/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/mmass/ana_Lambda/Lambda_small_H_0910.root"));
  //TFile* f1=new TFile(Form("../root/Lambda_small_H.root"));
  //  TFile* f2=new TFile(Form("/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/mmass/ana_Lambda/Lambda_small_H_0th.root"));
  //  TFile* f2=new TFile(Form("/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/mmass/ana_Lambda/Lambda_small_H_0910.root"));
  //  TH1D* hmm=(TH1D*)f1->Get("h_mm");
  //  TH1D* hmm_acc=(TH1D*)f1->Get("h_mm_acc");
  TH1D* hmm=(TH1D*)f1->Get("h_mm_L");
  TH1D* hmm_0=(TH1D*)f2->Get("h_mm_L");
  TH1D* hmm_acc=(TH1D*)f1->Get("h_acc_L");
  TH1D* hmm_p=new TH1D("hmm_p","",500,0.5,1.5);
   hmm_p->Add(hmm,hmm_acc,1.,-2./40.);

   TCanvas* c0=new TCanvas("c0","c0");
   c0->cd();
   hmm_p->Draw();   

   double pbg_L[2],pL[3];
   
   TF1* fbg_L=new TF1("fbg_L","[0]*x+[1]",1.1,1.15);
   fbg_L->SetParameters(-1496.17 , 1361.23);
   hmm_p->Fit("fbg_L","","",1.1,1.14);
   pbg_L[0]=fbg_L->GetParameter(0);
   pbg_L[1]=fbg_L->GetParameter(1);
   TF1* fL=new TF1("fL","[0]*x+[1] + gausn(2)",1.1,1.3);
   fL->SetParameters(pbg_L[0] , pbg_L[1],  1.10100e+00 ,  1.11304e+00, 2.00810e-03);
   hmm_p->Fit("fL","","",1.1,1.13);
   pL[0]=fL->GetParameter(2);
   pL[1]=fL->GetParameter(3);
   pL[2]=fL->GetParameter(4);

   
   double pbg_S[2],pS[3];
   TF1* fbg_S=new TF1("fbg_S","[0]*x+[1]",1.2,1.3);
   fbg_S->SetParameters(-7.35954e+02 , 9.69466e+02);
   hmm_p->Fit("fbg_S","","",1.15,1.23);
   pbg_S[0]=fbg_S->GetParameter(0);
   pbg_S[1]=fbg_S->GetParameter(1);
   TF1* fS=new TF1("fS","[0]*x+[1] + gausn(2)",1.1,1.3);
   fS->SetParameters(pbg_S[0] , pbg_S[1], 4.74407e-01, 1.18962e+00 , 4.26818e-03);
   //   fS->FixParameter(0,pbg_S[0]);
   //   fS->FixParameter(1,pbg_S[1]);
   hmm_p->Fit("fS","","",1.18,1.20);
   pS[0]=fS->GetParameter(2);
   pS[1]=fS->GetParameter(3);
   pS[2]=fS->GetParameter(4);
   
   TF1* fLp=new TF1("fLp","gausn(0)",1.0,1.3);
   TF1* fSp=new TF1("fSp","gausn(0)",1.0,1.3);
   
   fLp->SetParameters(pL[0],pL[1],pL[2]);
   fSp->SetParameters(pS[0],pS[1],pS[2]);

   
   TCanvas* c1=new TCanvas("c1","c1");
   c1->cd();

   fLp->SetLineColor(2);
   //   fLp->SetFillColor(2);
   //   fLp->SetFillStyle(3002);
   fSp->SetLineColor(kBlue);
   //   fSp->SetFillColor(kBlue);
   //   fSp->SetFillStyle(3002);
   fSp->SetNpx(2000);
   fLp->SetNpx(2000);
   hmm_acc->SetLineColor(kGreen);
   hmm_acc->SetFillColor(kGreen);
   hmm_acc->SetFillStyle(3002);

   hmm_0->Draw();
   hmm_acc->Draw("same");
   hmm->Draw("same");
   fLp->Draw("same");
   fSp->Draw("same");

}
