

void mm_fit(){




  //  TFile* f1=new TFile(Form("/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/mmass/ana_Lambda/Lambda_small_H_0916.root"));
    TFile* f1=new TFile(Form("/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/mmass/ana_Lambda/2020-01-06/Lambda_small_OleH.root"));
  //  TFile* f1=new TFile(Form("/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/mmass/ana_Lambda/Lambda_small_H_0912.root"));
  //TFile* f1=new TFile(Form("../root/Lambda_small_H.root"));
  //  TFile* f2=new TFile(Form("/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/mmass/ana_Lambda/Lambda_small_H_0th.root"));
  //  TFile* f2=new TFile(Form("/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/mmass/ana_Lambda/Lambda_small_H_0910.root"));
  //  TH1D* hmm=(TH1D*)f1->Get("h_mm");
  //  TH1D* hmm_acc=(TH1D*)f1->Get("h_mm_acc");
  TH1D* hmm=(TH1D*)f1->Get("h_mm_L");
  TH1D* hmm_0=(TH1D*)f1->Get("h_mm_L");
  TH1D* hmm_acc=(TH1D*)f1->Get("h_acc_L");
  //  TH1D* hmm_p=new TH1D("hmm_p","",400,-0.2,0.2);
  TH1D* hmm_p=(TH1D*)f1->Get("h_peak_L");
  //   hmm_p->Add(hmm,hmm_acc,1.,-1.);

   TCanvas* c0=new TCanvas("c0","c0");
   c0->cd();
   hmm_p->Draw();   

   double pbg_L[2],pL[3];

   double min_L,max_L,min_S,max_S;
   min_L=-0.01;
   max_L=0.01;
   min_S=0.067;
   max_S=0.0825;
   TF1* fbg_L=new TF1("fbg_L","[0]*x+[1]",-0.2,0.2);
   fbg_L->SetParameters(-1496.17 , 1361.23);
   hmm_p->Fit("fbg_L","","",min_L,max_L);
   pbg_L[0]=fbg_L->GetParameter(0);
   pbg_L[1]=fbg_L->GetParameter(1);
   TF1* fL=new TF1("fL","[0]*x+[1] + gausn(2)",-0.2,0.2);
   fL->SetParameters(pbg_L[0] , pbg_L[1],  1.22944e+00,3.67807e-04,3.13158e-03);
   hmm_p->Fit("fL","","",min_L,max_L);
   pL[0]=fL->GetParameter(2);
   pL[1]=fL->GetParameter(3);
   pL[2]=fL->GetParameter(4);


   
   double pbg_S[2],pS[3];
   TF1* fbg_S=new TF1("fbg_S","[0]*x+[1]",-0.2,0.2);
   fbg_S->SetParameters(19.5663,127.912 );
   hmm_p->Fit("fbg_S","","",0.05,0.1);
   pbg_S[0]=fbg_S->GetParameter(0);
   pbg_S[1]=fbg_S->GetParameter(1);
   TF1* fS=new TF1("fS","[0]*x+[1] + gausn(2)",-0.2,0.2);
   fS->SetParameters(pbg_S[0] , pbg_S[1],6.10877e-01,7.60154e-02,6.36921e-03);
   //   fS->SetParameters(pbg_S[0] , pbg_S[1], 4.06005e+01,7.63064e-02,5.73203e-03);
   //   fS->SetParameters(0.0 , 0.0, 4.06005e+01,7.63064e-02,5.73203e-03);
   fS->FixParameter(0,pbg_S[0]);
   fS->FixParameter(1,pbg_S[1]);
   hmm_p->Fit("fS","","",min_S,max_S);
   pS[0]=fS->GetParameter(2);
   pS[1]=fS->GetParameter(3);
   pS[2]=fS->GetParameter(4);
   
   TF1* fLp=new TF1("fLp","gausn(0)",-0.2,0.2);
   TF1* fSp=new TF1("fSp","gausn(0)",-0.2,0.2);
   
   fLp->SetParameters(pL[0],pL[1],pL[2]);
   fSp->SetParameters(pS[0],pS[1],pS[2]);

   
   TCanvas* c1=new TCanvas("c1","c1");
   c1->cd();

   fLp->SetLineColor(2);
   fLp->SetFillColor(2);
   fLp->SetFillStyle(3002);
   fSp->SetLineColor(kBlue);
   fSp->SetFillColor(kBlue);
   fSp->SetFillStyle(3002);
   fSp->SetNpx(2000);
   fLp->SetNpx(2000);
   hmm_acc->SetLineColor(kGreen);
   hmm_acc->SetFillColor(kGreen);
   hmm_acc->SetFillStyle(3002);

   //hmm_0->Draw();
   hmm->Draw();
   hmm_acc->Draw("same");
   hmm_acc->SetLineColor(8);
   hmm_acc->SetFillColor(8);
   hmm_acc->SetFillStyle(3002);
   fLp->Draw("same");
   fSp->Draw("same");

   
   cout<<"======================================"<<endl;
   cout<<"======== Fitting Result =============="<<endl;
   cout<<"======================================"<<endl;

   cout<<" -----< Lambda > -------"<<endl;
   cout<<"Events : "<<pL[0]/0.002<< " Counts "<<endl;
   cout<<"Mean : "<<pL[1]*1.0e3<< " MeV "<<endl;
   cout<<"Sigma : "<<pL[2]*1.0e3<< " MeV "<<endl;
   cout<<" -----< Sigma  > -------"<<endl;
   cout<<"Events : "<<pS[0]/0.002<< " Counts "<<endl;
   cout<<"Mean : "<<pS[1]*1.0e3<< " MeV "<<endl;
   cout<<"Sigma : "<<pS[2]*1.0e3<< " MeV "<<endl;
   
}
