const double PI = 4.*atan(1.);
string dir = "/usr/local/apps/strangecalc/share/strangecalc-wrapper/data/";

void SetTH2(TH2F *h2, TString hname, TString xname, TString yname, double offset, double Tsize){
  h2->SetTitle(hname);
  h2->GetXaxis()->SetTitle(xname);
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->SetTitle(yname);
  h2->GetYaxis()->CenterTitle();
  h2->GetYaxis()->SetTitleOffset(offset);
  h2->GetXaxis()->SetLabelOffset(offset);
  h2->GetXaxis()->SetLabelSize(Tsize);
  h2->GetXaxis()->SetTitleSize(Tsize);
  h2->GetYaxis()->SetLabelSize(Tsize);
  h2->GetYaxis()->SetTitleSize(Tsize);
}

void SetTH2_forDraw(TH2F *h2, TString hname, TString xname, TString yname){
  h2->SetTitle(hname);
  h2->GetXaxis()->SetTitle(xname);
  h2->GetXaxis()->CenterTitle();
  h2->GetXaxis()->SetTitleFont(42);
  h2->GetXaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetTitleOffset(1.10);
  h2->GetXaxis()->SetLabelSize(0.04);
  h2->GetXaxis()->SetTitleSize(0.05);

  h2->GetYaxis()->SetTitle(yname);
  h2->GetYaxis()->CenterTitle();
  h2->GetYaxis()->SetTitleFont(42);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetTitleOffset(1.00);
  h2->GetYaxis()->SetLabelSize(0.04);
  h2->GetYaxis()->SetTitleSize(0.05);
}

void SetGr(TGraph *gr, TString hname, TString xname, TString yname, int LColor, int MColor, int MStyle, double MSize){
  gr->SetTitle(hname);
  gr->SetName(hname);
  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->SetLineColor(LColor);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
//  gr->GetYaxis()->SetTitleOffset(Yoffset);
}

void SetGr(TGraphErrors *gr, TString hname, TString xname, TString yname, int LColor, int MColor, int MStyle, double MSize){
  gr->SetTitle(hname);
  gr->SetName(hname);
  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->SetLineColor(LColor);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
//  gr->GetYaxis()->SetTitleOffset(Yoffset);
}

void all(){
  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso2.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("SAPHIR_2004__Glander__EPJ_A19_251__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<"  "<<dataset->GetEntries()<<endl;
  ifp->GetObject("CLAS_2004__McNabb__PRC69_042201R__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<"  "<<dataset->GetEntries()<<endl;
  ifp->GetObject("CLAS_2006__Bradford__PRC73_035202__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<"  "<<dataset->GetEntries()<<endl;
  ifp->GetObject("CLAS_2010__Dey__PRC_82_025202__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<"  "<<dataset->GetEntries()<<endl;
  ifp->GetObject("LEPS_2006__Sumihama__PRC73_035214__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<"  "<<dataset->GetEntries()<<endl;
  ifp->GetObject("LEPS_2006__Kohri__PRL_082003__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<"  "<<dataset->GetEntries()<<endl;

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void safhir_2004(){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso2.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("SAPHIR_2004__Glander__EPJ_A19_251__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<endl;
//  dataset->ViewSelections();

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,10);

  TGraph *gr[33];
  TLatex *t1[33], *t2[33];
  for(int i=0;i<33;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.35,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.31,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08);
    t2[i]->SetTextSize(0.08);
#if 0
    double *x = gr[i]->GetX(), *y = gr[i]->GetY();
    cout<<Form("  // W= %.3lf GeV",dataset->GetKinematics()->GetW()/1000.)<<endl;
    for(int j=gr[i]->GetN()-1;j>=0;j--){ cout<<Form("%+.2lf, ",x[j])<<flush; }
    cout<<endl<<"  { "<<flush;
    for(int j=gr[i]->GetN()-1;j>=1;j--){ cout<<Form("%.5lf, ",y[j])<<flush; }
    cout<<Form("%.5lf },",y[0])<<endl;
#endif
//  dataset->Scan("observable:qsquared:w:costhkcm:amplitude:error");
  }

  TCanvas *c1 = new TCanvas("c1","c1",1800,950);
  c1->Divide(10,5);
  for(int i=0;i<33;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void clas_2004(){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso2.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("CLAS_2004__McNabb__PRC69_042201R__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<endl;
//  dataset->ViewSelections();

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,10);

  TGraph *gr[50];
  TLatex *t1[50], *t2[50];
  for(int i=0;i<50;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.35,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.31,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08);
    t2[i]->SetTextSize(0.08);
#if 0
    double *x = gr[i]->GetX(), *y = gr[i]->GetY();
    cout<<Form("  // W= %.3lf GeV",dataset->GetKinematics()->GetW()/1000.)<<endl;
    for(int j=gr[i]->GetN()-1;j>=0;j--){ cout<<Form("%+.2lf, ",x[j])<<flush; }
    cout<<endl<<"  { "<<flush;
    for(int j=gr[i]->GetN()-1;j>=1;j--){ cout<<Form("%.5lf, ",y[j])<<flush; }
    cout<<Form("%.5lf },",y[0])<<endl;
#endif
//  dataset->Scan("observable:qsquared:w:costhkcm:amplitude:error");
  }

  TCanvas *c1 = new TCanvas("c1","c1",1800,950);
  c1->Divide(10,5);
  for(int i=0;i<50;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void clas_2006(){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso2.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("CLAS_2006__Bradford__PRC73_035202__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<endl;
//  dataset->ViewSelections();

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,10);

  TGraph *gr[73];
  TLatex *t1[73], *t2[73];
  for(int i=0;i<73;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.35,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.31,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08);
    t2[i]->SetTextSize(0.08);
#if 0
    double *x = gr[i]->GetX(), *y = gr[i]->GetY();
    cout<<Form("  // W= %.3lf GeV",dataset->GetKinematics()->GetW()/1000.)<<endl;
    for(int j=gr[i]->GetN()-1;j>=0;j--){ cout<<Form("%+.2lf, ",x[j])<<flush; }
    cout<<endl<<"  { "<<flush;
    for(int j=gr[i]->GetN()-1;j>=1;j--){ cout<<Form("%.5lf, ",y[j])<<flush; }
    cout<<Form("%.5lf },",y[0])<<endl;
#endif
//  dataset->Scan("observable:qsquared:w:costhkcm:amplitude:error");
  }

  TCanvas *c1 = new TCanvas("c1","c1",1800,950);
  c1->Divide(10,5);
  for(int i=0;i<50;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
  TCanvas *c2 = new TCanvas("c2","c2",1800,950);
  c2->Divide(10,5);
  for(int i=0;i<23;i++){ c2->cd(i+1); h2->Draw(); gr[i+50]->Draw("Psame"); t1[i+50]->Draw("same"); t2[i+50]->Draw("same"); }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void clas_2010(){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso2.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("CLAS_2010__Dey__PRC_82_025202__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<endl;
//  dataset->ViewSelections();

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,10);

  TGraph *gr[112];
  TLatex *t1[112], *t2[112];
  for(int i=0;i<112;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.35,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.31,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08);
    t2[i]->SetTextSize(0.08);
#if 0
    double *x = gr[i]->GetX(), *y = gr[i]->GetY();
    cout<<Form("  // W= %.3lf GeV",dataset->GetKinematics()->GetW()/1000.)<<endl;
    for(int j=gr[i]->GetN()-1;j>=0;j--){ cout<<Form("%+.2lf, ",x[j])<<flush; }
    cout<<endl<<"  { "<<flush;
    for(int j=gr[i]->GetN()-1;j>=1;j--){ cout<<Form("%.5lf, ",y[j])<<flush; }
    cout<<Form("%.5lf },",y[0])<<endl;
#endif
//  dataset->Scan("observable:qsquared:w:costhkcm:amplitude:error");
  }

  TCanvas *c1 = new TCanvas("c1","c1",1800,950);
  c1->Divide(10,5);
  for(int i=0;i<50;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
  TCanvas *c2 = new TCanvas("c2","c2",1800,950);
  c2->Divide(10,5);
  for(int i=0;i<50;i++){ c2->cd(i+1); h2->Draw(); gr[i+50]->Draw("Psame"); t1[i+50]->Draw("same"); t2[i+50]->Draw("same"); }
  TCanvas *c3 = new TCanvas("c3","c3",1800,950);
  c3->Divide(10,5);
  for(int i=0;i<12;i++){ c3->cd(i+1); h2->Draw(); gr[i+100]->Draw("Psame"); t1[i+100]->Draw("same"); t2[i+100]->Draw("same"); }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void leps_2006_1(){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso2.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("LEPS_2006__Sumihama__PRC73_035214__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<endl;
//  dataset->ViewSelections();

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,10);

  TGraph *gr[18];
  TLatex *t1[18], *t2[18];
  for(int i=0;i<18;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.35,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.31,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08);
    t2[i]->SetTextSize(0.08);
#if 0
    double *x = gr[i]->GetX(), *y = gr[i]->GetY();
    cout<<Form("  // W= %.3lf GeV",dataset->GetKinematics()->GetW()/1000.)<<endl;
    for(int j=gr[i]->GetN()-1;j>=0;j--){ cout<<Form("%+.2lf, ",x[j])<<flush; }
    cout<<endl<<"  { "<<flush;
    for(int j=gr[i]->GetN()-1;j>=1;j--){ cout<<Form("%.5lf, ",y[j])<<flush; }
    cout<<Form("%.5lf },",y[0])<<endl;
#endif
//  dataset->Scan("observable:qsquared:w:costhkcm:amplitude:error");
  }

  TCanvas *c1 = new TCanvas("c1","c1",1800,950);
  c1->Divide(10,5);
  for(int i=0;i<18;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void leps_2006_2(){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso2.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("LEPS_2006__Kohri__PRL_082003__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<endl;
//  dataset->ViewSelections();

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,10);

  TGraph *gr[18];
  TLatex *t1[18], *t2[18];
  for(int i=0;i<18;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.35,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.31,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08);
    t2[i]->SetTextSize(0.08);
#if 0
    double *x = gr[i]->GetX(), *y = gr[i]->GetY();
    cout<<Form("  // W= %.3lf GeV",dataset->GetKinematics()->GetW()/1000.)<<endl;
    for(int j=gr[i]->GetN()-1;j>=0;j--){ cout<<Form("%+.2lf, ",x[j])<<flush; }
    cout<<endl<<"  { "<<flush;
    for(int j=gr[i]->GetN()-1;j>=1;j--){ cout<<Form("%.5lf, ",y[j])<<flush; }
    cout<<Form("%.5lf },",y[0])<<endl;
#endif
//  dataset->Scan("observable:qsquared:w:costhkcm:amplitude:error");
  }

  TCanvas *c1 = new TCanvas("c1","c1",1800,950);
  c1->Divide(10,5);
  for(int i=0;i<18;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void all_pdf(){
  gStyle->SetPadGridX(0);  gStyle->SetPadGridY(1);
  gStyle->SetGridWidth(0);  gStyle->SetFrameLineWidth(0);  gStyle->SetLineWidth(0);
  gStyle->SetOptStat("");
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first
  gROOT->SetBatch(1);

  TFile *ifp;
  TDataset *dataset;
  ifp = new TFile(Form("%sdata_iso2.root",dir.c_str()));

  string ofname = "Exp_KpSigmaZ.pdf";
  TCanvas *c1 = new TCanvas("c1","c1",1800,950);
  c1->Print(Form("%s[",ofname.c_str()));

  TPaveText *p1 = new TPaveText(0.2,0.5,0.8,0.7,"NDC");
  p1->SetTextSize(0.05);
  p1->SetFillColor(10);
  p1->SetBorderSize(1);
  TText *title;

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.5);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,10);
  TGraph *gr[200];
  TLatex *t1[200], *t2[200];

/////////
  ifp->GetObject("SAPHIR_2004__Glander__EPJ_A19_251__dcs",dataset);
  for(int i=0;i<33;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.45,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.40,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08); t1[i]->SetTextColor(4);
    t2[i]->SetTextSize(0.08); t2[i]->SetTextColor(4);
  }
  title = p1->AddText("SAPHIR 2004"); p1->Draw();
  c1->Print(Form("%s",ofname.c_str()));  p1->Clear();  c1->Clear();
  TCanvas *c1 = new TCanvas("c1","c1",1800,950);
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<13;i++){ c1->cd(i+1); h2->Draw(); gr[i+20]->Draw("Psame"); t1[i+20]->Draw("same"); t2[i+20]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();

/////////
  ifp->GetObject("CLAS_2004__McNabb__PRC69_042201R__dcs",dataset);
  for(int i=0;i<50;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.45,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.40,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08); t1[i]->SetTextColor(4);
    t2[i]->SetTextSize(0.08); t2[i]->SetTextColor(4);
  }
  title = p1->AddText("CLAS 2004"); p1->Draw();
  c1->Print(Form("%s",ofname.c_str()));  p1->Clear();  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame");    t1[i]->Draw("same");    t2[i]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i+20]->Draw("Psame"); t1[i+20]->Draw("same"); t2[i+20]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<10;i++){ c1->cd(i+1); h2->Draw(); gr[i+40]->Draw("Psame"); t1[i+40]->Draw("same"); t2[i+40]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();

/////////
  ifp->GetObject("CLAS_2006__Bradford__PRC73_035202__dcs",dataset);
  for(int i=0;i<73;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.45,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.40,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08); t1[i]->SetTextColor(4);
    t2[i]->SetTextSize(0.08); t2[i]->SetTextColor(4);
  }
  title = p1->AddText("CLAS 2006"); p1->Draw();
  c1->Print(Form("%s",ofname.c_str()));  p1->Clear();  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame");     t1[i]->Draw("same");    t2[i]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i+20]->Draw("Psame");  t1[i+20]->Draw("same"); t2[i+20]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i+40]->Draw("Psame");  t1[i+40]->Draw("same"); t2[i+40]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<13;i++){ c1->cd(i+1); h2->Draw(); gr[i+60]->Draw("Psame");  t1[i+60]->Draw("same"); t2[i+60]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();

/////////
  ifp->GetObject("CLAS_2010__Dey__PRC_82_025202__dcs",dataset);
  for(int i=0;i<112;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.45,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.40,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08); t1[i]->SetTextColor(4);
    t2[i]->SetTextSize(0.08); t2[i]->SetTextColor(4);
  }
  title = p1->AddText("CLAS 2010"); p1->Draw();
  c1->Print(Form("%s",ofname.c_str()));  p1->Clear();  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i+20]->Draw("Psame");  t1[i+20]->Draw("same"); t2[i+20]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i+40]->Draw("Psame");  t1[i+40]->Draw("same"); t2[i+40]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i+60]->Draw("Psame");  t1[i+60]->Draw("same"); t2[i+60]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i+80]->Draw("Psame");  t1[i+80]->Draw("same"); t2[i+80]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<12;i++){ c1->cd(i+1); h2->Draw(); gr[i+100]->Draw("Psame"); t1[i+100]->Draw("same");t2[i+100]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();

/////////
  ifp->GetObject("LEPS_2006__Kohri__PRL_082003__dcs",dataset);
  for(int i=0;i<18;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.45,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.40,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08); t1[i]->SetTextColor(4);
    t2[i]->SetTextSize(0.08); t2[i]->SetTextColor(4);
  }
  title = p1->AddText("LEPS 2006 1"); p1->Draw();
  c1->Print(Form("%s",ofname.c_str()));  p1->Clear();  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<18;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();

/////////
  ifp->GetObject("LEPS_2006__Sumihama__PRC73_035214__dcs",dataset);
  for(int i=0;i<18;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.45,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.40,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08); t1[i]->SetTextColor(4);
    t2[i]->SetTextSize(0.08); t2[i]->SetTextColor(4);
  }
  title = p1->AddText("LEPS 2006 2"); p1->Draw();
  c1->Print(Form("%s",ofname.c_str()));  p1->Clear();  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<18;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();

  c1->Print(Form("%s]",ofname.c_str()));
  gSystem->Exit(1);

}

////////////////////////////////////////////////////////////////////////////////////////////
void draw_KpSigmaZ(){
  gStyle->SetOptDate(0);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLineWidth(0);
  gStyle->SetOptDate(0);
  gStyle->SetTextFont(42);
  gStyle->SetGridWidth(0);
  gStyle->SetFrameLineWidth(0);

  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso2.root",dir.c_str()));
  TDataset *dataset1, *dataset2, *dataset3, *dataset4;
  ifp->GetObject("SAPHIR_2004__Glander__EPJ_A19_251__dcs",dataset1);
  ifp->GetObject("CLAS_2004__McNabb__PRC69_042201R__dcs" ,dataset2);
  ifp->GetObject("CLAS_2006__Bradford__PRC73_035202__dcs",dataset3);
  ifp->GetObject("CLAS_2010__Dey__PRC_82_025202__dcs"    ,dataset4);

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.5);
  h2->SetStats(kFALSE);
//  SetTH2_forDraw(h2,"","cos(#theta_{K^{+}}^{c.m.})","d#sigma/d#Omega_{K}^{c.m.} [#mub/sr]");
  SetTH2_forDraw(h2,"","","");
  h2->GetXaxis()->SetLabelOffset(1.);
  h2->GetXaxis()->SetTickLength(0);
  h2->GetYaxis()->SetLabelOffset(1.);
  h2->GetYaxis()->SetLabelSize(0.);
  h2->GetYaxis()->SetTickLength(0);

// SAPHIR 2004
  TGraphErrors *gr1;
  dataset1->SetSelection(8);
  gr1 = dataset1->MakeGraph("costhkcm");
  double *x1 = gr1->GetX(), *y1 = gr1->GetY(), ey1[100];
    for(int i=0;i<gr1->GetN();i++){ x1[i] *= -1; ey1[i] = gr1->GetErrorY(i); }
  gr1->Clear();
  gr1 = new TGraphErrors(gr1->GetN(),x1,y1,0,ey1);
  SetGr(gr1,"","","",4,4,24,1.5);

// CLAS 2004
  TGraphErrors *gr2;
  dataset2->SetSelection(11);
  gr2 = dataset2->MakeGraph("costhkcm");
  double *x2 = gr2->GetX(), *y2 = gr2->GetY(), ey2[100];
    for(int i=0;i<gr2->GetN();i++){ x2[i] *= -1; ey2[i] = gr2->GetErrorY(i); }
  gr2->Clear();
  gr2 = new TGraphErrors(gr2->GetN(),x2,y2,0,ey2);
  SetGr(gr2,"","","",1,1,22,1.5);

// CLAS 2006
  TGraphErrors *gr3;
  dataset3->SetSelection(10);
  gr3 = dataset3->MakeGraph("costhkcm");
  double *x3 = gr3->GetX(), *y3 = gr3->GetY(), ey3[100];
    for(int i=0;i<gr3->GetN();i++){ x3[i] *= -1; ey3[i] = gr3->GetErrorY(i); }
  gr3->Clear();
  gr3 = new TGraphErrors(gr3->GetN(),x3,y3,0,ey3);
  SetGr(gr3,"","","",3,3,21,1.5);

// CLAS 2010
  TGraphErrors *gr4;
  dataset4->SetSelection(15);
  gr4 = dataset4->MakeGraph("costhkcm");
  double *x4 = gr4->GetX(), *y4 = gr4->GetY(), ey4[100];
    for(int i=0;i<gr4->GetN();i++){ x4[i] *= -1; ey4[i] = gr4->GetErrorY(i); }
  gr4->Clear();
  gr4 = new TGraphErrors(gr4->GetN(),x4,y4,0,ey4);
  SetGr(gr4,"","","",2,2,20,1.5);

// KMaid W=1.835 GeV
  double x5[37] = {
      0.,  5., 10., 15., 20., 25., 30., 35., 40., 45.,
     50., 55., 60., 65., 70., 75., 80., 85., 90., 95.,
    100.,105.,110.,115.,120.,125.,130.,135.,140.,145.,
    150.,155.,160.,165.,170.,175.,180. };
  double y5[37] = {
    0.1565E+03,0.1576E+03,0.1606E+03,0.1654E+03,0.1716E+03,0.1787E+03,0.1862E+03,0.1934E+03,0.2000E+03,0.2053E+03,
    0.2091E+03,0.2111E+03,0.2110E+03,0.2089E+03,0.2048E+03,0.1988E+03,0.1910E+03,0.1817E+03,0.1712E+03,0.1596E+03,
    0.1472E+03,0.1344E+03,0.1214E+03,0.1084E+03,0.9581E+02,0.8370E+02,0.7233E+02,0.6184E+02,0.5237E+02,0.4401E+02,
    0.3680E+02,0.3077E+02,0.2591E+02,0.2219E+02,0.1957E+02,0.1801E+02,0.1749E+02 };
  for(int i=0;i<37;i++){ x5[i] = -cos(x5[i]/180.*PI); y5[i] /= 1000.; }
  TGraphErrors *gr5 = new TGraphErrors(37,x5,y5,0,0);
  SetGr(gr5,"","","",1,1,1,1.0);
  gr5->SetLineWidth(2);
  gr5->SetLineStyle(1);


  TF1 *fx = new TF1("fx","-x",-1,1);
  TGaxis *ax = new TGaxis(-1,0, 1,0, "fx", 505 );
  ax->SetLabelFont(42); ax->SetTitleFont(42); ax->CenterTitle();
  ax->SetLabelOffset(0.01);
  ax->SetTitleSize(0.06);  ax->SetTitleOffset(0.80);
  ax->SetMaxDigits(1);
  ax->SetTitle("cos(#theta_{K}^{c.m.})");

  TGaxis *ay = new TGaxis(-1,0, -1,0.5, 0,0.5, 505 );
  ay->SetLabelFont(42); ay->SetTitleFont(42); ay->CenterTitle();
  ay->SetLabelOffset(0.01);
  ay->SetTitleSize(0.06);  ay->SetTitleOffset(0.80);
  ay->SetTitle("d#sigma/d#Omega_{K}^{c.m.} (#mub/sr)");

  TLegend *leg = new TLegend(0.40,0.65,0.98,0.98,"","NDC");
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(3);
  leg->SetLineWidth(1);
  leg->SetTextSize(0.050);
  leg->SetMargin(0.15);
  leg->SetNColumns(1);
  leg->AddEntry(gr1     ,Form("SAPHIR #scale[0.6]{EPJA19(2004) 251.}  #scale[0.5]{(#sqrt{s}=%.3lf GeV)}"
                        ,dataset1->GetKinematics()->GetW()/1000.),"pl");
  leg->AddEntry(gr2     ,Form("CLAS #scale[0.6]{PRC69(2004) 042201R.} #scale[0.5]{(#sqrt{s}=%.3lf GeV)}"
                        ,dataset2->GetKinematics()->GetW()/1000.),"pl");
  leg->AddEntry(gr3     ,Form("CLAS #scale[0.6]{PRC73(2006) 035202.}  #scale[0.5]{(#sqrt{s}=%.3lf GeV)}"
                        ,dataset3->GetKinematics()->GetW()/1000.),"pl");
  leg->AddEntry(gr4     ,Form("CLAS #scale[0.6]{PRC82(2010) 025202.}  #scale[0.5]{(#sqrt{s}=%.3lf GeV)}"
                        ,dataset4->GetKinematics()->GetW()/1000.),"pl");
  leg->AddEntry(gr5     ,"KMaid #scale[0.6]{PRC61(1999) 012201.}  #scale[0.5]{(#sqrt{s}=1.835 GeV)}","pl");

  TCanvas *c1 = new TCanvas("c1","c1",1100,800);
  c1->Divide(1,1,0.0001,0.0001);
  c1->cd(1);
  h2->Draw();
  gr1->Draw("Psame");
  gr2->Draw("Psame");
  gr3->Draw("Psame");
  gr4->Draw("Psame");
  //gr5->Draw("Lsame");
  ax->Draw();
  ay->Draw();
  leg->Draw();
}
