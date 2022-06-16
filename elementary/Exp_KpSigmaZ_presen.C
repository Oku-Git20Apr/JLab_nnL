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

  TH2F *h2 = new TH2F("h2","h2",1000,0.,180.,1000,0,0.5);
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
  //dataset1->SetSelection(8);
  dataset1->SetSelection(20);
  gr1 = dataset1->MakeGraph("costhkcm");
  double *x1 = gr1->GetX(), *y1 = gr1->GetY(), ey1[100];
    for(int i=0;i<gr1->GetN();i++){ x1[i] = acos(x1[i])*180./PI; ey1[i] = gr1->GetErrorY(i); }
  gr1->Clear();
  gr1 = new TGraphErrors(gr1->GetN(),x1,y1,0,ey1);
  SetGr(gr1,"","","",4,4,24,1.5);

// CLAS 2004
  TGraphErrors *gr2;
  //dataset2->SetSelection(11);
  dataset2->SetSelection(35);
  gr2 = dataset2->MakeGraph("costhkcm");
  double *x2 = gr2->GetX(), *y2 = gr2->GetY(), ey2[100];
    for(int i=0;i<gr2->GetN();i++){ x2[i] = acos(x2[i])*180./PI; ey2[i] = gr2->GetErrorY(i); }
  gr2->Clear();
  gr2 = new TGraphErrors(gr2->GetN(),x2,y2,0,ey2);
  SetGr(gr2,"","","",30,30,22,1.5);//Dark Green?

// CLAS 2006
  TGraphErrors *gr3;
  //dataset3->SetSelection(10);
  dataset3->SetSelection(35);
  gr3 = dataset3->MakeGraph("costhkcm");
  double *x3 = gr3->GetX(), *y3 = gr3->GetY(), ey3[100];
    for(int i=0;i<gr3->GetN();i++){ x3[i] = acos(x3[i])*180./PI; ey3[i] = gr3->GetErrorY(i); }
  gr3->Clear();
  gr3 = new TGraphErrors(gr3->GetN(),x3,y3,0,ey3);
  SetGr(gr3,"","","",17,17,21,1.5);//Grey

// CLAS 2010
  TGraphErrors *gr4;
  //dataset4->SetSelection(15);
  dataset4->SetSelection(43);
  gr4 = dataset4->MakeGraph("costhkcm");
  double *x4 = gr4->GetX(), *y4 = gr4->GetY(), ey4[100];
    for(int i=0;i<gr4->GetN();i++){ x4[i] = acos(x4[i])*180./PI; ey4[i] = gr4->GetErrorY(i); }
  gr4->Clear();
  gr4 = new TGraphErrors(gr4->GetN(),x4,y4,0,ey4);
  SetGr(gr4,"","","",1,1,20,1.5);//Black

// KMaid W=2.110 GeV
  double x5[37] = {
  //    0.,  5., 10., 15., 20., 25., 30., 35., 40., 45.,
  //   50., 55., 60., 65., 70., 75., 80., 85., 90., 95.,
  //  100.,105.,110.,115.,120.,125.,130.,135.,140.,145.,
  //  150.,155.,160.,165.,170.,175.,180. };// KMaid W=1.835 GeV
  0.0000,
  5.0000,
 10.0000,
 15.0000,
 20.0000,
 25.0000,
 30.0000,
 35.0000,
 40.0000,
 45.0000,
 50.0000,
 55.0000,
 60.0000,
 65.0000,
 70.0000,
 75.0000,
 80.0000,
 85.0000,
 90.0000,
 95.0000,
100.0000,
105.0000,
110.0000,
115.0000,
120.0000,
125.0000,
130.0000,
135.0000,
140.0000,
145.0000,
150.0000,
155.0000,
160.0000,
165.0000,
170.0000,
175.0000,
180.0000};

  double y5[37] = {
  //  0.1565E+03,0.1576E+03,0.1606E+03,0.1654E+03,0.1716E+03,0.1787E+03,0.1862E+03,0.1934E+03,0.2000E+03,0.2053E+03,
  //  0.2091E+03,0.2111E+03,0.2110E+03,0.2089E+03,0.2048E+03,0.1988E+03,0.1910E+03,0.1817E+03,0.1712E+03,0.1596E+03,
  //  0.1472E+03,0.1344E+03,0.1214E+03,0.1084E+03,0.9581E+02,0.8370E+02,0.7233E+02,0.6184E+02,0.5237E+02,0.4401E+02,
  //  0.3680E+02,0.3077E+02,0.2591E+02,0.2219E+02,0.1957E+02,0.1801E+02,0.1749E+02 };
0.1224E+02,
0.1694E+02,
0.3038E+02,
0.5070E+02,
0.7521E+02,
0.1009E+03,
0.1248E+03,
0.1448E+03,
0.1593E+03,
0.1678E+03,
0.1705E+03,
0.1681E+03,
0.1615E+03,
0.1519E+03,
0.1403E+03,
0.1276E+03,
0.1144E+03,
0.1014E+03,
0.8881E+02,
0.7701E+02,
0.6617E+02,
0.5641E+02,
0.4786E+02,
0.4060E+02,
0.3474E+02,
0.3038E+02,
0.2759E+02,
0.2640E+02,
0.2679E+02,
0.2864E+02,
0.3171E+02,
0.3562E+02,
0.3992E+02,
0.4406E+02,
0.4750E+02,
0.4977E+02,
0.5057E+02};

  for(int i=0;i<37;i++){  y5[i] /= 1000.; }
  TGraphErrors *gr5 = new TGraphErrors(37,x5,y5,0,0);
  SetGr(gr5,"","","",3,3,1,1.0);//Green
  gr5->SetLineWidth(2);
  gr5->SetLineStyle(1);

// SLA W=2.110 GeV
  double x6[37] = {
  0.0000,
  5.0000,
 10.0000,
 15.0000,
 20.0000,
 25.0000,
 30.0000,
 35.0000,
 40.0000,
 45.0000,
 50.0000,
 55.0000,
 60.0000,
 65.0000,
 70.0000,
 75.0000,
 80.0000,
 85.0000,
 90.0000,
 95.0000,
100.0000,
105.0000,
110.0000,
115.0000,
120.0000,
125.0000,
130.0000,
135.0000,
140.0000,
145.0000,
150.0000,
155.0000,
160.0000,
165.0000,
170.0000,
175.0000,
180.0000};
  double y6[37] = {
0.2126E+03,
0.2118E+03,
0.2096E+03,
0.2065E+03,
0.2030E+03,
0.1995E+03,
0.1958E+03,
0.1917E+03,
0.1867E+03,
0.1802E+03,
0.1718E+03,
0.1612E+03,
0.1484E+03,
0.1334E+03,
0.1169E+03,
0.9930E+02,
0.8155E+02,
0.6453E+02,
0.4915E+02,
0.3620E+02,
0.2636E+02,
0.2002E+02,
0.1735E+02,
0.1819E+02,
0.2213E+02,
0.2853E+02,
0.3657E+02,
0.4540E+02,
0.5417E+02,
0.6217E+02,
0.6889E+02,
0.7407E+02,
0.7769E+02,
0.7996E+02,
0.8121E+02,
0.8179E+02,
0.8195E+02};

  for(int i=0;i<37;i++){  y6[i] /= 1000.; }
  TGraphErrors *gr6 = new TGraphErrors(37,x6,y6,0,0);
  SetGr(gr6,"","","",51,51,1,1.0);//Purple
  gr6->SetLineWidth(2);
  gr6->SetLineStyle(2);


double x7[4] = {//LEPS W=2.163 GeV
0.75,
0.85,
0.925,
0.975
};
double y7[4] = {//LEPS W=2.163 GeV
1.154,
1.127,
0.883,
0.761
};
double xe7[4] = {//LEPS W=2.163 GeV
0.,
0.,
0.,
0.
};
double ye7[4] = {//LEPS W=2.163 GeV
0.0522,
0.0414,
0.0406,
0.041
};
  for(int i=0;i<4;i++){ x7[i] = acos(x7[i])*180./PI; y7[i] /= 2*PI; ye7[i] /= 2*PI;}
  TGraphErrors *gr7 = new TGraphErrors(4,x7,y7,xe7,ye7);
  SetGr(gr7,"","","",6,6,24,1.5);//Purple






  TF1 *fx = new TF1("fx","x",0.,180.);
  TGaxis *ax = new TGaxis(0.,0, 180.,0, "fx", 505 );
  ax->SetLabelFont(42); ax->SetTitleFont(42); ax->CenterTitle();
  ax->SetLabelOffset(0.01);
  ax->SetTitleSize(0.06);  ax->SetTitleOffset(0.80);
  ax->SetMaxDigits(3);
  ax->SetTitle("#theta_{#gammaK}^{c.m.} [deg]");

  //TGaxis *ay = new TGaxis(-1,0, -1,0.5, 0,0.5, 505 );
  TGaxis *ay = new TGaxis(0,0, 0,0.5, 0,0.5, 505 );
  ay->SetLabelFont(42); ay->SetTitleFont(42); ay->CenterTitle();
  ay->SetLabelOffset(0.01);
  ay->SetTitleSize(0.06);  ay->SetTitleOffset(0.80);
  ay->SetTitle("d#sigma/d#Omega_{K}^{c.m.} [#mub/sr]");

  TLegend *leg = new TLegend(0.45,0.65,0.98,0.98,"","NDC");
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(3);
  leg->SetLineWidth(1);
  leg->SetTextSize(0.050);
  leg->SetMargin(0.15);
  leg->SetNColumns(1);
  leg->AddEntry(gr1     ,Form("SAPHIR #scale[0.6]{EPJA19(2004) 251.}  #scale[0.5]{(W=%.3lf GeV)}"
                        ,dataset1->GetKinematics()->GetW()/1000.),"pl");
  leg->AddEntry(gr2     ,Form("CLAS #scale[0.6]{PRC69(2004) 042201R.}  #scale[0.5]{(W=%.3lf GeV)}"
                        ,dataset2->GetKinematics()->GetW()/1000.),"pl");
  leg->AddEntry(gr3     ,Form("CLAS #scale[0.6]{PRC73(2006) 035202.}   #scale[0.5]{(W=%.3lf GeV)}"
                        ,dataset3->GetKinematics()->GetW()/1000.),"pl");
  leg->AddEntry(gr4     ,Form("CLAS #scale[0.6]{PRC82(2010) 025202.}   #scale[0.5]{(W=%.3lf GeV)}"
                        ,dataset4->GetKinematics()->GetW()/1000.),"pl");
  leg->AddEntry(gr7     ,"LEPS #scale[0.6]{PRC73(2006) 035214}    #scale[0.5]{(W=2.163 GeV)}","pl");
  leg->AddEntry(gr5     ,"KM #scale[0.6]{PRC61(1999) 012201.}       #scale[0.5]{(W=2.110 GeV)}","pl");
  leg->AddEntry(gr6     ,"SLA #scale[0.6]{PRC58(1998) 75.}          #scale[0.5]{(W=2.110 GeV)}","pl");


	double x_result[1], y_result[1], xel_result[1], xeh_result[1], yel_result[1], yeh_result[1], xe_result[1], ye_result[1];
	double x2_result[2], y2_result[2], xel2_result[2], xeh2_result[2], yel2_result[2], yeh2_result[2], xe2_result[2], ye2_result[2];
	double x3_result[3], y3_result[3], xe3_result[3], yel3_result[3], yeh3_result[3];
//	string result_in = "./result_input.dat";//D.C.S. Result
//	string buf;
//	int npoint = 0;
//	int npoint2 = 0;
//	int npoint3 = 0;
//	int flag_LS = -1;//Lambda(1) or Sigma0(2)
//	int flag_dep = -1;//Angle(1) or Q2(2)
//	int flag_div = -1;//1-Div.(1) or 2-Div.(2) or 3-Div.(3)
//	double xval, xerr, yval, yerr;
//	ifstream ifp(result_in.c_str(),ios::in);
//	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
//cout << "Param file : " << result_in.c_str() << endl;
//	while(1){
//		getline(ifp,buf);
//		if(buf[0]=='#'){continue;}
//		if(ifp.eof())break;
//		stringstream sbuf(buf);
//		sbuf >> flag_LS >> flag_dep >> flag_div >> xval >> yval >> xerr >> yerr;
//		cout << flag_LS << ", " << flag_dep << ", " << flag_div << ", " << xval << ", " << yval << ", " << xerr << ", " << yerr <<endl;
//		
//		if(flag_LS==2&&flag_dep==1&&flag_div==1){//1-Div.
//			x_result[npoint] = xval; 
//			y_result[npoint] = yval;
//			xe_result[npoint]= xerr; 
//			ye_result[npoint]= yerr; 
//			npoint++;
//		}
//		if(flag_LS==2&&flag_dep==1&&flag_div==2){//2-Div.
//			x2_result[npoint2] = xval; 
//			y2_result[npoint2] = yval;
//			xe2_result[npoint2]= xerr; 
//			ye2_result[npoint2]= yerr; 
//			npoint2++;
//		}
//		if(flag_LS==2&&flag_dep==1&&flag_div==3){//3-Div.
//			x3_result[npoint3] = xval; 
//			y3_result[npoint3] = yval;
//			xe3_result[npoint3]= xerr; 
//			ye3_result[npoint3]= yerr; 
//			npoint3++;
//		}
//	}
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%2021. 1. 11. ver., M-thesis		%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	x_result[0]=8.;//deg
	y_result[0]=0.0635;//ub/sr
	xe_result[0]=8.;
	xel_result[0]=8.;
	xeh_result[0]=8.;
	ye_result[0]=0.009;
	yel_result[0]=0.0183;
	yeh_result[0]=0.0626;

	x2_result[0]=4.;//deg
	y2_result[0]=0.0525;//ub/sr
	xe2_result[0]=4.;
	xel2_result[0]=4.;
	xeh2_result[0]=4.;
	ye2_result[0]=0.010;//Stat.
	yel2_result[0]=0.0177;//Stat. & Syst.
	yeh2_result[0]=0.0562;//Stat. & Syst.

	x2_result[1]=12.;//deg
	y2_result[1]=0.0743;//ub/sr
	xe2_result[1]=4.;
	xel2_result[1]=4.;
	xeh2_result[1]=4.;
	ye2_result[1]=0.0123;
	yel2_result[1]=0.0239;
	yeh2_result[1]=0.0774;
//
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%2021. 1. 8. ver., VP Flux was wrong%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//	x_result[0]=8.;//deg
//	y_result[0]=0.0676;//ub/sr
//	xe_result[0]=8.;
//	xel_result[0]=8.;
//	xeh_result[0]=8.;
//	ye_result[0]=0.010;
//	yel_result[0]=0.0195;
//	yeh_result[0]=0.0666;
//
//	x2_result[0]=4.;//deg
//	y2_result[0]=0.0546;//ub/sr
//	ye2_result[0]=0.00998;
//	xe2_result[0]=4.;
//	xel2_result[0]=4.;
//	xeh2_result[0]=4.;
//	yel2_result[0]=0.0184;
//	yeh2_result[0]=0.0585;
//
//	x2_result[1]=12.;//deg
//	y2_result[1]=0.0805;//ub/sr
//	xe2_result[1]=4.;
//	xel2_result[1]=4.;
//	xeh2_result[1]=4.;
//	ye2_result[1]=0.0133;
//	yel2_result[1]=0.0259;
//	yeh2_result[1]=0.0838;
//	x3_result[0]=3.;//deg
//	y3_result[0]=0.0687;//ub/sr
//	xe3_result[0]=3.;
//	ye3_result[0]=0.008;
//	x3_result[1]=8.;//deg
//	y3_result[1]=0.1169;//ub/sr
//	xe3_result[1]=2.;
//	ye3_result[1]=0.011;
//	x3_result[2]=13.;//deg
//	y3_result[2]=0.1205;//ub/sr
//	xe3_result[2]=3.;
//	ye3_result[2]=0.011;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

  //TGraphErrors *gr_result = new TGraphErrors(1, x_result, y_result, xe_result, ye_result);
  //gr_result->SetMarkerSize(1.5);
  //gr_result->SetMarkerStyle(21);
  //gr_result->SetFillStyle(3004);
  ////gr_result->SetLineColor(2);
  ////gr_result->SetFillColor(2);
  //gr_result->SetMarkerColor(2);
  //TGraphErrors *gr2_result = new TGraphErrors(2, x2_result, y2_result, xe2_result, ye2_result);
  //gr2_result->SetMarkerSize(1.5);
  //gr2_result->SetMarkerStyle(21);
  //gr2_result->SetFillStyle(3005);
  ////gr2_result->SetLineColor(kViolet);
  ////gr2_result->SetFillColor(103);
  //gr2_result->SetMarkerColor(3);
  //TGraphAsymmErrors *gr_result = new TGraphAsymmErrors(1, x_result, y_result, xel_result, xeh_result, yel_result, yeh_result);
  TGraphErrors *gr_result = new TGraphErrors(1, x_result, y_result, xe_result, ye_result);
  gr_result->SetMarkerSize(1.5);
  gr_result->SetMarkerStyle(21);
  //gr_result->SetFillStyle(3004);
  gr_result->SetLineColor(2);
  gr_result->SetLineWidth(2);
  //gr_result->SetFillColor(2);
  gr_result->SetMarkerColor(2);
  //TGraphAsymmErrors *gr2_result = new TGraphAsymmErrors(2, x2_result, y2_result, xel2_result, xeh2_result, yel2_result, yeh2_result);
  TGraphErrors *gr2_result = new TGraphErrors(2, x2_result, y2_result, xe2_result, ye2_result);
  gr2_result->SetMarkerSize(1.5);
  gr2_result->SetMarkerStyle(21);
  //gr2_result->SetFillStyle(3005);
  gr2_result->SetLineColor(kAzure);
  gr2_result->SetLineWidth(2);
  //gr2_result->SetFillColor(103);
  gr2_result->SetMarkerColor(4);
//  TGraphErrors *gr3_result = new TGraphErrors(3, x3_result, y3_result, xe3_result, ye3_result);
//  gr3_result->SetMarkerSize(1.5);
//  gr3_result->SetMarkerStyle(21);
//  gr3_result->SetFillStyle(3004);
//  //gr3_result->SetLineColor(kOrange);
//  //gr3_result->SetFillColor(102);
//  gr3_result->SetMarkerColor(4);

  TCanvas *c1 = new TCanvas("c1","c1",1100,800);
  c1->Divide(1,1,0.0001,0.0001);
  c1->cd(1);
  h2->Draw();
  gr1->Draw("Psame");
  gr2->Draw("Psame");
  gr3->Draw("Psame");
  gr4->Draw("Psame");
  gr7->Draw("Psame");//LEPS
  gr5->Draw("Lsame");
  gr6->Draw("Lsame");
//TBox *b = new TBox(x_result[0]-xe_result[0],y_result[0]-ye_result[0],x_result[0]+xe_result[0],y_result[0]+ye_result[0]); 
//	b->SetFillColor(2); 
//	b->SetFillStyle(0);
//	b->SetLineColor(2);
//	b->SetLineWidth(1);
//TBox *b2_1 = new TBox(x2_result[0]-xe2_result[0],y2_result[0]-ye2_result[0],x2_result[0]+xe2_result[0],y2_result[0]+ye2_result[0]); 
////TBox *b2_1 = new TBox(0.,0.26527,8.,0.30327); 
//	b2_1->SetFillColor(3); 
//	b2_1->SetFillStyle(0);
//	b2_1->SetLineColor(3);
//	b2_1->SetLineWidth(1);
////TBox *b2_2 = new TBox(8.,0.31396,16.,0.35196); 
//TBox *b2_2 = new TBox(x2_result[1]-xe2_result[1],y2_result[1]-ye2_result[1],x2_result[1]+xe2_result[1],y2_result[1]+ye2_result[1]); 
//	b2_2->SetFillColor(3); 
//	b2_2->SetFillStyle(0);
//	b2_2->SetLineColor(3);
//	b2_2->SetLineWidth(1);
TBox *b = new TBox(x_result[0]-xel_result[0],y_result[0]-yel_result[0],x_result[0]+xeh_result[0],y_result[0]+yeh_result[0]); 
	b->SetFillColor(2); 
	b->SetFillStyle(0);
	b->SetLineColor(2);
	b->SetLineWidth(3);
	b->SetLineStyle(2);
TBox *b_stat = new TBox(x_result[0]-xe_result[0],y_result[0]-ye_result[0],x_result[0]+xe_result[0],y_result[0]+ye_result[0]); 
	b_stat->SetFillColor(2); 
	b_stat->SetFillStyle(0);
	b_stat->SetLineColor(2);
	b_stat->SetLineWidth(3);
TBox *b2_1 = new TBox(x2_result[0]-xel2_result[0],y2_result[0]-yel2_result[0],x2_result[0]+xeh2_result[0],y2_result[0]+yeh2_result[0]); 
//TBox *b2_1 = new TBox(0.,0.26527,8.,0.30327); 
	b2_1->SetFillColor(4); 
	b2_1->SetFillStyle(0);
	b2_1->SetLineColor(4);
	b2_1->SetLineWidth(3);
	b2_1->SetLineStyle(2);
TBox *b2_stat1 = new TBox(x2_result[0]-xe2_result[0],y2_result[0]-ye2_result[0],x2_result[0]+xe2_result[0],y2_result[0]+ye2_result[0]); 
//TBox *b2_1 = new TBox(0.,0.26527,8.,0.30327); 
	b2_stat1->SetFillColor(4); 
	b2_stat1->SetFillStyle(0);
	b2_stat1->SetLineColor(4);
	b2_stat1->SetLineWidth(3);
//TBox *b2_2 = new TBox(8.,0.31396,16.,0.35196); 
TBox *b2_2 = new TBox(x2_result[1]-xel2_result[1],y2_result[1]-yel2_result[1],x2_result[1]+xeh2_result[1],y2_result[1]+yeh2_result[1]); 
	b2_2->SetFillColor(4); 
	b2_2->SetFillStyle(0);
	b2_2->SetLineColor(4);
	b2_2->SetLineWidth(3);
	b2_2->SetLineStyle(2);
TBox *b2_stat2 = new TBox(x2_result[1]-xe2_result[1],y2_result[1]-ye2_result[1],x2_result[1]+xe2_result[1],y2_result[1]+ye2_result[1]); 
	b2_stat2->SetFillColor(4); 
	b2_stat2->SetFillStyle(0);
	b2_stat2->SetLineColor(4);
	b2_stat2->SetLineWidth(3);
//TBox *b3_1 = new TBox(0.,0.23115,6.,0.28115); 
//TBox *b3_1 = new TBox(x3_result[0]-xe3_result[0],y3_result[0]-ye3_result[0],x3_result[0]+xe3_result[0],y3_result[0]+ye3_result[0]); 
//	b3_1->SetFillColor(4); 
//	b3_1->SetFillStyle(0);
//	b3_1->SetLineColor(4);
//	b3_1->SetLineWidth(1);
////TBox *b3_2 = new TBox(6.,0.31818,10.,0.36818); 
//TBox *b3_2 = new TBox(x3_result[1]-xe3_result[1],y3_result[1]-ye3_result[1],x3_result[1]+xe3_result[1],y3_result[1]+ye3_result[1]); 
//	b3_2->SetFillColor(4); 
//	b3_2->SetFillStyle(0);
//	b3_2->SetLineColor(4);
//	b3_2->SetLineWidth(1);
////TBox *b3_3 = new TBox(10.,0.30168,16.,0.35168); 
//TBox *b3_3 = new TBox(x3_result[2]-xe3_result[2],y3_result[2]-ye3_result[2],x3_result[2]+xe3_result[2],y3_result[2]+ye3_result[2]); 
//	b3_3->SetFillColor(4); 
//	b3_3->SetFillStyle(0);
//	b3_3->SetLineColor(4);
//	b3_3->SetLineWidth(1);
//  gr3_result->Draw("P2same");
  //gr2_result->Draw("Psame");
//  gr_result->Draw("Psame");
//  b3_1->Draw();
//  b3_2->Draw();
//  b3_3->Draw();
 // b2_1->Draw();
 // b2_2->Draw();
 // b2_stat1->Draw();
 // b2_stat2->Draw();
 // b->Draw();
  //b_stat->Draw();
  ax->Draw();
  ay->Draw();
//  leg->Draw();
//for(int i=1;i<200;i++){
//dataset1->SetSelection(i);
//cout<<i<<": "<<dataset1->GetKinematics()->GetW()<<endl;
//dataset2->SetSelection(i);
//cout<<i<<": "<<dataset2->GetKinematics()->GetW()<<endl;
//dataset3->SetSelection(i);
//cout<<i<<": "<<dataset3->GetKinematics()->GetW()<<endl;
//dataset4->SetSelection(i);
//cout<<i<<": "<<dataset4->GetKinematics()->GetW()<<endl;
//}
//c1->Print("Elementary_gpKS.pdf");
//c1->Print("/data/41a/ELS/okuyama/JLab_nnL/okuya_macros/mthesis_Fig/pdf/CS_thetadepS2.pdf");
//c1->Print("/data/41a/ELS/okuyama/JLab_nnL/okuya_macros/mthesis_Fig/pdf/symp_gpKS.pdf");
}
