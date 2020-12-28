const double PI = 4.*atan(1.);
string dir = "/usr/local/apps/strangecalc/share/strangecalc-wrapper/data/";

void SetTH2(TH2F *h2, TString hname, TString xname, TString yname, double offset, double Tsize){
  h2->SetTitle(hname);
  h2->GetXaxis()->SetTitle(xname);
  h2->GetXaxis()->CenterTitle();
//  h2->GetXaxis()->SetTitleFont(42);
//  h2->GetXaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetLabelOffset(offset);
  h2->GetXaxis()->SetLabelSize(Tsize);
  h2->GetXaxis()->SetTitleSize(Tsize);

  h2->GetYaxis()->SetTitle(yname);
  h2->GetYaxis()->CenterTitle();
//  h2->GetYaxis()->SetTitleFont(42);
//  h2->GetYaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetTitleOffset(offset);
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
  ifp = new TFile(Form("%sdata_iso1.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("CLAS_2006__Bradford__PRC_73_035202",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<"  "<<dataset->GetEntries()<<endl;
  ifp->GetObject("LEPS_2007__Hicks__PRC_76_042201__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<"  "<<dataset->GetEntries()<<endl;
  ifp->GetObject("SAPHIR_2004__Glander__EPJA_19_251",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<"  "<<dataset->GetEntries()<<endl;
  ifp->GetObject("CLAS_2009__McCracken__PRC_81_025201__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<"  "<<dataset->GetEntries()<<endl;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void clas_2006(){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso1.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("CLAS_2006__Bradford__PRC_73_035202",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<endl;
//  dataset->ViewSelections();

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,10);

  TGraph *gr[79];
  TLatex *t1[79], *t2[79];
  for(int i=0;i<79;i++){
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
  for(int i=0;i<29;i++){ c2->cd(i+1); h2->Draw(); gr[i+50]->Draw("Psame"); t1[i+50]->Draw("same"); t2[i]->Draw("same"); }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void leps_2007(){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso1.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("LEPS_2007__Hicks__PRC_76_042201__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<endl;
//  dataset->ViewSelections();

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,20);

  TGraph *gr[6];
  TLatex *t1[6], *t2[6];
  for(int i=0;i<6;i++){
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
  for(int i=0;i<6;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void safhir_2004(){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso1.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("SAPHIR_2004__Glander__EPJA_19_251",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<endl;
//  dataset->ViewSelections();

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,10);

  TGraph *gr[36];
  TLatex *t1[36], *t2[36];
  for(int i=0;i<36;i++){
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
  for(int i=0;i<36;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void clas_2009(){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso1.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("CLAS_2009__McCracken__PRC_81_025201__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<endl;
//  dataset->ViewSelections();

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,10);

  TGraphErrors *gr[119];
  TLatex *t1[119], *t2[119];
  for(int i=0;i<119;i++){
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

//    for(int j=gr[i]->GetN()-1;j>=0;j--){ cout<<Form("%.5lf, ",gr[i]->GetErrorY(j))<<flush; }
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
  for(int i=0;i<19;i++){ c3->cd(i+1); h2->Draw(); gr[i+100]->Draw("Psame"); t1[i+100]->Draw("same"); t2[i+100]->Draw("same"); }
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
  ifp = new TFile(Form("%sdata_iso1.root",dir.c_str()));

  string ofname = "Exp_KpLambda.pdf";
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
  ifp->GetObject("SAPHIR_2004__Glander__EPJA_19_251",dataset);
  for(int i=0;i<36;i++){
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
  for(int i=0;i<16;i++){ c1->cd(i+1); h2->Draw(); gr[i+20]->Draw("Psame"); t1[i+20]->Draw("same"); t2[i+20]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();

/////////
  ifp->GetObject("CLAS_2006__Bradford__PRC_73_035202",dataset);
  for(int i=0;i<79;i++){
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
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame");    t1[i]->Draw("same");    t2[i]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i+20]->Draw("Psame"); t1[i+20]->Draw("same"); t2[i+20]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i+40]->Draw("Psame"); t1[i+40]->Draw("same"); t2[i+40]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<19;i++){ c1->cd(i+1); h2->Draw(); gr[i+60]->Draw("Psame"); t1[i+60]->Draw("same"); t2[i+60]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();

/////////
  ifp->GetObject("CLAS_2009__McCracken__PRC_81_025201__dcs",dataset);
  for(int i=0;i<119;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.45,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.40,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08); t1[i]->SetTextColor(4);
    t2[i]->SetTextSize(0.08); t2[i]->SetTextColor(4);
  }
  title = p1->AddText("CLAS 2009"); p1->Draw();
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
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i+60]->Draw("Psame");  t1[i+60]->Draw("same"); t2[i+60]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i+80]->Draw("Psame");  t1[i+80]->Draw("same"); t2[i+80]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<19;i++){ c1->cd(i+1); h2->Draw(); gr[i+100]->Draw("Psame"); t1[i+100]->Draw("same"); t2[i+100]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();

/////////
  ifp->GetObject("LEPS_2007__Hicks__PRC_76_042201__dcs",dataset);
  for(int i=0;i<6;i++){
    dataset->SetSelection(i+1);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.45,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.40,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08); t1[i]->SetTextColor(4);
    t2[i]->SetTextSize(0.08); t2[i]->SetTextColor(4);
  }
  title = p1->AddText("LEPS 2007"); p1->Draw();
  c1->Print(Form("%s",ofname.c_str()));  p1->Clear();  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<6;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();

  c1->Print(Form("%s]",ofname.c_str()));
  gSystem->Exit(1);

}

////////////////////////////////////////////////////////////////////////////////////////////
void draw_KpLambda(){
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
  ifp = new TFile(Form("%sdata_iso1.root",dir.c_str()));
  TDataset *dataset1, *dataset2, *dataset3;
  ifp->GetObject("CLAS_2009__McCracken__PRC_81_025201__dcs",dataset1);
  ifp->GetObject("CLAS_2006__Bradford__PRC_73_035202"      ,dataset2);
  ifp->GetObject("SAPHIR_2004__Glander__EPJA_19_251"       ,dataset3);

  //TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.5);
  TH2F *h2 = new TH2F("h2","h2",1000,0.,180.,1000,0,0.5);
  h2->SetStats(kFALSE);
//  SetTH2_forDraw(h2,"","cos(#theta_{K^{+}}^{c.m.})","d#sigma/d#Omega_{K}^{c.m.} [#mub/sr]");
  SetTH2_forDraw(h2,"","","");
  h2->GetXaxis()->SetLabelOffset(1.);
  h2->GetXaxis()->SetTickLength(0);
  h2->GetYaxis()->SetLabelOffset(1.);
  h2->GetYaxis()->SetLabelSize(0.);
  h2->GetYaxis()->SetTickLength(0);

// CLAS 2009
  TGraphErrors *gr1;
  //dataset1->SetSelection(22);//W=1.835
  dataset1->SetSelection(50);
  gr1 = dataset1->MakeGraph("costhkcm");
  double *x1 = gr1->GetX(), *y1 = gr1->GetY(), ey1[100];
    //for(int i=0;i<gr1->GetN();i++){ x1[i] *= -1; ey1[i] = gr1->GetErrorY(i); }
    for(int i=0;i<gr1->GetN();i++){ x1[i] = acos(x1[i])*180./PI; ey1[i] = gr1->GetErrorY(i); }
  gr1->Clear();
  gr1 = new TGraphErrors(gr1->GetN(),x1,y1,0,ey1);
  SetGr(gr1,"","","",1,1,20,1.5);//Black

  TGraphErrors *gr1_1;
  gr1_1 = dataset1->MakeGraph("wlab");
  double *x1_1 = gr1_1->GetX(), *y1_1 = gr1_1->GetY(), ey1_1[100];
    for(int i=0;i<gr1_1->GetN();i++){ cout<<i<<"  "<<x1_1[i]<<"  "<<y1_1[i]<<endl; x1_1[i] *= +1; ey1_1[i] = gr1_1->GetErrorY(i); }
  gr1_1->Clear();
  gr1_1 = new TGraphErrors(gr1_1->GetN(),x1_1,y1_1,0,ey1_1);
  SetGr(gr1_1,"","","",1,1,20,1.5);//Black

// CLAS 2006
  TGraphErrors *gr2;
  //dataset2->SetSelection(16);//W=1.835
  dataset2->SetSelection(41);
  gr2 = dataset2->MakeGraph("costhkcm");
  double *x2 = gr2->GetX(), *y2 = gr2->GetY(), ey2[100];
    for(int i=0;i<gr2->GetN();i++){ x2[i] = acos(x2[i])*180./PI; ey2[i] = gr2->GetErrorY(i); }
  gr2->Clear();
  gr2 = new TGraphErrors(gr2->GetN(),x2,y2,0,ey2);
  SetGr(gr2,"","","",17,17,21,1.5);//Grey

// SAPHIR 2004
  TGraphErrors *gr3;
  //dataset3->SetSelection(11);//W=1.835
  dataset3->SetSelection(23);
  gr3 = dataset3->MakeGraph("costhkcm");
  double *x3 = gr3->GetX(), *y3 = gr3->GetY(), ey3[100];
    for(int i=0;i<gr3->GetN();i++){ x3[i] = acos(x3[i])*180./PI; ey3[i] = gr3->GetErrorY(i); }
  gr3->Clear();
  gr3 = new TGraphErrors(gr3->GetN(),x3,y3,0,ey3);
  SetGr(gr3,"","","",4,4,24,1.5);//Blue

// KMaid W=2.110 GeV
  double x4[37] = {
//      0.,  5., 10., 15., 20., 25., 30., 35., 40., 45.,
//     50., 55., 60., 65., 70., 75., 80., 85., 90., 95.,
//    100.,105.,110.,115.,120.,125.,130.,135.,140.,145.,
//    150.,155.,160.,165.,170.,175.,180. };// KMaid W=1.835 GeV
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
  double y4[37] = {
//    0.2806E+03,0.2819E+03,0.2855E+03,0.2905E+03,0.2955E+03,0.2990E+03,0.2998E+03,0.2968E+03,0.2897E+03,0.2786E+03,
//    0.2640E+03,0.2467E+03,0.2279E+03,0.2084E+03,0.1894E+03,0.1715E+03,0.1553E+03,0.1412E+03,0.1292E+03,0.1193E+03,
//    0.1114E+03,0.1050E+03,0.9979E+02,0.9537E+02,0.9132E+02,0.8727E+02,0.8294E+02,0.7814E+02,0.7279E+02,0.6694E+02,
//    0.6079E+02,0.5460E+02,0.4874E+02,0.4363E+02,0.3964E+02,0.3710E+02,0.3623E+02 };//W=1.835 GeV
0.5361E+02,
0.6634E+02,
0.1012E+03,
0.1492E+03,
0.1990E+03,
0.2399E+03,
0.2649E+03,
0.2716E+03,
0.2618E+03,
0.2395E+03,
0.2101E+03,
0.1781E+03,
0.1471E+03,
0.1195E+03,
0.9634E+02,
0.7785E+02,
0.6375E+02,
0.5341E+02,
0.4610E+02,
0.4110E+02,
0.3774E+02,
0.3549E+02,
0.3393E+02,
0.3276E+02,
0.3181E+02,
0.3102E+02,
0.3040E+02,
0.3003E+02,
0.3001E+02,
0.3044E+02,
0.3137E+02,
0.3275E+02,
0.3445E+02,
0.3623E+02,
0.3781E+02,
0.3889E+02,
0.3928E+02};
  //for(int i=0;i<37;i++){ x4[i] = -cos(x4[i]/180.*PI); y4[i] /= 1000.; }
  for(int i=0;i<37;i++){ y4[i] /= 1000.; }
  TGraphErrors *gr4 = new TGraphErrors(37,x4,y4,0,0);
  SetGr(gr4,"","","",3,3,1,1.0);//Green
  gr4->SetLineWidth(2);
  gr4->SetLineStyle(1);

// SLA W=2.110 GeV
  double x5[37] = {
//      0.,  5., 10., 15., 20., 25., 30., 35., 40., 45.,
//     50., 55., 60., 65., 70., 75., 80., 85., 90., 95.,
//    100.,105.,110.,115.,120.,125.,130.,135.,140.,145.,
//    150.,155.,160.,165.,170.,175.,180. };// SLA W=1.835 GeV
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
//    0.3981E+03,0.3944E+03,0.3838E+03,0.3681E+03,0.3495E+03,0.3301E+03,0.3112E+03,0.2936E+03,0.2772E+03,0.2619E+03,
//    0.2473E+03,0.2331E+03,0.2192E+03,0.2054E+03,0.1919E+03,0.1788E+03,0.1662E+03,0.1543E+03,0.1431E+03,0.1326E+03,
//    0.1229E+03,0.1137E+03,0.1050E+03,0.9654E+02,0.8825E+02,0.7998E+02,0.7162E+02,0.6316E+02,0.5464E+02,0.4618E+02,
//    0.3798E+02,0.3029E+02,0.2341E+02,0.1765E+02,0.1329E+02,0.1057E+02,0.9650E+01 };
0.4193E+03,
0.4118E+03,
0.3925E+03,
0.3675E+03,
0.3424E+03,
0.3191E+03,
0.2964E+03,
0.2723E+03,
0.2456E+03,
0.2160E+03,
0.1848E+03,
0.1537E+03,
0.1246E+03,
0.9954E+02,
0.7980E+02,
0.6640E+02,
0.5980E+02,
0.5998E+02,
0.6652E+02,
0.7869E+02,
0.9548E+02,
0.1157E+03,
0.1381E+03,
0.1613E+03,
0.1839E+03,
0.2049E+03,
0.2231E+03,
0.2378E+03,
0.2484E+03,
0.2548E+03,
0.2573E+03,
0.2564E+03,
0.2533E+03,
0.2491E+03,
0.2449E+03,
0.2419E+03,
0.2409E+03};
  //for(int i=0;i<37;i++){ x5[i] = -cos(x5[i]/180.*PI); y5[i] /= 1000.; }
  for(int i=0;i<37;i++){ y5[i] /= 1000.; }
  TGraphErrors *gr5 = new TGraphErrors(37,x5,y5,0,0);
  SetGr(gr5,"","","",2,2,1,1.0);//Red
  gr5->SetLineWidth(2);
  gr5->SetLineStyle(2);

// RPR3 w=2.130 GeV
  double x6[38] = { // cos(theta_CM)
//  +1.00000, +0.99619, +0.98481, +0.96593, +0.93969, +0.90631, +0.86603, +0.81915, +0.76604, +0.70711,
//  +0.64279, +0.57358, +0.50000, +0.42262, +0.34202, +0.25882, +0.17365, +0.08716, +0.00000, -0.08716,
//  -0.17365, -0.25882, -0.34202, -0.42262, -0.50000, -0.57358, -0.64279, -0.70711, -0.76604, -0.81915,
//  -0.86603, -0.90631, -0.93969, -0.96593, -0.98481, -0.99619, -1.00000 };// RPR3 w=1.835 GeV
1		,	
0.996397,	
0.985616,	
0.967733,	
0.942877,	
0.911228,	
0.873014,	
0.82851	,	
0.778036,	
0.721956,	
0.660675,	
0.594633,	
0.524307,	
0.450204,	
0.372856,	
0.292823,	
0.210679,	
0.127018,	
0.0424412,	
-0.0424412,
-0.127018,	
-0.210679,	
-0.292823,	
-0.372856,	
-0.450204,	
-0.524307,	
-0.594633,	
-0.660675,	
-0.721956,	
-0.778036,	
-0.82851,	
-0.873014,	
-0.911228,	
-0.942877,	
-0.967733,	
-0.985616,	
-0.996397,	
-1};

  double y6[38] = {
//    0.43543, 0.39347, 0.36429, 0.34223, 0.32417, 0.30840, 0.29399, 0.28045, 0.26759, 0.25533,
//    0.24371, 0.23276, 0.22257, 0.21319, 0.20468, 0.19708, 0.19040, 0.18464, 0.17977, 0.17575,
//    0.17252, 0.16999, 0.16807, 0.16664, 0.16559, 0.16478, 0.16407, 0.16332, 0.16237, 0.16108,
//    0.15930, 0.15688, 0.15368, 0.14955, 0.14437, 0.13803, 0.13040 };
0.218105,
0.216723,
0.214651,
0.216062,
0.223337,
0.235061,
0.246988,
0.254359,
0.253793,
0.244075,
0.226018,
0.201864,
0.174575,
0.147218,
0.122506,
0.102502,
0.0884764,
0.080866,
0.0793447,
0.0829619,
0.0903386,
0.0998901,
0.110046,
0.119435,
0.127022,
0.132167,
0.134627,
0.134502,
0.132153,
0.128099,
0.12293,
0.117232,
0.111537,
0.106299,
0.101875,
0.098534,
0.096459,
0.0957557};
  //for(int i=0;i<37;i++){ x6[i] = -x6[i]; }//cos
  for(int i=0;i<38;i++){ x6[i] = acos(x6[i])*180./PI; }
  TGraphErrors *gr6 = new TGraphErrors(38,x6,y6,0,0);
  SetGr(gr6,"","","",1,1,1,1.0); gr6->SetLineWidth(2);
  gr6->SetLineStyle(3);

// RPR2011 w=2.13 GeV
  double x7[38] = { // cos(theta_CM)
  //+1.00000, +0.99619, +0.98481, +0.96593, +0.93969, +0.90631, +0.86603, +0.81915, +0.76604, +0.70711,
  //+0.64279, +0.57358, +0.50000, +0.42262, +0.34202, +0.25882, +0.17365, +0.08716, +0.00000, -0.08716,
  //-0.17365, -0.25882, -0.34202, -0.42262, -0.50000, -0.57358, -0.64279, -0.70711, -0.76604, -0.81915,
  //-0.86603, -0.90631, -0.93969, -0.96593, -0.98481, -0.99619, -1.00000 };// RPR2011 w=1.835 GeV
1		,
0.996397,
0.985616,
0.967733,
0.942877,
0.911228,
0.873014,
0.82851	,
0.778036,
0.721956,
0.660675,
0.594633,
0.524307,
0.450204,
0.372856,
0.292823,
0.210679,
0.127018,
0.0424412,
-0.0424412,
-0.127018,
-0.210679,
-0.292823,
-0.372856,
-0.450204,
-0.524307,
-0.594633,
-0.660675,
-0.721956,
-0.778036,
-0.82851 ,
-0.873014,
-0.911228, 
-0.942877,
-0.967733,
-0.985616,
-0.996397,
-1};


  double y7[38] = {
  //  0.53054, 0.45499, 0.40082, 0.36112, 0.33152, 0.30910, 0.29176, 0.27798, 0.26660, 0.25674,
  //  0.24773, 0.23905, 0.23037, 0.22144, 0.21216, 0.20250, 0.19253, 0.18236, 0.17216, 0.16216,
  //  0.15256, 0.14359, 0.13547, 0.12836, 0.12239, 0.11762, 0.11404, 0.11154, 0.10988, 0.10873,
  //  0.10762, 0.10592, 0.10286, 0.09751, 0.08876, 0.07535, 0.05581 };
0.465884 ,	
0.455459 ,
0.427822 ,
0.391372 ,
0.354215 ,
0.32089  ,
0.292313 ,
0.267452 ,   
0.244896 ,
0.223594 ,
0.202926 ,
0.182519 ,
0.162098 ,
0.14148  ,
0.120677 ,
0.100002 ,
0.0801058,
0.0618982,
0.0463706,
0.0343699,
0.0263968,
0.022483 ,
0.0221773,
0.0246333,
0.0287656,
0.0334326,
0.0376024,
0.0404756,
0.0415551,
0.0406632,
0.0379187,
0.0336857,
0.0285048,
0.0230174,
0.0178863,
0.0137208,
0.01101  ,
0.0100699};

  //for(int i=0;i<37;i++){ x7[i] = -x7[i]; }
  for(int i=0;i<38;i++){ x7[i] = acos(x7[i])*180./PI; }
  TGraphErrors *gr7 = new TGraphErrors(38,x7,y7,0,0);
  SetGr(gr7,"","","",4,4,1,1.0);//Blue
  gr7->SetLineWidth(2);
  gr7->SetLineStyle(4);


double x8[4] = {//LEPS W=2.163 GeV
0.75,
0.85,
0.925,
0.975
};
double y8[4] = {//LEPS W=2.163 GeV
1.426,
1.679,
1.728,
1.945
};
double xe8[4] = {//LEPS W=2.163 GeV
0.,
0.,
0.,
0.
};
double ye8[4] = {//LEPS W=2.163 GeV
0.0578,
0.0506,
0.0565,
0.066
};
  for(int i=0;i<4;i++){ x8[i] = acos(x8[i])*180./PI; y8[i] /= 2*PI; ye8[i] /= 2*PI;}
  TGraphErrors *gr8 = new TGraphErrors(4,x8,y8,xe8,ye8);
  SetGr(gr8,"","","",6,6,24,1.5);//Purple



  //TF1 *fx = new TF1("fx","-x",-1,1);
  TF1 *fx = new TF1("fx","x",0.,180.);
  TGaxis *ax = new TGaxis(0.,0, 180.,0, "fx", 505 );
  //TGaxis *ax = new TGaxis(-1,0, 1,0, "fx", 505 );
  ax->SetLabelFont(42); ax->SetTitleFont(42); ax->CenterTitle();
  ax->SetLabelOffset(0.01);
  ax->SetTitleSize(0.06);  ax->SetTitleOffset(0.80);
  ax->SetMaxDigits(3);
  //ax->SetTitle("cos(#theta_{K}^{c.m.})");
  ax->SetTitle("#theta_{#gammaK}^{CM} [deg]");

  TGaxis *ay = new TGaxis(0,0, 0,0.5, 0,0.5, 505 );
  //TGaxis *ay = new TGaxis(-1,0, -1,0.5, 0,0.5, 505 );
  ay->SetLabelFont(42); ay->SetTitleFont(42); ay->CenterTitle();
  ay->SetLabelOffset(0.01);
  ay->SetTitleSize(0.06);  ay->SetTitleOffset(0.80);
  ay->SetTitle("d#sigma/d#Omega_{K}^{CM} [#mub/sr]");

  TLegend *leg = new TLegend(0.40,0.60,0.98,0.98,"","NDC");
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(3);
  leg->SetLineWidth(1);
  leg->SetTextSize(0.050);
  leg->SetNColumns(1);
  leg->SetMargin(0.15);
  leg->AddEntry(gr3     ,Form("SAPHIR #scale[0.6]{EPJA19(2004) 251.}      #scale[0.5]{(#sqrt{s}=%.3lf GeV)}"
                        ,dataset3->GetKinematics()->GetW()/1000.),"pl");
  leg->AddEntry(gr2     ,Form("CLAS  #scale[0.6]{PRC73(2006) 035202.}      #scale[0.5]{(#sqrt{s}=%.3lf GeV)}"
                        ,dataset2->GetKinematics()->GetW()/1000.),"pl");
  leg->AddEntry(gr1     ,Form("CLAS  #scale[0.6]{PRC81(2010) 025201.}      #scale[0.5]{(#sqrt{s}=%.3lf GeV)}"
                        ,dataset1->GetKinematics()->GetW()/1000.),"pl");
  leg->AddEntry(gr8     ,"LEPS #scale[0.6]{PRC73(2006) 035214}        #scale[0.5]{(#sqrt{s}=2.163 GeV)}","pl");
  leg->AddEntry(gr5     ,"SLA    #scale[0.6]{PRC58(1998) 75.}           #scale[0.5]{(#sqrt{s}=2.110 GeV)}","pl");
  leg->AddEntry(gr4     ,"MAID #scale[0.6]{PRC61(1999) 012201.}       #scale[0.5]{(#sqrt{s}=2.110 GeV)}","pl");
  leg->AddEntry(gr6     ,"RPR3  #scale[0.6]{PRC73(2006) 045207.}      #scale[0.5]{(#sqrt{s}=2.130 GeV)}","pl");
  leg->AddEntry(gr7     ,"RPR2011 #scale[0.6]{PRL108(2012) 182002.}  #scale[0.5]{(#sqrt{s}=2.130 GeV)}","pl");


	double x_result[1], y_result[1], xe_result[1], ye_result[1];
	double x2_result[2], y2_result[2], xe2_result[2], ye2_result[2];
	double x3_result[3], y3_result[3], xe3_result[3], ye3_result[3];
	x_result[0]=8.;//deg
	y_result[0]=0.29291;//ub/sr
	xe_result[0]=8.;
	ye_result[0]=0.0085;
	x2_result[0]=4.;//deg
	y2_result[0]=0.28427;//ub/sr
	xe2_result[0]=4.;
	ye2_result[0]=0.019;
	x2_result[1]=12.;//deg
	y2_result[1]=0.33296;//ub/sr
	xe2_result[1]=4.;
	ye2_result[1]=0.019;
	x3_result[0]=3.;//deg
	y3_result[0]=0.25615;//ub/sr
	xe3_result[0]=3.;
	ye3_result[0]=0.025;
	x3_result[1]=8.;//deg
	y3_result[1]=0.34318;//ub/sr
	xe3_result[1]=2.;
	ye3_result[1]=0.025;
	x3_result[2]=13.;//deg
	y3_result[2]=0.32668;//ub/sr
	xe3_result[2]=3.;
	ye3_result[2]=0.025;
  TGraphErrors *gr_result = new TGraphErrors(1, x_result, y_result, xe_result, ye_result);
  gr_result->SetMarkerSize(1.5);
  gr_result->SetMarkerStyle(21);
  gr_result->SetFillStyle(3004);
  //gr_result->SetLineColor(2);
  //gr_result->SetFillColor(2);
  gr_result->SetMarkerColor(2);
  TGraphErrors *gr2_result = new TGraphErrors(2, x2_result, y2_result, xe2_result, ye2_result);
  gr2_result->SetMarkerSize(1.5);
  gr2_result->SetMarkerStyle(21);
  gr2_result->SetFillStyle(3005);
  //gr2_result->SetLineColor(kViolet);
  //gr2_result->SetFillColor(103);
  gr2_result->SetMarkerColor(3);
  TGraphErrors *gr3_result = new TGraphErrors(3, x3_result, y3_result, xe3_result, ye3_result);
  gr3_result->SetMarkerSize(1.5);
  gr3_result->SetMarkerStyle(21);
  gr3_result->SetFillStyle(3004);
  //gr3_result->SetLineColor(kOrange);
  //gr3_result->SetFillColor(102);
  gr3_result->SetMarkerColor(4);

  TCanvas *c1 = new TCanvas("c1","c1",1100,800);
  c1->Divide(1,1,0.0001,0.0001);
  c1->cd(1);
  h2->Draw();
  gr2->Draw("Psame");
  gr3->Draw("Psame");
  gr1->Draw("Psame");//CLAS2010
  gr8->Draw("Psame");//LEPS
  gr_result->Draw("P2same");
  gr2_result->Draw("P2same");
  gr3_result->Draw("P2same");
  gr4->Draw("Lsame");
  gr5->Draw("Lsame");
  gr6->Draw("Lsame");
  gr7->Draw("Lsame");
TBox *b = new TBox(0.,0.28441,16.,0.30141); 
	b->SetFillColor(2); 
	b->SetFillStyle(0);
	b->SetLineColor(2);
	b->SetLineWidth(1);
TBox *b2_1 = new TBox(0.,0.26527,8.,0.30327); 
	b2_1->SetFillColor(3); 
	b2_1->SetFillStyle(0);
	b2_1->SetLineColor(3);
	b2_1->SetLineWidth(1);
TBox *b2_2 = new TBox(8.,0.31396,16.,0.35196); 
	b2_2->SetFillColor(3); 
	b2_2->SetFillStyle(0);
	b2_2->SetLineColor(3);
	b2_2->SetLineWidth(1);
TBox *b3_1 = new TBox(0.,0.23115,6.,0.28115); 
	b3_1->SetFillColor(4); 
	b3_1->SetFillStyle(0);
	b3_1->SetLineColor(4);
	b3_1->SetLineWidth(1);
TBox *b3_2 = new TBox(6.,0.31818,10.,0.36818); 
	b3_2->SetFillColor(4); 
	b3_2->SetFillStyle(0);
	b3_2->SetLineColor(4);
	b3_2->SetLineWidth(1);
TBox *b3_3 = new TBox(10.,0.30168,16.,0.35168); 
	b3_3->SetFillColor(4); 
	b3_3->SetFillStyle(0);
	b3_3->SetLineColor(4);
	b3_3->SetLineWidth(1);
  b->Draw();
  b2_1->Draw();
  b2_2->Draw();
  b3_1->Draw();
  b3_2->Draw();
  b3_3->Draw();
  ax->Draw();
  ay->Draw();
  leg->Draw();

//W dependence
//  TCanvas *c2 = new TCanvas("c2","c2",1100,800);
//  c2->Divide(1,1,0.0001,0.0001);
//  c2->cd(1);
//  gr1_1->Draw("AP");



//for(int i=1;i<200;i++){
//dataset1->SetSelection(i);
//cout<<i<<": "<<dataset1->GetKinematics()->GetW()<<endl;
//dataset2->SetSelection(i);
//cout<<i<<": "<<dataset2->GetKinematics()->GetW()<<endl;
//dataset3->SetSelection(i);
//cout<<i<<": "<<dataset3->GetKinematics()->GetW()<<endl;
//}
c1->Print("Elementary_gpKL.pdf");
}
