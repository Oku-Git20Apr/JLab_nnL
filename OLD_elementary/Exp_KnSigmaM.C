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

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void leps_2006(){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso6.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("LEPS_2006__Kohri__PRL_97_082003__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<"  "<<dataset->GetEntries()<<endl;
//  dataset->ViewSelections();

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,10);

  TGraph *gr[18];
  TLatex *t1[18], *t2[18];
  for(int i=0;i<18;i++){
    dataset->SetSelection(i+5);
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
  for(int i=0;i<18;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same");}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void clas_2010(){
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TFile *ifp;
  ifp = new TFile(Form("%sdata_iso6.root",dir.c_str()));
  TDataset *dataset;
  ifp->GetObject("CLAS_2010__Anefalos__arXiv-0912-4833__dcs",dataset);
  cout<<"DATA NAME = "<<dataset->GetName()<<"  "<<dataset->GetEntries()<<endl;
//  dataset->ViewSelections();

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  SetTH2(h2,"","","",0.020,10);

  TGraph *gr[25];
  TLatex *t1[25], *t2[25];
  for(int i=0;i<25;i++){
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
  for(int i=0;i<25;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
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
  ifp = new TFile(Form("%sdata_iso6.root",dir.c_str()));

  string ofname = "Exp_KnSigmaM.pdf";
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
  ifp->GetObject("LEPS_2006__Kohri__PRL_97_082003__dcs",dataset);
  for(int i=0;i<18;i++){
    dataset->SetSelection(i+5);
    gr[i] = dataset->MakeGraph("costhkcm");
    SetGr(gr[i],"","","",1,1,20,1.0);
    t1[i] = new TLatex(-0.9,0.45,Form("w = %.3lf GeV",dataset->GetKinematics()->GetW()/1000.));
    t2[i] = new TLatex(-0.9,0.40,Form("E_{#gamma} = %.3lf GeV",dataset->GetKinematics()->GetWlab()/1000.));
    t1[i]->SetTextSize(0.08); t1[i]->SetTextColor(4);
    t2[i]->SetTextSize(0.08); t2[i]->SetTextColor(4);
  }
  title = p1->AddText("LEPS 2006"); p1->Draw();
  c1->Print(Form("%s",ofname.c_str()));  p1->Clear();  c1->Clear();
  TCanvas *c1 = new TCanvas("c1","c1",1800,950);
  c1->Divide(5,4);
  for(int i=0;i<18;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t1[i]->Draw("same"); t2[i]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();

/////////
  ifp->GetObject("CLAS_2010__Anefalos__arXiv-0912-4833__dcs",dataset);
  for(int i=0;i<25;i++){
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
  for(int i=0;i<20;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame");    t1[i]->Draw("same");    t2[i]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();
  c1->Divide(5,4);
  for(int i=0;i<5 ;i++){ c1->cd(i+1); h2->Draw(); gr[i+20]->Draw("Psame"); t1[i+20]->Draw("same"); t2[i+20]->Draw("same"); }
  c1->Print(Form("%s",ofname.c_str()));  c1->Clear();

  c1->Print(Form("%s]",ofname.c_str()));
  gSystem->Exit(1);

}

////////////////////////////////////////////////////////////////////////////////////////////
void draw_KnSigmaM(){
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
  ifp = new TFile(Form("%sdata_iso6.root",dir.c_str()));
  TDataset *dataset1;
  ifp->GetObject("CLAS_2010__Anefalos__arXiv-0912-4833__dcs",dataset1);

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.5);
  h2->SetStats(kFALSE);
//  SetTH2_forDraw(h2,"","cos(#theta_{K^{+}}^{c.m.})","d#sigma/d#Omega_{K}^{c.m.} [#mub/sr]");
  SetTH2_forDraw(h2,"","","");
  h2->GetXaxis()->SetLabelOffset(1.);
  h2->GetXaxis()->SetTickLength(0);
  h2->GetYaxis()->SetLabelOffset(1.);
  h2->GetYaxis()->SetLabelSize(0.);
  h2->GetYaxis()->SetTickLength(0);

// CLAS 2010
  TGraphErrors *gr1;
  dataset1->SetSelection(3);
  gr1 = dataset1->MakeGraph("costhkcm");
  double *x1 = gr1->GetX(), *y1 = gr1->GetY(), ey1[100];
    for(int i=0;i<gr1->GetN();i++){ x1[i] *= -1; ey1[i] = gr1->GetErrorY(i); }
  gr1->Clear();
  gr1 = new TGraphErrors(gr1->GetN(),x1,y1,0,ey1);
  SetGr(gr1,"","","",2,2,20,1.5);

// KMaid W=1.835 GeV
  double x2[37] = {
      0.,  5., 10., 15., 20., 25., 30., 35., 40., 45.,
     50., 55., 60., 65., 70., 75., 80., 85., 90., 95.,
    100.,105.,110.,115.,120.,125.,130.,135.,140.,145.,
    150.,155.,160.,165.,170.,175.,180. };
  double y2[37] = {
    0.1668E+03,0.1682E+03,0.1721E+03,0.1783E+03,0.1862E+03,0.1951E+03,0.2042E+03,0.2130E+03,0.2206E+03,0.2266E+03,
    0.2306E+03,0.2324E+03,0.2318E+03,0.2289E+03,0.2238E+03,0.2168E+03,0.2082E+03,0.1981E+03,0.1870E+03,0.1750E+03,
    0.1625E+03,0.1497E+03,0.1369E+03,0.1243E+03,0.1120E+03,0.1003E+03,0.8931E+02,0.7915E+02,0.6993E+02,0.6174E+02,
    0.5462E+02,0.4861E+02,0.4372E+02,0.3994E+02,0.3726E+02,0.3565E+02,0.3512E+02 };
  for(int i=0;i<37;i++){ x2[i] = -cos(x2[i]/180.*PI); y2[i] /= 1000.; }
  TGraphErrors *gr2 = new TGraphErrors(37,x2,y2,0,0);
  SetGr(gr2,"","","",1,1,1,1.0);
  gr2->SetLineWidth(2);
  gr2->SetLineStyle(1);


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

  TLegend *leg = new TLegend(0.40,0.75,0.94,0.94,"","NDC");
  leg->SetFillColor(0);
  leg->SetBorderSize(3);
  leg->SetLineWidth(1);
  leg->SetTextSize(0.050);
  leg->SetMargin(0.15);
  leg->SetNColumns(1);
  leg->AddEntry(gr1     ,Form("CLAS #scale[0.6]{PLB688(2010) 025202.} #scale[0.5]{(#sqrt{s}=%.3lf GeV)}"
                        ,dataset1->GetKinematics()->GetW()/1000.),"pl");
  leg->AddEntry(gr2     ,"KMaid #scale[0.6]{PRC61(1999) 012201.} #scale[0.5]{(#sqrt{s}=1.849 GeV)}","pl");

  TCanvas *c1 = new TCanvas("c1","c1",1100,800);
  c1->Divide(1,1,0.0001,0.0001);
  c1->cd(1);
  h2->Draw();
  gr1->Draw("Psame");
  gr2->Draw("Lsame");
  ax->Draw();
  ay->Draw();
  leg->Draw();
}
