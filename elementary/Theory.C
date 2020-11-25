void Theory(){
  int isospin = 1;
// (1)   ===>   gamma + p   ->  K+ + Lambda
// (2)   ===>   gamma + p   ->  K+ + Sigma0
// (3)   ===>   gamma + p   ->  K0 + Sigma+
// (4)   ===>   gamma + n   ->  K0 + Lambda
// (5)   ===>   gamma + n   ->  K0 + Sigma0
// (6)   ===>   gamma + n   ->  K+ + Sigma-
// (7)   ===>   gamma + p   ->  pi0 + p
// (8)   ===>   gamma + p   ->  pi+ + n
// (9)   ===>   gamma + n   ->  pi0 + n
// (10)  ===>   gamma + n   ->  pi- + p
// (11)  ===>   K- + p  ->  gamma  + Lambda
// (12)  ===>   K- + p  ->  gamma  + Sigma0
  double W[33] = { // invariant energy [MeV]
   1609, 1620, 1640, 1650, 1660, 1670, 1680, 1690, 1700, 1710,
   1720, 1730, 1740, 1750, 1760, 1770, 1780, 1790, 1800, 1820,
   1840, 1870, 1900, 1930, 1960, 1990, 2020, 2050, 2080, 2110,
   2140, 2170, 2200 };
  double Cost[37] = { // cos(theta_CM)
  +1.00000, +0.99619, +0.98481, +0.96593, +0.93969, +0.90631, +0.86603, +0.81915, +0.76604, +0.70711,
  +0.64279, +0.57358, +0.50000, +0.42262, +0.34202, +0.25882, +0.17365, +0.08716, +0.00000, -0.08716,
  -0.17365, -0.25882, -0.34202, -0.42262, -0.50000, -0.57358, -0.64279, -0.70711, -0.76604, -0.81915,
  -0.86603, -0.90631, -0.93969, -0.96593, -0.98481, -0.99619, -1.00000 };

  int minWbin = 0;
  if(isospin==1){ minWbin = 1; W[0] = 1609.; }
  if(isospin==2){ minWbin = 8; W[7] = 1686.; }
  if(isospin==6){ minWbin = 8; W[7] = 1691.; }

  TString strangecalcPath = "/home/wasou/binary/strangecalc/";

// published models
  TMultiModel model("model","Kproduction");
  model.SetStrangeModel(TMultiModel::kRPR2011); // iso=1+4, arXiv:1111.6511
//  model.SetStrangeModel(TMultiModel::kRPR2007); // iso=1-6, T.Corthals Phd thesis
  
//  TStrangeModel model("model","Kproduction");
//  model.SetStrangeModel(strangecalcPath + "share/strangecalc-wrapper/models/rpr-2011/iso1+4/init/fit_specification");
//  model.SetStrangeModel("/data1/Calc/RPR/strangecalc/strangecalc-models/trunk/iso1/rpr2011/init/fit_specification");
//  model.SetStrangeModel("/data1/Calc/RPR/strangecalc/strangecalc-models/trunk/iso2/rpr3-full/init/fit_specification");
//  model.SetStrangeModel("/data1/Calc/RPR/strangecalc/strangecalc-models/trunk/iso2/rpr4-full/init/fit_specification");
//  model.SetStrangeModel("/home/wasou/binary/strangecalc/share/strangecalc-wrapper/models/rpr-2007/iso2+6/init/fit_specification");

  TKinematics kin("kin","example kinematics",isospin,"costhkcm:w",0.,1800);
  kin.IsFixed();
  kin.SetVarRange(1,1.,-1.,37);
  kin.IsFixed();
  TCalcInfo crossSection(TCalcInfo::kPhoto,isospin,"dcs");
  double dcs[33][37];
  for(int w=0;w<33;w++){ for(int i=0;i<37;i++) dcs[w][i] = 0.; }

  for(int w=minWbin;w<33;w++){
    kin.SetVar(2,W[w]);
    for(int i=0;i<kin.GetNumberOfSteps();i++){
//    for(int i=0;i<37;i++){
//      kin.SetVar(1,Cost[i]);
      kin.GoTo(i);
      dcs[w][i] = model.GetCalcpoint(kin,&crossSection);
      Cost[i] = kin.GetCosthkcm();
//      cout<<Form("%.3lf  %+.5lf  %.5lf",kin.GetW()/1000.,kin.GetCosthkcm(),dcs)<<endl;
    }
  }

  for(int w=0;w<33;w++){
    cout<<Form("  // W= %.3lf GeV",W[w]/1000.)<<endl;
    cout<<"  {"<<flush;
    for(int i=0;i<10;i++){
      cout<<Form(" %.5lf,",dcs[w][i])<<flush;
    }
    cout<<endl<<"   "<<flush;
    for(int i=10;i<20;i++){
      cout<<Form(" %.5lf,",dcs[w][i])<<flush;
    }
    cout<<endl<<"   "<<flush;
    for(int i=20;i<30;i++){
      cout<<Form(" %.5lf,",dcs[w][i])<<flush;
    }
    cout<<endl<<"   "<<flush;
    for(int i=30;i<36;i++){
      cout<<Form(" %.5lf,",dcs[w][i])<<flush;
    }
    cout<<Form(" %.5lf },",dcs[w][36])<<endl;
  }

  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

  TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.4);
  h2->SetStats(kFALSE);
  h2->SetTitle("");
  h2->GetXaxis()->SetLabelOffset(0.020);
  h2->GetXaxis()->SetLabelSize(10);
  h2->GetXaxis()->SetTitleSize(10);
  h2->GetYaxis()->SetLabelSize(10);
  h2->GetYaxis()->SetTitleSize(10);

  TGraph *gr[33];
  TLatex *t[33];
  for(int i=0;i<33;i++){
    gr[i] = new TGraph(37,Cost,dcs[i]);
    gr[i]->SetLineColor(1);
    gr[i]->SetMarkerStyle(20);
    gr[i]->SetMarkerColor(1);
    gr[i]->SetMarkerSize(1);
    t[i] = new TLatex(-0.9,0.35,Form("w = %.3lf GeV",W[i]/1000.));
    t[i]->SetTextSize(0.08);
  }
  TCanvas *c1 = new TCanvas("c1","c1",1800,950);
  c1->Divide(10,5);
  for(int i=0;i<33;i++){ c1->cd(i+1); h2->Draw(); gr[i]->Draw("Psame"); t[i]->Draw("same"); }


//  ofstream ofp.open("output.dat");
//  output.close();
}
