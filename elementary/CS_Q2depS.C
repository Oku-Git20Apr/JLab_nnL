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

void CS_Q2depS(){
  //TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.5);
  TH2F *h2 = new TH2F("h2","h2",1000,0.,5.,1000,0,0.5);
  h2->SetStats(kFALSE);
//  SetTH2_forDraw(h2,"","cos(#theta_{K^{+}}^{c.m.})","d#sigma/d#Omega_{K}^{c.m.} [#mub/sr]");
  SetTH2_forDraw(h2,"","","");
  h2->GetXaxis()->SetLabelOffset(1.);
  h2->GetXaxis()->SetTickLength(0);
  h2->GetYaxis()->SetLabelOffset(1.);
  h2->GetYaxis()->SetLabelSize(0.);
  h2->GetYaxis()->SetTickLength(0);

  //TF1 *fx = new TF1("fx","-x",-1,1);
  TF1 *fx = new TF1("fx","x",0.,5.);
  TGaxis *ax = new TGaxis(0.,0, 5.,0, "fx", 505 );
  //TGaxis *ax = new TGaxis(-1,0, 1,0, "fx", 505 );
  ax->SetLabelFont(42); ax->SetTitleFont(42); ax->CenterTitle();
  ax->SetLabelOffset(0.01);
  ax->SetTitleSize(0.05);  ax->SetTitleOffset(0.90);
  ax->SetMaxDigits(3);
  //ax->SetTitle("cos(#theta_{K}^{c.m.})");
  ax->SetTitle("Q^{2} [(GeV/c)^{2}]");

  TGaxis *ay = new TGaxis(0,0, 0,0.5, 0,0.5, 505 );
  //TGaxis *ay = new TGaxis(-1,0, -1,0.5, 0,0.5, 505 );
  ay->SetDecimals(2);
  ay->SetLabelFont(42); ay->SetTitleFont(42); ay->CenterTitle();
  ay->SetLabelOffset(0.01);
  ay->SetTitleSize(0.05);  ay->SetTitleOffset(0.90);
  ay->SetTitle("d#sigma/d#Omega_{K}^{CM} [#mub/sr]");

  TGraphErrors *Ref = new TGraphErrors(4);
  Ref->SetPoint(0,0.52,0.091);
  Ref->SetPointError(0,0.,0.004);
  Ref->SetPoint(1,0.75,0.067);
  Ref->SetPointError(1,0.,0.003);
  Ref->SetPoint(2,1.,0.045);
  Ref->SetPointError(2,0.,0.003);
  Ref->SetPoint(3,2.,0.019);
  Ref->SetPointError(3,0.,0.002);
  TGraphErrors *Ref28 = new TGraphErrors(5);
  Ref28->SetPoint(0,1.18,0.074);
  Ref28->SetPointError(0,0.,0.030);
  Ref28->SetPoint(1,1.98,0.017);
  Ref28->SetPointError(1,0.,0.012);
  Ref28->SetPoint(2,3.98,0.010);
  Ref28->SetPointError(2,0.,0.010);
  Ref28->SetPoint(3,1.36,0.035);
  Ref28->SetPointError(3,0.,0.036);
  Ref28->SetPoint(4,3.46,0.006);
  Ref28->SetPointError(4,0.,0.010);
  TGraphErrors *Ref29 = new TGraphErrors(6);
  //Ref29->SetPoint(0,0.18,0.539);
  //Ref29->SetPointError(0,0.,0.089);
  Ref29->SetPoint(0,0.29,0.95);
  Ref29->SetPointError(0,0.,0.049);
  Ref29->SetPoint(1,0.40,0.060);
  Ref29->SetPointError(1,0.,0.040);
  Ref29->SetPoint(2,0.76,0.076);
  Ref29->SetPointError(2,0.,0.033);
  //Ref29->SetPoint(4,1.17,0.280);
  //Ref29->SetPointError(4,0.,0.048);
  Ref29->SetPoint(3,0.29,0.330);
  Ref29->SetPointError(3,0.,0.145);
  Ref29->SetPoint(4,0.29,0.189);
  Ref29->SetPointError(4,0.,0.083);
  Ref29->SetPoint(5,0.29,0.144);
  Ref29->SetPointError(5,0.,0.034);
  TGraphErrors *Ref30 = new TGraphErrors(4);
  Ref30->SetPoint(0,0.62,0.046);
  Ref30->SetPointError(0,0.,0.019);
  Ref30->SetPoint(1,1.20,0.049);
  Ref30->SetPointError(1,0.,0.010);
  Ref30->SetPoint(2,2.00,0.008);
  Ref30->SetPointError(2,0.,0.008);
  Ref30->SetPoint(3,1.18,0.036);
  Ref30->SetPointError(3,0.,0.012);
  Ref->SetMarkerStyle(20);
  Ref->SetMarkerSize(2);
  gStyle->SetEndErrorSize(10);
  Ref28->SetMarkerStyle(22);
  Ref28->SetMarkerColor(30);//Dark green
  Ref28->SetLineColor(30);
  Ref28->SetMarkerSize(2);
  Ref29->SetMarkerStyle(20);
  Ref29->SetMarkerColor(kGray);
  Ref29->SetLineColor(kGray);
  Ref29->SetMarkerSize(2);
  Ref30->SetMarkerStyle(22);
  Ref30->SetMarkerColor(9);//Dark purple 
  Ref30->SetLineColor(9);
  Ref30->SetMarkerSize(2);
  //Ref28->SetMarkerStyle(26);
  //Ref28->SetMarkerSize(2);
  //Ref29->SetMarkerStyle(24);
  //Ref29->SetMarkerSize(2);
  //Ref30->SetMarkerStyle(25);
  //Ref30->SetMarkerSize(2);
  TLegend *leg = new TLegend(0.40,0.60,0.90,0.90,"","NDC");
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(3);
  leg->SetLineWidth(1);
  leg->SetTextSize(0.040);
  leg->SetNColumns(1);
  leg->SetMargin(0.15);
  leg->AddEntry(Ref     ,"JLab E93-018 #scale[0.6]{PRC67 (2003) 055205.}","pl");
  leg->AddEntry(Ref28     ,"Cornell #scale[0.6]{PRD15 (1977) 594.}","pl");
  leg->AddEntry(Ref29     ,"Cambridge #scale[0.6]{PRL28 (1972) 1086.}","pl");
  leg->AddEntry(Ref30     ,"Harvard-Cornell #scale[0.6]{PRL32 (1974) 21.}","pl");
  //leg->AddEntry(gr3     ,Form("SAPHIR #scale[0.6]{EPJA19(2004) 251.}      #scale[0.5]{(#sqrt{s}=%.3lf GeV)}"
  //                      ,dataset3->GetKinematics()->GetW()/1000.),"pl");
  //leg->AddEntry(gr2     ,Form("CLAS  #scale[0.6]{PRC73(2006) 035202.}      #scale[0.5]{(#sqrt{s}=%.3lf GeV)}"
  //                      ,dataset2->GetKinematics()->GetW()/1000.),"pl");
  //leg->AddEntry(gr1     ,Form("CLAS  #scale[0.6]{PRC81(2010) 025201.}      #scale[0.5]{(#sqrt{s}=%.3lf GeV)}"
  //                      ,dataset1->GetKinematics()->GetW()/1000.),"pl");
  //leg->AddEntry(gr8     ,"LEPS #scale[0.6]{PRC73(2006) 035214}        #scale[0.5]{(#sqrt{s}=2.163 GeV)}","pl");
  //leg->AddEntry(gr5     ,"SLA    #scale[0.6]{PRC58(1998) 75.}           #scale[0.5]{(#sqrt{s}=2.110 GeV)}","pl");
  //leg->AddEntry(gr4     ,"MAID #scale[0.6]{PRC61(1999) 012201.}       #scale[0.5]{(#sqrt{s}=2.110 GeV)}","pl");
  //leg->AddEntry(gr6     ,"RPR3  #scale[0.6]{PRC73(2006) 045207.}      #scale[0.5]{(#sqrt{s}=2.130 GeV)}","pl");
  //leg->AddEntry(gr7     ,"RPR2011 #scale[0.6]{PRL108(2012) 182002.}  #scale[0.5]{(#sqrt{s}=2.130 GeV)}","pl");


	double x_result[1], y_result[1], xe_result[1], ye_result[1];
	double x2_result[2], y2_result[2], xe2_result[2], ye2_result[2];
	double x3_result[3], y3_result[3], xe3_result[3], ye3_result[3];
	x_result[0]=0.5;//deg
	y_result[0]=0.1032;//ub/sr
	xe_result[0]=0.5;
	ye_result[0]=0.0085;
	x2_result[0]=0.25;//deg
	y2_result[0]=0.1588;//ub/sr
	xe2_result[0]=0.25;
	ye2_result[0]=0.014;
	x2_result[1]=0.75;//deg
	y2_result[1]=0.0529;//ub/sr
	xe2_result[1]=0.25;
	ye2_result[1]=0.0055;
	x3_result[0]=0.225;//deg
	y3_result[0]=0.2291;//ub/sr
	xe3_result[0]=0.225;
	ye3_result[0]=0.021;
	x3_result[1]=0.5;//deg
	y3_result[1]=0.0783;//ub/sr
	xe3_result[1]=0.05;
	ye3_result[1]=0.008;
	x3_result[2]=0.775;//deg
	y3_result[2]=0.03304;//ub/sr
	xe3_result[2]=0.225;
	ye3_result[2]=0.005;
  TGraphErrors *gr_result = new TGraphErrors(1, x_result, y_result, xe_result, ye_result);
  gr_result->SetMarkerSize(1.5);
  gr_result->SetMarkerStyle(21);
  gr_result->SetFillStyle(3004);
  //gr_result->SetLineColor(2);
  gr_result->SetFillColor(1);
  gr_result->SetMarkerColor(2);
  TGraphErrors *gr2_result = new TGraphErrors(2, x2_result, y2_result, xe2_result, ye2_result);
  gr2_result->SetMarkerSize(1.5);
  gr2_result->SetMarkerStyle(21);
  gr2_result->SetFillStyle(3005);
  //gr2_result->SetLineColor(kViolet);
  gr2_result->SetFillColor(1);
  gr2_result->SetMarkerColor(3);
  TGraphErrors *gr3_result = new TGraphErrors(3, x3_result, y3_result, xe3_result, ye3_result);
  gr3_result->SetMarkerSize(1.5);
  gr3_result->SetMarkerStyle(21);
  gr3_result->SetFillStyle(3004);
  //gr3_result->SetLineColor(kOrange);
  gr3_result->SetFillColor(1);
  gr3_result->SetMarkerColor(4);

  TCanvas *c1 = new TCanvas("c1","c1",1400,950);
  c1->Divide(1,1,0.0001,0.0001);
  c1->cd(1);
  h2->Draw();
  Ref30->Draw("Psame");
  Ref29->Draw("Psame");
  Ref28->Draw("Psame");
  Ref->Draw("Psame");
TBox *b = new TBox(x_result[0]-xe_result[0],y_result[0]-ye_result[0],x_result[0]+xe_result[0],y_result[0]+ye_result[0]); 
	b->SetFillColor(2); 
	b->SetFillStyle(0);
	b->SetLineColor(2);
	b->SetLineWidth(1);
TBox *b2_1 = new TBox(x2_result[0]-xe2_result[0],y2_result[0]-ye2_result[0],x2_result[0]+xe2_result[0],y2_result[0]+ye2_result[0]); 
//TBox *b2_1 = new TBox(0.,0.26527,8.,0.30327); 
	b2_1->SetFillColor(3); 
	b2_1->SetFillStyle(0);
	b2_1->SetLineColor(3);
	b2_1->SetLineWidth(1);
//TBox *b2_2 = new TBox(8.,0.31396,16.,0.35196); 
TBox *b2_2 = new TBox(x2_result[1]-xe2_result[1],y2_result[1]-ye2_result[1],x2_result[1]+xe2_result[1],y2_result[1]+ye2_result[1]); 
	b2_2->SetFillColor(3); 
	b2_2->SetFillStyle(0);
	b2_2->SetLineColor(3);
	b2_2->SetLineWidth(1);
//TBox *b3_1 = new TBox(0.,0.23115,6.,0.28115); 
TBox *b3_1 = new TBox(x3_result[0]-xe3_result[0],y3_result[0]-ye3_result[0],x3_result[0]+xe3_result[0],y3_result[0]+ye3_result[0]); 
	b3_1->SetFillColor(4); 
	b3_1->SetFillStyle(0);
	b3_1->SetLineColor(4);
	b3_1->SetLineWidth(1);
//TBox *b3_2 = new TBox(6.,0.31818,10.,0.36818); 
TBox *b3_2 = new TBox(x3_result[1]-xe3_result[1],y3_result[1]-ye3_result[1],x3_result[1]+xe3_result[1],y3_result[1]+ye3_result[1]); 
	b3_2->SetFillColor(4); 
	b3_2->SetFillStyle(0);
	b3_2->SetLineColor(4);
	b3_2->SetLineWidth(1);
//TBox *b3_3 = new TBox(10.,0.30168,16.,0.35168); 
TBox *b3_3 = new TBox(x3_result[2]-xe3_result[2],y3_result[2]-ye3_result[2],x3_result[2]+xe3_result[2],y3_result[2]+ye3_result[2]); 
	b3_3->SetFillColor(4); 
	b3_3->SetFillStyle(0);
	b3_3->SetLineColor(4);
	b3_3->SetLineWidth(1);
  gr3_result->Draw("P2same");
  gr2_result->Draw("P2same");
  gr_result->Draw("P2same");
  b3_1->Draw();
  b3_2->Draw();
  b3_3->Draw();
  b2_1->Draw();
  b2_2->Draw();
  b->Draw();
  ax->Draw();
  ay->Draw();
  leg->Draw();

	c1->Print("/data/41a/ELS/okuyama/JLab_nnL/okuya_macros/mthesis_Fig/pdf/CS_Q2depS.pdf");
}
