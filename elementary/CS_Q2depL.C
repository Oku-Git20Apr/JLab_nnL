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

void CS_Q2depL(){


//%%%%%%%%%%%%%%%%%%//
//Electroproduction-
//Theory
//%%%%%%%%%%%%%%%%%%//
//
// KMaid W=2.130 GeV
  double x4[45] = {
  0.00,
  0.05,
  0.10,
  0.15,
  0.20,
  0.25,
  0.30,
  0.35,
  0.40,
  0.45,
  0.50,
  0.55,
  0.60,
  0.65,
  0.70,
  0.75,
  0.80,
  0.85,
  0.90,
  0.95,
  1.00,
  1.05,
  1.10,
  1.15,
  1.20,
  1.25,
  1.30,
  1.35,
  1.40,
  1.45,
  1.50,
  1.55,
  1.60,
  1.65,
  1.70,
  1.75,
  1.80,
  1.85,
  1.90,
  1.95,
  2.00,
  2.05,
  2.10,
  2.15,
  2.20};
  double y4[45] = {
0.0440,
0.4284,
0.4920,
0.4780,
0.4460,
0.4136,
0.3862,
0.3648,
0.3466,
0.3316,
0.3208,
0.3118,
0.3044,
0.2986,
0.2928,
0.2878,
0.2836,
0.2794,
0.2752,
0.2700,
0.2666,
0.2624,
0.2580,
0.2536,
0.2492,
0.2448,
0.2394,
0.2350,
0.2296,
0.2250,
0.2204,
0.2150,
0.2094,
0.2056,
0.2000,
0.1952,
0.1896,
0.1848,
0.1808,
0.1750,
0.1710,
0.1660,
0.1620,
0.1570,
0.1520};
  //for(int i=0;i<37;i++){ x4[i] = -cos(x4[i]/180.*PI); y4[i] /= 1000.; }
  TGraphErrors *gr4 = new TGraphErrors(45,x4,y4,0,0);
  SetGr(gr4,"","","",3,3,1,1.0);//Green
  gr4->SetLineWidth(2);
  gr4->SetLineStyle(1);


// BS1 W=2.130 GeV
  double x9[51] = {
  0.00,
  0.10,
  0.20,
  0.30,
  0.40,
  0.50,
  0.60,
  0.70,
  0.80,
  0.90,
  1.00,
  1.10,
  1.20,
  1.30,
  1.40,
  1.50,
  1.60,
  1.70,
  1.80,
  1.90,
  2.00,
  2.10,
  2.20,
  2.30,
  2.40,
  2.50,
  2.60,
  2.70,
  2.80,
  2.90,
  3.00,
  3.10,
  3.20,
  3.30,
  3.40,
  3.50,
  3.60,
  3.70,
  3.80,
  3.90,
  4.00,
  4.10,
  4.20,
  4.30,
  4.40,
  4.50,
  4.60,
  4.70,
  4.80,
  4.90,
  5.00};
  double y9[51] = {
3.72E+02,
4.30E+02,
4.46E+02,
4.42E+02,
4.27E+02,
4.08E+02,
3.87E+02,
3.65E+02,
3.44E+02,
3.23E+02,
3.03E+02,
2.85E+02,
2.68E+02,
2.52E+02,
2.37E+02,
2.23E+02,
2.10E+02,
1.98E+02,
1.87E+02,
1.77E+02,
1.67E+02,
1.59E+02,
1.51E+02,
1.43E+02,
1.36E+02,
1.29E+02,
1.23E+02,
1.18E+02,
1.12E+02,
1.07E+02,
1.03E+02,
9.81E+01,
9.40E+01,
9.01E+01,
8.64E+01,
8.30E+01,
7.97E+01,
7.67E+01,
7.38E+01,
7.10E+01,
6.84E+01,
6.59E+01,
6.36E+01,
6.14E+01,
5.93E+01,
5.72E+01,
5.53E+01,
5.35E+01,
5.18E+01,
5.01E+01,
4.85E+01};
  for(int i=0;i<51;i++){y9[i] /= 1000.; }//[nb/sr]-->[ub/sr]
  TGraphErrors *gr9 = new TGraphErrors(51,x9,y9,0,0);
  SetGr(gr9,"","","",6,6,1,1.0);//Magenta
  gr9->SetLineWidth(2);
  gr9->SetLineStyle(1);

// BS3 W=2.130 GeV
  double x10[51] = {
  0.00,
  0.10,
  0.20,
  0.30,
  0.40,
  0.50,
  0.60,
  0.70,
  0.80,
  0.90,
  1.00,
  1.10,
  1.20,
  1.30,
  1.40,
  1.50,
  1.60,
  1.70,
  1.80,
  1.90,
  2.00,
  2.10,
  2.20,
  2.30,
  2.40,
  2.50,
  2.60,
  2.70,
  2.80,
  2.90,
  3.00,
  3.10,
  3.20,
  3.30,
  3.40,
  3.50,
  3.60,
  3.70,
  3.80,
  3.90,
  4.00,
  4.10,
  4.20,
  4.30,
  4.40,
  4.50,
  4.60,
  4.70,
  4.80,
  4.90,
  5.00};
  double y10[51] = {
3.73E+02,
4.60E+02,
4.53E+02,
4.21E+02,
3.86E+02,
3.54E+02,
3.25E+02,
3.00E+02,
2.79E+02,
2.60E+02,
2.44E+02,
2.30E+02,
2.17E+02,
2.05E+02,
1.95E+02,
1.86E+02,
1.77E+02,
1.69E+02,
1.62E+02,
1.55E+02,
1.48E+02,
1.42E+02,
1.37E+02,
1.32E+02,
1.26E+02,
1.22E+02,
1.17E+02,
1.13E+02,
1.09E+02,
1.05E+02,
1.01E+02,
9.76E+01,
9.41E+01,
9.09E+01,
8.77E+01,
8.47E+01,
8.18E+01,
7.90E+01,
7.64E+01,
7.38E+01,
7.13E+01,
6.90E+01,
6.67E+01,
6.45E+01,
6.24E+01,
6.03E+01,
5.84E+01,
5.65E+01,
5.46E+01,
5.29E+01,
5.12E+01};
  for(int i=0;i<51;i++){y10[i] /= 1000.; }//[nb/sr]-->[ub/sr]
  TGraphErrors *gr10 = new TGraphErrors(51,x10,y10,0,0);
  SetGr(gr10,"","","",7,7,1,1.0);//Cyan
  gr10->SetLineWidth(2);
  gr10->SetLineStyle(1);

  //TH2F *h2 = new TH2F("h2","h2",1000,-1,1,1000,0,0.5);
  TH2F *h2 = new TH2F("h2","h2",1000,0.,5.,1000,0,1.0);
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

  TGaxis *ay = new TGaxis(0,0, 0,1.0, 0,1.0, 505 );
  //((TGaxis*)h2->GetYaxis())->SetMaxDigits(2);
  ay->SetDecimals(2);
  //TGaxis *ay = new TGaxis(-1,0, -1,0.5, 0,0.5, 505 );
  ay->SetLabelFont(42); ay->SetTitleFont(42); ay->CenterTitle();
  ay->SetLabelOffset(0.01);
  ay->SetTitleSize(0.05);  ay->SetTitleOffset(0.90);
  ay->SetTitle("d#sigma/d#Omega_{K}^{c.m.} [#mub/sr]");

  TGraphErrors *Ref = new TGraphErrors(4);
  Ref->SetPoint(0,0.52,0.378);
  Ref->SetPointError(0,0.,0.013);
  Ref->SetPoint(1,0.75,0.355);
  Ref->SetPointError(1,0.,0.011);
  Ref->SetPoint(2,1.,0.317);
  Ref->SetPointError(2,0.,0.012);
  Ref->SetPoint(3,2.,0.189);
  Ref->SetPointError(3,0.,0.006);
  TGraphErrors *Ref28 = new TGraphErrors(5);
  Ref28->SetPoint(0,1.18,0.253);
  Ref28->SetPointError(0,0.,0.046);
  Ref28->SetPoint(1,1.98,0.193);
  Ref28->SetPointError(1,0.,0.023);
  Ref28->SetPoint(2,3.98,0.079);
  Ref28->SetPointError(2,0.,0.016);
  Ref28->SetPoint(3,1.36,0.216);
  Ref28->SetPointError(3,0.,0.050);
  Ref28->SetPoint(4,3.46,0.120);
  Ref28->SetPointError(4,0.,0.038);
  TGraphErrors *Ref29 = new TGraphErrors(8);
  Ref29->SetPoint(0,0.18,0.539);
  Ref29->SetPointError(0,0.,0.089);
  Ref29->SetPoint(1,0.29,0.504);
  Ref29->SetPointError(1,0.,0.060);
  Ref29->SetPoint(2,0.40,0.402);
  Ref29->SetPointError(2,0.,0.049);
  Ref29->SetPoint(3,0.76,0.290);
  Ref29->SetPointError(3,0.,0.045);
  Ref29->SetPoint(4,1.17,0.280);
  Ref29->SetPointError(4,0.,0.048);
  Ref29->SetPoint(5,0.29,0.676);
  Ref29->SetPointError(5,0.,0.136);
  Ref29->SetPoint(6,0.29,0.360);
  Ref29->SetPointError(6,0.,0.092);
  Ref29->SetPoint(7,0.29,0.425);
  Ref29->SetPointError(7,0.,0.042);
  TGraphErrors *Ref30 = new TGraphErrors(4);
  Ref30->SetPoint(0,0.62,0.365);
  Ref30->SetPointError(0,0.,0.037);
  Ref30->SetPoint(1,1.20,0.229);
  Ref30->SetPointError(1,0.,0.017);
  Ref30->SetPoint(2,2.00,0.185);
  Ref30->SetPointError(2,0.,0.019);
  Ref30->SetPoint(3,1.18,0.279);
  Ref30->SetPointError(3,0.,0.020);
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
  TLegend *leg = new TLegend(0.45,0.60,0.90,0.90,"","NDC");
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
  leg->AddEntry(gr4     ,"KM #scale[0.6]{PRC61(1999) 012201.}","pl");
  leg->AddEntry(gr9     ,"BS1 #scale[0.6]{PRC93(2016) 025204., PRC97(2018) 025202.}","pl");
  leg->AddEntry(gr10     ,"BS3 #scale[0.6]{PRC97(2018) 025202.}","pl");
  //leg->AddEntry(gr3     ,Form("SAPHIR #scale[0.6]{EPJA19(2004) 251.}      #scale[0.5]{(W=%.3lf GeV)}"
  //                      ,dataset3->GetKinematics()->GetW()/1000.),"pl");
  //leg->AddEntry(gr2     ,Form("CLAS  #scale[0.6]{PRC73(2006) 035202.}      #scale[0.5]{(W=%.3lf GeV)}"
  //                      ,dataset2->GetKinematics()->GetW()/1000.),"pl");
  //leg->AddEntry(gr1     ,Form("CLAS  #scale[0.6]{PRC81(2010) 025201.}      #scale[0.5]{(W=%.3lf GeV)}"
  //                      ,dataset1->GetKinematics()->GetW()/1000.),"pl");
  //leg->AddEntry(gr8     ,"LEPS #scale[0.6]{PRC73(2006) 035214}        #scale[0.5]{(W=2.163 GeV)}","pl");
  //leg->AddEntry(gr5     ,"SLA    #scale[0.6]{PRC58(1998) 75.}           #scale[0.5]{(W=2.110 GeV)}","pl");
  //leg->AddEntry(gr4     ,"MAID #scale[0.6]{PRC61(1999) 012201.}       #scale[0.5]{(W=2.110 GeV)}","pl");
  //leg->AddEntry(gr6     ,"RPR3  #scale[0.6]{PRC73(2006) 045207.}      #scale[0.5]{(W=2.130 GeV)}","pl");
  //leg->AddEntry(gr7     ,"RPR2011 #scale[0.6]{PRL108(2012) 182002.}  #scale[0.5]{(W=2.130 GeV)}","pl");


	double x_result[1], y_result[1], xe_result[1], xel_result[1], xeh_result[1], ye_result[1], yel_result[1], yeh_result[1];
	double x2_result[2], y2_result[2], xe2_result[2], xel2_result[2], xeh2_result[2], ye2_result[2], yel2_result[2], yeh2_result[2];
	double x3_result[3], y3_result[3], xe3_result[3], ye3_result[3], yel3_result[3], yeh3_result[3];
	//string result_in = "./result_input.dat";//D.C.S. Result
	string result_in = "./result_input_2022.dat";//D.C.S. Result
	string buf;
	int npoint = 0;
	int npoint2 = 0;
	int npoint3 = 0;
	int flag_LS = -1;//Lambda(1) or Sigma0(2)
	int flag_dep = -1;//Angle(1) or Q2(2)
	int flag_div = -1;//1-Div.(1) or 2-Div.(2) or 3-Div.(3)
	double xval, xerr, yval, yerr, yerr_l, yerr_h;
	double x[0], y[0], xe[0], ye[0], yel[0], yeh[0];
	double x2[1], y2[1], xe2[1], ye2[1], yel2[1], yeh2[1];
	double x3[2], y3[2], xe3[2], ye3[2], yel3[1], yeh3[1];
	ifstream ifp(result_in.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << result_in.c_str() << endl;
	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> flag_LS >> flag_dep >> flag_div >> xval >> yval >> xerr >> yerr >> yerr_l >> yerr_h;
		cout << flag_LS << ", " << flag_dep << ", " << flag_div << ", " << xval << ", " << yval << ", " << xerr << ", " << yerr << ", " << yerr_l << ", " << yerr_h <<endl;
		
		if(flag_LS==1&&flag_dep==2&&flag_div==1){//1-Div.
			x_result[npoint] = xval; 
			y_result[npoint] = yval;
			xe_result[npoint]= xerr; 
			xel_result[npoint]= xerr; 
			xeh_result[npoint]= xerr; 
			ye_result[npoint]= yerr; 
			yel_result[npoint]= yerr_l; 
			yeh_result[npoint]= yerr_h; 
			npoint++;
		}
		if(flag_LS==1&&flag_dep==2&&flag_div==2){//2-Div.
			x2_result[npoint2] = xval; 
			y2_result[npoint2] = yval;
			xe2_result[npoint2]= xerr; 
			xel2_result[npoint2]= xerr; 
			xeh2_result[npoint2]= xerr; 
			ye2_result[npoint2]= yerr; 
			yel2_result[npoint2]= yerr_l; 
			yeh2_result[npoint2]= yerr_h; 
			npoint2++;
		}
		if(flag_LS==1&&flag_dep==2&&flag_div==3){//3-Div.
			x3_result[npoint3] = xval; 
			y3_result[npoint3] = yval;
			xe3_result[npoint3]= xerr; 
			ye3_result[npoint3]= yerr; 
			npoint3++;
		}
	}
	//x_result[0]=0.5;//deg
	//y_result[0]=0.3195;//ub/sr
	//xe_result[0]=0.5;
	//ye_result[0]=0.016;
	//x2_result[0]=0.25;//deg
	//y2_result[0]=0.4244;//ub/sr
	//xe2_result[0]=0.25;
	//ye2_result[0]=0.023;
	//x2_result[1]=0.75;//deg
	//y2_result[1]=0.2341;//ub/sr
	//xe2_result[1]=0.25;
	//ye2_result[1]=0.014;
	//x3_result[0]=0.225;//deg
	//y3_result[0]=0.5514;//ub/sr
	//xe3_result[0]=0.225;
	//ye3_result[0]=0.034;
	//x3_result[1]=0.5;//deg
	//y3_result[1]=0.2745;//ub/sr
	//xe3_result[1]=0.05;
	//ye3_result[1]=0.016;
	//x3_result[2]=0.775;//deg
	//y3_result[2]=0.2142;//ub/sr
	//xe3_result[2]=0.225;
	//ye3_result[2]=0.016;
	//
	//
  //TGraphAsymmErrors *gr_result = new TGraphAsymmErrors(1, x_result, y_result, xel_result, xeh_result, yel_result, yeh_result);
  TGraphErrors *gr_result = new TGraphErrors(1, x_result, y_result, xe_result, ye_result);
  gr_result->SetMarkerSize(1.5);
  gr_result->SetMarkerStyle(21);
  //gr_result->SetFillStyle(3004);
  gr_result->SetLineColor(2);
  gr_result->SetLineWidth(2);
  gr_result->SetFillColor(1);
  gr_result->SetMarkerColor(2);
  //TGraphAsymmErrors *gr2_result = new TGraphAsymmErrors(2, x2_result, y2_result, xel2_result, xeh2_result, yel2_result, yeh2_result);
  TGraphErrors *gr2_result = new TGraphErrors(2, x2_result, y2_result, xe2_result, ye2_result);
  gr2_result->SetMarkerSize(1.5);
  gr2_result->SetMarkerStyle(21);
  //gr2_result->SetFillStyle(3005);
  gr2_result->SetLineColor(kAzure);
  gr2_result->SetLineWidth(2);
  gr2_result->SetFillColor(1);
  gr2_result->SetMarkerColor(4);
	//
//  TGraphErrors *gr_result = new TGraphErrors(1, x_result, y_result, xe_result, ye_result);
//  gr_result->SetMarkerSize(1.5);
//  gr_result->SetMarkerStyle(21);
//  gr_result->SetFillStyle(3004);
//  //gr_result->SetLineColor(2);
//  gr_result->SetFillColor(1);
//  gr_result->SetMarkerColor(2);
//  TGraphErrors *gr2_result = new TGraphErrors(2, x2_result, y2_result, xe2_result, ye2_result);
//  gr2_result->SetMarkerSize(1.5);
//  gr2_result->SetMarkerStyle(21);
//  gr2_result->SetFillStyle(3005);
//  //gr2_result->SetLineColor(kViolet);
//  gr2_result->SetFillColor(1);
//  gr2_result->SetMarkerColor(3);
//  TGraphErrors *gr3_result = new TGraphErrors(3, x3_result, y3_result, xe3_result, ye3_result);
//  gr3_result->SetMarkerSize(1.5);
//  gr3_result->SetMarkerStyle(21);
//  gr3_result->SetFillStyle(3004);
//  //gr3_result->SetLineColor(kOrange);
//  gr3_result->SetFillColor(1);
//  gr3_result->SetMarkerColor(4);

  //TCanvas *c1 = new TCanvas("c1","c1",1400,950);
  TCanvas *c1 = new TCanvas("c1","c1",1300,950);
  c1->Divide(1,1,0.0001,0.0001);
  c1->cd(1);
  h2->Draw();
  gr4->Draw("Csame");//KMaid
  gr9->Draw("Csame");//BS1
  gr10->Draw("Csame");//BS3
  Ref30->Draw("Psame");
  Ref29->Draw("Psame");
  Ref28->Draw("Psame");
  Ref->Draw("Psame");


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
////TBox *b3_1 = new TBox(0.,0.23115,6.,0.28115); 
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
  //b3_1->Draw();
  //b3_2->Draw();
  //b3_3->Draw();
  b2_1->Draw();
  b2_2->Draw();
  //b2_stat1->Draw();
  //b2_stat2->Draw();
  //b->Draw();
  //b_stat->Draw();
  //gr3_result->Draw("P2same");
  gr2_result->Draw("Psame");
  //gr_result->Draw("Psame");
  ax->Draw();
  ay->Draw();
  leg->Draw();

//	c1->Print("/data/41a/ELS/okuyama/JLab_nnL/okuya_macros/mthesis_Fig/pdf/CS_Q2depL2_2022.pdf");
}
