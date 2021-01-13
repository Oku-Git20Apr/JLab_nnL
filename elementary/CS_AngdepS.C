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

void CS_AngdepS(){


//%%%%%%%%%%%%%%%%%%//
//Electroproduction-
//Theory
//%%%%%%%%%%%%%%%%%%//
//
// KMaid W=2.130 GeV
  double x4[61] = {
0,
3,
6,
9,
12,
15,
18,
21,
24,
27,
30,
33,
36,
39,
42,
45,
48,
51,
54,
57,
60,
63,
66,
69,
72,
75,
78,
81,
84,
87,
90,
93,
96,
99,
102,
105,
108,
111,
114,
117,
120,
123,
126,
129,
132,
135,
138,
141,
144,
147,
150,
153,
156,
159,
162,
165,
168,
171,
174,
177,
180};
  double y4[61] = {
0.0892,
0.0892,
0.0930,
0.0968,
0.1034,
0.1100,
0.1176,
0.1252,
0.1326,
0.1400,
0.1456,
0.1520,
0.1566,
0.1600,
0.1626,
0.1660,
0.1666,
0.1680,
0.1686,
0.1682,
0.1696,
0.1692,
0.1688,
0.1684,
0.1680,
0.1676,
0.1682,
0.1678,
0.1684,
0.1680,
0.1686,
0.1692,
0.1698,
0.1704,
0.1710,
0.1714,
0.1720,
0.1734,
0.1738,
0.1752,
0.1764,
0.1766,
0.1778,
0.1788,
0.1806,
0.1824,
0.1840,
0.1856,
0.1870,
0.1902,
0.1924,
0.1946,
0.1968,
0.1998,
0.2020,
0.2044,
0.2070,
0.2088,
0.2098,
0.2104,
0.2112};
  //for(int i=0;i<37;i++){ x4[i] = -cos(x4[i]/180.*PI); y4[i] /= 1000.; }
  TGraphErrors *gr4 = new TGraphErrors(61,x4,y4,0,0);
  SetGr(gr4,"","","",3,3,1,1.0);//Green
  gr4->SetLineWidth(2);
  gr4->SetLineStyle(1);


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

  //TF1 *fx = new TF1("fx","-x",-1,1);
  TF1 *fx = new TF1("fx","x",0.,180.);
  TGaxis *ax = new TGaxis(0.,0, 180.,0, "fx", 505 );
  //TGaxis *ax = new TGaxis(-1,0, 1,0, "fx", 505 );
  ax->SetLabelFont(42); ax->SetTitleFont(42); ax->CenterTitle();
  ax->SetLabelOffset(0.01);
  ax->SetTitleSize(0.05);  ax->SetTitleOffset(0.90);
  ax->SetMaxDigits(3);
  //ax->SetTitle("cos(#theta_{K}^{c.m.})");
  ax->SetTitle("#theta_{#gamma K}^{CM} [degree]");

  TGaxis *ay = new TGaxis(0,0, 0,0.5, 0,0.5, 505 );
  //((TGaxis*)h2->GetYaxis())->SetMaxDigits(2);
  ay->SetDecimals(2);
  //TGaxis *ay = new TGaxis(-1,0, -1,0.5, 0,0.5, 505 );
  ay->SetLabelFont(42); ay->SetTitleFont(42); ay->CenterTitle();
  ay->SetLabelOffset(0.01);
  ay->SetTitleSize(0.05);  ay->SetTitleOffset(0.90);
  ay->SetTitle("d#sigma/d#Omega_{K}^{CM} [#mub/sr]");

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
  TLegend *leg = new TLegend(0.45,0.80,0.90,0.90,"","NDC");
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(3);
  leg->SetLineWidth(1);
  leg->SetTextSize(0.040);
  leg->SetNColumns(1);
  leg->SetMargin(0.15);
  //leg->AddEntry(Ref     ,"JLab E93-018 #scale[0.6]{PRC67 (2003) 055205.}","pl");
  //leg->AddEntry(Ref28     ,"Cornell #scale[0.6]{PRD15 (1977) 594.}","pl");
  //leg->AddEntry(Ref29     ,"Cambridge #scale[0.6]{PRL28 (1972) 1086.}","pl");
  //leg->AddEntry(Ref30     ,"Harvard-Cornell #scale[0.6]{PRL32 (1974) 21.}","pl");
  leg->AddEntry(gr4     ,"KM #scale[0.6]{PRC61(1999) 012201.}","pl");
  //leg->AddEntry(gr9     ,"BS1 #scale[0.6]{PRC93(2016) 025204., PRC97(2018) 025202.}  #scale[0.5]{(W=2.130 GeV)}","pl");
  //leg->AddEntry(gr10     ,"BS3 #scale[0.6]{PRC97(2018) 025202.} #scale[0.5]{(W=2.130 GeV)}","pl");
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
	string result_in = "./result_input.dat";//D.C.S. Result
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
		
		if(flag_LS==2&&flag_dep==1&&flag_div==1){//1-Div.
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
		if(flag_LS==2&&flag_dep==1&&flag_div==2){//2-Div.
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
		if(flag_LS==2&&flag_dep==1&&flag_div==3){//3-Div.
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



  TGraphAsymmErrors *gr_result = new TGraphAsymmErrors(1, x_result, y_result, xel_result, xeh_result, yel_result, yeh_result);
  //TGraphErrors *gr_result = new TGraphErrors(1, x_result, y_result, xe_result, ye_result);
  gr_result->SetMarkerSize(1.5);
  gr_result->SetMarkerStyle(21);
  gr_result->SetFillStyle(3004);
  //gr_result->SetLineColor(2);
  gr_result->SetFillColor(1);
  gr_result->SetMarkerColor(2);
  TGraphAsymmErrors *gr2_result = new TGraphAsymmErrors(2, x2_result, y2_result, xel2_result, xeh2_result, yel2_result, yeh2_result);
  //TGraphErrors *gr2_result = new TGraphErrors(2, x2_result, y2_result, xe2_result, ye2_result);
  gr2_result->SetMarkerSize(1.5);
  gr2_result->SetMarkerStyle(21);
  gr2_result->SetFillStyle(3005);
  //gr2_result->SetLineColor(kViolet);
  gr2_result->SetFillColor(1);
  gr2_result->SetMarkerColor(4);

  TCanvas *c1 = new TCanvas("c1","c1",1400,950);
  c1->Divide(1,1,0.0001,0.0001);
  c1->cd(1);
  h2->Draw();
  gr4->Draw("Csame");//KMaid
  //gr9->Draw("Csame");//BS1
  //gr10->Draw("Csame");//BS3
  //Ref30->Draw("Psame");
  //Ref29->Draw("Psame");
  //Ref28->Draw("Psame");
  //Ref->Draw("Psame");
  //
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
  //
  //
  //
  //
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
  //gr3_result->Draw("P2same");
  gr2_result->Draw("P2same");
  //gr_result->Draw("P2same");
  //b3_1->Draw();
  //b3_2->Draw();
  //b3_3->Draw();
  b2_1->Draw();
  b2_2->Draw();
  b2_stat1->Draw();
  b2_stat2->Draw();
  //b->Draw();
  //b_stat->Draw();
  ax->Draw();
  ay->Draw();
  leg->Draw();

	//c1->Print("/data/41a/ELS/okuyama/JLab_nnL/okuya_macros/mthesis_Fig/pdf/CS_AngdepS2.pdf");
}
