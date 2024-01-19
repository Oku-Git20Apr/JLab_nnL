//-- Mixed Event Analysis  --//
//Systematic Error
//
//Appendix
//
//K. Okuyama (Dec. 29, 2020)
//
//This is taken over from MEA.C
//No array branch mode 
//
void SetTH1(TH1 *h, TString name, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h->SetTitle(name);
  h->SetLineColor(LColor);
  h->SetLineWidth(1);
  h->SetFillStyle(FStyle);
  h->SetFillColor(FColor);

  h->SetTitleFont(42,"");
  h->SetTitleSize(0.04,"");

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetNdivisions(505);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}
void SetTH2(TH2 *h, TString name, TString xname, TString yname, double min=0.8){
  h->SetTitle(name);
  h->SetMinimum(min);
  h->SetLineWidth(0);
  h->SetTitleSize(0.05,"");
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.5);
  h->SetMarkerColor(1);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.40);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}
double F_Voigt( double *x, double *par )
  {
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
    // par[3] : lorentz fwhm
    double val = par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4);
    return val;
  }

double pol2gaus(double *x, double *par) {
   //par[0]=x^2
   //par[1]=x^1
   //par[2]=x^0
   //par[3]=Norm
   //par[4]=Width (sigma) of convoluted Gaussian function
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double mpshift  = -0.22278298;       // Landau maximum location
  double np = 500.0;      // number of convolution steps
  double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, mpc, fland, sum = 0.0, xlow,xupp, step, i;
  double val;

// Range of convolution integral
  xlow = x[0] - sc * par[4];
  xupp = x[0] + sc * par[4];
  step = (xupp-xlow) / np;
  for(i=1.0; i<=np/2; i++) {
     xx = xlow + (i-.5) * step;
     fland = TMath::Gaus(xx,x[0],par[4]);
if(xlow>-0.125&&xupp<0.125){
     sum += fland * (par[0]*x[0]*x[0]+par[1]*x[0]+par[2])*(1./(2.*par[0]*pow(0.125,3.)/3.+2.*par[2]*0.125));
}else sum += fland;

     xx = xupp - (i-.5) * step;
     fland = TMath::Gaus(xx,x[0],par[4]);
if(xlow>-0.125&&xupp<0.125){
     sum += fland * (par[0]*x[0]*x[0]+par[1]*x[0]+par[2])*(1./(2.*par[0]*pow(0.125,3.)/3.+2.*par[2]*0.125));
}else sum += fland;
  }
  val = par[3] * step * sum * invsq2pi / par[4];

  return val;
}

double expgaus2(double *x, double *par, int num) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;
// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[num+2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[num+3];
     fland = TMath::Gaus(xx,x[0],par[num+2]);
     sum += fland * TMath::Exp(-xx/par[num+1]);
     xx = xupp - (i-.5) * step - par[num+3];
     fland = TMath::Gaus(xx,x[0],par[num+2]);
     sum += fland * TMath::Exp(-xx/par[num+1]);
  }
  //val = par[2] * step * sum * invsq2pi / par[3];
  val = par[num] * step * sum * invsq2pi / (par[num+2]*par[num+1]*exp(par[num+3]/par[num+1]));
  return val;
}

double landaugaus2(double *x, double *par, int num) {
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double mpshift  = -0.22278298;       // Landau maximum location
  double np = 500.0;      // number of convolution steps
  double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, mpc, fland, sum = 0.0, xlow,xupp, step, i;
  double val;

// MP shift correction
  mpc = par[num+1] - mpshift * par[num];
// Range of convolution integral
  xlow = x[0] - sc * par[num+3];
  xupp = x[0] + sc * par[num+3];
  step = (xupp-xlow) / np;
// Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
     xx = xlow + (i-.5) * step;
     fland = TMath::Landau(xx,mpc,par[num]) / par[num];
     sum += fland * TMath::Gaus(x[0],xx,par[num+3]);

     xx = xupp - (i-.5) * step;
     fland = TMath::Landau(xx,mpc,par[num]) / par[num];
     sum += fland * TMath::Gaus(x[0],xx,par[num+3]);
  }
  val = par[num+2] * step * sum * invsq2pi / par[num+3];

  return val;
}


double FMM_Lambda_Sigma( double *x, double *par , int num)
  {
  double val = par[num] * TMath::Gaus(x[0],par[num+1],par[num+2]);//Lambda Gaussian
  val += par[num+3] * TMath::Gaus(x[0],par[num+4],par[num+5]);//Sigma Gaussian
  return val;
}

double F_VZ( double *x, double *par )
{
  return pol2gaus(x,par)+par[5]*TMath::Gaus(x[0],par[6],par[7],1)+par[8]*TMath::Gaus(x[0],par[9],par[10],1)+par[11]*TMath::Gaus(x[0],par[12],par[13],1)+par[14]*TMath::Gaus(x[0],par[15],par[16],1);
}

double F_VZ_cell( double *x, double *par)
{
double a1, a2;
double b1, b2;
a1=par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.));
a2=par[0]*par[3]*exp(-0.5*pow((x[0]-par[1]+par[4])/par[5],2.));
b1=par[6]*exp(-0.5*pow((x[0]-par[7])/par[2],2.));
b2=par[6]*par[3]*exp(-0.5*pow((x[0]-par[7]+par[4])/par[5],2.));
double a = a1+a2;
double b = b1+b2;
return a+b;
}

double F_VZ2( double *x, double *par)
{
//par[0]: Al front scale
//par[1]: Al front pos.
//par[2]: Al Gauss sigma
//par[3]: Al second gauss strength
//par[4]: Al(front) second gauss pos. 
//par[5]: Al second gauss sigma 
//par[6]: Al rear scale
//par[7]: Al(rear) second gauss pos. 
//par[8]: pol2 coeff.1
//par[9]: pol2 coeff.2
//par[10]: total scale 

double a1, a2;
double b1, b2;
a1=par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.));
a2=par[0]*par[3]*exp(-0.5*pow((x[0]-par[1]+par[4])/par[5],2.));
b1=par[6]*exp(-0.5*pow((x[0]-par[7])/par[2],2.));
b2=par[6]*par[3]*exp(-0.5*pow((x[0]-par[7]+par[4])/par[5],2.));
double a = a1+a2;
double b = b1+b2;

double c = 0.;
int np = 2000;
for(int i=0;i<np;i++){
 double d = par[8]*pow((x[0]-par[9]),2.)+1.;
 double step = -1.+(double)i/1000.;
 if(step<-0.125) d = 0.;
 if(step> 0.125) d = 0.;

 double aa;
 aa = exp(-0.5*pow((x[0]-step)/par[2],2.))+par[3]*exp(-0.5*pow((x[0]+par[4]-step)/par[5],2.));
 c = c + d*aa;
}
c = par[10]*c;
return a+b+c;
}

double FMM_Response( double *x, double *par ){

   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //par[4]=tau of exp function
   //par[5]=Shift of Function Peak
   //par[6]=Relative Strength
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double mpshift  = -0.22278298;       // Landau maximum location
  double np = 500.0;      // number of convolution steps
  double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, mpc, fland, sum = 0.0, xlow,xupp, step, i;
  double val1, val2;

// MP shift correction
  mpc = par[1] - mpshift * par[0];
// Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow) / np;
// Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
     xx = xlow + (i-.5) * step;
     fland = TMath::Landau(xx,mpc,par[0]) / par[0];
     sum += fland * TMath::Gaus(x[0],xx,par[3]);

     xx = xupp - (i-.5) * step;
     fland = TMath::Landau(xx,mpc,par[0]) / par[0];
     sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  val1 = step * sum * invsq2pi / par[3];

/*------Landau * Gauss convluted------*/

// Range of convolution integral
  sum  = 0.;
  xlow = 0.;
  xupp = x[0] + 1.6 * sc * par[3];
  step = (xupp-xlow) / np;
  if(step<0.)step = 0.;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[5];
     fland = TMath::Gaus(xx,x[0],par[3]);
     sum += fland * TMath::Exp(-xx/par[4]);
     xx = xupp - (i-.5) * step - par[5];
     fland = TMath::Gaus(xx,x[0],par[3]);
     sum += fland * TMath::Exp(-xx/par[4]);
  }
  //val = par[2] * step * sum * invsq2pi / par[3];
  val2 =  step * sum * invsq2pi / (par[3]*par[4]*exp(par[5]/par[4]));
  //val2 =  step * sum * invsq2pi / (par[3]*par[4]);
/*------Exp * Gauss convluted------*/

  return par[2]*(val1+par[6]*val2);//N x (Landau*Gauss) + N' x (Exp*Gauss)

}

double FMM_Res( double *x, double *par ){

	return FMM_Response(x,par)+FMM_Response(x,&par[7]);

}

void empty_cell(){
  
  TChain *chain_dummy = new TChain("T");
	//chain_dummy = (TChain*)file_dummy->Get("T");
  TChain *chain_true = new TChain("T");// = (TChain*)file_true->Get("T");
  //TChain *chain_dummy = new TChain("tree_dummy","");
  //TChain *chain_true = new TChain("tree_true","");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111323.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111324.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111324_1.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111325.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111325_1.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_1.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_2.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_3.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_4.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_5.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_6.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_7.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_8.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_9.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_10.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_11.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_12.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_13.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_14.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_15.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_16.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_17.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_18.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_19.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_20.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_21.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_22.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_23.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_24.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_25.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_26.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_27.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_28.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_29.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_30.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_31.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_32.root");
	chain_dummy->Add("/data/41a/ELS/JLab/E12-17-003/root/tritium_111326_33.root");

	//chain_dummy->SetBranchStatus("*",0);
  	//chain_dummy->SetBranchStatus("L.tr.vz",1);  chain_dummy->SetBranchAddress("L.tr.vz", &L_tr_vz2);
  	//chain_dummy->SetBranchStatus("R.tr.vz",1);  chain_dummy->SetBranchAddress("R.tr.vz", &R_tr_vz2);
	//chain_true->SetBranchStatus("*",0);
  	//chain_true->SetBranchStatus("L.tr.vz",1);  chain_true->SetBranchAddress("L.tr.vz", &L_tr_vz3);
  	//chain_true->SetBranchStatus("R.tr.vz",1);  chain_true->SetBranchAddress("R.tr.vz", &R_tr_vz3);

  TH1F* h_zave_dummy  = new TH1F("h_zave_dummy","h_zave_dummy",1000,-25.,25.);
  TH1F* h_zave_dummy_wHcut  = new TH1F("h_zave_dummy_wHcut","h_zave_dummy |Z|<10",1000,-25.,25.);
  TH2F* h_zz_dummy  = new TH2F("h_zz_dummy","h_zz_dummy",400,-25.,25.,400,-25.,25.);

  int ENum_dummy = chain_dummy->GetEntries();
cout<<"Entries(dummy): "<<ENum_dummy<<endl;
 // for(int i=0;i<ENum_dummy;i++){
//	chain_dummy->GetEntry(i);
	//if(abs(R_tr_vz2-L_tr_vz2)<0.025)h_zave_dummy->Fill((R_tr_vz2+L_tr_vz2)/2.);
	chain_dummy->Project("h_zave_dummy","(R.tr.vz+L.tr.vz)/2.*100.","abs(R.tr.vz-L.tr.vz)<0.025","");
	chain_dummy->Project("h_zz_dummy","R.tr.vz*100.:L.tr.vz*100.","","");
	chain_dummy->Project("h_zave_dummy_wHcut","(R.tr.vz+L.tr.vz)/2.*100.","abs(R.tr.vz-L.tr.vz)<0.025&&abs((R.tr.vz+L.tr.vz)/2.)<0.1","");
//	}

	h_zave_dummy->SetStats(0);
	h_zave_dummy_wHcut->SetStats(1);
	h_zz_dummy->SetStats(1);
	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	h_zave_dummy->Draw();
	TCanvas *c1L = new TCanvas("c1L","c1L",800,600);
	h_zave_dummy->Draw();
	c1L->SetLogy(1);
	TCanvas *c1p = new TCanvas("c1p","c1p",800,800);
	h_zave_dummy_wHcut->Draw();
	TCanvas *c2 = new TCanvas("c2","c2",800,800);
	h_zz_dummy->Draw("colz");

	//c5->SetLeftMargin(0.14);
	//c5->SetRightMargin(0.14);
	//c5->SetTopMargin(0.14);
	//c5->SetBottomMargin(0.14);
	//c5->Modified();
	//c5->Update();
	//gPad->Modified();
	//gPad->Update();

	//c1->Print("./pdf/mea_lll.pdf");
	//c2->Print("./pdf/mea_lcr.pdf");
	//c3->Print("./pdf/mea_rrr.pdf");
	//c4->Print("./pdf/mea_compare.pdf");
	//c5->Print("./pdf/mea_deviation.pdf");
	
	//c1->Print("./pdf/empty_z.pdf");
	c1L->Print("./pdf/empty_z_log.pdf");
	//c1p->Print("./pdf/empty_z_wHcut.pdf");
	//c2->Print("./pdf/empty_zz.pdf");

cout << "Well done!" << endl;
}
