//-----------------------------//
//--  Path Length Calc.      --//
//--	from Cointime        --//
//-----------------------------//

//K. Okuyama (Jul. 18, 2021)
//-- Path Length (from Cointime) vs. Momentum
//K. Okuyama (Jun. 15, 2022)
//-- final result, systematic error
//K. Okuyama (Jan. 30, 2023)
//-- almost same as "kaonSR_2022.C"
//-- added SetTH1, etc.

double PI=4.*atan(1.);

//specified for this macro!
void SetGr(TGraph *gr, TString name, TString xname, TString yname, int LColor=1, int LWidth=1, int LStyle=1, int MColor=1, int MStyle=20, double MSize=1.){
  gr->SetTitle(name);
  gr->SetName(name);

  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetTitleOffset(0.90);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelOffset(0.01);

  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleOffset(1.20);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)gr->GetYaxis())->SetMaxDigits(4);

  gr->SetLineColor(LColor);
  gr->SetLineStyle(LStyle);
  gr->SetLineWidth(LWidth);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
}
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

double F_path(double *x, double *par){//2D func. (mom, t_diff)
	double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
	double Mp = 0.938272046;         // proton       mass (GeV/c2)
	double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
	double mom = x[0]/1000.;//GeV/c
	double tdiff = x[1];//ns
	double betapi = mom/sqrt(Mpi*Mpi+mom*mom);
	double betap  = mom/sqrt(Mp*Mp+mom*mom);
	double pathlen = betapi*betap*LightVelocity*tdiff/(betapi-betap);
	return pathlen;//m
}

double F1_path(double *x, double *par){
	double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
	double Mp = 0.938272046;         // proton       mass (GeV/c2)
	double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
	double mom = x[0]/1000.;//GeV/c
	double b = par[0];//intercept
	double a = par[1];//slope
	double tdiff = a*mom*1000.+b;//pol1
	double betapi = mom/sqrt(Mpi*Mpi+mom*mom);
	double betap  = mom/sqrt(Mp*Mp+mom*mom);
	double pathlen = betapi*betap*LightVelocity*tdiff/(betapi-betap);
	return pathlen;//m
}
double F1_sr(double *x, double *par){
	double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
	double Mp = 0.938272046;         // proton       mass (GeV/c2)
	double MK = 0.493677;            // charged Kaon mass (GeV/c2)
	double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
	double mom = x[0]/1000.;//GeV/c
	double b = par[0];//intercept
	double a = par[1];//slope
	double tdiff = a*mom*1000.+b;//pol1
	double betapi = mom/sqrt(Mpi*Mpi+mom*mom);
	double betap  = mom/sqrt(Mp*Mp+mom*mom);
	double pathlen = betapi*betap*LightVelocity*tdiff/(betapi-betap);
	return exp(-pathlen*MK/mom/3.71);//m
}

void func_shade(TCanvas *cc, TF1 *f1, TF1 *f2, TF1 *f3) {
//shade the area between f1 and f2 and draw f3 on top

//create a TGraph to store the function values
//shaded area is the fill/color/style of f1
TGraph *gr = new TGraph();
SetGr(gr, "", "Momentum [MeV/c]", "#varepsilon(Decay)", kAzure, kRed, kRed);
gr->SetFillColor(f1->GetFillColor());
gr->SetFillStyle(f1->GetFillStyle());
f1->Draw("l");
cc->Update();
Double_t ymin = cc->GetUymin();
f2->Draw("l");
cc->Update();
Double_t ymax = cc->GetUymax();
f3->Draw("l");
cc->Update();
//get picture range
Double_t xmin = cc->GetUxmin();
Double_t xmax = cc->GetUxmax();

//process first function
Int_t npx = f3->GetNpx();
Int_t npoints=0;
Double_t dx = (xmax-xmin)/npx;
Double_t x = xmin+0.5*dx;
while (x <= xmax) {
Double_t y = f1->Eval(x);
if (y < ymin) y = ymin;
if (y > ymax) y = ymax;
gr->SetPoint(npoints,x,y);
npoints++;
x += dx;
}
//process second function
x = xmax-0.5*dx;
while (x >= xmin) {
Double_t y = f2->Eval(x);
if (y < ymin) y = ymin;
if (y > ymax) y = ymax;
gr->SetPoint(npoints,x,y);
npoints++;
x -= dx;
}
gr->GetYaxis()->SetDecimals(3);
gr->GetXaxis()->SetRangeUser(1760.,1900.);
gr->Draw("fa"); //draw graph with fill area option
f3->Draw("lsame"); //superimpose function
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

double Dgauss( double *x, double *par )
  {
    // par[0] : scale1 
    // par[1] : location1
    // par[2] : sigma1
    // par[3] : scale2 
    // par[4] : location2
    // par[5] : sigma2
    double val = par[0] * TMath::Gaus(x[0],par[1],par[2]);
    val += par[3] * TMath::Gaus(x[0],par[4],par[5]);
    return val;
  }


double fcoin_template( double *x, double *par , int shift, int num)
{
  return par[num] * TMath::Gaus(x[0],par[num+1]-2.0*shift,par[num+2]);//Lambda Gaussian
}

double fcoin_acc( double *x, double *par)
{
	double val=0.;
	for(int i=0;i<=80;i++){
	//val += par[num]*TMath::Gaus(x[0],par[num+1]+par[num+12]*i-20,par[num+2]);
    val += par[0] * (1.-par[5]) * TMath::Voigt(x[0]-3.179-par[1]-par[4]*(double)i+30.,par[2],par[3],4);
    //val += par[0] * par[5] * TMath::Voigt(x[0]+7.925-par[1]-par[4]*i+30.,par[6],par[7],4);
    val += par[0] * par[5] * 0.6 * TMath::Gaus(x[0],-8.303+par[1]+par[4]*(double)i-30.,par[6])/(sqrt(2.*PI)*par[6]);
    val += par[0] * par[5] * 0.4 * TMath::Gaus(x[0],-7.463+par[1]+par[4]*(double)i-30.,par[7])/(sqrt(2.*PI)*par[7]);
	}
	return val;
}

//double fcoin_total( double *x, double *par ){
//
//	if(reject_flag&&x[0]>-12.&&x[0]<10.){TF1::RejectPoint();return 0;}
//	//return fcoin_template(x,par,-10,0)+expgaus2(x,par,6)+expgaus2(x,par,10);//+expgaus2(x,par,14);
//	else return fcoin_acc(x,par)+par[8]*TMath::Voigt(x[0]-par[9],par[10],par[11],4)+par[12]*TMath::Voigt(x[0]-par[13],par[14],par[15],4)+par[16]*TMath::Gaus(x[0],par[17],par[18])+par[19]*TMath::Gaus(x[0],par[20],par[21]);
//
//}

void kaonSR(){
//ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);
	string pdfname = "temp.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  //TFile *file = new TFile("../h2all5.root","read");//input file of all H2 run(default: h2all4.root)
  TFile *file = new TFile("../h2all_2020Nov.root","read");//input file of all H2 run(default: h2all4.root)
  TTree *tree = (TTree*)file->Get("tree_out");

    
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);



//---Physics Constant---//
 
 const double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
 const double Mp = 0.938272046;         // proton       mass (GeV/c2)
 const double MK = 0.493677;            // charged Kaon mass (GeV/c2)
 const double Me = 0.510998928e-3;      // electron     mass (GeV/c2)
 const double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
 const double PI=3.14159265359;
 const double ML = 1.115683;            // Lambda       mass (GeV/c2)
 const double MS0 = 1.192642;           // Sigma Zero   mass (GeV/c2)
 const double def_n_L=250.; 
 const double def_sig_L=0.003; 
 const double def_mean_L=0.0;
 const double def_n_S=70.; 
 const double def_sig_S=0.004;
 const double def_mean_S=MS0-ML;
 const double def_sig_p=0.852;
 const double def_mean_p=-8.0;
 const double def_sig_pi=0.443;
 const double def_mean_pi=3.0;
 const double def_sig_k=0.644;
 const double def_mean_k=0.0;
 const double def_acc=27.7;
 const double min_coin_c=-20.0;
 const double max_coin_c=20.0;
 const double tdc_time=0.056;//ns
 int bin_coin_c=(int)((max_coin_c-min_coin_c)/tdc_time);
 const double min_mm=-0.1;//GeV/c^2
 const double max_mm=0.2;//GeV/c^2
 int bin_mm=(max_mm-min_mm)/0.002; //Counts/2 MeV
 bin_mm=(int)bin_mm;


//---For Efficiency Result---//
 double center_pi, center_k, center_p, center_L, center_S;
 double range_pi, range_k, range_p, range_L, range_S;
 double n_pi_nocut, n_k_nocut, n_p_nocut, n_L_nocut, n_S_nocut;
 double const_pi_nocut, const_k_nocut, const_p_nocut, const_L_nocut, const_S_nocut;
 double n_pi_best, n_k_best, n_p_best, n_L_best, n_S_best;
 double const_pi_best, const_k_best, const_p_best, const_L_best, const_S_best;
 double n_pi[100], n_k[100],n_p[100], n_L[100], n_S[100];
 double mean_pi_nocut, mean_k_nocut, mean_p_nocut, mean_L_nocut, mean_S_nocut;
 double mean_pi_best, mean_k_best, mean_p_best, mean_L_best, mean_S_best;
 double mean_pi[100], mean_k[100],mean_p[100], mean_L[100], mean_S[100];
 double sig_pi_nocut, sig_k_nocut, sig_p_nocut, sig_L_nocut, sig_S_nocut;
 double sig_pi_best, sig_k_best, sig_p_best, sig_L_best, sig_S_best;
 double sig_pi[100], sig_k[100],sig_p[100], sig_L[100], sig_S[100];

//---------------------------------------//
//               Branch                  //
//---------------------------------------//

 int NLtr, NRtr, Ls2_pad[100], Rs2_pad[100];
 double ct, ct_eff;

	double L_tr_chi2;
	double L_tr_x, L_tr_y, L_tr_th, L_tr_ph;
	double L_tr_p;
	double L_tr_tg_th, L_tr_tg_ph;
	double L_tr_vz;
	double L_tr_vz_saved;
	double R_tr_chi2;
	double R_tr_x, R_tr_y, R_tr_th, R_tr_ph;
	double R_tr_p;
	double R_tr_tg_th, R_tr_tg_ph;
	double R_tr_vz;
	double L_mom, R_mom, B_mom; 
	double L_ene, R_ene, B_ene; 
	double ac1sum, ac2sum;//NPE SUM

	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("tr.ntrack_l",1);tree->SetBranchAddress("tr.ntrack_l",&NLtr);
	tree->SetBranchStatus("tr.ntrack_r",1);tree->SetBranchAddress("tr.ntrack_r",&NRtr);
	
  	tree->SetBranchStatus("ac1_npe_sum",1);  tree->SetBranchAddress("ac1_npe_sum", &ac1sum);
  	tree->SetBranchStatus("ac2_npe_sum",1);  tree->SetBranchAddress("ac2_npe_sum", &ac2sum);
  	tree->SetBranchStatus("Lp_c",1);  tree->SetBranchAddress("Lp_c", &L_mom);
  	tree->SetBranchStatus("Rp_c",1);  tree->SetBranchAddress("Rp_c", &R_mom);
  	tree->SetBranchStatus("Bp_c",1);  tree->SetBranchAddress("Bp_c", &B_mom);
  	tree->SetBranchStatus("ct_orig",1);  tree->SetBranchAddress("ct_orig", &ct);

  	tree->SetBranchStatus("L.tr.chi2",1);  tree->SetBranchAddress("L.tr.chi2", &L_tr_chi2);
  	tree->SetBranchStatus("L.tr.x",1);  tree->SetBranchAddress("L.tr.x", &L_tr_x);
  	tree->SetBranchStatus("L.tr.y",1);  tree->SetBranchAddress("L.tr.y", &L_tr_y);
  	tree->SetBranchStatus("L.tr.th",1);  tree->SetBranchAddress("L.tr.th", &L_tr_th);
  	tree->SetBranchStatus("L.tr.ph",1);  tree->SetBranchAddress("L.tr.ph", &L_tr_ph);
  	tree->SetBranchStatus("L.tr.p",1);  tree->SetBranchAddress("L.tr.p", &L_tr_p);
  	tree->SetBranchStatus("L.tr.tg_th",1);  tree->SetBranchAddress("L.tr.tg_th", &L_tr_tg_th );
  	tree->SetBranchStatus("L.tr.tg_ph",1);  tree->SetBranchAddress("L.tr.tg_ph", &L_tr_tg_ph );
  	tree->SetBranchStatus("L.tr.vz",1);  tree->SetBranchAddress("L.tr.vz", &L_tr_vz);

  	tree->SetBranchStatus("R.tr.chi2",1);  tree->SetBranchAddress("R.tr.chi2", &R_tr_chi2);
	tree->SetBranchStatus("R.tr.x" ,1);  tree->SetBranchAddress("R.tr.x" , &R_tr_x );
  	tree->SetBranchStatus("R.tr.y" ,1);  tree->SetBranchAddress("R.tr.y" , &R_tr_y );
  	tree->SetBranchStatus("R.tr.th",1);  tree->SetBranchAddress("R.tr.th", &R_tr_th);
  	tree->SetBranchStatus("R.tr.ph",1);  tree->SetBranchAddress("R.tr.ph", &R_tr_ph);
  	tree->SetBranchStatus("R.tr.p",1);  tree->SetBranchAddress("R.tr.p", &R_tr_p);
  	tree->SetBranchStatus("R.tr.tg_th",1);  tree->SetBranchAddress("R.tr.tg_th", &R_tr_tg_th);
  	tree->SetBranchStatus("R.tr.tg_ph",1);  tree->SetBranchAddress("R.tr.tg_ph", &R_tr_tg_ph);
  	tree->SetBranchStatus("R.tr.vz",1);  tree->SetBranchAddress("R.tr.vz", &R_tr_vz);


  const int slice=8;
  double Nbin_pos[slice+1];
	Nbin_pos[0] = 1730.;
	Nbin_pos[slice] = 1930.;
  for(int i=1;i<slice;i++){
	Nbin_pos[i] = 1730.+200.*(double)(i+1)/(slice+2);
  }
  TH1F* h_path_mom  = new TH1F("h_path_mom;Momentum [MeV/c];Path Length [m]","",slice,Nbin_pos);
  TH1F* h_tdiff_mom  = new TH1F("h_tdiff_mom;Momentum [MeV/c];t_{#pi}-t_{p} [m]","",slice,Nbin_pos);
  TH1F* h_pionpos_mom  = new TH1F("h_pionpos_mom;Momentum [MeV/c];Pion Pos. [ns]","",slice,Nbin_pos);
  TH1F* h_protonpos_mom  = new TH1F("h_protonpos_mom;Momentum [MeV/c];Proton Pos. [ns]","",slice,Nbin_pos);
  TH2F* h_coin_mom  = new TH2F("h_coin_mom;Cointime [ns];Momentum [MeV/c]","",2000,-20.,20.,1000,1730.,1930.);
  TH2F* h_mom_coin  = new TH2F("h_mom_coin;Momentum [MeV/c];Cointime [ns]","",1000,1730.,1930.,2000,-20.,20.);
  TH2F* h_mom_coin_pi  = new TH2F("h_mom_coin_pi;Momentum [MeV/c];Cointime [ns]","",1000,1730.,1930.,2000,-20.,20.);
  TH2F* h_mom_coin_p  = new TH2F("h_mom_coin_p;Momentum [MeV/c];Cointime [ns]","",1000,1730.,1930.,2000,-20.,20.);
  //TH1F* hcoin  = new TH1F("hcoin","",40000/56,-20.,20.);
  //TH1F* hcoin  = new TH1F("hcoin","",130000/56,-30.,100.);
  TH1F* hcoin[slice];
  TH1F* hcoin_pi[slice];
  TH1F* hcoin_p[slice];
	for(int i=0;i<slice;i++){
		//hcoin[i] = new TH1F(Form("hcoin[%d]",i),"",130000/56,-30.,100.);
		hcoin[i] = new TH1F(Form("hcoin[%d]",i),"",40000/56,-20.,20.);
		hcoin_pi[i] = new TH1F(Form("hcoin_pi[%d]",i),"",40000/56,-20.,20.);
		hcoin_p[i] = new TH1F(Form("hcoin_p[%d]",i),"",40000/56,-20.,20.);
	}
  

  TH1F* h_test  = new TH1F("h_test","",1000,1.8,2.4);

  bool L_Tr = false;
  bool L_FP = false;
  bool R_Tr = false;
  bool R_FP = false;
  bool ct_cut = false;
  bool event_selection = false;
  bool event_selection_nocut = false;
  double z_par[100], ac_par[100], ct_par[100];
  double z2_par[100][100], ac2_par[100][100];
  double rf_bunch=2.0;//ns (RF bunch structure)
  const double kcenter = 0.0;
  double mh = ML;//hypernuclei
  double mt = Mp;//target mass
  double B_p, L_p, R_p;//Momentum


  //tree->Draw(">>elist" , "fabs(ct_orig[0][0])<1.0");
  //tree->Draw(">>elist" , "fabs(ct_orig)<1.0");//ctsum (does NOT dintinguish #track)
  //TEventList *elist = (TEventList*)gROOT->FindObject("elist");
  //int ENum = elist->GetN(); 
  int ENum = tree->GetEntries(); 
cout<<"Entries: "<<ENum<<endl;
  int time_div=ENum/25;
  if(ENum<100000)time_div=10000;


	time_t start, end;
	start = time(NULL);
	time(&start);

  for(int i=0;i<ENum;i++){
	//tree->GetEntry(elist->GetEntry(i));
	tree->GetEntry(i);

    if(i%time_div==0){
      end = time(NULL);
      time(&end);
      double diff = difftime(end,start);
      double esttime = diff * ENum / (i+1) - diff;
      cout<<i<<" / "<<ENum<<" ("<<i*100/ENum<<"%) : "<<Form("%.0lf sec passed,  %.0lf sec left",diff,esttime)<<endl;
    }

      
        L_Tr = L_FP = false;
        if( L_tr_chi2<0.01 ) L_Tr = true;
        if( L_tr_th<0.17*L_tr_x+0.025
         && L_tr_th>0.17*L_tr_x-0.035
         && L_tr_th<0.40*L_tr_x+0.130 ) L_FP = true;
	
        R_Tr = R_FP = false;
        // FP and chi2 cuts
        if( R_tr_chi2<0.01 ) R_Tr = true;
        if( R_tr_th<0.17*R_tr_x+0.025
         && R_tr_th>0.17*R_tr_x-0.035
         && R_tr_th<0.40*R_tr_x+0.130 ) R_FP = true;


	

		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;
		
		//if(abs(R_mom-1.90)<0.02)hcoin->Fill(ct);
		for(int i=0;i<slice;i++){
			//if(abs(R_mom-1.76-0.2*i/(slice+2))<0.1/(slice+2))hcoin[i]->Fill(ct);
			if(abs(R_mom-1.73-0.2*(i+1.5)/(slice+2))<0.1/(slice+2))hcoin[i]->Fill(ct);
			if(i==0&&abs(R_mom-1.73-0.1/(slice+2))<0.1/(slice+2))hcoin[i]->Fill(ct);
			if(i==slice-1&&abs(R_mom-1.93+0.1/(slice+2))<0.1/(slice+2))hcoin[i]->Fill(ct);
		//---Pion (AC1 > 3 p.e. && AC2 > 3 p.e.)
			if(abs(R_mom-1.73-0.2*(i+1.5)/(slice+2))<0.1/(slice+2)&&ac1sum>3.&&ac2sum>3.)hcoin_pi[i]->Fill(ct);
			if(i==0&&abs(R_mom-1.73-0.1/(slice+2))<0.1/(slice+2)&&ac1sum>3.&&ac2sum>3.)hcoin_pi[i]->Fill(ct);
			if(i==slice-1&&abs(R_mom-1.93+0.1/(slice+2))<0.1/(slice+2)&&ac1sum>3.&&ac2sum>3.)hcoin_pi[i]->Fill(ct);
		//---Proton (AC1 < 3 p.e. && AC2 < 3 p.e.)
			if(abs(R_mom-1.73-0.2*(i+1.5)/(slice+2))<0.1/(slice+2)&&ac1sum<3.&&ac2sum<3.)hcoin_p[i]->Fill(ct);
			if(i==0&&abs(R_mom-1.73-0.1/(slice+2))<0.1/(slice+2)&&ac1sum<3.&&ac2sum<3.)hcoin_p[i]->Fill(ct);
			if(i==slice-1&&abs(R_mom-1.93+0.1/(slice+2))<0.1/(slice+2)&&ac1sum<3.&&ac2sum<3.)hcoin_p[i]->Fill(ct);
		}
		h_coin_mom->Fill(ct,R_mom*1000.);
		h_mom_coin->Fill(R_mom*1000.,ct);
		if(ac1sum>3.&&ac2sum>3.&&abs(ct-3.15)<1.)h_mom_coin_pi->Fill(R_mom*1000.,ct);
		if(ac1sum<3.&&ac2sum<3.&&ct>10.*R_mom-27.3&&ct<10.*R_mom-25.3)h_mom_coin_p->Fill(R_mom*1000.,ct);

}//ENum

	double chisq = 0.;
	double dof = 0.;
	double pion_par0[slice] = {0}; 
	double pion_par1[slice] = {0};
	double pion_par2[slice] = {0};
	double proton_par0[slice] = {0}; 
	double proton_par1[slice] = {0};
	double proton_par2[slice] = {0};
	double pion_parerr0[slice] = {0}; 
	double pion_parerr1[slice] = {0};
	double pion_parerr2[slice] = {0};
	double proton_parerr0[slice] = {0}; 
	double proton_parerr1[slice] = {0};
	double proton_parerr2[slice] = {0};
	double t_diff[slice] = {0};
	double pathlen[slice] = {0};
	double t_diff_sig[slice] = {0};
	double pathlen_sig[slice] = {0};

		TCanvas *c4 = new TCanvas("c4", "c4", 800, 800);
//Pion Fitting
cout<<"Pion Fitting is performed as a first step"<<endl;
cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
//	 TF1 *fcoin_pion= new TF1("fcoin_pion","gausn",1.,5.);
  TF1* fcoin_pion[slice]; 
  TF1 *fcoin_proton[slice];
	for(int i=0;i<slice;i++){
		fcoin_pion[i] = new TF1(Form("fcoin_pion[%d]",i),"gausn",1.,5.);
	 fcoin_pion[i]->SetNpx(20000);
	 fcoin_pion[i]->SetParameter(0,100000./slice);
	 fcoin_pion[i]->SetParameter(1,3.15);
	 fcoin_pion[i]->SetParameter(2,0.4);
	 //hcoin[i]->Fit(Form("fcoin_pion[%d]",i),"","",2.,4.);
	 hcoin_pi[i]->Fit(Form("fcoin_pion[%d]",i),"","",2.,4.);
	 pion_par0[i] = fcoin_pion[i]->GetParameter(0);
	 pion_par1[i] = fcoin_pion[i]->GetParameter(1);
	 pion_par2[i] = fcoin_pion[i]->GetParameter(2);
	 pion_parerr0[i] = fcoin_pion[i]->GetParError(0);
	 pion_parerr1[i] = fcoin_pion[i]->GetParError(1);
	 pion_parerr2[i] = fcoin_pion[i]->GetParError(2);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 //cout<<"Pion Fitting after BG subtraction (fcoin_pion) "<<endl;
	 cout<<"Pion Fitting w/o BG subtraction (fcoin_pion) "<<endl;
	 chisq = fcoin_pion[i]->GetChisquare();
	 dof  = fcoin_pion[i]->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

//Proton Fitting
cout<<"Proton Fitting is performed as a next step"<<endl;
cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
	 fcoin_proton[i]= new TF1(Form("fcoin_proton[%d]",i),"gausn",-10,-5.);
	 fcoin_proton[i]->SetNpx(20000);
	 fcoin_proton[i]->SetParameter(0,2000.);
	 if(i<slice/4){
		fcoin_proton[i]->SetParameter(1,-8.5);
	 }else if(i<slice*2/4){
		fcoin_proton[i]->SetParameter(1,-8.0);
	 }else if(i<slice*3/4){
		fcoin_proton[i]->SetParameter(1,-7.5);
	 }else{
		fcoin_proton[i]->SetParameter(1,-7.0);
	 }
	 fcoin_proton[i]->SetParameter(2,0.4);
	 if(i<slice/8){
		//hcoin[i]->Fit(Form("fcoin_proton[%d]",i),"","",-10.,-8.);
		hcoin_p[i]->Fit(Form("fcoin_proton[%d]",i),"","",-10.,-8.);
	 }else if(i<slice*2/4){
		//hcoin[i]->Fit(Form("fcoin_proton[%d]",i),"","",-9.,-7.);
		hcoin_p[i]->Fit(Form("fcoin_proton[%d]",i),"","",-9.,-7.);
	 }else if(i<slice*3/4){
		//hcoin[i]->Fit(Form("fcoin_proton[%d]",i),"","",-8.5,-6.5);
		hcoin_p[i]->Fit(Form("fcoin_proton[%d]",i),"","",-8.5,-6.5);
	 }else{
		//hcoin[i]->Fit(Form("fcoin_proton[%d]",i),"","",-8.,-6.);
		hcoin_p[i]->Fit(Form("fcoin_proton[%d]",i),"","",-8.,-6.);
	 }
	 proton_par0[i] = fcoin_proton[i]->GetParameter(0);
	 proton_par1[i] = fcoin_proton[i]->GetParameter(1);
	 proton_par2[i] = fcoin_proton[i]->GetParameter(2);
	 proton_parerr0[i] = fcoin_proton[i]->GetParError(0);
	 proton_parerr1[i] = fcoin_proton[i]->GetParError(1);
	 proton_parerr2[i] = fcoin_proton[i]->GetParError(2);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 //cout<<"Proton Fitting after BG subtraction (fcoin_proton) "<<endl;
	 cout<<"Proton Fitting w/o BG subtraction (fcoin_proton) "<<endl;
	 chisq = fcoin_proton[i]->GetChisquare();
	 dof  = fcoin_proton[i]->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;
	}//i loop for fitting


	TCanvas *c5 = new TCanvas("c5", "c5", 800, 800);
	c5->Divide(slice/2,2);
	for(int i=0;i<slice;i++){
		//c5->cd(i+1)->DrawFrame(-20.,0.,20.,8000.);
		c5->cd(i+1);
		hcoin[i]->SetLineColor(kAzure);
		hcoin[i]->Draw("same");
		fcoin_pion[i]->SetLineColor(kOrange);
		fcoin_pion[i]->Draw("same");
		fcoin_proton[i]->SetLineColor(kRed);
		fcoin_proton[i]->Draw("same");
	}

	double mom = -1.82;//GeV/c
	double betapi = mom/sqrt(Mpi*Mpi+mom*mom);
	double betap  = mom/sqrt(Mp*Mp+mom*mom);
	double betapi_sig = 0.02/sqrt(Mpi*Mpi+mom*mom);
	double betap_sig  = 0.02/sqrt(Mp*Mp+mom*mom);
	for(int i=0;i<slice;i++){
		if(i==0){
		mom = 1.73+0.2/(slice+2);
		}else if(i==slice-1){
		mom = 1.93-0.2/(slice+2);;
		}else{
		mom = 1.73+0.2*(i+1.5)/(slice+2);//GeV/c
		}
		cout<<"mom="<<mom<<endl;
		betapi = mom/sqrt(Mpi*Mpi+mom*mom);
		betap  = mom/sqrt(Mp*Mp+mom*mom);
		if(i==0||i==slice-1){
		betapi_sig = 0.4/(slice+2)/pow(sqrt(Mpi*Mpi+mom*mom),3.);
		betap_sig  = 0.4/(slice+2)/pow(sqrt(Mp*Mp+mom*mom),3.);
		}else{
		betapi_sig = 0.2/(slice+2)/pow(sqrt(Mpi*Mpi+mom*mom),3.);
		betap_sig  = 0.2/(slice+2)/pow(sqrt(Mp*Mp+mom*mom),3.);
		}
		cout<<"betapi="<<betapi<<endl;
		cout<<"betap="<<betap<<endl;
		//cout<<"betapi_sig="<<betapi_sig<<endl;
		//cout<<"betap_sig="<<betap_sig<<endl;
		t_diff[i] = proton_par1[i]-pion_par1[i];
		pathlen[i] = betapi*betap*LightVelocity*t_diff[i]/(betap-betapi);
		//t_diff_sig[i] = sqrt(proton_par2[i]*proton_par2[i]+pion_par2[i]*pion_par2[i]);
		//t_diff_sig[i] = sqrt(proton_par2[i]*proton_par2[i]/proton_par0[i]*0.056+pion_par2[i]*pion_par2[i]/pion_par0[i]*0.056);
		t_diff_sig[i] = sqrt(pion_parerr1[i]*pion_parerr1[i]+proton_parerr1[i]*proton_parerr1[i]);
		//pathlen_sig[i] = betapi*betap*LightVelocity*t_diff_sig[i]/(betap-betapi);
		pathlen_sig[i] = sqrt(betapi*betap*LightVelocity/(betap-betapi)*betapi*betap*LightVelocity/(betap-betapi)*t_diff_sig[i]*t_diff_sig[i]+pow(betap,4.)*LightVelocity*LightVelocity/(std::pow(betapi-betap,4))*(betapi_sig*betapi_sig)+pow(betapi,4.)*LightVelocity*LightVelocity/(std::pow(betapi-betap,4))*(betap_sig*betap_sig));
	cout<<"t_diff["<<i<<"]="<<t_diff[i]<<"+/-"<<abs(t_diff_sig[i])<<" [ns]"<<endl;
	cout<<"pathlen["<<i<<"]="<<pathlen[i]<<"+/-"<<abs(pathlen_sig[i])<<" [m]"<<endl;
	//cout<<"t_diff_clac="<<t_diff[i]<<"+/-"<<sqrt(pion_parerr1[i]*pion_parerr1[i]+proton_parerr1[i]*proton_parerr1[i])<<endl;
		h_path_mom->SetBinContent(i+1,pathlen[i]);
		h_path_mom->SetBinError(i+1,pathlen_sig[i]);
		h_tdiff_mom->SetBinContent(i+1,abs(t_diff[i]));
		h_tdiff_mom->SetBinError(i+1,t_diff_sig[i]);
		h_pionpos_mom->SetBinContent(i+1,pion_par1[i]);
		h_pionpos_mom->SetBinError(i+1,pion_parerr1[i]);
		h_protonpos_mom->SetBinContent(i+1,proton_par1[i]);
		h_protonpos_mom->SetBinError(i+1,proton_parerr1[i]);
	}

	TCanvas *c6 = new TCanvas("c6", "c6", 800, 800);
	//c6->DrawFrame(1730.,26.3,1930.,27.3);
	h_path_mom->Draw("esame");
	TCanvas *c7 = new TCanvas("c7", "c7", 800, 800);
	//c7->DrawFrame(1730.,26.3,1930.,27.3);
	h_tdiff_mom->Draw("esame");
cout<<"pion positon;"<<endl;
	for(int i=0;i<slice;i++){
cout<<pion_par1[i]<<endl;
	}
cout<<"proton positon;"<<endl;
	for(int i=0;i<slice;i++){
cout<<proton_par1[i]<<endl;
	}
	TCanvas *c8 = new TCanvas("c8", "c8", 800, 800);
	h_pionpos_mom->Draw("e");
	TCanvas *c9 = new TCanvas("c9", "c9", 800, 800);
	h_protonpos_mom->Draw("e");
	TCanvas *c10 = new TCanvas("c10", "c10", 800, 800);
	//h_coin_mom->Draw("colz");
	h_mom_coin->Draw("colz");
	//TF1* fpi_line=new TF1("fpi_line","pol1",1730.,1930.);
	TF1* fpi_line=new TF1("fpi_line","pol1",1760.,1900.);
	fpi_line->SetNpx(2000);
	TProfile* ppi = h_mom_coin_pi->ProfileX();
	TProfile* pp = h_mom_coin_p->ProfileX();
	fpi_line->SetParameter(0,3.15);//intercept
	//fpi_line->SetParLimits(0,3.0,3.20);
	fpi_line->SetParameter(1,0.00001);//slope
	//TF1* fp_line=new TF1("fp_line","pol1",1730.,1930.);
	TF1* fp_line=new TF1("fp_line","pol1",1760.,1900.);
	fp_line->SetNpx(2000);
	fp_line->SetParameter(0,-26.3);//intercept
	fp_line->SetParameter(1,0.01);//slope
	//h_mom_coin_pi->Fit("fpi_line");
	ppi->Fit("fpi_line");
	//h_mom_coin_pi->Fit("fp_line");
	pp->Fit("fp_line");
	h_mom_coin->Draw("colz");
	fpi_line->Draw("same");
	fp_line->Draw("same");
	TCanvas *c11 = new TCanvas("c11", "c11", 800, 800);
	ppi->Draw("");
	//h_mom_coin_pi->Draw("colz");
	TCanvas *c12 = new TCanvas("c12", "c12", 800, 800);
	pp->Draw("");
	//h_mom_coin_p->Draw("colz");
	
	TCanvas *c13 = new TCanvas("c13", "c13", 800, 800);
	TF1* f_tdiff=new TF1("f_tdiff","pol1",1730.,1930.);
	f_tdiff->SetParameter(0,fpi_line->GetParameter(0)-fp_line->GetParameter(0));
	f_tdiff->SetParameter(1,fpi_line->GetParameter(1)-fp_line->GetParameter(1));
	h_tdiff_mom->Draw("e");
	f_tdiff->Draw("same");

	TCanvas *c14 = new TCanvas("c14", "c14", 800, 800);
	TF2* f_path = new TF2("f_path", F_path, 1730.,1930.,10.,12.);
	f_path->Draw("colz");
	f_tdiff->Draw("same");

	TCanvas *c15 = new TCanvas("c15", "c15", 800, 800);
	TF1* f1_path = new TF1("f1_path", F1_path, 1730.,1930.,2);
	f1_path->SetParameter(0,fpi_line->GetParameter(0)-fp_line->GetParameter(0));
	f1_path->SetParameter(1,fpi_line->GetParameter(1)-fp_line->GetParameter(1));
	//f1_path->Draw("");
	//h_path_mom->Draw("esame");
	TF1* f1_path_uerr = new TF1("f1_path_uerr", F1_path, 1730.,1930.,2);
	TF1* f1_path_derr = new TF1("f1_path_derr", F1_path, 1730.,1930.,2);
	f1_path_derr->SetParameter(0,fpi_line->GetParameter(0)+fpi_line->GetParError(0)-fp_line->GetParameter(0)-fp_line->GetParError(0));
	f1_path_derr->SetParameter(1,fpi_line->GetParameter(1)+fpi_line->GetParError(1)-fp_line->GetParameter(1)-fp_line->GetParError(1));
	f1_path_uerr->SetParameter(0,fpi_line->GetParameter(0)-fpi_line->GetParError(0)-fp_line->GetParameter(0)+fp_line->GetParError(0));
	f1_path_uerr->SetParameter(1,fpi_line->GetParameter(1)-fpi_line->GetParError(1)-fp_line->GetParameter(1)+fp_line->GetParError(1));
//	f1_path_uerr->SetLineColor(kAzure);
//	f1_path_derr->SetLineColor(kAzure);
//	f1_path_uerr->Draw("same");
//	f1_path_derr->Draw("same");
	f1_path_derr->SetFillColor(kBlue);
	f1_path_derr->SetFillStyle(3001);
	func_shade(c15,f1_path_derr,f1_path_uerr,f1_path);
	h_path_mom->Draw("esame");

//SURVIVAL RATIO
	TCanvas *c19 = new TCanvas("c19", "c19", 800, 800);
	TF1* f1_sr = new TF1("f1_sr", F1_sr, 1730.,1930.,2);
	f1_sr->SetParameter(0,fpi_line->GetParameter(0)-fp_line->GetParameter(0));
	f1_sr->SetParameter(1,fpi_line->GetParameter(1)-fp_line->GetParameter(1));
	//f1_sr->Draw("");
	//h_path_mom->Draw("esame");
	TF1* f1_sr_uerr = new TF1("f1_sr_uerr", F1_sr, 1730.,1930.,2);
	TF1* f1_sr_derr = new TF1("f1_sr_derr", F1_sr, 1730.,1930.,2);
	f1_sr_derr->SetParameter(0,fpi_line->GetParameter(0)+fpi_line->GetParError(0)-fp_line->GetParameter(0)-fp_line->GetParError(0));
	f1_sr_derr->SetParameter(1,fpi_line->GetParameter(1)+fpi_line->GetParError(1)-fp_line->GetParameter(1)-fp_line->GetParError(1));
	f1_sr_uerr->SetParameter(0,fpi_line->GetParameter(0)-fpi_line->GetParError(0)-fp_line->GetParameter(0)+fp_line->GetParError(0));
	f1_sr_uerr->SetParameter(1,fpi_line->GetParameter(1)-fpi_line->GetParError(1)-fp_line->GetParameter(1)+fp_line->GetParError(1));
	f1_sr_derr->SetFillColor(kBlue);
	f1_sr_derr->SetFillStyle(3001);
	func_shade(c19,f1_sr_derr,f1_sr_uerr,f1_sr);

//
//void fshade() {
//TF1 *f1 = new TF1(¡Èf1¡É,¡È1/x +0.04¡É,1,5);
//TF1 *f2 = new TF1(¡Èf2¡É,¡È1/x -0.04¡É,1,5);
//TF1 *f3 = new TF1(¡Èf3¡É,¡È1/x +0¡É,1,5);
//f1->SetFillColor(kBlue);
//f1->SetFillStyle(3001);
//TCanvas *c1 = new TCanvas(¡Èc1¡É);
//shade(c1,f1,f2,f3);
//}
	//double t_diff_sig = sqrt(proton_par2*proton_par2+pion_par2*pion_par2);
	//double t_diff_sig2 = sqrt(proton_par2*proton_par2/proton_par0*0.056+pion_par2*pion_par2/pion_par0*0.056);
	//double pathlen_sig = betapi*betap*LightVelocity*t_diff_sig/(betap-betapi);
	//double pathlen_sig2 = betapi*betap*LightVelocity*t_diff_sig2/(betap-betapi);
	//cout<<"Npi="<<pion_par0/0.056<<endl;
	//cout<<"Np="<<proton_par0/0.056<<endl;
	//cout<<"betapi="<<betapi<<endl;
	//cout<<"betap="<<betap<<endl;
	//cout<<"t_diff="<<abs(t_diff)<<" [ns]"<<endl;
	//cout<<"pathlen="<<pathlen<<"+/-"<<abs(pathlen_sig)<<" [m]"<<endl;
	//cout<<"pathlen2="<<pathlen<<"+/-"<<abs(pathlen_sig2)<<" [m]"<<endl;
	
  TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/test0811_1.root","read");//
  TTree *tree_simc = (TTree*)file_simc->Get("SNT");

//---------------------------------------//
//               Branch                  //
//---------------------------------------//


	float R_mom_simc; 
	float R_tr_pathl;
	float R_x;//h_xfp
	float R_y;//h_yfp
	float R_dx;//h_xpfp
	float R_dy;//h_ypfp
	float R_ksr;

	tree_simc->SetBranchStatus("*",0);
  	tree_simc->SetBranchStatus("Rp_rec",1);  tree_simc->SetBranchAddress("Rp_rec", &R_mom_simc);
    tree_simc->SetBranchStatus("h_xfp",1);  tree_simc->SetBranchAddress("h_xfp", &R_x);
    tree_simc->SetBranchStatus("h_yfp",1);  tree_simc->SetBranchAddress("h_yfp", &R_y);
    tree_simc->SetBranchStatus("h_xpfp",1);  tree_simc->SetBranchAddress("h_xpfp", &R_dx);
    tree_simc->SetBranchStatus("h_ypfp",1);  tree_simc->SetBranchAddress("h_ypfp", &R_dy);
    tree_simc->SetBranchStatus("PathLength",1);  tree_simc->SetBranchAddress("PathLength", &R_tr_pathl);
    tree_simc->SetBranchStatus("SurvivalRatio",1);  tree_simc->SetBranchAddress("SurvivalRatio", &R_ksr);




  //TH1F* h1  = new TH1F("h1","",400,-20.,20.0);
  //h1->GetXaxis()->SetTitle("coin time (ns)");
  //h1->GetYaxis()->SetTitle("Counts / 100 ps");
  //h1->GetXaxis()->SetRangeUser(-14.0,17.);
  //h1 ->SetLineColor(2);
  //h1->SetLineWidth(2);
  TH1F* hR_mom  = new TH1F("hR_mom","hR_mom",1000,1.7,2.0);
  TH1F* hR2_mom  = new TH1F("hR2_mom","hR2_mom",1000,1.7,2.0);
  TH1F* hR_T2S2  = new TH1F("hR_T2S2","hR_T2S2",1000,26.,29.);
  TH2F* hR_T2S2_mom  = new TH2F("hR_T2S2_mom","hR_T2S2_mom",100,25.,27.,100,1.7,2.);
  TH2F* hR2_T2S2_mom  = new TH2F("hR2_T2S2_mom","hR2_T2S2_mom",100,25.,27.,100,1.76,1.9);

  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  TH2F* hR_mom_x  = new TH2F("hR_mom_x","hR_mom_x;Momentum [GeV/c];X(FP) [m]",100,1.76,1.9,100,-1.0,1.0);
  TH2F* hR_mom_th  = new TH2F("hR_mom_th","hR_mom_th;Momentum [GeV/c];X'(FP) [rad]",100,1.76,1.9,100,-0.15,0.15);
  TH2F* hR_len_x  = new TH2F("hR_len_x","hR_len_x;Path Length (T to S2) [m];X(FP) [m]",100,25.,27.,100,-1.0,1.0);
  TH2F* hR_len_th  = new TH2F("hR_len_th","hR_len_th;Path Length (T to S2) [m];X'(FP) [rad]",100,25.,27.,100,-0.15,0.15);
  TH2F* hR_lenfp_x  = new TH2F("hR_lenfp_x","hR_lenfp_x;Path Length (T to FP) [m];X(FP) [m]",100,21.,25.,100,-1.0,1.0);
  TH2F* hR_lenfp_th  = new TH2F("hR_lenfp_th","hR_lenfp_th;Path Length (T to FP) [m];X'(FP) [rad]",100,21.,25.,100,-0.15,0.15);
  TH2F* hR_mom_y  = new TH2F("hR_mom_y","hR_mom_y;Momentum [GeV/c];Y(FP) [m]",100,1.76,1.9,100,-0.1,0.1);
  TH2F* hR_mom_ph  = new TH2F("hR_mom_ph","hR_mom_ph;Momentum [GeV/c];Y'(FP) [rad]",100,1.76,1.9,100,-0.25,0.25);
  TH2F* hR_len_y  = new TH2F("hR_len_y","hR_len_y;Path Length (T to S2) [m];Y(FP) [m]",100,25.,27.,100,-0.1,0.1);
  TH2F* hR_len_ph  = new TH2F("hR_len_ph","hR_len_ph;Path Length (T to S2) [m];Y'(FP) [rad]",100,25.,27.,100,-0.25,0.25);

  TH2F* h2_kSR  = new TH2F("h2_kSR","h2_kSR",100,25.,27.,100,1.7,2.);
  int srbin = 80;
  TH1D* h_kSR  = new TH1D("h_kSR","h_kSR",srbin,1.65,2.);
  TH1D* h_kSR2  = new TH1D("h_kSR2","h_kSR2 from SIMC",srbin,1.65,2.);
  double kaonsr[srbin];
  double ave[srbin];//Ensemble average: \bar{X} 
  double dif[srbin];// (X_i - \bar{X})^2
  //double sigma[srbin];//sigma of ave. \sigma = \sqrt{\frac{\sum_i^n (X_i-\bar{X})^2}{n(n-1)}}
  int counts[srbin];
  for(int i=0;i<srbin;i++){//Initialization
	kaonsr[i]=0.;
	ave[i]=0.;
	counts[i]=0;
	dif[i]=0.;
  }

//---  Loading... ---//average of ksr
	string daq_file = "./simc_ksr_tmp.dat";//input file
	int binnum;//running bin
	int nbin;//total # of bin
	double ave_tmp;
	double sig_tmp;
	for(int i=0;i<srbin;i++){
		ave[i]=0.;
	}

	string buf;
	ifstream ifp(daq_file.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed. Please make a input file." << endl; exit(1);}
cout << "Param file : " << daq_file.c_str() << endl;
	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(buf=="Number of Bin"){
			getline(ifp,buf);
			stringstream sbinbuf(buf);
			sbinbuf >> nbin;
			//cout << nbin <<endl;
			if(nbin==srbin)continue;
			else{
				cout<<"Number of Bins is different!"<<endl;
				cout<<"srbin should be set "<<nbin<<endl;
				break;
			}
		}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> binnum >> ave_tmp;
		ave[binnum] = ave_tmp;
		cout << binnum<< ", " << ave[binnum] <<endl;
	}
	cout<<"Average of SR was successfully loaded from external file."<<endl;
	


  int ENum_simc = tree_simc->GetEntries(); 
cout<<"Entries: "<<ENum_simc<<endl;
  int time_div_simc=ENum_simc/25;
  if(ENum_simc<100000)time_div_simc=10000;


	time_t start_simc, end_simc;
	start_simc = time(NULL);
	time(&start_simc);

  for(int i=0;i<ENum_simc;i++){
	//tree->GetEntry(elist->GetEntry(i));
	tree_simc->GetEntry(i);

	R_mom_simc /= 1000.;// MeV/c-->GeV/c
	R_tr_pathl /= 100.;// cm-->m

    if(i%time_div_simc==0){
      end_simc = time(NULL);
      time(&end_simc);
      double diff_simc = difftime(end_simc,start_simc);
      double esttime = diff_simc * ENum_simc / (i+1) - diff_simc;
      cout<<i<<" / "<<ENum_simc<<" ("<<i*100/ENum_simc<<"%) : "<<Form("%.0lf sec passed,  %.0lf sec left",diff_simc,esttime)<<endl;
    }

			//R_tr_pathl += -1.83*(1.+R_dx*R_dx+R_dy*R_dy);// 183cm (Z_aerogel)
			//R_tr_pathl += -0.5845*(1.+R_dx*R_dx+R_dy*R_dy);// (Final -> S2)
			hR_mom->Fill(R_mom_simc);
			hR_T2S2->Fill(R_tr_pathl);
			hR_T2S2_mom->Fill(R_tr_pathl,R_mom_simc);

			hR_mom_x->Fill(R_mom_simc,R_x*0.01);
			hR_mom_th->Fill(R_mom_simc,R_dx);
			hR_len_x->Fill(R_tr_pathl,R_x*0.01);
			hR_len_th->Fill(R_tr_pathl,R_dx);
			hR_mom_y->Fill(R_mom_simc,R_y*0.01);
			hR_mom_ph->Fill(R_mom_simc,R_dy);
			hR_len_y->Fill(R_tr_pathl,R_y*0.01);
			hR_len_ph->Fill(R_tr_pathl,R_dy);

			double ksr = exp(-1.*R_tr_pathl*MK/R_mom_simc/3.713);
//cout<<"survival ratio = "<<ksr<<endl;
//cout<<"R_tr_pathl = "<<R_tr_pathl<<endl;
//cout<<"R_mom = "<<R_mom<<endl;
			//ksr = R_ksr;
			if(R_mom_simc>1.7&&R_mom_simc<2.0&&R_tr_pathl>25.&&R_tr_pathl<29.){	
			kaonsr[(int)((R_mom_simc-1.65)*srbin/0.35)] += ksr;
			counts[(int)((R_mom_simc-1.65)*srbin/0.35)] ++;
			double ksrave = ave[(int)((R_mom_simc-1.65)*srbin/0.35)];
			dif[(int)((R_mom_simc-1.65)*srbin/0.35)] += (ksr - ksrave)*(ksr - ksrave);
			}
			h2_kSR->SetBinContent(h2_kSR->GetXaxis()->FindBin(R_tr_pathl),h2_kSR->GetYaxis()->FindBin(R_mom_simc),ksr);
			h_kSR2->Fill(R_ksr);
//cout<<"SR= "<<ksr<<endl;
	


}//ENum

double tmp,tmp2;
double ave_new[srbin];
  for(int i=0;i<srbin;i++){
	if(counts[i]!=0){
		tmp = kaonsr[i]/counts[i];//average
		ave_new[i] = tmp;
	}
	else tmp = 0.;
	if(counts[i]!=0&&counts[i]!=1) tmp2 = sqrt(dif[i]/counts[i]/(counts[i]-1));//sigma of ave.
	else tmp2 = 0.;
//cout<<"survival ratio is "<<tmp<<endl;
//cout<<"kaonsr["<<i<<"] = "<<kaonsr[i]<<endl;
//cout<<"counts["<<i<<"] = "<<counts[i]<<endl;
	h_kSR->SetBinContent(i+1,tmp);
	h_kSR->SetBinError(i+1,tmp2);
  }

	TCanvas* c60 = new TCanvas("c60","c60");
	TH1F *frame = c60->DrawFrame(1.65,0.,2.,0.20);
	h_kSR->SetLineColor(kGreen);
	h_kSR->SetLineWidth(2);
	h_kSR->Draw("esame");

	TCanvas* c70 = new TCanvas("c70","c70");
	TH1F *frame2 = c70->DrawFrame(1.65,0.,2.,0.20);
	h_kSR2->SetLineColor(kRed);
	h_kSR2->SetLineWidth(2);
	h_kSR2->Draw("esame");

	TCanvas *c80 = new TCanvas("c80", "c80", 800, 800);
	TH1D* h_kSR_simc  = new TH1D("h_kSR_simc","h_kSR_simc",srbin,1650.,2000.);
    SetTH1(h_kSR_simc, "", "Momentum [MeV/c]", "#varepsilon^{SR}", kAzure, kRed, kRed);
	func_shade(c80,f1_sr_derr,f1_sr_uerr,f1_sr);
	double syst_max1 = 0.;//lower
	double syst_max2 = 0.;//higher
	double df_temp = 0.;
	double dh_temp = 0.;
	for(int i=0;i<srbin;i++){
		h_kSR_simc->SetBinContent(i+1,h_kSR->GetBinContent(i+1));
		h_kSR_simc->SetBinError(i+1,h_kSR->GetBinError(i+1));
		dh_temp = (h_kSR_simc->GetBinContent(i+1))-(f1_sr->Eval(h_kSR_simc->GetBinCenter(i+1)));
		dh_temp = 100.*dh_temp/(f1_sr->Eval(h_kSR_simc->GetBinCenter(i+1)));
		if(syst_max1<dh_temp)syst_max1 = dh_temp;
		df_temp = (f1_sr->Eval(h_kSR_simc->GetBinCenter(i+1)))-(f1_sr_uerr->Eval(h_kSR_simc->GetBinCenter(i+1)));
		df_temp = 100.*df_temp/(f1_sr->Eval(h_kSR_simc->GetBinCenter(i+1)));
		if(syst_max2<df_temp)syst_max2 = df_temp;
	//cout << "f1_sr_derr = " << f1_sr_derr->Eval(h_kSR_simc->GetBinCenter(i+1)) <<endl;
	//cout << "f1_sr_uerr = " << f1_sr_uerr->Eval(h_kSR_simc->GetBinCenter(i+1)) <<endl;
	}
	cout << "SR(SIMC) - SR(Data) = " << syst_max1  << "%" <<endl;
	cout << "SR(Data) - min.SR(Data) = " << syst_max2 << "%" <<endl;
	h_kSR_simc->Draw("esame");
	c80->SetLeftMargin(0.18);
	c80->SetRightMargin(0.10);
	c80->SetTopMargin(0.14);
	c80->SetBottomMargin(0.14);
	c80->Modified();
	c80->Update();
	gPad->Modified();
	gPad->Update();
	c80->Print("./pdf/kaonSR.pdf");

	ofstream fout("./simc_ksr_tmp.dat");
		fout<<"#from PathLength/kaonSR_2022.C"<<endl;
		fout<<"Number of Bin"<<endl;
		fout<<srbin<<endl;
		fout<<"#"<<endl;
		fout<<"#"<<endl;
		fout<<"#"<<endl;
		double vpflux_temp;
		double vpflux_total=0.;
	for(int i=0; i<srbin; i++){
		fout<<i<<" "<<ave_new[i]<<endl;
	}
	cout<<"Average of each bin was successfully written in ./simc_ksr_tmp.dat."<<endl;

}//path_from_ct
