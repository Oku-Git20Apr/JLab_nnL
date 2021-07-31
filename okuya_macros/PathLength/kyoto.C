//-----------------------------//
//--  Path Length Calc.      --//
//--	from Cointime        --//
//-----------------------------//

//- taken from mom_calib.C
//K. Okuyama (Jul. 29, 2021)
// -- suggested by Nagao-san after ELS#143

double PI=4.*atan(1.);

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
gr->Draw("f"); //draw graph with fill area option
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

void kyoto(){
//ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);
	string pdfname = "temp.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  //TFile *file = new TFile("../h2all5.root","read");//input file of all H2 run(default: h2all4.root)
  //TFile *file = new TFile("../h2all_2020Nov.root","read");//input file of all H2 run(default: h2all4.root)
  TFile *file = new TFile("../Kyoto_root/h2_suzuki_nocut.root","read");//input file of all H2 run(default: h2all4.root)
  TTree *tree = (TTree*)file->Get("tree");

    
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
 double ct_wseg[16][16];
for(int i=0;i<16;i++){
	for(int j=0;j<16;j++){
		ct_wseg[i][j]=-1000.;
	}
}
 int LS2_seg;
 int RS2_seg;

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
	
  	tree->SetBranchStatus("AC1",1);  tree->SetBranchAddress("AC1", &ac1sum);
  	tree->SetBranchStatus("AC2",1);  tree->SetBranchAddress("AC2", &ac2sum);
  	tree->SetBranchStatus("Lp",1);  tree->SetBranchAddress("Lp", &L_mom);
  	tree->SetBranchStatus("Rp",1);  tree->SetBranchAddress("Rp", &R_mom);
  	tree->SetBranchStatus("ctime",1);  tree->SetBranchAddress("ctime", &ct);

  	tree->SetBranchStatus("Lth",1);  tree->SetBranchAddress("Lth", &L_tr_th);
  	tree->SetBranchStatus("Lph",1);  tree->SetBranchAddress("Lph", &L_tr_ph);
  	tree->SetBranchStatus("Lvz",1);  tree->SetBranchAddress("Lvz", &L_tr_vz);

  	tree->SetBranchStatus("Rth",1);  tree->SetBranchAddress("Rth", &R_tr_th);
  	tree->SetBranchStatus("Rph",1);  tree->SetBranchAddress("Rph", &R_tr_ph);
  	tree->SetBranchStatus("Rvz",1);  tree->SetBranchAddress("Rvz", &R_tr_vz);


  const int slice=6;
  double Nbin_pos[slice+1];
	Nbin_pos[0] = 1730.;
	Nbin_pos[slice] = 1930.;
  for(int i=1;i<slice;i++){
	Nbin_pos[i] = 1730.+200.*(double)(i+1)/(slice+2);
  }
  TH1F* h_path_mom  = new TH1F("h_path_mom;Momentum","h_path_mom;Momentum [MeV/c];Path Length [m]",slice,Nbin_pos);
  TH1F* h_tdiff_mom  = new TH1F("h_tdiff_mom;Momentum","h_tdiff_mom;Momentum [MeV/c];t_{#pi}-t_{p} [m]",slice,Nbin_pos);
  TH1F* h_mom_from_ct  = new TH1F("h_mom_from_ct","Mom_from_ct;Momentum [MeV/c];Momentum [MeV/c]",slice,Nbin_pos);
  TH1F* h_mom_from_cte1  = new TH1F("h_mom_from_cte1","Mom_from_ct;Momentum [MeV/c];Momentum [MeV/c]",slice,Nbin_pos);
  TH1F* h_mom_from_cte2  = new TH1F("h_mom_from_cte2","Mom_from_ct;Momentum [MeV/c];Momentum [MeV/c]",slice,Nbin_pos);
  TH1F* h_mom_from_cte3  = new TH1F("h_mom_from_cte3","Mom_from_ct;Momentum [MeV/c];Momentum [MeV/c]",slice,Nbin_pos);
  TH1F* h_mom_from_cte4  = new TH1F("h_mom_from_cte4","Mom_from_ct;Momentum [MeV/c];Momentum [MeV/c]",slice,Nbin_pos);
  h_mom_from_cte1->SetLineColor(kRed);
  h_mom_from_cte2->SetLineColor(kAzure);
  h_mom_from_cte3->SetLineColor(kGreen);
  h_mom_from_cte4->SetLineColor(kOrange);
  TH1F* h_pionpos_mom  = new TH1F("h_pionpos_mom","h_pionpos_mom;Momentum [MeV/c];Pion Pos. [ns]",slice,Nbin_pos);
  TH1F* h_protonpos_mom  = new TH1F("h_protonpos_mom","h_protonpos_mom;Momentum [MeV/c];Proton Pos. [ns]",slice,Nbin_pos);
  TH2F* h_coin_mom  = new TH2F("h_coin_mom","h_coin_mom;Cointime [ns];Momentum [MeV/c]",2000,-20.,20.,1000,1730.,1930.);
  TH2F* h_mom_coin  = new TH2F("h_mom_coin","h_mom_coin;Momentum [MeV/c];Cointime [ns]",1000,1730.,1930.,2000,-20.,20.);
  TH2F* h_mom_coin_pi  = new TH2F("h_mom_coin_pi","h_mom_coin_pi;Momentum [MeV/c];Cointime [ns]",1000,1730.,1930.,2000,-20.,20.);
  TH2F* h_mom_coin_p  = new TH2F("h_mom_coin_p","h_mom_coin_p;Momentum [MeV/c];Cointime [ns]",1000,1730.,1930.,2000,-20.,20.);
  TH2F* h_coin_seg  = new TH2F("h_coin_seg","#LS2seg=5: fixed;Cointime [ns];#Segment",2000,-20.,20.,16,1.,17.);
  TH2F* h_mom_seg  = new TH2F("h_mom_seg","Mom. vs #Seg;Momentum [GeV/c];#Segment",100,1730.,1930.,16,1.,17.);
  //h_coin_seg->GetXaxis()->SetTitle("Cointime [ns]");
  //h_coin_seg->GetYaxis()->SetTitle("#Segment");
  TH2F* h_coin_seg2  = new TH2F("h_coin_seg2","#LS2seg=#RS2seg;Cointime [ns];#Segment",2000,-20.,20.,16,1.,17.);
  //h_coin_seg2->GetXaxis()->SetTitle("Cointime [ns]");
  //h_coin_seg2->GetYaxis()->SetTitle("#Segment");
  TH1F* hcoin_test = new TH1F("hcoin_test","",40000/56,-20.,20.);
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
ct = -ct;//Kyoto ctime was inverse

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

	bool seg_flag = false;
	for(int k=0;k<16;k++){
		if(abs(ct_wseg[k][8])<20.)seg_flag=true;
	}
	seg_flag=true;
if(seg_flag){//&&fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.025){
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
		if(ac1sum>1.&&ac2sum>5.&&abs(ct-3.)<1.)h_mom_coin_pi->Fill(R_mom*1000.,ct);
		if(ac1sum<1.&&ac2sum<5.&&ct>10.*R_mom-27.3&&ct<10.*R_mom-25.3)h_mom_coin_p->Fill(R_mom*1000.,ct);
		for(int j=0;j<16;j++){
			//if(ct_wseg[j][j]>-30.)h_coin_seg->Fill(ct_wseg[j][j],h_coin_seg->GetYaxis()->GetBinCenter(j+1));
			if(ct_wseg[j][j]>-30.)h_coin_seg->Fill(ct_wseg[j][j],h_coin_seg->GetYaxis()->GetBinCenter(j+1));
			if(ct_wseg[8][j]>-30.)h_mom_seg->Fill(R_mom*1000.,h_mom_seg->GetYaxis()->GetBinCenter(j+1));
			//if(ct_wseg[5][j]>-30.)h_coin_seg2->Fill(ct_wseg[5][j],h_coin_seg2->GetYaxis()->GetBinCenter(j+1));
			
			hcoin_test->Fill(ct_wseg[j][j]);
		}
		//	h_coin_seg2->Fill(ct_wseg[LS2_seg][RS2_seg],LS2_seg);
	}//#LS2seg=#RS2seg=5

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
	double mom_from_ct[slice] = {0};
	double mom_from_cte1[slice] = {0};
	double mom_from_cte2[slice] = {0};
	double mom_from_cte3[slice] = {0};
	double mom_from_cte4[slice] = {0};

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
		double a = 27.3;//m; L
		double b = LightVelocity*abs(t_diff[i]);//ct_diff
		double c = Mpi;
		double d = Mp;
		double term1 = a*a*b*b*c*c/(4.*a*a*b*b-pow(b,4.));
		double term2 = a*a*b*b*d*d/(4.*a*a*b*b-pow(b,4.));
		double term3 = 2.*sqrt(pow(a,4.)*b*b*(a*a*pow(c,4.)-2.*pow(a*c*d,2.)+pow(a*d,2.)*d*d+2.*pow(b*c*d,2.)))/(4.*a*a*b*b-pow(b,4.)); 
//		mom_from_ct[i]=sqrt(pow(a*c,4.)-2.*pow(a,4.)*c*c*d*d+pow(a*d,4.)-2.*pow(a*b*c,2.)-2.*pow(a*b*d,2.)+pow(b,4.))/2/a/b;
		//mom_from_ct[i]=sqrt(-term1-term2-term3);
		mom_from_ct[i]=sqrt(-term1-term2+term3);
		double dpath = 0.025;
		double ae1 = 27.3-dpath;//m; L
		double ae2 = 27.3-dpath;//m; L
		double ae3 = 27.3+dpath;//m; L
		double ae4 = 27.3+dpath;//m; L
		double be1 = LightVelocity*abs(t_diff[i]-t_diff_sig[i]);//ct_diff
		double be2 = LightVelocity*abs(t_diff[i]+t_diff_sig[i]);//ct_diff
		double be3 = LightVelocity*abs(t_diff[i]-t_diff_sig[i]);//ct_diff
		double be4 = LightVelocity*abs(t_diff[i]+t_diff_sig[i]);//ct_diff
		double term1e1 = ae1*ae1*be1*be1*c*c/(4.*ae1*ae1*be1*be1-pow(be1,4.));
		double term1e2 = ae2*ae2*be2*be2*c*c/(4.*ae2*ae2*be2*be2-pow(be2,4.));
		double term1e3 = ae3*ae3*be3*be3*c*c/(4.*ae3*ae3*be3*be3-pow(be3,4.));
		double term1e4 = ae4*ae4*be4*be4*c*c/(4.*ae4*ae4*be4*be4-pow(be4,4.));
		double term2e1 = ae1*ae1*be1*be1*d*d/(4.*ae1*ae1*be1*be1-pow(be1,4.));
		double term2e2 = ae2*ae2*be2*be2*d*d/(4.*ae2*ae2*be2*be2-pow(be2,4.));
		double term2e3 = ae3*ae3*be3*be3*d*d/(4.*ae3*ae3*be3*be3-pow(be3,4.));
		double term2e4 = ae4*ae4*be4*be4*d*d/(4.*ae4*ae4*be4*be4-pow(be4,4.));
		double term3e1 = 2.*sqrt(pow(ae1,4.)*be1*be1*(ae1*ae1*pow(c,4.)-2.*pow(ae1*c*d,2.)+pow(ae1*d,2.)*d*d+2.*pow(be1*c*d,2.)))/(4.*ae1*ae1*be1*be1-pow(be1,4.)); 
		double term3e2 = 2.*sqrt(pow(ae2,4.)*be2*be2*(ae2*ae2*pow(c,4.)-2.*pow(ae2*c*d,2.)+pow(ae2*d,2.)*d*d+2.*pow(be2*c*d,2.)))/(4.*ae2*ae2*be2*be2-pow(be2,4.)); 
		double term3e3 = 2.*sqrt(pow(ae3,4.)*be3*be3*(ae3*ae3*pow(c,4.)-2.*pow(ae3*c*d,2.)+pow(ae3*d,2.)*d*d+2.*pow(be3*c*d,2.)))/(4.*ae3*ae3*be3*be3-pow(be3,4.)); 
		double term3e4 = 2.*sqrt(pow(ae4,4.)*be4*be4*(ae4*ae4*pow(c,4.)-2.*pow(ae4*c*d,2.)+pow(ae4*d,2.)*d*d+2.*pow(be4*c*d,2.)))/(4.*ae4*ae4*be4*be4-pow(be4,4.)); 
		mom_from_cte1[i]=sqrt(-term1e1-term2e1+term3e1);
		mom_from_cte2[i]=sqrt(-term1e2-term2e2+term3e2);
		mom_from_cte3[i]=sqrt(-term1e3-term2e3+term3e3);
		mom_from_cte4[i]=sqrt(-term1e4-term2e4+term3e4);
	cout<<"mom_from_ct["<<i<<"]="<<mom_from_ct[i]<<" [GeV/c]"<<endl;
		h_path_mom->SetBinContent(i+1,pathlen[i]);
		h_path_mom->SetBinError(i+1,pathlen_sig[i]);
		h_tdiff_mom->SetBinContent(i+1,abs(t_diff[i]));
		h_tdiff_mom->SetBinError(i+1,t_diff_sig[i]);
		h_pionpos_mom->SetBinContent(i+1,pion_par1[i]);
		h_pionpos_mom->SetBinError(i+1,pion_parerr1[i]);
		h_protonpos_mom->SetBinContent(i+1,proton_par1[i]);
		h_protonpos_mom->SetBinError(i+1,proton_parerr1[i]);
		h_mom_from_ct->SetBinContent(i+1,mom_from_ct[i]*1000.);
		h_mom_from_cte1->SetBinContent(i+1,mom_from_cte1[i]*1000.);
		h_mom_from_cte2->SetBinContent(i+1,mom_from_cte2[i]*1000.);
		h_mom_from_cte3->SetBinContent(i+1,mom_from_cte3[i]*1000.);
		h_mom_from_cte4->SetBinContent(i+1,mom_from_cte4[i]*1000.);
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
	TF1* fpi_line=new TF1("fpi_line","pol1",1730.,1930.);
	//TF1* fpi_line=new TF1("fpi_line","pol1",1760.,1900.);
	fpi_line->SetNpx(2000);
	TProfile* ppi = h_mom_coin_pi->ProfileX();
	TProfile* pp = h_mom_coin_p->ProfileX();
	fpi_line->SetParameter(0,3.15);//intercept
	//fpi_line->SetParLimits(0,3.0,3.20);
	fpi_line->SetParameter(1,-0.001);//slope
	TF1* fp_line=new TF1("fp_line","pol1",1730.,1930.);
	//TF1* fp_line=new TF1("fp_line","pol1",1760.,1900.);
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

	TCanvas *c16 = new TCanvas("c16", "c16", 800, 800);
	h_coin_seg->Draw("colz");
	TCanvas *c17 = new TCanvas("c17", "c17", 800, 800);
//	hcoin_test->Draw("");
	h_mom_seg->Draw("colz");
	TCanvas *c18 = new TCanvas("c18", "c18", 800, 800);
	h_mom_from_ct->Draw("");
	h_mom_from_cte1->Draw("same");
	h_mom_from_cte2->Draw("same");
	h_mom_from_cte3->Draw("same");
	h_mom_from_cte4->Draw("same");
	//h_coin_seg2->Draw("colz");
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

}//path_from_ct
