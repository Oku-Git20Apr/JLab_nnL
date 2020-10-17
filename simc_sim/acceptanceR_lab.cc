//------------------------------------//
//------------------------------------//
//---   HRS-R Acceptance           ---//
//------------------------------------//



#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <time.h>
#include <string>

#include "TStyle.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TCut.h"
#include "TLine.h"
#include "TRandom.h"
using namespace std;
void SetTitle( TH2D*, const char*, const char*, const char*);
void SetTitle( TH1D*, const char*, const char*, const char*);

static const double PI = 4.*atan(1.0);

double vpflux_lab(double *par){
	double Einc  = 4.3;//[GeV]
	//double Escat = 2.1;//[GeV]
	double Escat = par[0];//[GeV]
	double theta = par[1];	

	double Me=pow(511,-6.);//[GeV/c^2]
	double Mp=0.9382720;//[GeV/c^2]
	double Qsq=2*Einc*Escat*(1-cos(theta));
	double omega = Einc - Escat;
	double q2=Qsq+omega*omega;
	double kg=omega-Qsq/(2*Mp);
	double eps=1/(1+2*(q2/Qsq)*tan(theta/2)*tan(theta/2));
//if(theta>0.230&&theta<0.231){
//cout<<"Qsq="<<Qsq<<endl;
//cout<<"q2="<<q2<<endl;
//cout<<"kg="<<kg<<endl;
//cout<<"eps="<<eps<<endl;
//}

	double vpflux=Escat*kg/(137*2*PI*PI*Einc*Qsq*(1-eps));
	return vpflux;
}
//TF1* func3 = new TF1("func3",vpflux_lab, 0.001, 0.3,1);
//func3->SetNpx(600);
//func3->SetParameter(0,2.1);
//func3->SetLineColor(kGreen);
//func3->SetLineWidth(4);
//func3->Draw("same");

int main(int argc, char** argv){
	//////////////////////
	// command argument //
	//////////////////////
	int option;
	//string filename = "1H_kaon";
	string filename = "RHRS";
	bool BatchFlag = true;
	bool PDFFlag = true;
	while((option=getopt(argc, argv, "f:bp"))!=-1){
		switch(option){
			case 'f':
				filename = optarg;
				break;
			case 'b':
				cout << "batch mode" << endl;
				BatchFlag = true;
				break;
			case 'p':
				PDFFlag = true;
				break;
			case 'h':
				cout<<"-f : input root filename"<<endl;
				cout<<"-b : execute in batch mode"<<endl;
				cout<<"-p : print pdf file"<<endl;
				return 0;
				break;
			case '?':
				cout << "unknown option: " << option << endl;
				return 0;
				break;
			default:
				cout<<"type -h to see help!!"<<endl;
				return 0;
				break;
		}
	}
  TApplication theApp("App", &argc, argv);
	ostringstream RootFile;
	ostringstream PDFFile;
	RootFile << "../" << filename << ".root";	
	PDFFile << "./vpflux.pdf";	

	///////////////////////
	// General condition //
	///////////////////////
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleSize(0.05, "X");
  gStyle->SetTitleSize(0.05, "Y");
  gStyle->SetTitleOffset(0.9, "X");
  gStyle->SetTitleOffset(0.9, "Y");
	gStyle->SetOptStat(0);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);

	FILE *fp;
	const int SizeOfBuffer = 32;
    char str[SizeOfBuffer];
	double centralmom, centraltheta, thetawidth;
	double centralphi, phiwidth;


	
	centraltheta=13.2*PI/180.;
	thetawidth=0.07;
	centralphi=0.;
	phiwidth=2*PI;
	cout<< "central theta = " << centraltheta*180/PI << " [deg]" << endl;
	cout<< "theta width= " << thetawidth*180/PI << " [deg]" << endl;
	centralphi = 0.;

	TFile *f = new TFile(RootFile.str().c_str());
	TTree *t = (TTree*)f->Get("SNT");
	TTree *tree = (TTree*)f->Get("SNT");

	double thetamin = 1.*(1.*centraltheta - thetawidth);
	double thetamax = 1.*(1.*centraltheta + thetawidth);
	double phimin = 1.*(centralphi - phiwidth);
	double phimax = 1.*(centralphi + phiwidth);
	double omega = 2.*PI* (1-cos(thetawidth))*1000.; // [msr]
	double momwidth = 0.5;
	double MAX = 0.00001;
	double volume = omega*MAX/1000.; 
	cout<< "thetamin = " << thetamin*180/PI << " [deg]" <<endl;
	cout<< "thetamax = " << thetamax*180/PI << " [deg]" <<endl;
	cout<< "phimin = " << phimin*180/PI << " [deg]" <<endl;
	cout<< "phimax = " << phimax*180/PI << " [deg]" <<endl;
	cout<< "omega = " << omega << " [msr]" <<endl;
	cout<< "volume = " << volume <<endl;


	int bin_mom = 150;//expanded?
	double min_mom = 1600;//nnL
	double max_mom = 2000;//nnL
	//double min_mom = 2.6;//K40
	//double max_mom = 3.4;//K40
	int bin_th = 150;
	//double min_th = 0.10;
	//double max_th = 0.35;
	double min_th = 0.95;//cos
	double max_th = 1.0;//cos
	double min_ph = -0.4292;//default=0.57
	double max_ph = 0.4292;
	int bin_2D_mom = 50;
	int bin_2D_th = 50;
	int bin_2D_ph = 50;

	TH1D *h_mom_gen = new TH1D( "h_mom_gen", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_result = new TH1D( "h_mom_result", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_result = new TH1D("h_sa_mom_result", "", bin_mom, min_mom, max_mom);
	t->Project("h_mom_result","Rp_orig");

	TH1D *h_vp_mom = new TH1D( "h_vp_mom", "", bin_mom, min_mom, max_mom);
	TH1D *h_vp_mom2 = new TH1D( "h_vp_mom2", "w/ HRS-R Acceptance", bin_mom, min_mom, max_mom);
	TH1D *h_vp_mom_result = new TH1D("h_vp_mom_result", "", bin_mom, min_mom, max_mom);
	TH1D *h_vp_mom_result2 = new TH1D("h_vp_mom_result2", "", bin_mom, min_mom, max_mom);//w/ RHRS Acceptance

	int n1, n2;
	double val;
	double err;
	const double min_mom_gen = h_mom_gen->GetBinCenter(h_mom_gen->FindFirstBinAbove()) - 0.05;
	const double max_mom_gen = h_mom_gen->GetBinCenter(h_mom_gen->FindLastBinAbove()) + 0.05;
	const double min_mom_bin = h_mom_gen->FindBin(min_mom_gen);
	const double max_mom_bin = h_mom_gen->FindBin(max_mom_gen);
	cout <<"mom range: "<< min_mom_gen << " " << max_mom_gen << endl;
	cout <<"mom_bin: "<< min_mom_bin << " " << max_mom_bin << endl;
//----------------------------//
//------------Fill------------//
//----------------------------//
	int ENum=0;
	float mom=0.;
	float th=0.;
	float ph=0.;
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Rp_gen",1);tree->SetBranchAddress("Rp_gen",&mom);
	tree->SetBranchStatus("Rth_gen",1);tree->SetBranchAddress("Rth_gen",&th);
    tree->SetBranchStatus("Rph_gen"    ,1);tree->SetBranchAddress("Rph_gen"    ,&ph     );
    ENum = tree->GetEntries();
  for(int i=0;i<ENum;i++){
    tree->GetEntry(i);
    if(i%100000==0)cout<<i<<" / "<<ENum<<endl;
	double cosine=cos(th);
	
	double Einc  = 4300.;//[MeV]
	//double Escat = 2.1;//[MeV]
	double Escat = mom;//[MeV]
	double theta = th+centraltheta;	

	double Me=pow(511,-3.);//[MeV/c^2]
	double Mp=938.2720;//[MeV/c^2]
	double Qsq=2*Einc*Escat*(1-cos(theta));
	double deltaE = Einc - Escat;
	double q2=Qsq+deltaE*deltaE;
	double kg=deltaE-Qsq/(2*Mp);
	double eps=1/(1+2*(q2/Qsq)*tan(theta/2)*tan(theta/2));
		//h_mom_gen->SetBinContent(h_mom_gen->FindBin(mom),deltaE);
		//h_mom_result->SetBinContent(h_mom_gen->FindBin(mom),mom);
		for(int i=0;i<bin_mom;i++){
		//h_mom_gen->SetBinContent(i,26147*(1./189.)*(max_mom-min_mom)/bin_mom);//first try
		h_mom_gen->SetBinContent(i,2549979*(1./189.)*(max_mom-min_mom)/bin_mom);//RHRS 1,000,000 (2020/10/4)
		}
	double vpflux=Escat*kg/(137*2*PI*PI*Einc*Qsq*(1-eps));
	double k = MAX*gRandom->Uniform();
		//h_mom_result->SetBinContent(h_mom_gen->FindBin(mom),mom);
//cout<<"Qsq="<<Qsq<<endl;
//cout<<"q2="<<q2<<endl;
//cout<<"kg="<<kg<<endl;
//cout<<"eps="<<eps<<endl;

	 if(vpflux>k){//Full
		h_vp_mom->Fill(mom);
		if(mom>2100)h_vp_mom2->Fill(mom);
	}
	//cout<<"vpflux vs k = "<<vpflux<<" : "<<k<<endl;
		//if(fabs(th-0.225)<0.0125&&fabs(ph-0.25)<0.125){
	 //if(vpflux>k&&fabs(th-0.225)<0.025&&fabs(ph)<0.25){//new top-quality (6msr)
	//}
	 //}
	}

	double vpflux_tot_mom=0.;
	double vpflux_tot_mom_err=0.;
	for(int i=0; i<bin_mom; i++){
		n1 = 0;
		n1 = (int)h_mom_gen->GetBinContent(i+1);

		// === all cuts ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_result->GetBinContent(i+1);
		//cout<<"nbin_mom:"<<i<<endl;
		//cout<<"n1="<<n1<<endl;
		//cout<<"n2="<<n2<<endl;
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_result->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_result->SetBinError(i+1, err);

		// === VP Flux ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_vp_mom->GetBinContent(i+1);
		//cout<<"nbin_mom:"<<i<<endl;
		//cout<<"n1="<<n1<<endl;
		//cout<<"n2="<<n2<<endl;
		if(n1!=0 && n2!=0)val = volume*(1.0*n2/n1)*1000.;//MeV->GeV
		h_vp_mom_result->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_vp_mom_result->SetBinError(i+1, err);

		// === VP Flux (w/ RHRS Acceptance)===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_vp_mom2->GetBinContent(i+1);
		//cout<<"nbin_mom:"<<i<<endl;
		//cout<<"n1="<<n1<<endl;
		//cout<<"n2="<<n2<<endl;
		if(n1!=0 && n2!=0)val = volume*(1.0*n2/n1)*1000;//MeV->GeV
		h_vp_mom_result2->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_vp_mom_result2->SetBinError(i+1, err);

		double bin_width=(max_mom-min_mom)/bin_mom/1000.;//GeV/c
		vpflux_tot_mom+=val*bin_width;
		vpflux_tot_mom_err+=sqrt(err*bin_width*err*bin_width);
	}
	cout<<"vpflux_tot_mom="<<vpflux_tot_mom<<endl;
	cout<<"vpflux_tot_mom_err="<<vpflux_tot_mom_err<<endl;

////////////////////
// Draw histgrams //
////////////////////
	TH1D *hframe;
	TLine *lmom = new TLine( centralmom, 0., centralmom, omega);
	TLine *lth = new TLine( centraltheta, 0., centraltheta, omega);
	lmom->SetLineColor(4);
	lth->SetLineColor(4);

	TCanvas *c1 = new TCanvas("c1", "Momentum Acceptance (Exit)");
	TCanvas *c2 = new TCanvas("c2", "Angular Acceptance (Exit)");
	TCanvas *c3 = new TCanvas("c3", "Momentum Acceptance Expansion");
	TCanvas *c4 = new TCanvas("c4", "Angular Acceptance Expansion");
	TCanvas *c5 = new TCanvas("c5", "Momentum vs. Theta");
	TCanvas *c6 = new TCanvas("c6", "Theta vs. Phi");
	TCanvas *c7 = new TCanvas("c7", "VP Flux");


	//======= Momentum Acceptance Expansion ======
	c3->Divide(2,2);
	c3->cd(1);
//	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	//SetTitle(h_sa_mom_result, "Solid Angle vs. Momentum at Reference Plane", "Momentum [GeV/c]", "Solid Angle [msr]");
	SetTitle(h_sa_mom_result, "Solid Angle vs. Momentum (w/ all Cuts)", "Momentum [MeV/c]", "Solid Angle [msr]");
//	h_sa_mom_result->GetXaxis()->SetNdivisions(506, kFALSE);
	h_sa_mom_result->GetXaxis()->SetNdivisions(505, kFALSE);
//	h_sa_mom_result->Draw("same");
	h_sa_mom_result->Draw();
//	lmom->Draw("same");
	c3->cd(2);
	SetTitle(h_mom_gen, "Solid Angle vs. Momentum (generated)", "Momentum [MeV/c]", "Solid Angle [msr]");
	h_mom_gen->GetXaxis()->SetNdivisions(505, kFALSE);
	h_mom_gen->Draw();
	c3->cd(3);
	SetTitle(h_mom_result, "Solid Angle vs. Momentum (cut)", "Momentum [MeV/c]", "Solid Angle [msr]");
	h_mom_result->GetXaxis()->SetNdivisions(505, kFALSE);
	h_mom_result->Draw();
	c3->cd(4);
	h_mom_gen->Draw();
	h_mom_result->SetLineColor(kAzure);
	h_mom_result->Draw("same");

	//======= VP Flux  ======
	c7->Divide(2,2);
	c7->cd(1);
	//SetTitle(h_vp_mom_result, "Integrated VP Flux vs. Momentum (w/ all Cuts)", "Momentum [GeV/c]", "Integrated VP Flux [/GeV]");
		h_vp_mom_result->Scale(bin_mom);
	h_vp_mom_result->GetXaxis()->SetNdivisions(505, kFALSE);
	h_vp_mom_result->Draw();
	//	h_vp_mom_result2->Scale(bin_mom);
	h_vp_mom_result2->GetXaxis()->SetNdivisions(505, kFALSE);
	h_vp_mom_result2->SetLineColor(kAzure);
	h_vp_mom_result2->Draw("same");
	double ymax = (h_vp_mom_result->GetBinContent(h_vp_mom_result->GetMaximumBin()));
	TLine *RHRS_min = new TLine( 2100., 0., 2100, 1.1*ymax);
	TLine *RHRS_max = new TLine( 2220, 0., 2220, 1.1*ymax);
	RHRS_min->SetLineColor(kRed);
	RHRS_max->SetLineColor(kRed);
	RHRS_min->Draw("same");
	RHRS_max->Draw("same");
	c7->cd(3);
	h_vp_mom->Draw();

	if(!BatchFlag){
		theApp.Run();
	}

	if(PDFFlag){
		cout << "Creating a PDF file ... " << endl;
		c1 ->Print(Form("%s[",PDFFile.str().c_str()) );
		c1 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c2 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c3 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c4 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c5 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c6 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c7 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c7 ->Print(Form("%s]" ,PDFFile.str().c_str()) );
	}
c1->Close();
c2->Close();
c3->Close();
c4->Close();
c5->Close();
c6->Close();
c7->Close();
	
	

	ofstream fout2("RHRS_SIMC.dat");//RHRS
		fout2<<"#Acceptance by SIMC"<<endl;
		fout2<<"#1.8<pe[GeV/c]<2.4, 150 partition --> 1bin=4MeV/c"<<endl;
		double Acceptance_temp;
		double Acceptance_total=0.;
	for(int i=0; i<bin_mom; i++){
		Acceptance_temp = h_sa_mom_result->GetBinContent(i+1);
		Acceptance_temp = Acceptance_temp;
		fout2<<i+1<<" "<<Acceptance_temp<<endl;
		Acceptance_total+=Acceptance_temp;
	}
	cout<<"Acceptance (average)="<<Acceptance_total/150.<<endl;


	cout<< "Finish!" << endl;
	return 0;
}

void SetTitle( TH2D *h, const char *title, const char *xaxis, const char *yaxis){
	h->GetXaxis()->SetTitle(xaxis);
	h->GetYaxis()->SetTitle(yaxis);
	h->SetTitle(title);
}

void SetTitle( TH1D *h, const char *title, const char *xaxis, const char *yaxis){
	h->GetXaxis()->SetTitle(xaxis);
	h->GetYaxis()->SetTitle(yaxis);
	h->SetTitle(title);
}
