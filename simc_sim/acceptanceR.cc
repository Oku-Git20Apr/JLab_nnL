//------------------------------------//
//------------------------------------//
//---   HRS-R Acceptance           ---//
//------------------------------------//
//K. Okuyama (Oct. 13, 2020)


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
#include "TLorentzVector.h"
using namespace std;
void SetTitle( TH2D*, const char*, const char*, const char*);
void SetTitle( TH1D*, const char*, const char*, const char*);

static const double PI = 4.*atan(1.0);

int main(int argc, char** argv){
	//////////////////////
	// command argument //
	//////////////////////
	int option;
	//string filename = "1H_kaon";
	//string filename = "RHRS_new";//0<theta_max<0.067
	//string filename = "RHRS_big";//0<theta_max<0.01
	string filename = "RHRS_cm";//0<theta_max<0.01
	//string filename = "RHRS2";//0<theta_max<0.01
	bool BatchFlag = false;
	bool PDFFlag = true;

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
	//thetawidth=0.067;//new
	//thetawidth=0.1;//big
	thetawidth=PI;//CM
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
	double Domega = 2.*PI* (1-cos(thetawidth))*1000.; // [msr]
	double momwidth = 0.5;
	double MAX = 0.00001;
	double volume = Domega*MAX/1000.; 
	cout<< "thetamin = " << thetamin*180/PI << " [deg]" <<endl;
	cout<< "thetamax = " << thetamax*180/PI << " [deg]" <<endl;
	cout<< "phimin = " << phimin*180/PI << " [deg]" <<endl;
	cout<< "phimax = " << phimax*180/PI << " [deg]" <<endl;
	cout<< "Domega = " << Domega << " [msr]" <<endl;
	cout<< "volume = " << volume <<endl;


	int bin_mom = 150;
	double min_mom = 1.6;//nnL
	double max_mom = 2.0;//nnL
	double min_mom_cm = 0.45;//nnL
	double max_mom_cm = 0.85;//nnL
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
	//t->Project("h_mom_result","Rp_orig");

	TH1D *h_mom_gen_cm = new TH1D( "h_mom_gen_cm", "", bin_mom, min_mom_cm, max_mom_cm);
	TH1D *h_mom_result_cm = new TH1D( "h_mom_result_cm", "", bin_mom, min_mom_cm, max_mom_cm);
	TH1D *h_sa_mom_result_cm = new TH1D("h_sa_mom_result_cm", "", bin_mom, min_mom_cm, max_mom_cm);

	TH1D *h_momR = new TH1D( "h_momR", "", 200,1.73,1.93);
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
	float R_mom=0.;
	float R_th=0.;
	float R_ph=0.;
	float L_mom=0.;
	float L_th=0.;
	float L_ph=0.;
	float vertex=0.;
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Rp_gen",1);tree->SetBranchAddress("Rp_gen",&R_mom);
	tree->SetBranchStatus("Rth_gen",1);tree->SetBranchAddress("Rth_gen",&R_th);
    tree->SetBranchStatus("Rph_gen"    ,1);tree->SetBranchAddress("Rph_gen"    ,&R_ph     );
	tree->SetBranchStatus("Lp_gen",1);tree->SetBranchAddress("Lp_gen",&L_mom);
	tree->SetBranchStatus("Lth_gen",1);tree->SetBranchAddress("Lth_gen",&L_th);
    tree->SetBranchStatus("Lph_gen"    ,1);tree->SetBranchAddress("Lph_gen"    ,&L_ph     );
    tree->SetBranchStatus("zposi"    ,1);tree->SetBranchAddress("zposi"    ,&vertex     );
    ENum = tree->GetEntries();
  for(int i=0;i<ENum;i++){
    tree->GetEntry(i);
    if(i%100000==0)cout<<i<<" / "<<ENum<<endl;
	
	L_mom /= 1000.;//MeV-->GeV
	R_mom /= 1000.;//MeV-->GeV
	double Einc  = 4.318;//[GeV]
	//double Escat = 2.1;//[GeV]
	double Escat = L_mom;//[GeV]
	double theta = L_th+centraltheta;	

	double Me=511.*pow(10.,-6.);//[GeV/c^2]
	double Mp=0.9382720;//[GeV/c^2]
	double MK=0.494677;//[GeV/c^2]
    double ML = 1.115683;//[GeV/c2]
    double mh = ML;//hypernuclei
    double mt = Mp;//target mass
	double Qsq=2*Einc*Escat*(1-cos(theta));
	double deltaE = Einc - Escat;
	double q2=Qsq+deltaE*deltaE;
	double kg=deltaE-Qsq/(2*Mp);
	double eps=1/(1+2*(q2/Qsq)*tan(theta/2)*tan(theta/2));
		//h_mom_gen->SetBinContent(h_mom_gen->FindBin(mom),deltaE);
		//h_mom_result->SetBinContent(h_mom_gen->FindBin(mom),mom);
		for(int i=0;i<bin_mom;i++){
		//h_mom_gen->SetBinContent(i,26147*(1./189.)*(max_mom-min_mom)/bin_mom);//first try
		//h_mom_gen->SetBinContent(i,2547932*(1./163.8)*1000.*(max_mom-min_mom)/bin_mom);//RHRS_new 1,000,000 (2020/10/18)
		//h_mom_gen->SetBinContent(i,5682429*(1./163.8)*1000.*(max_mom-min_mom)/bin_mom);//RHRS2 1,000,000 (2020/11/2)// true density
		h_mom_gen->SetBinContent(i,168647157*(1./163.8)*1000.*(max_mom-min_mom)/bin_mom);//RHRS_cm 1,000,000 (2020/11/21)// true density
		//h_mom_gen->SetBinContent(i,5710184*(1./163.8)*1000.*(max_mom-min_mom)/bin_mom);//RHRS_big 1,000,000 (2020/10/27)// true density
		//h_mom_gen_cm->SetBinContent(i,2547932*(4./25.)*1000.*(1./163.8)*(max_mom_cm-min_mom_cm)/bin_mom);//RHRS_new 1,000,000 (2020/10/18)
		h_mom_gen_cm->SetBinContent(i,168647157*1000.*(1./163.8)*(max_mom_cm-min_mom_cm)/bin_mom);//RHRS_cm 1,000,000 (2020/11/21)
		//h_mom_gen->SetBinContent(i,2549979*(1./163.8)*1000.*(max_mom-min_mom)/bin_mom);//RHRS 1,000,000 (2020/10/4), w/ z dep.
		//h_mom_gen_cm->SetBinContent(i,2549979*(4./25.)*1000.*(1./163.8)*(max_mom_cm-min_mom_cm)/bin_mom);//RHRS 1,000,000 (2020/10/4), w/ z dep.
		//h_mom_gen->SetBinContent(i,2549979*(1./163.8)*1000.*(max_mom-min_mom)/bin_mom);//RHRS 1,000,000 (2020/10/4)
		//h_mom_gen_cm->SetBinContent(i,2549979*(1./163.8)*1000.*(max_mom_cm-min_mom_cm)/bin_mom);//RHRS 1,000,000 (2020/10/4)
		//h_mom_gen->SetBinContent(i,5248172*(1./163.8)*1000.*(max_mom-min_mom)/bin_mom);//BOTH 1,000,000 (2020/10/4)
		//h_mom_gen_cm->SetBinContent(i,5248172*(1./163.8)*1000.*(max_mom_cm-min_mom_cm)/bin_mom);//BOTH 1,000,000 (2020/10/4)
		}
	double vpflux=Escat*kg/(137*2*PI*PI*Einc*Qsq*(1-eps));
	double k = MAX*gRandom->Uniform();
		//h_mom_result->SetBinContent(h_mom_gen->FindBin(mom),mom);
//cout<<"Qsq="<<Qsq<<endl;
//cout<<"q2="<<q2<<endl;
//cout<<"kg="<<kg<<endl;
//cout<<"eps="<<eps<<endl;

	
	    //===== Right Hand Coordinate ====//
	    //th and phi are originally meant tan(theta) and tan(phi),
	    //so, they should not be treated like tan(R_tr_tr_th) //2020.6.30 Okuyama

		double B_mom = Einc;
		
	    double R_pz = R_mom/sqrt(1.0*1.0 + pow(R_th, 2.0) + pow( R_ph,2.0));
	    double R_px = R_pz * ( R_th );
	    double R_py = R_pz * ( R_ph );

	    double L_pz = L_mom/sqrt(1.0*1.0 + pow(L_th, 2.0) + pow(L_ph,2.0));
	    double L_px = L_pz * ( L_th );
	    double L_py = L_pz * ( L_ph );


	    double B_E =sqrt(B_mom*B_mom + Me*Me);
	    double R_E =sqrt(R_mom*R_mom + MK*MK);
	    double L_E =sqrt(L_mom*L_mom + Me*Me);


		TLorentzVector L_4vec;//Left
		TLorentzVector R_4vec;//Right
		TLorentzVector B_4vec;//Beam
		TLorentzVector T_4vec;//Target
		TLorentzVector G_4vec;//Gamma (Virtual Photon)
		L_4vec.SetPxPyPzE(L_px, L_py, L_pz, L_E);
        R_4vec.SetPxPyPzE(R_px, R_py, R_pz, R_E);
        B_4vec.SetPxPyPzE(0.0 ,  0.0,B_mom, B_E);
        T_4vec.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);




	    double pL    = L_mom;//GeV
	    double pR    = R_mom;//GeV
		theta = L_th;
		double theta_R = R_th;
		double phi = L_ph;
		double phi_R = R_ph;
		double phi0=13.2*PI/180;//rad
		double phi_L = L_4vec.Phi();//LHRS frame
		double phi_RHRS = R_4vec.Phi();//RHRS frame
	    L_4vec.RotateX( -13.2/180.*PI );
	    R_4vec.RotateX(  13.2/180.*PI );
        double mass,mm;
		TLorentzVector Missing;
		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass = Missing.M();
        //mass = sqrt( (Ee + mt - L_E - R_E)*(Ee + mt - L_E - R_E)-(B_v - L_v - R_v)*(B_v - L_v - R_v) );
	    mm=mass - mh;//shift by ML

		double theta_ee = L_4vec.Theta();
		//double theta_ee = acos((-phi*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
		double theta_ek = R_4vec.Theta();
		//double theta_ek = acos((phi_R*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
		double phi_ee = L_4vec.Phi();//original frame
//cout<<"phi_ee="<<phi_ee<<endl;
		double phi_ek = R_4vec.Phi()+2*PI;//original frame
//cout<<"phi_ek="<<phi_ek<<endl;

		G_4vec = B_4vec - L_4vec;
		double mom_g=sqrt(G_4vec.Px()*G_4vec.Px()+G_4vec.Py()*G_4vec.Py()+G_4vec.Pz()*G_4vec.Pz());
		Qsq = G_4vec.M()*G_4vec.M();
		double phi_g = G_4vec.Phi()+2*PI;
		double theta_g = G_4vec.Theta();
//cout<<"theta_gk="<<(theta_ek-theta_g)*180./PI<<endl;
		double theta_gk_lab = G_4vec.Angle(R_4vec.Vect());
//cout<<"theta_gk_lab(TLorentz)="<<theta_gk_lab*180./PI<<endl;
		double omega=G_4vec.E();
		double pY = sqrt((omega+Mp-sqrt(pR*pR+MK*MK))*(omega+Mp-sqrt(pR*pR+MK*MK))-ML*ML);
		double W = sqrt((omega+Mp)*(omega+Mp)-mom_g*mom_g);
		//double theta_gk_lab_test = acos((mom_g*mom_g+pR*pR-pY*pY)/(2.*pR*mom_g));
		double theta_eg_lab = acos((B_E-L_E*(1.-Qsq/2./B_E/L_E))/mom_g);
//cout<<"theta_eg_lab="<<theta_eg_lab*180./PI<<endl;
		double theta_gk_lab_test = 13.2*PI/180.-theta_eg_lab;
//cout<<"theta_gk_lab="<<theta_gk_lab_test*180./PI<<endl;
		double beta=mom_g/(omega+Mp);
		double ER=sqrt(pR*pR+MK*MK);
		double gamma=1./sqrt(1-beta*beta);
		double theta_gk_cm_test = atan((pR*sin(theta_gk_lab_test))/(-1.*gamma*beta*ER+gamma*pR*cos(theta_gk_lab_test)));
//cout<<"theta_gk_cm="<<theta_gk_cm_test*180./PI<<endl;
	
		TVector3 boost;
		TLorentzVector GT_4vec;
		GT_4vec=G_4vec+T_4vec;
		boost=GT_4vec.BoostVector();
		R_4vec.Boost(-boost);
		L_4vec.Boost(-boost);
		B_4vec.Boost(-boost);
		double theta_gk_cm = G_4vec.Angle(R_4vec.Vect());
		double pR_cm=sqrt(R_4vec.Px()*R_4vec.Px()+R_4vec.Py()*R_4vec.Py()+R_4vec.Pz()*R_4vec.Pz());
		double pL_cm=sqrt(L_4vec.Px()*L_4vec.Px()+L_4vec.Py()*L_4vec.Py()+L_4vec.Pz()*L_4vec.Pz());
		double pB_cm=sqrt(B_4vec.Px()*B_4vec.Px()+B_4vec.Py()*B_4vec.Py()+B_4vec.Pz()*B_4vec.Pz());

//cout<<"theta_gk_cm(TLorentz)="<<theta_gk_cm*180./PI<<endl;
		double p_cm=sqrt(GT_4vec.Px()*GT_4vec.Px()+GT_4vec.Py()*GT_4vec.Py()+GT_4vec.Pz()*GT_4vec.Pz());
		double E_cm = GT_4vec.E();
		double ER_cm=sqrt(pR_cm*pR_cm+MK*MK);
//cout<<"beta="<<beta<<endl;
//cout<<"gamma="<<gamma<<endl;

		double labtocm = (gamma*pR_cm*pR_cm*(pR_cm*cos(theta_gk_cm)+beta*ER_cm))/(pow(sqrt(pR_cm*pR_cm*sin(theta_gk_cm)*sin(theta_gk_cm)+gamma*gamma*(pR_cm*cos(theta_gk_cm)+beta*ER_cm)*(pR_cm*cos(theta_gk_cm)+beta*ER_cm)),3.));
//cout<<"labtocm="<<labtocm<<endl;
		double tan_lab1 = sin(theta_gk_cm)/(gamma*(cos(theta_gk_cm)+beta*sqrt(MK*MK+pR_cm*pR_cm)/pR_cm));
		double tan_lab2 = sin(theta_gk_cm)/(gamma*(cos(theta_gk_cm)+(omega*Mp-Qsq*Qsq)/(omega*Mp+Mp*Mp)));

		h_mom_result_cm->Fill(pR_cm);
		
		//Z partition; SIMC_RHRS_z1(,z2,z3,z4,z5).dat
		//if(6.<=vertex&&vertex<10.)h_mom_result->Fill(R_mom);
		h_mom_result->Fill(R_mom);


		if(L_mom>2.&&L_mom<2.21)h_momR->Fill(R_mom);
	}

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
		if(n1!=0 && n2!=0)val = Domega*(1.0*n2/n1);
		h_sa_mom_result->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_result->SetBinError(i+1, err);

		// === all cuts in CM frame ===
		n1 = 0;
		n1 = (int)h_mom_gen_cm->GetBinContent(i+1);
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_result_cm->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = Domega*(1.0*n2/n1);
		h_sa_mom_result_cm->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_result_cm->SetBinError(i+1, err);
	}

////////////////////
// Draw histgrams //
////////////////////
	TH1D *hframe;
	TLine *lmom = new TLine( centralmom, 0., centralmom, Domega);
	TLine *lth = new TLine( centraltheta, 0., centraltheta, Domega);
	lmom->SetLineColor(4);
	lth->SetLineColor(4);

	TCanvas *c1 = new TCanvas("c1", "Momentum Acceptance (Lab)");
	TCanvas *c2 = new TCanvas("c2", "Angular Acceptance (CM)");


	//======= Momentum Acceptance (Lab) ======
	c1->Divide(2,2);
	c1->cd(1);
//	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, Domega + 0.5);
	//SetTitle(h_sa_mom_result, "Solid Angle vs. Momentum at Reference Plane", "Momentum [GeV/c]", "Solid Angle [msr]");
	SetTitle(h_sa_mom_result, "Solid Angle vs. Momentum (w/ all Cuts)", "Momentum [GeV/c]", "Solid Angle [msr]");
//	h_sa_mom_result->GetXaxis()->SetNdivisions(506, kFALSE);
	h_sa_mom_result->GetXaxis()->SetNdivisions(505, kFALSE);
	h_sa_mom_result->SetMarkerColor(kSpring-1);
	h_sa_mom_result->SetLineColor(kSpring-1);
	h_sa_mom_result->SetLineWidth(2);
//	h_sa_mom_result->Draw("same");
	h_sa_mom_result->Draw("e2");
	h_sa_mom_result->Draw("p same");
//	lmom->Draw("same");
	c1->cd(2);
	SetTitle(h_mom_gen, "Solid Angle vs. Momentum (generated)", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_mom_gen->GetXaxis()->SetNdivisions(505, kFALSE);
	h_mom_gen->Draw();
	c1->cd(3);
	SetTitle(h_mom_result, "Solid Angle vs. Momentum (cut)", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_mom_result->GetXaxis()->SetNdivisions(505, kFALSE);
	h_mom_result->Draw();
	c1->cd(4);
	h_mom_gen->Draw();
	h_mom_result->SetLineColor(kAzure);
	h_mom_result->Draw("same");

	//======= Momentum Acceptance (CM) ======
	c2->Divide(2,2);
	c2->cd(1);
//	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen_cm, 0., max_mom_gen_cm, Domega + 0.5);
	//SetTitle(h_sa_mom_result_cm, "Solid Angle vs. Momentum at Reference Plane (CM)", "Momentum [GeV/c]", "Solid Angle [msr]");
	SetTitle(h_sa_mom_result_cm, "Solid Angle vs. Momentum in CM (w/ all Cuts)", "Momentum [GeV/c]", "Solid Angle [msr]");
//	h_sa_mom_result_cm->GetXaxis()->SetNdivisions(506, kFALSE);
	h_sa_mom_result_cm->GetXaxis()->SetNdivisions(505, kFALSE);
//	h_sa_mom_result_cm->Draw("same");
	h_sa_mom_result_cm->Draw();
//	lmom->Draw("same");
	c2->cd(2);
	SetTitle(h_mom_gen_cm, "Solid Angle vs. Momentum in CM (generated)", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_mom_gen_cm->GetXaxis()->SetNdivisions(505, kFALSE);
	h_mom_gen_cm->Draw();
	c2->cd(3);
	SetTitle(h_mom_result_cm, "Solid Angle vs. Momentum in CM (cut)", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_mom_result_cm->GetXaxis()->SetNdivisions(505, kFALSE);
	h_mom_result_cm->Draw();
	c2->cd(4);
	h_mom_gen_cm->Draw();
	h_mom_result_cm->SetLineColor(kAzure);
	h_mom_result_cm->Draw("same");

	TCanvas *c3 = new TCanvas("c3", "c3");
	h_momR->Draw("");

	if(!BatchFlag){
		theApp.Run();
	}

	if(PDFFlag){
		cout << "Creating a PDF file ... " << endl;
		c1 ->Print(Form("%s[",PDFFile.str().c_str()) );
		c1 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c2 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c2 ->Print(Form("%s]" ,PDFFile.str().c_str()) );
	}
c1->Close();
c2->Close();
	
	

	ofstream fout2("RHRS_SIMC.dat");//RHRS
		fout2<<"#Acceptance by SIMC"<<endl;
		fout2<<"#1.8<pe[GeV/c]<2.4, 150 partition --> 1bin=4MeV/c"<<endl;
		double Acceptance_temp;
		double Acceptance_total=0.;
		int Acceptance_bin=0;
	for(int i=0; i<bin_mom; i++){
		Acceptance_temp = h_sa_mom_result->GetBinContent(i+1);
		Acceptance_temp = Acceptance_temp;
		fout2<<i+1<<" "<<Acceptance_temp<<endl;
		Acceptance_total+=Acceptance_temp;
		if(Acceptance_temp!=0)Acceptance_bin++;
	}
	cout<<"Acceptance (average)="<<Acceptance_total/(double)Acceptance_bin<<endl;


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
