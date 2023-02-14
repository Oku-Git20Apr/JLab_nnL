#include <TRandom.h>

bool RHRS=false;
const int nfoil = 10;
double fcent[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
		       0.00, 0.025, 0.05, 0.10, 0.125};
double fcent_real[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
			      0.000, 0.025, 0.050, 0.100, 0.125};
double selection_width = 0.01; // event selection width for z
const double l0 = 1.03; // [m]
//const double l0 = 103; // [cm]
double dth[nfoil];
const double step = 0.492 * 2.54;
const int nrow = 11; // the number of row in SS pattern
const int ncol = 8;  // the number of column in SS pattern
const int nsshole = nrow*ncol; // the number of holes to consider
const double hrs_ang = 13.2/180.*3.14;
double l[nfoil];
double refx[nsshole];
double refy[nsshole];
double ssx_cent_real[nsshole];
double ssy_cent_real[nsshole];
double conv_inch_to_cm =2.54;



void SieveCalc(){


  TRandom random;
  
  for(int i=0;i<nfoil;i++){
    l[i] = 0;
    l[i]=(l0-fcent_real[i]/cos(hrs_ang))*100.;
    dth[i] = atan(l0*sin(hrs_ang)/(l0*cos(hrs_ang) -fcent_real[i])); }
  

  double theta,phi,vz;
  double ssx,ssy;
  int zfoil;
  double theta_res = 0.006;
  double phi_res   = 0.006;
  double vz_res    = 0.006; // [m]
  double hole_size = 0.157; // cm
  double hole_res  = 0.4; // cm
  
  TH1D* hssx = new TH1D("hssx","",1000,-10,10);
  TH1D* hssy = new TH1D("hssy","",1000,-10,10);
  TH1D* hth = new TH1D("hth","",1000,-0.1,0.1);
  TH1D* hph = new TH1D("hph","",1000,-0.1,0.1);  
  TH2D* hss  = new TH2D("hss","Sieve Sliet; ssx [cm] ; ssy [cm]",1000,-10,10,1000,-10,10);
  hssx->SetLineColor(2);
  hssy->SetLineColor(4);
  
  int  ngen =0;
  int  nhole=0;
  // generate hole posi -> angle

  
  for(int i=0;i<ncol;i++){ //  hole postion ssy
    for(int j=0;j<nrow;j++){ // hole positon ssx
      ngen =0;
      ssy_cent_real[i] = -3.0*step + step*i;
      if(j%2==0)ssy_cent_real[i] = ssy_cent_real[i] - step/2.0;
      ssx_cent_real[j] = 5.0*step - step*j;
      refx[nhole] = ssx_cent_real[j];
      refy[nhole] = ssy_cent_real[i];
      
      for(int h=0; h<nfoil;h++){ // select z-foil
	
	while(ngen<1.0e5){
	  
	  theta = random.Gaus(0.0,0.001);
	  phi   = random.Gaus(0.0,0.001);
	  vz    = random.Gaus(fcent_real[h],vz_res);
	  
	  for(int k=0 ; k<nfoil ; k++){
	    if(fcent[k]-selection_width<vz  && vz<fcent[k]+selection_width){
	      zfoil=k;
	      
	      if(!RHRS ) ssy=l[k]*sin(atan(-phi))/cos(dth[k]-atan(-phi));
	      if(RHRS  ) ssy=l[k]*sin(atan(-phi))/cos(dth[k]+atan(-phi));	     
	      double lx;
	      if(ssy>0)lx=sqrt(pow(l[k],2.0) + pow(ssy,2.0) + 2.0*l[k]*ssy*sin(dth[k]));
	      else     lx=sqrt(pow(l[k],2.0) + pow(ssy,2.0) - 2.0*l[k]*ssy*sin(dth[k]));
	      ssx = - theta * lx;
	      
	      hssx->Fill(ssx);
	      hssy->Fill(ssy);
	      hss->Fill(ssx,ssy);
	      hth->Fill(theta);
	      hph->Fill(phi);
	      ngen++;	  
	      
	    }
	  }
	}// end ngen
      }
      nhole++;
    } // end nrow
  } // end ncol
  
  
  TCanvas* c0 =new TCanvas("c0","c0");
  c0->Divide(2,1);
  c0->cd(1);
  hssx->Draw();
  c0->cd(2);
  hssy->Draw("");

  
  TCanvas* c1 =new TCanvas("c1","c1");
  c1->Divide(2,1);
  c1->cd(1);
  hth->Draw();
  c1->cd(2);
  hph->Draw();

  TCanvas* c2 =new TCanvas("c2","c2");
  c2->cd();
  hss->Draw("colz");


  string   rname = "test.root";
  TFile* f =new TFile(rname.c_str(),"recreate");
  hssx->Write();
  hssy->Write();
  hth->Write();
  hph->Write();
  hss->Write();
}
