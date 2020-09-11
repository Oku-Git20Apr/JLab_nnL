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
  val = par[num] * step * sum * invsq2pi / (par[num+2]*par[num+1]*exp(-par[num+3]/par[num+1]));
  return val;
}

double FMM_Lambda_Sigma( double *x, double *par , int num)
  {
  double val = par[num] * TMath::Gaus(x[0],par[num+1],par[num+2]);//Lambda Gaussian
  val += par[num+3] * TMath::Gaus(x[0],par[num+4],par[num+5]);//Sigma Gaussian
  return val;
}
double FMM_Res( double *x, double *par ){

	return FMM_Lambda_Sigma(x,par,0)+expgaus2(x,par,6)+expgaus2(x,par,10)+expgaus2(x,par,14);

}

void resp(){

	TCanvas *c1 = new TCanvas("c1","c1",800,800);
 const double ML = 1.115683;            // Lambda       mass (GeV/c2)
 const double MS0 = 1.192642;           // Sigma Zero   mass (GeV/c2)
 const double def_n_L=250.; 
 const double def_sig_L=0.003; 
 const double def_mean_L=0.0;
 const double def_n_S=70.; 
 const double def_sig_S=0.004;
 const double def_mean_S=MS0-ML;
 const double const_L_best=258.305;
 const double mean_L_best=-0.000457363;
 const double sig_L_best=0.00190997;
 const double const_S_best=78.7173;
 const double mean_S_best=0.0766371;
 const double sig_S_best=0.00190055;

 TF1 *f4 = new TF1("f4", FMM_Res, -0.1, 0.2, 18);
	 f4->SetNpx(20000);
	 f4->SetParameter(0,281);
	 f4->SetParameter(1,-0.00064);
	 f4->SetParameter(2,0.00148);
	 f4->SetParameter(3,73.29);
	 f4->SetParameter(4,0.0766);
	 f4->SetParameter(5,0.0015);
	 f4->SetParameter(6,0.0766);
	 f4->SetParameter(7,0.0324);
	 f4->SetParameter(8,0.0235);
	 f4->SetParameter(9,-0.0165);
	 f4->SetParameter(10,20);
	 f4->SetParameter(11,0.03);
	 f4->SetParameter(12,0.00116);
	 f4->SetParameter(13,-0.0698);
	 f4->SetParameter(14,0.7097);
	 f4->SetParameter(15,0.00724);
	 f4->SetParameter(16,0.000175);
	 f4->SetParameter(17,-0.00150);
 f4->Draw("");
}
