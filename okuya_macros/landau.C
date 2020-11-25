double landaugaus( double *x, double *par ){

   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //par[4]=Scale
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
	return par[4]*val1;
}

/*------Landau * Gauss convluted------*/

void landau(){

	TCanvas *c1 = new TCanvas("c1","c1",800,800);
//c1->cd();
// TF1 *f1 = new TF1("f1", FMM_4Poly_wRes, -0.1, 0.2, 17);
//f1->SetNpx(2000);
//f1->SetParameter(0,242.);
//f1->SetParameter(1,-0.000510);
//f1->SetParameter(2,0.00163);
//f1->SetParameter(3,71.3);
//f1->SetParameter(4,0.0764);
//f1->SetParameter(5,0.00191);
//f1->SetParameter(6,-1.38);
//f1->SetParameter(7,-15.3);
//f1->SetParameter(8,94.5);
//f1->SetParameter(9,244.);
//f1->SetParameter(10,97.8);
//f1->SetParameter(11,73.8);
//f1->SetParameter(12,24.3);
//f1->SetParameter(13,-0.041);
//f1->SetParameter(14,6.3);
//f1->SetParameter(15,0.09);
//f1->SetParameter(16,0.135);
////f1->SetParameter(17,0.0107);
//f1->SetParLimits(0,0.,1000000);//positive

 TF1 *f4 = new TF1("f4", landaugaus, -10, 100, 5);
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //par[4]=Scale
 f4->SetNpx(2000);
 f4->SetParameter(0,0.5);
 f4->SetParameter(1,1.5);
 f4->SetParameter(2,0.5);
 f4->SetParameter(3,0.8);
 f4->SetParameter(4,10.);
 f4->Draw("");
cout<<"Integral="<<f4->Integral(-10.,100.)<<endl;
}
