double FMM_4Poly_wRes( double *x, double *par )
  {
    double Napier = 2.7182818;
    double val = par[0] * TMath::Gaus(x[0],par[1],par[2]);//Lambda Gaussian
    val += par[3] * TMath::Gaus(x[0],par[4],par[5]);//Sigma Gaussian
    val += par[6] + par[7]*x[0] + par[8]*x[0]*x[0] + par[9]*x[0]*x[0]*x[0] + par[10]*TMath::Power(x[0],4.);//4Poly BG
    if(x[0]-par[1]>0.)val += par[11] * TMath::Power(Napier,-par[12]*(x[0]-par[13]))*TMath::Gaus(x[0],par[13],par[16]);//Lambda Radiative tail
    //if((x[0]-par[4]-2*par[5])>0)val += par[14] * TMath::Power(Napier,-par[15]*(x[0]-par[16]));//Simga Radiative tail
    //if(x[0]-par[4]>0.)val += par[14] * TMath::Power(Napier,-par[12]*(x[0]-par[15]))*TMath::Gaus(x[0],par[15],par[17]);//Sim     ga Radiative tail//relevant to Lambda Radiative tail
    if(x[0]-par[4]>0.)val += par[14] * TMath::Power(Napier,-par[12]*(x[0]-par[15]))*TMath::Gaus(x[0],par[13],par[16]);//Sim     ga Radiative tail//relevant to Lambda Radiative tail
    //val += par[11] * TMath::Power(Napier,-par[12]*(x[0]-par[13]));//Lambda Radiative tail
    //val += par[14] * TMath::Power(Napier,-par[15]*(x[0]-par[16]));//Simga Radiative tail
    return val;
}//&&(x[0]-par[4]+2*par[5])<0

//double expgaus2(double *x, double *par) {
//  //par[0]=Total area
//  //par[1]=tau of exp function
//  //par[2]=Width (sigma) of convoluted Gaussian function
//  //par[3]=Shift of Function Peak
//  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
//  double np = 500.0;      // number of convolution steps
//  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
//  double xx, fland, sum = 0.0, xlow, xupp, step, i;
//  double val;
//// Range of convolution integral
//  xlow = 0.;
//  xupp = x[0] + sc * par[2];
//  step = (xupp-xlow) / np;
//// Convolution integral
//  for(i=1.0; i<=np/2; i++){
//     xx = xlow + (i-0.5) * step - par[3];
//     fland = TMath::Gaus(xx,x[0],par[2]);
//     sum += fland * TMath::Exp(-xx/par[1]);
//     xx = xupp - (i-.5) * step - par[3];
//     fland = TMath::Gaus(xx,x[0],par[2]);
//     sum += fland * TMath::Exp(-xx/par[1]);
//  }
//  //val = par[2] * step * sum * invsq2pi / par[3];
//  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));
//  return val;
//}


double expgaus2(double *x,double *par){
//  //par[0]=Total area
//  //par[1]=tau of exp function
//  //par[2]=Width (sigma) of convoluted Gaussian function
//  //par[3]=Shift of Function Peak
// Range of convolution integral
double  sc   = 5.;
double  np   = 500.;
double  sum  = 0.;
double  xlow = 0.;
double  val2 = 0.;
double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
double fland = 0.;
double    xx = 0.;
double     i = 0.;
double  xupp = x[0] + 1.6 * sc * par[2];
double  step = (xupp-xlow) / np;
  if(step<0.)step = 0.;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
  //val = par[2] * step * sum * invsq2pi / par[3];
  val2 =  par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));
  //val2 =  step * sum * invsq2pi / (par[3]*par[4]);
/*------Exp * Gauss convluted------*/
	return val2;
}


void exp(){

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

 TF1 *f4 = new TF1("f4", expgaus2, -1, 10, 4);
 f4->SetNpx(2000);
 f4->SetParameter(0,10.);//scale
 f4->SetParameter(1,0.5);//att.
 f4->SetParameter(2,0.003);//sigma
 f4->SetParameter(3,-0.08);//peak pos.
 f4->Draw("");
cout<<"Integral="<<f4->Integral(-100.,10.)<<endl;
}
