double FMM_4Poly_wRes( double *x, double *par )
  {
    double Napier = 2.7182818;
    double val = par[0] * TMath::Gaus(x[0],par[1],par[2]);//Lambda Gaussian
    val += par[3] * TMath::Gaus(x[0],par[4],par[5]);//Sigma Gaussian
    val +=-TMath::Floor(x[0]-0.007)*TMath::Ceil(x[0]-0.002)*( par[6] + par[7]*x[0] + par[8]*x[0]*x[0] + par[9]*x[0]*x[0]*x[0] + par[10]*TMath::Power(x[0],4.));//4Poly BG
    //val +=TMath::Ceil(x[0]-par[1]-par[2])* par[11] * TMath::Power(Napier,-par[12]*(x[0]-par[13]));//*TMath::Gaus(x[0],par[13],par[16]);//Lambda Radiative tail
    //val +=TMath::Ceil(x[0]-par[1]-par[2])* par[16] * TMath::Power(Napier,-par[17]*(x[0]-par[18]));//Simga Radiative tail
    //if(x[0]-par[4]>0.)val += par[14] * TMath::Power(Napier,-par[12]*(x[0]-par[15]))*TMath::Gaus(x[0],par[15],par[17]);//Sim     ga Radiative tail//relevant to Lambda Radiative tail
    val += TMath::Ceil(x[0]-par[4]-par[5])* par[14] * TMath::Power(Napier,-par[12]*(x[0]-par[15]));//*TMath::Gaus(x[0],par[13],par[16]);//Sim     ga Radiative tail//relevant to Lambda Radiative tail
    //val += par[11] * TMath::Power(Napier,-par[12]*(x[0]-par[13]));//Lambda Radiative tail
    //val += par[14] * TMath::Power(Napier,-par[15]*(x[0]-par[16]));//Simga Radiative tail
    return val;
}//&&(x[0]-par[4]+2*par[5])<0

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
double FMM_Res_test( double *x, double *par ){

	return FMM_Lambda_Sigma(x,par,0)+expgaus2(x,par,6)+expgaus2(x,par,10)+expgaus2(x,par,14);

}
double sqrt2pi = sqrt(2.*3.141592);
double F_Pois( double *x, double *par )
{
  /*
    par[0] : lambda, average of #photon
    par[1] : energy resolution factor
    par[2] : menseki
  */
  double val = 0;
  for( int np=1; np<200; np++ ){
    double pois;
    double sigma = par[1]*sqrt(np);
    if(np<50){
      pois = par[2]*pow( par[0], np )*exp(-par[0])/TMath::Gamma(np+1);
    }
    else{ // stirling's approximation
      pois = par[2]*pow( par[0]/np, np )*exp(-par[0]+np)/(sqrt2pi*pow(np,0.5));
    }
    val += pois/(sqrt2pi*sigma)*exp( -pow(x[0]-np,2)/2./sigma/sigma ); // adding gauss    ian distribution
  }
  return val;
}

double F_res( double *x, double *par )
{
    double Napier = 2.7182818;
	double val=0.;
    val += par[0] * TMath::Gaus(x[0],par[1],par[2]);//Lambda Gaussian
    val += par[3] * TMath::Gaus(x[0],par[4],par[5]);//Sigma Gaussian
   // if((x[0]-par[1]-par[2]>0.)&&(x[0]-par[1]-4*par[2]<0.))
    val += (par[6]+par[7]*x[0]+par[8]*x[0]*x[0]); 
    val +=TMath::Ceil(x[0]-par[1]-par[2])* par[9] * TMath::Power(Napier,-par[11]*(x[0]-par[13]));//Simga Radiative tail
    //if(x[0]-par[4]>0.)val += par[14] * TMath::Power(Napier,-par[12]*(x[0]-par[15]))*TMath::Gaus(x[0],par[15],par[17]);//Sim     ga Radiative tail//relevant to Lambda Radiative tail
    val += TMath::Ceil(x[0]-par[4]-par[5])* par[10] * TMath::Power(Napier,-par[12]*(x[0]-par[14]));//*TMath::Gaus(x[0],par[13],par[16]);//Sim     ga Radiative tail//relevant to Lambda Radiative tail
    //val += par[11] * TMath::Power(Napier,-par[12]*(x[0]-par[13]));//Lambda Radiative tail
    //val += par[14] * TMath::Power(Napier,-par[15]*(x[0]-par[16]));//Simga Radiative tail
    return val;
}
double F_Lambda( double *x, double *par )
{
    double Napier = 2.7182818;
	double val=0.;
    val += par[0] * TMath::Gaus(x[0],par[1],par[2]);//Lambda Gaussian
    val += par[3] * TMath::Ceil(x[0]-par[1])*TMath::Power(Napier,-par[4]*(x[0]-par[5]))*TMath::Gaus(x[0],par[6],par[7]);//Lambda Radiative tail
   // val += par[8] * TMath::Power(Napier,-par[9]*(x[0]-par[10]))*TMath::Gaus(x[0],par[11],par[12]);//Lambda Radiative tail
    return val;
}
void res(){

 TF1 *f1 = new TF1("f1", FMM_4Poly_wRes, -0.1, 0.2, 19);
f1->SetNpx(2000);
f1->SetParameter(0,259.);
f1->SetParameter(1,-0.000510);
f1->SetParameter(2,0.0012);
f1->SetParameter(3,45.3);
f1->SetParameter(4,0.077);
f1->SetParameter(5,0.00231);
f1->SetParameter(6,50.);
f1->SetParameter(7,-5000.);
f1->SetParameter(8,10.);
f1->SetParameter(9,0.);
f1->SetParameter(10,0.);
f1->SetParameter(11,22.9);
f1->SetParameter(12,16.);
f1->SetParameter(13,0.01);
f1->SetParameter(14,20.);
f1->SetParameter(15,0.03);
f1->SetParameter(16,20.);
f1->SetParameter(17,63.);
f1->SetParameter(18,0.0015);
f1->SetParLimits(0,0.,1000000);//positive

 TF1 *f2 = new TF1("f2", F_res, 0.005, 0.01, 15);
f2->SetNpx(2000);
//f2->SetParameter(0,259.);
f2->SetParameter(0,0.);
f2->SetParameter(1,-0.000510);
f2->SetParameter(2,0.0012);
//f2->SetParameter(3,45.3);
f2->SetParameter(3,0.);
f2->SetParameter(4,0.077);
f2->SetParameter(5,0.00231);
f2->SetParameter(6,50.);//poly
f2->SetParameter(7,-5000);//poly
f2->SetParameter(8,10.);//poly
//f2->SetParameter(9,20.);
f2->SetParameter(9,0.);
//f2->SetParameter(10,10.);
f2->SetParameter(10,0.);
//f2->SetParameter(11,16.);
f2->SetParameter(11,0.);
//f2->SetParameter(12,16.);
f2->SetParameter(12,0.);
//f2->SetParameter(13,0.);
f2->SetParameter(13,0.);
//f2->SetParameter(14,0.077);
f2->SetParameter(14,0.);
 //TF1 *f2 = new TF1("f2", "TMath::Ceil(x)", -0.1,0.2);
 
 TF1 *f3 = new TF1("f3", F_Lambda, -0.1, 0.2, 13);
 f3->SetNpx(2000);
 f3->SetParameter(0,259.);
 f3->SetParameter(1,-0.000510);
 f3->SetParameter(2,0.0012);
 f3->SetParameter(3,22.9);
 f3->SetParameter(4,16.);
 f3->SetParameter(5,0.01);
 f3->SetParameter(6,0.0012);
 f3->SetParameter(7,0.0052);
 f3->SetParameter(8,20.);
 f3->SetParameter(9,2.);
 f3->SetParameter(10,0.0012);
 f3->SetParameter(11,0.0012);
 f3->SetParameter(12,0.1);
// f3->Draw("");

 //TF1 *f4 = new TF1("func4", expgaus2, -0.1, 0.2,4);
 //f4->SetNpx(2000);
 //f4->SetParameter(0,1000.);//scale
 //f4->SetParameter(1,50.);//att.
 //f4->SetParameter(2,0.003);//sigma
 //f4->SetParameter(3,0.);//peak pos.
 //f4->Draw("");
// TF1 *f0 = new TF1("func0",FMM_Res_test, -0.1, 0.2,6);
// f0->SetNpx(2000);
 //f0->SetParameter(0,259.);
 //f0->SetParameter(1,-0.000510);
 //f0->SetParameter(2,0.0012);
 //f0->SetParameter(3,45.3);
 //f0->SetParameter(4,0.077);
 //f0->SetParameter(5,0.00231);
//f4->Draw();
 TF1 *f5 = new TF1("f5",FMM_Res_test, -0.1, 0.2, 18);
 f5->SetParameter(0,259.);
 //f5->SetParameter(0,0.);
 f5->SetParameter(1,-0.000510);
 f5->SetParameter(2,0.0012);
 f5->SetParameter(3,45.3);
 //f5->SetParameter(3,0.);
 f5->SetParameter(4,0.077);
 f5->SetParameter(5,0.00231);
 f5->SetParameter(6,0.3);//scale
 f5->SetParameter(7,0.02);//att.
 f5->SetParameter(8,0.002);//sigma
 f5->SetParameter(9,0.0);//peak pos.
 f5->SetParameter(10,10.);
 f5->SetParameter(11,0.04);
 f5->SetParameter(12,0.002);
 f5->SetParameter(13,-0.077);
 f5->SetParameter(14,0.3);
 f5->SetParameter(15,0.04);
 f5->SetParameter(16,0.002);
 f5->SetParameter(17,0.0);
// f5->SetNpx(2000);
// f5->SetParameter(6,10.);//scale
// f5->SetParameter(7,5000.);//att.
// f5->SetParameter(8,0.003);//sigma
// f5->SetParameter(9,0.);//peak pos.
 f5->Draw("");
}
