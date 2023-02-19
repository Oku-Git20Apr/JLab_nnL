//-------------------------------------------//
//--  Statistical Error					   --//
//--	Clopper-Peason's confidence level  --//
//-------------------------------------------//
//
//K. Okuyama (Oct. 6, 2020)
//K. Okuyama (Dec. 13, 2020)//Al Contamination
//
void pearson(){ 

	//string error_in = "./input/error_input.dat";//center value of error
	string error_in = "./input/error_input_2023.dat";//updated 2023.1.26
	//string ofname = "./input/error_out.pdf";

	string buf;
	int nofdata = 100;
	int npoint = 0;
	int npoint2 = 0;
	int npoint3 = 0;
	int flag = -1;
	double xval, xerr, yval, yerr;
	double x[nofdata], y[nofdata], xe[nofdata], ye[nofdata];
	double x2[nofdata], y2[nofdata], xe2[nofdata], ye2[nofdata];
	double x3[nofdata], y3[nofdata], xe3[nofdata], ye3[nofdata];



	ifstream ifp(error_in.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << error_in.c_str() << endl;
	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> flag >> xval >> yval;
		cout << flag << ", " << xval << ", " << yval <<endl;
		
		if(flag==1){//nominator
			x[npoint] = xval; 
			y[npoint] = yval;
			xe[npoint]= 0.; 
			ye[npoint]= 0.; 
			npoint++;
		}
		if(flag==2){//denominator
			x2[npoint2] = xval; 
			y2[npoint2] = yval;
			xe2[npoint2]= 0.; 
			ye2[npoint2]= 0.; 
			npoint2++;
		}

	}

	TH1F* he1_n = new TH1F("he1_n", "he1_n", 20, 0., 20.);//nominator
	TH1F* he1_d = new TH1F("he1_d", "he1_d", 20, 0., 20.);//denominator
	for(int i=0;i<npoint;i++){
	he1_n->SetBinContent(x[i]+1,y[i]);
	he1_d->SetBinContent(x2[i]+1,y2[i]);
	}
	he1_n->SetLineColor(kRed);
	he1_d->SetLineColor(kAzure);

cout << "TEfficiency!" << endl;
cout<<"pEff1:"<<endl;
TEfficiency *pEff1;
if(TEfficiency::CheckConsistency(*he1_n,*he1_d,"w")){
pEff1 = new TEfficiency(*he1_n,*he1_d);
}



cout << "Drawing Start " << endl;
TCanvas* c1 = new TCanvas("c1","c1",800.,800.);
he1_d->Draw("");
he1_n->Draw("same");
TCanvas* c2 = new TCanvas("c2","c2",800.,800.);
pEff1->Draw("");

cout<<"Lambda Z (best)"<<endl;
cout<<pEff1->GetEfficiency(2)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(2)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(2)<<endl;
cout<<"Lambda AC (best)"<<endl;
cout<<pEff1->GetEfficiency(3)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(3)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(3)<<endl;
cout<<"Lambda CT (best)"<<endl;
cout<<pEff1->GetEfficiency(4)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(4)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(4)<<endl;
cout<<"Lambda AC (strict)"<<endl;
cout<<pEff1->GetEfficiency(6)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(6)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(6)<<endl;
cout<<"Lambda CT (strict)"<<endl;
cout<<pEff1->GetEfficiency(7)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(7)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(7)<<endl;

cout<<"Sigma0 Z (best)"<<endl;
cout<<pEff1->GetEfficiency(8)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(8)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(8)<<endl;
cout<<"Sigma0 AC (best)"<<endl;
cout<<pEff1->GetEfficiency(9)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(9)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(9)<<endl;
cout<<"Sigma0 CT (best)"<<endl;
cout<<pEff1->GetEfficiency(10)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(10)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(10)<<endl;
cout<<"Sigma0 AC (strict)"<<endl;
cout<<pEff1->GetEfficiency(11)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(11)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(11)<<endl;
cout<<"Sigma0 CT (strict)"<<endl;
cout<<pEff1->GetEfficiency(12)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(12)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(12)<<endl;

cout<<"Al Contamination"<<endl;
cout<<pEff1->GetEfficiency(13)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(13)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(13)<<endl;

cout<<"TEST"<<endl;
cout<<pEff1->GetEfficiency(14)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(14)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(14)<<endl;
//  cout<<"Print is starting"<<endl;
//  cout<<"pdf name: "<<ofname<<endl;
//	 c1->Print(Form("%s[",ofname.c_str()));
//	 c1->Print(Form("%s",ofname.c_str()));
//	 c2->Print(Form("%s",ofname.c_str()));
//	 c2->Print(Form("%s]",ofname.c_str()));

cout<<"Well Done!"<<endl;

}
