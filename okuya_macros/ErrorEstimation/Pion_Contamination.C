//-------------------------------------------//
//--  Statistical Error					   --//
//--	Clopper-Peason's confidence level  --//
//-------------------------------------------//
//
//K. Okuyama (Oct. 6, 2020)
//K. Okuyama (Dec. 12, 2020)
//
void Pion_Contamination(){ 

	string error_in = "./input/Pion_Contamination.dat";//center value of error
	//string error_in = "./input/error_input.dat";//center value of error
	//string ofname = "./input/error_out.pdf";

	string buf;
	int nofdata = 100;
	int npoint = 0;
	int npoint2 = 0;
	int npoint3 = 0;
	int flag = -1;
	double tru, kaon, pion;
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
		sbuf >> flag >> tru >> kaon >> pion;
		cout << flag << ", " << tru << ", " << kaon << ", " << pion <<endl;
		
			x[npoint] = flag;
			y[npoint] = pion;
			xe[npoint]= 0.; 
			ye[npoint]= 0.; 
			x2[npoint] = flag;
			y2[npoint] = tru;
			xe2[npoint]= 0.; 
			ye2[npoint]= 0.; 
			x3[npoint] = flag;
			y3[npoint] = pion+kaon;
			xe3[npoint]= 0.; 
			ye3[npoint]= 0.; 
			npoint++;

	}

	TH1F* he1_n = new TH1F("he1_n", "he1_n", 20, 0., 20.);//nominator
	TH1F* he1_d = new TH1F("he1_d", "he1_d", 20, 0., 20.);//denominator
	TH1F* he1_d2 = new TH1F("he1_d2", "he1_d2", 20, 0., 20.);//denominator2
	for(int i=0;i<npoint;i++){
	he1_n->SetBinContent(x[i],y[i]);
	he1_d->SetBinContent(x2[i],y2[i]);
	he1_d2->SetBinContent(x3[i],y3[i]);
	}
	he1_n->SetLineColor(kRed);
	he1_d->SetLineColor(kAzure);
	he1_d2->SetLineColor(kGreen);

cout << "TEfficiency!" << endl;
cout<<"pEff1:"<<endl;
TEfficiency *pEff1;
if(TEfficiency::CheckConsistency(*he1_n,*he1_d,"w")){
pEff1 = new TEfficiency(*he1_n,*he1_d);
}
TEfficiency *pEff2;
if(TEfficiency::CheckConsistency(*he1_n,*he1_d2,"w")){
pEff2 = new TEfficiency(*he1_n,*he1_d2);
}



cout << "Drawing Start " << endl;
TCanvas* c1 = new TCanvas("c1","c1",800.,800.);
he1_d->Draw("");
he1_n->Draw("same");
TCanvas* c2 = new TCanvas("c2","c2",800.,800.);
he1_d2->Draw("");
he1_n->Draw("same");
TCanvas* c3 = new TCanvas("c3","c3",800.,800.);
pEff1->Draw("");
pEff2->Draw("same");

cout<<"Method 1: Pion/(All-BG)"<<endl;
cout<<"Double Gauss Fixed"<<endl;
cout<<pEff1->GetEfficiency(1)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(1)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(1)<<endl;
cout<<"Double Gauss Free"<<endl;
cout<<pEff1->GetEfficiency(2)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(2)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(2)<<endl;
cout<<"Voigt Fixed"<<endl;
cout<<pEff1->GetEfficiency(3)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(3)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(3)<<endl;
cout<<"Voigt Free"<<endl;
cout<<pEff1->GetEfficiency(4)<<endl;
cout<<pEff1->GetEfficiencyErrorLow(4)<<endl;
cout<<pEff1->GetEfficiencyErrorUp(4)<<endl;

cout<<"Method 2: Pion/(Pion+Kaon)"<<endl;
cout<<"Double Gauss Fixed"<<endl;
cout<<pEff2->GetEfficiency(1)<<endl;
cout<<pEff2->GetEfficiencyErrorLow(1)<<endl;
cout<<pEff2->GetEfficiencyErrorUp(1)<<endl;
cout<<"Double Gauss Free"<<endl;
cout<<pEff2->GetEfficiency(2)<<endl;
cout<<pEff2->GetEfficiencyErrorLow(2)<<endl;
cout<<pEff2->GetEfficiencyErrorUp(2)<<endl;
cout<<"Voigt Fixed"<<endl;
cout<<pEff2->GetEfficiency(3)<<endl;
cout<<pEff2->GetEfficiencyErrorLow(3)<<endl;
cout<<pEff2->GetEfficiencyErrorUp(3)<<endl;
cout<<"Voigt Free"<<endl;
cout<<pEff2->GetEfficiency(4)<<endl;
cout<<pEff2->GetEfficiencyErrorLow(4)<<endl;
cout<<pEff2->GetEfficiencyErrorUp(4)<<endl;
//  cout<<"Print is starting"<<endl;
//  cout<<"pdf name: "<<ofname<<endl;
//	 c1->Print(Form("%s[",ofname.c_str()));
//	 c1->Print(Form("%s",ofname.c_str()));
//	 c2->Print(Form("%s",ofname.c_str()));
//	 c2->Print(Form("%s]",ofname.c_str()));

cout<<"Well Done!"<<endl;

}
