//--  Beam Current vs #event  --//
//
//K. Okuyama (Jan. 7, 2024)
//
//This is taken over from ebeam.C
void SetTH1(TH1 *h, TString name, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h->SetTitle(name);
  h->SetLineColor(LColor);
  h->SetLineWidth(1);
  h->SetFillStyle(FStyle);
  h->SetFillColor(FColor);

  h->SetTitleFont(42,"");
  h->SetTitleSize(0.04,"");

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetNdivisions(505);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}
void SetTH2(TH2 *h, TString name, TString xname, TString yname, double min=0.8){
  h->SetTitle(name);
  h->SetMinimum(min);
  h->SetLineWidth(0);
  h->SetTitleSize(0.05,"");
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.005);
  h->SetMarkerColor(1);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(0.90);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetNoExponent();

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetDecimals(3);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

void current_monitor(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

	unsigned long long int nevent=1;
	int flag = 0;
	double I_e; 
	double time;
	double time_at_last;
	int strange_Ie_flag=0;
	

	TH1F* h  = new TH1F("h","Beam Current",1000,0.0,30.0);
	TH2F* hist  = new TH2F("hist","Beam Current vs #event",50000,0.,5000000.,1000,-2.0,30.0);
  	SetTH2(hist, "", "Num. of events", "Ie [#mu A]", 0.4);
	TH2F* hist2  = new TH2F("hist2","Beam Current vs Time",150,33330.,33480.,500,-2.0,30.0);
	//TH2F* hist2  = new TH2F("hist2","Beam Current vs Time",30000,0.,300000.,1000,-2.0,30.0);
  	SetTH2(hist2, "", "Time [s]", "Ie [#mu A]", 0.4);
    hist2->SetMarkerSize(0.2);

	string file_name;
	string files_list = "./current_files.list";//list of current records (after executing ./get_filename.sh)
	string buf;
	ifstream ifp(files_list.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
	while(1){
		strange_Ie_flag=0;
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> file_name;
		cout << file_name <<endl;
		ifstream ifp2(file_name.c_str(),ios::in);
		if (ifp2.fail()){ cout << "Failed" << endl; exit(1);}
		while(1){
		getline(ifp2,buf);
		if(buf[0]=='#'){continue;}
		if(ifp2.eof())break;
		stringstream sbuf2(buf);
		sbuf2 >> time >> I_e;
		time += time_at_last;
		if(I_e<15.)strange_Ie_flag++;
		//cout << I_e <<endl;
		hist->Fill((double)nevent, I_e);
		//if(flag<500)I_e=29.;
		hist2->Fill(time, I_e);
		h->Fill(I_e);
		nevent++;
		flag++;
	}
cout<<"flag="<<flag<<endl;
		flag=0;
		cout<<file_name<<": #strange Ie="<<strange_Ie_flag<<endl;
		time_at_last=time;
	}


	TCanvas* c_temp = new TCanvas("c_temp","c_temp",1800,1000);
	h->Draw("");
	

	TCanvas* c = new TCanvas("c","c",1800,1000);
	//hist->Draw("");
	c->SetLeftMargin(0.14);
	c->SetRightMargin(0.14);
	c->SetTopMargin(0.14);
	c->SetBottomMargin(0.14);
	c->Modified();
	c->Update();
	gPad->Modified();
	gPad->Update();
	//c->Print("./pdf/current_ev.pdf");
	TCanvas* c2 = new TCanvas("c2","c2",1800,1000);
	hist2->Draw("");
	c2->SetLeftMargin(0.14);
	c2->SetRightMargin(0.14);
	c2->SetTopMargin(0.14);
	c2->SetBottomMargin(0.14);
	c2->Modified();
	c2->Update();
	gPad->Modified();
	gPad->Update();
	//c2->Print("./pdf/current_time.pdf");

cout << "Well done!" << endl;
}//fit
