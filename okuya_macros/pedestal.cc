void pedestal(){

  string ofname = "nothing";
cout << " output file name is " << ofname << endl;
  
	double Ra1a[100], Ra2a[100], Ra1a_c[100], Ra2a_c[100], Ra1a_p[100], Ra2a_p[100];
	double Ra1sum, Ra2sum;

	double adc1[100], adc2[100];

    TChain *tr = new TChain("T");
	TCanvas *c = new TCanvas("c","Fitting Result",800.,600.);
   tr->Add("~/rootfile_nnL/tritium_111157_1.root");//or runlist
   tr->Add("~/rootfile_nnL/tritium_111157_2.root");
   tr->Add("~/rootfile_nnL/tritium_111157_3.root");
   tr->Add("~/rootfile_nnL/tritium_111157_4.root");
   tr->Add("~/rootfile_nnL/tritium_111157.root");


	int seg = 2;//a1 segment number

	TH1D *H, *Hp;
	H = new TH1D("H",Form("R.a1.a[%d]",seg),1000,5000.,7000.);
	Hp = new TH1D("Hp",Form("R.a1.a_p[%d]",seg),999,2.,2000.);
	tr->Project("H",Form("R.a1.a[%d]",seg),"","");
	tr->Project("Hp",Form("R.a1.a_p[%d]",seg),"","");
	
	TF1 *f1 = new TF1("f1","gausn",5000.,7000.);
	TF1 *f2 = new TF1("f2","gausn",2.,2000.);

	c->Divide(2,1);
	c->cd(1)->SetLogy(1);
	H->Draw("");
	H->Fit("f1","","",5175.,5250.);
	f1->SetLineStyle(2);
	f1->Draw("same");

	c->cd(2)->SetLogy(1);
	Hp->Draw("");
	Hp->Fit("f2","","",0.,50.);
	f2->SetLineStyle(2);
	f2->Draw("same");

	H->SetStats(0);
	Hp->SetStats(0);
   
//	TCanvas *c1;//AC1 viewer (Linear)
//	TCanvas *c2;//AC1 viewer (Log)
//	c1 = new TCanvas("c1","A1 Linear",800.,800.);
//	c2 = new TCanvas("c2","A1 Log",800.,800.);
//	c1->Divide(6,4);
//	c2->Divide(6,4);
//	for(int i=0;i<24;i++){
//	c1->cd(i+1);
//	h_a1[i]->Draw();
//	h_a1[i]->SetStats(0);
//	c2->cd(i+1)->SetLogy();
//	h_a1[i]->Draw();
//	h_a1[i]->SetStats(0);
//	}
// 
//	TCanvas *c3;//AC2 viewer (Linear)
//	TCanvas *c4;//AC2 viewer (Log)
//	c3 = new TCanvas("c3","A2 Linear",800.,800.);
//	c4 = new TCanvas("c4","A2 Log",800.,800.);
//	c3->Divide(7,4);
//	c4->Divide(7,4);
//	for(int i=0;i<26;i++){
//	c3->cd(i+1);
//	h_a2[i]->Draw();
//	h_a2[i]->SetStats(0);
//	c4->cd(i+1)->SetLogy();
//	h_a2[i]->Draw();
//	h_a2[i]->SetStats(0);
//	}
//
//cout << "Print is starting" << endl;
//	c1->Print(Form("%s[",ofname.c_str()));
//	c1->Print(Form("%s",ofname.c_str()));
//	c2->Print(Form("%s",ofname.c_str()));
//	c3->Print(Form("%s",ofname.c_str()));
//	c4->Print(Form("%s",ofname.c_str()));
//	c4->Print(Form("%s]",ofname.c_str()));
//
	

 return 0;

}


