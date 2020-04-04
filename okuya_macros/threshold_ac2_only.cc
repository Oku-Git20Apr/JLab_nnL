void threshold_ac2_only(){

  string ofname = "./threshold_ac2_only.pdf";
cout << " output file name is " << ofname << endl;
  string pname = "../ac/param/offset_ac.dat";
     ifstream ifp(pname.c_str(),ios::in);
     if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
 cout << "Param file : " << pname.c_str() << endl;
 
     string buf;
     int nofsegment = 24;//a1
     int nofsegment2 = 26;//a2
     int npoint = 0;
     int npoint2 = 0;
     int ac, seg;
     double ped_val, ope_val, gain_val;
     double gain[nofsegment], gain2[nofsegment2], ped[nofsegment], ped2[nofsegment2];

	double Ra1a[100], Ra2a[100], Ra1a_c[100], Ra2a_c[100], Ra1a_p[100], Ra2a_p[100];
	double Ra1sum, Ra2sum;

	double adc1[100], adc2[100];

   TChain *tr = new TChain("T");
   tr->SetBranchStatus("*",0);
   tr->Add("~/rootfile_nnL/tritium_111157_1.root");//or runlist
   tr->Add("~/rootfile_nnL/tritium_111157_2.root");
   tr->Add("~/rootfile_nnL/tritium_111157_3.root");
   tr->Add("~/rootfile_nnL/tritium_111157_4.root");
   tr->Add("~/rootfile_nnL/tritium_111157.root");
   
//	tr->SetBranchStatus("R.a1.a",1);
//	tr->SetBranchAddress("R.a1.a",Ra1a);
//	tr->SetBranchStatus("R.a2.a",1);
//	tr->SetBranchAddress("R.a2.a",Ra2a);
	tr->SetBranchStatus("R.a1.a_c",1);
	tr->SetBranchAddress("R.a1.a_c",Ra1a_c);
	tr->SetBranchStatus("R.a2.a_c",1);
	tr->SetBranchAddress("R.a2.a_c",Ra2a_c);
	tr->SetBranchStatus("R.a1.a_p",1);
	tr->SetBranchAddress("R.a1.a_p",Ra1a);
	tr->SetBranchStatus("R.a2.a_p",1);
	tr->SetBranchAddress("R.a2.a_p",Ra2a);
	tr->SetBranchStatus("R.a1.asum_c",1);
	tr->SetBranchAddress("R.a1.asum_c",&Ra1sum);
	tr->SetBranchStatus("R.a2.asum_c",1);
	tr->SetBranchAddress("R.a2.asum_c",&Ra2sum);
 // tr->SetBranchStatus("R.a1.asum_p",1);
 // tr->SetBranchAddress("R.a1.asum_p",&Ra1sum);
 // tr->SetBranchStatus("R.a2.asum_p",1);
 // tr->SetBranchAddress("R.a2.asum_p",&Ra2sum);
	
	double ENum = tr->GetEntries();
	cout << "ENum = " << ENum << endl;

	TH1D *h_a1[24], *h_a2[26];
	for(int i=0;i<24;i++){
	h_a1[i] = new TH1D(Form("h_a1[%d]",i),Form("h_a1[%d]_myNPE",i),500,-1.,10.);
	}
	for(int i=0;i<26;i++){
	h_a2[i] = new TH1D(Form("h_a2[%d]",i),Form("h_a2[%d]_myNPE",i),1000,-1.,10.);
	}

while(1){
         getline(ifp,buf);
         if(buf[0]=='#'){continue;}
         if(ifp.eof())break;
         stringstream sbuf(buf);
         sbuf >> ac >> seg >> ped_val >> ope_val;
 
         gain_val = ope_val - ped_val;
		if(ac==1){
         ped[seg]  = ped_val;
         gain[seg] = gain_val;
cout << ped[seg] << ":" << gain[seg] << endl;
		}
		if(ac==2){
         ped2[seg]  = ped_val;
         gain2[seg] = gain_val;
		}
     }

	for(int n=0;n<ENum;n++){
		tr->GetEntry(n);
	for(int i=0;i<24;i++){
		double adc1 = Ra1a[i] / gain[i];
cout << adc1 << endl;
		h_a1[i] -> Fill( adc1 );
	}
	for(int i=0;i<26;i++){
		double adc2 = Ra2a[i] / gain2[i];
		h_a2[i] -> Fill( adc2 );
	}
	}
	TCanvas *c1;//AC1 viewer (Linear)
	TCanvas *c2;//AC1 viewer (Log)
	c1 = new TCanvas("c1","A1 Linear",800.,800.);
	c2 = new TCanvas("c2","A1 Log",800.,800.);
	c1->Divide(6,4);
	c2->Divide(6,4);
	for(int i=0;i<24;i++){
	c1->cd(i+1);
	h_a1[i]->Draw();
	h_a1[i]->SetStats(0);
	c2->cd(i+1)->SetLogy();
	h_a1[i]->Draw();
	h_a1[i]->SetStats(0);
  }
 
	TCanvas *c3;//AC2 viewer (Linear)
	TCanvas *c4;//AC2 viewer (Log)
	c3 = new TCanvas("c3","A2 Linear",800.,800.);
	c4 = new TCanvas("c4","A2 Log",800.,800.);
	c3->Divide(7,4);
	c4->Divide(7,4);
	for(int i=0;i<26;i++){
	c3->cd(i+1);
	h_a2[i]->Draw();
	h_a2[i]->SetStats(0);
	c4->cd(i+1)->SetLogy();
	h_a2[i]->Draw();
	h_a2[i]->SetStats(0);
	}

cout << "Print is starting" << endl;
	c1->Print(Form("%s[",ofname.c_str()));
	c1->Print(Form("%s",ofname.c_str()));
	c2->Print(Form("%s",ofname.c_str()));
	c3->Print(Form("%s",ofname.c_str()));
	c4->Print(Form("%s",ofname.c_str()));
	c4->Print(Form("%s]",ofname.c_str()));

	

 return 0;

}


