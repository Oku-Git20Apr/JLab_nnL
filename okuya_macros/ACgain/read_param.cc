void read_param(){

	string pname = "../ac/param/offset_ac.dat";
	ifstream ifp(pname.c_str(),ios::in);

cout << "Parameter file: " << pname.c_str() << endl;

	string buf;
	int ac, seg;
	double ped, ope;

	while(1){
		getline(ifp,buf);
		if(buf[0]=='#') continue;// skip comment out
		if(ifp.eof()) break;
		if(buf.empty()) {cout << "buf.empty()"<< endl;break;}
	
		stringstream sbuf(buf);
		sbuf >> ac >> seg >> ped >> ope;
			if(ac==1){
//cout << "a1| seg: " << seg << ", pedestal: " << ped << ", 1 p.e.: " << ope << endl;
cout << "a1| seg: " << seg << ", gain: " << ope - ped << endl;

			}
			else if(ac==2){
//cout << "a2| seg: " << seg << ", pedestal: " << ped << ", 1 p.e.: " << ope << endl;
cout << "a2| seg: " << seg << ", gain: " << ope - ped << endl;
			}
			else{
				cout << "Error has occured." << endl; exit(1);
			}
	}

	return 0;

}

