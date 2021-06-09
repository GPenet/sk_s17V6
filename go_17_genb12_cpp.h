void GEN_BANDES_12::SGUA2::Debug(const char * lib) {
	cout << lib << " gua2 \t" << i_81 << " cols " << col1 + 1 << col2 + 1
		<< " digs " << dig1 + 1 << dig2 + 1 
		<< " digs 0"<<oct<<digs<<dec<< endl;
}

void GEN_BANDES_12::InitialSockets2Setup() {//load permanent data
	for (int i = 0; i < 81; i++) {// initial socket 2
		SGUA2 & w= tsgua2[i];
		w.i_81 = i;
		register int i9 = i/9;
		int tpcol[3][2] = { {1,2},{0,2},{0,1} };
		int rdcol = i9 % 3,dcol=3*(i9/3),*p=tpcol[rdcol];
		w.i9 = i9;
		w.col1 = dcol + p[0];
		w.col2 = dcol + p[1];
		int rd1 = i - 9 * w.i9;//relative d1;d2		
		w.id1 = rd1 / 3 + 3 * w.col1;
		w.id2 = rd1 % 3 + 3 * w.col2;
		if (0) {
			cout << "sua2=" << w.i_81 << " i9=" << w.i9 + 1
				<< " cols12 " << w.col1 + 1 << w.col2 + 1
				<< " id1;id2 " << w.id1 << ";" << w.id2 << endl;
		}
	}
}
void GEN_BANDES_12::SecondSockets2Setup() {
	Build_CheckUAs_Subsets_Table();
	for (i81 = 0; i81 < 81; i81++) {// initial socket 2
		SGUA2 & w = tsgua2[i81];
		w.dig1 = gang27[w.id1];
		w.dig2 = gang27[w.id2];
		w.digs = (1 << w.dig1) | (1 << w.dig2);
		w.valid = 0;
		//skip if no band3 uses it gua2 gua4 gua6 / build guas 
		for (int iband3 = 0; iband3 < nband3; iband3++) {
			if(bands3[iband3].IsGua(i81))// setup guas in band
				w.valid=1;
		}
		if (!w.valid)continue;
		// build revised gangster
		memcpy(w.gangcols, gangb12, sizeof gangb12);
		w.gangcols[w.col1] ^= w.digs;
		w.gangcols[w.col2] ^= w.digs;
		// find guas of the socket
		zh2b_g.InitGangster(gangb12, w.gangcols);
		zh2b5_g.sizef5 = GUALIMSIZE;
		zh2b5_g.modevalid = 0;
		w.nua = 0;
		//w.nua_start = ntua2;
		w.tua= guaw.tua;
		ptua2 = w.tua;
		nua2 = 0;
		//================== GUA collector 2 bands 
		GuaCollect(w.digs );
		for (int i = 0; i < 84; i++) {// find UAs 3 digits
			int fl = floors_3d[i];
			if ((fl & w.digs) == w.digs)	GuaCollect(fl);
		}
		for (int i = 0; i < 126; i++) {// find UAs 4 digits
			int fl = floors_4d[i];
			if ((fl & w.digs) == w.digs) GuaCollect(fl);
		}
		for (int i = 0; i < 126; i++) {// find UAs 5digits
			int fl = 0x1ff ^ floors_4d[i];
			if ((fl & w.digs) == w.digs) GuaCollect(fl);
		}
		SecondSockets2MoreUAs();
		if (nua2 > 30)nua2 = 30;
		guaw.Init(nua2, 0, i81);
		tguas.AddStart(guaw); // can be empty
	}
}

void ZH2B::Init_2digits_banda(BF64  cellsbf) {
	int txcells[50], nxcells = 0;
	BitsInTable64(txcells,nxcells, cellsbf.bf.u64);
	memset(zh2b_g.Digit_cell_Assigned_init, 0, sizeof zh2b_g.Digit_cell_Assigned_init);
	for (int i = 0; i < nxcells; i++) {
		int xcell=txcells[i],cell=From_128_To_81[xcell],
		digit = zh2b_g.puz0[cell];
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_init[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_init[i];
}
//void ZH2B::EndInit_2digits_bandb(int fl, int ibandb) {
//}
/*
extract from bands 1 2 and UAs bands1+2 all uas (uint64_t)
having a chance to be a subset of a gua
Note : this is to save time. A UA having a subset is not an error.
Note 2 : All digits of a subset are in the UA to check  
*/
void GEN_BANDES_12::Build_CheckUAs_Subsets_Table() {
	nua_for_check = 0;
	for (uint32_t i = 0; i < genuasb12.nuab1b2; i++) {// uAs bands
		register uint64_t R = genuasb12.tuab1b2[i]& BIT_SET_2X;
		uint64_t cc = _popcnt64(R);
		if (cc > 10) continue;// nearly no chance to be subset
		tua_for_check[nua_for_check++] = R| (cc<<59);
	}
	for (uint32_t i = 0; i < genuasb12.nua; i++){
		register uint64_t R = genuasb12.tua[i];
		if ((R >> 59) > 10) continue;// nearly no chance to be subset
		tua_for_check[nua_for_check++] = R;
	}
	// and build the digit pattern of the uas in tua_for_check
	for (uint32_t i = 0; i < nua_for_check; i++) {
		uint32_t digs = 0,res;
		register uint64_t R = tua_for_check[i]&BIT_SET_2X;
		while (bitscanforward64(res, R)) {
			R ^= (uint64_t)1 << res;
			int cell = From_128_To_81[res],digit=grid0[cell];
			digs |= 1 << digit;
		}
		uadigs_for_check[i] = digs;
		if (0) {
			cout << Char2Xout(tua_for_check[i]) << " ua for check digs 0"
				<< oct << digs << dec <<" ndigs "<<_popcnt32(digs)<< endl;
		}
	}

}
void GEN_BANDES_12::Build_tuacheck(int fl) {
	nua_check=0;
	int ndigs = 0x1ff ^ fl;
	for (uint32_t i = 0; i < nua_for_check; i++) {
		uint32_t digs = uadigs_for_check[i];
		if (digs & ndigs) continue;
		tuacheck[nua_check++]= tua_for_check[i];
	}
}

int GEN_BANDES_12::Have_tuacheck_subset(  ) {// ua 2x27  + 5 bit length
	uint64_t * t = tuacheck;
	register uint64_t ua2x = ua & BIT_SET_2X;
	for (uint32_t iua = 0; iua < nua_check; iua++) {
		register uint64_t R = t[iua];
		R &= BIT_SET_2X;
		if ((R&ua2x) == R)		return 1;// we have a subset
	}
	return 0;
}

void GEN_BANDES_12::SecondSockets2MoreUAs() {
	diagmore = 0;
	//if (i81 == 29)diagmore = 1;
	// this is for a given GUA2 socket i81
	SGUA2 s = tsgua2[i81];
/*	struct SGUA2 {// 81 possible UA2 sockets
uint64_t tua[SIZETGUA];
int digs, gangcols[9];// revised gangster
uint32_t nua;// nua_start, nua_end;
*/
	int rx0[2], bx0[2], dx3[2];
	rx0[0] = zh2b_g.fd_sols[0][s.dig1].SolRow(s.col2);
	dx3[0] = genb12.grid0[9 * rx0[0] + s.col1];
	rx0[1] = zh2b_g.fd_sols[0][s.dig2].SolRow(s.col1);
	dx3[1] = genb12.grid0[9 * rx0[1] + s.col2];
	bx0[0] = rx0[0] / 3;
	bx0[1] = rx0[1] / 3;

	if (bx0[0] == bx0[1]) return; // want digitsin 2 bands
	if (dx3[0] == dx3[1]) return; // ua 6 cells
	for (int idig = 0; idig < 2; idig++) {
		int dx[2], cx[2],row,dig3;
		if (idig) {
			dx[0] = s.dig1; dx[1] = s.dig2;
			cx[0] = s.col1; cx[1] = s.col2;
			row = rx0[0]; dig3 = dx3[0];
		}
		else {
			dx[0] = s.dig2; dx[1] = s.dig1;
			cx[0] = s.col2; cx[1] = s.col1;
			row = rx0[1]; dig3 = dx3[1];
		}
		band=row/3;
		// dig3 must be possible in 
		if (!((1 << dig3) & s.gangcols[cx[1]])) {
			//cout << " not a valid pattern dig3 col" << endl;
			continue;
		}
		if (diagmore) {
			cout << "SecondSockets2MoreUAs" << endl;
			cout << "digs " << s.dig1 + 1 << s.dig2 + 1 << "\tcols " << s.col1 + 1 << s.col2 + 1 << endl;
			cout << "rows " << rx0[0] + 1 << rx0[1] + 1 << "\tdigs3 " << dx3[0] + 1 << dx3[1] + 1 << endl;
		}
		// now fill band bx[0] except the 2 cells
		zh2b[0].Init_gang();
		BF64 cells_to_fill;
		cells_to_fill.bf.u64=BIT_SET_2X;// all cells
		cells_to_fill.bf.u32[1-band] = 0;// nothing in the other band
		c1 = row * 9 + s.col1;
		c2 = row * 9 + s.col2;
		if (row >= 3) { c1 += 5; c2 += 5; }// must be 2X mode
		cells_to_fill.bf.u64 ^= (uint64_t)1 << c1;
		cells_to_fill.bf.u64 ^= (uint64_t)1 << c2;
		zh2b[0].Init_2digits_banda(cells_to_fill);
		zh2b[0].FullUpdate();
		if(diagmore)zh2b[0].ImageCandidats();
		// at this point,we have only one band with a gangster 
		//not equal to the start and a partial GUA 2cells in the other band
		//we  prepare ZHONE to find partial uas in one band
		uint32_t tua[300]; // supply a tua to ZHONE
		zh1b_g.tua = tua;
		for (int i = 0; i < 9; i++) {// load the gangster pm
			zh1b_g.fd_sols[1][i] = zh2b[0].FD[i].bf.u32[1 - band];
		}
		STD_B416 * myb = &myband1;
		if(!band) myb= &myband2; // to catch the solution
		memcpy(zh1b_g.fd_sols[0], myb->fd_sols[0], sizeof zh1b_g.fd_sols[0]);
		zh1b_g.band0 = myb->band0;
		digp = s.digs | (1 << dig3);// must be in the floor
		//==============  now find more uas 6/7 digits
		ngua6_7 = 6;
		wua0 = ((uint64_t)1 << c1) | ((uint64_t)1 << c2);
		for (int i6 = 0; i6 < 84; i6++) {
			int fl3 = floors_3d[i6];// , fl6 = 0x1ff ^ fl3;
			if (fl3&digp) continue;// digits must be in the multi floors
			floors = 0x1ff ^ fl3;
			if (diagmore)cout << "try 6 floors 0" << oct << floors << dec << endl;
			GuaCollectMore();
		}
		ngua6_7 = 7;
		for (int i7 = 0; i7 < 36; i7++) {
			int fl2 = floors_2d[i7];// , fl7 = 0x1ff ^ fl2;
			if (fl2&digp) continue;// digits must be in the multi floors
			floors = 0x1ff ^ fl2;
			if (diagmore)cout << "try 7 floors 0" << oct << floors << dec << endl;
			GuaCollectMore();
		}
	}
}

void GEN_BANDES_12::GuaCollectMore() {
	//if (diagmore)if (floors == 0663) { diagmore = 2; zh1b_g.diag = 1; }
	//else { diagmore = 1; zh1b_g.diag = 0; }
	zh1b_g.nua = 0;
	zhone[0].Start_nFloors(floors);
	zhone[0].InitGuess();
	if (diagmore > 1) zhone[0].ImageCandidats();
	else zh1b_g.diag = 0;
	if (ngua6_7 == 6)zhone[0].Guess6();
	else zhone[0].Guess7();
	if (zh1b_g.nua) {
		Build_tuacheck(floors);
		SGUA2 s = tsgua2[i81];
		for (uint32_t i = 0; i < zh1b_g.nua; i++) {
			uint32_t w = zh1b_g.tua[i]&BIT_SET_27,	cc=_popcnt32(w),cell;
			//cout << Char27out(w) << endl;;
			if (cc < 8) continue;
			if (cc < 12 && ngua6_7 == 7) continue;
			if (cc > GUALIMSIZE - 2) continue;
			int digpw = digp;// assigned outside the band
			register uint32_t R = w;
			while (bitscanforward(cell, R)) {
				R ^= 1 << cell;
				int dig = zh1b_g.band0[cell]; // digit changed
				digpw |= 1 << dig;
			}
			if (_popcnt32(digpw) != ngua6_7) continue;
			//cout << "valid ua to add to the table" << endl;
			ua = wua0;
			ua |= (uint64_t)w << 32 * (1 - band);
			if (Have_tuacheck_subset()) {
				//cout << Char2Xout(ua) << " ua has subset" << endl;
				continue;
			}
			ua |= (uint64_t)(cc + 2) << 59;
			if (nua2 >= SIZETGUA)nua2 = SIZETGUA - 1; // guess it will be a smaller
			int ir = 		AddUA64(ptua2, nua2, ua);
			if(ir && diagmore)		cout << Char2Xout(ua) << " wua added cc=" << cc + 2 <<" nua2="<<nua2<< endl;

		}
	}
	// try to add the ua to the file
}

int GEN_BANDES_12::Debug17(SGUA2 & w) {// check validity of guas
	STD_B3 & myb3 = bands3[0];
	if( myb3.guas.isguasocket2.Off_c(w.i_81)) return 0;
	//this is a valid pair if not hit in band 3by the known 17 
	//must be hit in bands 1-2 
	if (myb3.guas.ua_pair[w.i_81] & g17b.p17diag.bf.u32[2]) return 0;//hit in band 3
	//cout << Char27out(myb3.guas.ua_pair[w.i_81]) << " i81=" << w.i_81 << endl;

	for (uint32_t i = 0; i < w.nua; i++) {// exit false if not hit 
		register uint64_t R = w.tua[i];
		if (0)
			cout << Char2Xout(R) << " cc=" << (R >> 59) << " i=" << i << endl;
		if (R&g17b.p17diag.bf.u64[0])continue;
		cout << "not valid gua2 for i_81=" << w.i_81 << endl;
		cout << Char2Xout(R) << " cc=" << (R >> 59) <<" i="<<i<< endl;
		return 1;
	}
	return 0; 
	
}


void GEN_BANDES_12::GuaCollect(int fl,int diag) {//use revised gangster
	uint64_t solved_cells = zh2b5_g.FindUAsInit(fl, 1);
	if (!solved_cells) return;// one digit solved true
	if (diag) cout << "return from zh2b5_g.FindUAsInit(fl, 1);" << endl;
	zh2b5_g.CollectUas5();// collect uas for this set of floors
	if (diag) cout << "zh2b5_g.CollectUas5();" << endl;
	if (!zh2b5_g.nuaf5) return;
	Build_tuacheck( fl);
	// check subsets and add to main table
	for (uint32_t i = 0; i < zh2b5_g.nuaf5; i++) {
		ua = zh2b5_g.tuaf5[i].bf.u64;
		//genuasb12.ua = ua;
		if (diag)cout << Char2Xout(ua) << " to add if new2 nua=" << nua2 << endl;
		if (Have_tuacheck_subset()) continue;// superset of a previous ua
		//		if (genuasb12.CheckOld()) continue;// superset of a previous ua
		uint64_t cc = _popcnt64(ua&BIT_SET_2X);
		ua |= cc << 59;
		//if (genuasb12.CheckMain(ua)) continue;// superset of a previous ua
		// see why, seems to be of no effect
		// protect against table limit
		if (nua2 >= SIZETGUA)nua2 = SIZETGUA - 1; // guess it will be a smaller
		int ir = AddUA64(ptua2, nua2,ua);
		if (ir) {
			if (diag)cout << Char2Xout(genuasb12.ua) <<" ir="<<ir << " added nua2="<<nua2 << endl;
		}
	}
}

void GEN_BANDES_12::BuildGang9x3() {
	gang27 = gang[0];
	for (int i = 0; i < 9; i++) {// 9 cols to build out of gangcols
		int istack = C_box[i];
		int * d = gang[i], c = gangcols[i];
		uint32_t bit;
		bitscanforward(bit, c);
		c ^= (1 << bit);
		d[0] = bit;
		gang_digits_cols[bit][istack] = i;
		bitscanforward(bit, c);
		c ^= (1 << bit);
		d[1] = bit;
		gang_digits_cols[bit][istack] = i;
		bitscanforward(bit, c);
		d[2] = bit;
		gang_digits_cols[bit][istack] = i;
	}

}

void GEN_BANDES_12::InitialSockets3Setup() {//load permanent data
	for (int i = 0; i < 81; i++) {// initial socket 3
		SGUA3 & w = tsgua3[i];
		int minir(i / 27);
		w.i_81 = i;
		w.imini = minir;
		w.col1 =3*minir;// minirow first column in gangster
		int dp = i - 27 * minir,// the 9 perms of the gangster minirow
			dg = 9 * minir;// gangster 27 start
		 // pointers to gangster digits
		w.id1 = dp / 9;
		dp -= 9 * w.id1;
		w.id1 += dg;
		dg += 3;// next gangster column inthe minirow
		w.id2 = dp / 3;
		dp -= 3 * w.id2;
		w.id2 += dg;
		dg += 3;// last gangster column in the minirow
		w.id3 = dp +dg;
		if (0) {
			cout << "sua2=" << w.i_81 << " col1=" << w.col1 + 1
				<< " id1;id2,id3 " << w.id1 
				<< ";" << w.id2 << ";" << w.id3 << endl;
		}
	}
}
void GEN_BANDES_12::SecondSockets3Setup() {
	for (int i81 = 0; i81 < 81; i81++) {// initial socket 2
		SGUA3 &w = tsgua3[i81];
		w.dig1 = gang27[w.id1];
		w.dig2 = gang27[w.id2];
		w.dig3 = gang27[w.id3];
		w.digs = 1 << w.dig1 | 1 << w.dig2 | 1 << w.dig3;
		w.valid = 0;

		//skip if no band3 uses it   
		for (int iband3 = 0; iband3 < nband3; iband3++) {
			if (bands3[iband3].IsGua3(i81))// setup guas in band
				w.valid = 1;
		}
		if (!w.valid)continue;
		// Setup the perms for gangsters in minirow
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs=bita|bitb|bitc;
		int triplet_perms[2][3];

		int * p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		w.nua = 0;
		w.tua = guaw.tua;
		nua2 = 0;
		ptua2 = w.tua;

		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			// build revised gangster
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, gangb12, sizeof gangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];

			// find guas of the socket
			zh2b_g.InitGangster(gangb12, rgangcols);
			zh2b5_g.sizef5 = UALIMSIZE-1;//; - 3;
			zh2b5_g.modevalid = 0;
			//w.nua_start = ntua3;

			//================== GUA collector 2 bands 
			GuaCollect(digs);// find UAs 3 digits
			for (int i = 0; i < 126; i++) {// find UAs 4 digits
				int fl =  floors_4d[i];
				if((fl& digs) == digs)	GuaCollect(fl);
			}
			for (int i = 0; i < 126; i++) {// find UAs 5digits
				int fl = 0x1ff ^ floors_4d[i];
				if ((fl & digs) == digs) GuaCollect(fl);
			}
		}
		if (nua2 > 30)nua2 = 30;
		guaw.Init(nua2, 1, i81);
		tguas.AddStart(guaw);// can be empty
	}
	//cout << "endSecondSockets3Setup ntua3=" << ntua3
	//	<< " nactive i81=" << nactive3 << endl;


}

void GEN_BANDES_12::GetStartB2(int ip) {//set  rows 3_9 column 1
	char const *tp[20] = {// 3 out of 6 ordered 
		"012345", "345012", "013245", "245013", "014235", "235014", "015234", "234015",
		"023145", "145023", "024135", "135024", "025134", "134025",
		"034125", "125034", "035124", "124035", "045123", "123045"
	};
		char tpw[7];	strcpy( tpw , tp[ip]);
	for (int i = 0; i < 6; i++)boxd[i] = 0x1ff;
	for (int j = 0, jc = 27; j < 6; j++, jc += 9) {
		int ic = tpw[j] - '0', c0 = tc[ic], bit = 1 << c0;
		grid0[jc] = c0;
		zsol[jc] = c0+'1';
		rowd[j] = 0x1ff ^ bit;
		if (j < 3)boxd[0] ^= bit; else  boxd[3] ^= bit;
	}
	if (sgo.vx[8]) cout << "start b2 ip= " << ip << " "
		<< tp[ip] <<" sol  "
		<< zsol[27] << zsol[36] << zsol[45] << endl;
}
void GEN_BANDES_12::Start(int mode) {
	modeb12 = mode;
	myband1.Initstd();
	zsol[81] = 0;
	nb12 = 0;
}
void GEN_BANDES_12::NewBand1(int iw) {
	go_back = 0;
	i1t16 = iw;
	it16 = tn6_to_416[iw];
	myband1.InitG12(it16);
	memcpy(grid0, myband1.band0, sizeof myband1.band0);
	memcpy(gcheck, myband1.band0, sizeof myband1.band0);

	strcpy(zsol, myband1.band);
	n_auto_b1 = bandminlex.GetAutoMorphs(it16, t_auto_b1);
	for (int i = 0; i < 9; i++) // init columns status
		cold[i] = 0x1ff ^ myband1.gangster[i];
	memcpy(coldf, cold, sizeof coldf);
	zsol[27] = 0;
	g17b.ExpandB1( );// expand and set index

	//myband1.ExpandBand( );// expand and set index

	//tulock.LockExpand(myband1.nmyi3, myband1.nmybv5, myband1.nmybv6);
	//tulock.Store1();
	cout << "i1t16=" << i1t16 << " it16=" << it16 
		<< " n auto morphs=" << n_auto_b1 << endl;
	ntc = 0;
	BitsInTable32(tc, ntc, cold[0]);// first col 6 digits in table
	for (ip20 = 0; ip20 < 20; ip20++) {//0;ip<20 setup initial values for rows columns
		GetStartB2(ip20);
		Find_band2B();
		if (go_back)return;
	}
}
int GEN_BANDES_12::Band2Check() {
	int * zs0 = &gcheck[27];
	n_auto_b1b2 = 0;
	if (n_auto_b1) {
		for (int imorph = 0; imorph < n_auto_b1; imorph++) {
			BANDMINLEX::PERM &p = t_auto_b1[imorph];
			int band[27];// morph the band
			for (int i = 0; i < 9; i++) {
				band[i] = p.map[zs0[p.cols[i]]];
				band[i + 9] = p.map[zs0[p.cols[i] + 9]];
				band[i + 18] = p.map[zs0[p.cols[i] + 18]];
			}
			int ir = G17ComparedOrderedBand(zs0, band);
			if (ir == 1)				return 1;
			else if (!ir) {// auto morph b1 b2 store it for later
				t_auto_b1b2[n_auto_b1b2++] = p;
			}
		}
	}
	n_auto_b2b1 = 0;// possible automorph after perm b1b2
	if (i1t16 == ib2check) {// must try perm bands 12 auto morphs
		int b23[3][9];
		for (int i = 0; i < 3; i++) {// morph band1 to band2 minlex
			register int * rrd = b23[i], *rro = &gcheck[9 * i];
			for (int j = 0; j < 9; j++)
				rrd[j] = pcheck2.map[rro[pcheck2.cols[j]]];
		}
		int ir = G17ComparedOrderedBand(zs0, b23[0]);// is it same as base
		if (ir == 1) 			return 1;
		else if (!ir)// auto morph b1 b2 store it for later
			t_auto_b2b1[n_auto_b2b1++].InitBase(ib2check);
		// must also test all auto morphs b2b1
		for (int imorph = 0; imorph < n_auto_b1; imorph++) {// same automorphs b1 b2
			BANDMINLEX::PERM &pp = t_auto_b1[imorph];
			int b23_a[3][9];
			for (int i = 0; i < 3; i++) {
				register int * rrd = b23_a[i], *rro = b23[i];
				for (int j = 0; j < 9; j++)		rrd[j] = pp.map[rro[pp.cols[j]]];
			}
			int ir = G17ComparedOrderedBand(zs0, b23_a[0]);
			if (ir == 1)return 1;
			else if (!ir)// auto morph b1 b2 store it for later
				t_auto_b2b1[n_auto_b2b1++] = pp;
		}
	}
	return 0;
}
int GEN_BANDES_12::Band3Check() {
	if (i1t16 == ib3check && ib3check == ib2check) {// 3 bands equal use diagonal test 
		BANDMINLEX::PERM * p = minlexusingbands.pout;
		p[0].InitBase(i1t16);
		p[1] = pcheck2;
		p[2] = pcheck3;
		if (minlexusingbands.IsLexMinDirect(gcheck, i1t16, t_auto_b1, n_auto_b1))
			return 1;
		return 0; 
	}
	{
		//========================== morphs on b1b2 base test
		if (n_auto_b1b2) {// still direct automorphism b1b2
			for (int imorph = 0; imorph < n_auto_b1b2; imorph++) {
				BANDMINLEX::PERM &p = t_auto_b1b2[imorph];
				int b23[3][9];
				// direct
				for (int i = 0; i < 3; i++) {// band 3 only
					register int * rrd = b23[i], *rro = &gcheck[54 + 9 * i];
					for (int j = 0; j < 9; j++)		rrd[j] = p.map[rro[p.cols[j]]];
				}
				if (G17ComparedOrderedBand(&gcheck[54], b23[0]) == 1)
					return 1;
			}
		}
		//=========================== perm b1b2 and base test (b1=b2)
		if (n_auto_b2b1) {// possible lower band3 with a perm band1 band2
			int b23[3][9];//first morph to band 2 min lexical
			for (int i = 0; i < 3; i++) {// rows 4 to 9 as of band 2 perm
				register int * rrd = b23[i], *rro = &gcheck[54 + 9 * i];
				for (int j = 0; j < 9; j++)
					rrd[j] = pcheck2.map[rro[pcheck2.cols[j]]];
			}
			for (int imorph = 0; imorph < n_auto_b2b1; imorph++) {// then apply auto morphs 
				BANDMINLEX::PERM &pp = t_auto_b2b1[imorph];
				int b23_a[3][9];
				for (int i = 0; i < 3; i++) {
					register int * rrd = b23_a[i], *rro = b23[i];
					for (int j = 0; j < 9; j++)		rrd[j] = pp.map[rro[pp.cols[j]]];
				}
				if (G17ComparedOrderedBand(&gcheck[54], b23_a[0]) == 1)
					return 1;
			}
		}
		//========================= (b2=b3)#b1  perm b2b3 to consider (direct done)
		if (ib3check == ib2check) {// check b3b2 on  auto morphs b1
			if (gcheck[27] - 1) return 1; // must be '2' in r4c1
			for (int imorph = 0; imorph < n_auto_b1; imorph++) {
				BANDMINLEX::PERM &p = t_auto_b1[imorph];
				int b23[6][9];
				for (int i = 0; i < 6; i++) {// rows 4 to 9 from band 2
					register int * rrd = b23[i], *rro = &gcheck[27 + 9 * i];
					for (int j = 0; j < 9; j++)		rrd[j] = p.map[rro[p.cols[j]]];
				}
				int ir = G17ComparedOrderedBand(&gcheck[27], b23[3]);
				if (ir == 1)return 1;
				if (ir < 1 && G17ComparedOrderedBand(&gcheck[54], b23[0]) == 1)return 1;
			}
		}
		//============================= b1=b3 #b2 
		if (minlexusingbands.IsLexMinDiagB(gcheck, i1t16, ib2check, ib3check, t_auto_b1, n_auto_b1))
			return 1;
	}
	return 0;
}
void GEN_BANDES_12::Find_band2B() {
	int * zs0= &grid0[27];
	register int  *crcb, bit;
	int cd[9], rd[6], bd[6];
	memcpy(rd, rowd, sizeof rd);
	memcpy(cd, coldf, sizeof cd);
	memcpy(bd, boxd, sizeof bd);
	char * zs = zsol;
	// now loop over the 24 cells not yet assigned in the band to fill the band (use relative cell) 
	int ii = -1, free[24];
	uint32_t d;
nextii:
	ii++;
	{	crcb = tgen_band_cat[ii];//cell_r_c_b  24 cells to fill
		register int fr0 = cd[crcb[2]] & bd[crcb[3]], fr = rd[crcb[1]] & fr0;
		if (crcb[4])if (_popcnt32(fr0) < 3) goto back; // 3 clues needed here
		if (!fr)goto back;
		free[ii] = fr;
	}
	goto next_first;
next:// erase previous fill and look for next
	crcb = tgen_band_cat[ii];
	d = zs0[crcb[0]];
	bit = 1 << d;
	rd[crcb[1]] ^= bit; cd[crcb[2]] ^= bit; bd[crcb[3]] ^= bit;
	if (!free[ii])goto back;
	{
	next_first:
		crcb = tgen_band_cat[ii];// be sure to have the good one
		bitscanforward(d, free[ii]);
		bit = 1 << d;
		free[ii] ^= bit;
		zs[crcb[0] + 27] = (char)(d + '1');
		zs0[crcb[0]] = d;
		rd[crcb[1]] ^= bit; cd[crcb[2]] ^= bit; bd[crcb[3]] ^= bit;
		if (ii < 23) goto nextii;
		// this is a valid band, check if lexically minimale 
		int ir = bandminlex.Getmin(zs0, &pband2, 0);
		if (ir < 0) {//would be bug  did not come in enumeration
			cerr << "Find B2 invalid return  Getmin" << endl;
			cout << "Find B2 invalid return  Getmin" << endl;
			return;
		}
		pcheck2 = pband2;
		it16_2 = pband2.i416;
		ib2check = i2t16 = t416_to_n6[it16_2];
		if (i2t16 < i1t16)goto next;// not canonical
		//if (i2t16 == i1t16)if (it16_2 < it16)goto next;// not canonical
		if (sgo.vx[4] == 1)memcpy(&gcheck[54], zs0, 27 * sizeof gcheck[0]);
		else {
			memcpy(&gcheck[27], zs0, 27 * sizeof gcheck[0]);
			if (Band2Check())goto next;// do nothing if p2b
		}
		nb12++;
		if (ValidBand2()) { cout << "stop b2" << endl;	go_back = 1; return; }
		goto next;
	}

back:
	if (--ii >= 0) goto next;
}
void Go_c17_91_go();
void GEN_BANDES_12::ValidInitGang() {
	for (int i = 0; i < 9; i++) {// init columns status
		cold[i] = 0x1ff;
		for (int j = 0; j < 6; j++)	cold[i] ^= 1 << grid0[i + 9 * j];
		gangb12[i] = 0x1ff ^ cold[i];
	}
	memcpy(gangcols, cold, sizeof gangcols);
}
int GEN_BANDES_12::ValidBand2() {
	myband2.InitBand2_3(it16_2, &zsol[27], pband2);
	//_______________________ std process
	if (modeb12 < 10) {
		nband3 = 0;
		//if (p_cpt2g[0]) return 1;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<	
		if ((nb12 >> 6) < skip) return 0;// here restart value, kept untouched if no band 3 found
		{// print a restart point every 64 bands 1+2 seen
			uint64_t w = genb12.nb12, w1 = w >> 6;
			w &= 63;
			if (w == 0) {
				long tfin = GetTimeMillis();
				cout << "next skip value to use=\t" << w1 << "\tmil=" << (tfin - sgo.tdeb) / 1000 << "\tnb2=" << p_cpt2g[0] << endl;
			}
		}
		ValidInitGang();// also called from existing 17 in test
		if (sgo.vx[7]) {
			if (sgo.vx[7] == i2t16) {
				Find_band3B(0);
			}
		}
		else Find_band3B();
		//cout << "band12=" << nb12 << "\tnband3=" << nband3 << " it16_2="<< it16_2 << endl;
		if ((nb12 >> 6) >= last)return 1;
		return 0;
	}
	//______________________ testing options 
	switch (modeb12) {
	case 10: {// test ua collection
		if (++p_cptg[0] > sgo.vx[1]) return 1;
		cout << "deux bandes à tester" << endl;
		Go_c17_91_go();
		return 0;
	}
	case 11: {// enumeration test
		for (int i = 0; i < 9; i++) {// init columns status
			cold[i] = 0x1ff;
			for (int j = 0; j < 6; j++)	cold[i] ^= 1 << grid0[i + 9 * j];
		}
		memcpy(gangcols, cold, sizeof gangcols);
		Find_band3B(11);
		if (nband3) {
			p_cpt[0]++;
			p_cpt[1] += nband3;
		}
		return 0;
	}
	}
	return 0;
}
void GEN_BANDES_12::Find_band3B(int m10) {
	//BANDMINLEX::PERM pout;
	register int  *crcb, bit;
	nband3 = 0;
	int *rd = rowdb3, *cd = cold, *bd = boxdb3; // to updates rows cols boxes
	char * zs = zsol;
	int * zs0 = &grid0[54];
	memcpy(boxdb3, &boxd[3], sizeof boxdb3);
	memcpy(rowdb3, &rowd[3], sizeof rowdb3);
	// now loop over the 24 cells not yet assigned in the band to fill the band use relative cell 
	int ii = -1, free[24];
	uint32_t d;
nextii:
	ii++;
	{
		crcb = tgen_band_cat[ii];//cell_row_col_box one of the 24 cells to fill
		register int fr = cd[crcb[2]] & bd[crcb[3]] & rd[crcb[1]];
		if (!fr)goto back;
		free[ii] = fr;
	}
	goto next_first;

next:// erase previous fill and look for next
	crcb = tgen_band_cat[ii];
	d = zs0[crcb[0]];
	bit = 1 << d;
	rd[crcb[1]] ^= bit; cd[crcb[2]] ^= bit; bd[crcb[3]] ^= bit;
	if (!free[ii])goto back;
	{
	next_first:
		crcb = tgen_band_cat[ii];// be sure to have the good one
		bitscanforward(d, free[ii]);
		bit = 1 << d;
		free[ii] ^= bit;
		zs[crcb[0] + 54] = (char)(d + '1');
		zs0[crcb[0]] = d;
		rd[crcb[1]] ^= bit; cd[crcb[2]] ^= bit; bd[crcb[3]] ^= bit;
		if (ii < 23) goto nextii;
		// this is a valid band, check if canonical 
		int ir = bandminlex.Getmin(zs0, &pband3, 0);
		if (ir < 0) {//would be bug  did not come in enumeration
			cerr << "gen band 3 invalid return Getmin" << endl;
			return;
		}
		int it16_3 = pband3.i416;
		ib3check=i3t16 = t416_to_n6[it16_3];
		if (sgo.vx[8] && i3t16 != sgo.vx[8]) goto next;
		if (i3t16 < i1t16)goto next;// not canonical
		if (sgo.vx[4] != 1) {
			if (i3t16 < i2t16)goto next;// not canonical (must be in this case
		}
		else if (i3t16 > i2t16)goto next;// not canonical (must be in this case
		//==============================  b1=b2=b3 use minlex check (simplest code, not common)
		if (sgo.vx[4] < 2) {// do nothing in p1
			if (sgo.vx[4] < 1) {// p2a
				pcheck3 = pband3;
				memcpy(&gcheck[54], zs0, 27 * sizeof gcheck[0]);
				if ( Band3Check())goto next;
			}
			else {// p2b exchanging band 2 band 3
				memcpy(&gcheck[27], zs0, 27 * sizeof gcheck[0]);				
				pcheck3 = pband2;
				pcheck2 = pband3;
				ib2check = i3t16;
				ib3check = i2t16;
				if (Band2Check())goto next;// band 3 must be a valid "band2"
				if (Band3Check())goto next;// then band 2 a valid band3
			}
		}
		bands3[nband3++].InitBand3(it16_3, &zs[54], pband3);
		if (sgo.vx[8]) {
			//fout1 << zs << endl;
			fout1 << zs << ";" << nb12<< " i2t16="<<i2t16<<" i3t16="<<i3t16 << endl;
		}
		goto next;
	}
back:
	if (--ii >= 0) goto next;
	if (m10 != 1)return;
	if (nband3)		g17b.GoM10();// call the process for that entry
}

int GEN_BANDES_12::DebugFreshUA(uint64_t ua) {
	// purely debugging code, a fresh UA <= limsize must have >5 digits
	uint32_t cell, digits = 0;
	while (bitscanforward64(cell, ua)) {
		ua ^= (uint64_t)1 << cell;
		digits |= 1 << grid0[From_128_To_81[cell]];
	}
	if (_popcnt32(digits) >= 6) return 0;
	cout << "bug more ua bands 1+2 not enough digits " << endl;
	return 1;
}
