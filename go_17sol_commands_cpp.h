
//========================================
const char * zh_g_cpt[10] = { "npuz", "guess", "close_d ", "upd1 ", "upd2 ",
"fupd ", "hpair ", "htripl ", " ", " " };

const char * libs_c17_00_cpt2g[60] = {
	"0 bands1+2 processed entries M10",//0
	"1 total bands 3",//1
	"2 steps external loop ",//2
	"3 loop band1 ",//3
	"4 loop band 2  ",//4
	"5 clean entries    ",//5
	"6 pass all filters",//6
	"7 calls brute force ",//7
	"8 valid brute force",//8
	"9 3B nb3 tot miss",//9
	"10 miss 0",//10
	"11 miss 1",//11
	"12 miss 2",//12
	"13 expand",//13
	"14 ",//14
	"15 3A nb ok",//15
	"16 b2 minplus 5",//16
	"17 3A blocs nb3",//17
	"18 3A blocs nb stack",//18
	"19 3A blocs nb sock2",//19
	"20 3A blocs nb ok",//20
	"21 entry 3a",//21
	"22 3a mincount <=6",//22
	"23 ",//23
	"24  ",//24
	"25 ",
	"26  ",
	"27  entry crit",	
	"28  need to assign", 
	"29  b3 other call check",
	"30  b3_expand call check",	
	"31 ajouts uas",
	"32  b3 new ua <=16 ",
	"33  b3 new ua b3=4",
	"34 TU_GUAN::DoStepb2()",
	"35 max ngua stepb2",	
	"36 ng2_2 + ng3_2",	
	"37 setv nguan",	
	"38 setv nguas", 
	"39 all steps",
	"40 b3 new size 2",
	"41 b3 new size 3",
	"42 ",
	"43 ",
	"44 ",
	"45 step dead  b2",
	"46 ",
	"47 ",
	"48 ",
	"49 ",
	"50 count clean",
	"51  <= switch bandes3 bandes",
	"52  > switch bandes3 bandes",
	"53  calls band3b",
	"54 call IsMultiple(int bf)",
	"55 zhou[1].CallMultipleB3",
	"56 nmiss1 ua of found",
	"57 nmiss1 go subcritical",
	"58",
	"59 test",
};
void Go_c17_00( ) {// p2 process
	cout << "Go_c17_00 search batch 17 clues 656 566 " << endl;
	cout << sgo.vx[0] << " -v0- band 0_415" << endl;
	cout << sgo.vx[2] << " -v2- skip first nnn restart after batch failure" << endl;
	cout << sgo.vx[3] << " -v3- last entry number for this batch must be > vx[2]" << endl;
	cout << sgo.vx[4] << " -v4- 0 if p2a 1 if p2b" << endl;
	if(sgo.vx[5])cout << sgo.vx[5] << " -v5- band filter" << endl;
	if (sgo.vx[6])cout << sgo.vx[6] << " -v6- diag option" << endl;

	int it16_start = sgo.vx[0];
	g17b.debug17 = g17b.debug17_check =   g17b.aigstop=0;
	//g17b.diag = sgo.vx[6];
	genb12.skip = sgo.vx[2];
	genb12.last = sgo.vx[3];
	if (sgo.vx[2] < 0) {
		cerr << "invalid value for skip" << endl;
		return;
	}
	if (sgo.vx[3] < sgo.vx[2]) {
		cerr << "invalid value for last to process" << endl;
		return;
	}
	zh_g.modevalid = 1;
	zh_g2.grid0 = genb12.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	memset(p_cptg, 0, sizeof p_cptg);// used in debugging sequences only
	memset(p_cpt1g, 0, sizeof p_cpt1g);// used in debugging sequences only
	memset(p_cpt2g, 0, sizeof p_cpt2g);// used in debugging sequences only
	genb12.Start(0);
	genb12.NewBand1(sgo.vx[0]);
	cout << "print final stats" << endl;
	for (int i = 0; i < 60; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
	cout << "exit final stats" << endl;
}

void Go_c17_09() {// p2 process locate 
	cout << "Go_c17_09 search locate batch 17 clues 656 566 " << endl;
	cout << sgo.vx[0] << " -v0- band 0_415" << endl;
	cout << sgo.vx[2] << " -v2- skip first nnn restart after batch failure" << endl;
	cout << sgo.vx[3] << " -v3- last entry number for this batch must be > vx[2]" << endl;
	cout << sgo.vx[4] << " -v4- 0 if p2a 1 if p2b" << endl;
	if (!sgo.vx[7] || !sgo.vx[8]) sgo.vx[7] = sgo.vx[8] = 0;
	cout << sgo.vx[7] << " -v7- band2 searched" << endl;
	cout << sgo.vx[8] << " -v8- band3 searched" << endl;
	int it16_start = sgo.vx[0];
	g17b.debug17 = g17b.debug17_check=0;
	g17b.diag =  sgo.vx[6];
	genb12.skip = sgo.vx[2];
	genb12.last = sgo.vx[3];
	if (sgo.vx[2] < 0) {
		cerr << "invalid value for skip" << endl;
		return;
	}
	if (sgo.vx[3] < sgo.vx[2]) {
		cerr << "invalid value for last to process" << endl;
		return;
	}
	zh_g.modevalid = 1;
	zh_g2.grid0 = genb12.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	memset(p_cpt2g, 0, sizeof p_cpt2g);// used in debugging sequences only
	genb12.Start(0);
	genb12.NewBand1(sgo.vx[0]);
}
//========================= known s17 file 10/19
void Go_c17_10( ) {
	zh_g.modevalid = 1;
	zh_g2.grid0 = genb12.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	char * ze = finput.ze;
	int * zs0 = genb12.grid0, npuz = 0;
	cout << "Go_c17_10() search 17 using a file having known 17 656 " << endl;
	while (finput.GetLigne()) {
		if(strlen(ze)<160) continue;// skip blank lines
		npuz++;
		g17b.npuz = npuz;
		g17b.a_17_found_here = 0;
		g17b.aigstop= 0;
		if (npuz <= (int)sgo.vx[2]) continue;
		if (npuz > (int)sgo.vx[3]) break;
		g17b.debug17 = sgo.vx[0];
		g17b.debug17_check = 0;
		//if (npuz >5) return;
		cout << "\n\nto process  n="<<dec << npuz <<" debug="<< g17b.debug17 << endl;
		long tdeb = GetTimeMillis();
		//================================ to avoid the 665 case
		int ncb3 = 0;
		for (int i = 0; i < 27; i++) 
			if (ze[i + 136] - '.')ncb3++;
		if (ncb3 == 5) {// change band3 <-> band2
			for (int i = 0; i < 27; i++) {
				char temp = ze[i + 27];	ze[i + 27] = ze[i + 54];	ze[i + 54] = temp;
				temp = ze[i + 109];	ze[i + 109] = ze[i + 136];	ze[i + 136] = temp;
			}
		}

		if (g17b.debug17)
			cout << ze <<  " to process  n="  << npuz << endl;

		// =======================morph entry 
		for (int i = 0; i < 81; i++)zs0[i] = ze[i] - '1';
		BANDMINLEX::PERM perm_ret;
		bandminlex.Getmin(zs0, &perm_ret);
		int ib1 = perm_ret.i416, ib1a = t416_to_n6[ib1];
		myband1.InitBand2_3(ib1, ze, perm_ret, 0);
		bandminlex.Getmin(&zs0[27], &perm_ret);
		int ib2 = perm_ret.i416, ib2a = t416_to_n6[ib2];
		myband2.InitBand2_3(ib2, &ze[27], perm_ret, 1);
		bandminlex.Getmin(&zs0[54], &perm_ret);
		int ib3 = perm_ret.i416, ib3a = t416_to_n6[ib3];
		genb12.bands3[0].InitBand3(ib3, &ze[54], perm_ret);
		genb12.nband3 = 1;
		ze[81] = 0;
		char * ze2 = &ze[82];
		g17b.p17diag.SetAll_0();
		for (int i = 0; i < 81; i++) if (ze2[i] != '.')
			g17b.p17diag.Set_c(i);
		if (g17b.debug17)
			cout << Char2Xout(g17b.p17diag.bf.u64[0]) << " b12 pattern for the 17" << endl;
		genb12.ValidInitGang();
		g17b.npuz = npuz;
		g17b.GoM10Known();
		//g17b.a_17_found_here = 1;
		if (!g17b.a_17_found_here) {
			cout << "puz="<<npuz << " failed to find the searched 17" << endl;
			cerr << "puz=" << npuz << " failed to find the searched 17" << endl;
			break;
		}
	}
	cout << "print final stats" << endl;
	for (int i = 0; i < 60; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
}


void Go_c17_11() {// extract cout xx5 file1 xx6
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	char * ze = finput.ze;
	int * zs0 = genb12.grid0, npuz = 0;
	cout << "Go_c17_11() split known 17 b3_5 b3_6 " << endl;
	while (finput.GetLigne()) {
		// =======================morph entry to have min n6 count in first
		for (int i = 0; i < 81; i++)zs0[i] = ze[i] - '1';
		BANDMINLEX::PERM perm_ret;
		bandminlex.Getmin(zs0, &perm_ret);
		int ib1 = perm_ret.i416, ib1a = t416_to_n6[ib1];
		bandminlex.Getmin(&zs0[27], &perm_ret);
		int ib2 = perm_ret.i416, ib2a = t416_to_n6[ib2];
		bandminlex.Getmin(&zs0[54], &perm_ret);
		int ib3 = perm_ret.i416, ib3a = t416_to_n6[ib3];
		if (ib1a > ib2a) {// change band2 <-> band1
			for (int i = 0; i < 27; i++) {
				char temp = ze[i];	ze[i] = ze[i + 27];	ze[i + 27] = temp;
				temp = ze[i + 82];	ze[i + 82] = ze[i + 109]; ze[i + 109] = temp;
				int w = ib1a;	ib1a = ib2a;	ib2a = w;
			}
		}
		if (ib1a > ib3a) {// change band3 <-> band1
			for (int i = 0; i < 27; i++) {
				char temp = ze[i];	ze[i] = ze[i + 54];		ze[i + 54] = temp;
				temp = ze[i + 82];	ze[i + 82] = ze[i + 136];	ze[i + 136] = temp;
				int w = ib1a;	ib1a = ib3a;	ib3a = w;
			}
		}
		int ncb3 = 0;
		for (int i = 0; i < 27; i++) {
			if (ze[i + 136] - '.')ncb3++;
		}
		if (ncb3 == 6) {
			fout1 <<&ze[82] << ";" << ib2a 
				<< ";" << ib3a <<endl;
				continue;
		}
		// now 5 clues in band 3 exchange bands
		for (int i = 0; i < 27; i++) {
			char temp = ze[i + 27];	ze[i + 27] = ze[i + 54];	ze[i + 54] = temp;
				temp = ze[i + 109];	ze[i + 109] = ze[i + 136];	ze[i + 136] = temp;
		}
		int w = ib2a;	ib2a = ib3a;	ib3a = w;
		cout << &ze[82]  << ";" << ib2a
			<< ";" << ib3a << endl;

	}

}


void Go_c17_12() {// check diagonal status in a 665
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	char * ze = finput.ze;
	int  zs0 [81],zs0d[81],
		npuz = 0;
	cout << "Go_c17_11() split known 17 b3_5 b3_6 " << endl;
	while (finput.GetLigne()) {
		// =======================morph entry to have min n6 count in first
		for (int i = 0; i < 81; i++) {
			zs0[i] = ze[i] - '1';
			zs0d[i]=ze[C_transpose_d[i]] - '1';
		}
		BANDMINLEX::PERM perm_ret;
		bandminlex.Getmin(zs0, &perm_ret);
		int ib1 = perm_ret.i416, ib1a = t416_to_n6[ib1];
		bandminlex.Getmin(&zs0[27], &perm_ret);
		int ib2 = perm_ret.i416, ib2a = t416_to_n6[ib2];
		bandminlex.Getmin(&zs0[54], &perm_ret);
		int ib3 = perm_ret.i416, ib3a = t416_to_n6[ib3];
		bandminlex.Getmin(zs0d, &perm_ret);
		int ib1d = perm_ret.i416, ib1ad = t416_to_n6[ib1d];
		bandminlex.Getmin(&zs0d[27], &perm_ret);
		int ib2d = perm_ret.i416, ib2ad = t416_to_n6[ib2d];
		bandminlex.Getmin(&zs0d[54], &perm_ret);
		int ib3d = perm_ret.i416, ib3ad = t416_to_n6[ib3d];
		cout << ze << ";\t" << ib1 << ";" << ib2 << ";" << ib3
			<<";\t" << ib1a << ";" << ib2a << ";" << ib3a
			<< ";\t" << ib1d << ";" << ib2d << ";" << ib3d 
			<< ";\t" << ib1ad << ";" << ib2ad << ";" << ib3ad
			<< endl;
	}

}

//=========================== regressive tests

void Go_c17_91_go() {
	if (p_cptg[0] > 1) return;
	myband1.PrintStatus();
	myband2.PrintStatus();
	genuasb12.Initgen();
}
void Go_c17_80() {// enumeration test
	cout << "Go_c17_20 phase 2a enumeration test " << endl;
	cout << sgo.vx[0] << " -v0- first id 0_415" << endl;
	cout << sgo.vx[1] << " -v1- second id 0_415" << endl;
	cout << sgo.vx[2] << " -v2- if 1 printout asked" << endl;
	cout << sgo.vx[4] << " -v4- 0 mode p2a 1 mode p2b 2 mode p1" << endl;
	int it16_start = sgo.vx[0], it16_end = sgo.vx[1];
	genb12.skip = 0;
	genb12.last = 100000;

	if (it16_start > 415 || it16_end > 415 || it16_start > it16_end) {
		cerr << "invalid it16_start it16_end" << endl;
		return;
	}
	memset(p_cptg, 0, sizeof p_cptg);// used in debugging sequences only
	genb12.Start(11);// mode p2a or p2b or p1
	for (int i1t16 = it16_start; i1t16 <= it16_end; i1t16++) {
		memset(p_cpt, 0, sizeof p_cpt);// band2 and band 3 count
		genb12.nb12 = 0;
		genb12.NewBand1(i1t16);
		cout << genb12.i1t16 << "\t" << genb12.it16 << "\t" << p_cpt[0]
			<< "\t" << p_cpt[1] << endl;
		p_cptg[0] += p_cpt[0];
		p_cptg[1] += p_cpt[1];
	}
	cout << "total\t\t" << p_cptg[0]
		<< "\t" << p_cptg[1] << endl;
}
void Go_c17_90() {// creating band UAs
	//the table exists, this is somehow a regressive test  
	//or to rebuild the table
	cout << "Go_10 entry  creating band uas test" << endl;
	STD_B416 wband;
	wband.Initstd();
	int maxua = 0;
	zh1b_g.diag = 0;
	for (int i416 = 0; i416 < 416; i416++) {
		//for (int i416 = 64; i416 <= 64; i416++) {	
		wband.InitC10(i416);
		int nua = t16_nua[i416]; // reference known
		zhone[0].InitOne_std_band();
		cout << wband.band << " i416=" << i416 << endl;

		for (int i = 0; i < 36; i++) {// find UAs 2 digits
			if (zhone[0].Start_nFloors(floors_2d[i])) continue;
			zhone[0].InitGuess();
			zhone[0].Guess2();
		}
		for (int i = 0; i < 84; i++) {// find UAs 3 digits
			if (zhone[0].Start_nFloors(floors_3d[i])) continue;
			zhone[0].InitGuess();
			zhone[0].Guess3();
		}
		if (t16_nua[i416] == zh1b_g.nua)goto ok;
		for (int i = 0; i < 126; i++) {// find UAs 4 digits
			if (zhone[0].Start_nFloors(floors_4d[i])) continue;
			zhone[0].InitGuess();
			zhone[0].Guess4();
		}
		for (int i = 0; i < 126; i++) {// find UAs 5 digits
			if (zhone[0].Start_nFloors(0x1ff ^ floors_4d[i])) continue;
			zhone[0].InitGuess();
			if (0 && floors_4d[i] == 0113) {
				zhone[0].ImageCandidats();
				zh1b_g.diag = 1;
			}
			else zh1b_g.diag = 0;
			zhone[0].Guess5();
		}
		for (int i = 0; i < 84; i++) {// find UAs 6 digits
			if (zhone[0].Start_nFloors(0x1ff ^ floors_3d[i])) continue;
			zhone[0].InitGuess();
			zhone[0].Guess6();
		}
		if (t16_nua[i416] == zh1b_g.nua)goto ok;

		//if(1)zh1b_g.PrintTua();

		zh1b_g.FindMissingUAs();

		if (t16_nua[i416] == zh1b_g.nua)goto ok;
		{
			cout << "uas final table" << endl;
			zh1b_g.PrintTua();
			cout << "uas known table" << endl;
			int id = t16_indua[i416], iff = id + t16_nua[i416];
			for (int i = id; i < iff; i++)
				cout << i - id << "\t"
				<< Char27out(t16_UAs[i]) << "\t" << _popcnt32(t16_UAs[i]) << endl;
			continue;
		}
	ok:
		cout << "ua count equal" << endl;
	}
}
void Go_c17_91() {// UAs collector 2 bands
	cout << "Go_11 test UAs collector 2 ba,ds" << endl;
	cout << sgo.vx[0] << " -v0- id bande low order 0_415" << endl;
	cout << sgo.vx[0] << " -v1- number of band2" << endl;
	memset(p_cpt, 0, sizeof p_cpt);// used in debugging sequences only
	memset(p_cpt1, 0, sizeof p_cpt1);// used in debugging sequences only
	genb12.Start(10);
	genb12.NewBand1(sgo.vx[0]);
}
/* sample used for Go_c17_92()
1 0     0       0       0       0       729     18225
28 27   0       38      2025    26285   157845  610977
29 28   9       567     7398    51516   237762  803574
412 411 0       0       0       324     10503   112491
413 412 0       0       0       0       729     18225
32 31   0       0       0       546     12291   101597
		*/
void Go_c17_92() {// test UAs 5 6 expand
	uint32_t maxstep = 0, maxvalid=0,
		max2 = 0, max3 = 0, max4 = 0, max5 = 0, max6 = 0;
	for (int ib = 0; ib < 416; ib++) {
		myband1.Initstd();
		myband1.InitG12(ib);
		g17b.ExpandB1();// expand and set index
		if (g17b.nmybi2t > maxstep) maxstep = g17b.nmybi2t;
		uint32_t last = g17b.mybi2t[g17b.nmybi2t - 1].iend;
		if (last > maxvalid) maxvalid = last;
		cout << "last=" << last << endl;

		// find max values in step
	}
	cout << endl << endl
		<< "maxstep=" << maxstep << " maxvalid=" << maxvalid;
}


//====== 15  extraction des 17 6 6 5
void Go_c17_15() {
	char * ze = finput.ze;

	while (finput.GetLigne()) {
		int count[6];
		memset(count, 0, sizeof count);
		for (int i = 0; i < 81; i++) {
			if (ze[i] - '.') {
				int band = i / 27, stack = C_stack[i] + 3;
				++count[band];
				++count[stack];
			}
		}
		if (count[0] > 6 || count[1] > 6 ||
			count[2] > 6 || count[3] > 6 ||
			count[4] > 6 || count[5] > 6)
			cout << ze << endl;
		else fout1 << ze << endl;
	}
}

void Msp_ReorderBand(char * ze,char * zep=0)
{
	char temp[9];
	if (ze[0] > ze[9]) {
		memmove(temp, ze, sizeof temp);
		memmove(ze, &ze[9], sizeof temp);
		memmove(&ze[9], temp,sizeof temp);
		if (zep) {
			memmove(temp, zep, sizeof temp);
			memmove(zep, &zep[9], sizeof temp);
			memmove(&zep[9], temp, sizeof temp);

		}
	}
	if (ze[0] > ze[18]) {
		memmove(temp, ze, sizeof temp);
		memmove(ze, &ze[18], sizeof temp);
		memmove(&ze[18], temp, sizeof temp);
		if (zep) {
			memmove(temp, zep, sizeof temp);
			memmove(zep, &zep[18], sizeof temp);
			memmove(&zep[18], temp, sizeof temp);
		}
	}
	if (ze[9] > ze[18]) {
		memmove(temp, &ze[9], sizeof temp);
		memmove(&ze[9], &ze[18], sizeof temp);
		memmove(&ze[18], temp, sizeof temp);
		if (zep) {
			memmove(temp, &zep[9], sizeof temp);
			memmove(&zep[9], &zep[18], sizeof temp);
			memmove(&zep[18], temp, sizeof temp);
		}
	}
}

void Go_c17_16() {
	char * ze = finput.ze,zout[82],zdiag[82],zew[82],zes[82], zesf[82];
	zout[81] = zdiag[81] = zes[81] = 0;
	zh_g2.zsol = zout;
	while (finput.GetLigne()) {
		if (zhou[0].CheckValidityQuick(ze) == 1) {
			int count[6];
			memset(count, 0, sizeof count);
			for (int i = 0; i < 81; i++) {
				if (ze[i] - '.') {
					int band = i / 27, stack = C_stack[i] + 3;
					++count[band];
					++count[stack];
				}
			}
			for (int i = 0; i < 81; i++)
				zdiag[C_transpose_d[i]] = zout[i];
			int zw[27], indexr[6], indexmin = 500;
			BANDMINLEX::PERM pr;
			char * zz = zout;
			for (int i = 0; i < 27; i++)
				zw[i] = zz[i] - '1';
			bandminlex.Getmin(zw, &pr);
			indexmin = indexr[0] = t416_to_n6[pr.i416];
			zz = &zout[27];
			for (int i = 0; i < 27; i++)
				zw[i] = zz[i] - '1';
			bandminlex.Getmin(zw, &pr);
			indexr[1] = t416_to_n6[pr.i416];
			if (indexmin > indexr[1])indexmin = indexr[1];
			zz = &zout[54];
			for (int i = 0; i < 27; i++)
				zw[i] = zz[i] - '1';
			bandminlex.Getmin(zw, &pr);
			indexr[2] = t416_to_n6[pr.i416];
			if (indexmin > indexr[2])indexmin = indexr[2];
			zz = zdiag;
			for (int i = 0; i < 27; i++)
				zw[i] = zz[i] - '1';
			bandminlex.Getmin(zw, &pr);
			indexr[3] = t416_to_n6[pr.i416];
			if (indexmin > indexr[3])indexmin = indexr[3];
			zz = &zdiag[27];
			for (int i = 0; i < 27; i++)
				zw[i] = zz[i] - '1';
			bandminlex.Getmin(zw, &pr);
			indexr[4] = t416_to_n6[pr.i416];
			if (indexmin > indexr[4])indexmin = indexr[4];
			zz = &zdiag[54];
			for (int i = 0; i < 27; i++)
				zw[i] = zz[i] - '1';
			bandminlex.Getmin(zw, &pr);
			indexr[5] = t416_to_n6[pr.i416];
			if (indexmin > indexr[5])indexmin = indexr[5];
			// sort and keep track of the count or split the file


			if (sgo.vx[6] && indexmin != sgo.vx[6]) continue;
			uint32_t tx[3], tdx[3], temp;
			tx[0] = (indexr[0] << 8) | count[0];
			tx[1] = (indexr[1] << 8) | count[1]|010;
			tx[2] = (indexr[2] << 8) | count[2]|020;
			tdx[0] = (indexr[3] << 8) | count[3];
			tdx[1] = (indexr[4] << 8) | count[4]|010;
			tdx[2] = (indexr[5] << 8) | count[5]|020;
			if (tx[0] > tx[1]) { temp = tx[0]; tx[0] = tx[1]; tx[1] = temp; }
			if (tx[0] > tx[2]) { temp = tx[0]; tx[0] = tx[2]; tx[2] = temp; }
			if (tx[1] > tx[2]) { temp = tx[1]; tx[1] = tx[2]; tx[2] = temp; }
			if (tdx[0] > tdx[1]) { temp = tdx[0]; tdx[0] = tdx[1]; tdx[1] = temp; }
			if (tdx[0] > tdx[2]) { temp = tdx[0]; tdx[0] = tdx[2]; tdx[2] = temp; }
			if (tdx[1] > tdx[2]) { temp = tdx[1]; tdx[1] = tdx[2]; tdx[2] = temp; }
			BF64 va, vb; va.bf.u64 = vb.bf.u64 = 0;
			va.bf.u16[2] = tx[0] >> 8; va.bf.u16[1] = tx[1] >> 8; va.bf.u16[0] = tx[2] >> 8;
			vb.bf.u16[2] = tdx[0] >> 8; vb.bf.u16[1] = tdx[1] >> 8; vb.bf.u16[0] = tdx[2] >> 8;
			if (sgo.vx[6]) {// check bands 2 3 for test
				cout << "va1=" << va.bf.u16[2] << " va2=" << va.bf.u16[1]
					<< " va3=" << va.bf.u16[0] << "\tvb1=" << vb.bf.u16[2] << " vb2=" << vb.bf.u16[1]
					<< " vb3=" << vb.bf.u16[0];
				if (va.bf.u64 <= vb.bf.u64)			cout << "\tdirect"<<endl;
				else {
					cout << "\tdiag"<<endl;
				}
			}
			strcpy(zew, ze);// to morph the expected 17 puzzle
			if (va.bf.u64 > vb.bf.u64)	{
				va = vb;
				memcpy(tx, tdx, sizeof tx);
				strcpy(zout, zdiag);
				for (int i = 0; i < 81; i++)zew[C_transpose_d[i]] = ze[i];
			}
			if (sgo.vx[6]) {// check bands 2 3 for test
				if (va.bf.u16[2] != sgo.vx[6]) continue;
				if (va.bf.u16[1] != sgo.vx[7]) continue;
			}
			char zs[82], zsa[82]; zsa[81] = zs[81] = 0;
			int i1 = (tx[0] & 070) >> 3, i2 = (tx[1] & 070) >> 3, i3 = (tx[2] & 070) >> 3;
			memcpy(zsa, &zout[27 * i1], 27);
			memcpy(&zsa[27], &zout[27 * i2], 27);
			memcpy(&zsa[54], &zout[27 * i3], 27);
			memcpy(zes, &zew[27 * i1], 27);
			memcpy(&zes[27], &zew[27 * i2], 27);
			memcpy(&zes[54], &zew[27 * i3], 27);
			int zs0[81];
			// redo id to build tables
			int ib1, ib1a, ib2, ib2a, ib3, ib3a;
			for (int i = 0; i < 81; i++)zs0[i] = zsa[i] - '1';
			BANDMINLEX::PERM pb1, prw;
			bandminlex.Getmin(zs0, &pb1);
			ib1 = pb1.i416;
			ib1a = t416_to_n6[ib1];
			//myband1.InitBand2_3(ib1, zsa, perm_ret, 0);
			bandminlex.Getmin(&zs0[27], &prw);
			ib2 = prw.i416;
			ib2a = t416_to_n6[ib2];
			//myband2.InitBand2_3(ib2, &zsa[27], perm_ret, 1);
			bandminlex.Getmin(&zs0[54], &prw);
			ib3 = prw.i416;
			ib3a = t416_to_n6[ib3];
			strcpy(zs, "12345678945");
			strncpy(&zs[11], t416[ib1], 16);
			strcpy(zesf, empty_puzzle);
			// relabel all
			for (int i = 3, ij = 27; i < 9; i++)
				for (int j = 0; j < 9; j++, ij++) {
					int is = 9 * i + pb1.cols[j];
					int c = zsa[is] - '1';
					zs[ij] = pb1.map[c] + '1';
					if (zes[is] != '.')zesf[ij] = zs[ij];
				}
			if (sgo.vx[6]) {
				cout << "\ti1,2,3 " << i1 << i2 << i3 << endl;
				cout << ze << " ze" << endl;
				cout << zes << " zes" << endl;

				cout << zout << endl;
				cout << zsa << endl;
				//genb12.bands3[0].InitBand3(ib3, &zsa[54], perm_ret);
				
			}
			Msp_ReorderBand(&zs[27], &zesf[27]);
			Msp_ReorderBand(&zs[54], &zesf[54]);

			// reorder band1 p17
			for (int i = 0, ij = 0; i < 3; i++)
				for (int j = 0; j < 9; j++, ij++) {
					int is = 9 * pb1.rows[i] + pb1.cols[j],
						c = zes[is];
					if (c != '.')zesf[ij] = zs[ij];
				}
			if (sgo.vx[6]) {
				cout << "ib1a=" << ib1a << " ib2a=" << ib2a << " ib3a=" << ib3a << endl;
				cout << zs << "zsa morphed relabeled" << endl;
				cout << zs << "zsa reordered" << endl;
				cout << zesf << "zes reordered relabelled" << endl;

			}
			fout1 << zs << ";" << zesf << ";" << indexmin << ";";
			fout1  << va.bf.u16[1] << ";"
				<< va.bf.u16[0] << ";" << (tx[2] & 07)  << endl;


		}
		else cout << ze << "non resolu" << endl;
	}
}



/*

void GEN_BANDES_12::Morph_ED_B12_known17(char * z, int ib1){
	i1t16 = ib1;  it16 = tn6_to_416[ib1];
	char zw[164], *zb1p = &z[82], *zb2p = &z[27 + 82], *zb3p = &z[54 + 82];
	strcpy(zw, z);
	unsigned long  g0[81], *z02 = &g0[27], *z03 = &g0[54],zcomp[27];
	for (int i = 0; i < 81; i++)g0[i] = z[i] - '1';
	__movsd(zcomp,z02,27);
	n_auto_b1 = bandminlex.GetAutoMorphs(it16, t_auto_b1);
	cout << "processing band rank=" << i1t16 << " N° 1_416=" << it16 + 1
		<< " n auto morphs=" << n_auto_b1 << endl;
	if (n_auto_b1){
		for (int imorph = 0; imorph < n_auto_b1; imorph++){
			BANDMINLEX::PERM &p = t_auto_b1[imorph];
			unsigned long band[27];// morph the band
			for (int i = 0; i < 9; i++){
				band[i] = p.map[z02[p.cols[i]]];
				band[i + 9] = p.map[z02[p.cols[i] + 9]];
				band[i + 18] = p.map[z02[p.cols[i] + 18]];
			}
			if( G17ComparedOrderedBand((int *)zcomp,(int *) band)==1){//morph to p orderd
				int tsort[3];
				G17BuildSort((int *)band, tsort);
				cout << p.rows[0] << p.rows[1] << p.rows[2] << " rows"<<endl;
				for (int i = 0; i < 3; i++){
					int ir = tsort[i] &= 3,irb1=p.rows[i];
					__movsd(&zcomp[9*i], &band[9 * ir], 9);// new zs0
					for (int j = 0; j < 9; j++){// new entry
						int wc = band[9 * ir + j]+ '1',wc2=wc;
						if (zb2p[ 9 * ir + p.cols[j]] == '.') wc2 = '.';
						zw[27 + 9 * i + j] = wc;
						zw[27+82 + 9 * i + j] = wc2;
						wc = p.map[z03[9 * i + p.cols[j]] ] + '1';
						wc2 = wc;
						if (zb3p[ 9 * i + p.cols[j]] == '.') wc2 = '.';
						zw[54 + 9 * i + j] = wc;
						zw[54 + 82 + 9 * i + j] = wc2;
						wc2 = p.map[g0[9 * irb1 + p.cols[j]]] + '1';
						if (z[82+9 * irb1 + p.cols[j]] == '.') wc2 = '.';
						zw[82 + 9 * i + j] = wc2;

					}
				}
				cout << zw << "après morph une etape" << endl;
			}
		}
	}
	cout << zw << "après morph ED band 12" << endl;
}
*/

void msp_change12(char * ze, int &ib1a, int & ib2a) {
	for (int i = 0; i < 27; i++) {
		char temp = ze[i];	ze[i] = ze[i + 27];	ze[i + 27] = temp;
		temp = ze[i + 82];	ze[i + 82] = ze[i + 109]; ze[i + 109] = temp;
		int w = ib1a;	ib1a = ib2a;	ib2a = w;
	}
}
void msp_change13(char * ze, int &ib1a, int & ib3a) {
	for (int i = 0; i < 27; i++) {
		char temp = ze[i];	ze[i] = ze[i + 54];		ze[i + 54] = temp;
		temp = ze[i + 82];	ze[i + 82] = ze[i + 136];	ze[i + 136] = temp;
		int w = ib1a;	ib1a = ib3a;	ib3a = w;
	}
}
void msp_change23(char * ze, int &ib2a, int & ib3a) {
	msp_change12(&ze[27], ib2a, ib3a);
}
void M1_S17_Morph(char * z, BANDMINLEX::PERM & p1) {
	char zw[164]; strcpy(zw, z); // pose 0 et ;
	// map band 1
	for (int ir = 0, i = 0; ir < 3; ir++)for (int ic = 0; ic < 9; ic++, i++) {
		int iw = 9 * p1.rows[ir] + p1.cols[ic];
		char cw1 = p1.map[z[iw] - '1'] + '1', cw2 = cw1;
		if (z[iw + 82] == '.') cw2 = '.';
		zw[i] = cw1; zw[i + 82] = cw2;
	}
	{// map band2 to band 1 reorder columns
		int iw4 = 27 + p1.cols[0], iw5 = 36 + p1.cols[0], iw6 = 45 + p1.cols[0];
		int cw4 = p1.map[z[iw4] - '1'], cw5 = p1.map[z[iw5] - '1'], cw6 = p1.map[z[iw6] - '1'];
		int tsort[3], w;// sort bands / stacks increasing order of the  id
		tsort[0] = cw4 << 8;
		tsort[1] = 1 | (cw5 << 8);
		tsort[2] = 2 | (cw6 << 8);
		for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
			if (tsort[i] > tsort[j]) { w = tsort[i]; tsort[i] = tsort[j]; tsort[j] = w; }
		tsort[0] &= 3;  tsort[1] &= 3; tsort[2] &= 3;
		for (int ir = 0, i = 27; ir < 3; ir++)for (int ic = 0; ic < 9; ic++, i++) {
			int iw = 9 * tsort[ir] + p1.cols[ic] + 27;
			char cw1 = p1.map[z[iw] - '1'] + '1', cw2 = cw1;
			if (z[iw + 82] == '.') cw2 = '.';
			zw[i] = cw1; zw[i + 82] = cw2;
		}

	}

	{// map band3 to band 1 reorder columns
		int iw5 = 54 + p1.cols[0], iw6 = 63 + p1.cols[0], iw7 = 72 + p1.cols[0];
		int cw5 = p1.map[z[iw5] - '1'], cw6 = p1.map[z[iw6] - '1'], cw7 = p1.map[z[iw7] - '1'];
		int tsort[3], w;// sort bands / stacks increasing order of the  id
		tsort[0] = cw5 << 8;
		tsort[1] = 1 | (cw6 << 8);
		tsort[2] = 2 | (cw7 << 8);
		for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
			if (tsort[i] > tsort[j]) { w = tsort[i]; tsort[i] = tsort[j]; tsort[j] = w; }
		tsort[0] &= 3;  tsort[1] &= 3; tsort[2] &= 3;
		for (int ir = 0, i = 54; ir < 3; ir++)for (int ic = 0; ic < 9; ic++, i++) {
			int iw = 9 * tsort[ir] + p1.cols[ic] + 54;
			char cw1 = p1.map[z[iw] - '1'] + '1', cw2 = cw1;
			if (z[iw + 82] == '.') cw2 = '.';
			zw[i] = cw1; zw[i + 82] = cw2;
		}

	}
	strcpy(z, zw);
	//zw[54] = 0; zw[54 + 82] = 0;
	//fout1 << zw  << endl;
}
void M1_S17(char * finput_name, char * foutput_name, uint32_t * vx, uint32_t * bx) {
	// temporary code extract solution with a known 17 6 6 5
	finput.open(finput_name);
	if (!finput.is_open()) { cerr << "error open file " << finput_name << endl; return; }
	fout1.open(foutput_name);
	char * ze = finput.ze, *ze2 = &ze[82];
	char zdiag[164], zw[164];
	zdiag[81] = ';'; zdiag[163] = 0;
	int count[6], zwi[81], zei[81], zdiagi[81];
	while (finput.GetLigne()) {
		if (strlen(ze) < 163)continue;
		ze[163] = 0; // forget more info
		memset(count, 0, sizeof count);
		for (int i = 0; i < 81; i++) {
			if (ze2[i] - '.') {
				int band = i / 27, stack = C_stack[i] + 3;
				++count[band];
				++count[stack];
			}
			zdiag[C_transpose_d[i]] = ze[i];
			zdiag[C_transpose_d[i] + 82] = ze[i + 82];
		}
		BANDMINLEX::PERM pr[6];
		for (int ib = 0; ib < 6; ib++)if (count[ib] > 6) {
			if (0) goto next;
			int x;
			// morph the puzzle to the high band in band
			if (ib > 2)memmove(zw, zdiag, sizeof zw);
			else memmove(zw, ze, sizeof zw);
			// put the >6 clues band as band3
			if (ib == 0 || ib == 3) msp_change13(zw, x, x);
			else if (ib == 1 || ib == 4) msp_change23(zw, x, x);
			for (int i = 0; i < 54; i++) zwi[i] = zw[i] - '1';
			bandminlex.Getmin(zwi, &pr[0]);
			int ib1 = pr[0].i416, ib1a = t416_to_n6[ib1];
			bandminlex.Getmin(&zwi[27], &pr[1]);
			int ib2 = pr[1].i416, ib2a = t416_to_n6[ib2],
				i2 = ib2a, i1 = ib1a;
			if (ib1a > ib2a) {
				msp_change12(zw, i1, i2);
				M1_S17_Morph(zw, pr[1]);
			}
			else 	M1_S17_Morph(zw, pr[0]);

			fout1 << zw << ";" << i1 << ";" << i2 << endl;
			cout << zw << endl;
			//genb12.Morph_ED_B12_known17(zw, i1);
			if (0) {// mode look for 566 656 no 7 in stack
				char zout[164];
				for (int i = 0; i < 81; i++) {
					zei[i] = ze[i] - '1';
					zdiagi[i] = zdiag[i] - '1';
				}
				BANDMINLEX::PERM perm_reti[6];
				int ibi[6], ibri[6];
				for (int i = 0; i < 3; i++) {
					bandminlex.Getmin(&zei[27 * i], &perm_reti[i]);
					ibi[i] = perm_reti[i].i416;
					ibri[i] = t416_to_n6[ibi[i]];
					bandminlex.Getmin(&zdiagi[27 * i], &perm_reti[i + 3]);
					ibi[i + 3] = perm_reti[i + 3].i416;
					ibri[i + 3] = t416_to_n6[ibi[i + 3]];
				}
				int tsort[3], w;// sort bands / stacks increasing order of the  id
				tsort[0] = ibri[0] << 8;
				tsort[1] = 1 | (ibri[1] << 8);
				tsort[2] = 2 | (ibri[2] << 8);
				for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
					if (tsort[i] > tsort[j]) { w = tsort[i]; tsort[i] = tsort[j]; tsort[j] = w; }
				int ib1 = tsort[0] & 3, ib2 = tsort[1] & 3, ib3 = tsort[2] & 3;
				tsort[0] = 3 | ibri[3] << 8;
				tsort[1] = 4 | (ibri[4] << 8);
				tsort[2] = 5 | (ibri[5] << 8);
				for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
					if (tsort[i] > tsort[j]) { w = tsort[i]; tsort[i] = tsort[j]; tsort[j] = w; }
				int ibd1 = tsort[0] & 3, ibd2 = tsort[1] & 3, ibd3 = tsort[2] & 3;
				char *z = ze;
				int * cpt = count, iout = ibri[ib1];
				if (ibri[ib1] < ibri[ibd1]) goto bands;
				if (ibri[ib1] > ibri[ibd1]) goto stacks;
				if (ibri[ib2] < ibri[ibd2]) goto bands;
				if (ibri[ib2] > ibri[ibd2]) goto stacks;
				goto bands;// not a key point, use band
			stacks:
				z = zdiag;
				cpt = &count[3];
				ib1 = ibd1 - 3; ib2 = ibd2 - 3; ib3 = ibd3 - 3;
				iout = ibri[ibd1];
			bands: {
				strcpy(zout, z);
				if (ib1) {
					memcpy(zout, &ze[27 * ib1], 27);
					memcpy(&zout[82], &ze[27 * ib1 + 82], 27);
				}
				if (count[ib3] == 5) {//  must be pass2b
					if (ib3 != 2) {
						memcpy(&zout[27], &z[27 * ib3], 27);
						memcpy(&zout[82 + 27], &z[27 * ib3 + 82], 27);
					}
					if (ib2 != 3) {
						memcpy(&zout[54], &z[27 * ib2], 27);
						memcpy(&zout[82 + 54], &z[27 * ib2 + 82], 27);
					}
					fout1 << z << ";" << iout << endl;
				}
				else {//  pass2a
					if (ib2 != 2) {
						memcpy(&zout[27], &z[27 * ib2], 27);
						memcpy(&zout[82 + 27], &z[27 * ib2 + 82], 27);
					}
					if (ib3 != 3) {
						memcpy(&zout[54], &z[27 * ib3], 27);
						memcpy(&zout[82 + 54], &z[27 * ib3 + 82], 27);
					}
					cout << z << ";" << iout << endl;

				}
				}
			}
		}
	next:;
	}
}


