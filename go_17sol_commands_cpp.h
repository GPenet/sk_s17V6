
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
	"14 miss 2 expand",//14
	"15 3A nb ok",//15
	"16 b2 minplus 5",//16
	"17 3A blocs nb3",//17
	"18 3A blocs nb stack",//18
	"19 3A blocs nb sock2",//19
	"20 3A blocs nb ok",//20
	"21 pass mincount",//21
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
	"46 n cleamore",
	"47 n cleanplus",
	"48 n clean b",
	"49 ",
	"50 count clean",
	"51 clean0 ",
	"52 cleanmore  ",
	"53 cleanplus  ",
	"54 clean b ",
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
	g17b.aigstop=0;
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
