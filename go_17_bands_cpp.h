

// standard first band (or unique band)

void STD_B416::Initstd() {
	strcpy(band, "12345678945");
	for (int i = 0; i < 11; i++) band0[i] = band[i] - '1';
	for (int i = 0; i < 27; i++)map[i] = i;// no mapping
	dband = 0;
}
void STD_B416::GetBandTable(int i) {
	i416 = i;
	strncpy(&band[11], t416[i], 16);
	band[27] = 0;
	for (int i = 11; i < 27; i++) band0[i] = band[i] - '1';
}
void STD_B416::SetGangster() {
	memset(gangster, 0, sizeof gangster);
	for (int ir = 0, i = 0; ir < 3; ir++)
		for (int ic = 0; ic < 9; ic++, i++)
			gangster[ic] |= 1 << band0[i];
	// build sol per digit and pm per digit at start
	memset(fd_sols, 0, sizeof fd_sols);
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) {
			int cell = 9 * j + i, dig = band0[cell];
			fd_sols[1][dig] |= Zhoucol << i; // add candidates in the column
			fd_sols[0][dig] |= 1 << cell;
		}
	}


}
void STD_B416::InitC10(int i) {
	GetBandTable(i); SetGangster();
	zh1b_g.GetBand(band0, tua);// set zhone_i
}
void STD_B416::InitG12(int i) {
	GetBandTable(i); SetGangster(); GetUAs();
	MorphUas(); // to add the count
}
void STD_B416::MorphUas() {
	// morph all uas
	for (uint32_t i = 0; i < nua; i++) {
		register int uao = tua[i]&BIT_SET_27, ua = 0;
		register uint32_t cc;
		while (bitscanforward(cc, uao)) {
			uao ^= 1 << cc;
			ua |= 1 << map[cc];
		}
		tua[i] = ua | _popcnt32(ua) << 27;// store with count
	}
}

void STD_B416::InitBand2_3(int i16, char * ze, BANDMINLEX::PERM & p
,int iband) {
	i416 = i16;
	dband = 27*iband;
	GetUAs();
	strncpy(band, ze, 27);
	for (int i = 0; i < 27; i++) band0[i] = band[i] - '1';
	// create the cell map in output
	for (int i = 0; i < 3; i++) {
		int vr = 9 * p.rows[i], vr0 = 9 * i;
		for (int j = 0; j < 9; j++)
		map[vr0 + j] = vr + p.cols[j];
	}
	MorphUas();// morph all uas

	SetGangster();
}
void STD_B416::PrintStatus() {
	cout << "band status i=" << i416 << "\tstart=" << dband<<endl<<"map ";
	for (int i = 0; i < 27; i++)cout << map[i] << " ";
	cout <<endl;
	cout << band << endl<<"gangster status"<<endl;;
	zh1b_g.GetBand(band0, tua);// set zhone_i
	zhone[0].InitOne_std_band();
	zh1b_g.ndigits = 9;
	zhone[0].ImageCandidats(); // gangster status
	cout << "UAs table" << endl;
	for (uint32_t i = 0; i < nua; i++)
		cout << Char27out(tua[i]) << endl;
}
void STD_B1_2::FillMiniDigsMiniPairs(STD_B1_2 & bb) {
	if (0) {
		cout << "gangster other band" << oct << endl;
		for (int i = 0; i < 9; i++)cout << bb.gangster[i] << " ";
		cout << endl << "my gangster  " << endl;
		for (int i = 0; i < 9; i++)cout << gangster[i] << " ";
		cout << dec << endl;
		cout << band << " my band" << endl;
	}

	nvpairs = 0;
	for (int i = 0, j = 0; i < 9; i++, j += 3) {
		int a = (1 << band0[j]), b = (1 << band0[j + 1]), c = (1 << band0[j + 2]);
		mini_digs[i] = a | b | c;
		mini_pairs[j] = b | c;// missing a  relative columns 2,3
		mini_pairs[j+1] = a | c;// missing b
		mini_pairs[j+2] = a | b;// missing c 
		int jcol = j % 9;// start col for the mini row
		int * gg = bb.gangster;
		if ((gg[jcol + 1] & c) && (gg[jcol + 2] & b))
			tv_pairs[nvpairs++] = j;
		if ((gg[jcol ] & c) && (gg[jcol + 2] & a))
			tv_pairs[nvpairs++] = j+1;
		if ((gg[jcol] & b) && (gg[jcol + 1] & a))
			tv_pairs[nvpairs++] = j + 2;
	}
	if (0) {
		cout << " valid pairs table ";
		for (int i = 0; i < nvpairs; i++)
			cout << tv_pairs[i] << " ";
		cout << endl;
		cout << " valid pairs pairs " << oct;
		for (int i = 0; i < nvpairs; i++)
			cout << mini_pairs[tv_pairs[i]] << " ";
		cout << dec << endl;
	}
}
int STD_B1_2::ReviseG_triplet(int imini, int ip, STD_B1_2 * bb) {
	int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
	int dcell=3*imini,dcol= dcell % 9,
		*cols = &bb->gangster[dcol],
		*myp= tp3f[ip];
	int digit[3], digit_bit[3], pdigit[3], pdigit_bit[3], digp = 0;
	for (int i = 0; i < 3; i++) {// collect digits
		digit[i] = band0[dcell + i];
		int bit = 1 << digit[i];
		digit_bit[i] = bit;
		digp |= bit;
	}
	for (int i = 0; i < 3; i++) {// collect digits perm
		pdigit[i] = digit[myp[i]];
		pdigit_bit[i] = 1 << pdigit[i];
	}
	int *g12 = &zh2b_g.gangster[dcol];
	if (!(pdigit_bit[0] & g12[0]) ||
		!(pdigit_bit[1] & g12[1]) ||
		!(pdigit_bit[2] & g12[2])) return 0;
	for (int ic = 0; ic < 3; ic++)
		bb->revised_g[dcol + ic] ^= (pdigit_bit[ic] | digit_bit[ic]);
	return digp;
}

uint32_t  STD_B1_2::GetMiniData(int index, uint32_t & bcells, STD_B1_2 *bb) {
	//index is cell 0-26 in the band assumed free in the mini-row
	int tpcol[3][2] = { {1,2},{0,2},{0,1} };
	uint32_t tcells[3] = { 6,5,3 };// corresponding pairs
	int imini= index / 3,	dmini = 3*imini,dcol=dmini%9,
		 perm = index % 3,	*pcol = tpcol[perm];
	bcells = tcells[perm] << (3 * imini);
	uint32_t digs= mini_pairs[index];
	//bb->InitRevisedg();// must be done by the caller
	bb->revised_g[dcol + pcol[0]] ^= digs;
	bb->revised_g[dcol + pcol[1]] ^= digs;
	return digs;
}


void STD_B1_2::PrintShortStatus() {
	cout  << band << "\t\tband main data i0-415="<<i416 << endl;
	cout << "nua    \t" << nua << endl;

}

void STD_B3::InitBand3(int i16, char * ze, BANDMINLEX::PERM & p) {
	InitBand2_3(i16, ze, p);
	memset(&guas, 0, sizeof guas);
	// setup minirows bit fields
	for (int i = 0; i < 9; i++) {
		minirows_bf[i] = 0;
		int * p = &band0[3 * i];
		for (int j = 0; j < 3; j++)
			minirows_bf[i] |= 1 << p[j];
	}
	//cout << band << " " << oct << minirows_bf[0]
	//	<< " " << minirows_bf[1] << dec << endl;

}
/* GUA4 GUA6
.x. .y.
..y .x.

.x. ... .y.   c2a=c2b c3a=c3b
..y .x. ...   3 rows
... .y. .x.

.x. .y. .z.   2 rowssame third
..y .z. .x.
*/
int STD_B3::IsGua(int i81) {
	GEN_BANDES_12::SGUA2 w81 = genb12.tsgua2[i81];
	int * d10 = &band0[w81.col1], *d20 = &band0[w81.col2];
	int d1 = w81.dig1, d2 = w81.dig2,r1,r2,r3,ua,cell1,cell2;
	for (int irow = 0; irow < 3; irow++) {
		int c1 = d10[9 * irow], c2 = d20[9 * irow];
		if (c1 == d1 && c2 == d2) {// gua2
			r1 = r2 = irow;		cell1 = 9 * irow+ w81.col1;
			cell2 = 9 * irow + w81.col2;			break;
		}
		if (c1 == d1) {	r1 = irow; cell1 = 9 * irow + w81.col1;		}
		else if (c2 == d2) {r2 = irow; cell2 = 9 * irow + w81.col2;	}
		else r3 = irow;
	}
	ua = (1 << cell1) | (1 << cell2);
	if (r1 == r2) {// gua2
		guas.ua_pair[i81] = ua;
		int i27 = 9 * r1 + w81.i9;// index 0-26 of the pair
		i_27_to_81[i27] = i81;
		guas.isguasocket2.Set_c(i81);
		int imini= 3 * r1 + w81.i9 / 3;
		guas.ua2_imini[i81] = imini;
		guas.ua2pair27[i81] = (7 << (3 * imini))^ua;
		guas.ua2bit[i81] = 1 << imini;
		//cout << "set gua2 i81=" << i81 << " imini=" << 3 * r1 + w81.i9 / 3 
		//	<<	" guas.ua2_i27="<< guas.ua2_i27[i81] << endl;;
		return 1;
	}
	// is it a gua4 gua6 catch data to use
	int tc[2], ntc = 0,digs=w81.digs;//colums with the 2 digits
	int col1, col2, *p1 = &band0[9 * r1], *p2 = &band0[9 * r2];
	for (int i = 0; i < 9; i++) if ((gangster[i] & digs) == digs)
		tc[ntc++] = i;
	for (int icol = 0; icol < 9; icol++, p1++, p2++) {
		if (*p1 == d2)col1 = icol;
		if (*p2 == d1)col2 = icol;
	}


	if (ntc == 2) {// gua6 first type find and store ua
		int cella = 9 * r1 + col1, cellb = 9 * r2 + col2,
			cellc = 9 * r3 + col1, celld = 9 * r3 + col2;
		ua |= (1 << cella) | (1 << cellb) | (1 << cellc) | (1 << celld);
		guas.ua_pair[i81] = ua;
		guas.isguasocket2_46.Set_c(i81);
		return 4;
	}
	if (ntc) {
		int c = band0[9 * r3 + tc[0]];
		if (c != d1 && c != d2) {// gua4 
			int cella = 9 * r1 + col1, cellb = 9 * r2 + col2;
			ua |= (1 << cella) | (1 << cellb);
			guas.ua_pair[i81] = ua;
			guas.isguasocket2_46.Set_c(i81);
			return 2;
		}
	}
	// last  is gua6 with d1,r2 ; d2,r1
	if (band0[9 * r1 + col2] == band0[9 * r2 + col1]) {// gua6 second type
		int cella = 9 * r1 + col1, cellb = 9 * r2 + col2,
			cellc = 9 * r2 + col1, celld = 9 * r1 + col2;
		ua |= (1 << cella) | (1 << cellb) | (1 << cellc) | (1 << celld);
		guas.ua_pair[i81] = ua;
		guas.isguasocket2_46.Set_c(i81);
		return 8;
	}
	return 0;
}
int STD_B3::IsGua3(int i81) {
	GEN_BANDES_12::SGUA3 w81 = genb12.tsgua3[i81];
	int *g = genb12.gang27,//the gangster 
		d1= g[w81.id1],d2= g[w81.id2],d3= g[w81.id3];
	//catch the minirow pattern needed
	int bita = 1 << d1, bitb = 1 << d2, bitc = 1 << d3;
	int mrpat = bita | bitb | bitc;
	int stack = w81.col1 / 3;
	// minirow of the stack must fit
	for (int irow = 0; irow < 3; irow++) {
		int imini = stack + 3 * irow, 
			*pmini = &band0[9 * irow + 3 * stack];
		i_9_to_81[imini] = i81;
		if (mrpat != minirows_bf[imini])continue;
		// possible triplet, must be right digit in right place
		if (d1 != pmini[0] || d2 != pmini[1])continue;
		guas.triplet[imini] = i81;// valid triplet
		guas.triplet_imini[i81] = imini;// valid triplet
		guas.isguasocket3.Set_c(i81);
		guas.ua_triplet[i81] = 7 << (3 * imini);
		guas.ua3_imini[i81] = imini;
		guas.ua3bit[i81] =1<< imini;
		return 1;
	}
	return 0;
}

void STD_B3::PrintB3Status() {
	cout << "band3 status" << endl;
	for (int i = 0, ij = 0; i < 3; i++) {
		for (int j = 0; j < 9; j++, ij++) cout << band0[ij] + 1;
		cout << endl;
	}
	cout << endl<<"gua2 gua4 gua6s" << endl;
	for (int i = 0; i < 81; i++){
		int w = guas.ua_pair[i];
		if (w) cout << Char27out(w) << " i81=" << i << endl;
	}
	cout  << "gua3s" << endl;
	for (int i = 0; i < 81; i++) {
		int w = guas.ua_triplet[i];
		if (w) cout << Char27out(w) << " i81=" << i << endl;
	}
	char ws[129];
	const char* w3 = "123456789...---...123456789";
	cout << w3 << w3 << w3 << endl;;
	cout << guas.isguasocket2.String128(ws) << " sock2" << endl;
	cout << guas.isguasocket2_46.String128(ws) << " sock2_46" << endl;
	cout << guas.isguasocket3.String128(ws) << " sock3" << endl;
}

//==================== sockets UA2s UA3s control

int GENUAS_B12::Initgen() {// buil start myband1 myband2
	limsize = UALIMSIZE;
	zh2b5_g.sizef5 = UALIMSIZE;
	zh2b5_g.modevalid = 0;
	// prepare zh2b_g___________________________________________
	memcpy(zh2b_g.puz0, myband1.band0, sizeof myband1.band0);
	memcpy(&zh2b_g.puz0[27], myband2.band0, sizeof myband2.band0);
	for (int i = 0; i < 9; i++)
		zh2b_g.gangster[i] = myband1.gangster[i] | myband2.gangster[i];
	zh2b_g.GetBands(myband1.gangster, myband2.gangster);// set sol/pm
	zh2b_i.Init_std_bands();
	//_______________________________________________________
	nua = 0;// final table of uas bands 12 empty at start

/*
	for (uint32_t i = 0; i < myband1.nua; i++) {// collect band 1
		register uint64_t  ua = myband1.tua[i] & BIT_SET_27;
		uint64_t cc = _popcnt64(ua);
		ua |= cc << 59;
		AddUA64(tua, nua, ua);
	}
	for (uint32_t i = 0; i < myband2.nua; i++) {// collect band 2
		register uint64_t  ua = myband2.tua[i] & BIT_SET_27;
		uint64_t cc = _popcnt64(ua);
		ua <<= 32;
		ua |= cc << 59;
		AddUA64(tua, nua, ua);
	}*/
	nuab1b2 = 0;// switch uas band1 and uas band2 to 2X mode
	for (uint32_t i = 0; i < myband1.nua; i++) // collect band 1
		tuab1b2[nuab1b2++] = myband1.tua[i] & BIT_SET_27;
	for (uint32_t i = 0; i < myband2.nua; i++) {// collect band 2
		register uint64_t  R = myband2.tua[i] & BIT_SET_27;
		R <<= 32;
		tuab1b2[nuab1b2++] = R;
	}

	//___________________________ Start collection of uas
	zh2b_g.nua = 0;// new uas 
	for (int i = 0; i < 36; i++) BuildFloorsAndCollectOlds(floors_2d[i]);
	for (int i = 0; i < 84; i++) BuildFloorsAndCollectOlds(floors_3d[i]);
	for (int i = 0; i < 126; i++)BuildFloorsAndCollectOlds(floors_4d[i]);
	for (int i = 0; i < 126; i++) BuildFloorsAndCollectOlds(0x1ff ^ floors_4d[i]);
	//if(g17b.debug17&&g17b.GodebugCheckUas("standars uas")) return 1;

	//==================== collect more uas 6/7 digit 
	//cout << "initial status for UAS bands 1+2 nua="<<nua << endl;
	CollectMore2digits();
	myband1.FillMiniDigsMiniPairs(myband2);
	myband2.FillMiniDigsMiniPairs(myband1);
	//cout << "after collect more nua=" << nua << endl;
	CollectTriplets();
	//cout << "after collect triplets nua=" << nua << endl;
	CollectMore2minirows();
	//cout << "after collect 2 mini rows nua=" << nua << endl;
	//if (g17b.debug17&&g17b.GodebugCheckUas(" all us")) return 1;
	//DebugUas();
	return 0; // ok
}
int GENUAS_B12::DebugUas() {
	cout << "  debug uas" << endl;

	for (uint32_t i = 0; i < nua; i++) {
		uint64_t w = tua[i];
		int cc = (w >> 59);
		//if(cc==13)
		cout << Char2Xout(w) << " " << cc << " i=" << i << endl;
	}
	cout << " end debug uas" << endl;
	return 0;
}

void GENUAS_B12::BuildFloorsAndCollectOlds(int fl) {
	int diag = 0;
	floors = fl;// for debugging only
	uint64_t solved_cells = zh2b5_g.FindUAsInit(fl, 0);
	if (!solved_cells) return;// one digit solved true
	// now collect UAs not hit by solved cells  
	nuaold = 0;
	{	register uint64_t R = solved_cells;
	for (uint32_t i = 0; i < nuab1b2; i++)
		if (!(R & tuab1b2[i])) tuaold[nuaold++] = tuab1b2[i];
	/*
	*/
	for (uint32_t i = 0; i < nua; i++)
		if (!(R & tua[i])) tuaold[nuaold++] = tua[i];
	}
	zh2b5_g.CollectUas5();// collect uas for this set of floors
	// check subsets and add to main table
	for (uint32_t i = 0; i < zh2b5_g.nuaf5; i++) {
		ua = zh2b5_g.tuaf5[i].bf.u64&BIT_SET_2X;// be sure to use only relevant bits
		uint64_t cc = _popcnt64(ua);
		if (cc > limsize) continue;
		if (CheckOld()) continue;// superset of a previous ua
		ua |= cc << 59;
		AddUA64(tua, nua, ua);
	}

}
int GENUAS_B12::BuilOldUAs(uint32_t r0) {
	//====extract uas not hit in the band where is the mini row
	nuaold = 0;
	register uint64_t R = BIT_SET_27 ^ r0,R2= BIT_SET_27;
	if (ib)R <<= 32;
	else R2 <<= 32;
	for (uint32_t i = 0; i < nua; i++) {
		register uint64_t U = tua[i];
		if (R&U) continue;// hitting another cell in band B
		tuaold[nuaold++] = U;
		if (_popcnt64(R2&U) < 4) return 1; // will be subset for all 
	}
	for (uint32_t i = 0; i < nuab1b2; i++)
		if (!(R & tuab1b2[i])) tuaold[nuaold++] = tuab1b2[i];
	//cout<<Char2Xout(R) << "BuilOldUAs nuas=" << nuaold << endl;
	return 0;
}
int GENUAS_B12::CheckOld() {// ua 2x27  + 5 bit length
	uint64_t * t = tuaold;
	register uint64_t ua2x = ua & BIT_SET_2X;
	for (uint32_t iua = 0; iua < nuaold; iua++) {
		register uint64_t R = t[iua];
		R &= BIT_SET_2X;
		if ((R&ua2x) == R)		return 1;// we have a subset
	}
	return 0;
}

int GENUAS_B12::CheckMain(uint64_t wua) {//subset in the main table
	register uint64_t ua2x = wua & BIT_SET_2X;
	for (uint32_t iua = 0; iua < nua; iua++) {
		register uint64_t R = tua[iua];
		if (R > wua) return 0; // no subset
		if (R < wua) {// is it subset
			R &= BIT_SET_2X;
			if ((R&ua2x) == R) return 1;// subset discard
		}
		else if (R == ua) return 1;// same discard
	}
	return 0;
}
void GENUAS_B12::CollectMore2digits() {// special 6 7 digits minirow
	BF64 * mypm = zh2b_g.fd_sols[0];
	//zh1b_g.tua = tuamore;
	for (int dig1 = 0; dig1 < 8; dig1++) { 
		for (int dig2 = dig1 + 1; dig2 < 9; dig2++) { 			
			digp = (1 << dig1) | (1 << dig2);
			p12 = (mypm[dig1] | mypm[dig2]).bf.u64;
			uint32_t pbx[2],colsx[2],mx[2];
			pbx[0] = p12 & BIT_SET_27;
			pbx[1] = (p12>>32) & BIT_SET_27;
			{// build columns and mini rows
				register uint32_t px = pbx[0], R = 0x1ff;
				colsx[0] = (px & R) | ((px >> 9)&R) | ((px >> 18)&R);
				mx[0]= TblShrinkMask[px & R] |
					TblShrinkMask[(px >> 9) & R] << 3 |
					TblShrinkMask[(px >> 18) & R] << 6;
				//______________________
				px = pbx[1];
				colsx[1] = (px & R) | ((px >> 9)&R) | ((px >> 18)&R);
				mx[1]= TblShrinkMask[px & R] |
					TblShrinkMask[(px >> 9) & R] << 3 |
					TblShrinkMask[(px >> 18) & R] << 6;
			}
			for (ib = 0; ib < 2; ib++) {
				cola = colsx[1-ib];
				uint32_t mpbx=pbx[ib],mcolsx=colsx[ib],mmini=mx[ib],
					nmini = _popcnt32(mmini), ncols = _popcnt32(mcolsx);
				//cout << Char9out(digp) << " digits ib="<<ib
				//	<<" nmini="<<nmini<<" ncols="<<ncols<<endl;
				
				switch (nmini) {
				case 3: {// try each mini and each pair
					for (int irow = 0; irow < 3; irow++) {
						patb = mpbx & (0777 << (9 * irow));
						colb = patb>>(9*irow);// just catch the row
						Collect2digits2_4_cols();// first call one minirow
						colb ^= mcolsx ;
						patb ^= mpbx;
						Collect2digits2_4_cols();// second call 2 minirows
					}
				} //m case 3
					break;
				case 5: {// several patterns
					switch (ncols) {
					case 4: {// one mini pair (one in band)
						for (int ist = 0; ist < 3; ist++) {
							colb = mcolsx & (7 <<(3*ist) );
							if (_popcnt32(colb) != 2)continue;
							patb = mpbx & (07007007 << (3 * ist));
							Collect2digits2_4_cols();
						}
					}// c case 4
							break;
					case 5: {// one mini pair one 2 cols 4 cells
						for (int ist = 0; ist < 3; ist++) {
							uint32_t  minis = mmini & (0111 << ist);
							if (_popcnt32(minis) != 1)continue;
							patb = mpbx & (07007007 << (3 * ist));
							colb = mcolsx & (7 << (3 * ist));
							Collect2digits2_4_cols();// first call one minirow
							uint32_t rpatb = patb;
							for (int ist2 = 0; ist2 < 3; ist2++) {// second stack
								if (ist == ist2)continue;
								colb = mcolsx & (7 << (3 * ist2));
								if (_popcnt32(colb) != 2)continue;
								patb =rpatb ^ mpbx;// other cells of the band
								Collect2digits2_4_cols();// second  call 2 cols four cells
							}
						}
					}// c case 5
							break;
					case 6: {// one mini pair one 4 cols
						for (int ist = 0; ist < 3; ist++) {
							uint32_t minis = mmini & (0111 <<  ist);
							if (_popcnt32(minis) != 1)continue;
							patb = mpbx & (07007007 << (3 * ist));
							colb = mcolsx & (7 << (3 * ist));
							Collect2digits2_4_cols();// first call one minirow
							colb ^= mcolsx;
							patb ^= mpbx;
							Collect2digits2_4_cols();// second call 2 minirows
						}
						break;
					}//c case 6
					}// c switch
				}// m case 5
					break;
				case 6: {if (ncols != 4) continue;
					// one 2 columns 6 cells
					for (int ist = 0; ist < 3; ist++) {
						colb = mcolsx & (7 << (3 * ist));
						if (_popcnt32(colb) != 2)continue;
						patb = mpbx;
						Collect2digits2_4_cols();// first call one minirow
					}

					}//m case6
				}//m switch
			}// ib
		}// dig2
	}// dig1
	/* possible maps 2 digits one band
	   ab|  |		ab| |		ab| |		ab|  |		a | |b		 
	     |ab|		  |a|b		  |a|b		  |a |b		 b|a|		 
		 |  |ab		  |b|a		  |b| a		  | b| a	  |b|a	 
	  6 cols		4 cols		5 cols		6 cols		4 cols		 
	  3 minis		5 minis		5 minis		5 minis		6 minis		 
	  3(2) 3 (4)	1 (2)   	1(2) 1(2-4)	1(2) 1(4)	1 (2-6) 
	  a |  |b		a |  |b		a |c [b   ignore also
	   b|a |		 b|a |		 b|a |c   3 columns
	    | b|a		  | b| a      | b[ a   
	  5 cols		6 cols		note : ignored here
	  6 minis		6 mini		3 digits in band
	 0 case		0 case
	*/
}
void GENUAS_B12::Collect2digits2_4_cols() {
	if ((cola & colb) != colb) return; // must have the same columns valid
	STD_B1_2 * mybx[2] = { &myband1 ,&myband2 };
	ba = mybx[ib];
	bb = mybx[1 - ib];
	if (BuilOldUAs(patb)) return;

	if (0 &&_popcnt32(patb)>2) {
		cout << Char9out(digp) << " digits";
		cout << Char9out(colb) << " socket";
		cout << Char27out(patb) << " pat in B ib="<<ib << endl;
		cout << Char2Xout(p12) << " pat12 for the digits" << endl;
	}
	w0 = patb;// w0 is the ua part located in band ba
	if (ib) w0 <<= 32;
	// can cancel the search if a small old ua is granted subset
	zh1b_g.GetBand(bb->band0, tuamore);
	bb->InitRevisedg();
	for(int icol=0,bit=1;icol<9;icol++, bit<<=1)
		if(colb&bit)bb->revised_g[icol] ^= digp;
	CollectMoreTry6_7();// then go 6/7
}

void GENUAS_B12::CollectMoreTry6_7() {
	nfloors = 6;// for debugging

	//____________ try 6 digits unsolved in band a
	for (int i6 = 0; i6 < 84; i6++) {
		int fl3 = floors_3d[i6];// , fl6 = 0x1ff ^ fl3;
		if (fl3&digp) continue;// digits must be in the multi floors
		floors = 0x1ff ^ fl3;
		zhone[0].Start_Uas_Mini(floors, digp);
		zhone[0].ApplyGangsterChanges(bb->gangster, bb->revised_g);
		zh1b_g.nua = 0;
		zhone[0].InitGuess();
		if(zhone[0].Update6())
		zhone[0].Guess6();
		if (zh1b_g.nua) EndCollectMoreStep();
	}

	nfloors = 7;// for debugging

	//____________ try now 7 digits unsolved in band a
	for (int i7 = 0; i7 < 36; i7++) {
		int fl2 = floors_2d[i7];// , fl7 = 0x1ff ^ fl2;
		if (fl2&digp) continue;// digits must be in the multi floors
		floors = 0x1ff ^ fl2;
		zhone[0].Start_Uas_Mini(floors, digp);
		zhone[0].ApplyGangsterChanges(bb->gangster, bb->revised_g);
		zh1b_g.nua = 0;
		zhone[0].InitGuess();
		if (zhone[0].Update7())		zhone[0].Guess7();
		if (zh1b_g.nua) EndCollectMoreStep();
	}
}
void GENUAS_B12::EndCollectMoreStep() {
	for (uint32_t i = 0; i < zh1b_g.nua; i++) {
		ua = zh1b_g.tua[i] &= BIT_SET_27;
		if (!ib) ua <<= 32;
		ua |= w0;
		ua &= BIT_SET_2X;// be sure to use only relevant bits
		uint64_t cc = _popcnt64(ua);
		if (cc > limsize) continue;
		if (CheckOld()) continue;
		cc <<= 59;
		ua |= cc;
		AddUA64(tua, nua,ua);
	}
}
void GENUAS_B12::CollectTriplets() {// special 6 7 digits full minirow
	modemore = 3;
	STD_B1_2 * mybx[2] = { &myband1 ,&myband2 };
	for ( ib = 0; ib < 2; ib++) {
		//check possible  triplets  must be three possible digits
		ba = mybx[ib];
		bb = mybx[1 - ib];
		// init the brute force 
		zh1b_g.GetBand(bb->band0, tuamore);
		for (int imini = 0, cell = 0; imini < 9; imini++, cell += 3) {
			for (int ip = 0; ip < 2; ip++) {// 2 false triplets in mini row
				bb->InitRevisedg();
				digp = ba->ReviseG_triplet(imini, ip, bb);
				if (!digp)continue;// not a possible perm
				// ______________________need the cells to assign in band bb
				uint32_t R0 = 7 << (3 * imini);
				w0 = R0;
				if (ib) w0 <<= 32; // 2 cells each ua
				if (BuilOldUAs(R0)) continue;
				CollectMoreTry6_7();// then go 6/7
			}
		}
	}
}
void GENUAS_B12::CollectMore2minirows() {
	modemore = 4;
	STD_B1_2 * mybx[2] = { &myband1 ,&myband2 };
	for (ib = 0; ib < 2; ib++) {
		//if(ib)	zh1b_g.diag = 1;
		ba = mybx[ib];
		bb = mybx[1 - ib];
		zh1b_g.GetBand(bb->band0, tuamore);
		if (zh1b_g.diag) {
			cout << "new 2 minis ib=" << ib << endl;
			zhone[0].CheckSolPerDigit(); 
		}
		int npairs = ba->nvpairs;
		//cout << " CollectMore2minirows vpairs status nvpairs=" << npairs <<" ib="<<ib<< endl;
		uint32_t bcells1,bcells2;
		for (int i1 = 0; i1 < npairs-1; i1++) {
			int cell1 = ba->tv_pairs[i1],
				box1=cellsFixedData[cell1].eb;
			for (int i2 = i1+1; i2 < npairs ; i2++) {
				int cell2 = ba->tv_pairs[i2],
					box2 = cellsFixedData[cell2].eb;
				if (box1 == box2) continue;
				// 2 mini rrows
				bb->InitRevisedg();
				digp = ba->GetMiniData(cell1, bcells1, bb);
				digp |= ba->GetMiniData(cell2, bcells2, bb);
				if (_popcnt32(digp) < 3) continue; // already done
				uint32_t  R0 = bcells1 | bcells2;
				w0 = R0;// w0 is tha ua part located in band ba
				if (ib) w0 <<= 32; // 2 cells each ua
				if (BuilOldUAs(R0)) continue;
				CollectMoreTry6_7();// then go 6/7
			}
			//_______________________ get 1  pair + 1 triplet
			for (int ibox2 = 0; ibox2 < 3; ibox2++) {
				if (box1 == ibox2) continue;
				for (int iminirow = 0; iminirow < 3; iminirow++) {
					for (int ip = 0; ip < 2; ip++) {// 2 false triplets in mini row
						int imini = 3 * iminirow + ibox2;
						bb->InitRevisedg();
						digp = ba->ReviseG_triplet(imini, ip, bb);
						if (!digp)continue;// not a possible perm
						digp |= ba->GetMiniData(cell1, bcells1, bb);
						uint32_t R0 = 7 << (3 * imini) | bcells1;
						w0 = R0;// w0 is the ua part located in band ba
						if (ib) w0 <<= 32; // 2 cells each ua
						if(BuilOldUAs(R0)) continue;
						CollectMoreTry6_7();
					}
				}
			}
		}
	}
	
}
//=============== ua collector socket 2

void GENUAS_B12::ProcessSocket2(int i81) {
	GEN_BANDES_12::SGUA2 &wi81 = genb12.tsgua2[i81];
	zh2b_g.InitGangster(genb12.gangcols, wi81.gangcols);

}
