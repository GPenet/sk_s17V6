//#define DODEBUG
//#define DO2X2
#define ZST6 800
#define ZST5 500
#define ZST4 300
#define ZST3 50
#define ZST2 5
struct MINCOUNT {
	uint32_t mini_bf1, mini_bf2, mini_bf3, mini_triplet,
		critbf, pairs27;
	uint32_t all_used_minis, mincount, minplus;
	inline void SetMincount() {// after direct setting minis
		all_used_minis = mini_bf1 | mini_triplet;
		mini_triplet &= ~mini_bf1;// count only triplets with no pair
		mincount = _popcnt32(all_used_minis) + _popcnt32(mini_bf3);
		mini_bf1 &= ~mini_bf2;// now pure one pair
		mini_bf2 &= ~mini_bf3;// now pure 2 pairs

		// set up pair + triplet bitfield
		if (mini_triplet) {// must add triplet minirow
			for (int i = 0, bit = 1, field = 7; i < 9; i++, bit <<= 1, field <<= 3)
				if (mini_triplet&bit)
					critbf |= field;
		}
		minplus = mincount;
	}
	GINT64_t  Count_per_stack() {
		GINT64 cc; cc.u64 = 0;
		for (int i = 0, st = 0111; i < 3; i++, st <<= 1) {
			cc.u16[i] = _popcnt32(all_used_minis&st) +
				_popcnt32(mini_bf3&st);
		}
		return cc;
	}
	void Status(const char * lib) {
		cout << lib << "critical Status mincount =" << mincount << " minplus=" << minplus << endl;
		cout << Char27out(critbf) << " critical bf" << endl;
		cout << Char27out(pairs27) << " pairs 27" << endl;
		cout << Char9out(mini_bf1) << "     minis bf1" << endl;
		cout << Char9out(mini_bf2) << "     minis bf2" << endl;
		cout << Char9out(mini_bf3) << "     minis bf3" << endl;
		cout << Char9out(mini_triplet) << " mini triplets" << endl << endl;

	}
};
struct BI2 {//Index 2 valid in band
	uint64_t bf, // 2 cells in bit fiekd
		active; // remaining possible cells
	uint32_t tval[2], // 2 cells in int mode
		istart, iend;// in the VALIB table

};

struct VALIDB {//  valid band 1+2
	uint64_t bf; // 2 cells in bit fiekd
	uint32_t tval[4]; // 3-4 cells in int mode
	uint64_t nval;// n cells over the 2
	inline void Enter(uint64_t ebf, uint32_t * tc2) {
		bf = ebf;
		nval = _popcnt64(bf) - 2;
		memcpy(tval, tc2, sizeof tval);
	}
	int DebugFindKnown17();
	void Print() {// debugging
		cout << Char2Xout(bf) << " valid n=" << nval << endl;
	}
};
struct VALIDB64 {// minimal valid band no index
	uint64_t bf; // 2 cells in bit fiekd
	uint32_t tval[6]; // 0-6 cells in int mode
	uint64_t nval;// n cells over the 2
	void Print() { // debugging
		cout << Char2Xout(bf) << " valid n=" << nval << endl;
	}
};
struct CTR162 {// bit field 162 bits for guas2_3
	uint64_t bf[3];
	inline void New(CTR162 & o) {
		bf[0] &=  ~o.bf[0];
		bf[1] &=  ~o.bf[1];
		bf[2] &=  ~o.bf[2];
	}
	inline void Or(CTR162 & o) {
		bf[0] |= o.bf[0];
		bf[1] |= o.bf[1];
		bf[2] |= o.bf[2];
	}
	inline uint64_t NotEmpty() { return (bf[0] | bf[1] | bf[2]); }

}bf162all, bf162new, bf162more;

// standard first band (or unique band)
struct STD_B416 {
	BI2 * my_bi2;
	VALIDB * my_validb;
	uint32_t nbi2, nvalidb; 
	char band[28];
	int band0[27], i416, gangster[9], map[27], dband;
	uint32_t tua[100], nua;//   maximum 81  
	uint32_t fd_sols[2][9];//start puzzle/ solution
	void Initstd();
	void GetBandTable(int i);
	void SetGangster();
	inline void GetUAs() {
		nua = t16_nua[i416];
		memcpy(tua, &t16_UAs[t16_indua[i416]], 4 * nua);
	}
	void MorphUas()	;
	void InitC10(int i);
	void InitG12(int i);
	void InitBand2_3(int i16, char * ze, BANDMINLEX::PERM & p
		, int iband = 1);


	void PrintStatus();
};
struct STD_B1_2 :STD_B416 {
	// row solution pattern in digit
	int mini_digs[9], mini_pairs[27],
		revised_g[9];// after false forced in 1/2 minirows
	int  tv_pairs[27], nvpairs; //9 or 27 bits 
	void FillMiniDigsMiniPairs(STD_B1_2 & bb);
	inline void InitRevisedg() {
		memcpy(revised_g, gangster, sizeof gangster);
	}
	int ReviseG_triplet(int imini, int ip, STD_B1_2 * bb);
	uint32_t GetMiniData(int index, uint32_t & bcells, STD_B1_2 *bb);
	void PrintShortStatus();
};
struct STD_B3 :STD_B416 {// data specific to bands 3
	BF128 tua128[1000], tua128_b2[1000];
	struct GUAs {
		BF128 isguasocket2, isguasocket3, isguasocket2_46;// active i81
		BF128 isguasocketc2, isguasocketc3, isguasocketc2_46;// active i81
		int triplet[9];//same gua3s
		int triplet_imini[81];
		int ua_pair[81], ua_triplet[81]; // storing ua bitfields
		int ua2_imini[81], ua3_imini[81],
			ua2pair27[81], ua2bit[81], ua3bit[81];
	}guas;
	BF128 isguasocketc246;//all linked to a socket 2
	CTR162 bf162b,bf162b46;// g2 g3  in 162 bits mode 

	uint32_t ntua128,ntua128_b2;
	int minirows_bf[9];
	int triplet_perms[9][2][3];
	uint32_t i_27_to_81[27], i_9_to_81[9]; //band i81 for guas guas3
	//_______________ handling mincount
	MINCOUNT smin;
	GINT64  stack_count;

	//_______________________

	void InitBand3(int i16, char * ze, BANDMINLEX::PERM & p);
	void InitStack(int i16, int  * z0, BANDMINLEX::PERM & p, int iband);
	int IsGua(int i81);
	int IsGua3(int i81);
	void Setup_bf162b() {
		bf162b.bf[0] = guas.isguasocketc2.bf.u64[0];
		register uint64_t R= guas.isguasocketc2.bf.u64[1],
			R2= guas.isguasocketc3.bf.u64[0];
		bf162b.bf[1] = R | (R2 << 17);//17g2 47 g3
		R2 >>= 47; //now 17 g3
		R = guas.isguasocketc3.bf.u64[1];// last 17 g3
		bf162b.bf[2] = R2 | (R << 17);//17+17=34 g3
	}
	inline int GetI81_2(int bf) {
		for (uint32_t i = 0; i < 27; i++) {
			register uint32_t i81 = i_27_to_81[i];
			if (guas.ua_pair[i81] == bf) return i81;
		}
		return -1;
	}
	inline int GetI81_3(int bf) {
		for (uint32_t i = 0; i < 9; i++) {
			register uint32_t i81 = i_9_to_81[i];
			if (guas.ua_triplet[i81] == bf) return i81;
		}
		return -1;
	}
	int Clean0(); //first 256 guas 
	int Cleanmore(); // all guas
	inline void Insert2(uint32_t i81) {
		//tua2_3[ntua2_3++] = guas.ua_pair[i81];
		register uint32_t	bit = guas.ua2bit[i81];
		smin.mini_bf3 |= smin.mini_bf2&bit;
		smin.mini_bf2 |= smin.mini_bf1&bit;
		smin.mini_bf1 |= bit;
		smin.critbf |= guas.ua_pair[i81];
		smin.pairs27 |= guas.ua2pair27[i81];
	}
	inline void Insert3(uint32_t i81) {
		//tua2_3[ntua2_3++] = guas.ua_triplet[i81];
		smin.mini_triplet |= guas.ua3bit[i81];
	}
	void PrintB3Status();
};

struct BINDEXN {
	BI2 * t2;
	VALIDB * tvb;
	uint32_t nt2, ntvb;
	inline void Attach(BI2 * t2e, VALIDB * tvbe) { t2 = t2e; tvb = tvbe; }
	void Copy(STD_B1_2 & b);
	void Copy(BINDEXN & b);
}bin_b1, bin_b2, bin_b1yes, bin_b2yes,
bin2_b1, bin2_b2, bin2_b1yes, bin2_b2yes;
struct INDEX_XY {
	uint64_t bf, and_g, or_g;// 2 cells common to the lot
	uint64_t ntotvb, ncluesmin;
	struct ITEM {// one of the tables per lot 
		VALIDB64 * tvb;
		uint32_t ntvb, sum_vb;
	}titem[2];
	void Debug() {
		cout << Char2Xout(and_g) << "\tnand=" << _popcnt64(and_g) << endl;;
		cout << Char2Xout(or_g) << "\tnor=" << _popcnt64(or_g) << endl;;
		for (int i = 0; i < 2; i++)
			cout << "\t" << titem[i].ntvb << ";" << titem[i].sum_vb;
	}

}index_xy_b1, index_xy_b2;

BF128 uas_buffer[20000], uasnn_buffer[30000];

struct G17TMORE {// FIFO table of more for bands 1+2
	uint64_t  t[G17MORESIZE];
	int nt, maxt, curt;
	inline void Init() { maxt = G17MORESIZE; nt = 0; }
	inline void Add(uint64_t v) {//add a new more in FIFO 
		if (nt < maxt) {// use next location
			curt = nt;
			t[nt++] = v;
		}
		else {// replace the oldest
			curt++;
			if (curt == maxt)curt = 0;
			t[curt] = v;
		}

	}
	inline void Add_If_New(uint64_t v) {// check oldest first
		register uint64_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// from curt to 0
		if ((*Rt) == V)return;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) goto add;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt > Rtl)if ((*Rt) == V)return;
	add:
		Add(v);
	}
	inline int Check(uint64_t v) {// check oldest first
		if (!nt) return 0;
		register uint64_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		if (!((*Rt) & V))return 1;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return 0;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt > Rtl)if (!((*Rt) & V))return 1;
		return 0;
	}
	void Print(int modegua) {
		register uint64_t *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		{
			register uint64_t w = *Rt;
			if (modegua)			cout << Char54out(w) << " " << (w >> 56) << endl;
			else cout << Char2Xout(w) << endl;
		}
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return;

		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt > Rtl) {
			register uint64_t w = *Rt;
			if (modegua)			cout << Char54out(w) << " " << (w >> 56) << endl;
			else cout << Char2Xout(w) << endl;
		}

	}


	void PrintUasDirect() {
		for (int i = 0; i < nt; i++) {
			register uint64_t w = t[i];
			cout << Char2Xout(w) << endl;
		}
	}

};
struct MORE32 {// FIFO table of more for band b
	uint32_t  t[32];
	int nt, maxt, curt;
	inline void Init() { maxt = 32; nt = 0; }
	inline void Add(uint32_t v) {//add a new more in FIFO 
		if (nt < maxt) {// use next location
			curt = nt;
			t[nt++] = v;
		}
		else {// replace the oldest
			curt++;
			if (curt == maxt)curt = 0;
			t[curt] = v;
		}

	}

	inline int Check(uint32_t v) {// check oldest first
		if (!nt) return 0;
		register uint32_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		if (!((*Rt) & V))return 1;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return 0;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt > Rtl)if (!((*Rt) & V))return 1;
		return 0;
	}

};
struct MOREALL {// FIFO table of more for bands 1+2
	BF128 ta[G17MORESIZE];
	uint64_t  t[G17MORESIZE];
	int nt, maxt, curt;
	inline void Init() { maxt = G17MORESIZE; nt = 0; }
	inline void Add(BF128 v) {//add a new more in FIFO 
		if (nt < maxt) {// use next location
			curt = nt;
			ta[nt++] = v;
		}
		else {// replace the oldest
			curt++;
			if (curt == maxt)curt = 0;
			ta[curt] = v;
		}

	}
	inline int Check(uint64_t v) {// check oldest first
		if (!nt) return 0;
		register uint64_t V = v;
		register BF128 *Rt = &ta[curt], *Rtl;
	loop1:// form curt to 0
		if (!(Rt->bf.u64[0] & V))return 1;
		if (--Rt >= ta)goto loop1;
		if (nt < maxt) return 0;
		Rtl = &ta[curt];
		Rt = &ta[maxt];
		while (--Rt > Rtl)if (!(Rt->bf.u64[0] & V))return 1;
		return 0;
	}
	void PrintUasDirect() {
		for (int i = 0; i < nt; i++) {
			BF128 wa = ta[i];
			cout << Char27out(wa.bf.u32[2]) << "\t";
			cout << Char2Xout(wa.bf.u64[0]) << endl;
		}
	}
};

struct GVCELLS {// vectors cells for guas
	BF128 v21[81][54]; // first 128 uas GUA2
	BF128 v31[81][54]; // first 128 uas GUA3
}gvcells;
struct GVSTEPS {// base vectors start,loopb1,loopb2
	BF128 v21[81];
	BF128 v31[81];
}gvs_start, gvs_b1, gvs_b2;
struct G4TABLE {
	struct TI {
		uint64_t ua64;
		uint32_t i81_1, i81_2;
	}ti[2000];;
	BF128 tua4[2000];
	uint32_t ntua4;
	inline void Init() { ntua4 = 0; }
	inline void Add(TI & a) {
		if (ntua4 < 2000)ti[ntua4++] = a;
	}
	void Shrink(G4TABLE & o, uint64_t bf) {
		ntua4 = 0;
		for (uint32_t i = 0; i < o.ntua4; i++)
			if (!(bf & o.tua4[i].bf.u64[0]))
				tua4[ntua4++] = o.tua4[i];
	}
	void Catch(uint64_t bf, uint32_t * td, uint32_t & ntd) {
		for (uint32_t i = 0; i < ntua4; i++)
			if (!(bf & tua4[i].bf.u64[0]))
				td[ntd++] = tua4[i].bf.u32[2];
	}
	void Debug(int nodet = 1) {
		cout << "g4t_.. status ntua4=" << ntua4 << endl;
		if (nodet)return;
		for (uint32_t i = 0; i < ntua4; i++) {
			cout << Char27out(tua4[i].bf.u32[2]) << "\t";
			cout << Char2Xout(tua4[i].bf.u64[0]) << " ";
			cout << _popcnt64(tua4[i].bf.u64[0]) << endl;
		}
	}

}g4t_start, g4t_b1, g4t_b2;

//================== UA collector 2 bands 

struct GENUAS_B12 {// for uas collection in bands 1+2 using brute force 
	int dig_cells[9][9],
		gangbf[9],// columns for band 3 in bit field
revised_gangbf[9],// same revised UA2s UA3s ***
mini_digs[9], mini_pairs[27], // UA2s UA3  ***
//valid_pairs, //  27 bits valid sockets UA2s ***
nfloors, limstep, map[9], cptdebug, modemore;
BF128 valid_sockets;

//=============== uas collector 
int limsize, floors;
uint64_t  tuaold[1000],// previous non hit uas infinal table of uas for bands 1+2
tua[TUA64_12SIZE]// 
, tuab1b2[200];// collecting bands uas in 2x mode
uint32_t nuaold, nua, nuab1b2,
tuamore[500];
//_______________uamore control
STD_B1_2 *ba, *bb;
uint32_t patb, ib, digp, colb, cola;
uint64_t w0, ua, p12;
//_____________________ functions collect UAs bands 1+2
int Initgen();
void BuildFloorsAndCollectOlds(int fl);
//int AddUA64(uint64_t * t, uint32_t & nt);
inline void AddUA(uint64_t v) {
	ua = v; AddUA64(tua, nua, ua);
}
inline void AddUACheck(uint64_t v) {
	if (nua >= TUA64_12SIZE) nua = TUA64_12SIZE - 1;
	ua = v; AddUA64(tua, nua, ua);
}
int BuilOldUAs(uint32_t r0);
int CheckOld();
int CheckMain(uint64_t wua);
void CollectMore2digits();
void Collect2digits2_4_cols();
void CollectMoreTry6_7();
void EndCollectMoreStep();
void CollectTriplets();
void CollectMore2minirows();
//_____________________ functions collect UA2s UAs3 socket 

void ProcessSocket2(int i81);
int DebugUas();
};
// gaus table size 64 but 30 as initial limit
#define SIZETGUA 35
//#define GUAREDSIZE 100
struct GEN_BANDES_12 {// encapsulating global data 
	STD_B3 bands3[512];
	int modeb12, go_back, diagmore, diagbug, ip20,
		it16, it16_2, imin16_1, imin16_2, imin16_3;
	int i1t16, i2t16, i3t16, maxnb3; // index 416 ordered in increasing size of valid clues 6
	char zsol[82], rband2[28];
	int grid0[82], tc[6], ntc;
	int gcheck[82], ib2check, ib3check;
	int skip, last;// restart point; last entry in the batch
	uint64_t   nb12;
	BANDMINLEX::PERM t_auto_b1[108], // maxi is 107excluding start position
		t_auto_b1b2[108], t_auto_b2b1[108],
		pband2, pband3, pcheck2, pcheck3;
	int n_auto_b1, n_auto_b1b2, n_auto_b2b1;
	int cold[9], coldf[9], rowd[6], boxd[6], rowdb3[3], boxdb3[3]; //free digits 
	//_________________ gangster 
	int gangcols[9];// 9 bits fields band3 for each colum (see valid band2)
	int gang[9][3]; // gangcols expanded (buildgang ) 3 digits
	int gangb12[9]; // digit bf bands 12 per column
	int   *gang27; // redefines gang[9][3] as 27 integer
	int   gang_digits_cols[9][3];// active cols for a given digit
	//____________structs hosting the 81 GUA entries
	struct SGUA2 {// 81 possible UA2 sockets
		// permanent data
		uint64_t * tua;
		int col1, col2;// columns of the socket
		int i_81; // index 0_80 for this 
		int i9;// complementary column in minirow
		int id1, id2; // index of digits in gang 27 
		// Current band1+2 data
		int digs, dig1, dig2;// depending on gang27 status
		int valid, // valid if guas 
			validuas,// gua2s found
			used;// if needed in bands3
		int gangcols[9];// revised gangster
		uint32_t nua;// nua_start, nua_end;
		void Debug(const char * lib);

	}tsgua2[81];

	int GET_I81_G2(int digs, uint32_t pat) {
		register  uint32_t A = pat,
			B = (A | (A >> 9) | (A >> 18)) & 0777; // all columns
		uint32_t dstack;
		bitscanforward(dstack, B);
		dstack = (dstack / 3) * 27;// now stack0 0, 27 57 chunks of 27 i81
		for (uint32_t i = dstack; i < dstack + 27; i++)
			if (digs == tsgua2[i].digs) return i;
		return 0; // would be a bug
	}
	int GET_I81_G2_4(int digs, uint32_t pat) {
		register  uint32_t A = pat,
			B = (A | (A >> 9) | (A >> 18)) & 0777; // all columns
		// keep the stack with 2 columns
		for (int i = 0, mask = 7; i < 3; i++, mask <<= 3) {
			if (_popcnt32(B&mask) == 2) {
				B = mask;
				break;
			}
		}
		uint32_t dstack; // now same as pat 2 cells
		bitscanforward(dstack, B);
		dstack = (dstack / 3) * 27;// now stack0 0, 27 57 chunks of 27 i81
		for (uint32_t i = dstack; i < dstack + 27; i++)
			if (digs == tsgua2[i].digs) return i;
		return 0; // would be a bug
	}
	struct SGUA3 {// 81 possible UA3 sockets
		// permanent data
		uint64_t * tua;// , killer;
		int col1;// first columns 0-9 
		int i_81, imini;// , iguan; // index 0_80 for this 
		int id1, id2, id3; // index of digits in gang 27 
		// Current band1+2 data
		int  dig1, dig2, dig3,digs;// depending on gang27 status
		int valid, // valid if guas 
			validuas,// gua2s found
			used;// if needed in bands3
		uint32_t nua;// nua_start, nua_end, nua;
		void Debug(const char * lib);
	}tsgua3[81];

	int GET_I81_G3(int digs, uint32_t pat) {
		register  uint32_t A = pat,
			B = (A | (A >> 9) | (A >> 18)) & 0777; // all columns
		uint32_t dstack;
		bitscanforward(dstack, B);
		dstack = (dstack / 3) * 27;// now stack0 0, 27 57 chunks of 27 i81
		for (uint32_t i = dstack; i < dstack + 27; i++)
			if (digs == tsgua3[i].digs) return i;
		return 0; // would be a bug
	}

	// __________________________  primary UAs tables and creation of such tables
	uint64_t  *ptua2;// pointer to current table cycle search 2/3
	uint32_t   nua2; // nua2 for cycle search 2/3  
	//================== bands 3 and gangster band 3 analysis
	int nband3;
	int   tcolok[2], ncolok;

	int ngua6_7, c1, c2, band, floors, digp, i81;
	uint64_t wua0, ua;// partial gua ua to check
	uint64_t tuacheck[100], tua_for_check[500];
	uint32_t uadigs_for_check[500], nua_for_check, nua_check;
	//================ A creating a catalogue for the 17 search 
	//sorted increasing number of valid bands 6 clues

	GEN_BANDES_12() {
		gang27 = gang[0];
		InitialSockets2Setup();
		InitialSockets3Setup();
	}
	void InitialSockets2Setup();// batch level
	void InitialSockets3Setup();// batch level
	void BuildGang9x3();
	void Build_CheckUAs_Subsets_Table();
	void Build_tuacheck(int fl);
	int Have_tuacheck_subset();
	void SecondSockets2Setup();// band1+2 level
	void SecondSockets2MoreUAs();// band1+2 level
	void GuaCollectMore();
	void SecondSockets3Setup();// band1+2 level
	void GuaCollect(int fl, int diag = 0);
	//================================= functions
	void GetStartB2(int i); // one of the 20 starts 
	void Start(int mode = 0);
	void NewBand1(int iw);
	int Band2Check();
	int Band3Check();
	void Find_band2B();
	int ValidBand2();
	void ValidInitGang();
	void Find_band3B(int m10 = 1);
	int DebugFreshUA(uint64_t ua);
	int Debug17(SGUA2 & w);
	//int FindBand3Unique();//test or  debugging code see the corresponding file
	//================ B creating a catalogue for the 17 search 
	//same as A exchanging bands 2/3


	//============= loops control for UAs 5;6;7 digits collection (collect more=
	int iband, ibox, iminirow, ibox2, iminirow2, pat1, pat2, ncells;
	int tcells[6], tcols[6];
	int bcols[2][9], mycols[9], myfloors;
	uint64_t mybf;
	// debugging code special call
	//int Afterb1b2(int option = 0);
};

/*
entry 92
maxindex= 983
maxn5= 51516
maxn6= 237762
maxdet5= 261
maxdet6= 2004
*/

/*
a GUA can be 
type 
0 one of the 81 gangster pairs
1 one of the 81 gangster triplets

for a given band we can have a BF128
UA produced in band 3 process not gua2 gua3
a stack ua for the band 
could be also expansion to 6 cells of the type 0

in bit fields "active" 
all 81 types are encoded in 3x27 mode in a 128 bits field



At the start, a working GUA struct is open
if less than 23uas are collected, the correponding uas 
 are stored in a separate table of GUAR.

for each 2 clues step in band 1 or band 2 
the GUARs table is first reduced to still active GUARs
then the GUAs are reduced and if the count is <3 
 the result is sent in the GUAR table


the search is done directly in the GUAs GUARs tables, 
but an index of active gangsters is setup at each step.

For each bands 1+2 passing the filters, 
 the final status of the gangsters sockets is built in the same way
 for a GUA, the first not hit forces the ganster true,
 for a GUAR, the same happens, but the process continues with the next entry

*/

struct GUA {
	uint64_t tua[64];
	uint32_t nua,type, i81 ;
	inline void Add(uint64_t ua) {
		if (nua < 64)tua[nua++] = ua;
	}
	inline void Adduacheck(uint64_t ua) {
		if (nua > 63) nua = 63;
		AddUA64(tua, nua, ua);
	}
	inline void Init(uint32_t n, uint32_t t, uint32_t i) {
		nua = n; type = t; i81 = i; ;
	}
	inline void Init(GUA & g) {
		nua = 0; type = g.type; i81 = g.i81; 
	}
	void Debug(int nodet=1) {
		cout << "gua type=" << type << " i81=" 	<< i81
			<< "\tnua=" << nua << endl;
		if (nodet) return;
		for (uint32_t i = 0; i < nua; i++) {
			cout << "\t" << Char2Xout(tua[i]) 
				<<" "<< _popcnt64(tua[i]& BIT_SET_2X)
				<< " " << i << endl;
		}
	}
}guaw;
struct GUARx {
	uint64_t ua;
	uint32_t type, i81;
	void Init(GUA & o) {
		type = o.type;
		i81 = o.i81;
	}
	inline void Enter(uint64_t eua, uint32_t etype, uint32_t ei81) {
		i81 = ei81; type = etype; ua = eua;
	}
	void Debug() {
		cout << "gr t=" << type << " i81=" << i81
			 << "\t" << Char2Xout(ua)
			<<" "<<_popcnt64(ua &BIT_SET_2X)<< endl;
		
	}
}guarw;

struct TGUAS {
	struct INDEX_TO_GUA {
		int g2[81], g3[81];
		inline void Set(uint32_t i_gua, GUA & gua) {
			// set the relevant index pointer  to the gua
			if (gua.type) 		g3[gua.i81] = i_gua;
			else g2[gua.i81] = i_gua;  
		}
	}ix_start;
	struct CTRL {// check for redundancy
		BF128 g2, g3;
		void Debug() {	g2.Print("g2");	g3.Print("g3");	}
	}ctrla;


	GUA tgua_start[162];// max is 81+81
	//GUAR tguar_start[2000];// , tguar_b1[2000], tguar_b2[4096]; // 4096=128*32
	uint32_t ntgua_start, 	nvect;

	uint32_t  tvti[4096]; // index 0_161 of the vector item
	BF128 vect[32], cellsv[54][32],vect_b2[32], vect_xy[32],
		cur_vect;
	void InitStart() {
		memset(&ix_start, 255, sizeof ix_start);
		ntgua_start = 0;
	}
	inline void InitApply() {
		memset(&bf162all, 0, sizeof bf162all);
	}
	inline void AddStart(GUA & g) {
		ix_start.Set(ntgua_start, g);
		tgua_start[ntgua_start++] = guaw;
	}
	inline void SetVect54(uint64_t ua, uint32_t i) {
		long long * w = (long long*)&vect[0].bf.u64[0];
		_bittestandset64(w, i);  // set the bit to 1
		register uint64_t W = ua ;
		uint32_t cc54;
		while (bitscanforward64(cc54, W)) {// look for  possible cells
			W ^= (uint64_t)1 << cc54;// clear bit
			tguas.cellsv[cc54][0].Clear(i);
		}
	}
	inline void SetVect(uint64_t ua, uint32_t i) {
		long long * w = (long long*)&vect[0].bf.u64[0];
		_bittestandset64(w, i);  // set the bit to 1
		register uint64_t W = ua & BIT_SET_2X;
		uint32_t cc64;
		while (bitscanforward64(cc64, W)) {// look for  possible cells
			W ^= (uint64_t)1 << cc64;// clear bit
			tguas.cellsv[From_128_To_81[cc64]][0].Clear(i);
		}
	}
	inline void AddVect(uint64_t eua, uint32_t etype, uint32_t ei81) {
		if (nvect > 4090) return; // limit for vectors is 4096
		SetVect(eua, nvect);
		tvti[nvect++] = ei81+81* etype;
	}
	inline uint32_t Get_nblocs128() {
		return (nvect + 127) / 128;
	}
	void ApplyB2();
	void ApplyXY(uint64_t bf);
	void Apply128(uint32_t i);
	void ApplyFirst128();
	void ApplyFirst256();

	void DebugStart(int nodet = 1) {
		cout << "gua tables n=" << ntgua_start << endl;
		for (uint32_t i = 0; i < ntgua_start; i++)
			tgua_start[i].Debug(nodet);
	}


	void Debugxy() {	ctrla.Debug();	}
}tguas;

/*
struct GUAN {// one occurrence of the active guan
	uint64_t *tua, killer;
	uint32_t colbf, digsbf, ncol, // 2-4 cols
		nua, nfree;
	int i81;// i81 or pattern (if 4 columns)
	void Enter(uint64_t *t, uint32_t n, uint32_t cbf,
		int32_t dbf, uint32_t ind) {
		i81 = ind;// index 0-80 if gua2 gua3
		tua = t; nua = n; colbf = cbf; digsbf = dbf;
		killer = BIT_SET_2X;
		for (uint32_t i = 0; i < nua; i++)killer &= tua[i];
		ncol = _popcnt32(colbf);
	}
	void Debug1Guan(int i, int all = 0) {
		cout << Char2Xout(killer);
		cout << "kill  i=" << i << "\tcols" << Char9out(colbf);
		cout << "\tdigs" << Char9out(digsbf) << " ncol=" << ncol
			<< " nua=" << nua << " i81=" << i81
			<< endl;
		if (all) {
			for (uint32_t i = 0; i < nua; i++)
				cout << Char2Xout(tua[i]) << endl;
		}
	}

}guanw;


struct TU_GUAN {// GUAN process (used GUA all kinds)
	GUAN tguan[256], guanw, tguanb1[160],g2[160];
	uint64_t  guabuf[15000], *pguabuf;
	uint32_t //i81socket2[80], i81socket3[80],
		tsockets2[80], tsockets3[80],// tsockets4[80],
		ntsockets2, ntsockets3, //ntsockets4,
		ntsockets2_2, ntsockets3_2;// , ntsockets4_2;
	//______________________ valid handler
	//uint64_t vv2, vv3, vv2_2, vv3_2;
	uint32_t nguan, nguanb1, nguastepb2, ng2 ;
	void AddGuan(uint64_t *t, uint32_t n, uint32_t cbf,
		int32_t dbf, uint32_t ind) {
		if (nguan < 256) {
			pguabuf = &pguabuf[n];// lock space
			tguan[nguan++].Enter(t, n, cbf, dbf, ind);
		}
		//tguan[nguan - 1].Debug1Guan(nguan - 1);
	}

	void Init(){
		pguabuf = guabuf;// reinit gua buffer use
		nguan = 0;//and guan table
	}


};*/



struct G17B3HANDLER {
	MINCOUNT smin;
	int known_b3, rknown_b3, active_b3, ib3, nb3,
		active_sub, ndead, wactive0, nmiss, //ncritical,
		irloop, wua, stack;
	uint32_t *uasb3if, nuasb3if, *uasb3of, nuasb3of, andoutf;
	GINT64 stack_count;
	int diagh;
	// ================== entry in the proces
	void GoMiss0(STD_B3 & b3);
	void GoMiss1(STD_B3 & b3);
	void Do_miss1();
	void GoMiss2Init(STD_B3 & b3);
	void GoMiss2(STD_B3 & b3, uint32_t uamin);
	inline void AddCellMiss1(uint32_t cell, int bit) {
		stack_count.u16[C_stack[cell]]++;
		known_b3 |= bit;
	}
	inline void AddCell_Miss2(uint32_t * t);
	inline int AddCell_Of(uint32_t cell, int bit) {
		register int s = C_stack[cell];
		if (stack_count.u16[s] > 5) return 0;
		stack_count.u16[s]++;
		if (stack_count.u16[s] > 5) {
			s = ~(07007007 << (3 * s));// mask
			wua &= s;
			wactive0 &= s;
		}
		nmiss--;
		known_b3 |= bit;
		return 1;
	}
	uint32_t IsMultiple(int bf);
	int ShrinkUas1();
	//=============== process critical
	void CriticalAssignCell(int Ru);
	void Critical2pairs();
	void Go_Critical();
	void CriticalLoop();
	//==================== process subcritical no cell added outside the GUAs field
	void SubMini(int M, int mask);
	void Go_Subcritical();
	void Go_SubcriticalMiniRow();

};
struct ZS64 {// final loop 64 bits 
	uint64_t bf,  v;
};


struct G17B {// hosting the search in 6 6 5 mode combining bands solutions
	G17B();// initial tasks all commands

	BF128 p17diag;// known 17 pattern for tests
	int b3lim, debug17,debug17_check,
		diag, diagbug,debugb3, aigstop, aigstopxy,
		iretb1,doloopnotok,
		npuz, a_17_found_here;
	uint32_t	iband1,iband2, step1count;
	G17B3HANDLER hh0;
	//______sockets common to  all bands 3  
	BF128 isguasocket2all, isguasocket3all;
	BF128 socks2x2;// sockets 2x2  4 cells  2 digits in band3
	int ts2x2[81], nts2x2; // storing active sockets 2x2
	uint64_t tuas2x2[81]; // first ua (usually one) for an active socket 2x2
	int ts2x2_clean[81], nts2x2_clean; // same status for a given XY 

	//====== data for band expansion
	uint32_t nexp, bnua;
	uint64_t btua[300],start_active, b1cpt[8], b2cpt[8],b1cptdiag;
	BI2 * mybi2t,wbi_1,wbi_2;
	VALIDB * myokt, validb_known,validb_known_b1, validb_known_b2;
	uint32_t nmybi2t, nmyokt,nbi2_1,nbi2_2,
		nvb1,nvb1steps[10],
		nzs1_6,nzs2_6, nzs1_5, nzs2_5;
	//======= status after step 2 in band 2 then 2 uas in band 1
	uint64_t tusb1[2000], tusb2_12[1000], tusb2[1000], temptyb1[500];
	uint32_t ntusb1, ntusb2, ntusb2_12, nemptyb1;
	VALIDB64 * tvb1go;// part A or part B
	uint32_t ntvb1go;
	uint64_t  temptyb2[500];
	uint32_t  ntua_128, ntua_256, nemptyb2;
	uint64_t fb1,acb1, fb2, fb12, acb2a,  acb2, acb12;


	//_____ data and functions of the external loop

	int loopb1;
	uint64_t breakcountx, b2count, totb2;
	struct EXTL {
		uint64_t noxyes;
		uint32_t bfx, tbfy[20], ntbfy, mode, ratio;
		void Init(uint32_t bf, int mod) { bfx = bf; ntbfy = 0; mode = mod; }
		void Debug() {
			uint32_t orw = 0;
			if (mode == 1) {
				cout << Char27out(bfx) << " bfx b1 ntbfy= " << ntbfy << endl;
				for (uint32_t i1 = 0; i1 < ntbfy; i1++) {
					orw |= tbfy[i1];
					cout << "\t\t" << Char27out(tbfy[i1]) << " bfy" << endl;
				}
				cout << "\t\t" << Char27out(orw) << " orw" << endl;
			}
			else {
				cout << "\t\t" << Char27out(bfx) << " bfx b2 ntbfy= " << ntbfy << endl;
				for (uint32_t i1 = 0; i1 < ntbfy; i1++) {
					cout << Char27out(tbfy[i1]) << " bfy" << endl;
					orw |= tbfy[i1];
				}
				cout << Char27out(orw) << " orw" << endl;

			}

		}
	}extl1[10], extl2[10], extlw, extlr;
	uint32_t nextl1, nextl2;
	uint64_t n_yesb1, n_yesb2, n_nob1, n_nob2,minratio;
	uint64_t FindSockets(uint64_t active, uint64_t lim);
	void ExtractMin(uint64_t active, BINDEXN & bin1, BINDEXN & bin2);
	void ExtSplitY(BINDEXN & binw, uint32_t *tbf, uint32_t ntbf,
		uint32_t & activer,int bande=1);
	void ExtSplitX(BINDEXN & bin1no, BINDEXN & bin1yes,
		uint32_t bf, uint32_t & activer,int bande=1);
	void Go2_Ext_Loop();
	void Go2b_Ext_Loop(uint64_t activeloop, uint32_t mode2);

	//___ process sub lots 
	uint32_t kn_ir1, kn_ir2; // to ckeck known
	int G3_SplitBi2( BINDEXN & binw, uint32_t ibi2,
		INDEX_XY & ixyw, VALIDB64 * pvb);
	
	void Go3(BINDEXN & bin1, BINDEXN & bin2);



	//============ vectors 64 128 bits 
	uint64_t v64uas, vc64[54];
	BF128 v64_192uas, vc64_192[54], bv192;
	BF128 v192_320uas, vc192_320[54], bv320;
	BF128 v320_448uas, vc320_448[54], bv448;
	BF128 v448_576uas, vc448_576[54], bv576;
	BF128 v576_704uas, vc576_704[54], bv704;
	BF128 v704_832uas, vc704_832[54], bv832;
	BF128 v832_960uas, vc832_960[54], bv960;
	uint64_t  n_to_clean;
	ZS64 zs64b1[MAX_56], zs64b2[MAX_56];

	//============================ b12 no more uas to test
	G17TMORE moreuas_12_13, moreuas_14, moreuas_15,
		 moreuas_AB_small, moreuas_AB, moreuas_AB_big;
	uint64_t wb12bf, wb12active,myua;
	GINT64  stack_count_step, stack_count, stack_countf;
	uint32_t tclues[40],tb3[256],*tcluesxy; 
	int nclues_step, nclues,ntb3;
	BF128 bands_active_pairs, bands_active_triplets,
		valid_vect;
	//============  go band3
	STD_B3* myband3;
	uint32_t fstk, andmiss1, noutmiss1, wactive0;
	uint32_t free1, free2, free3;
	uint32_t tua128_b3[1000], ntua128_b3;
	BF128 wg46;
	uint32_t cur_ib;
	uint32_t tcluesb12[20], ncluesb3x;
	uint32_t   nmiss;
	uint32_t uasb3_1[2000], uasb3_2[2000], uas_in[2000],
		nuasb3_1, nuasb3_2, nuas_in, b3_andout;
	MINCOUNT smin;
	MORE32 moreuas_b3, moreuas_b3_small;

	
	inline void AddB3_2_to_B3_1() {
		memcpy(&uasb3_1[nuasb3_1], uasb3_2, nuasb3_2 * 4);
		nuasb3_1 += nuasb3_2;
	}
	
	
	
	
	//=====================process for a new band 2 / set of bands 3
	void GoM10();// standard entry
	void GoM10Known();// entry for a known 17
	void GoM10Uas();// collect uas guas
	void StackUas();
	void ExpandB1();
	void ExpandB2();
	void ExpandB3();
	void ExpandOneBand(int ib);// find 2-6 valid bands no redundant clue
	void Apply_Band1_Step();// shrink tables
	int Apply_Band2_Step();// shrink tables again
			//________ end of step band 1

					  
					  //___________ outer loop on steps band2 band1

		//___________ extract potential valid bands 1+2 (no more uas)
	void Do64uas();
	void DoChunk64(ZS64 * a, ZS64 * b, uint64_t na, uint64_t nb);
	void Do64uas_11(ZS64 * a, ZS64 * b, uint64_t na, uint64_t nb);


	//_______ processing potential valid bands 1+2
	void CleanAll();
	inline void AddXClue(uint32_t * t, int & n, uint32_t xcell) {
		uint32_t cell = From_128_To_81[xcell];
		stack_count.u16[C_stack[cell]]++;
		t[n++] = cell;
	}
	int Is_B12_Not_Unique();
	void NewUaB12();
	void DebugAdd12();
	void GoB3(STD_B3 & b);
	void GoB3Expand();//mincount too low
	void GoB3_6();//mincount 6 no room for outfield
	void GoB3_5();//mincount 5 max 1 outfield
	void GoB3_4();//mincount 4 max 2 outfield
	void GoB3_4Expand();//mincount 4 max 2 outfield

	void FinalCheckB3(uint32_t bfb3);
	void Out17(uint32_t bfb3);
	void NewUaB3();
	void NewUaB3_g2(uint32_t i81, uint64_t ua12);
	void NewUaB3_g3(uint32_t i81, uint64_t ua12);

	void Debug_If_Of_b3();




	//================ debugging code
	inline uint32_t k17x(int ix) {	return p17diag.bf.u32[ix];	}
	uint64_t DebugFindKnown17_valids_2bands();
	void Debug_b1b2cpt();
	void DebugGetPuz(const char * p) {
		p17diag.SetAll_0();
		for (int i = 0; i < 81; i++)
			if (p[i] != '.')p17diag.Set_c(i);

		cout <<"this is a "
			<<_popcnt32(p17diag.bf.u32[0])
			<< _popcnt32(p17diag.bf.u32[1])
			<< _popcnt32(p17diag.bf.u32[2]) 
			<<" pattern for the expected puzzle"<< endl;
	}
	void GodebugInit(int mode);
	int GodebugCheckUas(const char * lib);
};