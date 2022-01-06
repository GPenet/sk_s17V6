
// initial taks all commands
G17B::G17B() {
	bin_b1.Attach(bi2_b1w, vab1w);
	bin_b2.Attach(bi2_b2w, vab2w);
	bin_b1yes.Attach(bi2_b1yes, vab1yes);
	bin_b2yes.Attach(bi2_b2yes, vab2yes);
	bin2_b1.Attach(bi2_b1w2, vab1w2);
	bin2_b2.Attach(bi2_b2w2, vab2w2);
	bin2_b1yes.Attach(bi2_b1yes2, vab1yes2);
	bin2_b2yes.Attach(bi2_b2yes2, vab2yes2);
	aigstop = 0;
}
void G17B::GoM10() {// processing an entry 656 566 with the relevant table of ba,ds3
	if (aigstop)return;
	p_cpt2g[0] ++;
	p_cpt2g[1] += genb12.nband3;
	ExpandB2();// expand band2 272 ms
	if (!myband2.nvalidb)return ; // mode 656 only no 5
	b1cptdiag = b1cpt[1];
	GoM10Uas();//expand bands 3  collect UAs 350 ms
	StackUas();
	bin_b1.Copy(myband1);
	bin_b2.Copy(myband2);
	Go2_Ext_Loop();//next is external (outer) loop 
}
 
//__________________ bands expansion collect valid 5/6

void G17B::ExpandB2() {
	STD_B1_2 &b = myband2;
	bnua = b.nua;
	for (uint32_t i = 0; i < bnua; i++) btua[i] = (uint64_t)b.tua[i] << 32;
	start_active = (uint64_t)BIT_SET_27 << 32;
	nexp = 6;
	mybi2t = bi2_2;
	memset(bi2_2, 0, sizeof bi2_2);
	myokt = vab_2;
	ExpandOneBand(2);
	nbi2_2 = nmybi2t;
	memcpy(b2cpt, p_cpt, sizeof b2cpt);
	b2count = totb2 = b2cpt[0]; // debugging for max expected
	b.my_bi2 = bi2_2;
	b.nbi2 = nbi2_2;
	b.my_validb = vab_2;
	b.nvalidb = (uint32_t)(p_cpt[0] + p_cpt[1]);
}
void G17B::ExpandB1() {
	STD_B1_2 &b = myband1;
	bnua = b.nua;
	for (uint32_t i = 0; i < bnua; i++) btua[i] = b.tua[i];// just expand to 64
	start_active = BIT_SET_27;
	nexp = 6;
	mybi2t = bi2_1;
	memset(bi2_1, 0, sizeof bi2_1);
	myokt = vab_1;
	ExpandOneBand(1);
	nbi2_1 = nmybi2t;
	memcpy(b1cpt, p_cpt, sizeof b1cpt);
	b.my_bi2 = bi2_1;
	b.nbi2 = nbi2_1;
	b.my_validb = vab_1;
	b.nvalidb = (uint32_t)(p_cpt[0] + p_cpt[1]);
}
struct SPOT_E64 {// spots to find band12 valid solutions n clues
	SPOT_E64 * sp;
	uint64_t  all_previous_cells, active_cells;
	uint32_t * start_possibles, n_possibles, ipos, ispot;
	uint64_t * tua;
	uint32_t stack[3], bands[2], missing_clues, nua;
	inline void Copy(SPOT_E64 * old) {
		*this = *old;
		start_possibles += n_possibles;
		ispot++;
		missing_clues--;
		ipos = 0;
		tua += nua;
	}
	inline void AddCellBandStack(int cell, uint32_t * ncb) {
		// if the stack is limit update sn active
		int st = C_stack[cell];
		stack[st]++;
		if (stack[st] > 5) {
			//cout << "stack pleine" << st << endl;
			active_cells &= ~band3xBM[st + 3].u64[0];

		}
		// if the band is limit update sn active
		int b = C_div27[cell];
		bands[b]++;
		if (bands[b] > 5) {// more in mode b 656 only
			//cout << "bande pleine" << b<< endl;
			active_cells &= ~band3xBM[b].u64[0];
		}
	}
	inline void GetOne(uint64_t v) {
		n_possibles = 1;
		bitscanforward64(start_possibles[0], v);
		start_possibles[0] = From_128_To_81[start_possibles[0]];
	}
	inline void AddPossibles(uint64_t v) {
		uint32_t cc;
		while (bitscanforward64(cc, v)) {// look for  possible cells
			register uint64_t bit2 = (uint64_t)1 << cc;
			v ^= (uint64_t)1 << cc;// clear bit
			start_possibles[n_possibles++] = From_128_To_81[cc];
		}
	}
	inline int GetUa(uint64_t v) {
		n_possibles = 0;
		AddPossibles(v);
		return n_possibles;
	}
	uint32_t GetPossibles() {
		// UAs are limited to active cells no empty or single ua
		if (missing_clues < 2) return 0; // minimum to call this process
		uint32_t cells_count[64],
			min_elims = (nua + missing_clues - 1) / missing_clues;
		memset(cells_count, 0, sizeof cells_count);
		uint32_t cc;
		for (uint32_t iua = 0; iua < nua; iua++) {
			register uint64_t Rw = tua[iua] & BIT_SET_2X;
			while (bitscanforward64(cc, Rw)) {// look for  possible cells
				register uint64_t bit2 = (uint64_t)1 << cc;
				Rw ^= (uint64_t)1 << cc;// clear bit
				cells_count[cc]++;
				//if (!cc) cout << "cell 0 pour i=" << iua << " cpt="<< cells_count[cc] << endl;
			}
		}
		//cout << "compte brut critical "<< min_elims << endl;
		//for (int i = 0; i < 64; i++)if (cells_count[i])
		//	cout << From_128_To_81[i] << "\t" << cells_count[i] << endl;

		// collect cells over critical count
		GINT64 tcells[64], temp;
		uint32_t ntcells = 0;
		for (int i = 0; i < 54; i++) {
			register uint32_t my_cell_count = cells_count[C_To128[i]];
			if (my_cell_count >= min_elims) {
				GINT64 & w = tcells[ntcells++];
				w.u32[0] = i;
				w.u32[1] = my_cell_count;
			}
		}
		if (!ntcells) return 0;
		if (ntcells > 1) {// sort in decreasing order
			for (uint32_t i = 0; i < ntcells - 1; i++) {
				for (uint32_t j = i + 1; j < ntcells; j++) {
					if (tcells[i].u64 < tcells[j].u64) {
						temp.u64 = tcells[i].u64;
						tcells[i].u64 = tcells[j].u64;
						tcells[j].u64 = temp.u64;
					}
				}
			}
			//if (ntcells > 64 - missing_clues) ntcells = 64 - missing_clues;
		}
		// load the final table of cells to consider
		//cout << "final count" << endl;
		for (uint32_t i = 0; i < ntcells; i++) {
			start_possibles[i] = tcells[i].u32[0];
			//cout <<i<<"\t"<< tcells[i].u32[0]<<"\t"<< tcells[i].u32[1] <<endl;
		}
		n_possibles = ntcells;
		return ntcells;
	}
	inline int GetLast() {
		n_possibles = 0;
		uint64_t andx = (uint64_t)BIT_SET_2X;
		for (uint32_t iua = 0; iua < nua; iua++) {
			andx &= tua[iua];
			if (!andx)return 0;
		}
		uint32_t cc;
		while (bitscanforward64(cc, andx)) {// look for  possible cells
			register uint64_t bit2 = (uint64_t)1 << cc;
			andx ^= (uint64_t)1 << cc;// clear bit
			start_possibles[n_possibles++] = From_128_To_81[cc];
		}
		return n_possibles;
	}
	inline uint32_t GetPossiblesP() {// using parallel
		uint32_t lim = (nua + missing_clues - 1) / missing_clues,
			lim2 = lim + 5;
		//uint64_t tval[17];
		uint64_t tval[14];
		memset(tval, 0, sizeof tval);
		tval[1] = tua[0] & BIT_SET_2X;
		uint32_t tend = 1; // last used so far
		for (uint32_t iua = 1; iua < nua; iua++) {
			register uint64_t Rw = tua[iua] & BIT_SET_2X,
				Rmore = tval[tend] & Rw;
			if (tend < lim2 && Rmore) tend++;
			switch (tend) {
				//case 15:tval[15]|=tval[14] & Rw;
				//case 14:tval[14]|=tval[13] & Rw;
			case 13:tval[13] |= tval[12] & Rw;
			case 12:tval[12] |= tval[11] & Rw;
			case 11:tval[11] |= tval[10] & Rw;
			case 10:tval[10] |= tval[9] & Rw;
			case 9:tval[9] |= tval[8] & Rw;
			case 8:tval[8] |= tval[7] & Rw;
			case 7:tval[7] |= tval[6] & Rw;
			case 6:tval[6] |= tval[5] & Rw;
			case 5:tval[5] |= tval[4] & Rw;
			case 4:tval[4] |= tval[3] & Rw;
			case 3:tval[3] |= tval[2] & Rw;
			case 2:tval[2] |= tval[1] & Rw;
			}
			tval[1] |= Rw;
		}
		if (tend < lim)return 0;// nothing to do
		n_possibles = 0;
		for (uint32_t i = lim; i < tend; i++) tval[i] &= ~tval[i + 1];
		for (uint32_t i = tend; i >= lim; i--)AddPossibles(tval[i]);
		return n_possibles;
	}
	inline uint32_t GetPossibles2() {// using parallel
		register uint64_t R5 = 0, R4 = 0, R3 = 0, R2 = 0, R1 = 0, R;
		for (uint32_t iua = 1; iua < nua; iua++) {
			R = tua[iua] & BIT_SET_2X;
			R5 |= R4 & R; R4 |= R3 & R; R3 |= R2 & R;	R2 |= R1 & R;
			R1 |= R;
		}
		//if (!R2) return 0;
		R2 &= ~R3; R3 &= ~R4; R4 &= ~R5;
		n_possibles = 0;
		if (R5)					AddPossibles(R5);
		if (R4)					AddPossibles(R4);
		if (R3)					AddPossibles(R3);
		if (R2)					AddPossibles(R2);
		return n_possibles;
	}
	inline uint32_t GetPossibles3() {// using parallel
		register uint64_t R6 = 0, R5 = 0, R4 = 0, R3 = 0, R2 = 0, R1 = 0, R;
		for (uint32_t iua = 1; iua < nua; iua++) {
			R = tua[iua] & BIT_SET_2X;
			R6 |= R5 & R; R5 |= R4 & R; R4 |= R3 & R; R3 |= R2 & R;	R2 |= R1 & R;
			R1 |= R;
		}
		//if (!R3) return 0;
		R3 &= ~R4; R4 &= ~R5; R5 &= ~R6;
		n_possibles = 0;
		if (R6)	AddPossibles(R6);
		if (R5)	AddPossibles(R5);
		if (R4)	AddPossibles(R4);
		if (R3)	AddPossibles(R3);
		return n_possibles;
	}
	inline uint32_t GetPossibles4() {// using parallel
		register uint64_t R6 = 0, R5 = 0, R4 = 0, R3 = 0, R2 = 0, R1 = 0, R;
		for (uint32_t iua = 1; iua < nua; iua++) {
			R = tua[iua] & BIT_SET_2X;
			R6 |= R5 & R; R5 |= R4 & R; R4 |= R3 & R; R3 |= R2 & R;	R2 |= R1 & R;
			R1 |= R;
		}
		//if (!R4) return 0;
		R4 &= ~R5; R5 &= ~R6;
		n_possibles = 0;
		if (R6)	AddPossibles(R6);
		if (R5)	AddPossibles(R5);
		if (R4)	AddPossibles(R4);
		return n_possibles;
	}
	inline uint32_t GetPossibles5() {// using parallel
		register uint64_t R7 = 0, R6 = 0, R5 = 0, R4 = 0, R3 = 0, R2 = 0, R1 = 0, R;
		for (uint32_t iua = 1; iua < nua; iua++) {
			R = tua[iua] & BIT_SET_2X;
			R7 |= R6 & R; R6 |= R5 & R; R5 |= R4 & R; R4 |= R3 & R; R3 |= R2 & R;	R2 |= R1 & R;
			R1 |= R;
		}
		//if (!R5) return 0;
		R4 &= ~R5; R5 &= ~R6; R6 &= ~R7;
		n_possibles = 0;
		if (R7)	AddPossibles(R7);
		if (R6)	AddPossibles(R6);
		if (R5)	AddPossibles(R5);
		return n_possibles;
	}
	void D1() {
		cout << "D1\t" << ispot << "\t" << ipos << endl;
		//cout << Char2Xout(active_cells) << " spot=" << ispot << " pos=" << ipos
		//	<< " npos=" << n_possibles << endl;
	}
	void D2(int all = 0) {
		//cout << Char2Xout(active_cells) << " Shrink active"  << endl;
		cout << Char2Xout(all_previous_cells) << " known nuas=" << nua << endl;
		if (!all) return;
		for (uint32_t iua = 0; iua < nua; iua++)
			cout << Char2Xout(tua[iua]) << endl;
	}

	void D3() {
		cout << Char2Xout(all_previous_cells) << " known" << endl;
		cout << "get possibles  nua=" << nua << "  missing_clues " << missing_clues
			<< "n_possibles" << n_possibles << endl;
	}
	inline void Ddead(uint32_t iua) {
		//			cout <<"\t\t"<< ispot <<" "<<ipos<<" dead branch iua=" << iua << endl;
	}
	inline void Dass(uint32_t iua, uint64_t Ru) {
		//			cout<<"\t\t" << ispot << " " << ipos << " assign iua=" << iua
		//				<< " cell=" << start_possibles[0] << endl;
		//			cout << ispot << " " << ipos << " assign iua=" << iua <<endl
		//				<< Char2Xout(Ru) <<" cell="<< start_possibles[0] << endl;
	}
	inline void DNoMore() {
		cout << Char2Xout(all_previous_cells) << "no more uas "
			<< bands[0] << bands[1] << " "
			<< stack[0] << stack[1] << stack[2] << endl;

	}
};
void G17B::ExpandOneBand(int ib) {
	int notfirst = 0,
		c56 = 011;// collect 5_6
	if (sgo.vx[4]) {// this is pass b only 656
		if (ib == 1)c56 = 010; // must be 6
		else c56 = 1;// must be 5
	}
	memset(p_cpt, 0, sizeof p_cpt);
	uint32_t tclues[12], bufp[1000], lastspot = nexp - 1;
	uint64_t bufua12[30000];
	SPOT_E64 spt[8], *s, *sn;
	nmybi2t = nmyokt = 0;
	s = spt;
	memset(s, 0, sizeof spt[0]);// init the stack status ...
	s->missing_clues = nexp;
	s->active_cells = start_active;// all cells active
	s->start_possibles = bufp;
	s->tua = bufua12;// copy init table to the buffer
	s->nua = bnua;
	memcpy(s->tua, btua, 8 * bnua);
	s->GetUa(s->tua[0] & BIT_SET_2X);
next:
	{// catch and apply cell in bitfields
		uint32_t iw = s->ipos++;
		if (iw >= s->n_possibles)goto back;
		register uint32_t cell = s->start_possibles[iw];
		tclues[s->ispot] = cell;
		register uint64_t bit = (uint64_t)1 << C_To128[cell];
		register  uint64_t filter = s->all_previous_cells | bit,
			ac = s->active_cells ^ bit;
		sn = s + 1; sn->Copy(s); // prepare next spot
		sn->all_previous_cells = filter;
		sn->active_cells = s->active_cells = ac;
		{//level>0 shrink the ua table in new
			register uint32_t nua1 = s->nua, iua;
			sn->nua = 0;
			register uint64_t Ra = sn->active_cells;
			for (iua = 0; iua < nua1; iua++) {
				register uint64_t Ru = s->tua[iua];
				if (Ru&filter)continue;
				Ru &= Ra;
				if (!Ru) goto next;// dead branch  
				register uint64_t cc = _popcnt64(Ru);
				Ru |= (cc << 59);
				AddUA64(sn->tua, sn->nua, Ru);
			}
		}
		if (s->ispot == 1) {// open a new index2
			if (notfirst) {// save previous if active
				BI2 & pr = mybi2t[nmybi2t];
				if (pr.istart != pr.iend) {
					nmybi2t++;
					BI2 & pn = mybi2t[nmybi2t];
					pn.istart = pn.iend = pr.iend;
				}
			}
			notfirst = 1;
			BI2 & pn = mybi2t[nmybi2t];// init the ne status
			pn.bf = sn->all_previous_cells;
			pn.active = sn->active_cells;
			memcpy(pn.tval, tclues, sizeof pn.tval);
		}
		if (!sn->nua)goto no_more_uas;
		else if (sn->missing_clues == 1) { if (!sn->GetLast())goto next; }
		else sn->GetUa(sn->tua[0] & BIT_SET_2X);
		if (s->ispot < lastspot)s++;// {	s++; s->D3(); }

		goto next;
	}
no_more_uas:
	{	BI2 & pi = mybi2t[nmybi2t];
	register uint64_t R0 = sn->all_previous_cells;// ^pi.bf;
	if (s->ispot == 5) {
		if (c56 & 010) {
			myokt[pi.iend++].Enter(R0, &tclues[2]);
			p_cpt[1]++;
		}
	}
	else {//if below 6 loop  for redundant clues
		int tc[64], nt = 0;
		uint64_t tbit[64];
		{	uint32_t register xcell;
		register uint64_t ac = s->active_cells&BIT_SET_2X;
		while (bitscanforward64(xcell, ac)) {// put active cells in table
			uint64_t bit = (uint64_t)1 << xcell;
			ac ^= bit;
			tc[nt] = From_128_To_81[xcell];
			tbit[nt++] = bit;
		}
		if (s->ispot == 4) { //5 clues
			if (c56 & 1) {
				myokt[pi.iend++].Enter(R0, &tclues[2]);
				p_cpt[0]++;
			}
			for (int i6 = 0; i6 < nt; i6++) {
				tclues[5] = tc[i6];
				if (c56 & 010) {
					myokt[pi.iend++].Enter(tbit[i6] | R0, &tclues[2]);
					p_cpt[1]++;
				}
			}
		}
		else if (s->ispot == 3) { // valid 4 clues
			for (int i5 = 0; i5 < nt; i5++) {
				register uint64_t R5 = tbit[i5] | R0;
				tclues[4] = tc[i5];
				if (c56 & 1) {
					myokt[pi.iend++].Enter(R5, &tclues[2]);
					p_cpt[0]++;
				}
				if (c56 & 010)for (int i6 = i5 + 1; i6 < nt; i6++) {
					tclues[5] = tc[i6];
					myokt[pi.iend++].Enter(tbit[i6] | R5, &tclues[2]);
					p_cpt[1]++;
				}
			}
		}
		else if (s->ispot == 2) { // valid 3 clues
			for (int i4 = 0; i4 < nt; i4++) {
				register uint64_t R4 = tbit[i4] | R0;
				tclues[3] = tc[i4];
				for (int i5 = i4 + 1; i5 < nt; i5++) {
					register uint64_t R5 = tbit[i5] | R4;
					tclues[4] = tc[i5];
					if (c56 & 1) {
						myokt[pi.iend++].Enter(R5, &tclues[2]);
						p_cpt[0]++;
					}
					if (c56 & 010) for (int i6 = i5 + 1; i6 < nt; i6++) {
						tclues[5] = tc[i6];
						myokt[pi.iend++].Enter(tbit[i6] | R5, &tclues[2]);
						p_cpt[1]++;
					}
				}
			}
		}
		else if (s->ispot == 1) { // valid 2 clues
			for (int i3 = 0; i3 < nt; i3++) {
				register uint64_t R3 = tbit[i3] | R0;
				tclues[2] = tc[i3];
				for (int i4 = i3 + 1; i4 < nt; i4++) {
					register uint64_t R4 = tbit[i4] | R3;
					tclues[3] = tc[i4];
					for (int i5 = i4 + 1; i5 < nt; i5++) {
						register uint64_t R5 = tbit[i5] | R4;
						tclues[4] = tc[i5];
						if (c56 & 1) {
							myokt[pi.iend++].Enter(R5, &tclues[2]);
							p_cpt[0]++;
						}
						if (c56 & 010)	for (int i6 = i5 + 1; i6 < nt; i6++) {
							tclues[5] = tc[i6];
							myokt[pi.iend++].Enter(tbit[i6] | R5, &tclues[2]);
							p_cpt[1]++;
						}
					}
				}

			}
		}
		}
	}
	}
	goto next;

back:
	if (--s >= spt)goto next;
	// save the last index if 
	BI2 & pr = mybi2t[nmybi2t];
	if (pr.istart != pr.iend) 	nmybi2t++;
}

//____ ___________________________________initial tasks collect uas 12;  guas
void G17B::GoM10Uas() {
	//=========================== collect UAs  old process 
	zh1b_g.modegua = 0;//must be to activate filter in UAs b12 more
	if (genuasb12.Initgen()) return;
	genb12.BuildGang9x3();
	// _____ GUAs 
	tguas.InitStart();
	zh1b_g.modegua = 1;//must be to kill  filter in GUAs 6_7 more
	genb12.SecondSockets2Setup();// collect GUA2s 
	genb12.SecondSockets3Setup();// collect GUA3s 

	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3 & b = genb12.bands3[ib3];
		b.guas.isguasocketc2.Convert3X27to81(b.guas.isguasocket2);
		b.guas.isguasocketc3.Convert3X27to81(b.guas.isguasocket3);
		b.guas.isguasocketc2_46.Convert3X27to81(b.guas.isguasocket2_46);
		b.isguasocketc246 = b.guas.isguasocketc2 | b.guas.isguasocketc2_46;
		b.Setup_bf162b();
	}

}
void G17B::StackUas() {// band 3 uas used as gangsters via  C_transpose_d[81]
	STD_B3 wbs;

	// transpose bands 1+2  transpose bands 3
	for (int ib = 0; ib < genb12.nband3; ib++) {// bands 3
		STD_B3 &myb = genb12.bands3[ib];
		wbs.ntua128 = 0;
		memcpy(&genb12.grid0[54], myb.band0, 4 * 27);
		int zt[81];
		for (int i = 0; i < 81; i++) {
			zt[i] = genb12.grid0[C_transpose_d[i]];
			//cout << genb12.grid0[i] +1;
		}
		//cout << endl;
		BANDMINLEX::PERM perm_ret;
		//___ stack 1
		bandminlex.Getmin(zt, &perm_ret);
		int ib1 = perm_ret.i416;
		wbs.InitStack(ib1, zt, perm_ret,0);

		//___ stack 2
		bandminlex.Getmin(&zt[27], &perm_ret);
		int ib2 = perm_ret.i416;
		wbs.InitStack(ib2, &zt[27], perm_ret,1);	
		//___ stack 3
		bandminlex.Getmin(&zt[54], &perm_ret);
		int ib3 = perm_ret.i416;
		wbs.InitStack(ib3, &zt[54], perm_ret,2);

		// check the final status and copy to band3
		myb.ntua128 = 0;
		for (uint32_t i = 0; i < wbs.ntua128; i++) {
			BF128 w = wbs.tua128[i];
			int cc3 = _popcnt32(w.bf.u32[2]);
			if(cc3>3)	myb.tua128[myb.ntua128++] = wbs.tua128[i];
			else if (cc3 ==2) {// be sure to have gua 2 clues in table
				uint64_t cc12 = _popcnt64(w.bf.u64[0]);
				if (cc12 < 12) continue; // must be there

				// find the digits pattern from the current band 3
				int * cur_b3 = myb.band0, wdigs = 0, c27;
				{
					register uint32_t wua = w.bf.u32[2];
					while (bitscanforward(c27, wua)) {// look for  possible cells
						wua ^= 1 << c27;// clear bit
						wdigs |= 1 << cur_b3[c27];
					}
				}
				uint32_t my_i81 = genb12.GET_I81_G2(wdigs, w.bf.u32[2]);

				int ix_start = tguas.ix_start.g2[my_i81];
				if (ix_start >= 0) {
					GUA & g =tguas.tgua_start[ix_start];
					AddUA64(g.tua, g.nua, w.bf.u64[0]);
				}
				else myb.tua128[myb.ntua128++] = wbs.tua128[i];
			}
		}
	}
}
void STD_B3::InitStack(int i16, int * z0, BANDMINLEX::PERM & p, int iband) {
	i416 = i16;
	int dstack = 3 * iband;
	GetUAs();
	// create the cell map in output
	for (int i = 0; i < 3; i++) {
		int vr = 9 * p.rows[i], vr0 = 9 * i;
		for (int j = 0; j < 9; j++)
			map[vr0 + j] = vr + p.cols[j];
	}
	// morph all uas
	for (uint32_t i = 0; i < nua; i++) {
		register int uao = tua[i] & BIT_SET_27;
		BF128  ua; ua.SetAll_0();
		for (int i=0, bit = 1; i < 27; i++, bit <<= 1)
			if (bit & uao) {// set the bit in the band form
				int cell = C_transpose_d[map[i]];// +dstack;
				ua.Set_c(cell+dstack);
			}
		tua128[ntua128++] = ua;
	}
}

//_________________ external loop on bands 1+2 small uas

/*external loop control
 the external loop cut the process using small UAs bands 1+2
 if 'yes' is the part hitting a ua in a band
 yes1 hit in band 1 ; yes2 hit in band2
 the process can then be cut in 2 chunks
 yes1 * all2		all1-yes1 * yes2
 this is of interest if the count is significantly lower than
		all1 * all 2
 and if all1*all2 is big enough

 the process can be split several times with different uas

 here, the first chunk yes * all is split again
 the main loop continues to split    all1-yes1 * yes2
 Note the first band can be band 1 or band 2

 Control of the main loop

 EXLRATIO is the minimal wanted reduction of the count
 EXLNLOOP1 the maximum number of steps in the main loop
 EXLNLOOP2 the maximum number of steps in the first chunk split
 EXLBREAK the count minimum to loop
 */
#define EXLRATIO 80
#define EXLNLOOP1 4
#define EXLNLOOP2 3
#define EXLBREAK 4000000


void BINDEXN::Copy(STD_B1_2 & b) {
	ntvb = b.nvalidb;
	for (uint32_t i = 0; i < ntvb; i++) tvb[i] = b.my_validb[i];
	nt2 = b.nbi2;
	memcpy(t2, b.my_bi2, nt2 * sizeof t2[0]);
}
void BINDEXN::Copy(BINDEXN & b) {
	ntvb = b.ntvb;
	memcpy(tvb, b.tvb, ntvb * sizeof tvb[0]);
	nt2 = b.nt2;
	memcpy(t2, b.t2, (nt2 + 1) * sizeof t2[0]);
}

uint64_t G17B::FindSockets(uint64_t active, uint64_t lim) {
	nextl1 = 0, nextl2 = 0;
	uint64_t *t = genuasb12.tua;
	uint32_t n = genuasb12.nua;
	for (uint32_t i = 0; i < n; i++) {
		register uint64_t U = t[i] & active,
			U1 = U & BIT_SET_27, U2 = U & BIT_SET_B2;
		if ((!U1) || (!U2))continue;
		uint64_t n = _popcnt64(U& BIT_SET_27);
		if ((lim == 2 && n < lim) || n == lim) {
			uint32_t i1, aig = 1;

			for (i1 = 0; i1 < nextl1; i1++)
				if (U1 == extl1[i1].bfx) { aig = 0; break; }
			if (aig) {
				if (nextl1 < 6) { // open a new window (limit 10
					i1 = nextl1++; aig = 0;
					extl1[i1].Init((uint32_t)U1, 1);
				}
			}
			if (!aig)if (extl1[i1].ntbfy < 20)
				extl1[i1].tbfy[extl1[i1].ntbfy++] = (uint32_t)(U2 >> 32);

		}
		n = _popcnt64(U& BIT_SET_B2);
		if ((lim == 2 && n < lim) || n == lim) {
			register uint32_t U2a = (uint32_t)(U2 >> 32);
			//if (nt2sk < 5)t2sk[nt2sk++] = (uint32_t)(U2 >> 32);
			uint32_t i2, aig = 1;
			for (i2 = 0; i2 < nextl2; i2++)
				if (U2a == extl2[i2].bfx) { aig = 0; break; }
			if (aig) {
				if (nextl2 < 6) { // open a new window (limit 10
					i2 = nextl2++; aig = 0;
					extl2[i2].Init(U2a, 2);
				}
			}
			if (!aig)if (extl2[i2].ntbfy < 20)
				extl2[i2].tbfy[extl2[i2].ntbfy++] = (uint32_t)U1;
		}
	}
	if (nextl1 | nextl2) {	return lim - 1;	}
	else	return 0;
}

uint64_t GetYes_b1( VALIDB *vb, uint32_t nvb, uint32_t bf) {
	uint64_t n = 0;
	register uint64_t F = bf;
	for (uint32_t i = 0; i < nvb; i++) 
		if (F& vb[i].bf)n++;	
	return n;
}
uint64_t GetYes_b1(VALIDB *vb, uint32_t nvb, uint32_t *tbf, uint32_t ntbf) {
	uint64_t n = 0;
	for (uint32_t i = 0; i < nvb; i++) {
		VALIDB &w = vb[i];
		uint32_t aig = 0;
		register uint32_t U =(uint32_t) w.bf;
		for (uint32_t j = 0; j < ntbf; j++) {
			if ((tbf[j] & U)) { aig = 1; break; }
		}
		n += aig;
	}
	return n;
}
uint64_t GetYes_b2(VALIDB *vb, uint32_t nvb, uint32_t bf) {
	uint64_t n = 0;
	register uint64_t F = (uint64_t)bf << 32;
	for (uint32_t i = 0; i < nvb; i++) 
		if (F& vb[i].bf)n++;	
	return n;
}
uint64_t GetYes_b2(VALIDB *vb, uint32_t nvb, uint32_t *tbf, uint32_t ntbf) {
	uint64_t n = 0;
	for (uint32_t i = 0; i < nvb; i++) {
		VALIDB &w = vb[i];
		uint32_t aig = 0;
		register uint32_t U = (uint32_t)(w.bf >> 32);
		for (uint32_t j = 0; j < ntbf; j++) {
			if ((tbf[j] & U)) { aig = 1; break; }
		}
		n += aig;
	}
	return n;
}

void G17B::ExtractMin(uint64_t active, BINDEXN & bin1, BINDEXN & bin2) {
	uint64_t limb1 = _popcnt64(active & BIT_SET_27) - 3,
		limb2 = _popcnt64(active & BIT_SET_B2) - 3;
	for (uint32_t i = 0; i < nextl1; i++) {
		uint32_t orw = 0, *t = extl1[i].tbfy;
		for (uint32_t j = 0; j < extl1[i].ntbfy; j++) {
			orw |= t[j];
		}
		if (_popcnt32(orw) > (uint32_t)limb1) continue;//min 3 clues will end as 100% ratio
		n_yesb1= GetYes_b1(bin1.tvb, bin1.ntvb, extl1[i].bfx);
		if (extl1[i].ntbfy == 1)
			n_yesb2 = GetYes_b2(bin2.tvb, bin2.ntvb, extl1[i].tbfy[0]);
		else n_yesb2 = GetYes_b2(bin2.tvb, bin2.ntvb, extl1[i].tbfy, extl1[i].ntbfy);
		n_nob1 = bin1.ntvb - n_yesb1;
		n_nob2 = bin2.ntvb - n_yesb2;

		uint64_t ratio1 = n_yesb1 * bin2.ntvb + n_yesb2 * n_nob1,
			ratio = (uint64_t)100 * ratio1 / bin1.ntvb / bin2.ntvb;

		if (ratio < minratio) {
			minratio = ratio;
			extl1[i].ratio = (uint32_t)ratio;
			extlw = extl1[i];
			extlw.noxyes = n_yesb2 * n_nob1;
		}

	}
	for (uint32_t i = 0; i < nextl2; i++) {
		uint32_t orw = 0, *t = extl2[i].tbfy;
		for (uint32_t j = 0; j < extl2[i].ntbfy; j++) {
			orw |= t[j];
		}
		if (_popcnt32(orw) > (uint32_t)limb2) continue;// will end as 100% ratio
		n_yesb2 = GetYes_b2(bin2.tvb, bin2.ntvb, extl2[i].bfx);

		if (extl2[i].ntbfy == 1)
			n_yesb1 = GetYes_b1(bin1.tvb, bin1.ntvb, extl2[i].tbfy[0]);
		else n_yesb1 = GetYes_b1(bin1.tvb, bin1.ntvb, extl2[i].tbfy, extl2[i].ntbfy);
		n_nob1 = bin1.ntvb - n_yesb1;
		n_nob2 = bin2.ntvb - n_yesb2;
		uint64_t ratio1 = n_yesb2 * bin1.ntvb + n_yesb1 * n_nob2,
			ratio = (uint64_t)100 * ratio1 / bin1.ntvb / bin2.ntvb;
		if (ratio < minratio) {
			minratio = ratio;
			extl2[i].ratio = (uint32_t)ratio;
			extlw = extl2[i];
			extlw.noxyes = n_yesb1 * n_nob2;
		}
	}
}


void G17B::ExtSplitY(BINDEXN & binw, uint32_t *tbf, uint32_t ntbf,
	uint32_t & activer, int bande) {// source binw exit yes in binw
	uint32_t lim_source = binw.nt2;// table no is also source
	binw.nt2 = binw.ntvb = 0;
	VALIDB *vb = binw.tvb;
	activer = 0; // or status in new vb1 (no) start null
	for (uint32_t i = 0; i < lim_source; i++) { // all index 2
		BI2 wi = binw.t2[i], win = wi;
		win.istart = binw.ntvb;
		for (uint32_t j = wi.istart; j < wi.iend; j++) {
			VALIDB wj = binw.tvb[j];
			register uint64_t Bf = wj.bf;

			for (uint32_t k = 0; k < ntbf; k++) {
				register uint64_t F = tbf[k]; // UA to hit to be yes
				if (bande == 2)F <<= 32;
				if (Bf & F) {// a new "yes" 
					vb[binw.ntvb++] = wj;
					activer |= wj.bf;
					break;// one hit is enough 
				}
			}
		}
		if (binw.ntvb != win.istart) {// new group in "yes"
			win.iend = binw.ntvb;
			binw.t2[binw.nt2++] = win;
		}
	}
}
void  G17B::ExtSplitX(BINDEXN & bin1no, BINDEXN & bin1yes,
	uint32_t bf, uint32_t & activer, int bande) {
	uint32_t lim_source = bin1no.nt2;// table no is also source
	// initial status empty for table yes and no
	bin1no.nt2 = bin1no.ntvb = 0;
	bin1yes.nt2 = bin1yes.ntvb = 0;
	activer = 0; // or status in new vb1 (no) start null
	VALIDB *vb1 = bin1no.tvb, *vb2 = bin1yes.tvb;

	register uint64_t F = bf; // UA to hit to be yes
	if (bande == 2)F <<= 32;
	for (uint32_t i = 0; i < lim_source; i++) { // all index 2
		BI2 wi = bin1no.t2[i], win = wi, win1;
		win.istart = bin1yes.ntvb;
		if (F&wi.bf) {// all the group is "yes"
			uint32_t n = wi.iend - wi.istart;
			win.iend = bin1yes.ntvb + n;
			memcpy(&vb2[bin1yes.ntvb], &vb1[wi.istart], n * sizeof  vb2[0]);
			bin1yes.t2[bin1yes.nt2++] = win;
			bin1yes.ntvb += n;
			continue;
		}
		// must check each validb of the group
		win1 = wi;
		win1.istart = bin1no.ntvb;
		for (uint32_t j = wi.istart; j < wi.iend; j++) {
			VALIDB wj = vb1[j];
			if (F&wj.bf) 	vb2[bin1yes.ntvb++] = wj;
			else { vb1[bin1no.ntvb++] = wj; activer |= wj.bf; }
		}
		if (bin1no.ntvb != win1.istart) {// new group in "no"
			win1.iend = bin1no.ntvb;
			bin1no.t2[bin1no.nt2++] = win1;
		}
		if (bin1yes.ntvb != win.istart) {// new group in "yes"
			win.iend = bin1yes.ntvb;
			bin1yes.t2[bin1yes.nt2++] = win;
		}
	}

}

int CheckBf(BINDEXN & binw, uint64_t bfw) {
	register VALIDB * tvb = binw.tvb;
	for (uint32_t i = 0; i < binw.ntvb; i++) {
		if (tvb[i].bf == bfw)return i;
	}
	return -1;
}

void G17B::Go2_Ext_Loop() {	//_____________ outer loop
	loopb1 = 0;
	uint32_t activerb1, activerb2;
	uint64_t activeloop = BIT_SET_2X;
	activerb1= activerb2 = BIT_SET_27;
	//_________________ external loop
	while (++loopb1 << EXLNLOOP1) {
		if (aigstop)return;
		minratio = extlr.ratio=1000;
		uint64_t ir = FindSockets(activeloop,2);
		if (ir) 	ExtractMin(activeloop, bin_b1, bin_b2);
		if (!ir || minratio > EXLRATIO) {
			ir = FindSockets(activeloop,3);
			if (ir)ExtractMin(activeloop, bin_b1, bin_b2);
			if (!ir || minratio > EXLRATIO) {
				ir = FindSockets(activeloop,4);
				if (ir)ExtractMin(activeloop, bin_b1, bin_b2);
			}
		}
		if (minratio > EXLRATIO) break;
		else {
			extlr = extlw;
			if (extlw.mode == 1) {// this is a band1 X band2 Y
				ExtSplitX(bin_b1, bin_b1yes, extlw.bfx, activerb1);
				if (loopb1 == 1)
					Go2b_Ext_Loop(BIT_SET_2X | activerb1, 1);
				else  Go3(bin_b1yes, bin_b2);
				ExtSplitY(bin_b2, extlr.tbfy, extlr.ntbfy, activerb2, 2);
			}
			else {// this is a band2 X band1 Y
				ExtSplitX(bin_b2, bin_b2yes, extlw.bfx, activerb2, 2);
				if (loopb1 == 1)
					Go2b_Ext_Loop(BIT_SET_2X, 2);
				else  Go3(bin_b1, bin_b2yes);
				ExtSplitY(bin_b1,	extlr.tbfy, extlr.ntbfy,  activerb1);
			}
		}
		activeloop = activerb2; activeloop <<= 32; activeloop |= activerb1;
		if (extlr.noxyes < EXLBREAK) break;
	}
	Go3(bin_b1, bin_b2);// last call
}

void G17B::Go2b_Ext_Loop(uint64_t activeloop, uint32_t mode2) {
		//___________init the working areas bin2_b1, bin2_b2
		if (mode2 == 1) {// b1 yes b2 all
			bin2_b1.Copy(bin_b1yes);
			bin2_b2.Copy(bin_b2);
		}
		else {// b1 all b2 yes
			bin2_b1.Copy(bin_b1);
			bin2_b2.Copy(bin_b2yes);
		}
		uint32_t loopb2 = 0;
		uint32_t activerb1, activerb2;
		while (++loopb2 <= EXLNLOOP2) {
			if (aigstop) return;
			minratio = 1000;
			uint64_t ir = FindSockets(activeloop, 2);
			if (ir)  ExtractMin(activeloop, bin2_b1, bin2_b2);
			if (!ir || minratio > EXLRATIO) {
				ir = FindSockets(activeloop, 3);
				if (ir)ExtractMin(activeloop, bin2_b1, bin2_b2);
				if (!ir || minratio > EXLRATIO) {
					ir = FindSockets(activeloop, 4);
					if (ir)ExtractMin(activeloop, bin2_b1, bin2_b2);
				}
			}
			if (minratio >= EXLRATIO) break;
			else {
				if (extlw.mode == 1) {// this is a band1 X band2 Y
					ExtSplitX(bin2_b1, bin2_b1yes, extlw.bfx, activerb1);
					Go3(bin2_b1yes, bin2_b2);
					ExtSplitY(bin2_b2, extlw.tbfy, extlw.ntbfy, activerb2,2);
				}
				else {// this is a band2 X band1 Y
					ExtSplitX(bin2_b2, bin2_b2yes, extlw.bfx, activerb2,2);
					Go3(bin2_b1, bin2_b2yes);
					ExtSplitY(bin2_b1, extlw.tbfy, extlw.ntbfy, activerb1);
				}
			}
			activeloop = activerb2; activeloop <<= 32; activeloop |= activerb1;
			if (extlw.noxyes < EXLBREAK) break;
		}
		Go3(bin2_b1, bin2_b2);// last call
	}


//_________"step" _________ process a subset of valid {band1;band2}

	/* the sub lot is made of 2 pieces of the expansion
		usually the smaller piece is in band 1

		here 2 loops, outer band1 inner band 2
		for each loop, the process in cut in sub lots
		  depending on the expansion index (2 common cells)
		  valid bands are split by size
		  a pair of 2 cells index gives a step

		 UAs and GUAs tables are reduced to reach a "step size"
		 common cells and possible cells of the step are identified to optimize the process
	*/
	
	int G17B::G3_SplitBi2( BINDEXN & binw, uint32_t ibi2,
		INDEX_XY & ixyw, VALIDB64 * pvb) {// build tables per size 

			//___ sort  the lot per size of valid bands
		uint32_t buffer[MAX_56],
			nv[2] = { 0,0 }, // count per size
			stv[2] = { 0,MAXN5 },// start in buffer per size
			*psv[2]; // pointers to starts
		for (int i = 0; i < 2; i++) psv[i] = &buffer[stv[i]];

		BI2 biw = binw.t2[ibi2];
		ixyw.bf = biw.bf ;
		VALIDB * tvb = binw.tvb;
		uint32_t id = biw.istart, iend = biw.iend;
		for (uint32_t iv = id; iv < iend; iv++) {
			VALIDB & wv = tvb[iv];
			uint32_t nc =(uint32_t) wv.nval;
			nc -= 3;// 3=0 for 5 4=>1 for 6
			psv[nc][nv[nc]++] = iv;
		}
		ixyw.ntotvb = 0;
		for (int isize = 1; isize >= 0; isize--) {
			if (nv[isize]) {
				ixyw.ntotvb += nv[isize];
				ixyw.ncluesmin = isize;
			}
		}
		if (!ixyw.ntotvb)return 1;//  should not be !!!

		//_________ build the final tables in 64 mode in ixyw
		uint64_t	wand = BIT_SET_2X, wor = 0;
		VALIDB64 * pvb64 = pvb;
		uint32_t sumpvb = 0;
		for (int isize = 0; isize < 2; isize++) {//5/6
			uint32_t *psvw = psv[isize];
			INDEX_XY::ITEM &itemw = ixyw.titem[isize];
			itemw.ntvb = nv[isize];
			sumpvb += itemw.ntvb;
			itemw.sum_vb = sumpvb;
			itemw.tvb = pvb64;
			for (uint32_t j = 0; j < nv[isize]; j++) {
				VALIDB vb = tvb[psvw[j]];// source valid band
				VALIDB64 & vb64 = *pvb64++;
				uint64_t bf = vb.bf; 
				wand &= bf;
				wor |= bf;
				vb64.bf = bf;
				vb64.nval = vb.nval + 2;
				memcpy(vb64.tval, biw.tval, sizeof biw.tval);
				memcpy(&vb64.tval[2], vb.tval, sizeof vb.tval);
			}
		}
		ixyw.and_g = wand;
		ixyw.or_g = wor;
		return 0;
	}


	void G17B::Go3(BINDEXN & bin1, BINDEXN & bin2) {
		p_cpt2g[2]++;
		if (aigstop) return;
		if ((!bin1.ntvb) || (!bin2.ntvb)) return;

		//__________ loop on B1
		for (uint32_t ib1 = 0; ib1 < bin1.nt2; ib1++) {
			if (aigstop) return;
			if (G3_SplitBi2( bin1, ib1, index_xy_b1, vab64b1))
				continue;// X should never be killed
			fb1 = index_xy_b1.and_g;// including more common cells
			acb1 = index_xy_b1.or_g;// all these cells can be active in the step
			Apply_Band1_Step();// reduce uas build cells vectors

			// reset more uas tables 				   
			moreuas_12_13.Init();			moreuas_14.Init();
			moreuas_15.Init();			moreuas_AB_small.Init();
			moreuas_AB.Init();			moreuas_AB_big.Init();
			p_cpt2g[3]++;
					   
		//___________ loop on B2
			for (uint32_t ib2 = 0; ib2 < bin2.nt2; ib2++) {
				if (aigstop) return;
				p_cpt2g[4]++;
				if (G3_SplitBi2( bin2, ib2, index_xy_b2, vab64b2))		continue;
				if (Apply_Band2_Step())continue;// check dead branch set vectors 64
				tguas.ApplyB2();

				// and process the step filtering the first 4 uas
				n_to_clean = 0;// count "to clean" to 0
				Do64uas();
				CleanAll();
			}
		}
	}

	//___________ small functions to handle the uas band 1+2

	void BuildBaseAndCellsVector(uint32_t nuas, BF128 & bv, BF128 * cellsv, uint64_t * tu) {
		bv = maskLSB[nuas];// Uas vector
		memset(cellsv, 255, sizeof g17b.vc64_192);// all bits to 1
		uint32_t cc64;// build cells vectors A
		for (uint32_t i = 0; i < nuas; i++) {
			register uint64_t Rw = tu[i] & BIT_SET_2X;
			while (bitscanforward64(cc64, Rw)) {// look for  possible cells
				Rw ^= (uint64_t)1 << cc64;// clear bit
				cellsv[From_128_To_81[cc64]].clearBit(i);
			}
		}
	}
	void AddUaToVector(uint64_t ua12, BF128 * cellsv, uint32_t iloc) {
		register uint64_t Rw = ua12 & BIT_SET_2X;
		uint32_t cc64;// build cells vectors A
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			cellsv[From_128_To_81[cc64]].clearBit(iloc);
		}
	}
	inline void SetUpBaseVector(uint32_t ncl, uint32_t *tcl,
		BF128 & bv, BF128 * cellsv, BF128 & bvf) {
		bvf = bv;// start with valid uas
		for (uint32_t i = 0; i < ncl; i++) {//all cells common to the step
			bvf &= cellsv[tcl[i]];
		}
	}
	inline int SetUpStepV(uint32_t * tc, uint32_t ntc, BF128 & vb, BF128 & vd, BF128 * tvc) {
		vd = vb;
		if (vb.isEmpty())return 1;
		for (uint32_t i = 0; i < ntc; i++)vd &= tvc[tc[i]];
		return vd.isEmpty();
	}


void G17B::Apply_Band1_Step() {// shrink the uas table	
	{//shrink the ua table  no dead here ???
		register uint64_t filter=fb1,
			Ra = acb1 | BIT_SET_B2;
		register uint64_t * tua = genuasb12.tua;
		register uint32_t nua = genuasb12.nua;
		ntusb1 =  0;
		for (uint32_t iua = 0; iua < nua; iua++) {
			register uint64_t Ru = tua[iua];
			if (Ru&filter) continue;
			Ru &= Ra;
			Ru |= ((uint64_t )_popcnt64(Ru) << 59);
			if (ntusb1 >= 960)ntusb1 = 959;// limit for vectors
			AddUA64(tusb1, ntusb1, Ru);
		}
	}
	uint32_t ntua_64 = ntusb1;
	if (ntua_64 > 64) ntua_64 = 64;
	v64uas = maskLSB[ntua_64].u64[0];// Uas vector
	memset(vc64, 255, sizeof vc64);// all bits to 1
	uint32_t cc64;// build cells vectors A
	uint64_t biti = 1;
	for (uint32_t i = 0; i < ntua_64; i++, biti <<= 1) {
		register uint64_t Rw = tusb1[i] & BIT_SET_2X;
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			vc64[From_128_To_81[cc64]] ^= biti;
		}
	}
	// build other cells  vectors as needed 192 320 448 576 704 832 960

	if (ntusb1 > 64) {// 64_192
		uint32_t ntua_192 = (ntusb1 > 192) ? 128 : ntusb1 - 64;
		BuildBaseAndCellsVector(ntua_192, v64_192uas, vc64_192, &tusb1[64]);
	}
	else {
		memset(&v64_192uas, 0, sizeof v64_192uas);
		memset(vc64_192, 255, sizeof vc64_192);
	}
	if (ntusb1 > 192) {// 192_320
		uint32_t ntua_320 = (ntusb1 > 320) ? 128 : ntusb1 - 192;
		BuildBaseAndCellsVector(ntua_320, v192_320uas, vc192_320, &tusb1[192]);
	}
	else {
		memset(&v192_320uas, 0, sizeof v64_192uas);
		memset(vc192_320,255, sizeof vc64_192);
	}
	if (ntusb1 > 320) {// 320 448
		uint32_t ntua_448 = (ntusb1 > 448) ? 128 : ntusb1 - 320;
		BuildBaseAndCellsVector(ntua_448, v320_448uas, vc320_448, &tusb1[320]);
	}
	else {
		memset(&v320_448uas, 0, sizeof v64_192uas);
		memset(vc320_448, 255, sizeof vc64_192);
	}
	if (ntusb1 > 448) {// 448_576
		uint32_t ntua_576 = (ntusb1 > 576) ? 128 : ntusb1 - 448;
		BuildBaseAndCellsVector(ntua_576, v448_576uas, vc448_576, &tusb1[448]);
	}
	else {
		memset(&v448_576uas, 0, sizeof v64_192uas);
		memset(vc448_576, 255, sizeof vc64_192);
	}
	if (ntusb1 > 576) {// 576_704
		uint32_t ntua_704 = (ntusb1 > 704) ? 128 : ntusb1 - 576;
		BuildBaseAndCellsVector(ntua_704, v576_704uas, vc576_704, &tusb1[576]);
	}
	else {
		memset(&v576_704uas, 0, sizeof v64_192uas);
		memset(vc576_704, 255, sizeof vc64_192);
	}
	if (ntusb1 > 704) {// 704_832
		uint32_t ntua_832 = (ntusb1 > 832) ? 128 : ntusb1 - 704;
		BuildBaseAndCellsVector(ntua_832, v704_832uas, vc704_832, &tusb1[704]);
	}
	else {
		memset(&v704_832uas, 0, sizeof v64_192uas);
		memset(vc704_832, 255, sizeof vc64_192);
	}
	if (ntusb1 > 832) {// 832_960
		uint32_t ntua_960 = (ntusb1 > 960) ? 128 : ntusb1 - 832;
		BuildBaseAndCellsVector(ntua_960, v832_960uas, vc832_960, &tusb1[832]);
	}
	else {
		memset(&v832_960uas, 0, sizeof v64_192uas);
		memset(vc832_960, 255, sizeof vc64_192);
	}
	// _____ apply cells vectors to  band 1 step 3_6 or 3_7
	ntvb1go = (uint32_t)index_xy_b1.ntotvb;
	for (uint32_t iv = 0; iv < ntvb1go; iv++) {
		VALIDB64 vb = vab64b1[iv];
		register uint64_t bf = vb.bf;// must hit all "empty"
		ZS64 & w = zs64b1[iv];
		w.bf = bf;
		register uint64_t V = v64uas;// apply clues 2_5/6
		for (uint64_t j = 0; j < vb.nval; j++)
			V &= vc64[vb.tval[j]];
		w.v = V;
	}

}
int G17B::Apply_Band2_Step() {// prepare the main loop
	fb2 = index_xy_b2.and_g;// including more common cells if any
	acb2 = index_xy_b2.or_g;// all these cells can be active in the step
	fb12 = fb1 | fb2;
	acb12 = acb1 | acb2;

	// check dead branch

	register uint64_t Ra = acb12, filter = fb2;
	register uint64_t * tua = tusb1;
	uint32_t nua = ntusb1;
	for (uint32_t iua = 0; iua < nua; iua++) {
		register uint64_t Ru = tua[iua];
		if (!(Ru&filter)) {
			Ru &= Ra;
			if (!Ru) 	return 1; 
		}
	}

	//_________ setup the brute force start

	nclues_step = 0;
	uint64_t w = fb12;
	uint32_t xcell;
	stack_count_step.u64 = 0;
	while (bitscanforward64(xcell, w)) {
		w ^= (uint64_t)1 << xcell;
		uint32_t cell = From_128_To_81[xcell];
		tclues[nclues_step++] = cell;
		stack_count_step.u16[C_stack[cell]]++;
	}
	tcluesxy = &tclues[nclues_step];
	zh2b_i1.ValidXY_Step(tclues, nclues_step);

	// _____ apply cells vectors to  band 2 step no filter empty

	for (uint32_t iv = 0; iv < index_xy_b2.ntotvb; iv++) {
		VALIDB64 vb = vab64b2[iv];
		register uint64_t bf = vb.bf;
		ZS64 & w = zs64b2[iv];
		w.bf = bf;
		register uint64_t V = v64uas;
		for (uint64_t j = 0; j < vb.nval; j++)
			V &= vc64[vb.tval[j]];
		w.v = V;
	}

	// setup base vectors for more uas
	SetUpBaseVector(nclues_step, tclues, v64_192uas, vc64_192, bv192);
	SetUpBaseVector(nclues_step, tclues, v192_320uas, vc192_320, bv320);
	SetUpBaseVector(nclues_step, tclues, v320_448uas, vc320_448, bv448);
	SetUpBaseVector(nclues_step, tclues, v448_576uas, vc448_576, bv576);
	SetUpBaseVector(nclues_step, tclues, v576_704uas, vc576_704, bv704);
	SetUpBaseVector(nclues_step, tclues, v704_832uas, vc704_832, bv832);
	SetUpBaseVector(nclues_step, tclues, v832_960uas, vc832_960, bv960);

	// reduce the tua128 table in bands 3

	for (int ib = 0; ib < genb12.nband3; ib++) {// bands 3
		STD_B3 &myb = genb12.bands3[ib];
		myb.ntua128_b2 = 0;
		for (uint32_t i = 0; i < myb.ntua128; i++)
			if (!(myb.tua128[i].bf.u64[0] & filter))
				myb.tua128_b2[myb.ntua128_b2++] = myb.tua128[i]; 
	}

	return 0;
}

//___________________ gangster specific

void TGUAS::ApplyB2() {
	// relay table per size
	uint32_t ntt[21]; // count for tt
	nguasb2 = 0;
	BF128 tt[21][500];// guas 54 ua12 plus code 0_161 for the gangster
	memset(ntt, 0, sizeof ntt);

	register uint64_t Bf = g17b.fb12;
	register uint64_t Ac = g17b.acb12;

	for (uint32_t i = 0; i < ntgua_start; i++) {
		GUA & wg = tgua_start[i];
		guaw.Init(wg);
		for (uint32_t j = 0; j < wg.nua; j++) {// apply new subsets
			register uint64_t Ua = wg.tua[j];
			if (Ua & Bf) continue;
			Ua &= Ac;// could be empty
			uint64_t cc = _popcnt64(Ua);
			guaw.Adduacheck(Ua | (cc << 59)); // no redundancy
		}
		if (guaw.nua) {
			nguasb2++;
			// split per size in 54 + index 0_161 mode
			BF128 w;
			w.bf.u64[1] = wg.i81 + 81 * wg.type;
			for (uint32_t j = 0; j < guaw.nua; j++) {
				register uint64_t Ua = guaw.tua[j],
					cc = Ua >> 59,
					Ua1 = Ua & BIT_SET_27,
					Ua2 = Ua & BIT_SET_B2;
				if (cc > 16) continue; // should not be
				w.bf.u64[0] = Ua1 | (Ua2 >> 5);

				tt[cc][ntt[cc]++] = w;
			}
		}
	}

	nb64_1 = ((nguasb2 * 4) >> 6) + 1;
	//_____________ create vectors
	nvect = 0;
	for (uint32_t i = 0; i <= 16; i++) {
		uint32_t n = ntt[i];
		BF128 * tw = tt[i];
		for (uint32_t j = 0; j < n; j++) {
			BF128 w = tw[j];
			AddVect54(w.bf.u64[0], w.bf.u32[2]);
			if (nvect >= 784) break;
		}
		if (nvect >= 784) break;
	}
}
void TGUAS::ApplyFirst384() {
	memset(&bf162all, 0, sizeof bf162all);
	for (uint32_t iv = 0; iv < nb64_1; iv++) {
		TVG64 &vv = tvg64[iv];
		register uint64_t V = vv.v;
		for (int j = 0; j < g17b.nclues; j++)
			V &= vv.cells[g17b.tcluesxy[j]];
		register uint32_t cc64;
		while (bitscanforward64(cc64, V)) {
			V ^= (uint64_t)1 << cc64;// clear bit
			uint64_t i = vv.ti162[cc64], ibloc = i >> 6,
				bit = (uint64_t)1 << (i - 64 * ibloc);
			bf162all.bf[ibloc] |= bit;
		}
	}
}
void TGUAS::ApplyMore() {
	uint32_t nmore = (nvect + 63) >> 6;
	memset(&bf162more, 0, sizeof bf162more);
	for (uint32_t iv = nb64_1; iv < nmore; iv++) {
		TVG64 &vv = tvg64[iv];
		register uint64_t V = vv.v;
		for (int j = 0; j < g17b.nclues; j++)
			V &= vv.cells[g17b.tcluesxy[j]];
		register uint32_t cc64;
		while (bitscanforward64(cc64, V)) {
			V ^= (uint64_t)1 << cc64;// clear bit
			uint64_t i = vv.ti162[cc64], ibloc = i >> 6,
				bit = (uint64_t)1 << (i - 64 * ibloc);
			bf162more.bf[ibloc] |= bit;
		}
	}

	bf162more.New(bf162all);
	bf162all.Or(bf162more);
}
//______________main loop 64 filter on the first 64 uas
void G17B::Do64uas() {

	// 65 always if available
	if (index_xy_b1.titem[1].ntvb&& index_xy_b2.titem[0].ntvb) 
		DoChunk64(&zs64b1[index_xy_b1.titem[0].sum_vb], zs64b2,
			index_xy_b1.titem[1].ntvb, index_xy_b2.titem[0].ntvb);
	// 56 filtered at generation if not needed
	if (index_xy_b1.titem[0].ntvb&& index_xy_b2.titem[1].ntvb)
		DoChunk64(zs64b1, &zs64b2[index_xy_b2.titem[0].sum_vb],
			index_xy_b1.titem[0].ntvb, index_xy_b2.titem[1].ntvb);


}
void  G17B::DoChunk64(ZS64 * a, ZS64 * b, uint64_t na, uint64_t nb) {
	ZS64 * z1, *z2;
	uint64_t n1, n2;
	if (na < nb) { z1 = a; z2 = b; n1 = na; n2 = nb; }
	else { z1 = b; z2 = a; n1 = nb; n2 = na; }
	if ((nb * na) < 5000) Do64uas_11(z1, z2, n1, n2);
	else {// cut in chunks max Xchunk Ychunk
		uint64_t  ideb2 = 0, iend2 = YCHUNK64;
		if (iend2 > n2)iend2 = n2;
		while (ideb2 < n2) { //Y chunk
			uint64_t ny = iend2 - ideb2;
			uint64_t ideb1 = 0, iend1 = XCHUNK64;
			if (iend1 > n1)iend1 = n1;

			while (ideb1 < n1) {// X chunk
				uint64_t nx = iend1 - ideb1;
				Do64uas_11(&z1[ideb1], &z2[ideb2], nx, ny);
				ideb1 = iend1; iend1 += XCHUNK64;
				if (iend1 > n1)iend1 = n1;
			}
			ideb2 = iend2; iend2 += YCHUNK64;
			if (iend2 > n2)iend2 = n2;
		}
	}
}
inline void G17B::Do64uas_11(ZS64 * a, ZS64 * b, uint64_t na, uint64_t nb) {
	//check a matrix band 1 band2 for potential 2 bands valid 11 clues
	register ZS64 * Ra = &a[na - 1];
	register uint64_t * Rs = &to_clean[n_to_clean];
	for (; Ra >= a; Ra--) {
		register ZS64 * Rb = &b[nb - 1];
		register uint64_t va = (*Ra).v,  bfa = (*Ra).bf | fb12;
		for (; Rb >= b; Rb--) 	
			if (!(Rb->v&va))		*Rs++ = bfa | Rb->bf;
	}
	n_to_clean = Rs - to_clean;
	if (n_to_clean > 10000)CleanAll();
}


//___________ process potential valid bands 1+2

inline int Check128uas(uint32_t ncl, uint32_t *tcl, BF128  bv, BF128 * cellsv) {
	for (uint32_t i = 0; i < ncl; i++) {//all cells common to the step
		bv &= cellsv[tcl[i]];
	}
	return bv.isNotEmpty();
}
void G17B::CleanAll() {
	uint64_t nw = n_to_clean;
	n_to_clean = 0;
	if (!nw) return;
	if (aigstop)  return;
	p_cpt2g[5] += nw;
	for (uint64_t i = 0; i < nw; i++) {// loop on this chunk of XY
		if (aigstop) return;
		aigstopxy = 0;
		register uint64_t bf = to_clean[i];
		// setup clues over common clues
		wb12bf = bf;
		nclues = 0;
		uint64_t w = wb12bf ^ fb12;
		uint32_t xcell;
		stack_count = stack_count_step;
		while (bitscanforward64(xcell, w)) {
			w ^= (uint64_t)1 << xcell;
			uint32_t cell = From_128_To_81[xcell];
			tcluesxy[nclues++] = cell;
			stack_count.u16[C_stack[cell]]++;
		}
		if (stack_count.u16[0] > 6 || stack_count.u16[1] > 6 ||
			stack_count.u16[2] > 6) continue;

		// check UA's over 64 using clues specific
		if (ntusb1 > 64) {
			if (Check128uas(nclues, tcluesxy, bv192, vc64_192))continue;
			if (ntusb1 > 192) {
				if (Check128uas(nclues, tcluesxy, bv320, vc192_320))continue;
				if (ntusb1 > 320) {
					if (Check128uas(nclues, tcluesxy, bv448, vc320_448))continue;
					if (ntusb1 > 448) {
						if (Check128uas(nclues, tcluesxy, bv576, vc448_576))continue;
						if (ntusb1 > 576 && Check128uas(nclues, tcluesxy, bv704, vc576_704))continue;
						if (ntusb1 > 704 && Check128uas(nclues, tcluesxy, bv832, vc704_832))continue;
						if (ntusb1 > 832 && Check128uas(nclues, tcluesxy, bv960, vc832_960))continue;
					}
				}
			}
		}

		// check fifo additional tables
		if (moreuas_12_13.Check(bf))continue;
		if (moreuas_14.Check(bf))continue;
		if (moreuas_15.Check(bf))continue;
		if (moreuas_AB_small.Check(bf))continue;
		if (moreuas_AB.Check(bf)) continue;
		if (moreuas_AB_big.Check(bf)) continue;
		p_cpt2g[6]++;

		// no more filter to apply check each band 3

		uint32_t tvb3[256], nvb3 = 0; //Bands 3 still valid
		tguas.ApplyFirst384();// all short guas 
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
			if (genb12.bands3[ib3].Clean0())	tvb3[nvb3++] = ib3;
		if (!nvb3) continue;
		p_cpt2g[46]++;

		//____build mincount and filter on mincount in bands 3
		tguas.ApplyMore();
		if (bf162more.NotEmpty()) {
			uint32_t n = nvb3;
			nvb3 = 0;
			for (uint32_t iw = 0; iw < n; iw++) {
				uint32_t ib3 = tvb3[iw];
				if (genb12.bands3[ib3].Cleanmore())
					tvb3[nvb3++] = ib3;
			}
		}
		if (!nvb3) continue;
		p_cpt2g[47]++;

	//__________no more guas2 guas3 process bands 3	

		zhou[0].PartialInitSearch17(tclues, nclues + nclues_step);
		for (uint32_t iw = 0; iw < nvb3; iw++) {
			uint32_t ib3 = tvb3[iw];
			GoB3(genb12.bands3[ib3]);
			if (aigstopxy)break;// added ua in bands 1+2
		}
	}
}
/*
int G17B::Is_B12_Not_Unique() {
	p_cpt2g[7]++;
	myua = zh2b[0].Valid_XY(tcluesxy, nclues);
	if (myua) {//not unique do the best with the UA
		NewUaB12();
		return 1;
	}
	p_cpt2g[8]++;
	return 0;
}
*/


//_________________________________ apply guas in band3

//note : max per stack can be adjusted to the mode 566 or 656

int STD_B3::Clean0() {// critical code first 256 guas
	p_cpt2g[51]++;
	memset(&smin, 0, sizeof smin);
	register uint32_t cc64;
	register uint64_t V = bf162all.bf[0]& bf162b.bf[0];
	while (bitscanforward64(cc64, V)) {
		V ^= (uint64_t)1 << cc64;// clear bit
		Insert2(cc64);
	}
	V = bf162all.bf[1] &bf162b.bf[1];
	while (bitscanforward64(cc64, V)) {// 17 guas2;  47 guas 3
		V ^= (uint64_t)1 << cc64;// clear bit
		if(cc64<17) 		Insert2(cc64+64);
		else Insert3(cc64 -17);
	}
	V = bf162all.bf[2] & bf162b.bf[2];
	while (bitscanforward64(cc64, V)) {
		V ^= (uint64_t)1 << cc64;// clear bit
		Insert3(cc64 +47);
	}
	sminr = smin;
	smin.SetMincount();
	if (smin.mincount > 6)return 0;
	stack_count.u64 = g17b.stack_count.u64 + smin.Count_per_stack().u64;
	if (stack_count.u16[0] > 6 || stack_count.u16[1] > 6 ||
		stack_count.u16[2] > 6) return 0; // not ok
	return 1;
}
int STD_B3::Cleanmore() {// after  first 256 guas
	smin = sminr;
	register uint32_t cc64;
	register uint64_t V = bf162more.bf[0] & bf162b.bf[0];
	while (bitscanforward64(cc64, V)) {
		V ^= (uint64_t)1 << cc64;// clear bit
		Insert2(cc64);
	}
	V = bf162more.bf[1] & bf162b.bf[1];
	while (bitscanforward64(cc64, V)) {// 17 guas2;  47 guas 3
		V ^= (uint64_t)1 << cc64;// clear bit
		if (cc64 < 17) 		Insert2(cc64 + 64);
		else Insert3(cc64 - 17);
	}
	V = bf162more.bf[2] & bf162b.bf[2];
	while (bitscanforward64(cc64, V)) {
		V ^= (uint64_t)1 << cc64;// clear bit
		Insert3(cc64 + 47);
	}
	smin.SetMincount();
	if (smin.mincount > 6)return 0;
	stack_count.u64 = g17b.stack_count.u64 + smin.Count_per_stack().u64;
	if (stack_count.u16[0] > 6 || stack_count.u16[1] > 6 ||
		stack_count.u16[2] > 6) return 0; // not ok

	return 1;
}

//______________________ start process final b3
void G17B::GoB3(STD_B3 & b) {
	myband3 = &b;
	moreuas_b3.Init();
	memcpy(&genb12.grid0[54], b.band0, 4 * 27);
	memset(&hh0, 0, sizeof hh0);
	stack_countf = b.stack_count;
	smin = b.smin;
	register uint32_t Fstk = BIT_SET_27;
	if (stack_countf.u16[0] == 6) Fstk ^= 07007007;
	if (stack_countf.u16[1] == 6) Fstk ^= 070070070;
	if (stack_countf.u16[2] == 6) Fstk ^= 0700700700;
	wactive0 = fstk = Fstk;
	nuasb3_1 = nuasb3_2 = 0;
	{// load uas guas 2 guas 3 in field 
		register uint32_t cc64;
		register uint64_t V = bf162all.bf[0] & b.bf162b.bf[0];
		while (bitscanforward64(cc64, V)) {
			V ^= (uint64_t)1 << cc64;// clear bit
			uasb3_1[nuasb3_1++] = b.guas.ua_pair[cc64];
		}
		V = bf162all.bf[1] & b.bf162b.bf[1];
		while (bitscanforward64(cc64, V)) {// 17 guas2;  47 guas 3
			V ^= (uint64_t)1 << cc64;// clear bit
			if (cc64 < 17) 		uasb3_1[nuasb3_1++] = b.guas.ua_pair[cc64 + 64];
			else uasb3_1[nuasb3_1++] = b.guas.ua_triplet[cc64 - 17];
		}
		V = bf162all.bf[2] & b.bf162b.bf[2];
		while (bitscanforward64(cc64, V)) {
			V ^= (uint64_t)1 << cc64;// clear bit
			uasb3_1[nuasb3_1++] = b.guas.ua_triplet[cc64 + 47];
		}
	}

	{//____ reduce the ua128 table to active (no hint in bands 1+2)
		ntua128_b3 = 0;
		register uint64_t F = wb12bf;
		register uint32_t active_all = fstk | b.smin.critbf;
		for (uint32_t i = 0; i < b.ntua128_b2; i++)
			if (!(b.tua128_b2[i].bf.u64[0] & F)) {// forget uas hit in bands 1+2
				register uint32_t U = b.tua128_b2[i].bf.u32[2] & active_all;
				if (!U) return ; // This Ua will never be hit
				tua128_b3[ntua128_b3++] = U;// no can be hit by critical
			}
	}
	wg46.bf.u64[0] = bf162all.bf[0];
	wg46.bf.u64[1] = bf162all.bf[1];
	wg46 &= b.guas.isguasocketc2_46;// now mode 81 for guas 4_6
	if (b.smin.mincount == 6) {	GoB3_6(); 	return;	}
	if (b.smin.mincount == 5) {	GoB3_5();  	return;	}
	if (b.smin.mincount == 4) {	GoB3_4(); 	return;	}
	GoB3Expand();  
}

void G17B::GoB3Expand() {//mincount too low
	// fill "in field" with other uas and call expand
	int i81;
	register int  Rfilt = myband3->smin.critbf;
	for (uint32_t i = 0; i < ntua128_b3; i++) 
		uasb3_1[nuasb3_1++] = tua128_b3[i];	

	for (uint32_t i = 0; i < myband3->nua; i++) 
		uasb3_1[nuasb3_1++] = myband3->tua[i];
	
	while ((i81 = wg46.getFirst128()) >= 0) {
		wg46.Clear(i81);
		uasb3_1[nuasb3_1++] = myband3->guas.ua_pair[i81];
	}
	ExpandB3();
}

void G17B::GoB3_6() {//no room for outfield
	int i81;
	register int  Rfilt = myband3->smin.critbf;
	for (uint32_t i = 0; i < ntua128_b3; i++) {
		register uint32_t Ru = tua128_b3[i];
		if (Ru & Rfilt) uasb3_1[nuasb3_1++] = Ru;
		else return;
	}

	for (uint32_t i = 0; i < myband3->nua; i++) {
		register uint32_t Ru = myband3->tua[i];
		if (Ru & Rfilt) uasb3_1[nuasb3_1++] = Ru;
		else return;
	}
	while ((i81 = wg46.getFirst128()) >= 0) {
		wg46.Clear(i81);
		register uint32_t	Ru = myband3->guas.ua_pair[i81];
		if (Ru & Rfilt) uasb3_1[nuasb3_1++] = Ru;
		else return;
	}
	p_cpt2g[10]++;
	hh0.GoMiss0((*myband3));
}
void G17B::GoB3_5() {
	{// collect uas max 1 clue out field
		int i81;
		register int  Rfilt = myband3->smin.critbf;
		register uint32_t andout = fstk;
		noutmiss1 = 0;
		for (uint32_t i = 0; i < ntua128_b3; i++) {
			register uint32_t Ru = tua128_b3[i];
			if (Ru & Rfilt) uasb3_1[nuasb3_1++] = Ru;
			else { andout &= Ru;	noutmiss1 = 1; if (!andout) return; } //not ok		
		}

		for (uint32_t i = 0; i < myband3->nua; i++) {
			register uint32_t Ru = myband3->tua[i];
			if (Ru & Rfilt) uasb3_1[nuasb3_1++] = Ru;
			else { andout &= Ru;	noutmiss1 = 1; if (!andout) return; } //not ok		
		}
		while ((i81 = wg46.getFirst128()) >= 0) {
			wg46.Clear(i81);
			register uint32_t	Ru = myband3->guas.ua_pair[i81];
			if (Ru & Rfilt) uasb3_1[nuasb3_1++] = Ru;
			else { andout &= Ru;	noutmiss1 = 1; if (!andout) return; } //not ok		
		}
		andmiss1 = andout;
	}
	p_cpt2g[11]++;
	hh0.GoMiss1((*myband3));
}
void G17B::GoB3_4() {
	{// collect uas max 2 clues out field
		int i81;
		register int  Rfilt = myband3->smin.critbf;
		register uint32_t active_all = fstk | smin.critbf;
		register uint32_t andout = fstk;

		for (uint32_t i = 0; i < ntua128_b3; i++) {
			register uint32_t Ru = tua128_b3[i];
			if (Ru & Rfilt) uasb3_1[nuasb3_1++] = Ru;
			else { Ru&=active_all; 	andout &= Ru;	
			uasb3_2[nuasb3_2++] = Ru;	}
		}

		for (uint32_t i = 0; i < myband3->nua; i++) {
			register uint32_t Ru = myband3->tua[i];
			if (Ru & Rfilt) uasb3_1[nuasb3_1++] = Ru;
			else {	Ru &= active_all; 	andout &= Ru;
				uasb3_2[nuasb3_2++] = Ru;	}		
		}
		while ((i81 = wg46.getFirst128()) >= 0) {
			wg46.Clear(i81);
			register uint32_t	Ru = myband3->guas.ua_pair[i81];
			if (Ru & Rfilt) uasb3_1[nuasb3_1++] = Ru;
			else {	Ru &= active_all; 	andout &= Ru;
				uasb3_2[nuasb3_2++] = Ru;		}
		}
		if (nuasb3_2) {
			smin.minplus++;
			if (!andout) {//sminplus=6  must have outfield 2 clues valid -
				smin.minplus++;
				register uint32_t Ru = uasb3_2[0], bit, cc, andx = 0;
				while (bitscanforward(cc, Ru)) {// look for  possible cells
					bit = 1 << cc;
					Ru ^= bit;// clear bit
					andx = BIT_SET_27;
					for (uint32_t i = 1; i < nuasb3_2; i++) {
						register uint32_t Ru2 = uasb3_2[i];
						if (!(Ru2 & bit)) 	andx &= Ru2;
					}
					if (andx)break;
				}
				if (!andx) {//no 2 clues killing outfield
					// add outfield to infield and expand
					GoB3_4Expand();		return;		}
			}
			else {	GoB3_4Expand();	return;		}
		}
		else {	GoB3_4Expand();		return;	}
	}
	// 2 out field is possible assign first stacks locked 
	p_cpt2g[12]++;
	hh0.GoMiss2Init((*myband3));
	// start with the smallest ua next will be "and" of remaining uas
	uint32_t  uamin = uasb3_2[0];
	{
		register uint32_t min = _popcnt32(uamin);
		for (uint32_t i = 1; i < nuasb3_2; i++) {
			register uint32_t Ru = uasb3_2[i],cc= _popcnt32(Ru);
			if (cc < min) {		min = cc;		uamin = Ru;		}
		}
	}
	hh0.GoMiss2((*myband3), uamin);
}
void G17B::GoB3_4Expand() {
	p_cpt2g[14]++;
	memcpy(&uasb3_1[nuasb3_1], uasb3_2, nuasb3_2 * sizeof uasb3_1[0]);
	nuasb3_1 += nuasb3_2;
	ExpandB3();
}

//__________ phase 2___ find band 3 clues for one band 3

void G17B3HANDLER::GoMiss0(STD_B3 & b3) {
	smin = b3.smin;
	uasb3if = g17b.uasb3_1;
	nuasb3if = g17b.nuasb3_1;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	CriticalLoop();
}
void G17B3HANDLER::GoMiss1(STD_B3 & b3) {
	nmiss = 1;
	smin = b3.smin;
	stack_count = b3.stack_count;
	uasb3if = g17b.uasb3_1;
	nuasb3if = g17b.nuasb3_1;
	nuasb3of = g17b.noutmiss1;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	wactive0 = g17b.wactive0;
	wua = g17b.andmiss1;
	Do_miss1();
}
inline void G17B3HANDLER::AddCell_Miss2(uint32_t * t) {//uint32_t cell, int bit) {
	nuasb3of = t[2];
	wua = t[1]; 
	{
		register uint32_t cell = t[0];
		register int s = C_stack[cell];
		stack_count.u16[s]++;
		if (stack_count.u16[s] > 5) {
			s = ~(07007007 << (3 * s));// mask
			wua &= s;
			wactive0 &= s;
		}
	}
	nmiss--;
	known_b3 |= 1 << t[0];
	Do_miss1();
}
void G17B3HANDLER::Do_miss1(){
	if (!nuasb3of) {// subcritical in hn if solved
		int uabr = IsMultiple( active_b3);
		if (uabr) {// one ua outfield seen
			wua = uabr& wactive0;
		}
		else {// confirmed subcritical possible
			G17B3HANDLER hn = *this;
			hn.Go_Subcritical();
			wua &= ~active_b3; // don't re use this as first cell
		}
	}
	uint32_t res;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	while (bitscanforward(res, wua)) {
		int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
		G17B3HANDLER hn = *this;
		hn.AddCellMiss1(res, bit);
		hn.CriticalLoop();
	}
}
void G17B3HANDLER::GoMiss2Init(STD_B3 & b3) {
	smin = b3.smin;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	wactive0 = g17b.wactive0 & (BIT_SET_27 ^ active_b3);//  active cells out field
	stack_count = b3.stack_count;
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	for (int istack = 0, stp = 0111; istack < 3; istack++, stp <<= 1)
		if (stack_count.u16[istack] > 5) {// critical stack
			register int m2stack = stp & smin.mini_bf2, shrink = TblShrinkMask[m2stack];
			if (m2stack) {// common cell(s) to assign
				register int Mask = tbitsrows[shrink] << (3 * istack);
				//adjust count and known
				known_b3 |= Mask & (~smin.pairs27);// and set the common cell as assigned
				smin.mini_bf2 &= ~stp; // clear the 2pairs bit(s) in stack
				active_b3 &= (~Mask);// clear the  field bf
				smin.critbf &= (~Mask);
				smin.pairs27 &= (~Mask);
				smin.mincount -= _popcnt32(shrink);
			}
		}
}
void G17B3HANDLER::GoMiss2(STD_B3 & b3, uint32_t uamin) {
	nmiss = 2;
	uasb3if = g17b.uasb3_1;
	nuasb3if = g17b.nuasb3_1;
	uasb3of = g17b.uasb3_2;
	nuasb3of = g17b.nuasb3_2;
	if (!nuasb3of) {// subcritical in hn if solved
		wua = wactive0;
		int uabr = IsMultiple(active_b3| known_b3);
		if (uabr) {// one ua outfield seen
			wua = uabr;
		}
		else {// confirmed subcritical possible
			G17B3HANDLER hn = *this;
			hn.Go_Subcritical();
		}
	}
	else wua = uamin;
	// cells added must produce cells hitting all remaining uas
	uint32_t res, tcellsok[27][3], ntcellsok = 0 ;
	while (bitscanforward(res, wua)) {
		uint32_t   nout=0;
		register uint32_t  bit = 1 << res;
		wua ^= bit; 
		wactive0 ^= bit;
		register uint32_t andx = wactive0, s = C_stack[res];
		if (stack_count.u16[s] == 5) {
			s = ~(07007007 << (3 * s));// mask
			andx &= s;
		}
		for (uint32_t i = 0; i < nuasb3of; i++) {
			register uint32_t ua = uasb3of[i];
			if (!(ua&bit)) {
				nout = 1;	andx &= ua;	if (!andx) break;
			}
		}
		if (andx|| (!nuasb3of)) {
			tcellsok[ntcellsok][0] = res;
			tcellsok[ntcellsok][1] = andx;
			tcellsok[ntcellsok++][2] = nout; 

		}
	}

	for (uint32_t i = 0; i < ntcellsok; i++) {// now call 
		G17B3HANDLER hn = *this;
		hn.AddCell_Miss2(tcellsok[i]);
	}
}

void G17B3HANDLER::Critical2pairs() {// assign 2 pairs in minirow to common cell
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	if (smin.mini_bf2) {// and 2 pairs in minirow forced to common cell
		register int Rst = 07007007;// stack 0 pattern
		for (int ist = 0; ist < 3; ist++) {
			int shrink = TblShrinkMask[smin.mini_bf2 & (0111 << ist)];
			if (shrink) {// minirows 2 pairs in that stack
				register int Mask = tbitsrows[shrink] << (3 * ist);
				active_b3 &= (~Mask); // clear the minirow
				known_b3 |= Mask & (~smin.pairs27);// and set the common cell as assigned
			}
		}
		smin.mini_bf2 = 0;
	}
}
void G17B3HANDLER::CriticalLoop() {// after optional assignment
	while (1) {// first shrink uas in field
		irloop = 0;
		uint32_t * tn = &uasb3if[nuasb3if], n = 0;
		register uint32_t Ra = active_b3,
			Rfilt = known_b3;
		for (uint32_t iua = 0; iua < nuasb3if; iua++) {
			register int Ru = uasb3if[iua];
			if (Ru & known_b3) continue;// already hit, forget it
			Ru &= active_b3;
			if (!Ru) return;// dead branch
			if (_popcnt32(Ru) == 1) {// assign it and reduce the active cells
				CriticalAssignCell(Ru);
				Ra = active_b3; //can be  modified
				irloop = 1;// should loop for new singles
			}
			else tn[n++] = Ru;
		}
		uasb3if = tn;
		nuasb3if = n;
		if (!n) irloop = 0;// no need to loop again
		if (!irloop) break;
	}
	if (_popcnt32(known_b3) > 6) return;
	if (!active_b3) {// must be here expected number of clues
		if (nuasb3if) return; //can not be valid
		g17b.FinalCheckB3(known_b3);
		return; // branch closed
	}
	int wua = uasb3if[0] & active_b3, cell;
	while (bitscanforward(cell, wua)) {
		register int bit = 1 << cell;
		wua ^= bit;// clear bit

		// clean the bit in active_b3, this is now a dead cell downstream
		active_b3 ^= bit;
		G17B3HANDLER hn = *this;
		hn.CriticalAssignCell(bit);
		hn.CriticalLoop();
	}
}

void G17B::ExpandB3(){// uint32_t *tua, uint32_t nua) {// find all 5 and 6 clues solutions
	// Build tua
	p_cpt2g[13]++;
	uint32_t *tuaw=uasb3_1, nuaw = nuasb3_1;// use  in field pre loaded in the right way
	struct SPB3 {// spots to find band 3 minimum valid solutions

		// ====================== constant after initialization
		uint32_t  possible_cells, all_previous_cells, active_cells, iuab3,
			stack[3];
	}spb3[7], *s3, *sn3;
	s3 = spb3;
	s3->all_previous_cells = 0;
	s3->active_cells = BIT_SET_27;// all cells active
	// init the stack status
	for (int i = 0; i < 3; i++) {
		s3->stack[i] = stack_count.u16[i];// count before band 3 min count
		if (s3->stack[i] == 6)s3->active_cells &= ~(07007007 << (3 * i));
	}
	s3->iuab3 = 0; // copy the start table
	s3->possible_cells = tuaw[0] & s3->active_cells;
	int tcells[10];

	//____________ here start the search 6 clues
next:
	uint64_t ispot = s3 - spb3;
	{// catch and apply cell in bitfields
		register uint32_t cell, p = s3->possible_cells;
		if (!p)goto back;
		bitscanforward(cell, p);
		register int bit = 1 << cell;
		tcells[ispot] = cell;
		s3->possible_cells ^= bit;// clear bit
		register int filter = s3->all_previous_cells | bit,
			ac = s3->active_cells ^ bit;
		sn3 = s3 + 1; *sn3 = *s3; // (copy the stack count)
		sn3->all_previous_cells = filter;
		sn3->active_cells = s3->active_cells = ac;
		{// if the stack is limit update sn3 active
			int st = C_stack[cell];
			sn3->stack[st]++;
			if (sn3->stack[st] > 5)
				sn3->active_cells &= ~(07007007 << (3 * st));
		}
		// nextspot:take the next available ua to loop
		for (uint32_t i = s3->iuab3 + 1; i < nuaw; i++) {
			if (tuaw[i] & filter)continue;
			if (ispot >= 5) 	goto next;//passing the limit
			sn3->iuab3 = i;
			register uint32_t Ru = tuaw[i] & sn3->active_cells;
			if (!Ru)goto next;
			if (ispot == 4) {// last must hit all remaining uas
				for (uint32_t i2 = i + 1; i2 < nuaw; i2++) {
					if (tuaw[i2] & filter)continue;
					Ru &= tuaw[i2];
					if (!Ru)goto next;
				}
			}
			sn3->possible_cells = Ru;
			s3 = sn3; // switch to next spot
			goto next;
		}
	}
	if (ispot < 5) {// no more uas use active as possible
		sn3->possible_cells = sn3->active_cells;
		sn3->iuab3 = nuaw;
		s3 = sn3; // switch to next spot
		goto next;
	}
	p_cpt2g[30]++;	// this is a possible 17 do final check
	if (zhou[1].CallMultipleB3(zhou[0], sn3->all_previous_cells, 0)) {

		register uint32_t ua = zh_g2.cells_assigned.bf.u32[2];
		if (nuaw < 300)tuaw[nuaw++] = ua;
		NewUaB3();		
		if (!ua) return;// can now be empty (bands 1+2 not valid)
		s3->possible_cells &= ua;
	}
	else Out17(sn3->all_previous_cells);
	goto next;
	// going back, for a non empty index, count it back
back:
	if (--s3 >= spb3)goto next;
}



//=============== part 2  band 3 processing using guas2/3

int ZHOU::CallMultipleB3(ZHOU & o, uint32_t bf, int diagx) {
	*this = o;
	BF128 dca[9];
	int digitsbf = zh_g2.digitsbf;
	memcpy(dca, zh_g2.Digit_cell_Assigned, sizeof dca);
	{	
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = genb12.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc + 64; // the cell value in 3x32 of a 128 bits map
			if (FD[digit][0].Off(xcell))  return 0;// check not valid entry
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	if (_popcnt32(digitsbf < 8)) 	return 1;// can not be one solution
	
	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	//__________end assign last lot start solver
	zh_g.go_back = 0;	zh_g.nsol = 0; zh_g.lim = 1;// modevalid is set to  1
	int ir = Full17Update();
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0,0);

	return zh_g.nsol;  
}
int ZHOU::Apply17SingleOrEmptyCells() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched 
	BF128 R1 = FD[0][0], R2 = R1 & FD[1][0]; 	R1 |= FD[1][0];
	BF128 R3 = R2 & FD[2][0]; R2 |= R1 & FD[2][0]; R1 |= FD[2][0];
	BF128 R4= R3 & FD[3][0]; 
		R3 |= R2 & FD[3][0]; R2 |= R1 & FD[3][0]; R1 |= FD[3][0];
	BF128 R5 = R4 & FD[4][0]; R4 |= R3 & FD[4][0];
		R3 |= R2 & FD[4][0]; R2 |= R1 & FD[4][0]; R1 |= FD[4][0];
	R5 |= R4 & FD[5][0]; R4 |= R3 & FD[5][0];
		R3 |= R2 & FD[5][0]; R2 |= R1 & FD[5][0]; R1 |= FD[5][0];
	R5 |= R4 & FD[5][6]; R4 |= R3 & FD[6][0];
		R3 |= R2 & FD[6][0]; R2 |= R1 & FD[6][0]; R1 |= FD[6][0];
	R5 |= R4 & FD[7][0]; R4 |= R3 & FD[7][0];
		R3 |= R2 & FD[7][0]; R2 |= R1 & FD[7][0]; R1 |= FD[7][0];
	R5 |= R4 & FD[8][0]; R4 |= R3 & FD[8][0];
		R3 |= R2 & FD[8][0]; R2 |= R1 & FD[8][0]; R1 |= FD[8][0];
	if ((cells_unsolved - R1).isNotEmpty()) 	return 1; // empty cells
	R1 -= R2; // now true singles
	R1 &= cells_unsolved; // these are new singles
	if (R1.isEmpty()) {// no single store pairs and more
		zh_g.pairs = R2 - R3;
		zh_g2.triplets = R3 - R4;
		zh_g2.quads = R4 - R5;
		return 0;
	}
	int tcells[80], ntcells = R1.Table3X27(tcells);
	for (int i = 0; i < ntcells; i++) {
		int cell = tcells[i];
		for (int idig = 0; idig < 9; idig++) {
			if (FD[idig][0].On_c(cell)) {
				Assign(idig, cell, C_To128[cell]);
				goto nextr1;
			}
		}
		return 1; // conflict with previous assign within this lot
	nextr1:;
	}
	zh_g.single_applied = 1;
	return 0;
}
int ZHOU::Full17Update() {
	if (zh_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (Apply17SingleOrEmptyCells())	return 0; // locked empty cell or conflict singles in cells
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU::Guess17(int index, int diagx) {
	BF128 w = zh_g.pairs;
	if (w.isEmpty())w = zh_g2.triplets;
	if (w.isEmpty())w = zh_g2.quads;
	if (w.isEmpty())w = cells_unsolved;
	// here the target is to have ua band 3 as small as possible
	
	if (w.bf.u32[2]) w.bf.u32[0] = w.bf.u32[1] = 0;// select band 3 in priority
	int xcell = w.getFirst128(),
		cell = From_128_To_81[xcell],
		digit = zh_g2.grid0[cell];

	// true first if possible
	if (FD[digit][0].On(xcell)) {
		ZHOU * mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(digit, cell, xcell);
		mynext->Compute17Next(index + 1, 0);
		if (zh_g.go_back) return;

	}
	// if first step try first false
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (FD[idig][0].On(xcell)) {
			ZHOU * mynext = (this + 1);
			*mynext = *this;
			mynext->SetaCom(idig, cell, xcell);
			mynext->Compute17Next(index + 1, 0);
			if (zh_g.go_back) return;
		}
	}
}
void ZHOU::Compute17Next(int index, int diagx) {
	int ir = Full17Update();
	if (!ir) return;// locked 
	if (ir == 2) {//solved
		if (index) {// store false as ua
			BF128 & wua = zh_g2.cells_assigned;
			int * sol = genb12.grid0;
			wua.SetAll_0();;
			for (int i = 0; i < 81; i++) {
				int d = sol[i];
				if (FD[d][0].Off_c(i))	wua.Set_c(i);
			}
			if (wua.isNotEmpty()) {// ignore true solution
				zh_g.nsol++;
				zh_g.go_back = 1;// closed anyway
			}
		}
		return;
	}
	Guess17(index , 0);// continue the process
}

uint32_t G17B3HANDLER::IsMultiple(int bf) {
	if (bf == rknown_b3) return 0;
	uint32_t ua = 0;
	rknown_b3 = bf;
	G17B & bab = g17b;
	// check first if all tuab3 is hit
	p_cpt2g[55] ++;
	int ir = zhou[1].CallMultipleB3(zhou[0], bf, 0);
	if (ir) {
		ua = zh_g2.cells_assigned.bf.u32[2];
		g17b.moreuas_b3.Add(ua);
		g17b.NewUaB3();
	}
	return ua;
}

//================= critical process
void G17B3HANDLER::CriticalAssignCell(int Ru){// assign a cell within the critical cells
	// Ru is usually a regidster containing a 27 bits field with one bit on
	// 2 pairs in a miniriow have already been applied
	known_b3 |= Ru;
	uint32_t cell;
	bitscanforward(cell, Ru); // catch the cell
	register int mini = C_minirow[cell],// minirow to clear
		bit = 1 << mini,
		Mask = 7 << (3 * mini);
	if (bit & smin.mini_bf3){// the cell is in a minirow with 3 pairs active
		active_b3 &= ~Ru; //clear the cell
		smin.mini_bf3 ^= bit; // now only a pairto hit
		smin.mini_bf1 |= bit;
	}
	else{// either one pair or a triplet in the minirow
		active_b3 &= (~Mask); // kill the minirow as active
		smin.mini_bf1 &= ~bit;
		smin.mini_triplet &= ~bit;
	}
}

void G17B3HANDLER::Go_Critical(){// critical situation all clues in pairs tripl:ets
	active_b3 = smin.critbf;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	CriticalLoop();
}


//=============== sub critical process   missing(s)  in the critical area
void G17B3HANDLER::Go_SubcriticalMiniRow() {
	int c2[3] = { 3, 5, 6 };// 2 cells in a mini row
	int bit = 1 << ndead, mask = 7 << (3 * ndead);
	for (int i = ndead; i < 9; i++,  bit <<= 1, mask <<= 3) {
		stack = i % 3;
		register int M = active_sub & mask;
		if (!M)continue;
		ndead = i;
		if (bit & smin.mini_bf1) {// it was a gua2 pair assign both
			G17B3HANDLER hn = *this;
			hn.smin.mini_bf1 ^= bit;
			hn.SubMini( M, mask);
		}
		else if (bit & smin.mini_bf2)// it was 2 gua2 pair assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				G17B3HANDLER hn = *this;
				hn.smin.mini_bf2 ^= bit;
				hn.SubMini(M, mask);
			}
		else if (bit & smin.mini_bf3) {// it was 3 gua2 pair assign 3 out of 3
			G17B3HANDLER hn = *this;
			hn.smin.mini_bf3 ^= bit;
			hn.SubMini(M, mask);
		}
		else if (bit & smin.mini_triplet)// it was a gua3 triplet assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				G17B3HANDLER hn = *this;
				hn.smin.mini_triplet ^= bit;
				hn.SubMini(M, mask);
			}
		else {// second add in the mini row one residual cell take it
			G17B3HANDLER hn = *this;
			hn.SubMini(M, mask);
		}
	}
}
void G17B3HANDLER::SubMini( int M, int mask) {
	known_b3 |= M;// assign 1 or 2
	nmiss--;// one added
	active_b3 &= ~mask;
	active_sub ^= M;
	// now adjust the stack count
	stack_count.u16[stack]++;
	if (stack_count.u16[stack] > 5)active_sub &= ~(07007007 << (3 * stack));
	if (nmiss) Go_SubcriticalMiniRow();// continue till a"no missing clue condition"
	else {	// leave sub critical mode and enter the critical mode
		Critical2pairs();// assign 2 pairs in minirow to common cell
		CriticalLoop();
	}
}
void G17B3HANDLER::Go_Subcritical() {// nmiss to select in the critical field
	active_b3 = active_sub = smin.critbf;
	// check first if a global solution  is still possible
	for (int ist = 0; ist < 3; ist++) {// check stacks
		if (stack_count.u16[ist] > 5)active_sub &= ~(07007007 << (3 * ist));// kill the stack for more clues
	}
	ndead = 0;
	Go_SubcriticalMiniRow();// find the first miss
}




//________ final called by all branches
void G17B::FinalCheckB3(uint32_t bfb3) {
	p_cpt2g[29]++;
	if (moreuas_b3.Check(bfb3))return;
	register uint32_t ir = zhou[1].CallMultipleB3(zhou[0], bfb3, 0);
	if (ir) {
		register uint32_t ua = zh_g2.cells_assigned.bf.u32[2];
		NewUaB3();
		moreuas_b3.Add(ua);// if empty lock the call 
		return;
	}
	Out17(bfb3);
}
void G17B::Out17(uint32_t bfb3) {
	cout << Char27out(bfb3) << "\t\tone sol to print final check " << endl;
	char ws[82];
	strcpy(ws, empty_puzzle);
	for (int i = 0; i < (nclues_step+nclues); i++) {
		int cell = tclues[i];
		ws[cell] = genb12.grid0[cell] + '1';
	}
	for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)if (bfb3 & bit)
		ws[54 + i] = genb12.grid0[54 + i] + '1';
	fout1 << ws << ";" << genb12.nb12 / 64 << ";" << genb12.i1t16 << ";" << genb12.i2t16 << endl;
	a_17_found_here++;

}

void G17B::DebugAdd12() {
	aigstop = 1;
	cerr << "ua < 12 to add clean" << endl;
	cout << endl << endl << Char2Xout(myua) << " ua < 12 to add   clean" << endl;
	cout << "bug location band 2 id=" << genb12.nb12 << endl;
	cout << myband1.band << endl;
	cout << myband2.band << endl;
	cout << Char2Xout(wb12bf) << " b12 at call" << endl;
	cout << "ntusb1=" << ntusb1 << " n11=" << n_to_clean << endl;
	for (int i = 0; i < nclues_step; i++) cout << tclues[i] << " ";
	cout << "\t";
	for (int i = 0; i < nclues; i++) cout << tcluesxy[i] << " ";
	cout << endl;
	zh2b_i.ImageCandidats();
	zh2b_i1.ImageCandidats();
	zh2b[0].ImageCandidats();
	cout << "table uas" << endl;
	uint64_t *t = genuasb12.tua;
	uint32_t n = genuasb12.nua;
	for (uint32_t i = 0; i < n; i++) {
		uint64_t cc = _popcnt64(t[i] & BIT_SET_2X);
		if (cc > 12)break;
		cout << Char2Xout(t[i]) << " " << i << " " << cc << endl;

	}
	for (uint32_t i = 0; i < ntusb1; i++) {
		uint64_t cc = _popcnt64(tusb1[i] & BIT_SET_2X);
		if (cc > 12)break;
		cout << Char2Xout(tusb1[i]) << " b1 i=" << i << " " << cc << endl;

	}

}
void G17B::NewUaB12() {
	uint64_t cc64 = _popcnt64(myua&BIT_SET_2X);
	if (cc64 < 12) {// this should never be this is a check for a bug
		DebugAdd12();		return;
	}
	if (cc64 < 18) {
		if (cc64 < 14)moreuas_12_13.Add(myua);
		else if (cc64 == 14)moreuas_14.Add(myua);
		else if (cc64 == 15)moreuas_15.Add(myua);
		else moreuas_AB_small.Add(myua);
		register uint64_t ua_add = myua | (cc64 << 59);
		genuasb12.AddUACheck(ua_add);
		if(ntusb1<1000)tusb1[ntusb1++] = myua;
		p_cpt2g[31]++;
	}
	else if (cc64 < 21)			moreuas_AB.Add(myua);
	else moreuas_AB_big.Add(myua);
}

void G17B::NewUaB3_g2(uint32_t my_i81, uint64_t ua12) {
	int ix_start = tguas.ix_start.g2[my_i81];
	if (ix_start >= 0) { // should always be
		tguas.tgua_start[ix_start].Adduacheck(ua12);// for new steps
		uint64_t  ibloc = my_i81 >> 6, bit = (uint64_t)1 << (my_i81 - 64 * ibloc);
		bf162all.bf[ibloc] |= bit;
	}
	// add to vector 
	tguas.AddVect(ua12, 0, my_i81);
}
void G17B::NewUaB3_g3(uint32_t my_i81, uint64_t ua12) {
	int ix_start = tguas.ix_start.g3[my_i81];
	if (ix_start >= 0) {// should always be
		tguas.tgua_start[ix_start].Adduacheck(ua12);
		uint64_t  ibloc =( my_i81+81) >> 6, bit = (uint64_t)1 << (my_i81+81 - 64 * ibloc);
		bf162all.bf[ibloc] |= bit;
	}
	// add to vector 
	tguas.AddVect(ua12, 1, my_i81);
}

void G17B::NewUaB3() {// new ua from final check zh_g2.cells_assigned
	BF128 ua128 = zh_g2.cells_assigned;
	register uint64_t ua12 = ua128.bf.u64[0];
	register uint32_t ua = ua128.bf.u32[2],
		cc = _popcnt32(ua),
		cc0 = (uint32_t)_popcnt64(ua12);
	if (!cc) {// bands 1+2 not valid
		myua = ua12;
		NewUaB12();
		aigstopxy = 1;// don't process other bands 3
		return;
	}

	if (cc > 4) return; // see later if something of interest here
	if (cc0 > 16)return;
	p_cpt2g[32]++;// 2;3 or 4 cells in band 3

	// find the digits pattern from the current band 3
	int * cur_b3 = &genb12.grid0[54], wdigs = 0,c27;
	{
		register uint32_t wua = ua;
		while (bitscanforward(c27, wua)) {// look for  possible cells
			wua ^= 1 << c27;// clear bit
			wdigs |=1<<cur_b3[c27];
		}	
	}

	if (cc == 4) {// can be gua2
		register  uint32_t A = ua,
			B = (A | (A >> 9) | (A >> 18)) & 0777; // all columns
		uint32_t ncols=_popcnt32(B),ndigs= _popcnt32(wdigs);
		if (ncols == 3 && ndigs == 2) {
			p_cpt2g[40]++;
			uint32_t my_i81 = genb12.GET_I81_G2_4(wdigs, ua);
			NewUaB3_g2(my_i81, ua12);
			return;

		}
		else {
			if (cc0 > 12)return;// keep only small
			p_cpt2g[33]++;
			// could also be added later to other bands 3 where it is valid
			if (myband3->ntua128 < 1000) {
				myband3->tua128[myband3->ntua128++] = ua128;
				myband3->tua128_b2[myband3->ntua128_b2++] = ua128;
			}
			return;
		}
	}
	if (cc == 2) {// one of the 27 GUA2s add to the table
		p_cpt2g[40]++;
		uint32_t my_i81 = genb12.GET_I81_G2(wdigs, ua);
		NewUaB3_g2(my_i81, ua12);
		return;
	}
	if (cc == 3) {// one of the 3 GUA3s add to the table
		p_cpt2g[41]++;
		uint32_t my_i81 = genb12.GET_I81_G3(wdigs, ua);
		NewUaB3_g3(my_i81, ua12);
		return;
	}
}

