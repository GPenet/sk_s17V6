
#define _CRT_SECURE_NO_DEPRECATE
#define SEARCH17SOL


//#define DEFPHASE -4
#ifdef DEFPHASE
#endif

//#define DEBUGKNOWN 0
#ifdef DEBUGKNOWN
#else
#endif

//#define DEBUGONE 53
#ifdef DEBUGONE
#endif

//#define DEBUGINIT
#ifdef DEBUGINIT
#endif

//#define DEBUGEXL
#ifdef DEBUGEXL
#endif

//#define DEBUGL1L2 1
#ifdef DEBUGL1L2
#endif

//#define DEBUGSTEP 1
#ifdef DEBUGSTEP
#endif

//#define DEBUCLEANB 160656673
#ifdef DEBUCLEANB
#endif


/* program organisation
	main is the standard frame including the basic brute force 
*/

//#define MODE66_ON
#define CLEAN_SWITCH 8
#define XCHUNK64 200
#define YCHUNK64 200
#define XCHUNK128 200
#define YCHUNK128 200
#define XCHUNK256 100
#define YCHUNK256 100
#define GTEST17_ON 1
#define UALIMSIZE 20
#define GUALIMSIZE 18
#define UA32_10 0xffc00000
#define UA64_54 0x3fffffffffffff
#define TUA64_12SIZE 3000
/*entry 92
maxindex= 983
maxn5= 51516
maxn6= 237762
maxdet5= 261
maxdet6= 2004
*/
#define MAXN5 51520
#define MAXN6 237770 
#define MAX_56 300000 
#define MAXSTEP5 5000
#define MAXSTEP6 23000 
#define MAXNIND6 2100
#define MAXNIND5 300
// max for band 3 sockets 2,3,4,6
#define MAXSOCKB3 100

//============================================== 

#define G17MORESIZE 32

#define G17TESTUASGUASLIMITS 1

#include <sys/timeb.h>
#include "main.h"  // main and main tables and basic brute force
#include "go_17sol_tables.h"     
#include "Zh1b2b.h"  // brute force 2 bands  
extern SGO sgo;
//_________________ brute force handling 1 2 3 bands 
extern ZHOU    zhou[50],zhou_i;// , zhou_i, zhou_solve;
extern ZH_GLOBAL zh_g;
extern ZH_GLOBAL2 zh_g2;
extern ZH2B_GLOBAL   zh2b_g;
extern ZH2B5_GLOBAL   zh2b5_g;   // 2_5 digits 2 bands
extern ZH2B_1D_GLOBAL zh2b1d_g;  // one digit 2 bands
extern ZH2B zh2b[40], zh2b_i, zh2b_i1;
extern ZH2B5 zh2b5[10]; // solved digit per digit
extern ZH2B_1D zh2b1d[6]; // maxi 6 guesses, one per row
extern ZHONE_GLOBAL   zh1b_g;
extern ZHONE zhone[20];
extern ZHONE zhone_i;

FINPUT finput;
ofstream  fout1, fout2;

#include "sk_s17h.h"   //main classes of the project
#include "go_17sol_tables.h"
G17B g17b;
GENUAS_B12 genuasb12;
GEN_BANDES_12 genb12;
STD_B1_2 myband1, myband2;

//=== buffers to store valid bands and vectors
// max found in test all bands 243 steps 17063 valid bands
BI2 bi2_1[250], bi2_2[250];
VALIDB vab_1[MAX_56], vab_2[MAX_56];
VALIDB64 vab1_1[MAX_56],vab1a[MAX_56];

//uint64_t mainloop11[100000];
uint64_t to_clean[100000];

BI2 bi2_b1w[250], bi2_b1yes[250];
BI2 bi2_b2w[250], bi2_b2yes[250];
BI2 bi2_b1w2[250], bi2_b1yes2[250];
BI2 bi2_b2w2[250], bi2_b2yes2[250];
#define MAXEXP7 1200000
VALIDB vab1w[MAX_56], vab1yes[MAX_56];
VALIDB vab2w[MAX_56], vab2yes[MAX_56];
VALIDB vab1w2[MAX_56], vab1yes2[MAX_56];
VALIDB vab2w2[MAX_56], vab2yes2[MAX_56];

VALIDB64 vab64b1[MAX_56], vab64b2[MAX_56];
/*

*/


struct V3B128 {// for one bloc of 128 valid vand 3
	BF128 base;
	BF128 sockets [MAXSOCKB3];
	BF128 stacks[3][7];
};
struct VECTORBAND3 {// maximum of 1858 blocs of 128
	V3B128 b128[1858];
};

VECTORBAND3 tvectorband3[256];


uint32_t val3s[256][MAXN6];// expand up to 256 B3

GINT64 tempXY[30000];// limit chunkx * chunky here 100*200=20000
uint64_t valid_b12[30000];// in Clear tempxy




uint64_t p_cptg[40], p_cpt1g[20], p_cpt2g[60];
uint64_t p_cpt[40], p_cpt1[20];



#include "go_17_bands_cpp.h"  

#include "go_17_genb12_cpp.h"     
#include "go_17sol_bs_cpp.h"     
//#include "go_17sol_zx_cpp.h"  
#include "go_17sol_commands_cpp.h"



void Go_0() {
	// open  outputs files 1.txt
	if (sgo.foutput_name) {
		char zn[200];
		strcpy(zn, sgo.foutput_name);
		int ll = (int)strlen(zn);
		strcpy(&zn[ll], "_file1.txt");
		fout1.open(zn);
	}
	if (sgo.command >= 10
		&& sgo.command <20 ) {// input file expected
		if (!sgo.finput_name) {
			cerr << "missing input file name" << sgo.finput_name << endl; return;
		}
		finput.open(sgo.finput_name);
		if (!finput.is_open()) {
			cerr << "error open file " << sgo.finput_name << endl;
			return;
		}
	}
	cerr << "running command " << sgo.command << endl;
	switch (sgo.command) {
	case 0: Go_c17_00(); break; // search one band1
	case 9: Go_c17_09(); break; // search one band1 locate b2 b3
	case 10: Go_c17_10(); break; // search known 17s 
	case 11: Go_c17_11(); break; // split knwon 17s p2
	case 12: Go_c17_12(); break; // anal diag 665

	case 15: Go_c17_15(); break; // split knwon 17s 665 and others
	case 16: Go_c17_16(); break; // add solution+bands index


	case 80: Go_c17_80(); break; // enumeration test 
	case 90: Go_c17_90(); break; // regressive test uas one band	
	case 91: Go_c17_91(); break; // test uas collector 2 bands
	case 92: Go_c17_92(); break;// regressive test uas 5 6  expand
	}
	cerr << "go_0 return" << endl;
}

