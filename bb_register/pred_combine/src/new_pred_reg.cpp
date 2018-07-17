#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "ec_score_dict.cpp"

using namespace std;

enum StrandOrientation {PERI2EXTRA,EXTRA2PERI};
enum HBondPattern {SWV,VWS};

// load odds data
void load_odds(const string fn, vector<double>& arr, double default_val, bool needlog=true);
// calculate pairing energy of a hairpin
double pairing(HBondPattern hbp, int len, StrandOrientation o1, vector<int> &strand1, vector<int> &strand2) ;
// get pairid for a given aa pair
int aa_pair(int a, int b);
// determine hbond pattern of two neighboring strands
HBondPattern patterning(int pdbi, int extra, int orientation1);
int getres(int pdbi, int seqid);
// predict regsistration
int regsistrationing();

const static int AABoundary = 19;
const static int InvalidAA = 20;
const static double nref = 8.5; // avg loop length of ori 25 prot

int pdb_num = 0;
vector<int> res_len;
vector< vector<int> > res;

vector<double> strong_scores(210,0), vdw_scores(210,0), weak_scores(210,0);
vector<double> CoreIn(21,0), CoreOut(21,0), PeriIn(21,0), PeriOut(21,0), ExtraIn(21,0), ExtraOut(21,0);

// Todo add Level argument
int main(int argc,char*argv[]) {

	if(argc!=6){
		cerr << "# of arguments incorrect"<<endl;
		cerr << "Usage: " << argv[0] << " test_file ec_reg_file odds_folder input_folder level" << endl;
		exit(1);
	}

	char* testfn = argv[1];
	char* ecsfn = argv[2];
	string oddsdir = argv[3];
	string inputsdir = argv[4];
	int level_select = atoi(argv[5]);
	cerr << "Running level " << argv[5]<<endl;

	if (level_select<1 ||level_select > 5){
	    cerr << "Specified level must be between 1 and 5"<<endl;
	    exit(1);
	}

	// read ec score
	ECSDict dict;
	get_ec_dict(ecsfn,dict);

	// load test input data
	ifstream fin(testfn);
	if(fin.fail()) {
		cerr << "error opening file " << testfn <<endl;
		exit(1);
	}
	vector<int> strand_num;
	vector< vector<int> > ori_peris, ori_extras, ori_regs;
	vector< vector<bool> > iscorrect;
	vector<string> pdbs;
	string line, pdbname;
	int strandnum;
	while(getline(fin, line)){
		stringstream ss(line);
		ss >> pdbname >> strandnum;
		pdbs.push_back(pdbname);
		strand_num.push_back(strandnum);
		ori_peris.push_back(vector<int>(strandnum,0));
		ori_extras.push_back(vector<int>(strandnum,0));
		ori_regs.push_back(vector<int>(strandnum,0));
		iscorrect.push_back(vector<bool>(strandnum, false));
		for(int strandi=0; strandi < strandnum; strandi++) {
			ss >> ori_peris.back()[strandi] >> ori_extras.back()[strandi] >> ori_regs.back()[strandi];
		}
	}
	fin.close();
	// get total pdb num
	pdb_num = pdbs.size();

	// load res data
	res_len = vector<int>(pdb_num);
	res = vector< vector<int> >(pdb_num);
	for(int pdbi=0 ; pdbi<pdb_num; pdbi++) {
		res[pdbi].push_back(InvalidAA); // index 0
		string fn = inputsdir + "//" + pdbs[pdbi] + "//" + pdbs[pdbi] + ".res";
		fin.open(fn.c_str());
		if(fin.fail()) {
			cerr << "error  opening file " << fn <<endl;
			exit(1);
		}
		int dummy;
		while(fin >> dummy){
			if(dummy > AABoundary){
				dummy = InvalidAA;
			}
			res[pdbi].push_back(dummy);
		}
		fin.close();
		res_len[pdbi]=res[pdbi].size()-1;
	}

	// load pair odds
	vector< vector<double> > strongs(pdb_num), vdws(pdb_num), weaks(pdb_num);
	string oddsfn;
	for(int pdbi = 0 ; pdbi<pdb_num; pdbi++){
		strongs[pdbi] = vector<double>(210,0);
		weaks[pdbi]   = vector<double>(210,0);
		vdws[pdbi]    = vector<double>(210,0);
		load_odds(oddsdir+"//default//strong.odds", strongs[pdbi], -1.72);
		load_odds(oddsdir+"//default//vdw.odds", vdws[pdbi], -1.72);
		load_odds(oddsdir+"//default//weak.odds", weaks[pdbi], -1.72);
	}

	// load single body odds
	load_odds( oddsdir+"//ExtraOut.odds", ExtraOut, -3.9 );
	load_odds( oddsdir+"//ExtraIn.odds",  ExtraIn, -3.9 );
	load_odds( oddsdir+"//CoreOut.odds", CoreOut, -3.9 );
	load_odds( oddsdir+"//CoreIn.odds",  CoreIn, -3.9 );
	load_odds( oddsdir+"//PeriOut.odds", PeriOut, -3.9 );
	load_odds( oddsdir+"//PeriIn.odds",  PeriIn, -3.9 );

	double wec, penub, penneg, wstrong, wvdw, wweak; // weights
	double wecl, penubl, pennegl, wstrongl, wvdwl, wweakl; // lower bound
	double wecu, penubu, pennegu, wstrongu, wvdwu, wweaku; // upper bound
	double wecs, penubs, pennegs, wstrongs, wvdws, wweaks; // step size

	//// as input is fed as group, the first pdb is used to determine the group in order to reduce searching efforts
	//// TODO 3emn not added yet!!!!

	//Changed here for level 1-5
//	const char* tmpg01[] = { "1bxw", "1qj8", "1p4t", "2f1t", "1thq", "2erv", "2lhf", "2mlh", "3dzm", "1qd6", "2f1c", "1k24", "1i78", "2wjr", "4pr7", "unkn" };
//	const char* tmpg02[] = { "1t16", "1uyn", "1tly", "3aeh", "3bs0", "3dwo", "3fid", "3kvn", "4e1s" };
//	const char* tmpg03[] = { "2mpr", "1a0s", "2omf", "2por", "1prn", "1e54", "2o4v", "3vzt", "4k3c", "4k3b", "4c4v", "4n75" };
//	const char* tmpg04[] = { "2qdz", "2ynk", "3rbh", "3syb", "3szv", "4c00", "4gey" };
//	const char* tmpg05[] = { "1fep", "2fcp", "1kmo", "1nqe", "1xkw", "2vqi", "3csl", "3rfz", "3v8x", "4q35" };
//	vector<string> g01(tmpg01, tmpg01+16);
//	vector<string> g02(tmpg02, tmpg02+ 9);
//	vector<string> g03(tmpg03, tmpg03+12);
//	vector<string> g04(tmpg04, tmpg04+ 7);
//	vector<string> g05(tmpg05, tmpg05+10);
	if( level_select ==1 ){
		// default
		wecl = 1.000; penubl = 0.070; pennegl = 0.040; wstrongl = 0.000; wvdwl = 0.000; wweakl = 0.010;
		wecu = 1.001; penubu = 0.581; pennegu = 0.161; wstrongu = 0.131; wvdwu = 0.071; wweaku = 0.101;
		wecs = 0.100; penubs = 0.020; pennegs = 0.010; wstrongs = 0.010; wvdws = 0.005; wweaks = 0.005;
		wecs = 0.100; penubs = 0.040; pennegs = 0.020; wstrongs = 0.020; wvdws = 0.010; wweaks = 0.010;//rand

		#ifdef ECPARAM110010080// ec 110010080
		wecl = 1.000; penubl = 0.070; pennegl = 0.040; wstrongl = 0.010; wvdwl = 0.000; wweakl = 0.010;
		wecu = 1.001; penubu = 0.121; pennegu = 0.081; wstrongu = 0.031; wvdwu = 0.031; wweaku = 0.041;
		wecs = 0.100; penubs = 0.005; pennegs = 0.005; wstrongs = 0.002; wvdws = 0.002; wweaks = 0.002;
		#endif
		#ifdef ECPARAM110030080// ec 110030080
		wecl = 1.000; penubl = 0.230; pennegl = 0.040; wstrongl = 0.010; wvdwl = 0.030; wweakl = 0.020;
		wecu = 1.001; penubu = 0.261; pennegu = 0.061; wstrongu = 0.041; wvdwu = 0.051; wweaku = 0.041;
		wecs = 0.100; penubs = 0.005; pennegs = 0.005; wstrongs = 0.002; wvdws = 0.002; wweaks = 0.002;
		#endif
		#ifdef ECPARAM110020070// ec 110020070
		wecl = 1.000; penubl = 0.070; pennegl = 0.040; wstrongl = 0.000; wvdwl = 0.000; wweakl = 0.020;
		wecu = 1.001; penubu = 0.581; pennegu = 0.161; wstrongu = 0.131; wvdwu = 0.071; wweaku = 0.101;
		wecs = 0.100; penubs = 0.020; pennegs = 0.010; wstrongs = 0.010; wvdws = 0.005; wweaks = 0.005;
		#endif
	}
	else if(  level_select ==2){
		// default
		wecl = 1.000; penubl = 0.350; pennegl = 0.100; wstrongl = 0.020; wvdwl = 0.060; wweakl = 0.040;
		wecu = 1.001; penubu = 0.581; pennegu = 0.161; wstrongu = 0.091; wvdwu = 0.151; wweaku = 0.101;
		wecs = 0.100; penubs = 0.020; pennegs = 0.005; wstrongs = 0.005; wvdws = 0.005; wweaks = 0.005;
		wecs = 0.100; penubs = 0.040; pennegs = 0.010; wstrongs = 0.010; wvdws = 0.010; wweaks = 0.010;//rand
		
		#ifdef ECPARAM110010080// ec 110010080
		wecl = 1.000; penubl = 0.370; pennegl = 0.100; wstrongl = 0.020; wvdwl = 0.070; wweakl = 0.050;
		wecu = 1.001; penubu = 0.581; pennegu = 0.141; wstrongu = 0.081; wvdwu = 0.151; wweaku = 0.101;
		wecs = 0.100; penubs = 0.020; pennegs = 0.005; wstrongs = 0.005; wvdws = 0.005; wweaks = 0.005;
		#endif
		#ifdef ECPARAM110030080// ec 110030080
		wecl = 1.000; penubl = 0.350; pennegl = 0.100; wstrongl = 0.030; wvdwl = 0.060; wweakl = 0.040;
		wecu = 1.001; penubu = 0.581; pennegu = 0.161; wstrongu = 0.091; wvdwu = 0.151; wweaku = 0.101;
		wecs = 0.100; penubs = 0.020; pennegs = 0.005; wstrongs = 0.005; wvdws = 0.005; wweaks = 0.005;
		#endif
		#ifdef ECPARAM110020070// ec 110020070
		wecl = 1.000; penubl = 0.390; pennegl = 0.120; wstrongl = 0.030; wvdwl = 0.080; wweakl = 0.040;
		wecu = 1.001; penubu = 0.581; pennegu = 0.161; wstrongu = 0.091; wvdwu = 0.151; wweaku = 0.101;
		wecs = 0.100; penubs = 0.010; pennegs = 0.005; wstrongs = 0.005; wvdws = 0.005; wweaks = 0.005;
		#endif
	}
	else if( level_select ==3 ){
		// default
		wecl = 1.000; penubl = 0.050; pennegl = 0.060; wstrongl = 0.000; wvdwl = 0.030; wweakl = 0.000;
		wecu = 1.001; penubu = 0.181; pennegu = 0.181; wstrongu = 0.031; wvdwu = 0.121; wweaku = 0.051;
		wecs = 0.100; penubs = 0.010; pennegs = 0.010; wstrongs = 0.002; wvdws = 0.005; wweaks = 0.005;
		wecs = 0.100; penubs = 0.020; pennegs = 0.020; wstrongs = 0.004; wvdws = 0.010; wweaks = 0.010;//rand
		
		#ifdef ECPARAM110010080// ec 110010080
		wecl = 1.000; penubl = 0.050; pennegl = 0.080; wstrongl = 0.000; wvdwl = 0.030; wweakl = 0.000;
		wecu = 1.001; penubu = 0.181; pennegu = 0.181; wstrongu = 0.031; wvdwu = 0.121; wweaku = 0.051;
		wecs = 0.100; penubs = 0.010; pennegs = 0.010; wstrongs = 0.002; wvdws = 0.005; wweaks = 0.005;
		wecs = 0.100; penubs = 0.005; pennegs = 0.005; wstrongs = 0.002; wvdws = 0.003; wweaks = 0.003;
		#endif
		#ifdef ECPARAM110030080// ec 110030080
		wecl = 1.000; penubl = 0.050; pennegl = 0.060; wstrongl = 0.000; wvdwl = 0.070; wweakl = 0.000;
		wecu = 1.001; penubu = 0.081; pennegu = 0.101; wstrongu = 0.011; wvdwu = 0.101; wweaku = 0.021;
		wecs = 0.100; penubs = 0.005; pennegs = 0.005; wstrongs = 0.002; wvdws = 0.002; wweaks = 0.002;
		wecs = 0.100; penubs = 0.002; pennegs = 0.002; wstrongs = 0.002; wvdws = 0.002; wweaks = 0.002;
		#endif
		#ifdef ECPARAM110020070// ec 110020070
		wecl = 1.000; penubl = 0.090; pennegl = 0.060; wstrongl = 0.000; wvdwl = 0.080; wweakl = 0.000;
		wecu = 1.001; penubu = 0.141; pennegu = 0.101; wstrongu = 0.021; wvdwu = 0.121; wweaku = 0.021;
		wecs = 0.100; penubs = 0.005; pennegs = 0.005; wstrongs = 0.002; wvdws = 0.002; wweaks = 0.002;
		wecs = 0.100; penubs = 0.002; pennegs = 0.002; wstrongs = 0.002; wvdws = 0.002; wweaks = 0.002;
		#endif
	}
	else if(  level_select ==4 ){
		// default
		wecl = 1.000; penubl = 0.290; pennegl = 0.100; wstrongl = 0.030; wvdwl = 0.000; wweakl = 0.010;
		wecu = 1.001; penubu = 0.561; pennegu = 0.220; wstrongu = 0.141; wvdwu = 0.061; wweaku = 0.071;
		wecs = 0.100; penubs = 0.020; pennegs = 0.005; wstrongs = 0.005; wvdws = 0.005; wweaks = 0.002;
		wecs = 0.100; penubs = 0.040; pennegs = 0.010; wstrongs = 0.010; wvdws = 0.010; wweaks = 0.004;//rand
		
		#ifdef ECPARAM110010080// ec 110010080
		wecl = 1.000; penubl = 0.510; pennegl = 0.200; wstrongl = 0.110; wvdwl = 0.040; wweakl = 0.040;
		wecu = 1.001; penubu = 0.541; pennegu = 0.220; wstrongu = 0.141; wvdwu = 0.061; wweaku = 0.071;
		wecs = 0.100; penubs = 0.005; pennegs = 0.005; wstrongs = 0.002; wvdws = 0.002; wweaks = 0.002;
		#endif
		#ifdef ECPARAM110030080// ec 110030080
		wecl = 1.000; penubl = 0.290; pennegl = 0.100; wstrongl = 0.030; wvdwl = 0.000; wweakl = 0.010;
		wecu = 1.001; penubu = 0.501; pennegu = 0.181; wstrongu = 0.101; wvdwu = 0.061; wweaku = 0.051;
		wecs = 0.100; penubs = 0.020; pennegs = 0.005; wstrongs = 0.005; wvdws = 0.005; wweaks = 0.002;
		#endif
		#ifdef ECPARAM110020070// ec 110020070
		wecl = 1.000; penubl = 0.470; pennegl = 0.140; wstrongl = 0.090; wvdwl = 0.010; wweakl = 0.020;
		wecu = 1.001; penubu = 0.561; pennegu = 0.181; wstrongu = 0.131; wvdwu = 0.041; wweaku = 0.051;
		wecs = 0.100; penubs = 0.005; pennegs = 0.005; wstrongs = 0.002; wvdws = 0.002; wweaks = 0.002;
		#endif
	}
	else if(  level_select ==5 ){
		// default
		wecl = 1.000; penubl = 0.050; pennegl = 0.080; wstrongl = 0.000; wvdwl = 0.000; wweakl = 0.000;
		wecu = 1.001; penubu = 0.301; pennegu = 0.181; wstrongu = 0.081; wvdwu = 0.061; wweaku = 0.041;
		wecs = 0.100; penubs = 0.020; pennegs = 0.010; wstrongs = 0.005; wvdws = 0.005; wweaks = 0.002;
		wecs = 0.100; penubs = 0.040; pennegs = 0.020; wstrongs = 0.010; wvdws = 0.010; wweaks = 0.004;//rand
		
		#ifdef ECPARAM110010080// ec 110010080
		wecl = 1.000; penubl = 0.050; pennegl = 0.080; wstrongl = 0.000; wvdwl = 0.000; wweakl = 0.000;
		wecu = 1.001; penubu = 0.301; pennegu = 0.181; wstrongu = 0.081; wvdwu = 0.061; wweaku = 0.041;
		wecs = 0.100; penubs = 0.020; pennegs = 0.010; wstrongs = 0.005; wvdws = 0.005; wweaks = 0.002;
		#endif
		#ifdef ECPARAM110030080// ec 110030080
		wecl = 1.000; penubl = 0.070; pennegl = 0.120; wstrongl = 0.020; wvdwl = 0.000; wweakl = 0.000;
		wecu = 1.001; penubu = 0.181; pennegu = 0.161; wstrongu = 0.071; wvdwu = 0.041; wweaku = 0.031;
		wecs = 0.100; penubs = 0.010; pennegs = 0.005; wstrongs = 0.005; wvdws = 0.002; wweaks = 0.002;
		#endif
		#ifdef ECPARAM110020070// ec 110020070
		wecl = 1.000; penubl = 0.050; pennegl = 0.120; wstrongl = 0.000; wvdwl = 0.020; wweakl = 0.000;
		wecu = 1.001; penubu = 0.161; pennegu = 0.161; wstrongu = 0.031; wvdwu = 0.061; wweaku = 0.031;
		wecs = 0.100; penubs = 0.010; pennegs = 0.005; wstrongs = 0.002; wvdws = 0.002; wweaks = 0.002;
		#endif
	}
	else{
		cerr << "unknown pdb " << pdbs[0] <<endl;
		exit(1);
	}

	//wecl = 11.000; penubl = 0.000; pennegl = 0.000; wstrongl = 0.000; wvdwl = 0.000; wweakl = 0.000;
	//wecu = 52.000; penubu = 0.900; pennegu = 0.900; wstrongu = 0.900; wvdwu = 0.900; wweaku = 0.900;
	//wecs = 5.000; penubs = 0.050; pennegs = 0.050; wstrongs = 0.050; wvdws = 0.050; wweaks = 0.050;

	#ifdef L1O
	for(int leavepdbi=0; leavepdbi<pdb_num; leavepdbi++)
	{
	#endif//L1O



	int max_correct_num = 0, max_correct_num_2nd = 0;

	// search for the best weights
	#ifdef BENCH
	wec = 0;
	#endif//BENCH
	#if defined(SEARCH) && ! defined(BESTWEIGHT) // search for best weights
	#ifndef BENCH // using only statistical potential
	for(wec = wecl; wec < wecu; wec += wecs)// p
	{
	#endif//BENCH
	for(penub   = penubl;   penub   < penubu;   penub   += penubs   )
	{
	for(penneg  = pennegl;  penneg  < pennegu;  penneg  += pennegs  )
	{
	for(wstrong = wstrongl; wstrong < wstrongu; wstrong += wstrongs )
	{
	for(wvdw    = wvdwl;    wvdw    < wvdwu;    wvdw    += wvdws    )
	{
	for(wweak   = wweakl;   wweak   < wweaku;   wweak   += wweaks   )
	{
	#endif//SEARCH
	
		int tot_correct_num = 0;

		// loop for each pdb
		for(int pdbi=0 ; pdbi < pdb_num; pdbi++) {
			#ifdef L1O
			if(leavepdbi==pdbi){
				continue;
			}
			#endif//L1O

            // TODO : Change here for level 1-5

			#if defined(BESTWEIGHT) && ! defined(SEARCH) // use best weight directly
			if(strand_num[pdbi] < 10 || level_select ==1){
				#ifdef ECPARAM110010080// ec 110010080
				wec=1; penub=0.075; penneg=0.045; wstrong=0.018; wvdw=0.014; wweak=0.02;
				#endif
				#ifdef ECPARAM110030080// ec 110030080
				wec=1; penub=0.245; penneg=0.05; wstrong=0.026; wvdw=0.038; wweak=0.036;
				#endif
				#ifdef ECPARAM110020070// ec 110020070
				wec=1; penub=0.07; penneg=0.05; wstrong=0.01; wvdw=0.01; wweak=0.03;
				#endif
			}
			else if(level_select ==2){
				#ifdef ECPARAM110010080// ec 110010080
				wec=1; penub=0.43; penneg=0.115; wstrong=0.05; wvdw=0.1; wweak=0.07;
				#endif
				#ifdef ECPARAM110030080// ec 110030080
				wec=1; penub=0.45; penneg=0.12; wstrong=0.055; wvdw=0.1; wweak=0.075;
				#endif
				#ifdef ECPARAM110020070// ec 110020070
				wec=1; penub=0.46; penneg=0.125; wstrong=0.055; wvdw=0.11; wweak=0.085;
				#endif
			}
			else if(level_select ==3){
				#ifdef ECPARAM110010080// ec 110010080
				wec=1; penub=0.05; penneg=0.08; wstrong=0; wvdw=0.054; wweak=0.024;
				#endif
				#ifdef ECPARAM110030080// ec 110030080
				wec=1; penub=0.052; penneg=0.074; wstrong=0; wvdw=0.082; wweak=0.006;
				#endif
				#ifdef ECPARAM110020070// ec 110020070
				wec=1; penub=0.122; penneg=0.084; wstrong=0.004; wvdw=0.1; wweak=0;
				#endif
			}
			else if(level_select ==4){
				#ifdef ECPARAM110010080// ec 110010080
				wec=1; penub=0.535; penneg=0.215; wstrong=0.128; wvdw=0.058; wweak=0.052;
				#endif
				#ifdef ECPARAM110030080// ec 110030080
				wec=1; penub=0.29; penneg=0.1; wstrong=0.045; wvdw=0.02; wweak=0.024;
				#endif
				#ifdef ECPARAM110020070// ec 110020070
				wec=1; penub=0.47; penneg=0.15; wstrong=0.096; wvdw=0.018; wweak=0.034;
				#endif
			}
			else if(strand_num[pdbi] < 28 || level_select ==5){
				#ifdef ECPARAM110010080// ec 110010080
				wec=1; penub=0.05; penneg=0.08; wstrong=0; wvdw=0.035; wweak=0.018;
				#endif
				#ifdef ECPARAM110030080// ec 110030080
				wec=1; penub=0.11; penneg=0.135; wstrong=0.045; wvdw=0.024; wweak=0.014;
				#endif
				#ifdef ECPARAM110020070// ec 110020070
				wec=1; penub=0.07; penneg=0.155; wstrong=0.01; wvdw=0.044; wweak=0.008;
				#endif
			}
			#endif// ! SEARCH && BESTWEIGHT

			for(int i=0; i<210; i++) {
				strong_scores[i] = strongs[pdbi][i] * wstrong;
				vdw_scores[i] = vdws[pdbi][i] * wvdw;
				weak_scores[i] = weaks[pdbi][i] * wweak;
			}

			// useless, this part is for optimization
			vector<int> peris  = ori_peris [pdbi];
			vector<int> extras = ori_extras[pdbi];
			vector<int> regs   = ori_regs  [pdbi];

			int curr_correct_num = 0; 

			#if defined(OUTPUT_SCORE) && !defined(SEARCH) // output score details for shear adjustment 
			cout << "## " << pdbs[pdbi] << endl;
			#endif//OUTPUT_SCORE

			// loop for each strand for current pdb
			for(int strandi = 0; strandi < strand_num[pdbi]; strandi++) {
				int strandj = (strandi+1) % strand_num[pdbi];
				vector<int> newstrand1, newstrand2;
				StrandOrientation orientation1;
				if(strandi%2==0){
					orientation1 = PERI2EXTRA;
					// make strands
					newstrand1 = vector<int> (res[pdbi].begin()+peris[strandi], res[pdbi].begin()+extras[strandi]+1);
					newstrand2 = vector<int> (res[pdbi].begin()+extras[strandj], res[pdbi].begin()+peris[strandj]+1);
					reverse(newstrand1.begin(),newstrand1.end());
				}
				else{
					orientation1 = EXTRA2PERI;
					// make strands
					newstrand1 = vector<int> (res[pdbi].begin()+extras[strandi], res[pdbi].begin()+peris[strandi]+1);
					newstrand2 = vector<int> (res[pdbi].begin()+peris[strandj], res[pdbi].begin()+extras[strandj]+1);
					reverse(newstrand2.begin(),newstrand2.end());
				}
				// determine pattern
				HBondPattern hbond_pattern = patterning(pdbi, extras[strandi], orientation1);

////////////////////////////////////////
//////////////////////////////////////// core computation part
				// make tmp strand2
				vector<int> newstrand2tmp(newstrand1.size());

				double maxscore = -1000;
				int max_predreg  =0;
				// enumerate registration
				//for(int regoffset=-20; regoffset <= 16; regoffset++) { 
				for(int regoffset=-10; regoffset <= 6; regoffset++) { 
					for(int i=0; i < newstrand1.size(); i++) {
						if(i+regoffset < 0 || i+regoffset >= newstrand2.size())
							newstrand2tmp[i] = InvalidAA;
						else
							newstrand2tmp[i] = newstrand2[i+regoffset];
					}

					int curr_predreg = newstrand2.size()-newstrand1.size()-regoffset;
					double ecscore = dict[pdbs[pdbi]][strandi][curr_predreg];
					double negative_reg = curr_predreg>=0 ? 0 : curr_predreg;
					double currscore =	pairing(hbond_pattern, newstrand1.size(), orientation1, newstrand1, newstrand2tmp) // hbond
										- penub*log((abs(regoffset) + nref)/nref) // penalty for unbonded res
										+ penneg * negative_reg // panalty for neg reg
										+ wec * ecscore; // ec score

					#if defined(OUTPUT_SCORE) && !defined(SEARCH)
					cout << setprecision(3) << curr_predreg << ":" << currscore << " ";

					#endif//OUTPUT_SCORE

					if(currscore > maxscore) {
						maxscore = currscore;
						max_predreg = curr_predreg;
					}
				}
				#if defined(OUTPUT_SCORE) && !defined(SEARCH)
				cout << endl;
				#endif//OUTPUT_SCORE
//////////////////////////////////////// core computation part
////////////////////////////////////////

				if( max_predreg == regs[strandi]) {
					iscorrect[pdbi][strandi] = true;
					curr_correct_num++;
				}
				else{
					iscorrect[pdbi][strandi] = false;
				}
				
				//cout << pdbs[pdbi] << " " << strandi << " " << max_predreg << " " << regs[strandi] << endl;
			}// strand loop
			tot_correct_num += curr_correct_num;
		}// pdb loop

		//if( tot_correct_num >= max_correct_num_2nd) { //this is for parameter searching
		if( tot_correct_num >= max_correct_num) { //this is for parameter searching
			if(tot_correct_num > max_correct_num){
				max_correct_num_2nd = max_correct_num;
				max_correct_num = tot_correct_num;
			}
			else if(tot_correct_num < max_correct_num){
				max_correct_num_2nd =  tot_correct_num;
			}
			#ifndef OUTPUT_SCORE
			#ifdef L1O
			cout << "l1o: " << leavepdbi << " | " << pdbs[leavepdbi] << " | ";
			#endif//L1O
			cout << wec << " " << penub << " " << penneg << " " << wstrong << " " << wvdw << " " << wweak << " " << tot_correct_num << endl;
			#endif//OUTPUT_SCORE

			for(int pdbi=0 ; pdbi < pdb_num; pdbi++) {
				#ifdef SEARCH
				cout << pdbs[pdbi] << " " ;
				#endif//SEARCH
				int pdbicorrect = 0;
				for(int strandi=0; strandi < strand_num[pdbi]; strandi++) {
					if(iscorrect[pdbi][strandi]){
						pdbicorrect++;
						//cout << "o";
					}
					else{
						//cout << "x";
					}
				}
				#ifndef OUTPUT_SCORE
				cout << pdbicorrect << " ";
				//cout << endl;
				#endif//OUTPUT_SCORE
			}
			#ifndef OUTPUT_SCORE
			cout << endl;
			#endif//OUTPUT_SCORE
		}

	
	#ifdef SEARCH
	}
	}
	}
	}
	}
	#ifndef BENCH
	}
	#endif//BENCH
	#endif//SEARCH

	#ifdef L1O
	}
	#endif//L1O
	return 0;
}

/// 
int aa_pair(int aa1, int aa2) {
	if(aa1 > AABoundary || aa2 > AABoundary ) { return -1; }
	if(aa1 <= aa2) { return (aa1 * 20 + aa2 - (aa1*(aa1+1)/2)); }
	else { return (aa2 * 20 + aa1 - (aa2*(aa2+1)/2)); }
}

double pairing(HBondPattern hbp, int len, StrandOrientation o1, vector<int> &strand1, vector<int> &strand2) {
	double sum = 0.0;
	vector<int> &tmpstrand1 = (o1==PERI2EXTRA) ? strand2 : strand1;
	vector<int> &tmpstrand2 = (o1==PERI2EXTRA) ? strand1 : strand2;
	if(hbp == SWV) {
		for(int i = 0 ; i < len-1 ; i++) {
			if(aa_pair(tmpstrand1[i], tmpstrand2[i+1]) != -1) {
				sum += weak_scores[aa_pair(tmpstrand1[i], tmpstrand2[i+1])];
			}
		}
		for(int i = 0 ; i < len ; i = i+2) {
			if(aa_pair(tmpstrand1[i], tmpstrand2[i]) != -1) {
				sum += strong_scores[aa_pair(tmpstrand1[i], tmpstrand2[i])];
			}
		}
		for(int i = 1 ; i < len ; i = i+2) {
			if(aa_pair(tmpstrand1[i], tmpstrand2[i]) != -1) {
				// will: gly no vdw? diff from pairing()
				if( tmpstrand1[i] != 7 || tmpstrand2[i] != 7)
					sum += vdw_scores[aa_pair(tmpstrand1[i], tmpstrand2[i])]; 
			}
		}
	}
	else {
		for(int i = 0 ; i < len-1 ; i++) {
			if(aa_pair(tmpstrand1[i], tmpstrand2[i+1]) != -1) {
				sum += weak_scores[aa_pair(tmpstrand1[i], tmpstrand2[i+1])];
			}
		}
		for(int i = 0 ; i < len ; i = i+2) {
			if(aa_pair(tmpstrand1[i], tmpstrand2[i]) != -1) {
				if( tmpstrand1[i] != 7 || tmpstrand2[i] != 7)
					sum += vdw_scores[aa_pair(tmpstrand1[i], tmpstrand2[i])]; 
			}
		}
		for(int i = 1 ; i < len ; i = i+2) {
			if(aa_pair(tmpstrand1[i], tmpstrand2[i]) != -1) {
				sum += strong_scores[aa_pair(tmpstrand1[i], tmpstrand2[i])];
			}
		}	
	}
	return sum;
}



void load_odds(const string fn, vector<double>& arr, double default_val, bool needlog){
	ifstream fin(fn.c_str());
	if(fin.fail()) {
		cerr << "error  opening file " << fn <<endl;
		exit(1);
	}
	int i = 0;
	while(fin >> arr[i]){
		if(arr[i]==0){ arr[i]=default_val; }
		else{ arr[i] = needlog ? (double)log(arr[i]) : (double)arr[i]; }
		i++;
	}
	fin.close();
}


int getres(int pdbi, int seqid){
	if( seqid<1 || seqid>res_len[pdbi] ){
		return InvalidAA;
	}
	return res[pdbi][seqid];
}


HBondPattern patterning(int pdbi, int extra, int orientation1){
	double SumEven, SumOdd;
	if(orientation1 == EXTRA2PERI) {
		SumEven = ExtraIn[getres(pdbi,extra)]+ExtraOut[getres(pdbi,extra+1)]+CoreIn[getres(pdbi,extra+2)]+CoreOut[getres(pdbi,extra+3)]+CoreIn[getres(pdbi,extra+4)]+CoreOut[getres(pdbi,extra+5)]+CoreIn[getres(pdbi,extra+6)]+PeriOut[getres(pdbi,extra+7)]+PeriIn[getres(pdbi,extra+8)];
		SumOdd = ExtraOut[getres(pdbi,extra)]+ExtraIn[getres(pdbi,extra+1)]+CoreOut[getres(pdbi,extra+2)]+CoreIn[getres(pdbi,extra+3)]+CoreOut[getres(pdbi,extra+4)]+CoreIn[getres(pdbi,extra+5)]+CoreOut[getres(pdbi,extra+6)]+PeriIn[getres(pdbi,extra+7)]+PeriOut[getres(pdbi,extra+8)];
	}
	else {
		SumOdd = ExtraIn[getres(pdbi,extra)]+ExtraOut[getres(pdbi,extra-1)]+CoreIn[getres(pdbi,extra-2)]+CoreOut[getres(pdbi,extra-3)]+CoreIn[getres(pdbi,extra-4)]+CoreOut[getres(pdbi,extra-5)]+CoreIn[getres(pdbi,extra-6)]+PeriOut[getres(pdbi,extra-7)]+PeriIn[getres(pdbi,extra-8)];
		SumEven = ExtraOut[getres(pdbi,extra)]+ExtraIn[getres(pdbi,extra-1)]+CoreOut[getres(pdbi,extra-2)]+CoreIn[getres(pdbi,extra-3)]+CoreOut[getres(pdbi,extra-4)]+CoreIn[getres(pdbi,extra-5)]+CoreOut[getres(pdbi,extra-6)]+PeriIn[getres(pdbi,extra-7)]+PeriOut[getres(pdbi,extra-8)];
	}
	if(SumOdd > SumEven){
		return SWV;
	}
	else{
		return VWS;
	}
}

