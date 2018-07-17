#ifndef EC_SCORE_DICT_CPP
#define EC_SCORE_DICT_CPP

// This is code to read ec score file, and store all the data in an ECSDict
// ECSDict : a dictionary-like data structure


#include <fstream>
#include <sstream>
//#include <iostream>
#include <string>
#include <unordered_map>
using namespace std;

// ECSDict[pdb][strandid][registration] = score
typedef unordered_map< string, unordered_map<int, unordered_map<int, double> > > ECSDict;

void get_ec_dict(char* fn,ECSDict& dict){
	string line;
	string pdb;
	ifstream fin(fn);
	if(fin.fail()){
		cerr << "Error in opening file " << fn << endl;
		exit(1);
	}
	while(getline(fin,line)){
		int strandi;
		if(line.size()==4){
			pdb = line;
			strandi = 0;
			dict[pdb] = unordered_map<int, unordered_map<int, double> >();
		}
		else{
			dict[pdb][strandi] = unordered_map<int, double>();
			stringstream ss(line);
			int regnum;
			ss >> regnum;
			int reg;
			double score;
			for(int i=0; i<regnum; i++){
				ss >> reg >> score;
				dict[pdb][strandi][reg]=score;
			}
			strandi++;
		}
	}
}

#endif//EC_SCORE_DICT_CPP
