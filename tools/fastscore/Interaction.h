#ifndef INTERACTION_H_
#define INTERACTION_H_

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

struct rNCR {
	int ii;
	char c1;
	char hhaa;
	int i;
	rNCR(){}
	rNCR(int a, char b, char c, int i1){
		ii = a; c1 = b; hhaa = c; i = i1;
	}
};

struct triplet_list {
	int i,j,k;
	string triad;
	triplet_list(){}
	triplet_list(int ii, int jj, int kk, string tr){
		i = ii; j = jj; k = kk;
		triad = tr;
	}
};

struct duplet_list {
	int i1,j1;
	string dup;
	duplet_list(){}
	duplet_list(int i, int j, string str){
		i1 = i; j1 = j;
		dup = str;
	}
};

class Interaction {
	private:
		char hpos; //kaze covjek da ne mijenjam
		int len;



	public:
		Interaction(){
			hpos = 'f';
			len = 100;
			memset(hi,0,sizeof(hi));
		}
		void heptad_array();
		int hi[205];
		vector<char>ha;
		vector<triplet_list> triplets;
		vector<duplet_list>duplets;
		multimap<string, map<string,float> >weights;
		void init_complete_score(void);
		void get_heptad(void);
		void get_duplets(void);
		void get_triplets(void);
		float score_complete(string,string);
		float score_rfe();
		float score_fong_svm();
};

#endif
