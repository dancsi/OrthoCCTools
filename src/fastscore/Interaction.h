#ifndef INTERACTION_H_
#define INTERACTION_H_

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include "fast_map_hack.h"

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

struct first_level_hash
{
	size_t operator()(const string& s) const
	{
		const char* a = s.c_str();
		switch (a[0])
		{
		case 'A':
			goto A;
		case 'D':
			goto D;
		case 'E':
			goto E;
		case 'G':
			goto G;
		}
	A:
		switch (a[1])
		{
		case 'A':
			return a[2] ? 1 : 0;
		case 'D':
			return a[2] ? 3 : 2;
		}
	D:
		switch (a[1])
		{
		case 'A':
			return 4;
		case 'D':
			return a[2]?6:5;
		case 'E':
			return a[2] ? (a[2] == 'A' ? 8 : 9) : 7;
		case 'G':
			return 10;
		}
	E:
		return a[2] ? 12 : 11;
	G:
		switch (a[1])
		{
		case 'A':
			return a[2] ? (a[2]=='D'?14:15) : 13;
		case 'D':
			return 16;
		case 'E':
			return 17;
		}
		return -1;
	}
};

struct second_level_hash
{
	size_t operator()(const string& s)
	{
		const int m[] = {0, 0, 1, 2, 3, 4, 5, 6, 7, 0, 8, 9, 10, 11, 0, 12, 13, 14, 15, 16, 0, 17, 18, 0, 19};
		const int N = 20;
		const char* a = s.c_str();
		const int b[3] = { a[0] - 'A', a[1] - 'A', a[2] - 'A' };
		if (a[2])
		{
			return m[b[0]] + N * (m[b[1]] + N * m[b[2]]);
		}
		else
		{
			return 20 * 20 * 20 + m[b[0]] + N*m[b[1]];
		}
	}
};

class Interaction {
	private:
		static const char hpos='f'; //kaze covjek da ne mijenjam
		static const int len = 100;
		string readFile;

	public:
		Interaction(string rFile){
			memset(hi,0,sizeof(hi));
			readFile = rFile;
		}
		void heptad_array();
		int hi[205];
		vector<char>ha;
		vector<triplet_list> triplets;
		vector<duplet_list>duplets;
		fast_map< fast_map<float, 20 * 20 * 20 + 20 * 20, second_level_hash>, 18, first_level_hash > weights;
		void init_complete_score(void);
		void get_heptad(void);
		void get_duplets(void);
		void get_triplets(void);
		float score_complete(string&,string&, int = 0);
		float score_rfe();
		float score_fong_svm();
};

#endif
