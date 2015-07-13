%module fastscoreCC
 %{
#include "interaction.h"
 %}
 
 %include <cstring.i>
 
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
		//void heptad_array();
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
		//float score_rfe();
		//float score_fong_svm();
};

