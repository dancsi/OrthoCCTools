#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "Interaction.h"

using namespace std;

int main() {

	Interaction *inter = new Interaction();
	inter->init_complete_score(); 
	//cout << inter->weights["AA"]["AA"] << endl;
	//heptadi so ok
	/*for(int i=0;i<205;i++) {
		cout << inter->hi[i] << " ";
	}
	cout << endl;*/
	//cout << inter->ha.size() << " " << inter->duplets.size() << " " << inter->triplets.size() << endl;

	//DEIQALEEENAQLEQENAALEEEIAQLEYG
	cout << inter->score_complete("DEIQALEEENAQLEQENAALEEEIAQLEYG", "DKIAQLKQKIQALKQENQQLEEENAALEYG") << endl;
	cout << inter->score_complete("DEIQALEEENAQLEQENAALEEEIAQLEYG", "DENAALEEKIAQLKQKNAALKEEIQALEYG") << endl;
	cout << inter->score_complete("DEIQALEEENAQLEQENAALEEEIAQLEYG", "-KNAQLKEKIAALKEKIQQLKEENQALEYG") << endl;

	return 0;
}
