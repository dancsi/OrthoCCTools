#include "Interaction.h"

#include <algorithm>

bool cmp(rNCR a, rNCR b) {
	return a.ii == b.ii ? a.c1 < b.c1 : a.ii < b.ii;
}

string uc(string str) {
	string ret = "";
	for (int i = 0;i < str.length();i++) {
		ret.push_back(toupper(str[i]));
	}
	return ret;
}

void Interaction::get_heptad(void) {
	ha.push_back('a');ha.push_back('b');ha.push_back('c');ha.push_back('d');
	ha.push_back('e');ha.push_back('f');ha.push_back('g');
	int i = 0;
	//nadji hpos
	//dodaj elemente na kraj
	while (ha[0] != hpos) {
		char tmp = ha[0];
		ha.erase(ha.begin());
		ha.push_back(tmp);
	}
	int count = 0;
	for (int i = 0;i < 200;i++) {
		if (ha[i % 7] == 'a') {
			count++;
		}
		hi[i] = count;
	}
}

void Interaction::get_duplets() {
	int cnt = 0;
	//cout << len << endl;
	for (int i = 0;i < len;i++) {
		for (int j = 0;j < len;j++) {
			if (abs(i - j) < 7) {
				//cout << i << " " << j << endl;
				//cnt++;
				char h1 = ha[i % 7];
				char h2 = ha[j % 7];
				int i1 = i;
				int j1 = j + len;
				if (i > j) {
					//swap
					swap(h1, h2);
					swap(i1, j1);
					//stop here
				}
				string temp; temp.push_back(toupper(h1));temp.push_back(toupper(h2));
				auto it = weights.find(temp);
				if (it != weights.end()) {
					//cout << "Dodajem " << temp << endl;
					duplets.push_back(duplet_list(i1, j1, temp));
				}
				else {
					cnt++;
					//cout << "Nisam nasao: " << temp << endl;
				}
			}
		}
	}
	//cout << len << endl;
	//cout << cnt << endl;
}

void Interaction::get_triplets(void) {
	int ii, jj, kk;
	int cnt = 0;
	vector<char>chain; chain.push_back('A'); chain.push_back('B');
	map<string, pair<int, int> > triad_set;
	triad_set["gde"] = {0, 1};
	triad_set["gae"] = {1, 1};
	triad_set["deg"] = {1, 0};
	triad_set["ega"] = {1, 0};
	triad_set["gad"] = {1, 0};
	triad_set["ade"] = {1, 0};
	triad_set["dga"] = {0, 1};
	triad_set["dea"] = {1, 1};
	triad_set["aad"] = {1, 1};
	triad_set["dda"] = {1, 1};
	for (int i = 0;i < len;i++) {
		for (int j = i + 1;j < 2 * len;j++) {
			for (int k = j + 1;k < 2 * len;k++) {
				if (k < len) continue;
				//cnt++;
				ii = i%len;jj = j%len;kk = k%len;
				if (abs(ii - jj) >= 7 || abs(ii - kk) >= 7 || abs(jj - kk) >= 7) continue;
				char c1 = chain[(int)(i / len)];
				char c2 = chain[(int)(j / len)];
				char c3 = chain[(int)(k / len)];
				if (c1 == c2 && c2 == c3) continue;
				vector<rNCR> temp;
				temp.push_back(rNCR(ii, c1, ha[ii % 7], i));
				temp.push_back(rNCR(jj, c2, ha[jj % 7], j));
				temp.push_back(rNCR(kk, c3, ha[kk % 7], k));
				sort(temp.begin(), temp.end(), cmp);
				//cout << temp[0].ii << " " << temp[0].c1 << " " << temp[0].hhaa << " " << temp[0].i << endl;
				//sorting je ok
				if (temp[0].ii == temp[1].ii && temp[0].c1 == temp[2].c1) {
					swap(temp[0], temp[1]);
				}
				string triad = "";
				triad.push_back(temp[0].hhaa); triad.push_back(temp[1].hhaa); triad.push_back(temp[2].hhaa);
				//string triad({temp[0].hhaa, temp[1].hhaa, temp[2].hhaa});
				//cout << triad << endl;
				//triad je tudi ok
				int b1 = (temp[0].c1 != temp[1].c1) ? 1 : 0;
				int b2 = (temp[0].c1 != temp[2].c1) ? 1 : 0;

				pair<int, int> what_to_find = { b1, b2 };

				auto it1 = triad_set.find(triad);
				if (it1 != triad_set.end()) {
					if (it1->second == what_to_find && weights.find(uc(triad)) != weights.end()) {
						cnt++;
						triplets.push_back(triplet_list(temp[0].i, temp[1].i, temp[2].i, uc(triad)));
					}
				}
			}
		}
	}
	//cout << "Izvrsio sam se " << cnt << endl;
}

float Interaction::score_complete(string& p1, string& p2)
{
	static char seq[2 * len + 1] = 
		"----------------------------------------------------------------------------------------------------" \
		"----------------------------------------------------------------------------------------------------";
#pragma omp threadprivate(seq)

	memcpy(seq, p1.c_str(), p1.length());
	memcpy(seq + len, p2.c_str(), p2.length());

	float score = 0.0f;
	string ab(2, 0);
	//score duplets
	for (int x = 0;x < duplets.size();x++) {
		int i = duplets[x].i1;
		int j = duplets[x].j1;
		string& xy = duplets[x].dup;
		//char tmp[3]; 
		//critical section
		//cout << seq.at(i) << " " << seq.at(j) << endl;
		//tmp[0] = seq.at(i); tmp[1] = seq.at(j);
		ab[0] = seq[i];
		ab[1] = seq[j];


		//cout << xy << " " << ab << endl;
		auto it = weights.find(xy);
		if (it != weights.end()) {
			//cout << "Nasao xy: " << xy << " trazim: " << ab << endl;
			auto it1 = (*it).second.find(ab);
			if (it1 != (*it).second.end()) {
				//cout << "Nasao ab: " << ab << endl;
				score += (*it1).second;
			}
		}
	}

	//duplets su u redu
	//score triplets
	int cnt = 0;
	ab.resize(3);
	for (int x = 0;x < triplets.size();x++) {
		int i = triplets[x].i;
		int j = triplets[x].j;
		int k = triplets[x].k;
		string& xy = triplets[x].triad;

		//char tmp[4]; tmp[0] = seq.at(i); tmp[1] = seq.at(j); tmp[2] = seq.at(k);
		ab[0] = seq[i];
		ab[1] = seq[j];
		ab[2] = seq[k];
		//cout << ab << endl;
		auto it = weights.find(xy);
		if (it != weights.end()) {
			cnt++;
			auto it1 = (*it).second.find(ab);
			if (it1 != (*it).second.end()) {
				score += (*it1).second;
			}
		}
	}
	//imas neki konstant term
	score += -4.54197; // constant term w0
	//cout << "Izvrsio se " << cnt << endl;
	return score;
}

void Interaction::init_complete_score(void) {
	//ucitaj iz fajla
	ifstream in(this->readFile);
	if (!in) {
		cout << "Could not open file" << endl;
		return;
	}
	string in1, in2;
	float value;
	while (in >> in1 >> in2 >> value) {
		//cout << in1 << " " << in2 << endl;
		auto it = weights.find(in1);
		if (it != weights.end()) {
			(*it).second.insert(make_pair(in2, value));
		}
		else {
			weights.insert(make_pair(in1, fast_map<float, 20*20*20 + 20*20, second_level_hash>()));
			it = weights.find(in1);
			(*it).second.insert(make_pair(in2, value));
		}
	}
	in.close();
	this->get_heptad();
	this->get_duplets();
	this->get_triplets();
}
