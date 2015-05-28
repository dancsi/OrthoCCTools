#include <algorithm>
#include <cstring>
#include <ctime>
#include <iostream>
#include <sstream>
#include <thread>
#include "options.h"

#define TRACE_LEVEL 1
#define TRACE_MASK (TRACE_MASK_CLIQUE|TRACE_MASK_INITIAL)

#include "mcqd_para/MaximumCliqueBase.h"
#include "mcqd_para/ParallelMaximumClique.h"
#include "mcqd_para/BB_GreedyColorSort.h"
#include "mcqd_para/McrBB.h"

using namespace std;

float c1, c2;
bool search_homodimers, search_heterodimers;

string fname;
string out_name;
int n_peptides = 0;
float score[4100][4100];

map<string, int> peptide_id;
vector<string> reverse_id;

vector<pair<int, int> > vertices;
vector<int> degrees;
vector<vector<char> > conn;

int get_id(string peptide)
{
	auto it = peptide_id.find(peptide);
	if(it!=peptide_id.end())
	{
		return it->second;
	}
	else
	{
		peptide_id.insert(make_pair(peptide, n_peptides));
		reverse_id.push_back(peptide);
		return n_peptides++;
	}
}

bool parse_input(FILE* fin, char* id1, char* id2, float& score)
{
	char buf[50];
	char *ret = fgets(buf, 49, fin);
	if(feof(fin) || !ret) return false;

	int len = strlen(buf);
	if(len<5) return false;

	int i=0;
	while(i <len &&buf[i]!=',')
	{
		*id1++ = buf[i];
		i++;
	}
	if(i==len) return false;
	*id1 = 0; i++;
	while(i<len && buf[i]!=',')
	{
		*id2++ = buf[i];
		i++;
	}
	if(i==len) return false;
	*id2 = 0; i++;
	score = atof(buf+i);
	return true;
}

bool will_interact(pair<int, int> a, pair<int, int> b)
{
	int u1 = a.first, u2 = a.second, v1 = b.first, v2 = b.second;
	if(u1 == v1 || u1 == v2 || u2 == v1 || u2 == v2) return true;
	if(score[u1][v1]<c2 || score[u1][v2]<c2 || score[u2][v1]<c2 || score[u2][v2]<c2 ) return true;
	return false;
}

void dump_dimacs(const char *fname, int n, int m)
{
	FILE* fout = fopen(fname, "w");
	fprintf(fout, "p edge %d %d\n", n, m);
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<i;j++)
		{
			if(conn[i][j])
				fprintf(fout, "e %d %d\n", i+1, j+1);
		}
	}
	fclose(fout);
}


void print_clique(const BitstringSet& clique, Graph<BitstringSet>& graph)
{
	stringstream ss;
    bool comma = false;
	BitstringSet bs = clique;

	graph.remap(bs); //VERY IMPORTANT!!!

	while(bs.size()>0)
    {
        int i = bs.nextSetBit();
		bs.remove(i);
        {
            auto p = vertices[i];
            ss<<reverse_id[p.first]<<","<<reverse_id[p.second]<<"\n";
        }
    }

	printf("%s", ss.str().c_str());

	FILE* fout = fopen("output.txt", "w");
	fprintf(fout,"%s", ss.str().c_str());
	fclose(fout);
}

int main(int argc, char** argv)
{
	clock_t start_time = clock();

	if(argc>1)
	{
		fname = string(argv[1]);
		options::parse(argc-1, argv+1);
	}
	else
	{
		options::usage(argv);
		exit(0);
	}

	{
		c1 = options::get("binding-cutoff", -8.5f);
		c2 = options::get("nonbinding-cutoff", -7.f);
		search_homodimers = !options::get("hetero-only", false);
		search_heterodimers = !options::get("homo-only", false);
		//out_name = options::get("out-name", "output.txt");
		if(!(search_homodimers || search_heterodimers))
		{
			fprintf(stderr, "You can not disable both homo- and heterodimer search\n");
			exit(0);
		}
	}

	FILE* fin = fopen(fname.c_str(), "r");

	int id1, id2;
	float _score;
	char id1s[20], id2s[20];

	while(parse_input(fin, id1s, id2s, _score))
	{
		id1 = get_id(id1s);
		id2 = get_id(id2s);

		score[id1][id2] = score[id2][id1] = _score;
	}
	fclose(fin);

	cerr<<n_peptides<<" peptides\n";

	if(search_homodimers)
	{
		for(int i=0;i<n_peptides;i++)
		{
			if(score[i][i]<=c1) vertices.push_back(make_pair(i, i));
		}
	}
	if(search_heterodimers)
	{
		for(int i=0;i<n_peptides;i++)
		{
			for(int j=i+1;j<n_peptides;j++)
			{
				if(score[i][j]<=c1 &&
				   score[i][i] >= c2 &&
				   score[j][j] >= c2)
					vertices.push_back(make_pair(i, j));
			}
		}
	}


	size_t n = vertices.size(), m = 0;

	degrees.resize(n, 0);
	conn.resize(n);
  	for (int i=0; i < n; i++) {
    	conn[i].resize(n, 0);
  	}

	for(int i=0;i<n;i++)
	{
		for(int j=0;j<i;j++)
		{
			if(!will_interact(vertices[i], vertices[j]))
			{
				conn[i][j] = conn[j][i] = true;
				degrees[i]++;
				degrees[j]++;
				m++;
			}
		}
	}

	Graph<BitstringSet> graph;
	graph.init(conn, degrees);

	ParallelMaximumCliqueProblem<
            int,                        // vertex ID
            BitstringSet,               // vertex set
            Graph<BitstringSet>,        // graph
            BBGreedyColorSort<Graph<BitstringSet>>,        // color sort
            BBMcrSort        // initial sort
            > problem(graph);

    int n_threads = thread::hardware_concurrency(), n_jobs = 2*n_threads;
    std::vector<int> affinities;
    printf("Running on %d threads, %d jobs\n", n_threads, n_jobs);

    problem.search(n_threads, n_jobs, affinities);
    problem.outputStatistics(false); std::cout << "\n";
    std::cout << "Thread efficiency = " << std::setprecision(3) << problem.workerEfficiency() << "\n\n";
    print_clique(problem.getClique(), graph);

	clock_t stop_time = clock();
	printf( "Total elapsed time: %.2lfs\n", double(stop_time - start_time)/CLOCKS_PER_SEC);

	return 0;
}