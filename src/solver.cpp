#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#define TRACE_LEVEL 1
#define TRACE_MASK (TRACE_MASK_CLIQUE|TRACE_MASK_INITIAL)

#include "mcqd_para/MaximumCliqueBase.h"
#include "mcqd_para/ParallelMaximumClique.h"
#include "mcqd_para/BB_GreedyColorSort.h"
#include "mcqd_para/BB_ColorRSort.h"
#include "mcqd_para/McrBB.h"

using namespace std;

const float c1 = -8.5, c2 = -7;

string fname;
int n_peptides;
float score[4100][4100];

vector<pair<int, int> > vertices;
vector<int> degrees;
vector<vector<char> > conn;

string str(pair<int, int> p)
{
	stringstream ss;
	ss<<"("<<p.first<<", "<<p.second<<")";
	return ss.str();
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

void print_clique(int* clique, int clique_size)
{
	stringstream ss;
	ss<<"[ "<<str(vertices[clique[0]]);
	for(int i=1;i<clique_size;i++)
	{
		ss<<", "<<str(vertices[clique[i]]);
	}
	ss<<" ]\n";
	//printf("%s", ss.str().c_str());
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
		if(argc>2)
		{
			sscanf(argv[2], "%d", &n_peptides);
		}
	}
	else
	{	
		/*
		cerr<<"Input file not given\n";
		exit(EXIT_FAILURE);
		*/
		fname = "random.out";
		n_peptides = 4096;
	}

	FILE* fin = fopen(fname.c_str(), "r");

	int id1, id2;
	float _score;
	while(fscanf(fin, "P%d,P%d,%f\n", &id1, &id2, &_score)==3)
	{
		score[id1][id2] = score[id2][id1] = _score;
	}
	fclose(fin);

	for(int i=1;i<=n_peptides;i++)
	{
		for(int j=i;j<=n_peptides;j++)
		{
			if(score[i][j]<=c1) vertices.push_back(make_pair(i, j));
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
	/*
	Maxclique mm(conn, n, print_clique);
	int clique_size;
	int * clique = new int[n];
	mm.mcqdyn(clique, clique_size);
	*/

	clock_t stop_time = clock();
	printf( "Total elapsed time: %.2lfs\n", double(stop_time - start_time)/CLOCKS_PER_SEC);

	return 0;
}