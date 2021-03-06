#include <algorithm>
#include <cstring>
#include <ctime>
#include <iostream>
#include <sstream>
#include <thread>

#include "common.h"
#include "options.h"
#include "ioutil.h"

#include "mcqd_para/MaximumCliqueBase.h"
#include "mcqd_para/ParallelMaximumClique.h"
#include "mcqd_para/BB_GreedyColorSort.h"
#include "mcqd_para/McrBB.h"

using namespace std;

float c1, c2;
bool search_homodimers, search_heterodimers;

string fname, initial_set_fname;
string out_name, fasta_name;

int n_peptides = 0;
float **score;

map<string, int> peptide_id;
vector<string> reverse_id;

vector<pair<int, int> > vertices, initial_set;
vector<int> degrees;
vector<vector<char> > conn;

Graph<BitstringSet> graph;

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

bool will_interact(pair<int, int> a, pair<int, int> b)
{
	int u1 = a.first, u2 = a.second, v1 = b.first, v2 = b.second;
	if(u1 == v1 || u1 == v2 || u2 == v1 || u2 == v2) return true;
	if(score[u1][v1]<c2 || score[u1][v2]<c2 || score[u2][v1]<c2 || score[u2][v2]<c2 ) return true;
	return false;
}

bool will_interact_with_initial(pair<int, int> potential_pair)
{
	return any_of(initial_set.begin(), initial_set.end(), 
		[&potential_pair](pair<int, int>& initial_pair) {
			return will_interact(potential_pair, initial_pair);
		});
}

int max_count = 0;
struct ProgressReporter
{
	static string cmdline;
	static void set_cmdline(int argc, char** argv) {
		cmdline = "#";
		for (int i = 0; i < argc; i++) {
			cmdline += argv[i] + " "s;
		}
		cmdline.pop_back();
	}

	void operator()(const BitstringSet& clique) {
		bool overwrite = false;
		if (clique.size() > max_count)
		{
			max_count = clique.size();
			overwrite = true;
		}
		print_clique(out_name.c_str(), clique, graph, cmdline, overwrite);
	};
};
std::string ProgressReporter::cmdline;

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

		string default_fasta_name = basename(fname) + ".fasta";
		fasta_name = options::get("fasta-name", default_fasta_name);

		string default_output_name = basename(fname) + ".pairs";
		out_name = options::get("out-name", default_output_name);
		ProgressReporter::set_cmdline(argc, argv);

		if(!(search_homodimers || search_heterodimers))
		{
			fprintf(stderr, "You can not disable both homo- and heterodimer search\n");
			exit(0);
		}
		initial_set_fname = options::get("initial-set", string(""));
	}

	score = read_scores(fname, fasta_name);
	cerr<<n_peptides<<" peptides\n";

	read_initial_set();

	if(search_homodimers)
	{
		for(int i=0;i<n_peptides;i++)
		{
			if(score[i][i]<=c1 && !will_interact_with_initial(make_pair(i, i))) vertices.push_back(make_pair(i, i));
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
				   score[j][j] >= c2 && !will_interact_with_initial(make_pair(i, j)))
					vertices.push_back(make_pair(i, j));
			}
		}
	}


	size_t n = vertices.size(), m = 0;
	if (n == 0)
	{
		cerr << "Orthogonal set impossible with given constraints, aborting.\n";
		exit(1);
	}
	else
	{
		cerr << "Running max_clique on " << n << " vertices\n";
	}

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

	graph.init(conn, degrees);

	ParallelMaximumCliqueProblem<
            int,                        // vertex ID
            BitstringSet,               // vertex set
            Graph<BitstringSet>,        // graph
            BBGreedyColorSort<Graph<BitstringSet>>,        // color sort
            BBMcrSort ,       // initial sort
            ProgressReporter  //user callback
	> problem(graph);

    int n_threads = thread::hardware_concurrency(), n_jobs = 2*n_threads;
    std::vector<int> affinities;
    printf("Running on %d threads, %d jobs\n", n_threads, n_jobs);

    problem.search(n_threads, n_jobs, affinities);
    problem.outputStatistics(false); std::cout << "\n";
    std::cout << "Thread efficiency = " << std::setprecision(3) << problem.workerEfficiency() << "\n\n";
    //print_clique(out_name, problem.getClique(), graph, "");

	clock_t stop_time = clock();
	printf( "Total elapsed time: %.2lfs\n", double(stop_time - start_time)/CLOCKS_PER_SEC);

	return 0;
}