MaxCliquePara general information (version 2.2)

The following assumes you are running a console in a directory where you have copied/uncompressed the source files.

1. Compiling
    to compile sources run:
    make
    to clean (remove the binary) run:
    make clean
    
2. Running
    run the binary executable:
    ./maxClique <graph_file_name> [number_of_threads] [job_queue_size]
    where 
        graph_file_name : either relative or absolute path to the file containing an undirected graph in DIMACS file format (http://prolland.free.fr/works/research/dsat/dimacs.html)
        number_of_threads (optional) : the number of threads to run the algorithm on. (0 - execute serially using serial code; 
            n .. execute on n threads using parallel code). If the argument is not specified, 0 is selected as default
        job_queue_size (optional) : size of the job queue. This argument is a remnant of development cycle in which a more complex code was used. Currently the only reasonable numbers is 1 (the default).
    The results are reported in text form on the standard output, listed in the following form:
        Loading DIMACS graph <graph_file_name>
        <|V|> vertices <|E|> edges <graph_density>
        -- pMC[<T> threads, <J> jobs](int,bitstring based set,Graph,MCR sort for Bitstrings<greedy color sort on bitstrings>) 
        search took <t>; <s> steps
        Clique (<|Q|>): [Q1,Q2,...,Qq]
        Thread efficiency = <Eff>
    where the variables filled by the algorithm are marked in <...>. Text of the third line starts with pMC if parallel code is used and MC if serial code is used.
    Time reported is for the clique search only (reading file from disk and transforming the graph into adjacency matrix are not counted)
    
