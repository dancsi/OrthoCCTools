with open('transformed.txt', "r") as fin:
    with open('squared.fasta', 'w') as fout:
        peptides = set()
        for line in fin:
            line = line.strip()
            pair = line.split(',')
            peptides.add(pair[0][1:])
            peptides.add(pair[1][1:])
        peptides = list(peptides)
        print(peptides)
        pairs = {'-'+p+q for p in peptides for q in peptides}
        for idx, el in enumerate(pairs):
            print('>P'+str(idx), file=fout)
            print(el, file=fout)
            
