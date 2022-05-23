# Performance comparison with other tools
## Simulome
We tested the performance of Mutation-Simulator version 2 against Simulome,
simulating SNPs, insertions and deletions with a rate of 0.000825 per base each
on an [E.coli str. K-12](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3)
genome. In this test Mutation-Simulator performed 24 times faster than
Simulome.

Program | Parameters | Runtime \[s\] | Memory peak \[MB\]
--- | --- | --- | ---
Mutation-Simulator | ecoli_genome.fasta args -sn 0.000825353 -in 0.000825353 -de 0.000825353 | 1.7 | 106.34
Simulome | --genome=ecoli_genome.fasta --anno=ecoli_anno.gtf --output=output --snp=TRUE --num_snp=1 --whole_genome=TRUE --verbose=1 --indel=3 --ins_len=2 --num_ins=1 --del_len=2 --num_del=1 | 40.1 | 71.68

+ Testing Platform: Solus
+ Processor: Intel(R) Core(TM) i5-8265U CPU @ 1.60GHz
+ Mutation-Simulator version: 2
+ Measuring tool: [memory-profiler](https://pypi.org/project/memory-profiler/)

## SVsim
We also tested its performance in deletions against SVsim on the same E.coli
genome with a rate of 0.000098 deletions per base and it came out 207 times
faster.

Program | Parameters | Runtime \[s\] | Memory Peak \[MB\]
--- | --- | --- | ---
Mutation-Simulator | ecoli_genome.fasta args -de 0.000098241 -del 2 | 1.1 | 105.98
SVsim | -r ecoli_genome.fasta -o o/svsimfiles -i event | 228.0 | 22.34

+ Testing Platform: Solus
+ Processor: Intel(R) Core(TM) i5-8265U CPU @ 1.60GHz
+ Mutation-Simulator version: 2
+ Measuring tool: [memory-profiler](https://pypi.org/project/memory-profiler/)

## simuG
We also compared the performance of Mutation-Simulator version 2 against simuG
simulating SNPs and InDels at a rate of 0.01 and 0.1 on a 10 Mb test sequence
and Mutation-Simulator was 6785 times faster.

Program | Parameters | Runtime \[s\]
--- | --- | ---
Mutation-Simulator | 10Mb.fa args -sn 0.01 -in 0.005 -de 0.005 | 7
simuG | -refseq 10Mb.fa -snp_count 100000 -indel_count 100000 -prefix output_prefix | 47497
Mutation-Simulator | 10Mb.fa args -sn 0.1 -in 0.05 -de 0.05 | 47
simuG | -refseq 10Mb.fa -snp_count 1000000 -indel_count 1000000 -prefix output_prefix | aborted after several days

+ Testing Platform: Solus
+ Processor: Intel(R) Core(TM) i5-8265U CPU @ 1.60GHz
+ Mutation-Simulator version: 2
+ Measuring tool: [GNU time](https://www.gnu.org/software/time/)
