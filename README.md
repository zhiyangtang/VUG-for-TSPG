# VUG-for-TSPG

1.'Appendix.pdf' contains the appendix of our paper, which includes some proofs that were omitted from the paper due to space constraints.

2.'VUG_code' contains the source codes of our paper and sample datasets and queries. 
To run the codes, please read the following instructions.

## Datasets

We provide three datasets in the folder "Datasets/" for testing, D1 to D3.

In addition, "Datasets/TestGraph.txt" is the example that appears in our paper. 

If you want to use other datasets, please copy them into "Datasets/".

## Generate Queries

Usage of query generation program in "Datasets/GenQuery/":

```C++
cd Datasets/GenQuery/
make clean
make
./GenerateQueries <Graph File> <Number of queries>
cd ../../
```

- Graph File: input graph filename in "Datasets/"
- Number of queries: the number of random queries

After executions, generated query files are stored as "Datasets/GenQuery/{Graph File}.query"

## Execute VUG

Usage of VUG main program in "TSPG/":

```C++
cd TSPG/
make clean
make
./Run <Graph File> <Query File>
cd ../
```

- Graph File: input graph filename in "Datasets/"
- Query File: input query filename in "Datasets/"

After executions, the logs including running time are written in "Results/Logs.csv"

The output edges (number of edges,all edge ids) for input queries are stored in "Results/Answers/{Query File}.answer"

Statistics (space cost, number of quick upper-bound edges, number of tight upper-bound edges and number of exact edges) for answering each query are stored in "Results/Statistics/{Query File}.csv"
