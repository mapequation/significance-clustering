Significance clustering


Author:
--------------------------------------------------------
Martin Rosvall

For contact information, see 
http://www.mapequation.org/about


Getting started:
--------------------------------------------------------
In a terminal with the GNU Compiler Collection installed,
just run 'make' in the current directory to compile the
code with the included Makefile.

Call: ./sigclu [-s <seed>] [-c <confidencelevel>] [-w <weightsfile>] partitionsfile outfile
seed: Any positive integer.
confidencelevel: The confidence as a fraction, default is 0.95.
partitionsfile: Each column represents a partition, the first for the raw partition and the remaining for bootstrap partitions. Row number corresponds to node id.
outfile: 1 or 0 if a node does or does not belong to the significant core of its module. moduleId1 moduleId2 means that the significant core of moduleId1 cooccurs with moduleId2 more than a fraction 1 - conf of the samples.
weightsfile: One column for weights of each node.  Row number corresponds to node id.

Example:

./sigclu -s 34 -c 0.95 -w weightssample.txt bootsample.txt signodessample.txt

bootsample.txt
1 1 1 1 1 1 1 1 1 1 1

1 1 1 1 1 1 1 1 1 1 1
1 2 1 2 1 2 1 2 1 2 2
2 2 1 2 1 2 1 2 1 2 2
2 2 2 2 2 2 2 2 2 2 2
2 2 2 2 2 2 2 2 2 2 2
3 2 2 2 2 2 2 2 2 2 2
3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3

weightssample.txt
1
1
2
2
3
3
4
4
5

signodessample.txt
# 1 or 0 if a node does or does not belong to the significant core of its module.
0
0
1
0
1
1
0
1
1
# moduleId1 moduleId2 means that the significant core of moduleId1 cooccurs with moduleId2 more than a fraction 1 - conf of the samples.
1 2