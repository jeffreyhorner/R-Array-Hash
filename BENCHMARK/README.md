# List of Files

BENCHMARK.sh: driver file to run the tests
RUN.R: benchmark to run
results.RData: data file containing results data frame with 450 rows and allRs which contains info about each R tested.

The benchmark was set up similarly to the design in Nikolas Askitis and Justin Zobel's paper "Redesigning the String Hash Table, Burst Trie, and BST to Exploit Cache", easily found online. If you can't find it and want a copy, I'll send it.

Basically they created large data sets of words and tested hash table, Burst Trie, and Binary Search tree performance. In RUN.R we only test the performance of R environments as a proxy for a hash table.

There are two variants of the data sets, SKEW and DISTINCT. The SKEW data sets obey zipfs law, meaning there's tons of repeated words like 'a', 'an', and 'the'. Put another way the most frequent word in the set will appear twice as many times as the second most frequent word. I didn't know this beforehand, but apparently many corpi of spoken words obey this law. The DISTINCT data set contains only distinct words, no repeats, and they appear in the file in the order in which they were scraped, unordered.

So, what's measured? In the paper, they measured the total time it took to construct the hash table for each file, the total time it took to search the hash table, and the size in memory of the hash table. In RUN.R we measure those a little differently since we're dealing with an interpreted language and not C. Construct time is actually the amount of time it took call assign() for all words, search time is the amount of time it took to call get() for all words, and the run time measures the total amount of time to both construct and search the file. This includes the overhead of reading the files, gc() calls, hash table resizes etc. Also, we use calls to gc(reset=TRUE) and then another gc() to create a proxy for the size of the environment in memory.

The benchmark does 5 runs over 3 files ("SKEW.1mil","DISTINCT500thou","DISTINCT.1mil") over hash sizes of 2^10 to 2^19 over 3 different versions of R, so that's a total of 450 runs of RUN.R. Each file contains words scraped from Wikipedia articles, one word per line. The SKEW.1mil dataset has 1 million words with repeats, while the DISTINCT files contain no repeated words. Also none of the files are sorted; e.g they represent the position of the word in the articles as they were scraped. Mean string length of words in SKEW is 5.107329 characters and mean length of words in DISTINCT files is about 8.5.

Head of results (from results.RData):
                                        progname constructtime searchtime memory runtime hashsize cachemiss         datafile run
1 /home/jeffrey/Projects/R-Array-Hash/bin/exec/R         4.601      4.227   41.4  21.410     1024         0        SKEW.1mil   1
2        /home/jeffrey/Projects/devel/bin/exec/R         4.506      4.635   45.0  20.793     1024         0        SKEW.1mil   1
3 /usr/lib64/Revo-7.3/R-3.1.1/lib64/R/bin/exec/R         5.638      5.876   40.0  25.501     1024         0        SKEW.1mil   1
4 /home/jeffrey/Projects/R-Array-Hash/bin/exec/R         4.122      3.946  141.2  19.868     1024         0 DISTINCT.500thou   1
5        /home/jeffrey/Projects/devel/bin/exec/R         7.905      6.203  189.9  23.907     1024         0 DISTINCT.500thou   1
6 /usr/lib64/Revo-7.3/R-3.1.1/lib64/R/bin/exec/R        10.378      9.410  164.4  36.115     1024         0 DISTINCT.500thou   1

Ignore cachemiss as that's for measuring memory cache misses and that takes far too long to run the tests.

allRs (also from results.RData) is a list variable whose names() are:

[1] "/usr/lib64/Revo-7.3/R-3.1.1/lib64/R/bin/exec/R" "/home/jeffrey/Projects/devel/bin/exec/R"        "/home/jeffrey/Projects/R-Array-Hash/bin/exec/R"

and each list element has the Sys.info and R.version of the respective R's.
