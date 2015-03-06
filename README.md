# R with Array Hashes

This project replaces R's hash tables for environments and the global
variable cache with Array Hash Tables for Integer Keys[1]. Also, the
string cache and symbol table are implemented with an Array Hash Table
for Strings[2]. The intention is to evaluate run-time performance.

All changes to the R code base are licensed under the same terms as R using the GPL-2.

# Trying with Docker

Docker images are available with a builds of R
with Array Hashes. The dockerfile is based on the
[rocker:r-devel](https://registry.hub.docker.com/u/rocker/r-devel/)
image. It contains a copy of the latest stable version of R as well so
that we can compare performance.

To start an interactive shell:

```bash
docker run -t -i jeffreyhorner/r-array-hash sh -c '/bin/bash'
```

Within the image, use the `R` or `RScript` commands to run stable R. To use R-devel with Array Hashes, use `RD` or `RDScript` commands.

```r
RScript mybenchmarks.R
RDScript mybenchmarks.R
```

# Building from source

Configure and build with either of the following:

1. For debugging and development:

```
CC="gcc" \
CFLAGS="-ggdb3 -pipe -std=gnu99 -Wall -pedantic -DPROTECT_PARANOID" \
CXX="g++"                                \
CXXFLAGS="-ggdb -pipe -Wall -pedantic"   \
FC="gfortran" F77="gfortran"             \
LIBS="-lbfd" \
./configure --enable-memory-profiling --with-valgrind-instrumentation=2
time make
```

2. For speed (copied from ubuntu's R):

```
CC="gcc -std=gnu99" \
CFLAGS="-O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g" \
LDFLAGS="-Wl,-Bsymbolic-functions -Wl,-z,relro" \
F77=gfortran  \
FFLAGS="-g -O2 -fstack-protector --param=ssp-buffer-size=4"  \
CXX=g++  \
CXXFLAGS="-O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g" \
FC="gfortran" \
FCFLAGS="-g -O2 -fstack-protector --param=ssp-buffer-size=4" \
LIBnn=lib \
./configure  --with-cairo --with-jpeglib --with-readline --with-tcltk --with-system-bzlib \
--with-system-pcre --with-system-zlib  --enable-R-profiling --enable-R-shlib \
--enable-memory-profiling \
--build x86_64-linux-gnu build_alias=x86_64-linux-gnu
time make
```

# Implementation

An array hash is a cache-conscious data structure that takes advantage of
hardware prefetchers for improved performance on large hash tables, those
large enough to fit in main memory and larger than fast fixed size cpu caches.

However, their implementation is a radical departure from standard
chained hash tables. Rather than using chains of hash buckets for
collision resolution, array hashes use segements of contiguous memory
called dynamic arrays to store keys and values. Adding and deleting
items from the hash involve copying the entire segment to new areas
in memory. While this may seem wasteful and slow, it's surprisingly
efficient in both time and space[2].

In R, hashed environments are implemented using lists with each list
element (a CONS cell) acting as the hash bucket. The CONS cell is the
binding agent for a symbol and value. Hashed environments are searched
using the pointer address of the symbol rather than the symbol's
printed name.

R-Array-Hash takes advantage of this by implementing an integer array
hash[1] to store addresses of symbols and their associated values. Care is
also taken to account for whether or not a binding is locked, active, etc.

Similarly, R-Array-Hash reimplements R's string cache using a string
array hash. This introduces the most radical change to R's API: CHAR()
no longer returns an address that points to the area at the end of the
SEXP. Rather it returns an address located in one of the contiguous
dynamic arrays of the string hash table. Therefore, care must be taken
in C code to use the address immediately since additions and deletions
to the string hash could render the result of CHAR() useless. There are
many areas of the code that sidestep this by calling translateChar(),
which has been changed to *always* copy the string pointed by CHAR().

# Bibliography

[1] - Askitis, Nikolas. "Fast and compact hash tables for integer
keys." In Proceedings of the Thirty-Second Australasian Conference on
Computer Science-Volume 91, pp. 113-122. Australian Computer Society,
Inc., 2009.

[2] - Askitis, Nikolas, and Justin Zobel. "Redesigning the string hash
table, burst trie, and bst to exploit cache." Journal of Experimental
Algorithmics (JEA) 15 (2010): 1-7.
