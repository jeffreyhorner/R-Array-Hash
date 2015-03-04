# R with Array Hashes

This project replaces R's hash tables for environments and the global
variable cache with Array Hash Tables for Integer Keys[1]. Also, the string cache and symbol table are implemented with an Array Hash Table for Strings[2]. The intention is to evaluate
run-time performance.

All changes to the R code base are licenced under the same terms as R using the GPL-2.

[1] - Askitis, Nikolas. "Fast and compact hash tables for integer keys." In Proceedings of the Thirty-Second Australasian Conference on Computer Science-Volume 91, pp. 113-122. Australian Computer Society, Inc., 2009.

[2] - Askitis, Nikolas, and Justin Zobel. "Redesigning the string hash table, burst trie, and bst to exploit cache." Journal of Experimental Algorithmics (JEA) 15 (2010): 1-7.

# Configure and Build

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
