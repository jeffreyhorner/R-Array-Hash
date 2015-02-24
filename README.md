# R with Array Hashes

This project replaces R's hash tables for environments and the global
variable cache with Array Hashes[1]. The intention is to evaluate
run-time performance.

[1] - Fast and Compact Hash Tables for Integer Keys. Nicolas Askitis. ACSC '09
Proceedings of the Thirty-Second Australasian Conference on Computer Science -
Volume 91 Pages 113-122

# Configure and Build

1. You will want to download and install Google's Farmhash library, located here:
https://code.google.com/p/farmhash/. Step 2 below presumes you have installed it 
in /usr/local/.


2. Configure and build with the following:

For debugging and development:

  ```
  CC="ccache gcc"
  CFLAGS="-ggdb3 -pipe -std=gnu99 -Wall -pedantic -DPROTECT_PARANOID" \
  CXX="ccache g++"                                \
  CXXFLAGS="-ggdb -pipe -Wall -pedantic"          \
  FC="ccache gfortran"                            \
  F77="ccache gfortran"                           \
  LIBS="-lbfd" \
  ./configure --enable-memory-profiling --without-recommended-packages \
    --with-valgrind-instrumentation=2 --enable-R-shlib
  time make -j4
  ```

For speed (copied from ubuntu's R package):
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
    LIBS="-L/usr/local/lib -lfarmhash" \
    LIBnn=lib \
    ./configure  --with-cairo --with-jpeglib --with-readline --with-tcltk --with-system-bzlib \
    --with-system-pcre --with-system-zlib  --enable-R-profiling --enable-R-shlib \
    --enable-memory-profiling \
    --build x86_64-linux-gnu build_alias=x86_64-linux-gnu
  time make -j4
  ```
