# R with Judy Arrays

This project replaces R's hash tables in environments, global
variable cache, and global string hash table with Judy Arrays. The
intention is to evaluate run-time performance.  For an explanation of
Judy arrays, see http://judy.sourceforge.net/.

# Configure and Build

1. You will want to download and install Judy Arrays from the link above. The Judy
library can be compiled to support either 32 or 64 bit architectures. Be
sure to match its configuration with R's architecture.

(It's presumed that you installed Judy into /usr/local/, thus the following CFLAGS
and LDFALGS variables use that location. If you installed it somewhere else you
will want to update those.)

2. Configure and buid R with the following:

  ```
  CFLAGS="-I/usr/local/include" \
  LDFLAGS="-L/usr/local/lib" \
  LIBS="-lJudy" \
  ./configure
  make
  ```

An alternative way to configure for debugging and development is:

  ```
  CC="ccache gcc"
  CFLAGS="-ggdb3 -pipe -std=gnu99 -Wall -pedantic -I/usr/local/include -DPROTECT_PARANOID" \
  CXX="ccache g++"                                \
  CXXFLAGS="-ggdb -pipe -Wall -pedantic"          \
  FC="ccache gfortran"                            \
  F77="ccache gfortran"                           \
  LDFLAGS="-L/usr/local/lib" \
  LIBS="-lJudy" \
  ./configure --enable-memory-profiling --without-recommended-packages \
    --with-valgrind-instrumentation=2   
  make
  ```

For speed:
  ```
    CC="gcc -std=gnu99" \
    CFLAGS="-O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g" \
    LDFLAGS="-Wl,-Bsymbolic-functions -Wl,-z,relro -L/usr/local/lib" \
    F77=gfortran  \
    FFLAGS="-g -O2 -fstack-protector --param=ssp-buffer-size=4"  \
    CXX=g++  \
    CXXFLAGS="-O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g" \
    FC="gfortran" \
    FCFLAGS="-g -O2 -fstack-protector --param=ssp-buffer-size=4" \
    LIBnn=lib \
    LIBS="-lJudy" \
    ./configure  --with-cairo --with-jpeglib --with-readline --with-tcltk --with-system-bzlib \
    --with-system-pcre --with-system-zlib  --enable-R-profiling --enable-R-shlib \
    --enable-memory-profiling \
    --build x86_64-linux-gnu build_alias=x86_64-linux-gnu LIBnn=lib
  make -j4
  ```

