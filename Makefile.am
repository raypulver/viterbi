bin_PROGRAMS = voxelizer viterbi expo
AM_CXXFLAGS = $(libpng_CFLAGS) $(zlib_CFLAGS) $(CFLAGS) -std=c++0x
AM_FFLAGS = -static-libgfortran
voxelizer_SOURCES = voxelizer.cc core.cc pdbplugin.c hmm.c cache.c
voxelizer_LDADD = $(libpng_LIBS) $(zlib_LIBS) -ljson-c -lm $(LDFLAGS)
viterbi_SOURCES = viterbi.cc core.cc pdbplugin.c hmm.c cache.c
viterbi_LDADD = $(libpng_LIBS) $(zlib_LIBS) -ljson-c -lm $(LDFLAGS)
expo_SOURCES = exponentiate.c exp.c hmm.c cache.c pdbplugin.c core.cc expokit/fortran/dgpadm.f expokit/fortran/blas.f expokit/fortran/lapack.f expokit/fortran/
expo_LDADD = $(libpng_LIBS) $(zlib_LIBS) -ljson-c -lm $(LDFLAGS) -lgfortran
