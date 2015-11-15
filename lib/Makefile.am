#
# Path and flags for all files used in building the library
#

AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
                @CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@


#semble
nobase_include_HEADERS =	semble/semble_algebra.h			\
                          semble/semble_file_management.h		\
                          semble/semble_histogram.h		\
                          semble/semble_histogram_aux.h   \
                          semble/semble_linear_algebra.h		\
                          semble/semble_matrix.h			\
                          semble/semble_meta.h			\
                          semble/semble_multi_t0fit.h		\
                          semble/semble_semble.h			\
                          semble/semble_typedefs.h		\
                          semble/semble_vector.h			

#
#Library
#

lib_LIBRARIES = libsemble.a

#semble 
libsemble_a_SOURCES = 	semble/semble_file_management.cc	\
                        semble/semble_histogram.cc		\
                        semble/semble_linear_algebra.cc \
                        semble/semble_algebra.cc
