bin_PROGRAMS = LINKS
LINKS_CXXFLAGS = $(AM_CXXFLAGS) -pthread

LINKS_CPPFLAGS = -I$(top_srcdir)/btllib/include/ 

#arcs_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz

#arcs_LDFLAGS = $(OPENMP_CXXFLAGS)

LINKS_SOURCES = \
	main.cpp \
	InputParser.hpp \
	LINKS.hpp