CC=gcc
CXX=g++
#CC=icc
#CXX=icpc
#CFLAGS  = -parallel -par-report3 -par-threshold0 -O3

CFLAGS  = -O2 -MMD

BINDIR = ./bin
LIBDIR = ./lib
#PHCDIR = ./PHC

# for k0ana0
#ROOT	= /usr/local/cern/root_v5.22.00/bin
# for k0ana1 & k0ana2
#ROOT	= /usr/local/cern/root_v5.34.36/bin
# for farm
#ROOT	= /usr/local/cern/root_v5.32.04.x86_64_fc8_gcc4.1.2/bin
ROOT	= /usr/local/cern/root6.20.04/bin


ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

INCLUDE= -I./ -I$(PHCDIR) -I../../

TARGET1=    Acceptance 
OBJS1  =    Acceptance.o
TARGET2=    vpflux 
OBJS2  =    vpflux.o
TARGET3=    acceptanceR 
OBJS3  =    acceptanceR.o
TARGET4=    vpflux_mthesis 
OBJS4  =    vpflux_mthesis.o
TARGET5=    acceptanceR_mthesis 
OBJS5  =    acceptanceR_mthesis.o
TARGET6=    acceptanceR_2022 
OBJS6  =    acceptanceR_2022.o
TARGET7=    acceptanceL_2022 
OBJS7  =    acceptanceL_2022.o

all: $(TARGET2)\
		$(TARGET3)\
		$(TARGET4)\
		$(TARGET5)\
		$(TARGET6)\
		$(TARGET7)

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

$(LIBDIR)/%.o : %.C
	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

#$(LIBDIR)/%.o : $(PHCDIR)/%.C
#	$(CXX) $(CFLAGS) $(INCLUDE) -c -o $@ $< $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET2): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET3): $(patsubst %,$(LIBDIR)/%,$(OBJS3))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET4): $(patsubst %,$(LIBDIR)/%,$(OBJS4))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET5): $(patsubst %,$(LIBDIR)/%,$(OBJS5))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET6): $(patsubst %,$(LIBDIR)/%,$(OBJS6))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET7): $(patsubst %,$(LIBDIR)/%,$(OBJS7))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o core $(BINDIR)/* $(LIBDIR)/*.d

-include $(DEPENDS6)
