CC=gcc
CXX=g++

CFLAGS  = -O2

BINDIR = ./bin
LIBDIR = ./lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

# mp (Make Parameter)
#TARGET1= hrs_tof_twcorr
#OBJS1=   hrs_tof_twcorr.o 

#TARGET1= hrs_KID_F1_coin_run
#OBJS1=   hrs_KID_F1_coin_run.o  
#TARGET2= hrs_KID_F1_NOK
#OBJS2=   hrs_KID_F1_NOK.o  
#TARGET3= hrs_path_corr
#OBJS3=   hrs_path_corr.o  
#TARGET4= hrs_KID_F1_eff2
#OBJS4=   hrs_KID_F1_eff2.o Setting.o  
#TARGET5= hrs_KID_F1_check
#OBJS5=   hrs_KID_F1_check.o  Setting.o 
#TARGET6= hrs_MM_F1
#OBJS6=   hrs_MM_F1.o  
#TARGET7= hrs_KID_F1_acc
#OBJS7=   hrs_KID_F1_acc.o  
#TARGET8= hrs_KID_F1_eff3
#OBJS8=   hrs_KID_F1_eff3.o Setting.o 
#TARGET9= hrs_KID_F1_eff4
#OBJS9=   hrs_KID_F1_eff4.o Setting.o 
#TARGET10= hrs_KID_F1_eff5
#OBJS10=   hrs_KID_F1_eff5.o Setting.o 
#TARGET11= hrs_KID_F1_eff6
#OBJS11=   hrs_KID_F1_eff6.o Setting.o
#TARGET11= hrs_coin
#OBJS11=   hrs_coin.o Setting.o
#TARGET12= hrs_tuingAC
#OBJS12=   hrs_tuningAC.o Setting.o
TARGET13= tuning_ac_oku
OBJS13=   tuning_ac_oku.o Setting.o Tree.o ParamMan.o
TARGET14= eff_z
OBJS14=   eff_z.o Setting.o Tree.o ParamMan.o
TARGET15= eff_fp
OBJS15=   eff_fp.o Setting.o Tree.o ParamMan.o
TARGET16= eff_chi2
OBJS16=   eff_chi2.o Setting.o Tree.o ParamMan.o
TARGET17= eff_z2d
OBJS17=   eff_z2d.o Setting.o Tree.o ParamMan.o
TARGET18= eff_coin
OBJS18=   eff_coin.o Setting.o Tree.o ParamMan.o
TARGET19= eff_ac
OBJS19=   eff_ac.o Setting.o Tree.o ParamMan.o
TARGET20= z_analyzer 
OBJS20=   z_analyzer.o Setting.o Tree.o ParamMan.o
TARGET21= h2all 
OBJS21=   h2all.o Setting.o Tree.o ParamMan.o
#TARGET14= okuyana 
#OBJS14=   okuyana.o Setting.o Tree.o ParamMan.o

.PHONY: all
all:$(TARGET13)\
     $(TARGET20)\
     $(TARGET21)
#     $(TARGET11)
#$(TARGET1)\
#     $(TARGET2)\
#     $(TARGET3)\
#     $(TARGET4)\
#     $(TARGET5)\
#     $(TARGET6)\
#     $(TARGET7)\
#     $(TARGET8)\
#     $(TARGET9)\
#     $(TARGET10)\



$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< -MMD $(CXXFLAGS)

#$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
#$(TARGET2): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
#$(TARGET3): $(patsubst %,$(LIBDIR)/%,$(OBJS3))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
#$(TARGET4): $(patsubst %,$(LIBDIR)/%,$(OBJS4))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
#$(TARGET5): $(patsubst %,$(LIBDIR)/%,$(OBJS5))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
#$(TARGET6): $(patsubst %,$(LIBDIR)/%,$(OBJS6))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
#$(TARGET7): $(patsubst %,$(LIBDIR)/%,$(OBJS7))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
#$(TARGET8): $(patsubst %,$(LIBDIR)/%,$(OBJS8))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
#$(TARGET9): $(patsubst %,$(LIBDIR)/%,$(OBJS9))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
#$(TARGET10): $(patsubst %,$(LIBDIR)/%,$(OBJS10))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
#$(TARGET11): $(patsubst %,$(LIBDIR)/%,$(OBJS11))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
#
#$(TARGET12): $(patsubst %,$(LIBDIR)/%,$(OBJS12))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET13): $(patsubst %,$(LIBDIR)/%,$(OBJS13))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET14): $(patsubst %,$(LIBDIR)/%,$(OBJS14))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET15): $(patsubst %,$(LIBDIR)/%,$(OBJS15))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET16): $(patsubst %,$(LIBDIR)/%,$(OBJS16))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET17): $(patsubst %,$(LIBDIR)/%,$(OBJS17))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET18): $(patsubst %,$(LIBDIR)/%,$(OBJS18))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET19): $(patsubst %,$(LIBDIR)/%,$(OBJS19))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET20): $(patsubst %,$(LIBDIR)/%,$(OBJS20))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET21): $(patsubst %,$(LIBDIR)/%,$(OBJS21))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o $(LIBDIR)/*.d core $(BINDIR)/*

