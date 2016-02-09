C=g++
LD=g++
F77=g77
#O2 for optimization, g for debugging
SPECIALFLAGS=-g # -O2 -g
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs) -lMinuit

CFLAGS = $(SPECIALFLAGS) -I. -I./include -I$(SRT_PUBLIC_CONTEXT)/include -I$(ROOTSYS)/include -Wno-deprecated
LFLAGS = $(SPECIALFLAGS) -lz 


RCXX=$(CFLAGS) $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS)

SRC=SusyScan.cc ConfigFile.cc GeneratorMasses.cc simpleGenMasses.cc PlotTools.cc TheLimits.cc

%.o: %.cc
		$(C) $(RCXX) -c $<

all: plot_RA2 plot_SMS plot_mSUGRA smooth_SMS plot_SMS_fromCLsROOT plot_T2xx


ConfigFile.o: ConfigFile.cc ConfigFile.h
		$(C) $(RCXX) -c ConfigFile.cc 

SusyScan.o: SusyScan.cc SusyScan.h
		$(C) $(RCXX) -c SusyScan.cc 

GeneratorMasses.o: GeneratorMasses.cc GeneratorMasses.h
		$(C) $(RCXX) -c GeneratorMasses.cc 

simpleGenMasses.o: simpleGenMasses.cc simpleGenMasses.h
		$(C) $(RCXX) -c simpleGenMasses.cc 

PlotTools.o: PlotTools.cc PlotTools.h GlobalFunctions.h
		$(C) $(RCXX) -c PlotTools.cc 

TheLimits.o: TheLimits.cc TheLimits.h
		$(C) $(RCXX) -c TheLimits.cc 

plot_RA2.o: plot_RA2.cc plot.h GlobalFunctions.h simpleGenMasses.h StyleSettings.h
		$(C) $(RCXX) -c plot_RA2.cc 

plot_SMS.o: plot_SMS.cc plot.h GlobalFunctions.h simpleGenMasses.h StyleSettings.h
		$(C) $(RCXX) -c plot_SMS.cc 

plot_T2xx.o: plot_T2xx.cc plot.h GlobalFunctions.h simpleGenMasses.h StyleSettings.h
		$(C) $(RCXX) -c plot_T2xx.cc 

plot_RA2: $(SRC:.cc=.o) plot_RA2.o 
		$(LD) $(SRC:.cc=.o) plot_RA2.o $(RLXX) $(JCORR) -o plot_RA2
		@echo '-> plot_RA2 executable created.'

plot_SMS: $(SRC:.cc=.o) plot_SMS.o 
		$(LD) $(SRC:.cc=.o) plot_SMS.o $(RLXX) $(JCORR) -o plot_SMS
		@echo '-> plot_SMS executable created.'

plot_T2xx: $(SRC:.cc=.o) plot_T2xx.o 
		$(LD) $(SRC:.cc=.o) plot_T2xx.o $(RLXX) $(JCORR) -o plot_T2xx
		@echo '-> plot_T2xx executable created.'

plot_mSUGRA.o: plot_mSUGRA.cc plot.h GlobalFunctions.h simpleGenMasses.h StyleSettings.h
		$(C) $(RCXX) -c plot_mSUGRA.cc 

plot_mSUGRA: $(SRC:.cc=.o) plot_mSUGRA.o 
		$(LD) $(SRC:.cc=.o) plot_mSUGRA.o $(RLXX) $(JCORR) -o plot_mSUGRA
		@echo '-> plot_mSUGRA executable created.'

smooth_SMS.o: smooth_SMS.cc plot.h GlobalFunctions.h simpleGenMasses.h StyleSettings.h
		$(C) $(RCXX) -c smooth_SMS.cc

smooth_SMS: $(SRC:.cc=.o) smooth_SMS.o
		$(LD) $(SRC:.cc=.o) smooth_SMS.o $(RLXX) $(JCORR) -o smooth_SMS
		@echo '-> smooth_SMS executable created.'

plot_SMS_fromCLsROOT.o: plot_SMS_fromCLsROOT.cc plot.h GlobalFunctions.h simpleGenMasses.h StyleSettings.h
		$(C) $(RCXX) -c plot_SMS_fromCLsROOT.cc

plot_SMS_fromCLsROOT: $(SRC:.cc=.o) plot_SMS_fromCLsROOT.o
		$(LD) $(SRC:.cc=.o) plot_SMS_fromCLsROOT.o $(RLXX) $(JCORR) -o plot_SMS_fromCLsROOT
		@echo '-> plot_SMS_fromCLsROOT executable created.'

clean:
		@rm -f *.o 
		@rm -f plot plot_SMS plot_mSUGRA smooth_SMS plot_SMS_fromCLsROOT
