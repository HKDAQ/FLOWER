# Based on GNUmakefile for libWCSimBonsai 2015/08/12 T.Yano
# --------------------------------------------------------------


ROOTCFLAGS   := $(shell root-config --cflags) -DUSE_ROOT -fPIC
ROOTLIBS     := $(shell root-config --libs)
ifndef WCSIM_BUILD_DIR
$(error Environment variable WCSIM_BUILD_DIR is not set)
endif
WCSIMINCS     = -I$(WCSIM_BUILD_DIR)/include/WCSim
WCSIMLIBS     = -L$(WCSIM_BUILD_DIR)/lib -lWCSimRoot

CPPFLAGS  += -Wno-deprecated 
CPPFLAGS  += $(ROOTCFLAGS) $(WCSIMINCS)
EXTRALIBS += $(ROOTLIBS) $(WCSIMLIBS)
#CXXFLAGS  += -g

CXX = g++
#CXX = $(shell which g++)

WORKDIR := .
TMPDIR := $(WORKDIR)/tmp


FLOWERSO    := libWCSimFLOWER.so

FLOWERSRC := ./WCSimFLOWER.cpp \
			 ./WCSimFLOWER.h

FLOWEROBJS	:= $(TMPDIR)/WCSimFLOWER.o $(TMPDIR)/WCSimFLOWERDict.o


.PHONY: directories

all: directories libWCSimFLOWER.so

directories: $(TMPDIR)

$(TMPDIR) :
	mkdir -p $(TMPDIR)

libWCSimFLOWER.so : $(FLOWEROBJS) 
	@if [ ! -d $(TMPDIR) ] ; then mkdir $(TMPDIR) ; echo mkdir $(TMPDIR) ;fi
	$(CXX) -shared -O $^ -o $(FLOWERSO) $(FLOWERLIBS) $(EXTRALIBS)

WCSimFLOWERDict.cpp : $(FLOWERSRC)
	@echo Compiling rootcint ...
	rootcint  -f WCSimFLOWERDict.cpp -c -I$(shell root-config --incdir) $(WCSIMINCS) \
		WCSimFLOWER.h

rootcint: WCSimFLOWERDict.cpp

$(TMPDIR)/%.o : %.cpp
	@echo Compiling $*.cpp ...
	@if [ ! -d $(TMPDIR) ] ; then mkdir $(TMPDIR) ; echo mkdir $(TMPDIR) ;fi
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $(TMPDIR)/$(*F).o $<

clean :
	@rm -f $(TMPDIR)/*.o
	@rm -f libWCSimFLOWER.so
	@rm -f WCSimFLOWERDict*
