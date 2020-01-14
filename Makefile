# GNUmakefile for libWCSimBonsai 2015/08/12 T.Yano
# --------------------------------------------------------------


ROOTCFLAGS   := $(shell root-config --cflags) -DUSE_ROOT -fPIC
ROOTLIBS     := $(shell root-config --libs)
ifndef WCSIMDIR
$(error Environment variable WCSIMDIR is not set)
endif
WCSIMINCS     = -I$(WCSIMDIR)/include
WCSIMLIBS     = -L$(WCSIMDIR) -lWCSimRoot

CPPFLAGS  += -Wno-deprecated 
CPPFLAGS  += $(ROOTCFLAGS) $(WCSIMINCS)
EXTRALIBS += $(ROOTLIBS) $(WCSIMLIBS)
#CXXFLAGS  += ""

CXX = g++
#CXX = $(shell which g++)

WORKDIR := .
TMPDIR := $(WORKDIR)/tmp


EBONSAISO    := libWCSimEBonsai.so

EBONSAISRC := ./WCSimEBonsai.cpp \
			 ./WCSimEBonsai.h

EBONSAIOBJS	:= $(TMPDIR)/WCSimEBonsai.o $(TMPDIR)/WCSimEBonsaiDict.o


.PHONY: directories

all: directories libWCSimEBonsai.so

directories: $(TMPDIR)

$(TMPDIR) :
	mkdir -p $(TMPDIR)

libWCSimEBonsai.so : $(EBONSAIOBJS) 
	@if [ ! -d $(TMPDIR) ] ; then mkdir $(TMPDIR) ; echo mkdir $(TMPDIR) ;fi
	$(CXX) -shared -O $^ -o $(EBONSAISO) $(EBONSAILIBS) $(EXTRALIBS)

WCSimEBonsaiDict.cpp : $(EBONSAISRC)
	@echo Compiling rootcint ...
	rootcint  -f WCSimEBonsaiDict.cpp -c -I$(shell root-config --incdir) $(WCSIMINCS) \
		WCSimEBonsai.h

rootcint: WCSimEBonsaiDict.cpp

$(TMPDIR)/%.o : %.cpp
	@echo Compiling $*.cpp ...
	@if [ ! -d $(TMPDIR) ] ; then mkdir $(TMPDIR) ; echo mkdir $(TMPDIR) ;fi
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $(TMPDIR)/$(*F).o $<

clean :
	@rm -f $(TMPDIR)/*.o
	@rm -f libWCSimEBonsai.so
	@rm -f WCSimEBonsaiDict.*
