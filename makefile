CC= g++
INCLDIR= ./include/
SRCDIR= ./src/
OBJDIR= ./objs/
SRIMDIR = ./srim/
BADDIR = ./badDet/
ROOT= `root-config --cflags --glibs`
CFLAGS= -g -Wall $(ROOT)
CPPFLAGS= -I$(INCLDIR) -I./ 
LDFLAGS= -L$(INCLDIR) $(ROOT)
SRC= $(wildcard $(SRCDIR)*.cpp)
OBJS= $(SRC:$(SRCDIR)%.cpp=$(OBJDIR)%.o) 
DICT= $(SRCDIR)track_dict.cxx
LIB_BOOKS= $(INCLDIR)Track.h $(INCLDIR)SiHit.h $(INCLDIR)PCHit.h $(INCLDIR)CsIHit.h $(INCLDIR)LinkDef_Track.h 
LIB= $(OBJDIR)track_dict.o
EXE= analyzer
SRIM_EXE= srim_clean
BAD_EXE= badList


.PHONY: all clean

all: $(EXE) $(SRIM_EXE) $(BAD_EXE)

$(EXE): $(LIB) $(OBJS)
	$(CC) $(LDFLAGS) $^ -o $@

$(OBJDIR)%.o: $(SRCDIR)%.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@ 

$(LIB): $(DICT)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $^
	mv $(SRCDIR)*.pcm ./

$(DICT): $(LIB_BOOKS)
	rootcling -f $@ -c $^

$(SRIM_EXE): $(SRIMDIR)SRIMscript.cpp
	$(CC) -o $@ $^

$(BAD_EXE): $(DICT) $(BADDIR)BadDetector.cpp
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^

clean:
	$(RM) $(OBJS) $(DICT) $(EXE) *.pcm
