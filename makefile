CC= g++
INCLDIR= ./include/
SRCDIR= ./src/
OBJDIR= ./objs/
SRIMDIR = ./srim/
ROOT= `root-config --cflags --glibs`
CFLAGS= -g -Wall $(ROOT)
CPPFLAGS= -I$(INCLDIR) -I./ 
LDFLAGS= -L$(INCLDIR) $(ROOT)
SRC= $(wildcard $(SRCDIR)*.cpp)
OBJS= $(SRC:$(SRCDIR)%.cpp=$(OBJDIR)%.o)
DICT= track_dict.cxx
LIB_BOOKS= Track.h SiHit.h PCHit.h CsIHit.h LinkDef_Track.h 
EXE= analyzer
SRIM_EXE= srim_clean


.PHONY: all clean

all: $(EXE) $(SRIM_EXE)

$(EXE): $(DICT) $(OBJS)
	$(CC) $(LDFLAGS) $^ -o $@


$(OBJDIR)%.o: $(SRCDIR)%.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@ 

$(DICT): $(LIB_BOOKS)
	rootcint -f $@ -c $^

$(SRIM_EXE): $(SRIMDIR)SRIMscript.cpp
	$(CC) -o $@ $^

clean:
	$(RM) $(OBJS) $(DICT) $(EXE) *.pcm
