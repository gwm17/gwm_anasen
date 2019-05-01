CC= g++
INCLDIR= ./include/
SRCDIR= ./src/
LIBDIR= ./lib/
OBJDIR= ./objs/
ROOT= `root-config --cflags --glibs`
CFLAGS= -g -Wall $(ROOT)
CPPFLAGS= -I$(INCLDIR) -I$(LIBDIR)
LDFLAGS= -L$(LIBDIR) -L$(INCLDIR) $(ROOT)
SRC= $(wildcard $(SRCDIR)*.cpp)
OBJS= $(SRC:$(SRCDIR)%.cpp=$(OBJDIR)%.o)
DICT= track_dict.cxx
LIB_ITEMS= $(wildcard $(LIBDIR)*.h)
EXE= analyzer

.PHONY: all clean

all: $(EXE)

$(EXE): $(DICT) $(OBJS)
	$(CC) $(LDFLAGS) $^ -o $@

$(OBJDIR)%.o: $(SRCDIR)%.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@ 

$(DICT): $(LIB_ITEMS)
	rootcint -f $@ -c $^

clean:
	$(RM) $(OBJS) $(DICT) $(EXE) *.pcm
