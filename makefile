CC= g++
CFLAGS= -c -Wall `root-config --cflags --glibs`
LDFLAGS= `root-config --cflags --glibs`
OBJECTS= track_dict.cxx LookUp.o Reconstruct.o gwmAnalyzer.o main.o
EXECUTABLE= analyzer

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.cpp
	$(CC) -o $@ $(CFLAGS) $^

%.cxx: Track.h CsIHit.h SiHit.h PCHit.h LinkDef_Track.h
	rootcint -f $@ -c $^

.PHONY: clean
clean:
	rm *.o *.cxx *.pcm
