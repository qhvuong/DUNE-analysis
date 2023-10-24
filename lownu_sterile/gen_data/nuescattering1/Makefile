
CXX = g++
ROOTFLAGS = `root-config --libs --cflags`

# ROOT Flags are incomplete.
#LDLIBS += -L$(ROOTSYS)/lib
#INCLUDE += -I$(ROOTSYS)/include

all: nuescattering

nuescattering : nuescattering.cxx
	$(CXX) -o nuescattering nuescattering.cxx $(ROOTFLAGS) 

clean:
	-rm -f nuescattering
	-rm -f nuescattering.o

remake:
	make clean
	make all
