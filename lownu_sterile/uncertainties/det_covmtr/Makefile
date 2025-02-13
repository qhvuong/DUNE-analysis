
CXX = g++
ROOTFLAGS = `root-config --libs --cflags`

# ROOT Flags are incomplete.
#LDLIBS += -L$(ROOTSYS)/lib
#INCLUDE += -I$(ROOTSYS)/include

all: lownu

lownu : lownu.cxx
	$(CXX) -o lownu lownu.cxx $(ROOTFLAGS) 

clean:
	-rm -f lownu
	-rm -f lownu.o

remake:
	make clean
	make all
