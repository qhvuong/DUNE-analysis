
CXX = g++
ROOTFLAGS = `root-config --libs --cflags`
INCLUDE = -I/cvmfs/dune.opensciencegrid.org/products/dune/dune_plot_style/v01_00/NULL/include

# ROOT Flags are incomplete.
#LDLIBS += -L$(ROOTSYS)/lib
#INCLUDE += -I$(ROOTSYS)/include

all: DoTemplateFit

DoTemplateFit : DoTemplateFit.cxx
	$(CXX) -fopenmp DoTemplateFit.cxx -o DoTemplateFit $(ROOTFLAGS) $(INCLUDE) 

clean:
	-rm -f DoTemplateFit
	-rm -f DoTemplateFit.o

remake:
	make clean
	make all
