# Various flags

CXX  = g++
LINK = $(CXX)
# -O3 -funroll-loops -mmmx
#CXXFLAGS = -Wall -g 
CXXFLAGS = -std=c++11 -Wall -O3 -funroll-loops -pipe
LFLAGS = -lm


TARGET  = sigclu

HEADER  = significanceclustering.h
FILES = significanceclustering.cc

OBJECTS = $(FILES:.cc=.o)

$(TARGET): ${OBJECTS}
	$(LINK) $^ $(LFLAGS) -o $@

all: $(TARGET)

clean:
	rm -f $(OBJECTS)

distclean:
	rm -f $(OBJECTS) $(TARGET)

# Compile and dependency
$(OBJECTS): $(HEADER) Makefile




