CXX  = g++
LINK = $(CXX)
CXXFLAGS = -std=gnu++17 -Wall -O3 -Wall -Wextra
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
