CC	= g++
CPPFLAGS= -Wno-deprecated -std=c++11 -O3 -m64 -c -w #-Wall
LDFLAGS	= -O3 -m64 
SOURCES	= Run.cc
OBJECTS	= $(SOURCES:.cc=.o)
EXECUTABLE=Run

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(LDFLAGS) $@.o -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o
