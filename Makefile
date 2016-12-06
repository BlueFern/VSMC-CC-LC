
DOXY_DIR = doc

CC=g++
DOXYGEN = doxygen
CFLAGS =-c -Wall -I.

EXE = VSMCcoupledcells
DEPS = $(shell ls include/*.h)
SRCS = $(shell ls src/*.cpp) 
OBJS = $(SRCS:.cpp=.o)
DOXYF = $/doxyconfig

all: $(SRCS) $(EXE)

$(EXE): $(OBJS) 
	$(CC) $(OBJS) -o $@		
	$(DOXYGEN) $(DOXYF)

%.o: %.cpp $(DEPS)
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -fv $(OBJS)
	rm -fv $(EXE)
