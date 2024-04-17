# Compiler
CXX = g++

# Compile options
CXXFLAGS = -O2 -Wall

# Libraries
LIBS = -lfftw3 -llapack -lblas

# Target
TARGET = klitho_tcc

# Source files
SRCS = klitho_tcc.cpp file_io.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default target
.DEFAULT_GOAL := all

# Rules
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) $(TARGET)

.PHONY: all clean
