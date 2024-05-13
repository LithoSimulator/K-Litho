# Compiler
CXX = g++

# Compile options
CXXFLAGS = -O2 -Wall -I../common

# Libraries
LIBS = -lfftw3

# Target
TARGET = klitho_socs

# Source files
SRCS = klitho_socs.cpp ../common/mask.cpp ../common/file_io.cpp

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
