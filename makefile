CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -fopenmp
INCLUDES = -Isrc/include -Iheader
INCLUDE = $(INCLUDES)

SRC = $(wildcard src/*.cpp)
all:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o main main.cpp $(SRC)

%:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $@.cpp $(SRC)
	.\$@

clean:
	del /f *.exe main
