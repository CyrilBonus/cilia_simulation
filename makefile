# Makefile

CXX = g++
SRC = main.cpp wavobser_class_version.cpp reading.cpp
OUT = simulation

all: $(OUT)

$(OUT): $(SRC)
	$(CXX) $(SRC) -o $(OUT)

clean:
	rm -f $(OUT)