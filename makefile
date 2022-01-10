CC = mpicxx
 
all: nbody worker

nbody: main.cpp
	mpicxx -o nbody main.cpp -lglfw -framework Cocoa -framework OpenGL -framework IOKit -std=c++17

worker: worker.cpp
	mpicxx -o worker worker.cpp
