CC = g++
CFLAGS = -std=c++11 -O3 -march=native -mtune=native -fno-stack-protector -fomit-frame-pointer

all: CatForce.cpp LifeAPI.h
	$(CC) $(CFLAGS) -fopenmp -o CatForce CatForce.cpp
