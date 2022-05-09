CC = g++
CFLAGS = -std=c++11 -O3 -march=native -mtune=native -fno-stack-protector -fomit-frame-pointer

CatForce: CatForce.cpp LifeAPI.h
	$(CC) $(CFLAGS) -o CatForce CatForce.cpp
