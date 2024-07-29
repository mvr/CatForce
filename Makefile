CXX = clang++
CXXFLAGS = -std=c++20 -Wall -Wextra -pedantic -O3 -march=native -mtune=native -flto -fno-stack-protector -fomit-frame-pointer -g
# CXXFLAGS = -std=c++20 -Wall -Wextra -pedantic -g
LDFLAGS =

CXX = /usr/local/opt/llvm/bin/clang++
# CXX = /usr/local/bin/gcc-11
# CXXFLAGS = -O3 -std=c++11 -march=native -mtune=native -fno-stack-protector -fomit-frame-pointer
# LDFLAGS = -L /usr/local/opt/gcc/lib/gcc/11 -L /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -lSystem -lstdc++

.PHONY: all instrument test

all: CatForce

CatForce: CatForce.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(INSTRUMENTFLAGS) -o CatForce CatForce.cpp $(LDFLAGS)

CatEval: CatEval.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(INSTRUMENTFLAGS) -o CatEval CatEval.cpp $(LDFLAGS)

### Instrumeting

PROFDATAEXE = /Library/Developer/CommandLineTools/usr/bin/llvm-profdata
ifneq ($(wildcard instrumenting/pass2.profdata),)
	INSTRUMENTFLAGS = -fprofile-use=instrumenting/pass2.profdata
else
	INSTRUMENTFLAGS =
endif

instrument: CatForce.cpp *.hpp
	mkdir -p instrumenting
	$(CXX) $(CXXFLAGS) -fprofile-generate=instrumenting/pass1 -o instrumenting/pass1-CatForce CatForce.cpp
	instrumenting/pass1-CatForce instrumenting/farm.in
	$(PROFDATAEXE) merge instrumenting/pass1 -o instrumenting/pass1.profdata
	$(CXX) $(CXXFLAGS) -fno-lto -fprofile-use=instrumenting/pass1.profdata -fcs-profile-generate=instrumenting/pass2 -o instrumenting/pass2-CatForce CatForce.cpp
	instrumenting/pass2-CatForce instrumenting/farm.in
	$(PROFDATAEXE) merge instrumenting/pass1.profdata instrumenting/pass2 -o instrumenting/pass2.profdata
	touch CatForce.cpp

### Testing

GTEST_CFLAGS = `pkg-config --cflags gtest_main`
GTEST_LIBS = `pkg-config --libs gtest_main`

testapp: testapp.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@ $(GTEST_LIBS)

testapp.o: tests/*.cpp *.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -c -o $@ $(GTEST_CFLAGS)

test: testapp
	./testapp
