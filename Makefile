CC = clang++
CFLAGS = -std=c++11 -Wall -Wextra -pedantic -O3 -march=native -mtune=native -flto -fno-stack-protector -fomit-frame-pointer
# CFLAGS = -std=c++11 -Wall -Wextra -pedantic -g -fno-stack-protector -fomit-frame-pointer
LDFLAGS =

# CC = /usr/local/bin/gcc-11
# CFLAGS = -O3 -std=c++11 -march=native -mtune=native -fno-stack-protector -fomit-frame-pointer
# LDFLAGS = -L /usr/local/opt/gcc/lib/gcc/11 -L /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -lSystem -lstdc++

PROFDATAEXE = /Library/Developer/CommandLineTools/usr/bin/llvm-profdata
ifeq ($(wildcard $(instrumenting/profile.profdata)),)
	INSTRUMENTFLAGS = -fprofile-use=instrumenting/profile.profdata
else
	INSTRUMENTFLAGS =
endif

all: CatForce
CatForce: CatForce.cpp LifeAPI.h
	$(CC) $(CFLAGS) $(INSTRUMENTFLAGS) -o CatForce CatForce.cpp $(LDFLAGS)

instrument: CatForce.cpp LifeAPI.h
	$(CC) $(CFLAGS) -fprofile-instr-generate=instrumenting/pass1.profraw -o instrumenting/pass1-CatForce CatForce.cpp
	instrumenting/pass1-CatForce instrumenting/farm.in
	rm instrumenting/farm.rle
	$(PROFDATAEXE) merge instrumenting/pass1.profraw -o instrumenting/profile.profdata
	touch CatForce.cpp

# $(CC) $(CFLAGS) -fprofile-generate=instrumenting -o instrumenting/pass1-CatForce CatForce.cpp
# instrumenting/pass1-CatForce instrumenting/farm.in
# rm instrumenting/farm.rle
# $(PROFDATA) merge instrumenting/default*.profraw -o instrumenting/pass1.profdata
# $(CC) $(CFLAGS) -fprofile-use=instrumenting/pass1.profdata -fcs-profile-generate=instrumenting -o instrumenting/pass2-CatForce CatForce.cpp
# instrumenting/pass2-CatForce instrumenting/farm.in
# rm instrumenting/farm.rle
# $(PROFDATAEXE) merge instrumenting/default*.profraw instrumenting/pass1.profdata -o instrumenting/profile.profdata
