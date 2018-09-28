### make changes accordingly ###
CC       = gcc
CPP      = g++
CLINKER  = gcc
CCLINKER = g++
MAKE     = make --no-print-directory
SHELL    = /bin/sh
CFLAGS		= -Wall -pedantic 
OPTI            = -O3
#OPTI = -pg -g # for profiling
LDFLAGS	= 
INCLUDES	= 
LIBS	= -lm -lgsl -lgslcblas
DEFINES = -DVERBOSE 

default: choleramodel

choleramodel: Makefile GridCells.o Community.o Population.o R0Community.o driver.o
	$(CCLINKER) -o choleramodel driver.o GridCells.o Community.o Population.o R0Community.o $(LDFLAGS) $(LIBS)

grid: Makefile GridCells.o EpiModel.o sir.o
	$(CCLINKER) -o grid sir.o GridCells.o EpiModel.o $(LDFLAGS) $(LIBS)

GridCells.o: GridCells.cpp GridCells.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c GridCells.cpp

EpiModel.o: EpiModel.cpp EpiModel.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c EpiModel.cpp

Community.o: Community.cpp Community.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c Community.cpp

Population.o: Population.cpp Population.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c Population.cpp

R0Community.o: R0Community.cpp R0Community.h Community.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c R0Community.cpp

CholeraModel.o: CholeraModel.cpp CholeraModel.h EpiModel.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c CholeraModel.cpp

%.o: %.cpp GridCells.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c $<

emacs:
	emacs Makefile *.h *.cpp &

clean:
	rm -f *.o choleramodel

zip:
	cd ..; zip choleracode.zip code/Makefile code/README code/gpl.txt code/*.cpp code/*.h

