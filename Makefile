CCLIB=-lsdsl -ldivsufsort -ldivsufsort64 -Wno-comment -fopenmp
VLIB= -g -O0

LIB_DIR = ${HOME}/lib
INC_DIR = ${HOME}/include
MY_CXX_FLAGS= -std=c++11 -Wall -DNDEBUG -D__STDC_FORMAT_MACROS -fomit-frame-pointer -Wno-char-subscripts -Wno-sign-compare
#-D_FILE_OFFSET_BITS=64

MY_CXX_OPT_FLAGS= -O3 -m64 
#MY_CXX_OPT_FLAGS= $(VLIB)
MY_CXX=g++

##

LIBOBJ = \
	lib/file.o\
	lib/utils.o\
	external/malloc_count/malloc_count.o\
	external/gsacak.o
	
##

M64 = 0
DEBUG = 0
CHECK = 1

##

LFLAGS = -lm -lrt -ldl

DEFINES = -DDEBUG=$(DEBUG) -DM64=$(M64)

CXX_FLAGS=$(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -L$(LIB_DIR) $(LFLAGS) $(DEFINES)

CLAGS= -DSYMBOLBYTES=1


##

DIR = dataset/
INPUT = input.100.txt

K	= 1
MODE 	= 1
CHECK	= 1
OUTPUT	= 0

##

DEFINES = -DDEBUG=$(DEBUG) -DM64=$(M64) 

CFLAGS += $(DEFINES)

all: compile

clean:
	\rm -f *.o  external/*.o lib/*o external/malloc_count/*.o all-bwsd

##

lib: lib/file.c lib/utils.c external/gsacak.c external/malloc_count/malloc_count.o
	$(MY_CXX) $(CXX_FLAGS) $(DEFINES) -c lib/file.c -o lib/file.o 
	$(MY_CXX) $(CXX_FLAGS) $(DEFINES) -c lib/utils.c -o lib/utils.o 
	$(MY_CXX) $(CXX_FLAGS) $(DEFINES) -c external/gsacak.c -o external/gsacak.o 


compile: lib main.cpp ${LIBOBJ} 
	$(MY_CXX) $(CXX_FLAGS) -Wextra main.cpp $(CCLIB) -o all-bwsd ${LIBOBJ} 
#	$(CC) $(CFLAGS) $(LFLAGS) -o main main.c ${LIBOBJ}

run:
	./main $(DIR) $(INPUT) $(K) $(MODE) $(CHECK)

valgrind:
	valgrind --tool=memcheck --leak-check=full --track-origins=yes ./main $(DIR) $(INPUT) $(K) $(MODE) $(CHECK)
