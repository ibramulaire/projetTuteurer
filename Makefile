EXT = cpp
CXX=g++
CFLAGS=-Wall -I/usr/local/include
CFLAGS=
LDFLAGS= -lGL -lGLEW -lGLU -lglut  -larmadillo -lassimp -lpmp
#LDFLAGS= -lGL -lGLU -lGLEW -lglut -lm 
#LDFLAGS= -lGL -lGLU -lGLEW -lglut -lglui

#SRC=saisieinteractive.cpp courbe.cpp
SRC=$(wildcard *.$(EXT))
# $(wildcard ./utilstexture/*.$(EXT))
OBJ = $(SRC:.$(EXT)=.o)

DEBUBFLAG=-g

BIN=./
EXEC=main

all: $(EXEC)

$(EXEC): $(OBJ)
		$(CXX)  -o $(BIN)$@ $^ $(LDFLAGS)

%.o: %.cpp
		$(CXX)  -o $@ -c $< $(CFLAGS)


clean:
		rm -rf *.o
	    rm -rf $(BIN)$(EXEC)
mrproper: clean
		rm -rf $(BIN)$(EXEC)