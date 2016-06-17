

#
# Makefile for Root-based Applications
# Environment variable ROOTSYS must be set
#

# Compiler/Linker and Flags
CC = g++
CFLAGS = -D_REENTRANT -g -fPIC -w -I$(shell root-config --incdir)
LD = g++
LDFLAGS = -g -w

NAME = data2012_2DTagging_MC
ROOTLIBS = $(shell root-config --libs) 

all: $(NAME)

%.o: %.cpp
	@echo "compiling the classes..."
	$(CC) $(CFLAGS) -c -o $@ $<

$(NAME):$(NAME).o ptBinning.o alphaBinning.o QGSyst.o 
	@echo "linking..."
	$(LD) $(LDFLAGS) $(NAME).o -o $(NAME) ptBinning.o alphaBinning.o QGSyst.o $(ROOTLIBS) 

clean:
	@echo "cleaning up..."
	rm $(NAME).o ptBinning.o alphaBinning.o QGSyst.o *.d
