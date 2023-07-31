#TARG = common
CODE = mesh_$(TARG).a

# Set-up information ==========================================================

HOME = ..
include $(HOME)/Config/Make.cf.$(TARG)

# Files =======================================================================

SRC = baseObjects.cpp hexaMesh.cpp
TPL = 
OBJ = baseObjects.o hexaMesh.o
INC = baseObjects.h hexaMesh.h
C++C = g++

# Main Rules ==================================================================

.SUFFIXES:

.SUFFIXES: .o .cpp

.cpp.o: 
	$(C++C) $(C++FLAGS) -c $<

all: $(CODE)

$(CODE): $(INC) $(SRC) $(TPL)
	@echo 'LIBRARY: Mesh'
	@make $(OBJ)
	$(AR) ruv $(CODE) $(OBJ)
	@rm -f $(OBJ)

# Miscellaneous Rules =========================================================

clean: 
	rm -f *.o *~ $(CODE) testMesh testMesh.m testMesh.txt

allclean: clean
	rm -f *.a

