MAIN_SRC=main.f90
MAIN_FILE=main
CONST_SRC=kind_constants.f90
CONST_LIB=$(patsubst %.f90, %.o, $(CONST_SRC))
SRC_FILES=$(wildcard *.f90)
LIB_FILES=$(filter-out $(MAIN_SRC) $(CONST_SRC), $(wildcard *.f90))
OBJ_FILES=$(patsubst %.f90, %.o, *.o)

.PHONY : all
all : $(MAIN_FILE)

$(MAIN_FILE) : $(OBJ_FILES) $(CONST_LIB)
	gfortran -c -ffree-line-length-none $(MAIN_SRC)
	gfortran -o $(MAIN_SRC) $(OBJ_FILES)
	rm $(OBJ_FILES)


## constants :: compile constants module
.PHONY : constants
constants : $(CONST_LIB)

$(CONST_LIB) : $(CONST_SRC) 
	gfortran -c -ffree-line-length-none $(CONST_SRC)


## lib : compiles library files
.PHONY : lib
lib : $(OBJ_FILES) 

$(OBJ_FILES) : $(SRC_FILES) $(CONST_LIB)
	gfortran -c -ffree-line-length-none $(LIB_FILES)

## clean : remove files
.PHONY : clean 
clean : 
	rm -f *.o
	rm -f *.mod	

## variables : list variables
.PHONY : variables
variables :
	@echo MAIN_SRC : $(MAIN_SRC)
	@echo SRC_FILES : $(SRC_FILES)
	@echo LIB_FILES : $(LIB_FILES)
	@echo OBJ_FILES : $(OBJ_FILES)

.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<
