# Usage:
# make        # compile all binary
# make clean  # remove all binaries and objects

#* Variables

# Name of the project
PROJ_NAME = execGrupo5

# .cpp files
CPP_SOURCE = ${wildcard ./src/*.cpp}

# .hpp files
HPP_SOURCE = ${wildcard ./include/*.hpp}

# Object files
AUX_OBJ = ${subst .cpp,.o,${subst src,objects,${CPP_SOURCE}}}

# Adds main separately
CPP_SOURCE += main.cpp
AUX_OBJ += ./objects/main.o
OBJ = ${filter-out main.o,${AUX_OBJ}}

# Compiler
CC = g++

# Compiler flags
CC_FLAGS = -c -g


#* Compilation and linking

all: objectsFolder ${PROJ_NAME}

${PROJ_NAME}: ${OBJ}
	@${CC} $^ -o $@

./objects/%.o: ./src/%.cpp ./include/%.hpp
	@${CC} $< ${CC_FLAGS} -o $@

./objects/main.o: main.cpp ${H_SOURCE}
	@${CC} $< ${CC_FLAGS} -o $@

objectsFolder:
	@mkdir -p objects

clean:
	@rm -rf ./objects/*.o

print:
	@echo ${CPP_SOURCE}
	@echo ${HPP_SOURCE}
	@echo ${OBJ}

.PHONY: all clean