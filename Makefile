all: msg program
	@echo Rodando o programa
	@./program

msg:
	@echo Vou tentar compilar o programa.

main.o: main.cpp include/Vertex.hpp
	@g++ -c main.cpp 

vertex.o: src/Vertex.cpp include/Vertex.hpp
	@g++ -c src/Vertex.cpp

program: main.o vertex.o
	@g++ -o program main.o Vertex.o

limpa:
	@rm program *.o