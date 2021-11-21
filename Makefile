all: msg program
	@echo Rodando o programa
	@./program

msg:
	@echo Vou tentar compilar o programa.

main.o: main.cpp include/Vertex.hpp include/Edge.hpp
	@g++ -c main.cpp 

vertex.o: src/Vertex.cpp include/Vertex.hpp
	@g++ -c src/Vertex.cpp
	
edge.o: src/Edge.cpp include/Edge.hpp
	@g++ -c src/Edge.cpp

program: main.o vertex.o edge.o
	@g++ -o program main.o Vertex.o Edge.o

limpa:
	@rm program *.o