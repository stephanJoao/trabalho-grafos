all: msg program
	@echo Rodando o programa
	@./program

msg:
	@echo Vou tentar compilar o programa.

main.o: main.cpp include/Graph.hpp
	@g++ -c main.cpp 

graph.o: src/Graph.cpp include/Graph.hpp
	@g++ -c src/Graph.cpp

vertex.o: src/Vertex.cpp include/Vertex.hpp
	@g++ -c src/Vertex.cpp
	
edge.o: src/Edge.cpp include/Edge.hpp
	@g++ -c src/Edge.cpp

program: main.o edge.o vertex.o graph.o
	@g++ -o program main.o Edge.o Vertex.o Graph.o

limpa:
	@rm program *.o