all: program clear
	@echo Compilando

main.o: main.cpp include/Graph.hpp
	@g++ -c main.cpp -g

graph.o: src/Graph.cpp include/Graph.hpp
	@g++ -c src/Graph.cpp -g

vertex.o: src/Vertex.cpp include/Vertex.hpp
	@g++ -c src/Vertex.cpp -g
	
edge.o: src/Edge.cpp include/Edge.hpp
	@g++ -c src/Edge.cpp -g

program: main.o edge.o vertex.o graph.o
	@g++ -o execGrupo5 main.o Edge.o Vertex.o Graph.o -g

clear:
	@rm *.o