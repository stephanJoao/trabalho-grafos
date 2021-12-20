#include <iostream>
#include <fstream>
#include <unordered_map>

#include "include/Graph.hpp"

using namespace std;

Graph* readGraph(ifstream& input_file, bool directed, bool weighted_edge, bool weighted_vertex)
{    
    //Variáveis para auxiliar na criação dos nós no Grafo
    int idNodeSource;
    int idNodeTarget;
    int order;
    input_file >> order;

    //Criando objeto grafo
    Graph *graph = new Graph(order, directed, weighted_edge, weighted_vertex);
    //Leitura de arquivo

    if(!graph->isWeightedEdge() && !graph->isWeightedVertex())
    {
        while(input_file >> idNodeSource >> idNodeTarget)
        {
            graph->insertEdge(idNodeSource, idNodeTarget);
        }
    } 
    else if(graph->isWeightedEdge() && !graph->isWeightedVertex() )
    {
        float edgeWeight;

        while(input_file >> idNodeSource >> idNodeTarget >> edgeWeight)
        {
            graph->insertEdge(idNodeSource, idNodeTarget, edgeWeight);
        }
    } 
    else if(graph->isWeightedVertex() && !graph->isWeightedEdge())
    {
        float nodeSourceWeight, nodeTargetWeight;

        while(input_file >> idNodeSource >> nodeSourceWeight >> idNodeTarget >> nodeTargetWeight) 
        {
            graph->insertEdge(idNodeSource, idNodeTarget, 1, nodeSourceWeight, nodeTargetWeight);
        }
    } 
    else 
    {
        float nodeSourceWeight, nodeTargetWeight, edgeWeight;

        while(input_file >> idNodeSource >> nodeSourceWeight >> idNodeTarget >> nodeTargetWeight) 
        {
            graph->insertEdge(idNodeSource, idNodeTarget, edgeWeight, nodeSourceWeight, nodeTargetWeight);
        }
    }

    return graph;
}

int main(int argc, char const *argv[]) 
{
    // Verifies if all parameters have been provided
    if (argc != 6) 
    {
        cout << "ERROR: Expecting: ./<program_name> <input_file> <output_file> <directed> <weighted_edge> <weighted_vertex> " << endl;
        return 1;
    }

    string program_name(argv[0]);
    string input_file_name(argv[1]);
    std::cout << "Input file: " << input_file_name << std::endl;
    string output_file_name(argv[2]);
    std::cout << "Output file: " << output_file_name << std::endl;
    bool directed = atoi(argv[3]);
    std::cout << "Directed: " << directed << std::endl;
    bool weighted_edge = atoi(argv[4]);
    std::cout << "Weighted edge: " << weighted_edge << std::endl;
    bool weighted_vertex = atoi(argv[5]);
    std::cout << "Weighted vertex: " << weighted_vertex << std::endl;


    // // Read of input_file
    ifstream input_file;
    input_file.open(input_file_name, ios::in);

    cout << "Hello World!!!" << endl;
    Graph *g = new Graph(4);
    g->insertEdge(1, 2);
    g->insertEdge(1, 3);
    g->insertEdge(1, 4);
    g->insertEdge(2, 3);
    
    g = new Graph(8, false, true, false);
    // g->insertEdge(1, 2, 5);
    // g->insertEdge(1, 3, 8);
    // g->insertEdge(1, 6, 7);
    // g->insertEdge(1, 8, 6);
    // g->insertEdge(2, 3, 4);
    // g->insertEdge(2, 4, 8);
    // g->insertEdge(2, 5, 7);
    // g->insertEdge(2, 8, 9);
    // g->insertEdge(3, 5, 4);
    // g->insertEdge(3, 6, 5);
    // g->insertEdge(4, 5, 9);
    // g->insertEdge(4, 7, 4);
    // g->insertEdge(5, 6, 3);
    // g->insertEdge(5, 7, 10);
    // g->insertEdge(6, 7, 6);
    // Graph *g;
    // if(input_file.is_open()){

    //     g = readGraph(input_file, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    // } else {
    //     cout << "Unable to open " << input_file_name << std::endl;
    //     exit(1);
    // }
    
    // delete g;
    Graph* g = new Graph(7, true, true, false);
    g->insertEdge(0, 1, 4);
    g->insertEdge(0, 2, 2);
    g->insertEdge(1, 2, -3);
    g->insertEdge(1, 3, 1);
    g->insertEdge(1, 4, 3);
    g->insertEdge(2, 4, 2);
    g->insertEdge(2, 3, 3);
    g->insertEdge(3, 4, -2);
    g->insertEdge(3, 6, 4);
    g->insertEdge(4, 5, 3);
    g->insertEdge(4, 6, 3);
    g->insertEdge(5, 6, 1);

    // g->getInfo();

    // g->printAdjList();
    // g->saveToDot("graph1.dot");
    // g->Dijkstra(0, 1);

    std::set<std::pair<int, int>> *mst_edges = g->MST_Kruskal();
    g->saveToDot("mst_kruskal.dot", mst_edges);

    std::set<std::pair<int, int>> *back_edges = new std::set<std::pair<int, int>>;
    std::set<std::pair<int, int>> *tree_edges = g->BFS(1, back_edges);
    g->saveToDot("bfs_tree.dot", tree_edges, back_edges);

    g->topologicalSorting();

    delete g;

    input_file.close();
    // } else {
    //     cout << "Unable to open " << input_file_name << std::endl;
    //     exit(1);
    // }
    
    // delete g;
    Graph* g = new Graph(7, true, true, false);
    g->insertEdge(0, 1, 4);
    g->insertEdge(0, 2, 2);
    g->insertEdge(1, 2, -3);
    g->insertEdge(1, 3, 1);
    g->insertEdge(1, 4, 3);
    g->insertEdge(2, 4, 2);
    g->insertEdge(2, 3, 3);
    g->insertEdge(3, 4, -2);
    g->insertEdge(3, 6, 4);
    g->insertEdge(4, 5, 3);
    g->insertEdge(4, 6, 3);
    g->insertEdge(5, 6, 1);

    // g->getInfo();

    // g->printAdjList();
    // g->BFS(1);
    // g->saveToDot("graph1.dot");
    // g->saveToDot("graph1.dot");
    
    g->Dijkstra(1, 2);
    // g->Floyd(1, 2);
    delete g;

    return 0;
}
