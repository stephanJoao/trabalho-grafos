#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

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

    graph->insertMissingVertices();

    return graph;
}

void printOptions()
{
    cout << "MENU" << endl;
    cout << "----" << endl;
    cout << "[1] Subgrafo induzido por conjunto de vértices" << endl;
    cout << "[2] Subgrafo induzido por conjunto de vértices" << endl;
    cout << "[3] Caminho Mínimo entre dois vértices - Dijkstra" << endl;
    cout << "[4] Caminho Mínimo entre dois vértices - Floyd" << endl;
    cout << "[5] Árvore Geradora Mínima de Kruskal" << endl;
    cout << "[6] Imprimir caminhamento em largura" << endl;
    cout << "[7] Imprimir ordenacao topológica" << endl;
    cout << "[0] Sair" << endl;
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

    // Read of input_file
    
    ifstream input_file;
    input_file.open(input_file_name, ios::in);

    Graph *g;
    if(input_file.is_open())
    {
        g = readGraph(input_file, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
        g->setOutfileName(output_file_name);
    }
    
    // Menu

    int choice, a, b;
    mkdir("output", 0777);
    do 
    {
        cout << "Options..." << endl;
        printOptions();
        cin >> choice;

        switch (choice)
        {
            case 0:
                /* code */
                break;
            case 1:
                /* A */
                break;
            case 2:
                /* B */
                break;
            case 3:
                cout << "Type the vertices for the algorithm: ";
                cin >> a >> b;
                g->Dijkstra(a, b); 
                break;
            case 4:
                cout << "Type the vertices for the algorithm: ";
                cin >> a >> b;
                g->Floyd(a, b); 
                break;
            case 5:
                g->MST_Kruskal();
                break;
            case 6:
                cout << "Type the vertex for the algorithm: ";
                cin >> a;
                g->BFS(a);
                break;
            case 7:
                g->topologicalSorting();
                break;
            default:
                break;
        }
    }
    while(choice != 0);
    
    delete g;

    return 0;
}