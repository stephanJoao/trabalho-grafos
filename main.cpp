#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>
#include <chrono>

#include "include/Graph.hpp"

using namespace std;

Graph* readGraph(ifstream &input_file, bool directed, bool weighted_edge, bool weighted_vertex, int* clusters)
{
    Graph* graph;

    int order;

    int node_id, node_weight;
    int source, target;        
    
    string line;
    
    // Skips lines
    for(int i = 0; i < 3; i++)
        getline(input_file, line);
    for(int i = 0; i < 4; i++)
        getline(input_file, line, ' ');

    // Reads number of clusters
    *clusters = stoi(line);    
    cout << "Clusters: " << clusters << endl;

    // Skips lines
    for(int i = 0; i < 2; i++)
        getline(input_file, line);
    for(int i = 0; i < 2; i++)
        getline(input_file, line, ' ');
    
    // Reads order of the graph
    order = stoi(line);
    cout << "Order: " << order << endl;

    // Creates graph object
    graph = new Graph(order, directed, weighted_edge, weighted_vertex);

    // Reads nodes and node weight
    for(int i = 0; i < 6; i++)
        getline(input_file, line);
    for(int i = 0; i < order; i++) 
    {
        input_file >> node_id >> node_weight;
        // printf("Node:%4d", node_id);
        // printf(" Weight:%4d\n", node_weight);
        graph->insertVertex(node_id, node_weight);
    }   

    int cont = 0;

    // Reads edges
    for(int i = 0; i < 5; i++)
        getline(input_file, line);
    for(int i = 0; i < order; i++) 
    {
        int start, end, pos = 0;
        getline(input_file, line);
                
        // cout << "Line size: " << line.size() << endl;
        // cout << "Last char: " << line[line.size() - 1] << endl;
        while(pos < line.size()) 
        {
            start = line.find('(', pos);
            if(start == string::npos)
                break;
            end = line.find(',', start);
            
            // Gets first node of the edge
            source = stoi(line.substr(start + 1, end));
            // cout << "Source: " << source;
            
            // Fix positions
            start = end;
            end = line.find(')', start);
            
            // Gets second node of the edge
            target = stoi(line.substr(start + 1, end));
            // cout << " Target: " << target << endl;
            
            graph->insertEdge(source, target);
            cont++;

            // Fix position
            pos = end + 1;
            // cout << "Pos: " << pos << endl;
        }      
        //cout << line << endl;
    }
    // cout << "Edges: " << cont << endl;

    return graph;
}

void printOptions()
{
    cout << "MENU" << endl;
    cout << "----" << endl;
    cout << "[1] Vertex induced subgraph (Direct transitive closure)" << endl;
    cout << "[2] Vertex induced subgraph (Indirect transitive closure)" << endl;
    cout << "[3] Shortest path between two vertices (Dijkstra)" << endl;
    cout << "[4] Shortest path between two vertices (Floyd)" << endl;
    cout << "[5] Minimum spanning tree of a graph (Kruskal)" << endl;
    cout << "[6] Breadth-first search of a graph" << endl;
    cout << "[7] Topological sorting " << endl;
    cout << "[0] Leave menu" << endl;
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

    int* clusters = new int;
    Graph *g;
    if (input_file.is_open())
    {
        g = readGraph(input_file, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), clusters);
        g->setOutfileName(output_file_name);
    }
    else
    {
        std::cerr << "Input file not found!" << std::endl;
        exit(1);
    }

    std::clock_t c_start = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();

    g->Greedy(*clusters, 0);
    
    std::clock_t c_end = std::clock();
    auto t_end = std::chrono::high_resolution_clock::now();
    
    long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    
    std::cout << "CPU time: " << time_elapsed_ms << " ms\n";
    std::cout << "Wall clock time: " << std::chrono::duration<double, std::milli> (t_end - t_start).count() << " ms\n";
    
    delete clusters;

    // Menu

    /*
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
            // code 
            break;
        case 1:
            cout << "Insert the starting vertex: ";
            cin >> a;
            g->DirectTransitiveClosure(a);
            break;
        case 2:
            // B 
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
    } while (choice != 0);
    */

    delete g;

    return 0;
}