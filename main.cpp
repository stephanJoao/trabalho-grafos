#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>
#include <chrono>

#include "include/Graph.hpp"

using namespace std;

Graph* readGraph(ifstream &input_file, bool directed, bool weighted_edge, bool weighted_vertex)
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
    int clusters = stoi(line);    
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
    graph = new Graph(order, clusters, directed, weighted_edge, weighted_vertex);

    // Reads nodes and node weight
    for(int i = 0; i < 6; i++)
        getline(input_file, line);
    for(int i = 0; i < order; i++) 
    {
        input_file >> node_id >> node_weight;
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
                
        while(pos < line.size()) 
        {
            start = line.find('(', pos);
            if(start == string::npos)
                break;
            end = line.find(',', start);
            
            // Gets first node of the edge
            source = stoi(line.substr(start + 1, end));
            
            // Fix positions
            start = end;
            end = line.find(')', start);
            
            // Gets second node of the edge
            target = stoi(line.substr(start + 1, end));
            
            graph->insertEdge(source, target);
            cont++;

            // Fix position
            pos = end + 1;
        }      
    }

    return graph;
}

int main(int argc, char const *argv[])
{
    //* Verifies if all parameters have been provided
    if (argc != 4)
    {
        cout << "ERROR: Expecting: ./<program_name> <instance_directory_name> <instance_name> <output_file_name>" << endl;
        return 1;
    }

    //* Reads arguments of execution
    string program_name(argv[0]);
    string instance_directory_name(argv[1]);
    std::cout << "Directory name: " << instance_directory_name << std::endl;
    string instance_name(argv[2]);
    std::cout << "Instance name: " << instance_name << std::endl;
    string output_file_name(argv[3]);
    std::cout << "Output file: " << output_file_name << std::endl;

    //* Read of input_file
    ifstream input_file;
    input_file.open(instance_directory_name + "/" + instance_name, ios::in);

    Graph *g;
    if (input_file.is_open())
    {
        g = readGraph(input_file, false, false, true);
        // It's going to be greedy, greedyRA or greedyRAR
        g->setOutfileName(output_file_name);
    }
    else
    {
        std::cerr << "Input file not found!" << std::endl;
        exit(1);
    }

    g->saveToDot("teste");

    //* Experiments variables
    std::string instance_names[10] = {"n100d03p1i1.txt", "n100plap1i1.txt", 
                                      "n200d03p3i1.txt", "n200plap1i1.txt", 
                                      "n200plap1i3.txt", "n300d03p1i5.txt", 
                                      "n300plap1i1.txt", "n400plap1i5.txt", 
                                      "n500d06p3i3.txt", "n500plap3i5.txt"};
    
    //* Greedy randomized adaptative
    float alfas_gra[10]        = {0.10, 0.20, 0.30};
    const int iterations_gra   = 1000;
    const int experiments_gra  = 30;

    //* Greedy randomized adaptative reactive
    float alfas_grar[10]       = {0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50};
    const int iterations_grar  = 4000;
    const int experiments_grar = 30;


    // for(int a = 0; a < 3; a++) {
    //     for(int i = 0 ; i < 30; i++) {
    //         std::clock_t c_start = std::clock();
    //         auto t_start = std::chrono::high_resolution_clock::now();

    //         best_cost = g->GreedyRandomizedAdaptative(*clusters, alfa[a], seed, best_it, iterations);
            
    //         std::clock_t c_end = std::clock();
    //         auto t_end = std::chrono::high_resolution_clock::now();
            
    //         long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
            
    //         std::cout << "CPU time: " << time_elapsed_ms << " ms\n";
    //         std::cout << "Wall clock time: " << std::chrono::duration<double, std::milli> (t_end - t_start).count() << " ms\n";

    //         g->printGreedyRandomizedAdaptativeTxt("file_name", instance_names[2], iterations, alfa[a], *seed, best_cost, *best_it, time_elapsed_ms, std::chrono::duration<double, std::milli> (t_end - t_start).count());
    //     }
    // }

    //* Delete graph
    delete g;

    return 0;
}