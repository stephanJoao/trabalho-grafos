# Trabalho de Teoria do Grafos

Grafo implementado em lista de adjacência com auxílio de mapa, em que os IDs dos vértices são as chaves. 

O arquivo de entrada deve ser um arquivo de texto no formato especificado abaixo, onde os parênteses indicam argumentos obrigatórios e os colchetes indicam argumentos opcionais. Os parênteses e colchetes não devem ser incluídos.

```
(ordem)
(id_O) (id_D) [w]
(id_O) (id_D) [w]
        .
        .
        .
(id_O) (id_D) [w]
```
Em que 

* `ordem`: número de vértices no grafo

* `id_O`: ID do vértice origem na aresta

* `id_D`: ID do vértice destino na aresta

* `w`: peso da aresta

## Instruções de compilação
`make` no diretório raiz do projeto

## Instruções de execução
```./execGrupo5 "arquivoDeEntrada" "arquivoDeSaida" (direcionado) (pesoAresta) (pesoVertice)```

Em que 

* `arquivoDeEntrada`: nome do arquivo de entrada contendo o grafo no formato especificado

* `arquivoDeSaida`: nome do arquivo de saída que será impresso no formato .dot

* `direcionado`: 0 ou 1 (se o grafo é direcionado ou não)

* `pesoAresta`: 0 ou 1 (se as arestas têm peso)

* `pesoVertice`: 0 ou 1 (se os vértices têm peso)

## Funcionalidades
1. 

2. 

3. Dijkstra: imprime o menor caminho entre dois vértices segundo o algoritmo de Dijkstra

4. Floyd: imprime o menor caminho entre dois vértices segundo o algoritmo de Floyd

5. MST_Kruskal: imprime os vértices da árvore geradora mínima segundo o algoritmo de Kruskal

6. BFS: imprime o caminhamento em largura a partir de um vértice dado segundo o algoritmo Breadth First Search. 

7. topologicaoSorting: imprime uma ordenação topológica de um grafo acíclico e direcionado