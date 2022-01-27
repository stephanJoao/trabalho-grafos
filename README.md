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

Todos as funcionalidades, exceto a ordenação topológica apresentarão saídas no arquivo .dot (com nome especificado pelo comando de execução), os destaques dados no arquivo .dot às arestas de cada algoritmo serão explicados abaixo.

1. **DirectTransitiveClosure**: imprime na tela o fecho transitivo direto do vértice fornecido pelo usuário (o arquivo .dot destaca as arestas do subgrafo vértice induzido pelo conjunto com a cor vermelha)

2. **Dijkstra**: imprime na tela o menor caminho entre dois vértices utilizando o algoritmo de Dijkstra (arestas do caminho destacadas com a cor vermelha)

4. **Floyd**: imprime na tela o menor caminho entre dois vértices utilizando o algoritmo de Floyd (arestas do caminho destacadas com a cor vermelha)

5. **MST_Kruskal**: imprime na tela e no aquivo .dot os vértices da árvore geradora mínima utilizando o algoritmo de Kruskal (arestas da árvore em vermelho)

6. **BFS**: imprime em tela os vértices do caminhamento em largura a partir de um vértice dado segundo o algoritmo Breadth First Search e imprime no arquivo .dot (arestas de retorno em vermelho)

7. **topologicalSorting**: imprime na tela uma ordenação topológica de um grafo acíclico e direcionado (GAD)(sem saída .dot)