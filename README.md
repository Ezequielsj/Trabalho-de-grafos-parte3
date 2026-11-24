# Biblioteca de Grafos em C — Trabalho Prático 3 (Bellman-Ford & Floyd-Warshall)

Biblioteca de manipulação e análise de grafos ponderados em C, atualizada e otimizada para o **Trabalho Prático 3** com os algoritmos de **Bellman-Ford** (Aula 14) e **Floyd-Warshall** (Aula 13).

## Otimizações Implementadas

1. **Algoritmo de Bellman-Ford (Aula 14):**
   - **Otimização de memória:** utilização de um **vetor unidimensional de distâncias** ($\Theta(V)$) em vez de uma matriz bidimensional ($\Theta(V^2)$).
   - **Otimização de tempo (propagação ativa e parada antecipada):** atualiza apenas os vértices que tiveram algum vizinho atualizado na rodada anterior (via fila de vértices ativos com heurística *Small Label First* — SLF), **parando imediatamente quando nenhuma distância muda** (fila vazia), sem executar obrigatoriamente todas as $V - 1$ iterações.
   - **Detecção exata e imediata de ciclos negativos:** rastreia o número de arestas do caminho simples (`arestas_caminho[v] = arestas_caminho[u] + 1`); como qualquer caminho simples tem no máximo $V - 1$ arestas, `arestas_caminho[v] >= V` interrompe a execução imediatamente sem falsos positivos.
   - **Grafo transposto em CSR (*Compressed Sparse Row*):** resolve o problema *todos-para-um* (destino `100`) em uma única execução sem duplicar listas encadeadas na memória.
   - **Representação CSR direta e reutilização de workspace** na execução para todos os pares (`bellman_ford_todos_pares`).

2. **Algoritmo de Floyd-Warshall (Aula 13):**
   - **Otimização de memória e vetorização:** programação dinâmica *in-place* utilizando uma matriz contígua 2D (`float *distancias_1d`), consumindo **memória $\Theta(V^2)$ em vez de $\Theta(V^3)$**, enquanto a **complexidade de tempo permanece $\Theta(V^3)$**, com laço interno livre de desvios condicionais (`INFINITY` IEEE-754) para vetorização SIMD.
   - **Árvores de caminhos mínimos:** reconstrução completa de caminhos via matriz de predecessores `pred[i * V + j]`.
   - **Parada antecipada no laço `k`:** aborta imediatamente assim que `dist[i * V + i] < 0.0f` em qualquer passo intermediário.

3. **Representação em Memória e Cronometragem:**
   - Distâncias e pesos armazenados em **`float` (4 bytes)** e vértices/predecessores em **`int` (4 bytes)**, com alocação contígua 1D para vetorização (`-O3 -march=native`) e melhor localidade de cache.
   - **Arena Memory Pool** (`BlocoArenaNos`) na lista de adjacência para evitar milhões de chamadas individuais a `malloc`/`free`.
   - Cronometragem via `clock_gettime(CLOCK_MONOTONIC)` isolando apenas o tempo de processamento (sem I/O de disco).

4. **Controle de Viabilidade no Estudo de Caso 4 (`MAX_VERTICES_CASO4 = 5000`):**
   - O Estudo de Caso 4 mantém simultaneamente 3 matrizes $V \times V$ (`dist` do Bellman-Ford, `dist` e `pred` do Floyd-Warshall; a `pred` do Bellman-Ford é liberada antes de alocar o Floyd-Warshall), ou seja, pico de **$12 \cdot V^2$ bytes** (16·V² se as 4 matrizes coexistissem), além de tempo cúbico $\Theta(V^3)$ repetido por 10 rodadas.
   - Para $V > 5000$ o programa registra **`Não executado`**, informa a memória necessária e uma **estimativa por extrapolação** (não é medição direta):
     - *Floyd-Warshall:* 1 rodada em um grafo sintético conexo de 500 vértices com grau médio equivalente, escalada por $(V/500)^3$.
     - *Bellman-Ford todos os pares:* média de 10 execuções reais (origens espaçadas) no grafo completo, multiplicada por $V$.
   - *Nota metodológica:* em grafos esparsos, o Bellman-Ford para todos os pares ($O(V \cdot (V + E))$ no caso prático) pode ser mais rápido que o Floyd-Warshall ($\Theta(V^3)$). Use sempre os tempos efetivamente medidos na análise.

5. **Estudos de Caso:** o Caso 1 usa as origens {10, 20, 30} → destino 100 (conforme o enunciado); o Caso 2 usa os pares entre {10, 20, 30, 40, 50}. Se o grafo tiver ciclo negativo (ex.: `grafo_W_1`), os algoritmos reportam isso e as distâncias não estão definidas.

## Como Compilar e Executar no WSL (Ubuntu)

```bash
make
./programa
```

No menu inicial:
* **Representação:** `1` (Lista de Adjacência)
* **Direcionado:** `1` (Sim)
* **Arquivo do grafo:** informe o nome do arquivo (ex.: `grafo_1.txt`)
* **Opção:** `19` para executar todos os Estudos de Caso (1 a 4) e salvar as tabelas em `resultados/`.
