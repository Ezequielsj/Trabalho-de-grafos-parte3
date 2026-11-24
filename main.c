#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "biblioteca_grafos.h"

#define PASTA_RESULTADOS "resultados"
#define NUM_EXECUCOES_BENCHMARK 10

/* Limite pratico para Cases 2 e 3 (2 matrizes V x V = 8 * V^2 bytes; V=30000 -> ~6.71 GB) */
#define MAX_VERTICES_FLOYD 5000

/* Limite pratico para Caso 4 (4 matrizes V x V = 16 * V^2 bytes + tempo cubico x 10 rodadas) */
#define MAX_VERTICES_CASO4 5000
#define SUBGRAFO_CALIBRACAO_V 500

static double obter_tempo_segundos(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

void garantir_diretorio_resultados(void) {
#ifdef _WIN32
    _mkdir(PASTA_RESULTADOS);
#else
    mkdir(PASTA_RESULTADOS, 0777);
#endif
}

const char *extrair_nome_base(const char *caminho, char *buffer, size_t tam_buffer) {
    const char *barra1 = strrchr(caminho, '/');
    const char *barra2 = strrchr(caminho, '\\');
    const char *inicio = caminho;

    if (barra1 && barra2) {
        inicio = (barra1 > barra2) ? barra1 + 1 : barra2 + 1;
    } else if (barra1) {
        inicio = barra1 + 1;
    } else if (barra2) {
        inicio = barra2 + 1;
    }

    strncpy(buffer, inicio, tam_buffer - 1);
    buffer[tam_buffer - 1] = '\0';

    char *ponto = strrchr(buffer, '.');
    if (ponto && ponto != buffer) {
        *ponto = '\0';
    }
    return buffer;
}

/* Constroi grafo sintetico conexo de calibracao (V_sub vertices e grau medio equivalente)
 * para calibrar o tempo cubico do Floyd-Warshall sem vies de desconexao de subgrafo induzido. */
static Grafo *criar_grafo_calibracao_fw(int v_sub, int grau_medio_alvo) {
    if (v_sub <= 1) return NULL;
    Grafo *cal = criar_grafo_lista_adj(v_sub);
    if (!cal) return NULL;

    int deg = grau_medio_alvo;
    if (deg < 2) deg = 2;
    if (deg >= v_sub) deg = v_sub - 1;

    for (int u = 0; u < v_sub; u++) {
        for (int step = 1; step <= deg; step++) {
            int v = (u + step) % v_sub;
            float w = (float)(((u * 7 + step * 13) % 50) + 1);
            adicionar_aresta(cal, u, v, w, 1);
        }
    }
    return cal;
}

static bool arquivo_existe(const char *caminho) {
    FILE *f = fopen(caminho, "r");
    if (f) { fclose(f); return true; }
    return false;
}

static bool resolver_caminho_grafo(const char *entrada, char *resolvido, size_t tam) {
    const char *home = getenv("HOME");
    const char *pref[4] = {"", "../", "../../", home ? home : ""};
    int n_pref = (home && home[0]) ? 4 : 3;
    size_t len = strlen(entrada);
    bool tem_txt = (len >= 4 && strcmp(entrada + len - 4, ".txt") == 0);
    for (int p = 0; p < n_pref; p++) {
        snprintf(resolvido, tam, p == 3 ? "%s/%s" : "%s%s", pref[p], entrada);
        if (arquivo_existe(resolvido)) return true;
        if (!tem_txt) {
            snprintf(resolvido, tam, p == 3 ? "%s/%s.txt" : "%s%s.txt", pref[p], entrada);
            if (arquivo_existe(resolvido)) return true;
        }
    }
    snprintf(resolvido, tam, "%s", entrada);
    return false;
}

static double estimar_tempo_fw_1x(const Grafo *grafo) {
    int V = grafo->num_vertices;
    int v_sub = (V >= SUBGRAFO_CALIBRACAO_V) ? SUBGRAFO_CALIBRACAO_V : V;
    int grau = (V > 0) ? (int)((double)grafo->num_arestas / (double)V + 0.5) : 2;
    Grafo *cal = criar_grafo_calibracao_fw(v_sub, grau);
    double t_cal = 0.0;
    if (cal) {
        float *d = alocar_matriz_distancias_1d(v_sub);
        int *p = alocar_matriz_predecessores_1d(v_sub);
        if (d && p) {
            double t0 = obter_tempo_segundos();
            floyd_warshall(cal, d, p);
            t_cal = obter_tempo_segundos() - t0;
        }
        free(d); free(p); liberar_grafo(cal);
    }
    double r = (double)V / (double)v_sub;
    return t_cal * r * r * r;
}

void salvar_memoria_grafo(Grafo *grafo, int direcionado, const char *nome_arquivo_saida) {
    FILE *arquivo = fopen(nome_arquivo_saida, "w");
    if (!arquivo) {
        printf("Erro ao abrir o arquivo de memoria.\n");
        return;
    }

    long long memoria_bytes = 0;

    if (grafo->tipo_representacao == LISTA_ADJACENCIA) {
        memoria_bytes += (long long)grafo->num_vertices * sizeof(NoListaAdj *);
        long long total_nos = direcionado ? (long long)grafo->num_arestas
                                          : (long long)grafo->num_arestas * 2;
        memoria_bytes += total_nos * sizeof(NoListaAdj);
    } else {
        memoria_bytes += (long long)grafo->num_vertices * sizeof(float *);
        memoria_bytes += (long long)grafo->num_vertices * (long long)grafo->num_vertices * sizeof(float);
    }
    memoria_bytes += (long long)grafo->num_vertices * sizeof(int);
    memoria_bytes += sizeof(Grafo);

    double memoria_kb = memoria_bytes / 1024.0;
    double memoria_mb = memoria_kb / 1024.0;

    fprintf(arquivo, "Memoria estimada para a representacao do grafo (%s):\n",
            direcionado ? "Direcionado" : "Nao-Direcionado");
    fprintf(arquivo, "Bytes: %lld B\n", memoria_bytes);
    fprintf(arquivo, "Kilobytes: %.2f KB\n", memoria_kb);
    fprintf(arquivo, "Megabytes: %.2f MB\n", memoria_mb);

    fclose(arquivo);
}

void mostrar_menu(void) {
    printf("\n=== MENU DE OPCOES ===\n");
    printf("1. Salvar informacoes basicas do grafo\n");
    printf("2. Realizar Busca em Largura (BFS)\n");
    printf("3. Realizar Busca em Profundidade (DFS)\n");
    printf("4. Calcular distancia entre dois vertices (sem pesos)\n");
    printf("5. Calcular diametro do grafo (sem pesos)\n");
    printf("6. Encontrar componentes conexas\n");
    printf("7. Calcular memoria da representacao\n");
    printf("8. Benchmark de tempo BFS (100 vertices)\n");
    printf("9. Benchmark de tempo DFS (100 vertices)\n");
    printf("10. Caminho Minimo (Dijkstra - Todos os destinos a partir de uma origem)\n");
    printf("11. Caminho Minimo (Dijkstra - Entre dois vertices especificos)\n");
    printf("12. Caminho Minimo (Dijkstra - Destinos {10, 20, 30, 40, 50} a partir de 1)\n");
    printf("13. Caminho Minimo (Dijkstra SEM HEAP - Destinos {10, 20, 30, 40, 50} a partir de 1)\n");
    printf("14. Benchmark Dijkstra (Com Heap vs Sem Heap)\n");
    printf("--- TRABALHO PRATICO 3 (BELLMAN-FORD & FLOYD-WARSHALL) ---\n");
    printf("15. [TP3 - Caso 1] Bellman-Ford: origens {10, 20, 30, 40, 50} -> destino 100\n");
    printf("16. [TP3 - Caso 2] Floyd-Warshall: pares entre {10, 20, 30, 40, 50}\n");
    printf("17. [TP3 - Caso 3] Tempo medio (10x): Bellman-Ford vs Floyd-Warshall (destino 100)\n");
    printf("18. [TP3 - Caso 4] Tempo medio (10x): Todos os pares (Bellman-Ford vs Floyd-Warshall)\n");
    printf("19. [TP3 - COMPLETO] Executar todos os Estudos de Caso (1 a 4) em sequencia\n");
    printf("0. Sair\n");
    printf("Escolha uma opcao: ");
}

void imprimir_caminho_vetor_1_indexado(const int *pais, int origem_idx, int destino_idx, int V, FILE *saida) {
    if (destino_idx == origem_idx) {
        fprintf(saida, "%d", origem_idx + 1);
        return;
    }
    if (pais[destino_idx] == -1) {
        fprintf(saida, "Inalcancavel");
        return;
    }

    int *caminho = (int *)malloc((size_t)V * sizeof(int));
    if (!caminho) {
        fprintf(saida, "Erro de memoria");
        return;
    }

    int count = 0;
    int atual = destino_idx;
    while (atual != -1 && count < V) {
        caminho[count++] = atual;
        atual = pais[atual];
    }

    for (int i = count - 1; i >= 0; i--) {
        fprintf(saida, "%d", caminho[i] + 1);
        if (i > 0) fprintf(saida, " -> ");
    }
    free(caminho);
}

void imprimir_caminho_invertido_para_destino(const int *pais_inv, int origem_idx, int destino_idx, int V, FILE *saida) {
    if (origem_idx == destino_idx) {
        fprintf(saida, "%d", origem_idx + 1);
        return;
    }
    if (pais_inv[origem_idx] == -1) {
        fprintf(saida, "Inalcancavel");
        return;
    }

    int atual = origem_idx;
    int passos = 0;
    while (atual != -1 && passos < V) {
        fprintf(saida, "%d", atual + 1);
        if (atual == destino_idx) break;
        fprintf(saida, " -> ");
        atual = pais_inv[atual];
        passos++;
    }
}

/* ============================================================================
 * ESTUDO DE CASO 1: Bellman-Ford (10, 20, 30 -> 100)
 * ============================================================================ */
void executar_tp3_caso1(Grafo *grafo, const char *nome_base) {
    garantir_diretorio_resultados();
    char caminho_saida[512];
    snprintf(caminho_saida, sizeof(caminho_saida), "%s/%s_tp3_caso1_bellman_ford.txt", PASTA_RESULTADOS, nome_base);

    FILE *arq = fopen(caminho_saida, "w");
    if (!arq) {
        printf("Erro ao criar arquivo de saida '%s'.\n", caminho_saida);
        return;
    }

    int V = grafo->num_vertices;
    int origens[] = {10, 20, 30, 40, 50};
    int num_origens = 5;
    int destino = 100;

    printf("\n=== ESTUDO DE CASO 1: BELLMAN-FORD (Destinos -> 100) ===\n");
    fprintf(arq, "=== ESTUDO DE CASO 1: BELLMAN-FORD ===\n");
    fprintf(arq, "Grafo: %s (V = %d, E = %d)\n", nome_base, V, grafo->num_arestas);
    fprintf(arq, "Origens: {10, 20, 30, 40, 50} | Destino: 100\n\n");

    if (V < 100) {
        printf("Aviso: O grafo possui apenas %d vertices (menor que o destino 100).\n", V);
        fprintf(arq, "Grafo possui apenas %d vertices; vertice 100 nao existe.\n", V);
        fclose(arq);
        return;
    }

    float *distancias = (float *)malloc((size_t)V * sizeof(float));
    int *pais = (int *)malloc((size_t)V * sizeof(int));
    if (!distancias || !pais) {
        printf("Erro ao alocar memoria para Bellman-Ford.\n");
        free(distancias);
        free(pais);
        fclose(arq);
        return;
    }

    double t_ini = obter_tempo_segundos();
    bool sem_ciclo_negativo = bellman_ford_para_destino(grafo, destino, distancias, pais);
    double tempo_exec = obter_tempo_segundos() - t_ini;

    if (!sem_ciclo_negativo) {
        printf("Resultado: CICLO NEGATIVO DETECTADO! Nao existem caminhos minimos bem definidos.\n");
        fprintf(arq, "Status: CICLO NEGATIVO DETECTADO!\n");
        fprintf(arq, "Tempo de execucao: %.6f s\n", tempo_exec);
    } else {
        printf("Status: Concluido sem ciclo negativo (tempo: %.6f s)\n\n", tempo_exec);
        fprintf(arq, "Status: Concluido sem ciclo negativo\n");
        fprintf(arq, "Tempo de execucao (todos-para-100 via grafo transposto CSR): %.6f s\n\n", tempo_exec);
        fprintf(arq, "%-10s | %-10s | %-15s | %s\n", "Origem", "Destino", "Distancia", "Caminho Minimo");
        fprintf(arq, "-----------+------------+-----------------+------------------------------------------\n");

        for (int i = 0; i < num_origens; i++) {
            int u = origens[i];
            int u_idx = u - 1;
            int dest_idx = destino - 1;

            if (distancias[u_idx] == FLT_MAX) {
                printf("Origem %2d -> Destino %d: INF (Inalcancavel)\n", u, destino);
                fprintf(arq, "%-10d | %-10d | %-15s | Inalcancavel\n", u, destino, "INF");
            } else {
                printf("Origem %2d -> Destino %d: %.2f | Caminho: ", u, destino, distancias[u_idx]);
                imprimir_caminho_invertido_para_destino(pais, u_idx, dest_idx, V, stdout);
                printf("\n");

                fprintf(arq, "%-10d | %-10d | %-15.2f | ", u, destino, distancias[u_idx]);
                imprimir_caminho_invertido_para_destino(pais, u_idx, dest_idx, V, arq);
                fprintf(arq, "\n");
            }
        }
    }

    free(distancias);
    free(pais);
    fclose(arq);
    printf("Resultados salvos em: %s\n", caminho_saida);
}

/* ============================================================================
 * ESTUDO DE CASO 2: Floyd-Warshall (pares entre {10, 20, 30, 40, 50})
 * ============================================================================ */
void executar_tp3_caso2(Grafo *grafo, const char *nome_base) {
    garantir_diretorio_resultados();
    char caminho_saida[512];
    snprintf(caminho_saida, sizeof(caminho_saida), "%s/%s_tp3_caso2_floyd_warshall.txt", PASTA_RESULTADOS, nome_base);

    FILE *arq = fopen(caminho_saida, "w");
    if (!arq) {
        printf("Erro ao criar arquivo de saida '%s'.\n", caminho_saida);
        return;
    }

    int V = grafo->num_vertices;
    int vertices_alvo[] = {10, 20, 30, 40, 50};
    int k_alvo = 5;

    printf("\n=== ESTUDO DE CASO 2: FLOYD-WARSHALL (Pares entre {10, 20, 30, 40, 50}) ===\n");
    fprintf(arq, "=== ESTUDO DE CASO 2: FLOYD-WARSHALL ===\n");
    fprintf(arq, "Grafo: %s (V = %d, E = %d)\n", nome_base, V, grafo->num_arestas);

    if (V < 50) {
        printf("Aviso: O grafo possui apenas %d vertices (menor que 50).\n", V);
        fprintf(arq, "Grafo possui apenas %d vertices; conjunto {10..50} incompleto.\n", V);
        fclose(arq);
        return;
    }

    if (V > MAX_VERTICES_FLOYD) {
        /* 2 matrizes V x V (float + int) = 8 * V^2 bytes */
        double bytes_2_matrizes = (double)V * (double)V * (double)(sizeof(float) + sizeof(int));
        double gb_necessario = bytes_2_matrizes / (1024.0 * 1024.0 * 1024.0);
        printf("Status: Nao executado (V = %d excede MAX_VERTICES_FLOYD = %d; requer %.2f GB para 2 matrizes = 8*V^2 bytes).\n",
               V, MAX_VERTICES_FLOYD, gb_necessario);
        fprintf(arq, "Status: Nao executado\n");
        double est1 = estimar_tempo_fw_1x(grafo);
        fprintf(arq, "Motivo: V = %d excede o limite pratico MAX_VERTICES_FLOYD = %d (requer %.2f GB de RAM para 2 matrizes V x V = 8*V^2 bytes; tempo estimado ~%.2f s / %.2f h).\n",
                V, MAX_VERTICES_FLOYD, gb_necessario, est1, est1 / 3600.0);
        fclose(arq);
        return;
    }

    float *dist = alocar_matriz_distancias_1d(V);
    int *pred = alocar_matriz_predecessores_1d(V);
    if (!dist || !pred) {
        printf("Erro: Memoria insuficiente para matrizes %d x %d do Floyd-Warshall.\n", V, V);
        fprintf(arq, "Erro: Memoria insuficiente para matrizes %d x %d.\n", V, V);
        free(dist);
        free(pred);
        fclose(arq);
        return;
    }

    double t_ini = obter_tempo_segundos();
    bool sem_ciclo_negativo = floyd_warshall(grafo, dist, pred);
    double tempo_exec = obter_tempo_segundos() - t_ini;

    if (!sem_ciclo_negativo) {
        printf("Resultado: CICLO NEGATIVO DETECTADO! Nao existem caminhos minimos bem definidos.\n");
        fprintf(arq, "Status: CICLO NEGATIVO DETECTADO!\n");
        fprintf(arq, "Tempo de execucao: %.6f s\n", tempo_exec);
    } else {
        printf("Status: Concluido sem ciclo negativo (tempo: %.6f s)\n\n", tempo_exec);
        fprintf(arq, "Status: Concluido sem ciclo negativo\n");
        fprintf(arq, "Tempo de execucao: %.6f s\n\n", tempo_exec);

        fprintf(arq, "Tabela de Distancias Minimas (Origem x Destino):\n");
        fprintf(arq, "%-10s", "O \\ D");
        printf("%-10s", "O \\ D");
        for (int j = 0; j < k_alvo; j++) {
            fprintf(arq, " | %-12d", vertices_alvo[j]);
            printf(" | %-12d", vertices_alvo[j]);
        }
        fprintf(arq, "\n-----------+--------------+--------------+--------------+--------------+--------------\n");
        printf("\n-----------+--------------+--------------+--------------+--------------+--------------\n");

        for (int i = 0; i < k_alvo; i++) {
            int u = vertices_alvo[i];
            fprintf(arq, "%-10d", u);
            printf("%-10d", u);
            for (int j = 0; j < k_alvo; j++) {
                int v = vertices_alvo[j];
                float d = dist[(size_t)(u - 1) * (size_t)V + (size_t)(v - 1)];
                if (d == FLT_MAX) {
                    fprintf(arq, " | %-12s", "INF");
                    printf(" | %-12s", "INF");
                } else {
                    fprintf(arq, " | %-12.2f", d);
                    printf(" | %-12.2f", d);
                }
            }
            fprintf(arq, "\n");
            printf("\n");
        }

        fprintf(arq, "\nDetalhamento dos Caminhos Minimos:\n");
        for (int i = 0; i < k_alvo; i++) {
            for (int j = 0; j < k_alvo; j++) {
                int u = vertices_alvo[i];
                int v = vertices_alvo[j];
                size_t idx = (size_t)(u - 1) * (size_t)V + (size_t)(v - 1);
                float d = dist[idx];
                if (d == FLT_MAX) {
                    fprintf(arq, "%d -> %d: INF (Inalcancavel)\n", u, v);
                } else {
                    fprintf(arq, "%d -> %d (custo %.2f): ", u, v, d);
                    imprimir_caminho_vetor_1_indexado(&pred[(size_t)(u - 1) * (size_t)V], u - 1, v - 1, V, arq);
                    fprintf(arq, "\n");
                }
            }
        }
    }

    free(dist);
    free(pred);
    fclose(arq);
    printf("Resultados salvos em: %s\n", caminho_saida);
}

/* ============================================================================
 * ESTUDO DE CASO 3: Tempo medio (10 execucoes) - Caso 1 vs Caso 2
 * ============================================================================ */
void executar_tp3_caso3(Grafo *grafo, const char *nome_base) {
    garantir_diretorio_resultados();
    char caminho_saida[512];
    snprintf(caminho_saida, sizeof(caminho_saida), "%s/%s_tp3_caso3_tempo_destino100.txt", PASTA_RESULTADOS, nome_base);

    FILE *arq = fopen(caminho_saida, "w");
    if (!arq) {
        printf("Erro ao criar arquivo de saida '%s'.\n", caminho_saida);
        return;
    }

    int V = grafo->num_vertices;
    printf("\n=== ESTUDO DE CASO 3: TEMPO MEDIO (10 EXECUCOES) - DESTINO 100 ===\n");
    fprintf(arq, "=== ESTUDO DE CASO 3: TEMPO MEDIO DE EXECUCAO (10 RODADAS) ===\n");
    fprintf(arq, "Grafo: %s (V = %d, E = %d)\n", nome_base, V, grafo->num_arestas);
    fprintf(arq, "Nota: Tempos medidos sem leitura/escrita em disco (apenas processamento).\n\n");

    if (V < 100) {
        printf("Aviso: V = %d < 100. Destino 100 nao existe neste grafo.\n", V);
        fprintf(arq, "Nao executado: V = %d < 100.\n", V);
        fclose(arq);
        return;
    }

    /* 1. Benchmark Bellman-Ford (10 execucoes rumo ao destino 100 via CSR) */
    float *dist_bf = (float *)malloc((size_t)V * sizeof(float));
    int *pais_bf = (int *)malloc((size_t)V * sizeof(int));
    Fila *fila_bf = criar_fila(V);
    int *arestas_bf = (int *)malloc((size_t)V * sizeof(int));

    double soma_bf = 0.0;
    /* Nota: bf_ok guarda o retorno da ultima rodada; como o grafo e estatico
     * entre as 10 execucoes, a existencia de ciclo negativo e deterministica
     * e identica em todas as iteracoes. */
    bool bf_ok = true;

    double t_csr_ini = obter_tempo_segundos();
    GrafoCSR *csr = construir_grafo_transposto_csr(grafo);
    double tempo_csr = obter_tempo_segundos() - t_csr_ini;

    if (dist_bf && pais_bf && fila_bf && arestas_bf && csr) {
        for (int r = 0; r < NUM_EXECUCOES_BENCHMARK; r++) {
            double t0 = obter_tempo_segundos();
            bf_ok = bellman_ford_csr(csr, 100, dist_bf, pais_bf, fila_bf, arestas_bf);
            soma_bf += (obter_tempo_segundos() - t0);
        }
        double media_bf = soma_bf / NUM_EXECUCOES_BENCHMARK;
        printf("Bellman-Ford (10 execucoes -> destino 100): Total = %.6f s | Media = %.6f s (CSR build: %.6f s) [%s]\n",
               soma_bf, media_bf, tempo_csr, bf_ok ? "Sem ciclo negativo" : "Ciclo negativo");
        fprintf(arq, "1. Bellman-Ford (todos-para-100 via grafo transposto CSR):\n");
        fprintf(arq, "   - Status: %s\n", bf_ok ? "Sem ciclo negativo" : "Ciclo negativo detectado");
        fprintf(arq, "   - Tempo construcao CSR: %.6f s\n", tempo_csr);
        fprintf(arq, "   - Tempo total (10 execucoes): %.6f s\n", soma_bf);
        fprintf(arq, "   - Tempo medio por execucao: %.6f s\n\n", media_bf);
    }
    destruir_fila(fila_bf);
    free(arestas_bf);
    liberar_grafo_csr(csr);
    free(dist_bf);
    free(pais_bf);

    /* 2. Benchmark Floyd-Warshall (10 execucoes) */
    if (V > MAX_VERTICES_FLOYD) {
        /* 2 matrizes V x V (float + int) = 8 * V^2 bytes */
        double bytes_2_matrizes = (double)V * (double)V * (double)(sizeof(float) + sizeof(int));
        double gb_necessario = bytes_2_matrizes / (1024.0 * 1024.0 * 1024.0);
        printf("Floyd-Warshall: Nao executado (V = %d requer %.2f GB de RAM para 2 matrizes = 8*V^2 bytes).\n",
               V, gb_necessario);
        fprintf(arq, "2. Floyd-Warshall:\n");
        double est1 = estimar_tempo_fw_1x(grafo);
        fprintf(arq, "   - Status: Nao executado (V = %d excede MAX_VERTICES_FLOYD = %d; requer %.2f GB; tempo estimado ~%.2f s/rodada, ~%.2f h em 10 rodadas).\n",
                V, MAX_VERTICES_FLOYD, gb_necessario, est1, (est1 * NUM_EXECUCOES_BENCHMARK) / 3600.0);
    } else {
        float *dist_fw = alocar_matriz_distancias_1d(V);
        int *pred_fw = alocar_matriz_predecessores_1d(V);
        if (!dist_fw || !pred_fw) {
            printf("Floyd-Warshall: Memoria insuficiente para matrizes %d x %d.\n", V, V);
            fprintf(arq, "2. Floyd-Warshall: Memoria insuficiente.\n");
        } else {
            double soma_fw = 0.0;
            /* Mesma observacao: grafo estatico, retorno identico em todas as rodadas. */
            bool fw_ok = true;
            for (int r = 0; r < NUM_EXECUCOES_BENCHMARK; r++) {
                double t0 = obter_tempo_segundos();
                fw_ok = floyd_warshall(grafo, dist_fw, pred_fw);
                soma_fw += (obter_tempo_segundos() - t0);
            }
            double media_fw = soma_fw / NUM_EXECUCOES_BENCHMARK;
            printf("Floyd-Warshall (10 execucoes): Total = %.6f s | Media = %.6f s [%s]\n",
                   soma_fw, media_fw, fw_ok ? "Sem ciclo negativo" : "Ciclo negativo");
            fprintf(arq, "2. Floyd-Warshall (todos os pares):\n");
            fprintf(arq, "   - Status: %s\n", fw_ok ? "Sem ciclo negativo" : "Ciclo negativo detectado");
            fprintf(arq, "   - Tempo total (10 execucoes): %.6f s\n", soma_fw);
            fprintf(arq, "   - Tempo medio por execucao: %.6f s\n", media_fw);
        }
        free(dist_fw);
        free(pred_fw);
    }

    fclose(arq);
    printf("Resultados salvos em: %s\n", caminho_saida);
}

/* ============================================================================
 * ESTUDO DE CASO 4: Tempo medio (10 execucoes) - Todos os pares
 * ============================================================================ */
void executar_tp3_caso4(Grafo *grafo, const char *nome_base) {
    garantir_diretorio_resultados();
    char caminho_saida[512];
    snprintf(caminho_saida, sizeof(caminho_saida), "%s/%s_tp3_caso4_todos_pares.txt", PASTA_RESULTADOS, nome_base);

    FILE *arq = fopen(caminho_saida, "w");
    if (!arq) {
        printf("Erro ao criar arquivo de saida '%s'.\n", caminho_saida);
        return;
    }

    int V = grafo->num_vertices;
    printf("\n=== ESTUDO DE CASO 4: TODOS OS PARES (10 EXECUCOES) ===\n");
    fprintf(arq, "=== ESTUDO DE CASO 4: TODOS OS PARES DE VERTICES (10 EXECUCOES) ===\n");
    fprintf(arq, "Grafo: %s (V = %d, E = %d)\n\n", nome_base, V, grafo->num_arestas);

    /* Pico de RAM no Caso 4: 3 matrizes simultaneas V x V (dist_bf + dist_fw + pred_fw)
     * = 12 * V^2 bytes (pois pred_bf e liberado antes de alocar Floyd-Warshall),
     * alem de exigir tempo O(V^3) x 10 rodadas. */
    if (V > MAX_VERTICES_CASO4) {
        double bytes_pico_3_matrizes = (double)V * (double)V * (double)(2 * sizeof(float) + sizeof(int));
        double gb_pico = bytes_pico_3_matrizes / (1024.0 * 1024.0 * 1024.0);
        double gb_4_matrizes = (4.0 * 4.0 * (double)V * (double)V) / (1024.0 * 1024.0 * 1024.0);

        printf("Status: Nao executado (V = %d excede o limite pratico MAX_VERTICES_CASO4 = %d).\n",
               V, MAX_VERTICES_CASO4);
        printf("Memoria necessaria: %.2f GB de pico (3 matrizes = 12*V^2 bytes; %.2f GB se 4 matrizes = 16*V^2 bytes).\n",
               gb_pico, gb_4_matrizes);

        fprintf(arq, "Status: Nao executado\n");
        fprintf(arq, "Motivo: V = %d excede o limite pratico MAX_VERTICES_CASO4 = %d.\n", V, MAX_VERTICES_CASO4);
        fprintf(arq, "Memoria necessaria: %.2f GB de pico (3 matrizes simultaneas = 12*V^2 bytes; %.2f GB para 4 matrizes = 16*V^2 bytes).\n\n",
                gb_pico, gb_4_matrizes);

        /* 1. Extrapolacao do Bellman-Ford todos os pares a partir de 10 execucoes reais 1-para-todos no grafo completo */
        GrafoCSR *csr = construir_grafo_csr(grafo);
        float *d_sample = (float *)malloc((size_t)V * sizeof(float));
        int *p_sample = (int *)malloc((size_t)V * sizeof(int));
        Fila *f_sample = criar_fila(V);
        int *a_sample = (int *)malloc((size_t)V * sizeof(int));

        double t_bf_medio_1_origem = 0.0;
        bool bf_amostra_ok = true;
        if (csr && d_sample && p_sample && f_sample && a_sample) {
            for (int i = 0; i < NUM_EXECUCOES_BENCHMARK; i++) {
                int origem_amostra = 1 + (int)(((long long)i * V) / NUM_EXECUCOES_BENCHMARK);
                double t0 = obter_tempo_segundos();
                bool ok = bellman_ford_csr(csr, origem_amostra, d_sample, p_sample, f_sample, a_sample);
                t_bf_medio_1_origem += (obter_tempo_segundos() - t0);
                if (!ok) bf_amostra_ok = false;
            }
            t_bf_medio_1_origem /= NUM_EXECUCOES_BENCHMARK;
        }
        destruir_fila(f_sample);
        free(a_sample);
        free(d_sample);
        free(p_sample);
        liberar_grafo_csr(csr);

        double est_bf_1x = t_bf_medio_1_origem * (double)V;
        double est_bf_10x = est_bf_1x * NUM_EXECUCOES_BENCHMARK;

        /* 2. Extrapolacao cubica do Floyd-Warshall calibrada em grafo conexo de V_sub = 500 vertices */
        int v_sub = (V >= SUBGRAFO_CALIBRACAO_V) ? SUBGRAFO_CALIBRACAO_V : V;
        int grau_medio = (V > 0) ? (int)((double)grafo->num_arestas / (double)V + 0.5) : 2;
        Grafo *cal_fw = criar_grafo_calibracao_fw(v_sub, grau_medio);
        double t_cal_fw = 0.0;
        if (cal_fw) {
            float *s_dfw = alocar_matriz_distancias_1d(v_sub);
            int *s_pfw = alocar_matriz_predecessores_1d(v_sub);
            if (s_dfw && s_pfw) {
                double t0 = obter_tempo_segundos();
                floyd_warshall(cal_fw, s_dfw, s_pfw);
                t_cal_fw = obter_tempo_segundos() - t0;
            }
            free(s_dfw);
            free(s_pfw);
            liberar_grafo(cal_fw);
        }

        double razao_v = (double)V / (double)v_sub;
        double fator_cubico = razao_v * razao_v * razao_v;
        double est_fw_1x = t_cal_fw * fator_cubico;
        double est_fw_10x = est_fw_1x * NUM_EXECUCOES_BENCHMARK;

        printf("Estimativa de tempo por extrapolacao empirica nesta maquina (NAO e medicao direta do Caso 4 completo):\n");
        printf("  - Floyd-Warshall (calibracao conexa V_sub=%d -> escala cubica (V/V_sub)^3): ~%.2f s (%.2f h) por rodada | ~%.2f s (%.2f h) em 10 rodadas\n",
               v_sub, est_fw_1x, est_fw_1x / 3600.0, est_fw_10x, est_fw_10x / 3600.0);
        printf("  - Bellman-Ford todos os pares (10 origens reais no grafo completo x V):   ~%.2f s (%.2f h) por rodada | ~%.2f s (%.2f h) em 10 rodadas [%s]\n",
               est_bf_1x, est_bf_1x / 3600.0, est_bf_10x, est_bf_10x / 3600.0,
               bf_amostra_ok ? "sem ciclo neg. na amostra" : "ciclo neg. detectado na amostra");

        fprintf(arq, "Estimativa de tempo por extrapolacao empirica nesta maquina:\n");
        fprintf(arq, "AVISO: Os valores abaixo sao extrapolacoes matematicas calibradas nesta maquina, e NAO medicoes diretas do grafo completo.\n\n");
        fprintf(arq, "1. Floyd-Warshall (extrapolacao cubica Theta(V^3) a partir de calibracao conexa V_sub = %d, fator = (V/V_sub)^3 = %.1f):\n",
                v_sub, fator_cubico);
        fprintf(arq, "   - Tempo medido na calibracao (1 rodada, V_sub = %d): %.6f s\n", v_sub, t_cal_fw);
        fprintf(arq, "   - Tempo estimado (1 rodada): %.2f s (%.2f min / %.2f h)\n",
                est_fw_1x, est_fw_1x / 60.0, est_fw_1x / 3600.0);
        fprintf(arq, "   - Tempo estimado (10 rodadas): %.2f s (%.2f min / %.2f h)\n\n",
                est_fw_10x, est_fw_10x / 60.0, est_fw_10x / 3600.0);

        fprintf(arq, "2. Bellman-Ford para todos os pares (extrapolacao V * t_medio_1_origem a partir de 10 origens reais no grafo completo):\n");
        fprintf(arq, "   - Tempo medio real (1 origem no grafo completo): %.6f s (%s)\n",
                t_bf_medio_1_origem, bf_amostra_ok ? "sem ciclo negativo" : "ciclo negativo detectado");
        fprintf(arq, "   - Tempo estimado (1 rodada = V execucoes): %.2f s (%.2f min / %.2f h)\n",
                est_bf_1x, est_bf_1x / 60.0, est_bf_1x / 3600.0);
        fprintf(arq, "   - Tempo estimado (10 rodadas): %.2f s (%.2f min / %.2f h)\n",
                est_bf_10x, est_bf_10x / 60.0, est_bf_10x / 3600.0);

        fclose(arq);
        printf("Resultados salvos em: %s\n", caminho_saida);
        return;
    }

    /* 1. Bellman-Ford executado a partir de todos os V vertices (10 rodadas completas) */
    float *dist_bf = alocar_matriz_distancias_1d(V);
    int *pred_bf = alocar_matriz_predecessores_1d(V);
    if (!dist_bf || !pred_bf) {
        double gb_pico = ((double)V * (double)V * 12.0) / (1024.0 * 1024.0 * 1024.0);
        printf("Erro: Memoria insuficiente para matrizes %d x %d (%.2f GB necessarios).\n", V, V, gb_pico);
        fprintf(arq, "Status: Nao executado (memoria insuficiente para matrizes %d x %d = %.2f GB).\n", V, V, gb_pico);
        free(dist_bf);
        free(pred_bf);
        fclose(arq);
        return;
    }

    double soma_bf = 0.0;
    bool bf_ok = true;
    for (int r = 0; r < NUM_EXECUCOES_BENCHMARK; r++) {
        double t0 = obter_tempo_segundos();
        bf_ok = bellman_ford_todos_pares(grafo, dist_bf, pred_bf);
        soma_bf += (obter_tempo_segundos() - t0);
    }
    double media_bf = soma_bf / NUM_EXECUCOES_BENCHMARK;

    /* Libera pred_bf antes de alocar Floyd-Warshall, reduzindo o pico de RAM de 16*V^2 para 12*V^2 bytes */
    free(pred_bf);
    pred_bf = NULL;

    /* 2. Floyd-Warshall executado 10 vezes */
    float *dist_fw = alocar_matriz_distancias_1d(V);
    int *pred_fw = alocar_matriz_predecessores_1d(V);
    if (!dist_fw || !pred_fw) {
        double gb_pico = ((double)V * (double)V * 12.0) / (1024.0 * 1024.0 * 1024.0);
        printf("Erro: Memoria insuficiente para matrizes do Floyd-Warshall %d x %d (%.2f GB de pico).\n", V, V, gb_pico);
        fprintf(arq, "Status: Nao executado (memoria insuficiente para Floyd-Warshall %d x %d = %.2f GB).\n", V, V, gb_pico);
        free(dist_bf);
        free(dist_fw);
        free(pred_fw);
        fclose(arq);
        return;
    }

    double soma_fw = 0.0;
    bool fw_ok = true;
    for (int r = 0; r < NUM_EXECUCOES_BENCHMARK; r++) {
        double t0 = obter_tempo_segundos();
        fw_ok = floyd_warshall(grafo, dist_fw, pred_fw);
        soma_fw += (obter_tempo_segundos() - t0);
    }
    double media_fw = soma_fw / NUM_EXECUCOES_BENCHMARK;

    printf("Bellman-Ford (V execucoes x 10 rodadas): Total = %.6f s | Media = %.6f s [%s]\n",
           soma_bf, media_bf, bf_ok ? "Sem ciclo negativo" : "Ciclo negativo");
    printf("Floyd-Warshall (10 rodadas):             Total = %.6f s | Media = %.6f s [%s]\n",
           soma_fw, media_fw, fw_ok ? "Sem ciclo negativo" : "Ciclo negativo");

    fprintf(arq, "1. Bellman-Ford (executado V = %d vezes por rodada):\n", V);
    fprintf(arq, "   - Status: %s\n", bf_ok ? "Sem ciclo negativo" : "Ciclo negativo detectado");
    fprintf(arq, "   - Tempo total (10 rodadas): %.6f s\n", soma_bf);
    fprintf(arq, "   - Tempo medio por rodada: %.6f s\n\n", media_bf);

    fprintf(arq, "2. Floyd-Warshall:\n");
    fprintf(arq, "   - Status: %s\n", fw_ok ? "Sem ciclo negativo" : "Ciclo negativo detectado");
    fprintf(arq, "   - Tempo total (10 rodadas): %.6f s\n", soma_fw);
    fprintf(arq, "   - Tempo medio por rodada: %.6f s\n\n", media_fw);

    if (bf_ok && fw_ok) {
        int divergencias = 0;
        size_t total_pares = (size_t)V * (size_t)V;
        for (size_t idx = 0; idx < total_pares; idx++) {
            float d1 = dist_bf[idx];
            float d2 = dist_fw[idx];
            if (d1 == FLT_MAX && d2 == FLT_MAX) continue;
            float diff = (d1 > d2) ? (d1 - d2) : (d2 - d1);
            if (diff > 1e-3f) divergencias++;
        }
        printf("Verificacao de consistencia entre algoritmos: %d divergencias em %zu pares.\n", divergencias, total_pares);
        fprintf(arq, "Verificacao de consistencia (Bellman-Ford vs Floyd-Warshall): %d divergencias em %zu pares.\n",
                divergencias, total_pares);
    }

    free(dist_bf);
    free(dist_fw);
    free(pred_fw);
    fclose(arq);
    printf("Resultados salvos em: %s\n", caminho_saida);
}

static void gerar_resumo_consolidado_tp3(const char **bases, int qtd_bases) {
    garantir_diretorio_resultados();
    char caminho_resumo[512];
    snprintf(caminho_resumo, sizeof(caminho_resumo), "%s/RESUMO_GERAL_TP3.txt", PASTA_RESULTADOS);
    FILE *out = fopen(caminho_resumo, "w");
    if (!out) return;
    fprintf(out, "=== RESUMO CONSOLIDADO - TRABALHO PRATICO 3 ===\n\n");
    const char *suf[] = {"tp3_caso1_bellman_ford.txt", "tp3_caso2_floyd_warshall.txt",
                         "tp3_caso3_tempo_destino100.txt", "tp3_caso4_todos_pares.txt"};
    for (int b = 0; b < qtd_bases; b++) {
        fprintf(out, "==================== GRAFO: %s ====================\n\n", bases[b]);
        for (int s = 0; s < 4; s++) {
            char arq[512];
            snprintf(arq, sizeof(arq), "%s/%s_%s", PASTA_RESULTADOS, bases[b], suf[s]);
            FILE *in = fopen(arq, "r");
            if (in) {
                char linha[1024];
                while (fgets(linha, sizeof(linha), in)) fputs(linha, out);
                fprintf(out, "\n----------------------------------------------------\n\n");
                fclose(in);
            }
        }
    }
    fclose(out);
    printf("\n>>> RESUMO CONSOLIDADO SALVO EM: %s <<<\n", caminho_resumo);
}

static int executar_lote_tp3_automatico(void) {
    const char *cands[] = {"grafo_W_1.txt", "grafo_W_2.txt", "grafo_W_3.txt", "grafo_W_4.txt", "grafo_W_5.txt"};
    const char *bases[5];
    char storage[5][256];
    int total = 0;
    garantir_diretorio_resultados();
    for (int i = 0; i < 5; i++) {
        char cam[512];
        if (!resolver_caminho_grafo(cands[i], cam, sizeof(cam))) continue;
        printf("\n=== CARREGANDO: %s ===\n", cam);
        Grafo *g = ler_grafo_do_arquivo(cam, LISTA_ADJACENCIA, 1);
        if (!g) continue;
        extrair_nome_base(cam, storage[total], sizeof(storage[total]));
        bases[total] = storage[total];
        executar_tp3_caso1(g, bases[total]);
        executar_tp3_caso2(g, bases[total]);
        executar_tp3_caso3(g, bases[total]);
        executar_tp3_caso4(g, bases[total]);
        liberar_grafo(g);
        total++;
    }
    if (total == 0) { printf("Nenhum grafo_W_1..5 encontrado.\n"); return 1; }
    gerar_resumo_consolidado_tp3(bases, total);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "--tp3-todos") == 0) return executar_lote_tp3_automatico();
    if (argc >= 2) {
        char cam[512], base[256];
        resolver_caminho_grafo(argv[1], cam, sizeof(cam));
        Grafo *g = ler_grafo_do_arquivo(cam, LISTA_ADJACENCIA, 1);
        if (!g) { printf("Erro ao carregar '%s'!\n", cam); return 1; }
        extrair_nome_base(cam, base, sizeof(base));
        garantir_diretorio_resultados();
        executar_tp3_caso1(g, base); executar_tp3_caso2(g, base);
        executar_tp3_caso3(g, base); executar_tp3_caso4(g, base);
        liberar_grafo(g);
        return 0;
    }

    int opcao_rep, direcionado;
    char nome_arquivo[256], caminho_resolvido[512], nome_base[256];

    printf("Escolha a representacao do grafo (1 - Lista de Adjacencia, 2 - Matriz de Adjacencia): ");
    if (scanf("%d", &opcao_rep) != 1) return 1;

    RepresentacaoGrafo tipo_rep = (opcao_rep == 1) ? LISTA_ADJACENCIA : MATRIZ_ADJACENCIA;

    printf("O grafo e direcionado? (1 - Sim, 0 - Nao): ");
    if (scanf("%d", &direcionado) != 1) return 1;

    printf("Digite o nome do arquivo do grafo: ");
    if (scanf("%255s", nome_arquivo) != 1) return 1;

    resolver_caminho_grafo(nome_arquivo, caminho_resolvido, sizeof(caminho_resolvido));
    Grafo *grafo = ler_grafo_do_arquivo(caminho_resolvido, tipo_rep, direcionado);
    if (!grafo) {
        printf("Erro ao carregar o grafo '%s'!\n", caminho_resolvido);
        return 1;
    }

    extrair_nome_base(caminho_resolvido, nome_base, sizeof(nome_base));
    garantir_diretorio_resultados();

    int opcao;
    do {
        mostrar_menu();
        if (scanf("%d", &opcao) != 1) break;

        switch (opcao) {
            case 1: {
                char arquivo_saida[512];
                snprintf(arquivo_saida, sizeof(arquivo_saida), "%s/%s_informacoes_grafo.txt", PASTA_RESULTADOS, nome_base);
                salvar_informacoes_grafo(grafo, arquivo_saida);
                printf("Informacoes basicas salvas em '%s'.\n", arquivo_saida);
                break;
            }
            case 2: {
                int vertice_inicial;
                printf("Digite o vertice inicial para BFS: ");
                if (scanf("%d", &vertice_inicial) != 1) {
                    printf("Entrada invalida.\n");
                    break;
                }

                int *pais = (int *)malloc((size_t)grafo->num_vertices * sizeof(int));
                int *niveis = (int *)malloc((size_t)grafo->num_vertices * sizeof(int));
                busca_em_largura(grafo, vertice_inicial, pais, niveis);

                char arquivo_saida[512];
                snprintf(arquivo_saida, sizeof(arquivo_saida), "%s/%s_bfs_resultado.txt", PASTA_RESULTADOS, nome_base);
                FILE *arquivo = fopen(arquivo_saida, "w");
                if (arquivo) {
                    fprintf(arquivo, "Resultado da BFS a partir do vertice %d:\n", vertice_inicial);
                    for (int i = 0; i < grafo->num_vertices; i++) {
                        fprintf(arquivo, "Vertice %d: Pai = %d, Nivel = %d\n", i + 1, pais[i] != -1 ? pais[i] + 1 : -1, niveis[i]);
                    }
                    fclose(arquivo);
                }

                free(pais);
                free(niveis);
                printf("Resultado da BFS salvo em '%s'.\n", arquivo_saida);
                break;
            }
            case 3: {
                int vertice_inicial;
                printf("Digite o vertice inicial para DFS: ");
                if (scanf("%d", &vertice_inicial) != 1) {
                    printf("Entrada invalida.\n");
                    break;
                }

                int *pais = (int *)malloc((size_t)grafo->num_vertices * sizeof(int));
                int *niveis = (int *)malloc((size_t)grafo->num_vertices * sizeof(int));
                busca_em_profundidade(grafo, vertice_inicial, pais, niveis);

                char arquivo_saida[512];
                snprintf(arquivo_saida, sizeof(arquivo_saida), "%s/%s_dfs_resultado.txt", PASTA_RESULTADOS, nome_base);
                FILE *arquivo = fopen(arquivo_saida, "w");
                if (arquivo) {
                    fprintf(arquivo, "Resultado da DFS a partir do vertice %d:\n", vertice_inicial);
                    for (int i = 0; i < grafo->num_vertices; i++) {
                        fprintf(arquivo, "Vertice %d: Pai = %d, Nivel = %d\n", i + 1, pais[i] != -1 ? pais[i] + 1 : -1, niveis[i]);
                    }
                    fclose(arquivo);
                }

                free(pais);
                free(niveis);
                printf("Resultado da DFS salvo em '%s'.\n", arquivo_saida);
                break;
            }
            case 4: {
                int origem, destino;
                printf("Digite os vertices de origem e destino: ");
                if (scanf("%d %d", &origem, &destino) != 2) {
                    printf("Entrada invalida.\n");
                    break;
                }
                int distancia = calcular_distancia(grafo, origem, destino);
                if (distancia != -1) {
                    printf("A distancia entre %d e %d e: %d\n", origem, destino, distancia);
                } else {
                    printf("Nao ha caminho entre %d e %d.\n", origem, destino);
                }
                break;
            }
            case 5: {
                int diametro = calcular_diametro(grafo);
                printf("O diametro do grafo e: %d\n", diametro);
                break;
            }
            case 6: {
                char arquivo_saida[512];
                snprintf(arquivo_saida, sizeof(arquivo_saida), "%s/%s_componentes_conexas.txt", PASTA_RESULTADOS, nome_base);
                encontrar_componentes_conexas(grafo, arquivo_saida);
                printf("Componentes conexas salvas em '%s'.\n", arquivo_saida);
                break;
            }
            case 7: {
                char arquivo_saida[512];
                snprintf(arquivo_saida, sizeof(arquivo_saida), "%s/%s_memoria_grafo.txt", PASTA_RESULTADOS, nome_base);
                salvar_memoria_grafo(grafo, direcionado, arquivo_saida);
                printf("Memoria utilizada pela representacao salva em '%s'.\n", arquivo_saida);
                break;
            }
            case 8: {
                int *pais = (int *)malloc((size_t)grafo->num_vertices * sizeof(int));
                int *niveis = (int *)malloc((size_t)grafo->num_vertices * sizeof(int));
                double t0 = obter_tempo_segundos();
                srand((unsigned int)time(NULL));
                for (int i = 0; i < 100; i++) {
                    int vertice_inicial = (rand() % grafo->num_vertices) + 1;
                    busca_em_largura(grafo, vertice_inicial, pais, niveis);
                }
                double tempo_total = obter_tempo_segundos() - t0;
                printf("Tempo medio de 100 execucoes de BFS: %f segundos\n", tempo_total / 100.0);
                free(pais);
                free(niveis);
                break;
            }
            case 9: {
                int *pais = (int *)malloc((size_t)grafo->num_vertices * sizeof(int));
                int *niveis = (int *)malloc((size_t)grafo->num_vertices * sizeof(int));
                double t0 = obter_tempo_segundos();
                srand((unsigned int)time(NULL));
                for (int i = 0; i < 100; i++) {
                    int vertice_inicial = (rand() % grafo->num_vertices) + 1;
                    busca_em_profundidade(grafo, vertice_inicial, pais, niveis);
                }
                double tempo_total = obter_tempo_segundos() - t0;
                printf("Tempo medio de 100 execucoes de DFS: %f segundos\n", tempo_total / 100.0);
                free(pais);
                free(niveis);
                break;
            }
            case 10: {
                if (possui_arestas_negativas(grafo)) {
                    printf("Erro: O grafo possui arestas com pesos negativos. Dijkstra nao e aplicavel.\n");
                    break;
                }
                int vertice_inicial;
                printf("Digite o vertice de origem (1 a %d): ", grafo->num_vertices);
                if (scanf("%d", &vertice_inicial) != 1) {
                    printf("Entrada invalida.\n");
                    break;
                }

                if (vertice_inicial < 1 || vertice_inicial > grafo->num_vertices) {
                    printf("Vertice invalido.\n");
                    break;
                }

                int V = grafo->num_vertices;
                float *distancias = (float *)malloc((size_t)V * sizeof(float));
                int *pais = (int *)malloc((size_t)V * sizeof(int));

                dijkstra(grafo, vertice_inicial, distancias, pais);

                char arquivo_saida[512];
                snprintf(arquivo_saida, sizeof(arquivo_saida), "%s/%s_dijkstra_todos.txt", PASTA_RESULTADOS, nome_base);
                FILE *arquivo = fopen(arquivo_saida, "w");
                if (arquivo) {
                    fprintf(arquivo, "Resultados de Dijkstra a partir do vertice %d:\n\n", vertice_inicial);
                    for (int i = 0; i < V; i++) {
                        fprintf(arquivo, "Destino %d: ", i + 1);
                        if (distancias[i] == FLT_MAX) {
                            fprintf(arquivo, "Distancia = INFINITO | Caminho = Nenhum\n");
                        } else {
                            fprintf(arquivo, "Distancia = %.2f | Caminho = ", distancias[i]);
                            salvar_caminho_minimo(grafo, pais, vertice_inicial - 1, i, V, arquivo);
                        }
                    }
                    fclose(arquivo);
                    printf("Resultados salvos em '%s'.\n", arquivo_saida);
                }
                free(distancias);
                free(pais);
                break;
            }
            case 11: {
                if (possui_arestas_negativas(grafo)) {
                    printf("Erro: O grafo possui arestas com pesos negativos. Dijkstra nao e aplicavel.\n");
                    break;
                }
                int origem, destino;
                printf("Digite o vertice de origem (1 a %d): ", grafo->num_vertices);
                if (scanf("%d", &origem) != 1) {
                    printf("Entrada invalida.\n");
                    break;
                }
                printf("Digite o vertice de destino (1 a %d): ", grafo->num_vertices);
                if (scanf("%d", &destino) != 1) {
                    printf("Entrada invalida.\n");
                    break;
                }

                if (origem < 1 || origem > grafo->num_vertices || destino < 1 || destino > grafo->num_vertices) {
                    printf("Vertices invalidos.\n");
                    break;
                }

                int V = grafo->num_vertices;
                float *distancias = (float *)malloc((size_t)V * sizeof(float));
                int *pais = (int *)malloc((size_t)V * sizeof(int));

                dijkstra(grafo, origem, distancias, pais);

                int dest_idx = destino - 1;
                if (distancias[dest_idx] == FLT_MAX) {
                    printf("Nao existe caminho entre %d e %d.\n", origem, destino);
                } else {
                    printf("\nDistancia minima de %d para %d: %.2f\n", origem, destino, distancias[dest_idx]);
                    printf("Caminho: ");
                    salvar_caminho_minimo(grafo, pais, origem - 1, dest_idx, V, stdout);
                }
                free(distancias);
                free(pais);
                break;
            }
            case 12: {
                if (possui_arestas_negativas(grafo)) {
                    printf("Erro: O grafo possui arestas com pesos negativos. Dijkstra nao e aplicavel.\n");
                    break;
                }

                int V = grafo->num_vertices;
                int origem = 1;
                int destinos[] = {10, 20, 30, 40, 50};
                int num_destinos = 5;

                float *distancias = (float *)malloc((size_t)V * sizeof(float));
                int *pais = (int *)malloc((size_t)V * sizeof(int));

                dijkstra(grafo, origem, distancias, pais);

                char arquivo_saida[512];
                snprintf(arquivo_saida, sizeof(arquivo_saida), "%s/%s_dijkstra_alvos.txt", PASTA_RESULTADOS, nome_base);
                FILE *arquivo = fopen(arquivo_saida, "w");
                if (arquivo) {
                    fprintf(arquivo, "Caminhos minimos a partir do vertice %d (COM HEAP):\n\n", origem);
                    for (int i = 0; i < num_destinos; i++) {
                        int d = destinos[i];
                        if (d > V) continue;
                        int d_idx = d - 1;
                        fprintf(arquivo, "--- Destino %d ---\n", d);
                        if (distancias[d_idx] == FLT_MAX) {
                            fprintf(arquivo, "Distancia: INFINITO\nCaminho: Inalcancavel\n\n");
                        } else {
                            fprintf(arquivo, "Distancia: %.2f\nCaminho: ", distancias[d_idx]);
                            salvar_caminho_minimo(grafo, pais, origem - 1, d_idx, V, arquivo);
                            fprintf(arquivo, "\n");
                        }
                    }
                    fclose(arquivo);
                    printf("Resultados salvos em '%s'.\n", arquivo_saida);
                }
                free(distancias);
                free(pais);
                break;
            }
            case 13: {
                if (possui_arestas_negativas(grafo)) {
                    printf("Erro: O grafo possui arestas com pesos negativos. Dijkstra nao e aplicavel.\n");
                    break;
                }

                int V = grafo->num_vertices;
                int origem = 1;
                int destinos[] = {10, 20, 30, 40, 50};
                int num_destinos = 5;

                float *distancias = (float *)malloc((size_t)V * sizeof(float));
                int *pais = (int *)malloc((size_t)V * sizeof(int));

                dijkstra_sem_heap(grafo, origem, distancias, pais);

                char arquivo_saida[512];
                snprintf(arquivo_saida, sizeof(arquivo_saida), "%s/%s_dijkstra_sem_heap_alvos.txt", PASTA_RESULTADOS, nome_base);
                FILE *arquivo = fopen(arquivo_saida, "w");
                if (arquivo) {
                    fprintf(arquivo, "Caminhos minimos a partir do vertice %d (SEM HEAP):\n\n", origem);
                    for (int i = 0; i < num_destinos; i++) {
                        int d = destinos[i];
                        if (d > V) continue;
                        int d_idx = d - 1;
                        fprintf(arquivo, "--- Destino %d ---\n", d);
                        if (distancias[d_idx] == FLT_MAX) {
                            fprintf(arquivo, "Distancia: INFINITO\nCaminho: Inalcancavel\n\n");
                        } else {
                            fprintf(arquivo, "Distancia: %.2f\nCaminho: ", distancias[d_idx]);
                            salvar_caminho_minimo(grafo, pais, origem - 1, d_idx, V, arquivo);
                            fprintf(arquivo, "\n");
                        }
                    }
                    fclose(arquivo);
                    printf("Resultados salvos em '%s'.\n", arquivo_saida);
                }
                free(distancias);
                free(pais);
                break;
            }
            case 14: {
                if (possui_arestas_negativas(grafo)) {
                    printf("Erro: O grafo possui arestas com pesos negativos.\n");
                    break;
                }

                int k;
                printf("Digite o numero de vertices aleatorios (k) para o benchmark: ");
                if (scanf("%d", &k) != 1 || k <= 0) {
                    printf("Numero invalido.\n");
                    break;
                }

                int V = grafo->num_vertices;
                int *amostra = (int *)malloc((size_t)k * sizeof(int));
                srand((unsigned int)time(NULL));
                for (int i = 0; i < k; i++) {
                    amostra[i] = (rand() % V) + 1;
                }

                float *distancias = (float *)malloc((size_t)V * sizeof(float));
                int *pais = (int *)malloc((size_t)V * sizeof(int));

                double t0 = obter_tempo_segundos();
                for (int i = 0; i < k; i++) {
                    dijkstra(grafo, amostra[i], distancias, pais);
                }
                double tempo_com_heap = obter_tempo_segundos() - t0;

                t0 = obter_tempo_segundos();
                for (int i = 0; i < k; i++) {
                    dijkstra_sem_heap(grafo, amostra[i], distancias, pais);
                }
                double tempo_sem_heap = obter_tempo_segundos() - t0;

                printf("\n=== RESULTADO DO BENCHMARK (%d execucoes) ===\n", k);
                printf("Com Heap (MinHeap): Total = %.4fs | Medio = %.6fs\n", tempo_com_heap, tempo_com_heap / k);
                printf("Sem Heap (Vetor):   Total = %.4fs | Medio = %.6fs\n", tempo_sem_heap, tempo_sem_heap / k);
                if (tempo_com_heap > 0) {
                    printf("Speedup (Sem Heap / Com Heap): %.2fx\n", tempo_sem_heap / tempo_com_heap);
                }

                free(amostra);
                free(distancias);
                free(pais);
                break;
            }
            case 15:
                executar_tp3_caso1(grafo, nome_base);
                break;
            case 16:
                executar_tp3_caso2(grafo, nome_base);
                break;
            case 17:
                executar_tp3_caso3(grafo, nome_base);
                break;
            case 18:
                executar_tp3_caso4(grafo, nome_base);
                break;
            case 19:
                executar_tp3_caso1(grafo, nome_base);
                executar_tp3_caso2(grafo, nome_base);
                executar_tp3_caso3(grafo, nome_base);
                executar_tp3_caso4(grafo, nome_base);
                break;
            case 0:
                printf("Saindo...\n");
                break;
            default:
                printf("Opcao invalida!\n");
                break;
        }
    } while (opcao != 0);

    liberar_grafo(grafo);
    return 0;
}
