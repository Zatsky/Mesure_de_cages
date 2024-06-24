#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void read_graph_from_csv(char *filename, int ***adj_matrix, int *num_vertices);
void print_adj_matrix(int **adj_matrix, int num_vertices);
int is_clique(int **adj_matrix, int vertices[], int k);
int find_cliques_of_size_4(int **adj_matrix, int num_vertices);