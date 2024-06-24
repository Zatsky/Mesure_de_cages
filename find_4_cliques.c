#include "find_4_cliques.h"

// Fonction pour lire le fichier CSV et construire la matrice d'adjacence
void read_graph_from_csv(char *filename, int ***adj_matrix, int *num_vertices) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Unable to open file");
        exit(EXIT_FAILURE);
    }

    char line[1024];

    // Lire le nombre de sommets et d'arêtes
    fgets(line, sizeof(line), file);
    char *token = strtok(line, "_");
    *num_vertices = atoi(token);

    // Allouer dynamiquement la matrice d'adjacence
    *adj_matrix = (int **)malloc(*num_vertices * sizeof(int *));
    for (int i = 0; i < *num_vertices; i++) {
        (*adj_matrix)[i] = (int *)calloc(*num_vertices, sizeof(int));
    }

    // Ignorer les lignes décrivant les sommets
    for (int i = 0; i < *num_vertices; i++) {
        fgets(line, sizeof(line), file);
    }

    // Lire les arêtes et remplir la matrice d'adjacence
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '#') break; // fin du graphe
        token = strtok(line, "_");
        int vertex1 = atoi(token);
        token = strtok(NULL, "_");
        int vertex2 = atoi(token);
        // Ignorer le type de liaison
        (*adj_matrix)[vertex1][vertex2] = 1;
        (*adj_matrix)[vertex2][vertex1] = 1; // Graphe non orienté
    }

    fclose(file);
}

// Fonction pour imprimer la matrice d'adjacence
void print_adj_matrix(int **adj_matrix, int num_vertices) {
    printf("Adjacency Matrix:\n");
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < num_vertices; j++) {
            printf("%d ", adj_matrix[i][j]);
        }
        printf("\n");
    }
}

// Fonction pour vérifier si les sommets donnés forment une clique
int is_clique(int **adj_matrix, int vertices[], int k) {
    for (int i = 0; i < k; i++) {
        for (int j = i + 1; j < k; j++) {
            if (!adj_matrix[vertices[i]][vertices[j]]) {
                return 0;
            }
        }
    }
    return 1;
}

// Fonction pour trouver les cliques de taille 4
int find_cliques_of_size_4(int **adj_matrix, int num_vertices) {
    int vertices[4];
    int nb_cliques = 0;
    for (int i = 0; i < num_vertices - 3; i++) {
        for (int j = i + 1; j < num_vertices - 2; j++) {
            for (int k = j + 1; k < num_vertices - 1; k++) {
                for (int l = k + 1; l < num_vertices; l++) {
                    vertices[0] = i;
                    vertices[1] = j;
                    vertices[2] = k;
                    vertices[3] = l;
                    if (is_clique(adj_matrix, vertices, 4)) {
                        nb_cliques++;
                    }
                }
            }
        }
    }
    return nb_cliques;
}

int main(int argc, char *argv[]) {
    int **adj_matrix;
    int num_vertices;
    
    read_graph_from_csv(argv[1], &adj_matrix, &num_vertices);
    // print_adj_matrix(adj_matrix, num_vertices);
    int nb_cliques = find_cliques_of_size_4(adj_matrix, num_vertices);

    // Libérer la mémoire allouée dynamiquement
    for (int i = 0; i < num_vertices; i++) {
        free(adj_matrix[i]);
    }
    free(adj_matrix);
    printf("%d",nb_cliques);
    return 0;
}
