#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include "cagitude.h"

void init_sommet(SOMMET_COIN* sommet, int id ) {
    sommet->id = id;
    sommet->ids = (int*)malloc(3 * sizeof(int));
    sommet->poids = (int*)malloc(3 * sizeof(int)); // assuming a fixed size of 3 for the cycles
    sommet->liaisons_communs = (int*)malloc(3 * sizeof(int)); // assuming a fixed size of 3 for the common connections

    if (sommet->ids == NULL || sommet->poids == NULL || sommet->liaisons_communs == NULL) {
        // Handle memory allocation failure
        fprintf(stderr, "Failed to allocate memory\n");
        exit(1);
    }
    for (int i = 0; i < 3; i++) {
        sommet->ids[i] = -1;
        sommet->poids[i] = 0; // or some other default value
        sommet->liaisons_communs[i] = 0; // or some other default value
    }
}

void free_sommet(SOMMET_COIN* sommet) {
    free(sommet->ids);
    free(sommet->poids);
    free(sommet->liaisons_communs);
}

float calcul_mesure_coins_add(SOMMET_COIN* sommet, float alpha) {
    float total = 0.0;
    for (int i = 0; i < 3; i++) {
        if (sommet->liaisons_communs[i] == 1) {
            total += 1.0;
        } else if (sommet->liaisons_communs[i] == 2) {
            total += alpha;
        } else {
            return 0.0;
        }
    }
    return total;
}

float calcul_mesure_coins_mult(SOMMET_COIN* sommet, float alpha) {
    float calcul[3];
    for (int i = 0; i < 3; i++) {
        if (sommet->liaisons_communs[i] == 1) {
            calcul[i] = 1.0;
        } else if (sommet->liaisons_communs[i] == 2) {
            calcul[i] = alpha;
        } else {
            return 0.0;
        }
    }
    return (calcul[0] * calcul[1] * calcul[2]);
}

float lire_fichier_coins(char *nom_fichier, float alpha, char *arg2, char *arg3) {
    FILE *fichier = fopen(nom_fichier, "r");
    if (fichier == NULL) {
        perror("Erreur d'ouverture du fichier");
        return 0.0;
    }

    int nb_sommets, nb_aretes;
    if (fscanf(fichier, "%d_%d\n", &nb_sommets, &nb_aretes) != 2) {
        fprintf(stderr, "Erreur de lecture du nombre de sommets et d'arêtes\n");
        fclose(fichier);
        return 0.0;
    }

    float *liste_mesure = malloc(nb_sommets * sizeof(float));
    if (liste_mesure == NULL) {
        perror("Erreur d'allocation de mémoire");
        fclose(fichier);
        return 0.0;
    }

    float cagitude = 0.0;

    // Lire les sommets
    for (int i = 0; i < nb_sommets; i++) {
        char line[256];
        if (fgets(line, sizeof(line), fichier) != NULL) {
            SOMMET_COIN* sommet = (SOMMET_COIN*)malloc(sizeof(SOMMET_COIN));
            if (sommet == NULL) {
                perror("Erreur d'allocation de mémoire pour le sommet");
                free(liste_mesure);
                fclose(fichier);
                return 0.0;
            }
            init_sommet(sommet, i);
            if (sscanf(line, "%d_%d,%d,%d_%d,%d,%d_%d,%d,%d\n",
                       &sommet->id,
                       &sommet->poids[0], &sommet->poids[1], &sommet->poids[2],
                       &sommet->ids[0], &sommet->ids[1], &sommet->ids[2],
                       &sommet->liaisons_communs[0], &sommet->liaisons_communs[1], &sommet->liaisons_communs[2]) != 10) {
                fprintf(stderr, "Erreur de lecture des données du sommet\n");
                free_sommet(sommet);
                free(sommet);
                free(liste_mesure);
                fclose(fichier);
                return 0.0;
            }
            if (strcmp(arg2, "mult") == 0) {
                liste_mesure[sommet->id] = calcul_mesure_coins_mult(sommet, alpha);
            } else {
                liste_mesure[sommet->id] = calcul_mesure_coins_add(sommet, alpha);
            }
            cagitude += liste_mesure[sommet->id];
            free_sommet(sommet);
            free(sommet);
        }
    }

    if (strcmp(arg3, "connexe") == 0) {
        cagitude = 0.0;
        for (int i = 0; i < nb_aretes; i++) {
            char line[256];
            if (fgets(line, sizeof(line), fichier) != NULL) {
                ARETE_COIN arete;
                if (sscanf(line, "%d_%d_%d\n", &arete.id1, &arete.id2, &arete.poids) != 3) {
                    fprintf(stderr, "Erreur de lecture des données de l'arête\n");
                    free(liste_mesure);
                    fclose(fichier);
                    return 0.0;
                }
                cagitude += (liste_mesure[arete.id1] * liste_mesure[arete.id2]);
            }
        }
    }

    free(liste_mesure);
    fclose(fichier);
    return cagitude;
}


int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <directory> [<mode>] [<connexe>]\n", argv[0]);
        return EXIT_FAILURE;
    }
    char file[256];
    char* arg3 = "NONE";
    if (argc == 4) {
        arg3 = argv[3];
        snprintf(file, sizeof(file), "data/%s/results/liste_mesure_alpha_connexe.csv", argv[1]);
    }
    else{
        snprintf(file, sizeof(file), "data/%s/results/liste_mesure_alpha.csv", argv[1]);
    }

    char path[256];
    snprintf(path, sizeof(path), "data/%s/graphes_coins", argv[1]);
    DIR *rep = opendir(path);
    if (rep == NULL) {
        perror("opendir");
        return EXIT_FAILURE;
    }

    struct dirent *lecture;
    float alpha = 4.0;
    float temp = 0.0;
    

    FILE *f_mesure = fopen(file, "w");
    if (f_mesure == NULL) {
        printf("Impossible d'ouvrir le fichier %s\n", file);
        closedir(rep);
        return EXIT_FAILURE;
    }

    while ((lecture = readdir(rep)) != NULL) {
        if (strstr(lecture->d_name, ".csv")) {
            char file_path[512];
            snprintf(file_path, sizeof(file_path), "%s/%s", path, lecture->d_name);
            temp = (lire_fichier_coins(file_path, alpha, argv[2], arg3));
            if(temp !=0){
                temp =temp/3;
            }
            char nom_fichier_sans_extension[256];
            strcpy(nom_fichier_sans_extension, lecture->d_name);
            char *ptr = strrchr(nom_fichier_sans_extension, '.');
            if (ptr != NULL) {
                *ptr = '\0';
            }
            fprintf(f_mesure, "%s,%f\n", nom_fichier_sans_extension, temp);
        }
    }

    closedir(rep);
    fclose(f_mesure);
    return EXIT_SUCCESS;
}
