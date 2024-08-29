#ifndef MCIS_H
#define MCIS_H

#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "structure.h"
#include "utils_cage_moleculaire.h"
/**
 * @file similarite.h
 * @brief interface relative au calcul de similarité.
 */

void printMCIS();

GRAPHEPROD * produitModulaire_Cycle(GRAPHE_CYCLE *GC, GRAPHE_CYCLE *GC1);
GRAPHEPROD * produitModulaire_Coin(GRAPHE_COIN *GC, GRAPHE_COIN *GC1);
void cliqueRecursive(GRAPHEPROD *Gp, int *R, int *P, int *X, int **Cliques, int *nbCliques, int *currentMaxClique, double bronKerboschComplexite);
int * CliqueMax(GRAPHEPROD *Gp, int *tailleCliqueMaxtrouve, double bronKerboschComplexite);
int extractionMCIS(GRAPHEPROD *Gp, int **GC, int **GC1, int tailleMCIS, int *Clique);
float similarite(int nMCIS, int mMCIS, int nGraphe1, int mGraphe1, int nGraphe2, int mGraphe2);
float * MCIS(char * base, char * mol1, char* mol2);
int distAreteGC(GRAPHE_CYCLE *GC1, GRAPHE_CYCLE *GC2);
float McSplit(GRAPHE_CYCLE *G, GRAPHE_CYCLE *H);
int graphe_cycle_depuis_csv(const char *fichier, GRAPHE_CYCLE *graphe);
int graphe_coin_depuis_csv(const char *fichier, GRAPHE_COIN *graphe);
void liberer_graphe_cycle(GRAPHE_CYCLE *g);
void liberer_graphe_coin(GRAPHE_COIN *g);
void liberer_float(float * nums);
int * CliqueMaxCoin(GRAPHE_COIN *Gco, int *tailleCliqueMaxtrouve, double bronKerboschComplexite);
void cliqueRecursiveCoin(GRAPHE_COIN *Gco, int *R, int *P, int *X, int **Cliques, int *nbCliques, int *currentMaxClique, double bronKerboschComplexite);

void numerotationMCIS(graphe *G, graphe *H);
int verifNumerotationMCIS(int *num, int taille);
int initialisation(char * base, char * mol1, char * mol2, GRAPHE_CYCLE * gcy1, GRAPHE_CYCLE * gcy2, GRAPHE_COIN* gco1, GRAPHE_COIN * gco2);

#endif
