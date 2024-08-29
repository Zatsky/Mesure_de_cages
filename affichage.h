#ifndef AFFICHAGE_H
#define AFFICHAGE_H

#include <stdio.h>
#include <stdlib.h>

#include "structure.h"

/**
 * @file affichage.h
 * @brief interface relative aux fonctions d'affichage.
 */

void printAffichage();

void printGraphe(graphe *G);
void printDAG(graphe **DAG);
void printAretesGraphe(graphe *G);
void printAretesCycle(cycle *C, int m);
void printCycle(cycle *C, int m);
void printBase(cycle **Base, graphe *G, int tailleBase);
void printGrapheCycle(GRAPHE_CYCLE *GC);
void printGrapheProd(GRAPHEPROD *Gp);
void printList(LIST *liste, FILE *f);
void printListStats(LIST *liste, FILE *f);
void printTab(int *tab, int tailleTab, char *message);
void printChemin(int *Chemin, int IndiceX, int IndiceY);

#endif