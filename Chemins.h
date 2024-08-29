#ifndef CHEMINS_H
#define CHEMINS_H

#include <stdio.h>
#include <stdlib.h>

#include "structure.h"

/**
 * @file Chemins.h
 * @brief interface relative aux Chemins.
 */

void printChemins();

graphe ** GenererDAG(graphe *G);
int * plusCourtCheminDexAy(graphe **DAG, int IndiceX, int IndiceY, int tailleGraphe);
int distanceSoms(graphe *G, int som1, int som2,graphe **DAG);

#endif