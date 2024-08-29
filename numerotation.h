#ifndef NUMEROTATION_H
#define NUMEROTATION_H

#include <stdio.h>
#include <stdlib.h>

/**
 * @file numerotation.h
 * @brief interface relative à la gestion d'un numérotation d'un graphe.
 */

#include "structure.h"

void printNumerotation();

void tabInversionPlusUn(int *tabInversion, int tailleTab);
int * numToInversion(int *tabNum, int tailleTab);
int * inverseToNum(int *tabInverse, int tailleTab);
void conversionNum(graphe *G, int *tabNum);
long long int factoriel(int n);
int numerotationStandard(int *numerotation, int tailleTab);
int * intToInversion(long long int numInversion, int tailleTab);
long long int inversionToInt(int *inv, int tailleTab);
int distanceNumerotation(int *numUn, int *numDeux, int tailleTab);

int distanceDeg(graphe *G, graphe *H);
int degSom(graphe *G, int som);
int distanceClasseEquivalence(int *numUn, int *numDeux, int tailleTab, graphe *G);
int classeSommet(int pos, graphe *G);
int distanceEx(graphe *G, graphe *H, graphe **DAGG, graphe **DAGH);
int ExSom(graphe *G, int som, graphe **DAG);
#endif