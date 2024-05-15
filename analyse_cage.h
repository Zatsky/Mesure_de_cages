#include "utils_cage_moleculaire.h"
#include "graphe_cycles.h"
#include "lecture_molecule_sdf.h"
#include <dirent.h>
#include <stdio.h>
#include <math.h>

int nb_arete_base;
int taille_base;
ARETE *base_aretes;
int *arete_cycle;
int **arete_liste;
cycles *labase;

char *atom_name[NB_ATOM_NAMES];
#define NB_TAB 200

