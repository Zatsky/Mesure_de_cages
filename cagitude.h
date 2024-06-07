#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "analyse_cage.h"

float calcul_mesure_coins_add(SOMMET_COIN* sommet,float alpha);
float calcul_mesure_coins_mult(SOMMET_COIN* sommet,float alpha);
float lire_fichier_coins(char *nom_fichier,float alpha, char* arg2,char* arg3);