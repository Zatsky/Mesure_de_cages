#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

float calcul_mesure_coins_add(int length, char* line,float alpha);
float calcul_mesure_coins_mult(int length, char* line,float alpha);
float comparator(float* moy, FILE *f, float alpha);
void cage_ou_precage(float* moy, FILE *f, float alpha);
