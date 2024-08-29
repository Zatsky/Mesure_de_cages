#include "numerotation.h"
#include "affichage.h"
#include "Chemins.h"

/**
 * @file numerotation.c
 * @brief fichier permettant de gérer les numérotations des sommets graphes.
 */


void printNumerotation() {
  printf("\nfichier numerotation.c\n");
}

/**
 * Ajoute +1 à un tableau d'inversion donné en argument. (passe au tableau d'inversion suivant)
 * 
 * @param tabInversion : tableau d'inversion que l'on veut incrémenter
 * @param tailleTab : taille du tableau d'inversion
*/
void tabInversionPlusUn(int *tabInversion, int tailleTab) {

  int *tabRes;
  tabRes = malloc(sizeof(int) * tailleTab);
  
  for(int i = 0; i < tailleTab; i++) {
    tabRes[i] = tabInversion[i];
  }

  tabRes[0] = (tabInversion[0] + 1) % tailleTab;

  for(int i = 1; i < tailleTab; i++) {
    if((tabRes[i-1] == 0) && (tabInversion[i-1] == (tailleTab-i))) {
      tabRes[i] = (tabInversion[i] + 1) % (tailleTab-i);
    }
  }

  for(int i = 0; i < tailleTab; i++) {
    tabInversion[i] = tabRes[i];
  }

  free(tabRes);
}

/**
 * Prend un tableau de numérotation en argument, et renvoie le tableau d'inversion correspondant.
 * 
 * @param tabNum : tableau de numérotation.
 * @param tailleTab : taille du tableau
 * 
 * @return tableau d'inversion correspondant
*/
int * numToInversion(int *tabNum, int tailleTab) {

  int *tabInverse;
  tabInverse = malloc(sizeof(int) * tailleTab);

  for(int i = 0; i < tailleTab; i++) {
    int cpt = 0;
    for(int j = i; j < tailleTab; j++) {
      if(tabNum[j] < tabNum[i]) {
        cpt++;
      }
    }

    tabInverse[i] = cpt;
  }

  return tabInverse;
}

/**
 * Prend un tableau d'inversion en argument, et renvoie le tableau de numérotation correspondant.
 * 
 * @param tabNum : tableau d'inversion.
 * @param tailleTab : taille du tableau
 * 
 * @return tableau de numérotation correspondant
*/
int * inverseToNum(int *tabInverse, int tailleTab) {

  int tabSelec[tailleTab];
  for(int i = 0; i < tailleTab; i++) {
    tabSelec[i] = i;
  }  

  int *tabNum;
  tabNum = malloc(sizeof(int) * tailleTab);
  for(int i = 0; i < tailleTab; i++) {
    tabNum[i] = -1;
  }

  for(int i = 0; i < tailleTab; i++) {

    int indiceSelec = 0;

    for(int j = 0; j < tailleTab; j++) {
      if(tabSelec[j] != - 1) {
        if(indiceSelec == tabInverse[i]) {
          tabNum[i] = tabSelec[j];
          tabSelec[j] = -1;
        }
        indiceSelec++;
      }
    }
  }

  return tabNum;
}


/**
 * Prend en argument un graphe et une numérotation, et converti le graphe avec la nouvelle numérotation.
 * 
 * @param G : Graphe dont on veut modifier la numérotation.
 * @param tabNum : nouvelle numérotation à prendre en compte.
 * 
 * Step 1 : Pour les arêtes du graphe, change les indices des sommets en leurs nouvelles valeurs.
 * 
 * Step 2 : Remplace la nouvelle numérotation par la nouvelle.
 * 
 * Step 3 : Calcul la nouvelle matrice d'adjacence.
*/
void conversionNum(graphe *G, int *tabNum) {


/* --- Step 1 --- */

  for(int i = 0; i < G->nb_arete; i++) {
    for(int j = 0; j < G->nb_sommets; j++) {

      if(G->aretes[i].a1 == G->numerotation[j]) {

        G->aretes[i].a1 = tabNum[j];
        j=G->nb_sommets;
      }
    }
    for(int j = 0; j < G->nb_sommets; j++) {

      if(G->aretes[i].a2 == G->numerotation[j]) {
        G->aretes[i].a2 = tabNum[j];
        j=G->nb_sommets;
      }
    }
  }

/* --- Step 2 --- */

  for(int i = 0; i < G->nb_sommets; i++) {
    G->numerotation[i] = tabNum[i];
  }

/* --- Step 3 --- */

  for(int i = 0; i < G->nb_sommets; i++) {
    for(int j = i; j < G->nb_sommets; j++) {
      G->adjacence[i][j] = 0;
      G->adjacence[j][i] = 0;
    } 
  }

  for(int i = 0; i < G->nb_arete; i++) {
    G->adjacence[G->aretes[i].a1][G->aretes[i].a2] = 1;
    G->adjacence[G->aretes[i].a2][G->aretes[i].a1] = 1;
  }

/* --- Step 4 --- */

  for(int i = 0; i < G->nb_arete; i++) {
    if(G->aretes[i].a1 > G->aretes[i].a2) {
      int tmp;
      tmp = G->aretes[i].a1;
      G->aretes[i].a1 = G->aretes[i].a2;
      G->aretes[i].a2 = tmp;
    }
  }

  /* Re-tri les arêtes dans l'ordre lexicographique, en fonction de la nouvelle numérotation
  
  int changement = 1;
  while(changement) {

    changement = 0;

    for(int i = 1; i < G->nb_arete; i++) {

      if(G->aretes[i-1].a1 > G->aretes[i].a1) {
        ARETES a;
        a = G->aretes[i-1];
        G->aretes[i-1] = G->aretes[i];
        G->aretes[i] = a;
        changement = 1;

      } else if ((G->aretes[i-1].a1 == G->aretes[i].a1) && (G->aretes[i-1].a2 > G->aretes[i].a2)) {
        ARETES a;
        a = G->aretes[i-1];
        G->aretes[i-1] = G->aretes[i];
        G->aretes[i] = a;
        changement = 1;
      }
    }
  }*/

}

/**
 * Calcul la valeur de n! (n factoriel) d'un entier n donné en argument.
 *  
 * @param n : nombre entier
 * 
 * @return long long int contenant la valeur de n!
*/
long long int factoriel(int n) {
  long long int res = 1;

  for(int i = 0; i < n; i++) {
    res = res * (i+1);
  }

  return res;
}


/**
 * Vérifie qu'une numérotation donnée en argument est "Standard".
 * i.e : 0, 1, 2, ..., i-1, i, i+1, ..., n-1, n. 
 * 
 * @param numerotation : tableau de numérotation que l'on veut vérifier.
 * @param tailleTab : taille du tableau.
 * 
 * @return 1 si la numérotation est standard, 0 sinon.
*/
int numerotationStandard(int *numerotation, int tailleTab) {

  int vrai = 1;

  for(int i = 0; i < tailleTab; i++) {
    if(numerotation[i] != i) {
      vrai = 0;
      return vrai;
    }
  }

  return vrai;
}


/**
 * Converti un entier en un tableau d'inversion.
 * 
 * @param numInversion : entier que l'on veut convertir
 * @param tailleTab : taille du tableau que l'on veut obtenir
 * 
 * @return tableau d'inversion obtenu
 */
int * intToInversion(long long int numInversion, int tailleTab) {
  int *inv;
  inv = malloc(sizeof(int) * tailleTab);

  for(int i = 0; i < tailleTab; i++) {
    inv[i] = 0;
  }

  for(int i = 0; i < tailleTab; i++) {
    long long int val = factoriel(tailleTab) / factoriel(i+1);
    int ratio = 0;
    ratio = numInversion/val;
    inv[tailleTab - 1 - i] = ratio;
    numInversion = numInversion - (ratio * val);
  }

  return inv;
}


/**
 * Converti un tableau d'inversion en une valeur entière.
 * 
 * @param inv : tableau d'inversion que l'on veut convertir.
 * @param tailleTab : taille du tableau
 * 
 * @return valeur entière correspondante au tableau d'inversion.
 */
long long int inversionToInt(int *inv, int tailleTab) {

  long long int res = 0;

  for(int i = 0; i < tailleTab; i++) {
    res = res + inv[i] * (factoriel(tailleTab)/factoriel(tailleTab - i));
  }

  return res;
}


/**
 * Calcul la distance des tableau d'inversion de deux numérotation.
 * 
 * @param numUn : premier tableau de numérotation.
 * @param numDeux : deuxième tableau de numérotation.
 * @param tailleTab : taille des deux tableaux de numérotation.
 * 
 * @return distance en valeur entière des tableau d'inversion .
 */
int distanceNumerotation(int *numUn, int *numDeux, int tailleTab) {

  int *tabInvSuppUn;
  tabInvSuppUn = malloc(sizeof(int) * tailleTab);

  
  for(int i = tailleTab - 1; i >= 0; i--) {
    tabInvSuppUn[i] = 0;
    int nbSupp = 0;
    for(int j = 0; j < i; j++) {
      if(numUn[j] > numUn[i]) {
        nbSupp++;
      }
    }
    tabInvSuppUn[i] = nbSupp;
  }


  int *tabInvSuppDeux;
  tabInvSuppDeux = malloc(sizeof(int) * tailleTab);

  
  for(int i = tailleTab - 1; i >= 0; i--) {
    tabInvSuppDeux[i] = 0;
    int nbSupp = 0;
    for(int j = 0; j < i; j++) {
      if(numDeux[j] > numDeux[i]) {
        nbSupp++;
      }
    }
    tabInvSuppDeux[i] = nbSupp;
  }


  int distance = 0;

  for(int i = 0; i < tailleTab; i++) {

    int diff = tabInvSuppUn[i] - tabInvSuppDeux[i];
    if(diff < 0){
      diff = -diff;
    }
    distance = distance + diff;
  }

  free(tabInvSuppUn);
  free(tabInvSuppDeux);

  return distance;
}


/**
 * Calcul la différence des suites de degrés entre les numérotations de deux graphes. i.e : somme des différences de degrés des sommets ayant le même numéros. 
 * 
 * @param G : premier graphe.
 * @param H : deuxième graphe.
 * 
 */
int distanceDeg(graphe *G, graphe *H) {

  int dist = 0;

  for(int i = 0; i < G->nb_sommets; i++) {
    for(int j = 0; j < H->nb_sommets; j++) {

      if(G->numerotation[i] == H->numerotation[j]) {

        int diffDeg = 0;
        diffDeg = degSom(G, G->numerotation[i]) - degSom(H, H->numerotation[j]);

        if(diffDeg < 0) {
          diffDeg = -diffDeg;
        }

        dist += diffDeg; 
      }
    }
  }

  return dist;
}


/**
 * Calcul la différence des suites d'excentricités entre les numérotations de deux graphes. i.e : somme des différences d'excentricités des sommets ayant le même numéros. 
 * 
 * @param G : premier graphe.
 * @param H : deuxième graphe.
 * @param DAGG : tableau des DAGs de G.
 * @param DAGH : tableau des DAGs de H.
 * 
 */
int distanceEx(graphe *G, graphe *H, graphe **DAGG, graphe **DAGH) {


  int dist = 0;

  for(int i = 0; i < G->nb_sommets; i++) {
    for(int j = 0; j < H->nb_sommets; j++) {

      if(G->numerotation[i] == H->numerotation[j]) {
        int diffEx = 0;

        int exG = ExSom(G, G->numerotation[i], DAGG);
        int exH = ExSom(H, H->numerotation[j], DAGH);

        //printf("\nexcentricité  %d : %d\texcentricité  %d : %d\n",G->numerotation[i],exG, H->numerotation[j], exH);

        diffEx = exG - exH;

        if(diffEx < 0) {
          diffEx = -diffEx;
        }

        dist += diffEx; 
      }

    }
  }


  return dist;
}


/**
 * Calcul l'excentricité d'un sommet.
 * 
 * @param G : graphe G possédant le sommet en question.
 * @param som : indice du sommet en question.
 * @param DAG : tableau des DAGs du graphe G.
 */
int ExSom(graphe *G, int som, graphe **DAG) {
  

  int *dist;
  dist = malloc(sizeof(int) *G->nb_sommets);

  for(int i = 0; i < G->nb_sommets; i++) {
    dist[i] = distanceSoms(G, som, i, DAG);
  }

  int Ex = 0;

  for(int i = 0; i < G->nb_sommets; i++) {
    if(dist[i] > Ex) {
      Ex = dist[i];
    }
  }

  free(dist);

  return Ex;
}

/**
 * Calcul le degré d'un sommet.
 * 
 * @param G : graphe G possédant le sommet en question.
 * @param som : indice du sommet en question.
 */
int degSom(graphe *G, int som) {

  int deg = 0;

  for(int i = 0; i < G->nb_sommets; i++) {
    if((G->adjacence[som][i] == 1)) {
      deg++;
    }
  }

  return deg;
}

/**
 * @deprecated Fonction désuette (compliqué de définir les classes d'équivalence des sommets d'un graphe. Et surtout de les comparer, lorsque les graphes sont différents).
 * 
 * Calcul la différence des suites de classes d'équivalence entre les numérotations de deux graphes. i.e : somme des différences de classes d'équivalences des sommets ayant le même numéros. 
 * 
 * @param numUn : numérotation des sommets du premier graphe.
 * @param numDeux : numérotation des sommets du deuxième graphe.
 * @param tailleTab : taille des tableaux de numérotation
 * 
 */
int distanceClasseEquivalence(int *numUn, int *numDeux, int tailleTab, graphe *G) {

  int dist = 0;

  for(int i = 0; i < tailleTab; i++) {
    for(int j = 0; j < tailleTab; j++) {

      if(numUn[i] == numDeux[j]) {

        if(classeSommet(i, G) != classeSommet(j, G)) {
          dist++;
        }

      }
    }
  }

  return dist;
}


/**
 * @deprecated Fonction désuette (compliqué de définir les classes d'équivalence des sommets d'un graphe. Et surtout de les comparer, lorsque les graphes sont différents).
 * 
 * Calcul la classe d'équivalence d'un sommet d'un graphe.
 * 
 * @param pos : indice du sommet dont on veut la classe d'équivalence.
 * @param G : graphe dans lequel se trouve le sommet d'intérêt.
 * 
 */
int classeSommet(int pos, graphe *G) {

  return 0;
}

