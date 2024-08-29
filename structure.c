#include "structure.h"
#include "affichage.h"
#include "graphe_cycles.h"
#include "similarite.h"


/**
 * @file structure.c
 * @brief fichier permettant de gérer les différentes structure.
 */


void printStructure() {
  printf("\nfichier structure.c\n");
}

/**
 * Prend un nombre de sommets en argument, et renvoie un graphe intialiser avec ce nombre de sommets :
 * 
 * @param tailleGraphe : nombre de sommet du graphe que l'on veut initialiser.
 * 
 * @return Structure Graphe
 * 
 * NB : 
 *      - ce graphe n'a pas d'arête (le nombre maximale d'arête est alloué).
 *      - sa matrice d'adjacence est vide.
 *      - sa numérotation est  :  0, 1, 2, ..., i-1, i, i+1, ..., n-1, n.
*/
graphe * initGraphe(int tailleGraphe) {

  graphe *G;
  G = malloc(sizeof(graphe));
  if (G == NULL) {
    printf("PB malloc initGraphe \n");
    exit(0);
  }

  G->nb_sommets = tailleGraphe;
  //printf("Taille du graphe = %d\n",G->nb_sommets);

  G->adjacence = malloc(sizeof(int *) * G->nb_sommets);

  for(int i = 0; i < G->nb_sommets; i++) {
    G->adjacence[i] = malloc(sizeof(int) * G->nb_sommets);
  }

  for(int i = 0; i < G->nb_sommets; i++) {
    for(int j = 0; j < G->nb_sommets; j++) {
      if (&G->adjacence[0][0] == NULL) {
        printf("PB malloc initGraphe indice i=%d ; j=%d\n",i,j);
        exit(0);
      }
    }
  }

  for(int i = 0; i < G->nb_sommets; i++) {
    for(int j = 0; j < G->nb_sommets; j++) {
      G->adjacence[i][j] = 0;
    }
  }

  G->nb_arete = 0;
  G->aretes = malloc(sizeof(int *) * (G->nb_sommets*(G->nb_sommets - 1))/2);

  G->numerotation = malloc(sizeof(int) * G->nb_sommets);

  for(int i = 0; i < G->nb_sommets; i++) {
    G->numerotation[i] = i;
  }

  return G;
}

/**
 * Libère la mémoire allouée à un graphe.
 * 
 * @param G : graphe dont on veut libérer la mémoire.
*/
void freeGraphe(graphe *G) {
  for(int i = 0; i <G->nb_sommets; i++) {
    free(G->adjacence[i]);
  }
  free(G->adjacence);
  free(G->aretes);
  free(G->numerotation);
  free(G);
}

/**
 * Libère la mémoire allouée à un cycle.
 * 
 * @param C : cycle dont on veut libérer la mémoire.
*/
void freeCycle(cycle *C) {
  free(C->aretes);
  free(C->sommets);
  free(C);
}

/**
 * Libère la mémoire allouée à un graphe de cycle.
 * 
 * @param GC : graphe de cycle dont on veut libérer la mémoire.
*/
void freeGrapheCycle(GRAPHE_CYCLE *GC) {
  for(int i = 0; i < GC->nb_sommets; i++) {
    free(GC->adjacence[i]);
  }
  free(GC->adjacence);
  free(GC->liste_aretes);
  free(GC);
}

/**
 * Libère la mémoire allouée à un graphe produit modulaire.
 * 
 * @param Gp : Graphe produit modulaire dont on veut libérer la mémoire.
*/
void freeGrapheProd(GRAPHEPROD *Gp) {
  for(int i = 0; i <Gp->nbSommets; i++) {
    free(Gp->adjacence[i]);
  }
  free(Gp->adjacence);
  free(Gp->Sommets);
  free(Gp);
}


/**
 * Initialise une liste chainée de float avec une valeur donnée en argument.
 * 
 * @param valeur : valeur d'initialisation de la liste.
 * 
 * @return structure liste
*/
LIST *initList(float valeur) {
  LIST *l;
  l = malloc(sizeof(LIST));
  l->moySim = valeur;
  l->occurence = 1;
  l->distDeg = 0;
  l->distEx = 0;
  l->ecartType = 0.0;
  l->suivant = NULL;
  return l;
}

LISTGC *initListGC() {
  LISTGC *l;
  l = malloc(sizeof(LISTGC));
  l->GC = malloc(sizeof(GRAPHE_CYCLE));
  l->GC->nb_sommets = 0;
  l->GC->nb_aretes = 0;
  l->GC->indice = -1;
  l->GC->adjacence = NULL;
  l->GC->liste_aretes = NULL;
  l->GC->Cycles = NULL;
  l->suivant = NULL;
  return l;
}

/**
 * Libère la mémoire allouée à une liste chainée.
 * 
 * @param liste : liste dont on veut libérer la mémoire.
*/
void freeList(LIST *liste) {
    LIST *l = liste;
    while (l != NULL) {
        LIST *suivant = l->suivant;
        free(l);
        l = suivant;
    }
}

void freeListSim(LISTSIM *liste) {
    LISTSIM *l = liste;
    while (l != NULL) {
        LISTSIM *suivant = l->suivant;
        free(l);
        l = suivant;
    }
}

int insererGrapheCycles(LISTGC *liste, GRAPHE_CYCLE *GC, int nbGC, int m) {

  LISTGC *l = liste;

  while(l->suivant != NULL) {
    l = l->suivant;
    if(grapheCyclesIdentique(GC, l->GC, m)) {
      GC->indice = l->GC->indice;
      return nbGC;
    }
  }
 
  LISTGC *newGC = malloc(sizeof(LISTGC));
  GC->indice = nbGC;
  l->suivant = newGC;
  //newGC->GC = GC;
  newGC->GC = malloc(sizeof(GRAPHE_CYCLE));
  copierGrapheCycles(GC, newGC->GC, m);
  newGC->suivant = NULL;

  //printf("\nNouveau graphe de cycles trouvé ! %d\n",nbGC+1);
  //printGrapheCycle(GC);
  return (nbGC+1);
}

void copierGrapheCycles(GRAPHE_CYCLE *copy, GRAPHE_CYCLE *paste, int m) {
  paste->nb_aretes = copy->nb_aretes;
  paste->nb_sommets = copy->nb_sommets;
  paste->indice = copy->indice;

  paste->adjacence = malloc(sizeof(int *) * paste->nb_sommets);

  for(int i = 0; i < paste->nb_sommets; i++) {
    paste->adjacence[i] = malloc(sizeof(int) * paste->nb_sommets);
    for(int j = 0; j < paste->nb_sommets; j++) {
      paste->adjacence[i][j] = copy->adjacence[i][j];
    }
  }

  paste->liste_aretes = malloc(sizeof(couple) * paste->nb_aretes);
  for(int i = 0; i < paste->nb_aretes; i++) {
    paste->liste_aretes[i].id1 = copy->liste_aretes[i].id1;
    paste->liste_aretes[i].id2 = copy->liste_aretes[i].id2;
    paste->liste_aretes[i].poids = copy->liste_aretes[i].poids;
  }

  paste->Cycles = malloc(sizeof(cycle *) * paste->nb_sommets);
  for(int i = 0; i < paste->nb_sommets; i++) {
    paste->Cycles[i] = malloc(sizeof(cycle));
    copierCycles(copy->Cycles[i], paste->Cycles[i], m);
  }

  
}

void copierCycles(cycle *copy, cycle *paste, int m) {

  paste->taille = copy->taille;

  paste->aretes = malloc(sizeof(int)*m);
  for(int i = 0; i < m; i++) {
    paste->aretes[i] = copy->aretes[i];
  }

  paste->sommets = malloc(sizeof(int) *copy->taille);
  for(int i = 0; i < paste->taille; i++) {
    paste->sommets[i] = copy->sommets[i];
  }
}


/**
 * Insère une valeur flottante dans l'ordre croissant dans une liste chainée.
 * 
 * @param liste : liste dans laquelle on veut ajouter la valeur.
 * @param valeur : valeur que l'on souhaite ajouter à la liste.
*/
void insererValeur(LIST *liste, float valeur) {
  LIST *newVal = initList(valeur);
  
  // Si la liste est vide ou si la valeur est inférieure à celle du premier nœud
  if (liste == NULL || valeur < liste->moySim) {
    newVal->suivant = liste;
    liste = newVal;
    //return;
  } else {
    LIST *l = liste;
    // Trouver le bon emplacement pour insérer la nouvelle valeur
    while (l->suivant != NULL && l->suivant->moySim <= valeur) {
      
      l = l->suivant;
      if(l->moySim == valeur) {
        l->occurence++;
        freeList(newVal);
        return;
      }
      
    }
    
    newVal->suivant = l->suivant;
    l->suivant = newVal;
  }
}

void insererValeurStats(LIST *liste, float valeur) {
  LIST *newVal = initList(valeur);
  
  // Si la liste est vide
  if (liste == NULL) {
    newVal->suivant = liste;
    liste = newVal;
    //return;
  } else {
    LIST *l = liste;
    // Trouver le bon emplacement pour insérer la nouvelle valeur
    while (l->suivant != NULL) {
      
      l = l->suivant;      
    }
    
    newVal->suivant = l->suivant;
    l->suivant = newVal;
  }
}

/**
 * Insère une valeur flottante dans l'ordre croissant dans une liste chainée, si la valeur est déjà présente, calcule une moyenne et un écart-type.
 * 
 * @param liste : liste dans laquelle on veut ajouter la valeur.
 * @param valeur : valeur que l'on souhaite ajouter à la liste.
 * @param distDeg : distance des degrés des numérotations qui ont données la valeur.
 * @param distEx : distance des excentricités qui ont données la valeur.
 * @param exist : 1 si le couple distDeg/distEx existe, 0 sinon.
 * @param EcartType : valeur de l'écart-type à enregistrer, si le couple distDeg/distEx n'existe pas.
*/
void insererValeur2(LIST *liste, float valeur, int distDeg, int distEx, int exist, float EcartType) {
  LIST *newVal = initList(valeur);
  newVal->distDeg = distDeg;
  newVal->distEx = distEx;

  if(exist == 0) {
    newVal->moySim = valeur;
    newVal->ecartType = EcartType;
    newVal->occurence = -1.0;
  }
  
  // Si la liste est vide ou si la valeur est inférieure à celle du premier nœud
  if (liste == NULL || (newVal->distDeg < liste->distDeg && newVal->distEx < liste->distEx)) {
    newVal->suivant = liste;
    liste = newVal;
    //return;
  } else {
    LIST *l = liste;
    // Trouver le bon emplacement pour insérer la nouvelle valeur

     while ((l->suivant != NULL) && (l->suivant->distDeg < distDeg)) {
      l = l->suivant;     
    }

    while((l->suivant != NULL) && (l->suivant->distDeg == distDeg) && (l->suivant->distEx <= distEx)) {

            l = l->suivant;        

          if(l->distEx == distEx && l->distDeg == distDeg) {

            l->occurence++;
            l->moySim = l->moySim + ((valeur - l->moySim)/(float)l->occurence);

            float distMoy = l->moySim - valeur;
            if(distMoy < 0) distMoy = -distMoy;

            l->ecartType = l->ecartType + ((distMoy - l->ecartType)/(float)l->occurence);

            freeList(newVal);
            return;
          } 
  
      } 
    
    newVal->suivant = l->suivant;
    l->suivant = newVal;
  }
}

LISTSIM *initListSim(float valeur, int indiceGC1, int indiceGC2) {
  LISTSIM *l;
  l = malloc(sizeof(LISTSIM));
  l->sim = valeur;
  l->indiceGC1 = indiceGC1;
  l->indiceGC2 = indiceGC2;
  l->suivant = NULL;
  return l;
}

/*
float insererSim(LISTSIM *liste, GRAPHE_CYCLE *GC1, GRAPHE_CYCLE *GC2) {

  LISTSIM *newSim = initListSim(0.0, GC1->indice, GC2->indice);

  float res = 0.0;
  //printf("\nFLAG\n");
  
  // Si la liste est vide ou si la valeur est inférieure à celle du premier nœud
  if (liste == NULL) {

    newSim->suivant = liste;
    liste = newSim;
    newSim->sim = MCIS(GC1, GC2);
    res = newSim->sim;

  } else {

    if ((newSim->indiceGC1 < liste->indiceGC1 && newSim->indiceGC2 < liste->indiceGC2)) {

      newSim->suivant = liste;
      liste = newSim;
      newSim->sim = MCIS(GC1, GC2);
      res = newSim->sim;

    } else {

      LISTSIM *l = liste;
      // Trouver le bon emplacement pour insérer la nouvelle valeur

      while ((l->suivant != NULL) && (l->suivant->indiceGC1 < GC1->indice)) {
        l = l->suivant;     

      }

      while((l->suivant != NULL) && (l->suivant->indiceGC1 == GC1->indice) && (l->suivant->indiceGC2 <= GC2->indice)) {

              l = l->suivant;        

            if(l->indiceGC2 == GC2->indice && l->indiceGC1 == GC1->indice) {
                
              freeListSim(newSim);
              res = l->sim;
              return res;
            } 
    
        } 

      float val = MCIS(GC1, GC2);

      newSim->sim = val;

      newSim->suivant = l->suivant;
      l->suivant = newSim;
      res = newSim->sim;
    }
  }

  return res;
}
*/

/**
 * Calcul la moyenne des valeurs d'une liste chainée.
 * 
 * @param liste : liste dont on veut calculer la moyenne.
 * 
 * @return moyenne calculée, en valeur flottante, de la liste.
 */
float moyenneList(LIST *liste) {

  float mean = 0.0;
  long long int sum = 0;
  LIST *l = liste;

  while (l != NULL) {

    mean += l->occurence * l->moySim;
    sum += l->occurence;

    l = l->suivant;
  }

  mean = mean/(float)sum;

  return mean;
}

/**
 * Calcul la médiane des valeurs d'une liste chainée.
 * 
 * @param liste : liste dont on veut calculer la médiane.
 * 
 * @return médiane calculée, en valeur flottante, de la liste.
 */
float medianeList(LIST *liste) {

  long long int sum = 0;
  LIST *l = liste;

  while (l != NULL) {

    sum += l->occurence;
    l = l->suivant;
  }

  LIST *lprime = liste;

  long long int sumCourante = 0;


  while (lprime != NULL) {

    sumCourante += lprime->occurence;
    if(sumCourante > sum/2) {
      return (float)lprime->moySim;
    }
    lprime = lprime->suivant;
  }

  return 0.0;
}


/**
 * Prends un tableau de flottant en argument et le tri dans l'ordre croissant via l'algorithme de tri fusion (merged sort).
 * 
 * @param valeurs : tableau de valeurs flottantes.
 * @param tailleTab : taille du tableau.
 * 
*/
void triFusionFloat(float *valeurs, long long int tailleTab) {

  if(tailleTab <= 1) {
    return;
  }

  int tailleA = tailleTab/2;
  int tailleB = tailleTab/2 + tailleTab%2;

  float *valeursA;
  valeursA = malloc(sizeof(float) * tailleA);
  
  float *valeursB;
  valeursB = malloc(sizeof(float) * tailleB);

  for(int i = 0; i < tailleTab; i++) {
    if(i < tailleA) {
      valeursA[i] = valeurs[i];
    } else {
      valeursB[i - tailleA] = valeurs[i];
    }
  }

  triFusionFloat(valeursA, tailleA);
  triFusionFloat(valeursB, tailleB);

  float *fusion = fusionFloat(valeursA, tailleA, valeursB, tailleB);

  for(int i = 0; i < tailleTab; i++) {
    valeurs[i] = fusion[i];
  }

  free(valeursA);
  free(valeursB);
  free(fusion);

}

/**
 * Prends deux tableau de valeurs flottantes déjà trié et les fusionne en gardant l'ordre croissant effectif
 * 
 * @param tabA : premier tableau.
 * @param tailleA : taille du premier tableau.
 * @param tabB : deuxième tableau.
 * @param tailleB : taille du deuxième tableau.
 * 
*/
float * fusionFloat(float *tabA, long long int tailleA, float *tabB, long long int tailleB) {

  int tailleTot = tailleA + tailleB;
  float *Fusion;
  Fusion = malloc(sizeof(float) * tailleTot);

  int indiceA = 0;
  int indiceB = 0;
  for(int i = 0; i < tailleTot; i++) {
    if(indiceA >= tailleA) {
      Fusion[i] = tabB[indiceB];
      indiceB++;
    } else if (indiceB >= tailleB){
      Fusion[i] = tabA[indiceA];
      indiceA++;
    }else if(tabA[indiceA] <= tabB[indiceB]) {
      Fusion[i] = tabA[indiceA];
      indiceA++;
    } else {
      Fusion[i] = tabB[indiceB];
      indiceB++;
    }
  }

  return Fusion;
}




