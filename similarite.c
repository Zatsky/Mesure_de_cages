#include "similarite.h"
#include "affichage.h"
#include "numerotation.h"
#include "graphe_cycles.h"
#include "utils_cage_moleculaire.h"
#include "analyse_cage.h"

#define NB_MAX_CLIQUES 524288
/**
 * @file similarite.c
 * @brief fichier permettant de calculer un score de similarité entre 2 graphes.
 */

void printMCIS() {
  printf("\nfichier MCIS.c\n");
}

/**
 * Prends deux graphes de cycles en argument et renvoie le graphe produit modulaire correspondant.
 * 
 * @param GC : premier graphe de cycles.
 * @param GC1 : deuxième graphe de cycles.
 * 
 * @return graphe produit modulaire des deux graphes de cycles.
 * 
 * Step 1 : définie toutes les paires de sommets entre GC et GC1.
 * 
 * Step 2 : initialise la matrice d'adjacence à 0.
 * 
 * Step 3 : construit la matrice d'adjacence du graphe produit modulaire telle que :
 *            - pour chaque paire de sommets du graphe produit modulaire (où un sommet représente une paire entre un sommet de GC et un sommet de GC1 et où les cycles correspondant dans GC et GC1 ont la même taille) :
 *              - il y a une arête entre les deux sommets si les sommets correspondant dans GC sont reliés ET les sommets correspondant dans GC1 sont reliés.
 *              - il y a une arête entre les deux sommets si les sommets correspondant dans GC ne sont pas reliés ET les sommets correspondant dans GC1 ne sont pas reliés.
 * 
 * 
 * NB : 
 *      - le graphe produit modulaire contient un tableau représentant ses sommets (toutes les paires de sommets possibles des deux graphes de cycles i.e : produit cartésiens des sommets de GC et GC1).
 *      - le graphe produit modulaire contient la matrice d'adjacence 
 *      - le graphe produit modulaire ne contient rien d'autre que ces deux éléments.
*/
GRAPHEPROD * produitModulaire_Cycle(GRAPHE_CYCLE *GC, GRAPHE_CYCLE *GC1) {
  
  int nbSommets = 0;
  for(int i = 0; i < GC->nb_sommets; i++) {
    for(int j = 0; j < GC1->nb_sommets; j++) {
      if(GC->liste_sommets[i].poids == GC1->liste_sommets[j].poids) {
        nbSommets++;
      }
    }
  }

  //printf("\n\nTaille Graphe : %d\n\n",nbSommets);

/* --- Step 1 --- */

  couple *Sommets;
  Sommets = malloc(sizeof(couple)*nbSommets);

  int indiceSommets = 0;
  for(int i = 0; i < GC->nb_sommets; i++) {
    for(int j = 0; j < GC1->nb_sommets; j++) {
      if(GC->liste_sommets[i].poids == GC1->liste_sommets[j].poids) {
        Sommets[indiceSommets].a1 = i;
        Sommets[indiceSommets].a2 = j;
        indiceSommets++;
      }
    }
  }

/* --- Step 2 --- */

  int **prodAdjacence;
  prodAdjacence = malloc(sizeof(int *) * nbSommets);

  for(int i = 0; i < nbSommets; i++) {
    prodAdjacence[i] = malloc(sizeof(int) * nbSommets);
  }

  for(int i = 0; i < nbSommets; i++) {
    for(int j = 0; j < nbSommets; j++) {
      prodAdjacence[i][j] = 0;
    }
  }

/* --- Step 3 --- */

  int nb_sommetsEff = 0;

  for(int i = 0; i < nbSommets; i++) {
    for(int j = i+1; j < nbSommets; j++) {

      if((Sommets[i].a1 != Sommets[j].a1) && (Sommets[i].a2 != Sommets[j].a2)) {
 
          if((GC->adjacence[Sommets[i].a1][Sommets[j].a1] > 0) && (GC1->adjacence[Sommets[i].a2][Sommets[j].a2] > 0)) {
            prodAdjacence[i][j] = 2;
            prodAdjacence[j][i] = 2;
            nb_sommetsEff++;
          }

          if((GC->adjacence[Sommets[i].a1][Sommets[j].a1] == 0) && (GC1->adjacence[Sommets[i].a2][Sommets[j].a2] == 0)) {
            prodAdjacence[i][j] = 1;
            prodAdjacence[j][i] = 1;
            nb_sommetsEff++;
          }

          if((GC->adjacence[Sommets[i].a1][Sommets[j].a1] == -1) && (GC1->adjacence[Sommets[i].a2][Sommets[j].a2] == -1)) {
            prodAdjacence[i][j] = 2;
            prodAdjacence[j][i] = 2;
            nb_sommetsEff++;
          }

        
      }

    }
  }

  GRAPHEPROD *Gp;

  Gp = malloc(sizeof(GRAPHEPROD));
  //printf("%d\n\n",nbSommets);
  Gp->nbSommets = nbSommets;
  Gp->Sommets = Sommets;
  Gp->adjacence = prodAdjacence;

  //printGrapheProd(Gp);
  return Gp;  
}


GRAPHEPROD * produitModulaire_Coin(GRAPHE_COIN *GC, GRAPHE_COIN *GC1) {
  
  int nbSommets = 0;
  for(int i = 0; i < GC->nb_sommets; i++) {
    for(int j = 0; j < GC1->nb_sommets; j++) {
      if(meme_sommet(GC->liste_sommets[i], GC1->liste_sommets[j])) {
        nbSommets++;
      }
    }
  }
  //printf("\n\nTaille Graphe : %d\n\n",nbSommets);

/* --- Step 1 --- */

  couple *Sommets;
  Sommets = malloc(sizeof(couple)*nbSommets);

  int indiceSommets = 0;
  for(int i = 0; i < GC->nb_sommets; i++) {
    for(int j = 0; j < GC1->nb_sommets; j++) {
      if(meme_sommet(GC->liste_sommets[i], GC1->liste_sommets[j])) {
        Sommets[indiceSommets].a1 = i;
        Sommets[indiceSommets].a2 = j;
        indiceSommets++;
      }
    }
  }

/* --- Step 2 --- */

  int **prodAdjacence;
  prodAdjacence = malloc(sizeof(int *) * nbSommets);

  for(int i = 0; i < nbSommets; i++) {
    prodAdjacence[i] = malloc(sizeof(int) * nbSommets);
  }

  for(int i = 0; i < nbSommets; i++) {
    for(int j = 0; j < nbSommets; j++) {
      prodAdjacence[i][j] = 0;
    }
  }

/* --- Step 3 --- */

  int nb_sommetsEff = 0;

  for(int i = 0; i < nbSommets; i++) {
    for(int j = i+1; j < nbSommets; j++) {

      if((Sommets[i].a1 != Sommets[j].a1) && (Sommets[i].a2 != Sommets[j].a2)) {
 
          if((GC->adjacence[Sommets[i].a1][Sommets[j].a1] > 0) && (GC1->adjacence[Sommets[i].a2][Sommets[j].a2] > 0)) {
            prodAdjacence[i][j] = 2;
            prodAdjacence[j][i] = 2;
            nb_sommetsEff++;
          }

          if((GC->adjacence[Sommets[i].a1][Sommets[j].a1] == 0) && (GC1->adjacence[Sommets[i].a2][Sommets[j].a2] == 0)) {
            prodAdjacence[i][j] = 1;
            prodAdjacence[j][i] = 1;
            nb_sommetsEff++;
          }

          if((GC->adjacence[Sommets[i].a1][Sommets[j].a1] == -1) && (GC1->adjacence[Sommets[i].a2][Sommets[j].a2] == -1)) {
            prodAdjacence[i][j] = 2;
            prodAdjacence[j][i] = 2;
            nb_sommetsEff++;
          }

        
      }

    }
  }

  GRAPHEPROD *Gp;

  Gp = malloc(sizeof(GRAPHEPROD));

  Gp->nbSommets = nbSommets;
  Gp->Sommets = Sommets;
  Gp->adjacence = prodAdjacence;

  //printGrapheProd(Gp);
  return Gp;  
}


/**
 * Fonction récursive qui liste toutes les cliques maximales d'un graphe produit modulaire donné en argument, via l'algorithme de Bron-Kerbosch
 * 
 * @param Gp : graphe produit modulaire dont on veut récupérer les cliques.
 * @param R : tableau nécessaire à l'algorithme de Bron-Kerbosh
 * @param P : tableau nécessaire à l'algorithme de Bron-Kerbosh
 * @param X : tableau nécessaire à l'algorithme de Bron-Kerbosh
 * @param Cliques : tableau stockant toutes les cliques maximales trouvées
 * @param nbCliques : adresse d'un entier stockant le nombre de cliques trouvées.
 *  
 * 
*/
void cliqueRecursive(GRAPHEPROD *Gp, int *R, int *P, int *X, int **Cliques, int *nbCliques, int *currentMaxClique, double bronKerboschComplexite) {

  /*printf("\nnb cliques : %d\tnbSommets : %d\n",*nbCliques,Gp->nbSommets);
  for(int i = 0; i < Gp->nbSommets; i++) {
    for(int j = 0; j < Gp->nbSommets; j++) {
      printf(" %d ",Cliques[i][j]);
    }printf("\n");
  }
  */
  int Pvide = 1;

  for(int i = 0; i < Gp->nbSommets; i++) {
    if(P[i] == 1) {
      Pvide = 0;
      break;
    }
  }

  int Xvide = 1;

  for(int i = 0; i < Gp->nbSommets; i++) {
    if(X[i] == 1) {
      Xvide = 0;
      break;
    }
  }


  if((Pvide == 1) && (Xvide == 1)) {

    int nb_sommetsClique = 0;
    for(int i = 0; i < Gp->nbSommets; i++) {
      if(R[i] == 1){
        nb_sommetsClique++; 
      }
    }

    if(nb_sommetsClique >= *currentMaxClique) {

      if(nb_sommetsClique > *currentMaxClique) {

        for(int i = 0; i < NB_MAX_CLIQUES; i++) {
          for(int j = 0; j < Gp->nbSommets; j++) {
            Cliques[i][j] = 0;
          }
        }
        *nbCliques = 0;
      }

      for(int i = 0; i < Gp->nbSommets; i++) {
        Cliques[*nbCliques][i] = R[i];
      }
      *nbCliques = *nbCliques+1;
      *currentMaxClique = nb_sommetsClique;

    } else {
      return;
    }
    //if(nb_sommetsClique == 3) exit(0);

  } else {

    for(int i = 0; i < Gp->nbSommets; i++) {
      
      if(P[i] == 1) {

        int *NouveauR;
        int *NouveauP;
        NouveauR = malloc(sizeof(int) * Gp->nbSommets);
        NouveauP = malloc(sizeof(int) * Gp->nbSommets);

        for(int j = 0; j < Gp->nbSommets; j++) {
          NouveauR[j] = R[j];
          NouveauP[j] = P[j];
        }

        NouveauR[i] = 1;

        for(int j = 0; j < Gp->nbSommets; j++) {
          if(Gp->adjacence[i][j] == 0) {
            NouveauP[j] = 0;
            X[j] = 0;
          }
        }

        int borneMax = 0;

        for(int j = 0; j < Gp->nbSommets; j++) {
          if(NouveauR[j] == 1 || NouveauP[j] == 1) {
            borneMax++;
          }
        }

        if(borneMax < *currentMaxClique) {
          free(NouveauP);
          free(NouveauR);
          return;
        }

        cliqueRecursive(Gp, NouveauR, NouveauP, X, Cliques, nbCliques, currentMaxClique, bronKerboschComplexite);

        P[i] = 0;
        X[i] = 1;
        free(NouveauP);
        free(NouveauR);
      }
    }
  }

}

/**
 * Liste toutes les cliques maximales d'un graphe produit modulaire via l'algorithme de Bron-Kerbosch,
 * et renvoie celle (ou une de celle) de plus grande taille i.e : une clique maximum du graphe produit modulaire, possédant un nombre maximum d'arêtes de type 2.
 * 
 * @param Gp : graphe produit modulaire dont on veut la clique max.
 * @param tailleCliqueMaxtrouve : adresse d'un entier pour stocker la taille de la clique max que l'on trouvera.
 * 
 * @return tableau d'entier contenant les indices des sommets de la clique max.
 * 
 * Step 1 : Algorithme de Bron-Kerbosch
 * 
 * Step 2 : récupération d'une clique maximum
 * 
 * NB : 
 *      - Problème au niveau de l'allocation du tableau de Cliques (ligne 226 et 228) -> on ne sait pas combien de cliques maximales seront générées.
*/
int *   CliqueMax(GRAPHEPROD *Gp, int *tailleCliqueMaxtrouve, double bronKerboschComplexite) {


/* --- Step 1 --- */


  int **Cliques;

  Cliques = malloc(sizeof(int*) * NB_MAX_CLIQUES);

  for(int i = 0; i < NB_MAX_CLIQUES; i++) {
    Cliques[i] = malloc(sizeof(int) * Gp->nbSommets);

    for(int j = 0; j < Gp->nbSommets; j++) {

      Cliques[i][j] = 0;
    }
  }


  int nbCliques = 0;
  int currentMaxClique = 0;

  int *P;
  P = malloc(sizeof(int) * Gp->nbSommets);

  int *X;
  X = malloc(sizeof(int) * Gp->nbSommets);

  int *R;
  R = malloc(sizeof(int) * Gp->nbSommets);

  for(int i = 0; i < Gp->nbSommets; i++) {
    P[i] = 1;
    X[i] = 0;
    R[i] = 0;
  }

  cliqueRecursive(Gp, R, P, X, Cliques, &nbCliques, &currentMaxClique, bronKerboschComplexite);

  free(P);
  free(R);
  free(X);

  //printf("\nnb clique trouvée : %d\n",nbCliques);

/* --- Step 2 --- */


  int indiceCliqueMax = 0;
  int tailleCliqueMax = -10;
  int nbAreteMax = -10;

  int *nbMatchAretes;
  nbMatchAretes = malloc(sizeof(int) * nbCliques);

  for(int i = 0; i < nbCliques; i++) {
    nbMatchAretes[i] = 0;
  }

  for(int i = 0; i < nbCliques; i++) {

    for(int som1 = 0; som1 < Gp->nbSommets; som1++) {
      for(int som2 = 0; som2 < Gp->nbSommets; som2++) {

        if(Cliques[i][som1] == 1 && Cliques[i][som2] == 1) {
          if(Gp->adjacence[som1][som2] == 2) {
            nbMatchAretes[i]++;
          }
        }
      }
    }
  }

  for(int i = 0; i < nbCliques; i++) {
    int tailleClique = 0;

    for(int j = 0; j < Gp->nbSommets; j++) {
      if(Cliques[i][j] == 1) {
        tailleClique++;
      }
    }

    if (tailleClique > tailleCliqueMax) {
      tailleCliqueMax = tailleClique;
      nbAreteMax = nbMatchAretes[i];
      indiceCliqueMax = i;
    } else if ((tailleClique == tailleCliqueMax) && (nbMatchAretes[i] > nbAreteMax)) {
      tailleCliqueMax = tailleClique;
      nbAreteMax = nbMatchAretes[i];
      indiceCliqueMax = i;
    }
  }

  *tailleCliqueMaxtrouve = tailleCliqueMax;
  
  int *Res;
  Res = malloc(sizeof(int) * Gp->nbSommets);
  for(int i = 0; i < Gp->nbSommets; i++) {
    Res[i] = Cliques[indiceCliqueMax][i];
  }

  for(int i = 0; i < NB_MAX_CLIQUES; i++) {
    free(Cliques[i]);
  }
  free(Cliques);
  free(nbMatchAretes);

  return Res;

}


/**
 * Extrait le plus grand sous-graphe commun induit (MCIS) de deux graphes de cycles.
 * 
 * @param Gp : graphe produit modulaire des graphes de cycles GC et GC1.
 * @param GC : premier graphe de cycles.
 * @param GC1 : deuxième graphe de cycles.
 * @param tailleMCIS : taille du plus grand sous-graphe commun induit, correspondant à la taille de la clique max de Gp.
 * @param Clique : clique max de Gp correspondant au MCIS de GC et GC1.
 * 
 * @return nombre d'arête dans le MCIS.
*/
int extractionMCIS(GRAPHEPROD *Gp, int ** GC, int **GC1, int tailleMCIS, int *Clique) {

  int *sousGraphe1;
  int *sousGraphe2;
  sousGraphe1 = malloc(sizeof(int) * tailleMCIS);
  sousGraphe2 = malloc(sizeof(int) * tailleMCIS);

  int indiceSom = 0;

  for(int i = 0; i < Gp->nbSommets; i++) {
    
    if(Clique[i] == 1) {
  
      sousGraphe1[indiceSom] = Gp->Sommets[i].a1;
      sousGraphe2[indiceSom] = Gp->Sommets[i].a2;

      indiceSom++;
    }
  }
 
  /*printf("Pour le graphe 1 : \n\n");
  for(int i = 0; i < tailleMCIS; i++) {
    printf(" Cycle numéro %d de taille %d ;",sousGraphe1[i],GC->Cycles[sousGraphe1[i]]->taille);
  }printf("\n\n");

  printf("Pour le graphe 2 : \n\n");
  for(int i = 0; i < tailleMCIS; i++) {
    printf(" Cycle numéro %d de taille %d ;",sousGraphe2[i],GC1->Cycles[sousGraphe2[i]]->taille);
  }printf("\n\n");*/
  
  int nbArete = 0;

  for(int i = 0; i < tailleMCIS; i++) {
    for(int j = 0; j < tailleMCIS; j++) {
      if(GC[sousGraphe1[i]][sousGraphe1[j]] >0 || GC[sousGraphe1[i]][sousGraphe1[j]] == -1) {
        nbArete++;
      }
    }
  }

  free(sousGraphe1);
  free(sousGraphe2);

  //printf("\nNombre d'arêtes du MCIS : %d\n",nbArete/2);

  return nbArete/2;
}

/**
 * Calcul et renvoie un indice de similarité via la formule de Raymond, grâce au MCIS de deux graphes
 * 
 * @param nMCIS : nombre de sommets dans le MCIS.
 * @param mMCIS : nombre d'arêtes dans le MCIS.
 * @param nGraphe1 : nomdre de sommets dans le graphe 1.
 * @param mGraphe1 : nomdre d'arêtes dans le graphe 1.
 * @param nGraphe2 : nomdre de sommets dans le graphe 2.
 * @param mGraphe2 : nomdre d'arêtes dans le graphe 2.
 * 
 * @return score de similarité entre les deux graphes, en valeur flottante, via la formule de Raymond.
*/
float similarite(int nMCIS, int mMCIS, int nGraphe1, int mGraphe1, int nGraphe2, int mGraphe2) {

  float sim = 0.0;

  //printf("\n%d\t%d\n",nMCIS, mMCIS);
  //printf("%d\t%d\n",nGraphe1, mGraphe1);
  //printf("%d\t%d\n",nGraphe2, mGraphe2);

  sim = ( (float)(nMCIS + mMCIS) * (float)(nMCIS + mMCIS) ) / ( (float)(nGraphe1 + mGraphe1) * (float)(nGraphe2 + mGraphe2) );

  return sim;
}

/**
 * Prends en argument deux graphes de cycles et leur graphe produit modulaire et une valeur de similarité entre ces deux graphes de cycles.
 * 
 * @param Gp : Graphe produit modulaire des deux graphes de cycles.
 * @param GC : premier graphe de cycles.
 * @param GC1 : deuxième graphe de cycles.
 * 
 * Step 1 : recherche d'une clique maximum dans le graphe produit modulaire
 * 
 * Step 2 : extraction du MCIS des deux graphes de cycles, correspondant à la clique maximum de leur graphe produit modulaire.
 * 
 * Step 3 : calcul de la similarité
 * 
*/
float*  MCIS(char * base, char * mol1, char * mol2) {


  int nbAreteMCIS;
  int tailleCliqueMax1 = 0;
  float *sims = (float *)malloc(2 * sizeof(float));
  if (sims == NULL) {
      fprintf(stderr, "Erreur d'allocation de mémoire pour les similarités\n");
      return NULL;
  }
  

  GRAPHE_COIN gco1;
  GRAPHE_COIN gco2;
  double bronKerboschComplexite;
  int *CliquesMaximum2;
  int tailleCliqueMax2 = 0;
  char GCO1[256];
	snprintf(GCO1, sizeof(GCO1), "data/%s/graphes_coins/%s.csv",base,mol1);
  char GCO2[256];
	snprintf(GCO2, sizeof(GCO2), "data/%s/graphes_coins/%s.csv",base,mol2);

  if (graphe_coin_depuis_csv(GCO1, &gco1) != 0) {
    printf("erreur lecture 1");
  }
  if (graphe_coin_depuis_csv(GCO2, &gco2) != 0) {
    printf("erreur lecture 2");
  }

  if ((gco1.nb_sommets != gco2.nb_sommets) || (gco1.nb_aretes != gco2.nb_aretes)){
    sims[0] = 0;
    sims[1] = 0;

    liberer_graphe_coin(&gco1);
    liberer_graphe_coin(&gco2);
    return sims;
  }

  GRAPHEPROD *Gp2;
  Gp2 = produitModulaire_Coin(&gco1, &gco2);
  bronKerboschComplexite = pow(3, Gp2->nbSommets / 3.0);
  if (bronKerboschComplexite > 3*pow(10,10)){
    sims[1] = 0.0;
    freeGrapheProd(Gp2);

    liberer_graphe_coin(&gco1);
    liberer_graphe_coin(&gco2);

    return sims;
  }
/* --- Step 1 --- */
  CliquesMaximum2 = CliqueMax(Gp2, &tailleCliqueMax2, bronKerboschComplexite);


/* --- Step 2 --- */
  nbAreteMCIS = extractionMCIS(Gp2, gco1.adjacence, gco2.adjacence, tailleCliqueMax2, CliquesMaximum2);

/* --- Step 3 --- */

  sims[1] = 0.0;
  sims[1] = similarite(tailleCliqueMax2, nbAreteMCIS, gco1.nb_sommets, gco1.nb_aretes, gco2.nb_sommets, gco2.nb_aretes);

  free(CliquesMaximum2);
  freeGrapheProd(Gp2);
  
  liberer_graphe_coin(&gco1);
  liberer_graphe_coin(&gco2);

  GRAPHE_CYCLE gcy1;
  GRAPHE_CYCLE gcy2;

  char GCY1[256];
	snprintf(GCY1, sizeof(GCY1), "data/%s/graphes_cycles/%s.csv",base,mol1);
  char GCY2[256];
	snprintf(GCY2, sizeof(GCY2), "data/%s/graphes_cycles/%s.csv",base,mol2);
  
  
  if (graphe_cycle_depuis_csv(GCY1, &gcy1) != 0) {
    printf("erreur lecture 1");
  }
  if (graphe_cycle_depuis_csv(GCY2, &gcy2) != 0) {
    printf("erreur lecture 2");
  }

  int *CliquesMaximum1;

  if ((gcy1.nb_sommets != gcy2.nb_sommets) || (gcy1.nb_aretes != gcy2.nb_aretes)){
    sims[0] = 0;

    liberer_graphe_cycle(&gcy1);
    liberer_graphe_cycle(&gcy2);
    return sims;
  }
  GRAPHEPROD *Gp1;

  Gp1 = produitModulaire_Cycle(&gcy1, &gcy2);
  bronKerboschComplexite = pow(3, Gp1->nbSommets / 3.0);
  if (bronKerboschComplexite > 3*pow(10,10)){
    sims[0] = 0.0;
    sims[1] = 0.0;
    freeGrapheProd(Gp1);

    
    liberer_graphe_cycle(&gcy1);
    liberer_graphe_cycle(&gcy2);

    return sims;
  }
/* --- Step 1 --- */

  CliquesMaximum1 = CliqueMax(Gp1, &tailleCliqueMax1, bronKerboschComplexite);
  //printf("taille clique max = %d \n", tailleCliqueMax1);

/* --- Step 2 --- */

  nbAreteMCIS = extractionMCIS(Gp1, gcy1.adjacence, gcy2.adjacence, tailleCliqueMax1, CliquesMaximum1);

/* --- Step 3 --- */

  sims[0] = 0.0;
  sims[0] = similarite(tailleCliqueMax1, nbAreteMCIS, gcy1.nb_sommets, gcy1.nb_aretes, gcy2.nb_sommets, gcy2.nb_aretes);

  free(CliquesMaximum1);
  freeGrapheProd(Gp1);

  liberer_graphe_cycle(&gcy1);
  liberer_graphe_cycle(&gcy2);


  return sims;
}

float mesure_cagitude_clique(char * arg1, char * mol,double bronKerboschComplexite){
  GRAPHE_COIN gco;
  char GCO[256];
	snprintf(GCO, sizeof(GCO), "data/%s/graphes_coins/%s.csv",arg1,mol);
  if (graphe_coin_depuis_csv(GCO, &gco) != 0) {
    printf("erreur lecture 1");
  }
  int *CliquesMaximum1;
  int taille_clique_max;
  CliquesMaximum1 = CliqueMaxCoin(&gco,&taille_clique_max, bronKerboschComplexite);
  return 0.0;
}

int * CliqueMaxCoin(GRAPHE_COIN *Gco, int *tailleCliqueMaxtrouve, double bronKerboschComplexite) {


/* --- Step 1 --- */


  int **Cliques;

  Cliques = malloc(sizeof(int*) * NB_MAX_CLIQUES);

  for(int i = 0; i < NB_MAX_CLIQUES; i++) {
    Cliques[i] = malloc(sizeof(int) * Gco->nb_sommets);

    for(int j = 0; j < Gco->nb_sommets; j++) {

      Cliques[i][j] = 0;
    }
  }


  int nbCliques = 0;
  int currentMaxClique = 0;

  int *P;
  P = malloc(sizeof(int) * Gco->nb_sommets);

  int *X;
  X = malloc(sizeof(int) * Gco->nb_sommets);

  int *R;
  R = malloc(sizeof(int) * Gco->nb_sommets);

  for(int i = 0; i < Gco->nb_sommets; i++) {
    P[i] = 1;
    X[i] = 0;
    R[i] = 0;
  }

  cliqueRecursiveCoin(Gco, R, P, X, Cliques, &nbCliques, &currentMaxClique, bronKerboschComplexite);

  free(P);
  free(R);
  free(X);

  //printf("\nnb clique trouvée : %d\n",nbCliques);

/* --- Step 2 --- */


  int indiceCliqueMax = 0;
  int tailleCliqueMax = -10;
  int nbAreteMax = -10;

  int *nbMatchAretes;
  nbMatchAretes = malloc(sizeof(int) * nbCliques);

  for(int i = 0; i < nbCliques; i++) {
    nbMatchAretes[i] = 0;
  }

  for(int i = 0; i < nbCliques; i++) {

    for(int som1 = 0; som1 < Gco->nb_sommets; som1++) {
      for(int som2 = 0; som2 < Gco->nb_sommets; som2++) {

        if(Cliques[i][som1] == 1 && Cliques[i][som2] == 1) {
          if(Gco->adjacence[som1][som2] == 2) {
            nbMatchAretes[i]++;
          }
        }
      }
    }
  }

  for(int i = 0; i < nbCliques; i++) {
    int tailleClique = 0;

    for(int j = 0; j < Gco->nb_sommets; j++) {
      if(Cliques[i][j] == 1) {
        tailleClique++;
      }
    }

    if (tailleClique > tailleCliqueMax) {
      tailleCliqueMax = tailleClique;
      nbAreteMax = nbMatchAretes[i];
      indiceCliqueMax = i;
    } else if ((tailleClique == tailleCliqueMax) && (nbMatchAretes[i] > nbAreteMax)) {
      tailleCliqueMax = tailleClique;
      nbAreteMax = nbMatchAretes[i];
      indiceCliqueMax = i;
    }
  }

  *tailleCliqueMaxtrouve = tailleCliqueMax;
  
  int *Res;
  Res = malloc(sizeof(int) * Gco->nb_sommets);
  for(int i = 0; i < Gco->nb_sommets; i++) {
    Res[i] = Cliques[indiceCliqueMax][i];
  }

  for(int i = 0; i < NB_MAX_CLIQUES; i++) {
    free(Cliques[i]);
  }
  free(Cliques);
  free(nbMatchAretes);

  return Res;

}

void cliqueRecursiveCoin(GRAPHE_COIN *Gco, int *R, int *P, int *X, int **Cliques, int *nbCliques, int *currentMaxClique, double bronKerboschComplexite) {

  /*printf("\nnb cliques : %d\tnbSommets : %d\n",*nbCliques,Gp->nbSommets);
  for(int i = 0; i < Gp->nbSommets; i++) {
    for(int j = 0; j < Gp->nbSommets; j++) {
      printf(" %d ",Cliques[i][j]);
    }printf("\n");
  }
  */
  int Pvide = 1;

  for(int i = 0; i < Gco->nb_sommets; i++) {
    if(P[i] == 1) {
      Pvide = 0;
      break;
    }
  }

  int Xvide = 1;

  for(int i = 0; i < Gco->nb_sommets; i++) {
    if(X[i] == 1) {
      Xvide = 0;
      break;
    }
  }


  if((Pvide == 1) && (Xvide == 1)) {

    int nb_sommetsClique = 0;
    for(int i = 0; i < Gco->nb_sommets; i++) {
      if(R[i] == 1){
        nb_sommetsClique++; 
      }
    }

    if(nb_sommetsClique >= *currentMaxClique) {

      if(nb_sommetsClique > *currentMaxClique) {

        for(int i = 0; i < NB_MAX_CLIQUES; i++) {
          for(int j = 0; j < Gco->nb_sommets; j++) {
            Cliques[i][j] = 0;
          }
        }
        *nbCliques = 0;
      }

      for(int i = 0; i < Gco->nb_sommets; i++) {
        Cliques[*nbCliques][i] = R[i];
      }
      *nbCliques = *nbCliques+1;
      *currentMaxClique = nb_sommetsClique;

    } else {
      return;
    }
    //if(nb_sommetsClique == 3) exit(0);

  } else {

    for(int i = 0; i < Gco->nb_sommets; i++) {
      
      if(P[i] == 1) {

        int *NouveauR;
        int *NouveauP;
        NouveauR = malloc(sizeof(int) * Gco->nb_sommets);
        NouveauP = malloc(sizeof(int) * Gco->nb_sommets);

        for(int j = 0; j < Gco->nb_sommets; j++) {
          NouveauR[j] = R[j];
          NouveauP[j] = P[j];
        }

        NouveauR[i] = 1;

        for(int j = 0; j < Gco->nb_sommets; j++) {
          if(Gco->adjacence[i][j] == 0) {
            NouveauP[j] = 0;
            X[j] = 0;
          }
        }

        int borneMax = 0;

        for(int j = 0; j < Gco->nb_sommets; j++) {
          if(NouveauR[j] == 1 || NouveauP[j] == 1) {
            borneMax++;
          }
        }

        if(borneMax < *currentMaxClique) {
          free(NouveauP);
          free(NouveauR);
          return;
        }

        cliqueRecursiveCoin(Gco, NouveauR, NouveauP, X, Cliques, nbCliques, currentMaxClique,bronKerboschComplexite);

        P[i] = 0;
        X[i] = 1;
        free(NouveauP);
        free(NouveauR);
      }
    }
  }
}

int graphe_cycle_depuis_csv(const char *fichier, GRAPHE_CYCLE *graphe) {
    FILE *fp = fopen(fichier, "r");
    if (fp == NULL) {
        fprintf(stderr, "Erreur : impossible d'ouvrir le fichier %s\n", fichier);
        return -1;
    }

    // Lire le nombre de sommets et d'arêtes
    if (fscanf(fp, "%d_%d\n", &graphe->nb_sommets, &graphe->nb_aretes) != 2) {
        fprintf(stderr, "Erreur de lecture du nombre de sommets et d'arêtes\n");
        fclose(fp);
        return -1;
    }

    // Allouer de la mémoire pour les sommets et les arêtes
    graphe->liste_sommets = (SOMMET *)malloc(graphe->nb_sommets * sizeof(SOMMET));
    graphe->liste_aretes = (ARETE *)malloc(graphe->nb_aretes * sizeof(ARETE));
    if (graphe->liste_sommets == NULL || graphe->liste_aretes == NULL) {
        fprintf(stderr, "Erreur d'allocation de mémoire\n");
        fclose(fp);
        return -1;
    }

    // Lire les sommets
    for (int i = 0; i < graphe->nb_sommets; i++) {
        if (fscanf(fp, "%d_%d\n", &graphe->liste_sommets[i].id, &graphe->liste_sommets[i].poids) != 2) {
            fprintf(stderr, "Erreur de lecture des sommets\n");
            fclose(fp);
            return -1;
        }
    }

    // Lire les arêtes
    for (int i = 0; i < graphe->nb_aretes; i++) {
        if (fscanf(fp, "%d_%d_%d\n", &graphe->liste_aretes[i].id1, &graphe->liste_aretes[i].id2, &graphe->liste_aretes[i].poids) != 3) {
            fprintf(stderr, "Erreur de lecture des arêtes\n");
            fclose(fp);
            return -1;
        }
    }

    // Allouer et initialiser la matrice d'adjacence
    graphe->adjacence = (int **)malloc(graphe->nb_sommets * sizeof(int *));
    for (int i = 0; i < graphe->nb_sommets; i++) {
        graphe->adjacence[i] = (int *)malloc(graphe->nb_sommets * sizeof(int));
        for (int j = 0; j < graphe->nb_sommets; j++) {
            graphe->adjacence[i][j] = 0;
        }
    }

    // Remplir la matrice d'adjacence
    for (int i = 0; i < graphe->nb_aretes; i++) {
        int u = graphe->liste_aretes[i].id1;
        int v = graphe->liste_aretes[i].id2;
        graphe->adjacence[u][v] = graphe->liste_aretes[i].poids;
        graphe->adjacence[v][u] = graphe->liste_aretes[i].poids; // Si le graphe est non orienté
    }

    fclose(fp);
    return 0;
}

int graphe_coin_depuis_csv(const char *fichier, GRAPHE_COIN *graphe) {
    FILE *fp = fopen(fichier, "r");
    if (fp == NULL) {
        fprintf(stderr, "Erreur : impossible d'ouvrir le fichier %s\n", fichier);
        return -1;
    }

    // Lire le nombre de sommets et d'arêtes
    if (fscanf(fp, "%d_%d\n", &graphe->nb_sommets, &graphe->nb_aretes) != 2) {
        fprintf(stderr, "Erreur de lecture du nombre de sommets et d'arêtes\n");
        fclose(fp);
        return -1;
    }

    // Allouer de la mémoire pour les sommets et les arêtes
    graphe->liste_sommets = (SOMMET_COIN *)malloc(graphe->nb_sommets * sizeof(SOMMET_COIN));
    graphe->liste_aretes = (ARETE_COIN *)malloc(graphe->nb_aretes * sizeof(ARETE_COIN));
    if (graphe->liste_sommets == NULL || graphe->liste_aretes == NULL) {
        fprintf(stderr, "Erreur d'allocation de mémoire\n");
        fclose(fp);
        return -1;
    }
    graphe->liste_sommets->ids = (int *)malloc(3 * sizeof(int));
    // Lire les sommets
    for (int i = 0; i < graphe->nb_sommets; i++) {
      graphe->liste_sommets[i].ids = (int *)malloc(3 * sizeof(int));
      graphe->liste_sommets[i].poids = (int *)malloc(3 * sizeof(int));
      graphe->liste_sommets[i].liaisons_communs = (int *)malloc(3 * sizeof(int));
      if (fscanf(fp, "%d_%d,%d,%d_%d,%d,%d_%d,%d,%d\n",
        
        &graphe->liste_sommets[i].id, &graphe->liste_sommets[i].poids[0], &graphe->liste_sommets[i].poids[1], &graphe->liste_sommets[i].poids[2],
        &graphe->liste_sommets[i].ids[0], &graphe->liste_sommets[i].ids[1], &graphe->liste_sommets[i].ids[2],
        &graphe->liste_sommets[i].liaisons_communs[0], &graphe->liste_sommets[i].liaisons_communs[1], &graphe->liste_sommets[i].liaisons_communs[2]) != 10) {
            fprintf(stderr, "Erreur de lecture des sommets\n");
            fclose(fp);
            return -1;
        }
    }

    // Lire les arêtes
    for (int i = 0; i < graphe->nb_aretes; i++) {
        if (fscanf(fp, "%d_%d_%d\n", &graphe->liste_aretes[i].id1, &graphe->liste_aretes[i].id2, &graphe->liste_aretes[i].poids) != 3) {
            fprintf(stderr, "Erreur de lecture des arêtes\n");
            fclose(fp);
            return -1;
        }
    }

    // Allouer et initialiser la matrice d'adjacence
    graphe->adjacence = (int **)malloc(graphe->nb_sommets * sizeof(int *));
    for (int i = 0; i < graphe->nb_sommets; i++) {
        graphe->adjacence[i] = (int *)malloc(graphe->nb_sommets * sizeof(int));
        for (int j = 0; j < graphe->nb_sommets; j++) {
            graphe->adjacence[i][j] = 0;
        }
    }

    // Remplir la matrice d'adjacence
    for (int i = 0; i < graphe->nb_aretes; i++) {
        int u = graphe->liste_aretes[i].id1;
        int v = graphe->liste_aretes[i].id2;
        graphe->adjacence[u][v] = graphe->liste_aretes[i].poids;
        graphe->adjacence[v][u] = graphe->liste_aretes[i].poids; // Si le graphe est non orienté
    }

    fclose(fp);
    return 0;
}

void liberer_graphe_cycle(GRAPHE_CYCLE *g) {
    if (g == NULL) return;

    if (g->liste_sommets != NULL) {
        free(g->liste_sommets);
    }
    if (g->liste_aretes != NULL) {
        free(g->liste_aretes);
    }
    if (g->adjacence != NULL) {
        for (int i = 0; i < g->nb_sommets; i++) {
            if (g->adjacence[i] != NULL) {
                free(g->adjacence[i]);
            }
        }
        free(g->adjacence);
    }
}

void liberer_graphe_coin(GRAPHE_COIN* c)
{
	int i;
	if(c->liste_sommets != NULL){
		for(i=0;i<c->nb_sommets;i++)
		{
			free(c->liste_sommets[i].ids);
			free(c->liste_sommets[i].poids);
			free(c->liste_sommets[i].liaisons_communs);
		}
		free(c->liste_sommets);
	}
	if(c->liste_aretes != NULL){
		free(c->liste_aretes);
  }
	if (c->adjacence != NULL) {
        for (int i = 0; i < c->nb_sommets; i++) {
            if (c->adjacence[i] != NULL) {
                free(c->adjacence[i]);
            }
        }
        free(c->adjacence);
    }
}

void liberer_float(float * nums){
  free(nums);
}

int main(int argc, char *argv[]) {

  
  float* vals = MCIS("CHEBI","CHEBI_65526", "CHEBI_90211");
  printf("Les valeurs de MCIS sont %f %f",vals[0],vals[1]);
  //float* vals = MCIS(argv[1], argv[2], argv[3]);
  liberer_float(vals);
  return 0;
}