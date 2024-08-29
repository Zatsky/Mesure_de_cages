#include "affichage.h"

/**
 * @file affichage.c
 * @brief fichier contenant toutes les fonctions d'affichage.
 */

void printAffichage() {
  printf("\nfichier affichage.c\n");
}

/**
 * Affiche un graphe dans la console.
 * 
 * @param G : graphe que l'on veut afficher.
 * 
 * NB : 
 *      - affiche sa matrice d'adjacence.
 *      - puis sa numérotation actuelle.
 *      - puis ses arêtes.
*/
void printGraphe(graphe *G) {

  printf("\n|\t\t|");
  for (int i = 0; i < G->nb_sommets; i++) {
    printf(" %d |",i);
  }
  printf("\n");

  printf("-------------|");
  for (int i = 0; i < G->nb_sommets; i++) {
    printf(" ------ |");
  }
  printf("\n");

  for(int i = 0; i < G->nb_sommets; i++) {
    printf("Sommet %d :\t|",i);

    for(int j = 0; j < G->nb_sommets; j++) {
      printf(" %d |",G->adjacence[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  printf("\nNumérotion : ");
  for(int i = 0; i < G->nb_sommets; i++) {
    printf(" %d ",G->numerotation[i]);
  }printf("\n");

  printAretesGraphe(G);
}

/**
 * Affiche l'ensemble des DAG d'un tableau de DAG.
 * 
 * @param DAG : tableau de DAG que l'on veut afficher. 
*/
void printDAG(graphe **DAG) {
for (int i = 0; i < DAG[0]->nb_sommets;i++) {
    printf("\nD_%d\n",i+1);
    printGraphe(DAG[i]);
  }
}

/**
 * Affiche les arêtes d'un graphe.
 * 
 * @param G : Graphe dont on veut afficher les arêtes.
 * 
 * NB : 
 *      - affiche les arêtes dans le même ordre qu'elles sont lues dans le fichier, indépendemment d'un changement quelconque de numérotation.
*/
void printAretesGraphe(graphe *G) {

  for(int i = 0; i < G->nb_arete; i++) {
    printf("Arête numéro %d : Sommet 1 = %d ; Sommet 2 = %d\n",i,G->aretes[i].a1,G->aretes[i].a2);
  }

}


/**
 * Affiche les arêtes d'un cycle, sous forme d'un vecteur {0,1}^m.
 * 
 * @param C : Cycle dont on veut afficher les arêtes.
 * @param m : nombre d'arête du graphe.  
 * 
 * NB : 
 *      - Affiche les arêtes dans le même ordre que celui du Graphe (ordre de lecture du fichier correspondant).
 *      - Pour la i-ème arête, vaut 0 si le cycle ne la possède pas, et 1 sinon.
*/
void printAretesCycle(cycle *C, int m) {

  printf("\nArêtes du cycle : ");
  for(int i = 0; i < m; i++) {
    printf(" %d ;",C->aretes[i]);
  }

}

/**
 * Affiche un cycle.
 * 
 * @param C : Cycle que l'on veut afficher.
 * @param m : nombre d'arête du graphe.
 * 
 * NB : 
 *      - affiche la taille du cycle
 *      - puis les sommets qui le compose.
 *      - puis les arêtes qui le compose.
*/
void printCycle(cycle *C, int m) {

  printf("Taille du cycle : %d\n",C->taille);
  printf("Sommets du cycle : ");
  for (int i = 0; i < C->taille; i++) {
    printf(" %d ;",C->sommets[i]);
  }

  printAretesCycle(C, m);
  printf("\n\n\n");
}

/**
 * Prend un base de cycle en argument (tableau de cycle) et affiche tous les cycles qui la compose.
 * 
 * @param Base : base de cycle que l'on veut afficher.
 * @param G : graphe relatif à cette base de cycle.
 * @param tailleBase : nombre de cycle dans la base de cycle.
*/
void printBase(cycle **Base, graphe *G, int tailleBase) {

  printf("\n\nAFFICHAGE DES CYCLES DE LA BASE\n\n");
  for(int i = 0; i < tailleBase; i++) {
    printf("Cycle numéro %d\n\n",i+1);
    printCycle(Base[i], G->nb_arete);
  }
  printf("\nNombre de cycles dans la base : %d\n\n",tailleBase);

}

/**
 * Affiche un graphe de  cycle.
 * 
 * @param GC : Graphe de cycle que l'on veut afficher.
 * 
 * NB : 
 *      - affiche sa matrice d'adjacence (avec des valeurs valant 0, 1 ou -1).
 *        - 0 signifie que les cycles ne partagent rien
 *        - -1 signifie que les cycles partagent au moins un sommet, et aucune arêtes.
 *        - 1 signifie que les cycles partagent au moins une arête.
 *      - affiche le nombre d'arête (valant 1 ou -1).
*/
void printGrapheCycle(GRAPHE_CYCLE *GC) {
  printf("\nAFFICHAGE DU GRAPHE DE CYCLE : \n");
  printf("\n|\t\t|");
  for (int i = 0; i < GC->nb_sommets; i++) {
    printf(" %d |",i);
  }
  printf("\n");

  printf("------------|");
  for (int i = 0; i < GC->nb_sommets; i++) {
    printf(" ------ |");
  }
  printf("\n");

  for(int i = 0; i < GC->nb_sommets; i++) {
    printf("Sommet %d :\t|",i);

    for(int j = 0; j < GC->nb_sommets; j++) {
      printf(" %d |",GC->adjacence[i][j]);
    }
    printf("\n");
  }
  printf("\nNombre d'arêtes : %d\n",GC->nb_aretes);
  printf("\n");

}


/**
 * Affiche un graphe produit modulaire.
 * 
 * @param Gp : Graphe produit modulaire que l'on veut afficher.
 * 
 * NB :
 *      - n'affiche que sa matrice d'adjacence.
*/
void printGrapheProd(GRAPHEPROD *Gp) {

  printf("\n");
  for(int i = 0; i < Gp->nbSommets; i++) {
    printf("%d paires : %d ; %d\n",i,Gp->Sommets[i].a1,Gp->Sommets[i].a2);
  }
  printf("\n");

  printf("\nMatrice d'adjacence du graphe produit Modulaire : \n");

  printf("\n\t");
  for (int i = 0; i < Gp->nbSommets; i++) {
    printf("\t| %d |",i);
  }
  printf("\n");

  printf("-------------");
  for (int i = 0; i < Gp->nbSommets; i++) {
    printf("| ------ ");
  }
  printf("\n");

  for(int i = 0; i < Gp->nbSommets; i++) {
    printf("Sommet %d :\t",i);

    for(int j = 0; j < Gp->nbSommets; j++) {
      printf("| %d |\t",Gp->adjacence[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

/**
 * Écrit une liste chainée dans un fichier.
 * 
 * @param liste : liste que l'on veut écrire.
 * @param f : fichier dans lequel on veut écrire.
*/
void printList(LIST *liste, FILE *f) {
    LIST *l = liste;
    int indice = 0;
    l = l->suivant;
    while (l != NULL) {

        fprintf(f, "%d;%d;%d;%f;%f;\n", l->distDeg, l->distEx, l->occurence, l->moySim, l->ecartType);
        l = l->suivant;
        indice++;
    }

}

void printListStats(LIST *liste, FILE *f) {
    LIST *l = liste;
    int indice = 0;
    l = l->suivant;
    while (l != NULL) {

        fprintf(f, "%f;\n",l->moySim);
        l = l->suivant;
        indice++;
    }

}


/**
 * Affiche un tableau d'entier, ainsi qu'un message relatif à ce tableau.
 * 
 * @param tab : tableau que l'on veut afficher.
 * @param tailleTab : nomre d'éléments dans le tableau.
 * @param message : message que l'on veut afficher.
*/
void printTab(int *tab, int tailleTab, char *message) {
  
  printf("\n%s\n",message);

  for(int i = 0; i < tailleTab; i++) {
    printf(" %d ;",tab[i]);
  }
  printf("\n\n");
}

/**
 * Affiche un plus court chemin entre deux sommets x et y. 
 *  
 * @param Chemin : Plus court chemin que l'on veut afficher.
 * @param IndiceX : sommet de départ du plus court chemin.
 * @param IndiceY : sommet d'arrivé du plus court chemin.
 * 
 * NB :
 *      - affiche les indices des sommets qui compose ce plus court chemin.
*/
void printChemin(int *Chemin, int IndiceX, int IndiceY) {

  printf("Chemin de %d à %d : \t",IndiceX, IndiceY);

  int i = 0;
  while( Chemin[i] != -1) {
    printf("%d ; ",Chemin[i]);
    i++;
  }
  printf("\n\n");

}