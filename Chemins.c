#include "Chemins.h"
#include "affichage.h"

void printChemins() {
  printf("\nfichier Chemins.c\n");
}

/**
 * @file Chemins.c
 * @brief fichier contenant les fonctions liés aux plus courts chemins.
 */



/**
 * Pour chaque sommets v du graphe G, construit un DAG où v sera la racine d'une anti-arborescence.
 * Et où chaque parcours de x à v dans le DAG de v, correspond à un plus court chemin de x à v dans G.
 * 
 * @param G : graphe G duquel on extrait les DAGs.
 * 
 * @return Tableau de DAGs.
 * 
 * Step 1 : Effectue un parcours en largeur du graphe G à partir du sommet v.
 * 
 * Step 2 : Construit l'arc (y, z) dans le DAG de v, si l'arête (z,y) appartient à un plus court chemin de v à y dans G.
 * 
 * NB : 
 *      - le calcul des plus court chemins considère le nombre d'arête et non leurs poids (comme si chaque arête avait un poids égal à 1).
*/
graphe ** GenererDAG(graphe *G) {

  int n = G->nb_sommets;

  graphe **DAG;
  
  DAG = malloc(sizeof(graphe *) * n);
  if (DAG == NULL) {
    printf("PB malloc DAG\n");
    exit(0);
  }

  for(int i = 0; i < n; i++) {
    DAG[i] = initGraphe(n);
  }

  for(int i = 0; i < n; i++) {


    /* --- Step 1 --- */

    int ListeAttente[n];
    int tailleChemins[n];
    int somMarque[n];

    for(int j = 0; j < n; j++) {
      ListeAttente[j] = -1;
      somMarque[j] = 0;
      tailleChemins[j] = __INT_MAX__;
    }

    int indiceListe = 0;
    int indiceAjoute = 1;
    ListeAttente[indiceListe] = i;
    somMarque[i] = 1;
    tailleChemins[i] = 0;

    while (ListeAttente[indiceListe] != -1) {

      int somTraite = ListeAttente[indiceListe];

      for(int j = 0; j < n; j++) {


        /* --- Step 2 --- */

        if(G->adjacence[somTraite][j] == 1) {

          if(tailleChemins[somTraite] + 1 <= tailleChemins[j]) {
            tailleChemins[j] = tailleChemins[somTraite] + 1;
            DAG[i]->adjacence[j][somTraite] = 1;

            if(somMarque[j] == 0) {
              somMarque[j] = 1;
              ListeAttente[indiceAjoute] = j;
              indiceAjoute++;
              indiceAjoute = indiceAjoute % n;
            }

          }

        }
      }
      ListeAttente[indiceListe] = -1;
      indiceListe++;
      indiceListe = indiceListe % n;
    }
  }
  //printf("\nGénération des DAGs terminée\n");
  return DAG;
}


/**
 * Calcul le plus court chemins entre deux sommets d'un graphe.
 * 
 * @param DAG : tableau des DAGs.
 * @param IndiceX : indice du sommet de départ.
 * @param IndiceY : indice du sommet d'arrivé.
 * @param tailleGraphe : Nombre de sommets dans le graphe.
 * 
 * @return tableau d'entier de taille tailleGraphe. Contients les indices des sommets du plus court chemins dans l'ordre.
 * 
 * NB : 
 *      - le plus court chemin renvoyé est celui de plus petite numérotation (au sens de la numérotation des sommets) possible.
 *      - si l'algorithme à un moment un choix à faire entre deux sommets, pour définir un plus court chemin, il prendra le plus petit sommets des deux.
*/
int * plusCourtCheminDexAy(graphe **DAG, int IndiceX, int IndiceY, int tailleGraphe) {

  //printGraphe(DAG[IndiceY]);
  int *CheminXY;

  CheminXY = malloc(sizeof(int) * tailleGraphe);
  if(CheminXY == NULL) {
    printf("PB malloc Chemin\n");
    exit(0);
  }

  for(int i = 0; i < tailleGraphe; i++) {
    CheminXY[i] = -1;
  }

  int indiceParcours = IndiceX;
  int indiceChemin = 0;
  CheminXY[indiceChemin] = IndiceX;
  indiceChemin++;

  

  while(indiceParcours != IndiceY) {
    for(int i = 0; i < tailleGraphe; i++) {
      if (DAG[IndiceY]->adjacence[indiceParcours][i] == 1) {
        CheminXY[indiceChemin] = i;
        indiceChemin++;
        indiceParcours = i;
        i = tailleGraphe;
      }
    }
  }

  //printChemin(CheminXY, IndiceX, IndiceY);

  return CheminXY;
}

/**
 * Calcul la distance (taille d'un pcc) entre 2 sommets.
 * 
 * @param G : Graphe contenant les 2 sommets.
 * @param som1 : Premier sommet.
 * @param som2 : Deuxième sommet.
 * @param DAG : tableau des DAGs du graphe G.
 * 
 * @return renvoie la distance (en entier) entre les 2 sommets.
 */
int distanceSoms(graphe *G, int som1, int som2, graphe **DAG) {

  int *pcc;
  pcc = plusCourtCheminDexAy(DAG, som1, som2, G->nb_sommets);

  //printChemin(pcc, som1, som2);

  int taille = 0;

  for(int i = 0; i < G->nb_sommets; i++) {
    if(pcc[i] != -1) {
      taille++;
    }
  }

  free(pcc);

  return taille;
}

