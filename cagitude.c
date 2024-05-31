#include "cagitude.h"

float calcul_mesure_coins_add(int length, char* line,float alpha){
    float total = 0.0;
    char c;
    int x = 0;
    int i = length; // Utilisez une variable pour itérer à travers les caractères de la ligne

    while (line[i] != '\0') {
        c = line[i];
        if (c == ' ') {
            x = 1;
        }
        if (c == ',') {
            x = 0;
        }
        if (x == 1) { // Vous pouvez directement vérifier ici si c est '1' ou '2'
            if (c == '1') {
                total += 1.0;
            } else if (c == '2') {
                total += alpha;
            }
        }
        i++; // Passer au caractère suivant
    }
	return total;
}

float calcul_mesure_coins_mult(int length, char* line,float alpha){
	float mesure[3];
    float total = 0.0;
    char c;
    int x = 0;
    int i = length; // Utilisez une variable pour itérer à travers les caractères de la ligne
	int y = 0;
    while (line[i] != '\0') {
        c = line[i];
        if (c == ',') {
            x = 0;
        }
        if (x == 1) { // Vous pouvez directement vérifier ici si c est '1' ou '2'
            if (c == '1') {
                mesure[y] = 1.0;
				y++;
            } else if (c == '2') {
                mesure[y] = alpha;
				y++;
            }
			else if(c == '3'){
				mesure[y] = 1.0;
				y++;
			}
			if(y>2){
			y = 0;
			total += (mesure[0] * mesure[1] * mesure[2]);
			}
        }
        if (c == ' ') {
            x = 1;
        }
        i++; // Passer au caractère suivant
    }
	return total;
}


float comparator(float* moy, FILE *f, float alpha) {
    char line[1000];
    char *comma_ptr;
    float temp = 0.0;
    float rate = 0.0;
    int cage = 0;
	float diff1 = 0.0;
	float diff2 = 0.0;
    size_t length;
	float nb =0.0;
    while (fgets(line, sizeof(line), f) != NULL) {
        comma_ptr = strchr(line, ',');
        if (comma_ptr != NULL) {
            length = comma_ptr - line + 1;
            if (line[length] == 'c') {
                cage = 1;
				nb ++;
            }
			else if(line[length] == 'p'){
                cage = 2;
				nb ++;
			}
            comma_ptr = strchr(comma_ptr + 1, ',');
			
            if (comma_ptr != NULL) {
                // Calculer la longueur de la partie à écrire
                size_t start_index = comma_ptr - line + 1;
                size_t end_index = strlen(line) - 1; // Exclure le '\n' à la fin de la ligne
                // Extrait le texte après la deuxième virgule
                char extracted_text[1000];
                strncpy(extracted_text, line + start_index, end_index - start_index);
                extracted_text[end_index - start_index] = '\0'; // Ajouter la terminaison de la chaîne

                // Convertir le texte extrait en float
                sscanf(extracted_text, "%f", &temp);
                // Calcul du taux en fonction du critère de comparaison
				diff1 = fabs(temp - moy[0]);
				diff2 = fabs(temp - moy[1]);
				//printf("is cage = %d  val = %f , diff=  %f et %f\n",cage , temp ,diff1, diff2);
                if (((cage == 1) && diff1 <= diff2) || ((cage == 2) && diff2 <= diff1)) {
                    rate += 100;
                }
                cage = 0;
            }
        
        }
    }
    return rate / nb;
}

void cage_ou_precage(float* moy, FILE *f, float alpha) {
    char line[1000];
    char *comma_ptr;
    float temp = 0.0;
	float diff1 = 0.0;
	float diff2 = 0.0;
	float diff3 = 0.0;
    while (fgets(line, sizeof(line), f) != NULL) {
        comma_ptr = strchr(line, ',');
        if (comma_ptr != NULL) {
            comma_ptr = strchr(comma_ptr + 1, ',');
            if (comma_ptr != NULL) {
                // Calculer la longueur de la partie à écrire
                size_t start_index = comma_ptr - line + 1;
                size_t end_index = strlen(line) - 1; // Exclure le '\n' à la fin de la ligne
                // Extrait le texte après la deuxième virgule
                char extracted_text[1000];
                strncpy(extracted_text, line + start_index, end_index - start_index);
                extracted_text[end_index - start_index] = '\0'; // Ajouter la terminaison de la chaîne

                // Convertir le texte extrait en float
                sscanf(extracted_text, "%f", &temp);
                // Calcul du taux en fonction du critère de comparaison
				diff1 = fabs(temp - moy[0]);
				diff2 = fabs(temp - moy[1]);
				diff3 = fabs(temp - moy[2]);
				//printf("is cage = %d  val = %f , diff=  %f et %f\n",cage , temp ,diff1, diff2);
                if ((diff1 <= diff2) && diff1 <= diff3) {
                    fprintf(f,"cage,");
                }
				else if ((diff2 <= diff1) && diff2 <= diff3) {
                    fprintf(f,"pre cage,");
                }
				else {
					fprintf(f,"non cage,");
				}
            }
        }
    }
}


int main(int argc, char *argv[])
{	
	char * RESULTS_LISTE_COINS = "data/DEFAULT/results/liste_coins_reduit.csv";
	char * RESULTS_LISTE_MESURE = "data/DEFAULT/results/liste_mesure_alpha.csv";

	if (argc>2 && strcmp(argv[2], "CHEBI") == 0){ 
		RESULTS_LISTE_COINS = "data/CHEBI/results/liste_coins_reduit.csv";
		RESULTS_LISTE_MESURE = "data/CHEBI/results/liste_mesure_alpha.csv";
	}
	else if (argc>2 && strcmp(argv[2], "LOTUS") == 0){
		RESULTS_LISTE_COINS = "data/LOTUS/results/liste_coins_reduit.csv";
		RESULTS_LISTE_MESURE = "data/LOTUS/results/liste_mesure_alpha.csv";
	}

    if (argc > 1) {
		float alpha = 3.75; //à partir de 3.8 stabilisation à 87.5 pour mult
		float temp = 0.0;
        int nombreDeLignes = 0;
        char caractere;
		FILE *f_mesure,*f_liste;
		char *comma_ptr;

        f_liste = fopen(RESULTS_LISTE_COINS, "r"); // liste des coins par molécule (poids des sommets et poids des arêtes)
		if(f_liste == NULL)
		{
			printf("Impossible d'ouvrir le fichier %s\n", RESULTS_LISTE_COINS);
			exit(2);
		}

        while ((caractere = fgetc(f_liste)) != EOF) {
        // Si le caractère actuel est un retour à la ligne, incrémenter le compteur de lignes
            if (caractere == '\n') {
                nombreDeLignes++;
                }
        }
        char line[nombreDeLignes];
        fclose(f_liste);
        f_liste = fopen(RESULTS_LISTE_COINS, "r");
		
		f_mesure = fopen(RESULTS_LISTE_MESURE, "w"); // liste des coins par molécule (poids des sommets et poids des arêtes)
		if(f_mesure == NULL)
		{
			printf("Impossible d'ouvrir le fichier %s\n", RESULTS_LISTE_MESURE);
			exit(3);
		}
		int xz = 0;
		while (fgets(line, sizeof(line), f_liste) != NULL) {
			comma_ptr = strchr(line, ',');
			if (comma_ptr != NULL) {
				if (argc>2 && strcmp(argv[2], "DEFAULT") == 0){ 
					comma_ptr = strchr(comma_ptr+1, ',');
				}
				// Calculer la longueur de la partie à écrire
				size_t length = comma_ptr - line + 1;
				// Écrire la partie dans le fichier de sortie
				comma_ptr = strchr(comma_ptr + 1, ',');
				printf(comma_ptr);
				// Appeler calcul_mesure_coins pour traiter cette ligne
				if (argc > 1 && strcmp(argv[1], "add") == 0){
					temp =calcul_mesure_coins_add(length,line, alpha);
				}
				else if (argc > 1 && strcmp(argv[1], "mult") == 0){
					temp =calcul_mesure_coins_mult(length,line, alpha);
				}
				// Ajouter un retour à la ligne dans le fichier de sortie
                if (temp >0.0){
					fwrite(line, sizeof(char), length, f_mesure);
                    fprintf(f_mesure, "%f,\n", temp);
                }
				else{
					xz++;
				}
			}
		}
		fclose(f_mesure);
		//f_mesure = fopen(RESULTS_LISTE_MESURE, "r"); // liste des coins par molécule (poids des sommets et poids des arêtes)
		
		//rate = comparator(moy,f_mesure,alpha);
		fclose(f_liste);
		printf("%d à 0\n",xz);
		//fclose(f_mesure);
    
    /*
	if (argc > 2 && strcmp(argv[1], "mesurement") == 0) {
		float alpha = 3.75; //à partir de 3.8 stabilisation à 87.5 pour mult
		float cage = 0.0;
		float pcage = 0.0;
		float ncage = 0.0;
		float temp = 0.0;
		float rate = 0.0;
		float cmoy = 0.0;
		float pmoy = 0.0;
		float nmoy = 0.0;
        int nombreDeLignes = 0;
        char caractere;
		FILE *f_mesure,*f_liste;
		char *comma_ptr;

        f_liste = fopen(RESULTS_LISTE_COINS, "r"); // liste des coins par molécule (poids des sommets et poids des arêtes)
		if(f_liste == NULL)
		{
			printf("Impossible d'ouvrir le fichier %s\n", RESULTS_LISTE_COINS);
			exit(2);
		}

        while ((caractere = fgetc(f_liste)) != EOF) {
        // Si le caractère actuel est un retour à la ligne, incrémenter le compteur de lignes
        if (caractere == '\n') {
            nombreDeLignes++;
            }
        }
        char line[nombreDeLignes];
        fclose(f_liste);
        fopen(RESULTS_LISTE_COINS, "r");
		
		f_mesure = fopen(RESULTS_LISTE_MESURE, "w"); // liste des coins par molécule (poids des sommets et poids des arêtes)
		if(f_mesure == NULL)
		{
			printf("Impossible d'ouvrir le fichier %s\n", RESULTS_LISTE_MESURE);
			exit(3);
		}
		if (argc > 1 && strcmp(argv[1], "add") == 0){
			while (fgets(line, sizeof(line), f_liste) != NULL) {
				comma_ptr = strchr(line, ',');
				if (comma_ptr != NULL) {
					size_t length1 = comma_ptr - line + 1;
					comma_ptr = strchr(comma_ptr + 1, ',');
					if (comma_ptr != NULL) {
						// Calculer la longueur de la partie à écrire
						size_t length = comma_ptr - line + 1;
						// Écrire la partie dans le fichier de sortie
						comma_ptr = strchr(comma_ptr + 1, 
    # Compter le nombre de molécules pour chaque mesure',');
						fwrite(line, sizeof(char), length, f_mesure);
						// Appeler calcul_mesure_coins pour traiter cette ligne
						temp =calcul_mesure_coins_add(length,line, alpha);
						// Ajouter un retour à la ligne dans le fichier de sortie
						fprintf(f_mesure, "%f,\n", temp);
						if (line[length1] == 'c') {
							cage = cage +1.0;
							cmoy += temp;
						}
						else if(line[length1] == 'p'){
							pcage = pcage + 1.0; 
							pmoy += temp;
						}
					}
				}
			}
		}
		else if (argc > 2 && strcmp(argv[1], "mult") == 0){
			while (fgets(line, sizeof(line), f_liste) != NULL) {
				comma_ptr = strchr(line, ',');
				if (comma_ptr != NULL) {
					size_t length1 = comma_ptr - line + 1;
					comma_ptr = strchr(comma_ptr + 1, ',');
					if (comma_ptr != NULL) {
						// Calculer la longueur de la partie à écrire
						size_t length = comma_ptr - line + 1;
						// Écrire la partie dans le fichier de sortie
						comma_ptr = strchr(comma_ptr + 1, ',');
						fwrite(line, sizeof(char), length, f_mesure);
						// Appeler calcul_mesure_coins pour traiter cette ligne
						temp =calcul_mesure_coins_mult(length,line, alpha);
						// Ajouter un retour à la ligne dans le fichier de sortie
						fprintf(f_mesure, "%f,\n", temp);
						if (line[length1] == 'c') {
							cage = cage +1.0;
							cmoy += temp;
						}
						else if(line[length1] == 'p'){
							pcage = pcage + 1.0; 
							pmoy += temp;
						}
						else {
							ncage = ncage + 1.0; 
							nmoy += temp;
						}
					}
				}
			}
		}
		float moy[3];
		moy[0] = cmoy /cage;
		moy[1] = pmoy / pcage;
		moy[2] = nmoy / ncage;
		fclose(f_mesure);
		//f_mesure = fopen(RESULTS_LISTE_MESURE, "r"); // liste des coins par molécule (poids des sommets et poids des arêtes)
		
		//rate = comparator(moy,f_mesure,alpha);
		printf("%f ",rate);
		fclose(f_liste);
		//fclose(f_mesure);
		FILE* f_classification = fopen(CHEBI,"w");
		cage_ou_precage(moy,f_classification,alpha);
		fclose(f_classification);
        */
        

	}
    exit(0);
}