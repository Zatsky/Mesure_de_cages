CFLAGS=-g -Wall
CC=gcc
CXX=g++

run_cage_mesure: analyse_cage
	#valgrind -v --leak-check=full --show-leak-kinds=all ./analyse_cage mesure mult
	./analyse_cage mesure mult

run_cage: analyse_cage
	#valgrind -v --leak-check=full --show-leak-kinds=all ./analyse_cage
	./analyse_cage
	#gdb ./analyse_cage

analyse_cage: analyse_cage.o utils_cage_moleculaire.o lecture_molecule_sdf.o graphe_cycles.o
	$(CC) $(CFLAGS) -o $@ $^

analyse_cage.o: analyse_cage.c analyse_cage.h
	$(CC) $(CFLAGS) -c $< -o $@

utils_cage_moleculaire.o: utils_cage_moleculaire.c utils_cage_moleculaire.h structure.h
	$(CC) $(CFLAGS) -c $< -o $@

graphe_cycles.o: graphe_cycles.c graphe_cycles.h structure.h
	$(CC) $(CFLAGS) -c $< -o $@
	
lecture_molecule_sdf.o: lecture_molecule_sdf.c lecture_molecule_sdf.h 
	$(CC) $(CFLAGS) -c $< -o $@

clean: 
	rm -f analyse_cage
	rm -f sortie
	rm -f *.o

clean_results:
	rm -rf results/
	rm -rf data/smi_files_reduit/
	rm -rf data/dot_files_reduit/
	rm -rf data/png_files_reduit/
