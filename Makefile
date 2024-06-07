CFLAGS=-g -Wall
CC=gcc
CXX=g++

run: analyse_cage cagitude
	./scripts/script2.sh $(ARGS)
	./analyse_cage $(ARGS)
	./cagitude $(ARGS) add

mesure: cagitude
	#valgrind -v --leak-check=full --show-leak-kinds=all --track-origins=yes ./cagitude DEFAULT add
	./cagitude CHEBI add connexe


cagitude: cagitude.o
	$(CC) $(CFLAGS) -o $@ $^

cagitude.o: cagitude.c cagitude.h
	$(CC) $(CFLAGS) -c $< -o $@

run_cage: analyse_cage
	#valgrind -v --leak-check=full --show-leak-kinds=all ./analyse_cage $(ARGS)
	./analyse_cage $(ARGS)
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
	rm -rf data/DEFAULT/results/
	rm -rf data/CHEBI/results/
	rm -rf data/LOTUS/results/
	rm -rf data/DEFAULT/dot_files_reduit/
	rm -rf data/CHEBI/dot_files_reduit/
	rm -rf data/LOTUS/dot_files_reduit/
	rm -rf data/DEFAULT/png_files_reduit/
	rm -rf data/CHEBI/png_files_reduit/
	rm -rf data/LOTUS/png_files_reduit/


clean_pdf:
	rm -rf data/CHEBI/ID/*