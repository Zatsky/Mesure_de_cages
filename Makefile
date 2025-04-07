CFLAGS=-g -Wall -fPIC
CC=gcc
CXX=g++
OBJS = similarite.o numerotation.o affichage.o structure.o Chemins.o graphe_cycles.o utils_cage_moleculaire.o 
OBJS2 = analyse_cage.o utils_cage_moleculaire.o lecture_molecule_sdf.o graphe_cycles.o structure.o

run: analyse_cage cagitude
	./scripts/script2.sh $(ARGS)
	#valgrind -v --leak-check=full --show-leak-kinds=all ./analyse_cage $(ARGS)
	./analyse_cage $(ARGS)
	#./cagitude $(ARGS) add

python_lib: similarite.so

comparaison: similarite
	#valgrind -v --leak-check=full --show-leak-kinds=all --track-origins=yes ./similarite
	./similarite

mesure: cagitude
	#valgrind -v --leak-check=full --show-leak-kinds=all --track-origins=yes ./cagitude DEFAULT add
	./cagitude LOTUS add connexe

cliques: find_4_cliques
	./find_4_cliques

run_cage: analyse_cage
	#valgrind -v --leak-check=full --show-leak-kinds=all ./analyse_cage $(ARGS)
	./analyse_cage $(ARGS)
	#gdb ./analyse_cage 

analyse_cage: $(OBJS2)
	$(CC) -o $@ $^ -lm

similarite.so: $(OBJS)
	$(CC) -shared -o $@ $^

similarite: $(OBJS)
	$(CC) -o $@ $^ -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@


clean: 
	rm -f analyse_cage
	rm -f sortie
	rm -f *.o

clean_results:
	rm -rf data/CHEBI/results/*
	rm -rf data/LOTUS/results/*
	rm -rf data/CHEBI/dot_files_reduit/*
	rm -rf data/LOTUS/dot_files_reduit/*
	rm -rf data/CHEBI/png_files_reduit/*
	rm -rf data/LOTUS/png_files_reduit/*
	rm -rf data/CHIMISTE_1/results/*
	rm -rf data/CHIMISTE_2/results/*
	rm -rf data/CHIMISTE_1/dot_files_reduit/*
	rm -rf data/CHIMISTE_2/dot_files_reduit/*
	rm -rf data/CHIMISTE_1/png_files_reduit/*
	rm -rf data/CHIMISTE_2/png_files_reduit/*


clean_pdf:
	rm -rf data/CHEBI/ID/*
	rm -rf data/LOTUS/ID/*
	rm -rf data/CHIMISTE_1/ID/*
	rm -rf data/CHIMISTE_2/ID/*