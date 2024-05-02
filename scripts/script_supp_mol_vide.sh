# Parcourez tous les fichiers .mol dans le répertoire
for i in ./data/$1/mol_files/*
do
    # Vérifiez si le fichier est vide en utilisant la commande test -s
    if [ ! -s "$i" ]; then
        # Supprimez le fichier s'il est vide
        rm "$i"
        echo "Fichier vide supprimé : $i"
    fi
done