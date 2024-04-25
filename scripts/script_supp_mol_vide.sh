repertoire="data/chebi_smi/"

# Parcourez tous les fichiers .mol dans le répertoire
find "$repertoire" -type f -name "*.mol" | while read fichier; do
    # Vérifiez si le fichier est vide en utilisant la commande test -s
    if [ ! -s "$fichier" ]; then
        # Supprimez le fichier s'il est vide
        rm "$fichier"
        echo "Fichier vide supprimé : $fichier"
    fi
done