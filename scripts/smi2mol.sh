echo "SMI"
ls ./data/CHEBI/chebi_smi | wc -l

echo "MOL"
ls ./data/CHEBI/chebi_mol | wc -l

for i in ./data/CHEBI/chebi_smi/*
do
    inter=$(basename "$i" '.smi')
    nom="$inter.mol"
    if ! test -f "./data/CHEBI/chebi_mol/$nom"; then
        obgen "$i" > "./data/CHEBI/chebi_mol/$nom" 2> /dev/null
        echo $nom
    fi
done

