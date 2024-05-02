echo "SMI"
ls ./data/$1/smi_files | wc -l

echo "MOL"
ls ./data/$1/mol_files | wc -l

for i in ./data/$1/smi_files/*
do
    inter=$(basename "$i" '.smi')
    nom="$inter.mol"
    if ! test -f "./data/$1/mol_files/$nom"; then
        obgen "$i" > "./data/$1/mol_files/$nom" 2> /dev/null
    fi
done

