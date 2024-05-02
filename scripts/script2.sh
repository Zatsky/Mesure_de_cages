echo "creation des repertoires"
for i in "data/$1" "data/$1/smi_files/" "data/$1/mol_files/" "data/$1/dot_files_reduit/" "data/$1/dot_files_reduit/graphes_cycles/" "data/$1/dot_files_reduit/graphes_coins" "data/$1/png_files_reduit/" "data/$1/png_files_reduit/graphes_cycles" "data/$1/png_files_reduit/graphes_coins" "data/$1/results/" 
do 
     if [! -d $i ];then
	echo "Le dossier n'existe pas ($i). Il va être créé"
	mkdir $i
     fi
done

#bash scripts/smi2mol.sh  $1
bash scripts/script_supp_mol_vide.sh $1