echo "creation des repertoires"
for i in "data/$1" "data/$1/smi_files/" "data/$1/mol_files/" "data/$1/dot_files_reduit/" "data/$1/dot_files_reduit/graphes_cycles/" "data/$1/dot_files_reduit/graphes_coins" "data/$1/png_files_reduit/" "data/$1/png_files_reduit/graphes_cycles" "data/$1/png_files_reduit/graphes_coins" "data/$1/results/" "data/$1/ID/" "data/$1/png_files_reduit/pymol"
do 
     if [ ! -d $i ];then
	echo "Le dossier n'existe pas ($i). Il va être créé"
	mkdir $i
     fi
done
directory = "data/$1/mol_files/"
num_files=$(ls -p $directory | grep -v / | wc -l)
echo $num_files
#if [ "$num_files" -ne 276518 ]; then
     #python scripts/script_mol3.py  $1
#fi
bash scripts/script_supp_mol_vide.sh $1