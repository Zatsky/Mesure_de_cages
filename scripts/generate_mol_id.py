from pylatex import Document, Section, Command, Figure, SubFigure, Package, Itemize, Subsection, StandAloneGraphic
from pylatex.base_classes.command import Options
from pylatex.utils import fix_filename,NoEscape
import os
import subprocess
from pathlib import Path
import sys
from pymol import cmd
import csv
from rdkit import Chem
from PIL import Image
from PyPDF2 import PdfMerger
import pikepdf
import generation_graphiques
from pathlib import Path

def count_mol_files_with_more_than_250_atoms(arg1):
    count = 0
    dir_path = f"data/"+arg1+f"/mol_files"
    # Parcourir tous les fichiers dans le répertoire
    for file_name in os.listdir(dir_path):
        if file_name.endswith('.mol'):
            file_path = os.path.join(dir_path, file_name)
            
            # Lire le fichier .mol
            mol = Chem.MolFromMolFile(file_path)
            
            if mol is not None:
                num_atoms = mol.GetNumAtoms()
                if num_atoms > 250:
                    count += 1
            else:
                print(f"Erreur de lecture du fichier : {file_name}")

    return count

def count_files_in_directory(arg1):
    dir_path = f"data/"+arg1+f"/mol_files"
    try:
        # Créer un objet Path pour le répertoire
        path = Path(dir_path)
        # Compter les fichiers dans le répertoire
        return sum(1 for _ in path.iterdir() if _.is_file())
    except Exception as e:
        print(f"Erreur : {e}")
        return 0

def compress_pdf(input_pdf_path, output_pdf_path, quality='ebook'):
    """
    Compress a PDF file using Ghostscript.
    
    Args:
    - input_pdf_path: str, path to the input PDF file.
    - output_pdf_path: str, path to save the compressed PDF file.
    - quality: str, quality setting for Ghostscript (default is 'ebook').
               Options: 'screen', 'ebook', 'printer', 'prepress', 'default'.
    """
    try:
        # Define the Ghostscript command
        gs_command = [
            'gs',
            '-sDEVICE=pdfwrite',
            f'-dPDFSETTINGS=/{quality}',
            '-dNOPAUSE',
            '-dQUIET',
            '-dBATCH',
            f'-sOutputFile={output_pdf_path}',
            input_pdf_path
        ]
        
        # Run the Ghostscript command
        result = subprocess.run(gs_command, check=True, capture_output=True, text=True)
        
        # Check for errors in the output
        if result.returncode == 0:
            print(f"Compressed PDF saved as: {output_pdf_path}")
        else:
            print(f"Error compressing PDF: {result.stderr}")

    except subprocess.CalledProcessError as e:
        print(f"Error compressing PDF: {e}")
    except FileNotFoundError:
        print("Ghostscript not found. Please ensure it is installed and accessible in your PATH.")

def optimize_pdf(input_pdf_path, output_pdf_path):
    """
    Optimize PDF using pikepdf by linearizing and optimizing the structure.
    
    Args:
    - input_pdf_path: str, path to the input PDF file.
    - output_pdf_path: str, path to save the optimized PDF file.
    """
    with pikepdf.open(input_pdf_path,allow_overwriting_input=True) as pdf:
        # Linearize the PDF for faster web view
        pdf.save(output_pdf_path, linearize=True)
        print(f"Optimized PDF saved as: {output_pdf_path}")

class donnee_csv:
    def __init__(self,nom,cagitude):
        self.cagitude = cagitude
        self.nom = nom

    def other_cagitude(self,path):
        with open(path, 'r') as fichier:
            lecteur_csv = csv.reader(fichier)
            for i,ligne in enumerate(lecteur_csv):
                if ligne[0] == self.nom:
                    return int(float(ligne[1]))

def get_arete_coins(arg1,arg2):
    graphe_coin = 'data/' + arg1 + '/graphes_coins/'+ arg2 +'.csv'
    with open(graphe_coin, 'r') as file:
        # Lire la première ligne du fichier
        first_line = file.readline().strip()
        # Séparer la ligne en parties en utilisant '_' comme séparateur
        parts = first_line.split('_')
        if len(parts) > 1:
            try:
                # Retourner le deuxième élément converti en entier
                return int(parts[1])
            except ValueError:
                # Si la conversion échoue, afficher un message d'erreur
                print(f"Erreur : '{parts[1]}' n'est pas un entier valide.")
                return None
        else:
            print("Erreur : La première ligne ne contient pas de deuxième nombre.")
            return None

def lire_donnees_csv(nom_fichier,arg1):
    
    mesures = []
    nb_mesures = []
    cagitude = 0
    j = 0
    x = 0
    with open(nom_fichier, 'r') as fichier:
        lecteur_csv = csv.reader(fichier)
        for i,ligne in enumerate(lecteur_csv):
            try:
                cage = int(float(ligne[1]))
                #print(str(cage))
                if get_arete_coins(arg1,ligne[0]) <1000:
                    mesure = donnee_csv(ligne[0],cage)
                    mesures.append(mesure)
                else:
                    x+=1
            except ValueError:
                pass
    mesures.sort(key=lambda z: z.cagitude)
    cagitude = mesures[0].cagitude
    for i, mesure in enumerate(mesures):
        if cagitude != mesure.cagitude:
            nb_mesures.append(j)
            cagitude = mesure.cagitude
            j = 0
        j +=1
    nb_mesures.append(j)
    print("Il y'a "+ str(len(nb_mesures))+ " mesures différentes et "+ str(x) + "supérieures à 250")
    return mesures,nb_mesures,x

#Renvoie le smiles, le nombres d'atomes, ainsi que le nombres de liaisons d'une molécule
def get_molecule_info(arg1,arg2):
    mol_file = 'data/' + arg1 + '/mol_files/'+ arg2 +'.mol'
    try:
        # Lire la molécule à partir du fichier .mol
        mol = Chem.MolFromMolFile(mol_file)
        if mol is None:
            print(f"Erreur de lecture du fichier {mol_file}")
            return None, None, None
        
        # Récupérer le SMILES
        smiles = Chem.MolToSmiles(mol)
        
        # Récupérer le nombre d'atomes
        num_atoms = mol.GetNumAtoms()
        
        # Récupérer le nombre de liaisons
        num_bonds = mol.GetNumBonds()
        
        return smiles, num_atoms, num_bonds
    except Exception as e:
        print(f"Erreur: {e}")
        return None, None, None
    
#Donne la valeur de cagitude pour une molécule donnée
def get_cagitude(arg1,arg2):
    csv_file = 'data/' + arg1 + '/results/liste_mesure_alpha.csv'
    try:
        with open(csv_file, mode='r') as file:
            reader = csv.reader(file)
            for row in reader:
                if row[0] == arg2:
                    return row[1]
        return 0
    except FileNotFoundError:
        print(f"Le fichier {csv_file} n'existe pas.")
        return None
    
#Genere une image png représentant la molécule arg2
def generate_image_mol(arg1, arg2):
    cmd.load("data/" + arg1 +"/mol_files/"+ arg2 +".mol")
    path = "data/" + arg1 +"/png_files_reduit/pymol/"+ arg2 +".png"
    # Sauvegarde de l'image au format PNG
    cmd.png(path, width=200, height=200, dpi=300, ray=1)
    cmd.delete('all')

    convert_and_resize_image(path,path.replace('.png', '.jpg'))
    return path.replace('.png', '.jpg')

def escape_latex_special_chars(smiles):
    """Échappe les caractères spéciaux pour LaTeX dans une chaîne SMILES."""
    return smiles.replace('_', r'\_').replace('&', r'\&').replace('%', r'\%').replace('$', r'\$').replace('#', r'\#').replace('{', r'\{').replace('}', r'\}')

def format_smiles(smiles, line_length=70):
    """Formate une chaîne SMILES en insérant des sauts de ligne pour éviter qu'elle ne dépasse les marges."""
    escaped_smiles = escape_latex_special_chars(smiles)
    return '\n'.join([escaped_smiles[i:i+line_length] for i in range(0, len(escaped_smiles), line_length)])

def convert_and_resize_image(image_path, output_path, quality=100, size=(200, 200)):
    # Ouvrir l'image
    with Image.open(image_path) as img:
        # Convertir en JPEG avec la qualité spécifiée
        img = img.convert("RGB")
        img.save(output_path, quality=quality)
        print(f"Converted and saved image: {output_path}")

        # Redimensionner l'image
        img.thumbnail(size)
        img.save(output_path)

#Genere les images des graphes  de cycles et des graphes de coins de la molécule arg2
def generate_images(arg1, arg2):
    dot_cycle = f'data/{arg1}/dot_files_reduit/graphes_cycles/{arg2}.dot'
    dot_coin = f'data/{arg1}/dot_files_reduit/graphes_coins/{arg2}.dot'

    # Générer les images à partir des fichiers .dot
    image_path_coin = f'data/{arg1}/png_files_reduit/graphes_coins/{arg2}.png'
    command = ['dot', '-Tpng', dot_coin, '-o', image_path_coin]
    subprocess.run(command, check=True)
    print(f"Generated image: {image_path_coin}")
    
    image_path_cycle = f'data/{arg1}/png_files_reduit/graphes_cycles/{arg2}.png'
    command = ['dot', '-Tpng', dot_cycle, '-o', image_path_cycle]
    subprocess.run(command, check=True)
    print(f"Generated image: {image_path_cycle}")

    convert_and_resize_image(image_path_coin, image_path_coin.replace('.png', '.jpg'))
    convert_and_resize_image(image_path_cycle, image_path_cycle.replace('.png', '.jpg'))

    return image_path_coin.replace('.png', '.jpg'), image_path_cycle.replace('.png', '.jpg')

#Ecris le document latex contenant toute les données importantes de la molecule arg2
def generate_latex_document(arg1, arg2):
    output_file = f'data/{arg1}/ID/{arg2}'
    doc = Document(documentclass='article', document_options='a4paper')
    doc.packages.append(Package('titlesec'))
    doc.packages.append(Package('geometry'))
    doc.packages.append(Package('graphicx'))
    doc.packages.append(Package('float'))
    doc.packages.append(Package('hyperref'))    
    doc.packages.append(Command('usepackage', 'fancyhdr'))
    doc.preamble.append(Command('date', ''))
    doc.preamble.append(NoEscape(r'\pagestyle{fancy}'))
    doc.preamble.append(NoEscape(r'\fancyhf{}'))
    doc.preamble.append(NoEscape(r'\rfoot{\thepage}'))
    image_coin, image_cycle = generate_images(arg1, arg2)
    image_mol = generate_image_mol(arg1,arg2)
    images = [image_cycle, image_coin,image_mol]
    

    for image in images:
        if not Path(image).exists():
            print(f"Error: Image file {image} does not exist.")
            return

    with doc.create(Section(f"{arg2}'s ID")):
        with doc.create(Itemize()) as itemize:
            smiles,num_atoms,num_bonds = get_molecule_info(arg1,arg2)
            if smiles != None:
                smiles = format_smiles(smiles)
                itemize.add_item(NoEscape(r"SMILES: \begin{verbatim}" + smiles + r"\end{verbatim}"))
            else:
                 itemize.add_item(f"Nombres d'atomes: {smiles}")
            itemize.add_item(f"Nombres d'atomes: {num_atoms}")
            itemize.add_item(f"Nombres de liaisons: {num_bonds}")
            value = get_cagitude(arg1, arg2)
            itemize.add_item(f"Valeur de cagitude = {value}")

        with doc.create(Figure(position='H')) as fig:
            image_path_new = os.path.join('..', 'png_files_reduit', 'pymol', f'{arg2}.jpg')
            fig.add_image(image_path_new, width=NoEscape(r'0.4\textwidth'))
            fig.add_caption("Molécule en 3D")

        with doc.create(Figure(position='H')) as graph:
            with doc.create(SubFigure(position='b', width=NoEscape(r'0.4\textwidth'))) as left_image:
                image_path_cycles = os.path.join('..', 'png_files_reduit', 'graphes_cycles', f'{arg2}.jpg')
                left_image.add_image(image_path_cycles, width=NoEscape(r'\textwidth'))
                left_image.add_caption("Graphe des Cycles")

            with doc.create(SubFigure(position='b', width=NoEscape(r'0.4\textwidth'))) as right_image:
                image_path_coins = os.path.join('..', 'png_files_reduit', 'graphes_coins', f'{arg2}.jpg')
                right_image.add_image(image_path_coins, width=NoEscape(r'\textwidth'))
                right_image.add_caption("Graphe des Coins")
    if arg1 == "CHEBI":
        mol = arg2.replace("CHEBI_", "")
        link = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:' + mol
        doc.append(NoEscape(r"Lien: \href{" + link + r"}{" + link + r"}"))
    doc.generate_tex(output_file)
    doc.generate_pdf(output_file, clean_tex=True)
    return output_file + '.pdf'


def create_bestiaire(arg1):

    opti_path = f'data/{arg1}/ID/Bestiaire_des_coins_{arg1}.pdf'
    pdf_path = f'data/{arg1}/ID/Bestiaire_des_coins'
    compressed_path = f'data/{arg1}/ID/Bestiaire_des_coins_compressed.pdf'
    csv_file =f'data/{arg1}/results/liste_mesure_alpha_connexe.csv'
    csv_file2 =f'data/{arg1}/results/liste_mesure_alpha.csv'
    doc = Document(documentclass='article', document_options='a4paper')
    
    # Ajouter les packages nécessaires
    doc.packages.append(Command('usepackage', 'hyperref'))
    doc.packages.append(Command('usepackage', 'geometry'))
    doc.packages.append(Command('geometry', 'a4paper'))
    doc.packages.append(Command('usepackage', 'fancyhdr'))
    doc.packages.append(Package('xcolor', options='dvipsnames'))
    doc.packages.append(Package('titlesec'))
    doc.packages.append(Package('geometry'))
    doc.packages.append(Package('graphicx'))
    doc.packages.append(Package('float'))   
    doc.preamble.append(Command('date', ''))
    doc.preamble.append(NoEscape(r'\pagestyle{fancy}'))
    doc.preamble.append(NoEscape(r'\fancyhf{}'))
    doc.preamble.append(NoEscape(r'\rfoot{\thepage}'))
    doc.preamble.append(NoEscape(r'\setcounter{tocdepth}{1}'))
    doc.packages.append(Command('usepackage', 'fontenc', options='T1'))
    doc.packages.append(Command('usepackage', 'babel', options='french'))
    doc.preamble.append(NoEscape(r'\usepackage[french]{babel}'))
    doc.preamble.append(Command('hypersetup', NoEscape('colorlinks=true, linkcolor=blue, urlcolor=blue')))

    if arg1 == "CHEBI":
        seg_fault = 6
    elif arg1 == "LOTUS":
        seg_fault = 27
    else:
        seg_fault = 0
    all_files = count_files_in_directory(arg1)

    mesures,nb_mesure,trop_gros = lire_donnees_csv(csv_file,arg1)
    too_big = count_mol_files_with_more_than_250_atoms(arg1)

    no_corner = all_files - seg_fault - too_big -len(mesures) - trop_gros

    # Ajouter le titre principal
    doc.preamble.append(Command('title', 'Bestiaire des coins de '+arg1))
    doc.preamble.append(Command('date', NoEscape(r'\today')))
    doc.preamble.append(Command('author', 'Docherty Ronan'))
    doc.append(NoEscape(r'\maketitle'))
    if not Path(f'data/{arg1}/results/histogramme_discret.png').exists():
        generation_graphiques.gen_graphique(arg1)
    graph_path = f'../results/histogramme_discret.png'
    with doc.create(Section('Introduction', numbering=False)):
        doc.append("Ce document regroupe l'ensemble des molécules possédant au moins un coin et ayant moins de 250 atomes.")
        doc.append(NoEscape(r'\begin{figure}[h!]'))
        doc.append(NoEscape(r'\centering'))
        doc.append(NoEscape(r'\includegraphics[width=0.8\textwidth]{' + graph_path + r'}'))
        doc.append(NoEscape(r'\caption{Histogramme discret de la cagitude ChEBI}'))
        doc.append(NoEscape(r'\end{figure}'))
    doc.append(NoEscape(r'\newpage'))
    
    pluriel = ["cette","",""]

    conditions = [seg_fault,too_big,trop_gros,no_corner]
    if conditions.count(0) <4:
        if conditions.count(0) <3:
            pluriel =  ["ces","s"," dans cet ordre"]
        with doc.create(Section(('Contrainte'+pluriel[1]), numbering=False)):
                
            doc.append("Sur les " +str(all_files)+ " molécules intiales de " + arg1 + ", nous avons dû appliquer "+pluriel[0]+" restriction" + pluriel[1]+" excluant"+pluriel[2]+":")
            with doc.create(Itemize()) as itemize:
                if seg_fault>0:
                    if seg_fault == 1:
                        pluriel[1] = ""
                    else:
                        pluriel[1] = "s"
                    itemize.add_item(NoEscape(str(seg_fault)+r" molécule"+pluriel[1]+r" causant des segmentation fault."))
                if too_big>0:
                    if too_big == 1:
                        pluriel[1] = ""
                    else:
                        pluriel[1] = "s"
                    itemize.add_item(NoEscape(str(too_big) + r" molécule"+pluriel[1]+r" ayant plus de 250 atomes, car leurs graphes de cycles est trop long à génerer."))
                if no_corner>0:
                    if no_corner == 1:
                        pluriel[1] = ""
                    else:
                        pluriel[1] = "s"
                    itemize.add_item(NoEscape(str(no_corner) + r" molécule"+pluriel[1]+ r" n'ayant aucuns coins."))
                if trop_gros>0:
                    if trop_gros == 1:
                        pluriel[1] = ""
                    else:
                        pluriel[1] = "s"
                    itemize.add_item(NoEscape(str(trop_gros) + r" molécule"+pluriel[1]+r" avec plus de 1000 liaisons sur leurs graphes de coins, car l'image est trop grosse."))
            doc.append(r'Et nous finissons donc avec '+str(len(mesures))+r' molécules.')
            
    with doc.create(Section('Représentation', numbering=False)):
        doc.append("Chaque molécule possède une carte d'identité affichant les informations suivantes:")
        with doc.create(Itemize()) as itemize:
            itemize.add_item(NoEscape(r"Le SMILES."))
            itemize.add_item(NoEscape(r"Nombres d'atomes et de liaisons (une double liaison comptant comme une liaison simple)."))
            itemize.add_item(NoEscape(r"Sa valeur de cagitude."))
            itemize.add_item(NoEscape(r"Une représentation en 3D de la molécule grâce à Pymol."))
            itemize.add_item(NoEscape(r"Son graphe des cycles et son graphe des coins."))
            itemize.add_item(NoEscape(r"Un lien amenant sur la page de la molécule sur la base de donnée."))
            

    
    # Ajouter la table des matières
    doc.append(NoEscape(r'\newpage'))
    # Ajouter la table des matières
    doc.append(NoEscape(r'\hypertarget{toc}{}'))
    doc.append(NoEscape(r'\tableofcontents'))
    j = 0
    cagitude = -1
    for mesure in mesures:
        doc.append(NoEscape(r'\newpage'))
        image_coin, image_cycle = generate_images(arg1, mesure.nom)
        image_mol = generate_image_mol(arg1,mesure.nom)
        images = [image_cycle, image_coin,image_mol]
        

        for image in images:
            if not Path(image).exists():
                print(f"Error: Image file {image} does not exist.")
                return

        if mesure.cagitude != cagitude:
            if nb_mesure[j] >1:
                section_title = f"Cagitude = {mesure.cagitude} ({nb_mesure[j]} molécules)"
            else:
                section_title = f"Cagitude = {mesure.cagitude} ({nb_mesure[j]} molécule)"
            j+=1
            with doc.create(Section(section_title)):
                cagitude = mesure.cagitude

        with doc.create(Subsection(f"{mesure.nom}'s ID")):
            # Ajouter un lien hypertexte pour revenir à la table des matières
            doc.append(NoEscape(r'\hyperlink{toc}{Retour à la table des matières}'))
            with doc.create(Itemize()) as itemize:
                doc.append(NoEscape(r'\hypertarget{' + mesure.nom + r'}{}'))
                smiles, num_atoms, num_bonds = get_molecule_info(arg1, mesure.nom)
                if smiles != None:
                    smiles = format_smiles(smiles)
                    itemize.add_item(NoEscape(r"SMILES: \begin{verbatim}" + smiles + r"\end{verbatim}"))
                itemize.add_item(f"Nombres d'atomes: {num_atoms}")
                itemize.add_item(f"Nombres de liaisons: {num_bonds}")
                value = str(mesure.cagitude)
                value2 = str(mesure.other_cagitude(csv_file2))
                itemize.add_item(f"Mesure connexe de cagitude = {value}; autre mesure = {value2}")
                if arg1 == "CHEBI":
                    mol = mesure.nom.replace("CHEBI_", "")
                    link = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:' + mol
                    itemize.add_item(NoEscape(r"Lien: \href{" + link + r"}{" + link + r"}"))

            with doc.create(Figure(position='H')) as fig:
                image_path_new = os.path.join('..', 'png_files_reduit', 'pymol', f'{mesure.nom}.jpg')
                fig.add_image(image_path_new, width=NoEscape(r'0.4\textwidth'))
                fig.add_caption("Molécule en 3D")

            with doc.create(Figure(position='H')) as graph:
                with doc.create(SubFigure(position='c', width=NoEscape(r'0.5\textwidth'))) as left_image:
                    image_path_cycles = os.path.join('..', 'png_files_reduit', 'graphes_cycles', f'{mesure.nom}.jpg')
                    left_image.append(NoEscape(r'\scalebox{1.0}{\includegraphics{'))
                    left_image.append(image_path_cycles)
                    left_image.append(NoEscape(r'}}'))
                    left_image.add_caption("Graphe des Cycles")

                with doc.create(SubFigure(position='c', width=NoEscape(r'0.5\textwidth'))) as right_image:
                    image_path_coins = os.path.join('..', 'png_files_reduit', 'graphes_coins', f'{mesure.nom}.jpg')
                    right_image.append(NoEscape(r'\scalebox{1.0}{\includegraphics{'))
                    right_image.append(image_path_coins)
                    right_image.append(NoEscape(r'}}'))
                    right_image.add_caption(NoEscape(r'\begin{minipage}{\linewidth}Graphe des Coins\end{minipage}'))


    doc.generate_tex(pdf_path)
    doc.generate_pdf(pdf_path, clean_tex=True)
    compress_pdf(pdf_path + '.pdf', compressed_path)
    optimize_pdf(compressed_path, opti_path)


if __name__ == "__main__":
    if len(sys.argv) != 3 and len(sys.argv) !=2:
        print("Une molécule: python generate_mol_id.py <base de donnée> <molécule>")
        print("Bestiaire des coins d'une base: python generate_mol_id.py <base de donnée>")

    elif len(sys.argv) == 3:
        arg1 = sys.argv[1]
        arg2 = sys.argv[2]
        pdf_path = generate_latex_document(arg1, arg2,False)
        if os.path.exists(pdf_path):
            subprocess.run(['xdg-open', pdf_path]) 
    elif len(sys.argv) == 2:
        arg1 = sys.argv[1]
        create_bestiaire(arg1)