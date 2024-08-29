from pylatex import Document, Section, Command, Figure, SubFigure, Package, Itemize, Subsection, StandAloneGraphic, Tabular, Hyperref,LongTable
from pylatex.base_classes.command import Options
from pylatex.utils import fix_filename,NoEscape
import os
import subprocess
from pathlib import Path
import sys
from pymol import cmd
import csv
from rdkit import Chem
from PIL import Image, ImageOps
import pikepdf
import generation_graphiques
from pathlib import Path
from openbabel import pybel
import subprocess
import ctypes

def get_4_cliques(arg1, arg2):
    path = "data/" + arg1 + "/graphes_coins/" + arg2 + ".csv"
    result = subprocess.run(['./find_4_cliques', path], capture_output=True, text=True)
    output = result.stdout.strip()
    return output


def cleanup_temp_files(directory):
    extensions = ['.aux', '.log', '.out', '.toc']
    for file in os.listdir(directory):
        if any(file.endswith(ext) for ext in extensions):
            os.remove(os.path.join(directory, file))


def generate_pdf_safe(doc, pdf_path):
    try:
        doc.generate_pdf(pdf_path, compiler='latexmk', clean_tex=False)
    except subprocess.CalledProcessError as e:
        output = e.output
        with open('error_log.txt', 'wb') as f:
            f.write(output)
        try:
            print(output.decode('utf-8', errors='ignore'))
        except UnicodeDecodeError:
            print(output.decode('latin1', errors='ignore'))
        raise

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
    def __init__(self,nom,cagitude,size):
        self.cagitude = cagitude
        self.nom = nom
        self.size = size
        self.nb_sim = 0

    def other_cagitude(self,path):
        with open(path, 'r') as fichier:
            lecteur_csv = csv.reader(fichier)
            for i,ligne in enumerate(lecteur_csv):
                if ligne[0] == self.nom:
                    return int(float(ligne[1]))

class tableau:
    def __init__(self,taille,tab, cagitude):
        self.taille = taille
        self.nb_mol = taille
        self.tab = tab
        self.cagitude = cagitude
        self.classes_cycles = [[] for _ in range(taille)]
        self.classes_coins = [[] for _ in range(taille)]

    def tri_tableau_cycles(self,arg1):

        x = 0
        ignorer = 0
        if self.taille == 1:
            self.classes_cycles[0].append(self.tab[0])
        else:
            lib = ctypes.CDLL('./similarite.so')
            lib.liberer_float.argtypes = [ctypes.POINTER(ctypes.c_float)]
            lib.liberer_float.restype = None
            lib.MCIS.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
            lib.MCIS.restype = ctypes.POINTER(ctypes.c_float)
            c_str_base = ctypes.create_string_buffer(arg1.encode('utf-8'))
            print(f"Comparaisons de la cagitude: {self.cagitude} ({self.taille} molécules)")
        

            while x< self.taille and self.taille >1:
                #print(x)
                self.classes_cycles[x].append(self.tab[x])
                #too_much = 0
                #print(f"taille = {self.taille} et x = {x}")
                copy_tab = list(self.tab)
                c_str_mol1 = ctypes.create_string_buffer(self.tab[x].nom.encode('utf-8'))
                j = 0
                sup = []
                total =0
                for mesure in copy_tab:
                    if self.tab[x].nom != mesure.nom:
                        c_str_mol2 = ctypes.create_string_buffer(mesure.nom.encode('utf-8'))

                        #print(f"{self.tab[x].nom} et {mesure.nom}")
                        nums = lib.MCIS(c_str_base,c_str_mol1,c_str_mol2)
                        nums_array = [nums[i] for i in range(2)]
                        total += nums_array[0]
                        if nums_array[0] == 1.0:
                            sup.append(j)
                        #elif nums_array[0] == 0.0:
                        #    too_much += 1
                        lib.liberer_float(nums)
                    j = j+1
                if total ==0:
                    ignorer +=1
                #print(f"nombre de supp = {len(sup)} et nb à 0  = {too_much}")
                for index in sorted(sup, reverse=True):
                    self.classes_cycles[x].append(self.tab.pop(index))
                self.tab[x].nb_sim = len(sup) +1
                self.taille = self.taille - len(sup)
                x = x +1
            print(f"nombre de molécules ignorées = {ignorer}")

    def tri_tableau_coins(self,arg1):
        lib = ctypes.CDLL('./similarite.so')
        lib.liberer_float.argtypes = [ctypes.POINTER(ctypes.c_float)]
        lib.liberer_float.restype = None
        lib.MCIS.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
        lib.MCIS.restype = ctypes.POINTER(ctypes.c_float)
        c_str_base = ctypes.create_string_buffer(arg1.encode('utf-8'))
        x = 0
        ignorer = 0
        copy_tab1 = list(self.tab)
        copy_taille = self.taille
        if self.taille == 1:
            self.classes_coins[0].append(self.tab[0])
        else:
            while x< copy_taille:
                self.classes_coins[x].append(copy_tab1[x])
                #too_much = 0
                #print(f"taille = {self.taille} et x = {x}")
                copy_tab = list(copy_tab1)
                c_str_mol1 = ctypes.create_string_buffer(copy_tab1[x].nom.encode('utf-8'))
                j = 0
                sup = []
                total =0
                for mesure in copy_tab:
                    if copy_tab1[x].nom != mesure.nom:
                        c_str_mol2 = ctypes.create_string_buffer(mesure.nom.encode('utf-8'))

                        #print(f"{self.tab[x].nom} et {mesure.nom}")
                        nums = lib.MCIS(c_str_base,c_str_mol1,c_str_mol2)
                        nums_array = [nums[i] for i in range(2)]
                        total += nums_array[1]
                        if nums_array[1] == 1.0:
                            sup.append(j)
                        #elif nums_array[0] == 0.0:
                        #    too_much += 1
                        lib.liberer_float(nums)
                    j = j+1
                if total ==0:
                    ignorer +=1
                #print(f"nombre de supp = {len(sup)} et nb à 0  = {too_much}")
                for index in sorted(sup, reverse=True):
                    self.classes_coins[x].append(copy_tab1.pop(index))
                copy_taille = copy_taille - len(sup)
                x = x +1
            print(f"nombre de molécules ignorées = {ignorer}")
            z=0
            id = 0
            #print(f"longu = {len(self.classes_coins)}")
            while z < len(self.classes_coins):
                y = 0
                #print(f"longu = {len(self.classes_coins[z])}")
                while y < len(self.classes_coins[z]):
                    self.tab[id] = self.classes_coins[z][y]
                    id +=1
                    y+=1
                z+=1


def create_tableau(doc,mesures,arg1):
    total = 0
    id = 0
    coins = 0
    with doc.create(LongTable("|p{4cm}|p{4cm}|p{8cm}|")) as table:
        table.add_hline()
        table.add_row(("Représentant de la classe", "Nombres de molécules dans la classe", "Molécules dans la classe"))
        table.add_hline()
        table.end_table_header()
        for mesure in mesures.tab:
            if len(mesures.classes_coins[coins])>0:
                if mesure.nom == mesures.classes_coins[coins][0].nom:
                    table.add_hline()
                    coins +=1
            nom = NoEscape(r'\hyperlink{' + mesure.nom + r'}{'+mesure.nom.replace("_", r"\_")+'}')
            table.add_row((nom, f"{len(mesures.classes_cycles[id])}", liste_classe(mesures,arg1,id)))
            table.add_hline()
            total+=len(mesures.classes_cycles[id])
            id += 1
        table.add_row(("Total",f"{total}",""))
        table.add_hline()


def create_tableau_restrictions(doc,restrictions):
    restant = restrictions[0]
    with doc.create(Tabular('|c|c|c|')) as table:
        table.add_hline()
        table.add_row(("Restriction", "Nombre de molécules retirées", "Nombre de molécules restantes"))
        table.add_hline()
        table.add_hline()
        table.add_row(("Fichiers.mol", str(0), str(restrictions[0])))
        table.add_hline()
        restant -= restrictions[1]
        table.add_row(("Seg fault", str(restrictions[1]), str(restant)))
        table.add_hline()
        restant -= restrictions[2]
        table.add_row((">250 atomes", str(restrictions[2]), str(restant)))
        table.add_hline()
        restant -= restrictions[3]
        table.add_row(("Aucun coins", str(restrictions[3]), str(restant)))
        table.add_hline()
        table.add_row(("Total final",f"{str(restrictions[0]-restant)}",f"{restant}"))
        table.add_hline()

def liste_classe(mesures, arg1, id):
    noms = ""
    for mesure in mesures.classes_cycles[id]:
        if noms != "":
            noms += ", "
        if arg1 == "CHEBI":
            mol = mesure.nom.replace("CHEBI_", "")
            link = "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{}".format(mol)
            noms += r"\href{{{}}}{{{}}}".format(link, mesure.nom.replace("_", r"\_"))
        elif arg1 == "LOTUS":
            link = "https://lotus.naturalproducts.net/compound/lotus_id/{}".format(mesure.nom)
            noms += r"\href{{{}}}{{{}}}".format(link, mesure.nom.replace("_", r"\_"))
        else:
            noms += mesure.nom.replace("_", r"\_")

    return NoEscape(noms)

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
    mes = []
    cagitude = 0
    j = 0
    x = 0
    with open(nom_fichier, 'r') as fichier:
        lecteur_csv = csv.reader(fichier)
        for i,ligne in enumerate(lecteur_csv):
            try:
                cage = int(float(ligne[1]))
                #print(str(cage))
                if get_arete_coins(arg1,ligne[0]) <500:
                    mesure = donnee_csv(ligne[0],cage,True)
                    
                else:
                    mesure = donnee_csv(ligne[0],cage,False)
                    x+=1
                mesures.append(mesure)
            except ValueError:
                pass
    mesures.sort(key=lambda z: z.cagitude)
    cagitude = mesures[0].cagitude
    for i, mesure in enumerate(mesures):   
        if cagitude != mesure.cagitude:
            nb_mesures.append(tableau(j,list(mes),cagitude))
            mes.clear()
            cagitude = mesure.cagitude
            j = 0
        mes.append(mesure)
        j +=1
    nb_mesures.append(tableau(j,list(mes),cagitude))
    print("Il y'a "+ str(len(nb_mesures))+ " mesures différentes et "+ str(x) + "supérieures à 250")
    return mesures,nb_mesures,x

#Renvoie le smiles, le nombres d'atomes, ainsi que le nombres de liaisons d'une molécule
def get_molecule_info(arg1,arg2):
    mol_file = 'data/' + arg1 + '/mol_files/'+ arg2 +'.mol'
    mol = next(pybel.readfile("mol", mol_file))

    # Extraire les informations
    smiles = mol.write("can").strip()  # Utiliser le format canonique pour les SMILES
    num_atoms = mol.OBMol.NumAtoms()
    num_bonds = mol.OBMol.NumBonds()
    return smiles, num_atoms, num_bonds 
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
    cmd.png(path, width=930, height=650, dpi=300, ray=1)
    cmd.delete('all')

    convert_and_resize_image(path,path.replace('.png', '.jpg'))
    
    return path.replace('.png', '.jpg')

def escape_latex_special_chars(smiles):
    """Échappe les caractères spéciaux pour LaTeX dans une chaîne SMILES."""
    return smiles.replace('_', r'\_').replace('&', r'\&').replace('%', r'\%').replace('$', r'\$').replace('#', r'\#').replace('{', r'\{').replace('}', r'\}')

def format_smiles(smiles, line_length=118):
    """Formate une chaîne SMILES en insérant des sauts de ligne pour éviter qu'elle ne dépasse les marges."""
    escaped_smiles = escape_latex_special_chars(smiles)
    return '\n'.join([escaped_smiles[i:i+line_length] for i in range(0, len(escaped_smiles), line_length)])

def convert_and_resize_image(image_path, output_path, quality=100, target_size=(930, 650)):
    # Ouvrir l'image
    with Image.open(image_path) as img:
        # Convertir en JPEG avec la qualité spécifiée
        img = img.convert("RGB")
        
        # Calculer le ratio de l'image originale et du nouveau cadre
        original_ratio = img.width / img.height
        target_ratio = target_size[0] / target_size[1]

        if original_ratio > target_ratio:
            # L'image est plus large par rapport à la cible
            new_width = target_size[0]
            new_height = int(target_size[0] / original_ratio)
        else:
            # L'image est plus haute par rapport à la cible
            new_width = int(target_size[1] * original_ratio)
            new_height = target_size[1]

        # Redimensionner l'image tout en conservant le ratio
        img = img.resize((new_width, new_height), Image.LANCZOS)
        
        # Créer une nouvelle image avec un fond blanc
        new_img = Image.new("RGB", target_size, (255, 255, 255))
        
        # Calculer la position pour centrer l'image redimensionnée
        left = (target_size[0] - new_width) // 2
        top = (target_size[1] - new_height) // 2
        new_img.paste(img, (left, top))
        
        # Enregistrer l'image finale
        new_img.save(output_path, "JPEG", quality=quality)
        print(f"Converted and saved image: {output_path}")

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

def get_infos_csv(arg1, arg2):
    cycle_csv = f'data/{arg1}/graphes_cycles/{arg2}.csv'
    coin_csv = f'data/{arg1}/graphes_coins/{arg2}.csv'
    with open(cycle_csv, mode='r', newline='') as fichier:
        lecteur_csv = csv.reader(fichier)
        premiere_ligne = next(lecteur_csv)[0]  # Lire la première ligne
        nombre1, nombre2 = premiere_ligne.split('_')  # Séparer les nombres
        cycles = [nombre1,nombre2]
    with open(coin_csv, mode='r', newline='') as fichier:
        lecteur_csv = csv.reader(fichier)
        premiere_ligne = next(lecteur_csv)[0]  # Lire la première ligne
        nombre1, nombre2 = premiere_ligne.split('_')  # Séparer les nombres
        coins = [nombre1,nombre2]
    return cycles,coins

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
    doc.packages.append(Package('geometry', options=['left=2cm', 'right=2cm']))
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
    doc.packages.append(NoEscape(r'\usepackage{newunicodechar}'))
    doc.preamble.append(NoEscape(r'\newunicodechar{⁹}{\textsuperscript{9}}'))
    doc.preamble.append(NoEscape(r'\newunicodechar{²}{\textsuperscript{2}}'))
    doc.preamble.append(NoEscape(r'\newunicodechar{³}{\textsuperscript{3}}'))
    

    if arg1 == "CHEBI":
        seg_fault = 6
    elif arg1 == "LOTUS":
        seg_fault = 27
    else:
        seg_fault = 0
    all_files = count_files_in_directory(arg1)

    mesures,nb_mesure,trop_gros = lire_donnees_csv(csv_file,arg1)
    too_big = count_mol_files_with_more_than_250_atoms(arg1)
    total_final = 0
    for tab in nb_mesure:
        
        tab.tri_tableau_cycles(arg1)
        tab.tri_tableau_coins(arg1)
        total_final += len(tab.tab)
    no_corner = all_files - seg_fault - too_big -len(mesures) - trop_gros
    restrictions =[all_files,seg_fault,too_big,no_corner]
    nb_comparaisons = len(mesures) - total_final
    # Ajouter le titre principal
    doc.preamble.append(Command('title', 'Bestiaire des coins de '+arg1))
    doc.preamble.append(Command('date', NoEscape(r'\today')))
    doc.preamble.append(Command('author', 'Docherty Ronan'))
    doc.append(NoEscape(r'\maketitle'))
    if not Path(f'data/{arg1}/results/histogramme_discret.png').exists():
        generation_graphiques.gen_graphique(arg1)
    graph_path = f'../results/histogramme_discret.png'
    with doc.create(Section('Introduction', numbering=False)):
        doc.append("Ce document regroupe l'ensemble des molécules de la base " + arg1.replace("_", r"\_") + " possédant au moins un coin et ayant moins de 250 atomes.")
        doc.append(NoEscape(r'\begin{figure}[h!]'))
        doc.append(NoEscape(r'\centering'))
        doc.append(NoEscape(r'\includegraphics[width=0.8\textwidth]{' + graph_path + r'}'))
        doc.append(NoEscape(r'\caption{Histogramme discret de la cagitude ChEBI}'))
        doc.append(NoEscape(r'\end{figure}'))
    doc.append(NoEscape(r'\newpage'))
    
    pluriel = ["cette","",""]
    with doc.create(Section(('Traitement des données'), numbering=False)):
        with doc.create(Subsection(('Contraintes'), numbering=False)):
            doc.append(NoEscape(r"Après avoir récupéré les fichiers .mol depuis la base de donnée "+arg1.replace("_", r"\_") +r", nous devons appliquer des restrictions pour les molécules non-éligibles à ce bestiaire de coins. Tout d'abord certains des fichiers .mol extraits causent des problèmes pendant leur lecture et font planter le programme (seg fault). Ensuite, la complexité de l'algorithme de Vismara, nous permettant de générer les graphes des cycles, est dépendante du nombre d'arêtes dans le graphe moléculaire ce qui peut le rendre impossible à calculer pour des graphes trop denses. Pour éviter cela, nous avons décidé de retirer les molécules de plus de 250 atomes. Comme nous voulons lister seulement les molécules contenant des coins, nous avons retiré toutes celles n'en contenant pas. Le résultat de ces restrictions est affiché dans le tableau ci-dessous: "))
            doc.append(NoEscape(r"\\ \\"))
            create_tableau_restrictions(doc,restrictions)
        with doc.create(Section('Représentation des données', numbering=False)):
            doc.append(NoEscape(f"Afin d'éviter les doublons de molécules, nous effectuons l'algorithme du MCIS (Maximum Common Induced Subgraph) sur les graphes de cycles des ensembles de molécules possédant la même mesure de cagitude d'arête. Lorsque la mesure de MCIS entre les graphes de cycles est égale à 1 (et qu'elle possède donc le même graphe de cycles) nous plaçons ces molécules dans la même classe d'équivalence. Nous choisissons ensuite une molécule de cette classe pour représenter toutes les autres en indiquant combien de molécules elle représente. Cela nous a permis de représenter {len(mesures)} molécules avec {total_final} classes d'équivalences."))
            doc.append(NoEscape(r"\\"))
            doc.append(NoEscape("Nous avons généré un tableau pour chaque mesure de cagitude (présent à la première page des ensembles de molécules pour une mesure), détaillant les classes d'équivalences. Ce tableau indique quelle molécule représente sa classe d'équivalence (avec un lien permettant d'accéder à sa page), le nombre de molécules dans cette classe, ainsi que toutes les molécules lui appartenant (elles aussi ont un lieu amenant cette fois-ci au site de la base de données). Enfin, ces classes d'équivalences sont séparées par une ligne lorsque qu'on change de classe d'équivalence pour les graphes de cycles et par deux lignes lorsqu'on change de classe d'équivalence de graphes de coins."))
            doc.append(NoEscape(r"\\"))
            doc.append(NoEscape(r"\\"))
            doc.append("Chaque classe d'équivalence possède une carte d'identité moléculaire affichant les informations suivantes:")
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
    for mesured in nb_mesure:
        doc.append(NoEscape(r'\newpage'))
        if mesured.taille> 1:
            pluriel1 = "s"
        else:
            pluriel1=""
        if mesured.nb_mol >1:
            pluriel2 = "s"
        else:
            pluriel2 = ""
        section_title = f"Cagitude = {mesured.cagitude} ({mesured.taille} Classe{pluriel1} d'équivalence) ({mesured.nb_mol} molécule{pluriel2})"
        with doc.create(Section(section_title)):
            create_tableau(doc,mesured,arg1)
        for mesure in mesured.tab:
            doc.append(NoEscape(r'\newpage'))
                
            image_mol = generate_image_mol(arg1,mesure.nom)
            images = [image_mol]
            if mesure.size:
                image_coin, image_cycle = generate_images(arg1, mesure.nom)
                images.append(image_coin)
                images.append(image_cycle)
            

            for image in images:
                if not Path(image).exists():
                    print(f"Error: Image file {image} does not exist.")
                    return
            sub = f"{mesure.nom}"
            if mesure.nb_sim >0:
                sub = sub + f" ({str(mesure.nb_sim)} molécules)"
            with doc.create(Subsection(sub)):
                # Ajouter un lien hypertexte pour revenir à la table des matières
                doc.append(NoEscape(r'\hyperlink{toc}{Retour à la table des matières}'))
                doc.append(NoEscape(r'\hypertarget{' + mesure.nom + r'}{}'))
                with doc.create(Itemize()) as itemize:
                    smiles, num_atoms, num_bonds = get_molecule_info(arg1, mesure.nom)
                    cycle,coin=get_infos_csv(arg1,mesure.nom)
                    if smiles != None:
                        smiles = format_smiles(smiles)
                        itemize.add_item(NoEscape(r"{\footnotesize SMILES: \begin{verbatim}" + smiles + r"\end{verbatim}}"))
                    itemize.add_item(r"Molécules: "+ str(num_atoms)+r" atomes et "+str(num_bonds)+r" liaisons")
                    itemize.add_item(r"Graphes de cycles: "+cycle[0]+r" sommets et "+cycle[1]+r" arêtes")
                    som = r" sommets "
                    if int(coin[0])<2:
                        som = r" sommet "
                    are = r" arêtes"
                    if int(coin[1])<2:
                        are = r" arête"
                    itemize.add_item(r"Graphes de coins: "+coin[0]+som+r"et "+coin[1]+are)
                    value = str(mesure.cagitude)
                    value2 = str(mesure.other_cagitude(csv_file2))
                    #itemize.add_item(r"Nombres de cliques de taille 4: "+ str(get_4_cliques(arg1,mesure.nom)))
                    itemize.add_item(r"Mesure de cagitude de sommet = "+ str(value2) + r"; Mesure de cagitude d'arêtes= " + str(value))
                    if arg1 == "CHEBI":
                        mol = mesure.nom.replace("CHEBI_", "")
                        link = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:' + mol
                        itemize.add_item(NoEscape(r"Lien: \href{" + link + r"}{" + link + r"}"))
                    elif arg1 == "LOTUS":
                        link = 'https://lotus.naturalproducts.net/compound/lotus_id/' + mesure.nom
                    

                with doc.create(Figure(position='H')) as fig:
                    image_path_new = os.path.join('..', 'png_files_reduit', 'pymol', f'{mesure.nom}.jpg')
                    fig.add_image(image_path_new, width=NoEscape(r'0.6\textwidth'))
                    fig.add_caption("Molécule en 3D")

                if mesure.size:
                    with doc.create(Figure(position='H')) as graph:
                        with doc.create(SubFigure(position='b', width=NoEscape(r'0.5\textwidth'))) as left_image:
                            image_path_cycles = os.path.join('..', 'png_files_reduit', 'graphes_cycles', f'{mesure.nom}.jpg')
                            left_image.add_image(image_path_cycles, width=NoEscape(r'\textwidth'))
                            left_image.add_caption("Graphes des Cycles")

                        with doc.create(SubFigure(position='b', width=NoEscape(r'0.5\textwidth'))) as right_image:
                            image_path_coins = os.path.join('..', 'png_files_reduit', 'graphes_coins', f'{mesure.nom}.jpg')
                            right_image.add_image(image_path_coins, width=NoEscape(r'\textwidth'))
                            right_image.add_caption("Graphes des Coins")


    doc.generate_tex(pdf_path)
    #generate_pdf_safe(doc, pdf_path)
    doc.generate_pdf(pdf_path, clean_tex=False)
    compress_pdf(pdf_path + '.pdf', compressed_path)
    optimize_pdf(compressed_path, opti_path)

def create_doc_comparaison(arg1):

    
    csv_file =f'data/{arg1}/results/liste_mesure_alpha_connexe.csv'
    csv_file2 =f'data/{arg1}/results/liste_mesure_alpha.csv'
    lib = ctypes.CDLL('./similarite.so')

    # Déclarer le prototype de la fonction
    lib.MCIS.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
    lib.MCIS.restype = ctypes.POINTER(ctypes.c_float)
    input_str_base = arg1
    c_str_base = ctypes.create_string_buffer(input_str_base.encode('utf-8'))
    
    

    mesures,nb_mesure,ignore = lire_donnees_csv(csv_file,arg1)

    x = 0
    cagitude = -1
    for tab_mesure in nb_mesure:
        if tab_mesure.taille > 4 and tab_mesure.cagitude != 0: 
            for mesure in tab_mesure.tab:
                image_mol = generate_image_mol(arg1,mesure.nom)
                images = [image_mol]
                if mesure.size:
                    image_coin, image_cycle = generate_images(arg1, mesure.nom)
                    images.append(image_coin)
                    images.append(image_cycle)
                

                for image in images:
                    if not Path(image).exists():
                        print(f"Error: Image file {image} does not exist.")
                        return
            # Ajouter le titre principal
            opti_path = f'data/{arg1}/ID/Comparaisons/Comparaison_{arg1}_cagitude_{tab_mesure.cagitude}.pdf'
            pdf_path = f'data/{arg1}/ID/Comparaisons/Comparaison_{arg1}_cagitude_{tab_mesure.cagitude}'
            doc = Document(documentclass='article', document_options='a4paper')
            
            # Ajouter les packages nécessaires
            doc.packages.append(Command('usepackage', 'hyperref'))
            doc.packages.append(Package('geometry', options=['left=2cm', 'right=2cm']))
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
            doc.packages.append(NoEscape(r'\usepackage{newunicodechar}'))
            doc.preamble.append(NoEscape(r'\newunicodechar{⁹}{\textsuperscript{9}}'))
            doc.preamble.append(NoEscape(r'\newunicodechar{²}{\textsuperscript{2}}'))
            doc.preamble.append(NoEscape(r'\newunicodechar{³}{\textsuperscript{3}}'))
            doc.preamble.append(Command('title', 'Comparaison molécules de cagitude '+ str(tab_mesure.cagitude) +" " + arg1))
            doc.preamble.append(Command('date', NoEscape(r'\today')))
            doc.preamble.append(Command('author', 'Docherty Ronan'))
            doc.append(NoEscape(r'\maketitle'))
            with doc.create(Section('Introduction', numbering=False)):
                doc.append("Ce document regroupe les comparaisons par MCIS de l'ensemble des molécules de la base " + arg1 + " ayant une mesure de cagitude d'arête de "+ str(tab_mesure.cagitude) + ".")
            doc.append(NoEscape(r'\newpage'))
            
            
            # Ajouter la table des matières
            doc.append(NoEscape(r'\newpage'))
            # Ajouter la table des matières
            doc.append(NoEscape(r'\hypertarget{toc}{}'))
            doc.append(NoEscape(r'\tableofcontents'))
            j = 0
            cagitude = -1
            for mesure in tab_mesure.tab:

                doc.append(NoEscape(r'\newpage'))
                
                doc.append(NoEscape(r'\hypertarget{'+escape_latex_special_chars(mesure.nom)+r'}{}'))
                with doc.create(Section(f"{mesure.nom}")):
                    # Ajouter un lien hypertexte pour revenir à la table des matières
                    doc.append(NoEscape(r'\hyperlink{toc}{Retour à la table des matières}'))
                    with doc.create(Itemize()) as itemize:
                        smiles, num_atoms, num_bonds = get_molecule_info(arg1, mesure.nom)
                        if smiles != None:
                            smiles = format_smiles(smiles)
                            itemize.add_item(NoEscape(r"{\footnotesize SMILES: \begin{verbatim}" + smiles + r"\end{verbatim}}"))
                        itemize.add_item(r"Nombres d'atomes: "+ str(num_atoms))
                        itemize.add_item(r"Nombres de liaisons: "+ str(num_bonds))
                        value = str(mesure.cagitude)
                        value2 = str(mesure.other_cagitude(csv_file2))
                        #itemize.add_item(r"Nombres de cliques de taille 4: "+ str(get_4_cliques(arg1,mesure.nom)))
                        itemize.add_item(r"Mesure de cagitude de sommet = "+ str(value2) + r"; Mesure de cagitude d'arêtes= " + str(value))
                        if arg1 == "CHEBI":
                            mol = mesure.nom.replace("CHEBI_", "")
                            link = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:' + mol
                            itemize.add_item(NoEscape(r"Lien: \href{" + link + r"}{" + link + r"}"))
                        elif arg1 == "LOTUS":
                            link = 'https://lotus.naturalproducts.net/compound/lotus_id/' + mesure.nom
                            itemize.add_item(NoEscape(r"Lien: \href{" + link + r"}{" + link + r"}"))

                    with doc.create(Figure(position='H')) as fig:
                        image_path_new = os.path.join('..','..', 'png_files_reduit', 'pymol', f'{mesure.nom}.jpg')
                        fig.add_image(image_path_new, width=NoEscape(r'0.6\textwidth'))
                        fig.add_caption("Molécule en 3D")

                    
                    with doc.create(Figure(position='H')) as graph:
                            with doc.create(SubFigure(position='b', width=NoEscape(r'0.5\textwidth'))) as left_image:
                                image_path_cycles = os.path.join('..','..', 'png_files_reduit', 'graphes_cycles', f'{mesure.nom}.jpg')
                                left_image.add_image(image_path_cycles, width=NoEscape(r'\textwidth'))
                                left_image.add_caption("Graphes des Cycles")

                            with doc.create(SubFigure(position='b', width=NoEscape(r'0.5\textwidth'))) as right_image:
                                image_path_coins = os.path.join('..','..', 'png_files_reduit', 'graphes_coins', f'{mesure.nom}.jpg')
                                right_image.add_image(image_path_coins, width=NoEscape(r'\textwidth'))
                                right_image.add_caption("Graphes des Coins")

                            
                input_str_mol1 = mesure.nom
                c_str_mol1 = ctypes.create_string_buffer(input_str_mol1.encode('utf-8'))
                for mesure in nb_mesure[x].tab:
                    if nb_mesure[x].tab[j].nom != mesure.nom:
                        input_str_mol2 = mesure.nom
                        #print("MCIS entre " +str(nb_mesure[x].tab[j].nom) + " et " + str(mesure.nom))
                        c_str_mol2 = ctypes.create_string_buffer(input_str_mol2.encode('utf-8'))
                        nums = lib.MCIS(c_str_base,c_str_mol1,c_str_mol2)
                        doc.append(NoEscape(r'\newpage'))
                        hyperlink_1 = f"\\hyperlink{{{nb_mesure[x].tab[j].nom}}}{{{nb_mesure[x].tab[j].nom}}}"
                        hyperlink_2 = f"\\hyperlink{{{mesure.nom}}}{{{mesure.nom}}}"
                        hyperlink_1 = escape_latex_special_chars(hyperlink_1)
                        hyperlink_2 = escape_latex_special_chars(hyperlink_2)
                        section_title = rf"Comparaison de {hyperlink_1} avec {hyperlink_2}"
                        with doc.create(Subsection(NoEscape(section_title))):
                            doc.append(f"MCIS des graphes de cycles = {nums[0]}")
                            with doc.create(Figure(position='H')) as graph:
                                with doc.create(SubFigure(position='b', width=NoEscape(r'0.5\textwidth'))) as left_image:
                                    image_path_cycles = os.path.join('..','..', 'png_files_reduit', 'graphes_cycles', f'{nb_mesure[x].tab[j].nom}.jpg')
                                    left_image.add_image(image_path_cycles, width=NoEscape(r'\textwidth'))
                                    left_image.add_caption("Graphes des Cycles " + nb_mesure[x].tab[j].nom)

                                with doc.create(SubFigure(position='b', width=NoEscape(r'0.5\textwidth'))) as right_image:
                                    image_path_coins = os.path.join('..','..', 'png_files_reduit', 'graphes_cycles', f'{mesure.nom}.jpg')
                                    right_image.add_image(image_path_coins, width=NoEscape(r'\textwidth'))
                                    right_image.add_caption("Graphes des Cycles " + mesure.nom)
                            doc.append(f"MCIS des graphes de coins = {nums[1]}")
                            with doc.create(Figure(position='H')) as graph:
                                with doc.create(SubFigure(position='b', width=NoEscape(r'0.5\textwidth'))) as left_image:
                                    image_path_cycles = os.path.join('..','..', 'png_files_reduit', 'graphes_coins', f'{nb_mesure[x].tab[j].nom}.jpg')
                                    left_image.add_image(image_path_cycles, width=NoEscape(r'\textwidth'))
                                    left_image.add_caption("Graphes des Coins " + nb_mesure[x].tab[j].nom)
                                    
                                with doc.create(SubFigure(position='b', width=NoEscape(r'0.5\textwidth'))) as right_image:
                                    image_path_coins = os.path.join('..','..', 'png_files_reduit', 'graphes_coins', f'{mesure.nom}.jpg')
                                    right_image.add_image(image_path_coins, width=NoEscape(r'\textwidth'))
                                    right_image.add_caption("Graphes des Coins " + mesure.nom)
                            
                        
                j = j +1

            doc.generate_tex(pdf_path)
            #generate_pdf_safe(doc, pdf_path)
            doc.generate_pdf(pdf_path, clean_tex=True)
        x = x +1

if __name__ == "__main__":
    if len(sys.argv) != 3 and len(sys.argv) !=2:
        print("Une molécule: python generate_mol_id.py <base de donnée> <molécule>")
        print("Bestiaire des coins d'une base: python generate_mol_id.py <base de donnée>")
        #result = subprocess.run(['./similarite', "CHEBI" "CHEBI_140980", "CHEBI_40611"], capture_output=True, text=True)
        #output = result.stdout.strip()
        #print(str(output[0]) +""+ str(output[1]))
        create_doc_comparaison("CHEBI")

    elif len(sys.argv) == 3:
        arg1 = sys.argv[1]
        arg2 = sys.argv[2]
        pdf_path = generate_latex_document(arg1, arg2,False)
        if os.path.exists(pdf_path):
            subprocess.run(['xdg-open', pdf_path]) 
    elif len(sys.argv) == 2:
        arg1 = sys.argv[1]
        create_bestiaire(arg1)