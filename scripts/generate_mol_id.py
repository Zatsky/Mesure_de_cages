from pylatex import Document, Section, Subsection, Command, Figure,NewPage
from pylatex.utils import NoEscape, escape_latex
import os
import subprocess
from pathlib import Path
import sys

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

    return image_path_coin, image_path_cycle

def generate_latex_document(arg1, arg2):
    output_file = f'data/{arg1}/ID/{arg2}'
    doc = Document()

    image_coin, image_cycle = generate_images(arg1, arg2)

    images = [image_cycle, image_coin]

    # Vérifiez si les fichiers d'image existent
    for image in images:
        if not Path(image).exists():
            print(f"Error: Image file {image} does not exist.")
            return

    doc.preamble.append(Command('title', f"{arg2}'s ID"))
    doc.preamble.append(Command('author', arg1))
    doc.append(NoEscape(r'\maketitle'))

    with doc.create(Section('Introduction')):
        doc.append('This document is generated using Python and the pylatex library.')
        doc.append(' It includes text and images.')

    with doc.create(Section('Images')):
        # First image
        with doc.create(Subsection(f'{arg2}.png')):
            with doc.create(Figure(position='h!')) as graph:
                image_path_cycles = os.path.join('..', 'png_files_reduit', 'graphes_cycles', f'{arg2}.png')
                graph.add_image(image_path_cycles, width='0.8\\textwidth')
                graph.add_caption(f'{arg2}.png')
            doc.append(NewPage())

        # Second image
        with doc.create(Subsection(f'{arg2}.png')):
            with doc.create(Figure(position='h!')) as graph:
                image_path_coins = os.path.join('..', 'png_files_reduit', 'graphes_coins', f'{arg2}.png')
                graph.add_image(image_path_coins, width='0.8\\textwidth')
                graph.add_caption(f'{arg2}.png')
            doc.append(NewPage())

    doc.generate_tex(output_file)
    doc.generate_pdf(output_file, clean_tex=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_mol_id.py <arg1> <arg2>")
    else:
        arg1 = sys.argv[1]
        arg2 = sys.argv[2]
        generate_latex_document(arg1, arg2)
