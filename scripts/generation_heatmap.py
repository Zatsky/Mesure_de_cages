import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix

# Charger le fichier CSV
file_path = 'classement_diff_molecules.csv'
data = pd.read_csv(file_path)

# Extraire les deux classements
ranking1 = data['classement_1']
ranking2 = data['classement_2']

# Calculer la matrice de confusion
conf_matrix = confusion_matrix(ranking1, ranking2)

# Transformer la matrice de confusion en DataFrame pour une meilleure visualisation
conf_matrix_df = pd.DataFrame(conf_matrix, index=sorted(set(ranking1)), columns=sorted(set(ranking2)))

# Créer une heatmap à partir de la matrice de confusion
plt.figure(figsize=(10, 8))
sns.heatmap(conf_matrix_df, annot=True, cmap='viridis')
plt.title('Heatmap de la comparaison des classements')
plt.xlabel('Classement 2')
plt.ylabel('Classement 1')
plt.show()
