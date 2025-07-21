# Simulation

Ce script automatise le nettoyage, la compilation, l'exécution et la visualisation d'une simulation.

## Prérequis

- Python 3  
- bash  
- make  
- Dépendances Python : `matplotlib`, `numpy`, etc...

Un fichier "patch_centers.csv" contenant les coordonnées des particules.

## Installation

```bash

git clone https://github.com/CyrilBonus/cilia_simulation.git
cd cilia_simulation
conda create -n cil_sim python=3.10
conda activate cil_sim
pip install -r cil_sim_requirements.txt

```

## Utilisation

```bash

python3 simulation.py

```

Selectionner les particules pertinentes, puis fermez la fenêtre.

Entrez dans la console le nombre de temps souhaité pour la simulation, puis appuyez sur Entrée.
Cette valeur correspond au nombre de fichiers .dat générés, soit le nombre d’instantanés temporels simulés.

Le programme génère l’évolution des phases du système au cours du temps.

À la fin de la simulation, deux interfaces s’affichent automatiquement :

 - une visualisation 3D de la surface des phases,
 - et heatmap 2D synchronisée.
 - 
Dans l’interface 3D, un curseur vous permet de sélectionner l’instantané affiché dans les deux fenêtres.

