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
conda env create -f cil_sim_requirements.yml -n cil_sim
conda activate cil_sim

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

