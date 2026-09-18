# Optimisation du Réseau Vélib en Île-de-France

Projet de théorie des graphes et d'optimisation spatiale réalisé dans le cadre du BUT Informatique. 
L'objectif est d'analyser le réseau des 1 518 stations Vélib Métropole (Paris et Île-de-France), de construire un maillage spatial efficace avec la triangulation de Delaunay, d'optimiser l'interconnexion du réseau (MST) et d'offrir un calculateur d'itinéraire optimal (Dijkstra).

**Auteur** : Misbaou DIALLO (BUT 3 Informatique)

---

## Présentation et Objectifs

Le réseau Vélib d'Île-de-France comporte plus de 1 500 stations et 49 000 bornettes réparties sur 69 communes. Relier directement chaque station à toutes les autres créerait plus d'un million d'arêtes ($O(V^2)$), ce qui n'est ni réaliste ni efficace.

Ce projet applique plusieurs concepts d'algorithmique et de théorie des graphes :
1. **Triangulation de Delaunay** (`SciPy`) : Réduction des connexions candidates aux voisins géographiques directs (passage de ~1,15 million à **4 536 arêtes**).
2. **Coloration par Densité Spatiale** : Les triangles de la triangulation sont colorés selon leur superficie en échelle logarithmique. Plus la surface d'un triangle est petite (zone à forte densité de stations comme le centre de Paris), plus sa couleur est sombre.
3. **Arbre Couvrant Minimum (MST)** : Recherche de l'interconnexion minimale reliant toutes les stations sans cycle.
   - **Kruskal** (Union-Find)
   - **Prim** (Min-Heap)
4. **Calculateur d'Itinéraire (Dijkstra)** : Recherche du chemin le plus court entre deux stations sélectionnées dans la liste ou directement par clic sur la carte interactive.

---

## Aperçu du Dashboard & Graphiques

Le script génère deux livrables principaux :
- **`carte_velib_optimisee.html`** : Une carte interactive web (Leaflet / Folium) avec sélecteur d'itinéraire par clic, calque Delaunay activable/désactivable, et mini-dashboard Chart.js.
- **`data/graphiques_velib.png`** : Un tableau de bord analytique généré avec Matplotlib (distribution des capacités, répartition par commune, histogramme des distances inter-stations).

![Dashboard Analytique Vélib](data/graphiques_velib.png)

---

## Performances des Algorithmes

Les tests ont été effectués sur le jeu de données complet de **1 518 stations** :

| Algorithme | Usage | Complexité | Temps d'exécution |
|---|---|---|---|
| **Kruskal (Union-Find)** | Arbre Couvrant Minimum | $O(E \log E)$ | **~4.5 ms** |
| **Prim (Min-Heap)** | Arbre Couvrant Minimum | $O(E \log V)$ | **~52 ms** |
| **Dijkstra** | Plus court chemin entre 2 stations | $O((E + V) \log V)$ | **~4 ms** |

*Résultat MST* : Longueur totale minimale du réseau interconnecté = **502.09 km**.

---

## Installation et Exécution

### Prérequis
- Python 3.9+
- Dépendances : `scipy`, `pandas`, `numpy`, `folium`, `requests`, `matplotlib`

```bash
# 1. Cloner le dépôt
git clone https://github.com/misbaou672/Velib-Optim.git
cd Velib-Optim

# 2. Installer les bibliothèques requises
pip install -r requirements.txt

# 3. Lancer l'analyse et la génération de la carte
python optimisation_velib.py
```

L'application génèrera la carte interactive `carte_velib_optimisee.html` que vous pouvez ouvrir directement dans votre navigateur web.

---

## Structure du Projet

```text
Velib-Optim/
├── optimisation_velib.py          # Script principal (Data, Delaunay, MST, Dijkstra, Folium)
├── carte_velib_optimisee.html      # Application web interactive (Folium + Leaflet + Chart.js)
├── rapport_statistiques_velib.json # Rapport de synthèse au format JSON
├── plan.md                         # Roadmap et évolutions futures
├── anomalies_et_limites.md         # Rapport d'anomalies et points d'attention
├── data/
│   ├── stations_velib_idf_complete.json # Dataset des 1518 stations
│   └── graphiques_velib.png        # Tableau de bord analytique Matplotlib
└── README.md
```
