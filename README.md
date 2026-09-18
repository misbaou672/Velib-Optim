<div align="center">

  # 🚲 Optimisation Algorithmique du Réseau Vélib Métropole (1 518 Stations - Île-de-France)

  [![Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://python.org)
  [![SciPy](https://img.shields.io/badge/SciPy-Delaunay-8CAAE6?style=for-the-badge&logo=scipy&logoColor=white)](#)
  [![Pandas](https://img.shields.io/badge/Pandas-Data_Analysis-150458?style=for-the-badge&logo=pandas&logoColor=white)](#)
  [![Folium](https://img.shields.io/badge/Folium-Interactive_Maps-77B829?style=for-the-badge&logo=leaflet&logoColor=white)](#)
  [![OpenData](https://img.shields.io/badge/Data-Data.Gouv.fr_/_IDF-0055A5?style=for-the-badge)](#)
  [![License](https://img.shields.io/badge/License-MIT-blue.svg?style=for-the-badge)](#)

  <p align="center">
    <b>Étude algorithmique, spatialisation géodésique et optimisation par Arbre Couvrant Minimum (MST) de l'intégralité des 1 518 stations Vélib métropolitaines réparties sur 69 communes d'Île-de-France.</b>
  </p>

---

</div>

## 📌 Présentation du Projet

Ce projet a pour objectif d'optimiser l'infrastructure réseau et électrique reliant l'ensemble des **1 518 stations Vélib** réparties dans la région Île-de-France (Paris et 68 communes de la Petite et Grande Couronne : 92, 93, 94, 78, 91, 95).

### 📐 Approche Algorithmique & Data

1. **Données Ouvertes en Temps Réel (`Data.Gouv.fr / OpenData Île-de-France`)** : Extraction des coordonnées GPS, capacités et codes INSEE de 1 518 stations.
2. **Triangulation de Delaunay (`SciPy`)** : Réduction de la complexité spatiale des liaisons de $O(V^2)$ ($\approx 1.15 \text{ million}$) à un graphe planaire $O(V)$ de **4 536 arêtes candidates**.
3. **Distances Géodésiques (Haversine)** : Calcul de la distance réelle sur ellipsoïde terrestre.
4. **Arbre Couvrant Minimum (MST)** : 
   - **Algorithme de Kruskal** (*Union-Find / Ensemble Disjoint*, $O(E \log E)$).
   - **Algorithme de Prim** (*Min-Heap / Priority Queue*, $O(E \log V)$).
5. **Visualisation Cartographique Interactive (`Folium & MarkerCluster`)** : Carte vectorielle réactive avec clustering dynamique et contrôle de couches.

---

## ⚡ Résultats & Performances

| Métrique | Valeur |
|---|---|
| **Nombre de Stations** | **1 518 stations** |
| **Couverture Géographique** | **69 communes d'Île-de-France** |
| **Arêtes Candidates (Delaunay)** | **4 536 connexions** |
| **Longueur Totale du Réseau Optimisé (MST)** | **502.09 km** |
| **Temps d'exécution (Kruskal)** | **5.63 ms** |
| **Temps d'exécution (Prim)** | **5.30 ms** |

---

## 🛠️ Stack Technique

- **Langage** : Python 3.10+
- **Analyse Spatiale & Données** : SciPy, NumPy, Pandas, Requests
- **Visualisation Cartographique** : Folium, MarkerCluster, OpenStreetMap
- **Algorithmes de Graphe** : Kruskal (Disjoint Set), Prim (Min-Heap), Triangulation de Delaunay

---

## 🚀 Lancement Rapide

```bash
# 1. Cloner le projet
git clone https://github.com/misbaou672/optimisation-velib.git
cd optimisation-velib

# 2. Installer les dépendances
pip install scipy pandas numpy folium requests

# 3. Exécuter l'optimisation et générer la carte
python optimisation_velib.py
```

L'exécution génère le fichier `carte_velib_optimisee.html` consultable directement dans votre navigateur web.

---

<div align="center">
  <sub>Développé par <a href="https://github.com/misbaou672">Misbaou DIALLO</a></sub>
</div>
