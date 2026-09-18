<div align="center">

  # 🚲 Optimisation Algorithmique du Réseau Vélib

  [![Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://python.org)
  [![SciPy](https://img.shields.io/badge/SciPy-Delaunay-8CAAE6?style=for-the-badge&logo=scipy&logoColor=white)](#)
  [![Pandas](https://img.shields.io/badge/Pandas-Data_Analysis-150458?style=for-the-badge&logo=pandas&logoColor=white)](#)
  [![Folium](https://img.shields.io/badge/Folium-Interactive_Maps-77B829?style=for-the-badge&logo=leaflet&logoColor=white)](#)
  [![License](https://img.shields.io/badge/License-MIT-blue.svg?style=for-the-badge)](#)

  <p align="center">
    <b>Étude algorithmique, spatialisation géodésique et optimisation par Arbre Couvrant Minimum (MST) du réseau des stations Vélib Métropole (Paris & Petite Couronne).</b>
  </p>

---

</div>

## 📌 Présentation du Projet

Ce projet a pour objectif d'optimiser l'infrastructure réseau et électrique reliant les stations Vélib de la métropole parisienne en réduisant au minimum la longueur totale de câblage tout en garantissant la connectivité globale.

### 📐 Approche Algorithmique

1. **Triangulation de Delaunay (`SciPy`)** : Réduction du graphe complet $O(V^2)$ à un sous-graphe planaire $O(V)$ de connexions physiquement pertinentes.
2. **Distances Géodésiques (Haversine)** : Calcul de la distance sur sphère terrestre entre les coordonnées GPS ($\text{Lat}, \text{Lon}$).
3. **Arbre Couvrant Minimum (MST)** : 
   - **Algorithme de Kruskal** (Union-Find / Ensemble Disjoint, complexité $O(E \log E)$).
   - **Algorithme de Prim** (File de priorité Min-Heap, complexité $O(E \log V)$).
4. **Visualisation Cartographique Interactive (`Folium`)** : Génération d'une carte interactive HTML avec tracé vectoriel du réseau optimisé.

---

## ⚡ Résultats & Performance

| Algorithme | Longueur totale du réseau | Complexité Temporelle |
|---|---|---|
| **Kruskal (Union-Find)** | **23.45 km** | $O(E \log E)$ |
| **Prim (Min-Heap)** | **23.45 km** | $O(E \log V)$ |

---

## 🛠️ Stack Technique

- **Langage** : Python 3.10+
- **Analyse Spatiale & Données** : SciPy, NumPy, Pandas
- **Visualisation Cartographique** : Folium, OpenStreetMap
- **Algorithmes de Graphe** : Kruskal (Disjoint Set), Prim (Min-Heap), Triangulation de Delaunay

---

## 🚀 Lancement Rapide

```bash
# 1. Cloner le projet
git clone https://github.com/misbaou672/optimisation-velib.git
cd optimisation-velib

# 2. Installer les dépendances
pip install scipy pandas numpy folium

# 3. Exécuter l'optimisation et générer la carte
python optimisation_velib.py
```

L'exécution génère un fichier `carte_velib_optimisee.html` consultable directement dans votre navigateur web.

---

<div align="center">
  <sub>Développé par <a href="https://github.com/misbaou672">Misbaou DIALLO</a></sub>
</div>
