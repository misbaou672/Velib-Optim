# Cockpit Data & Optimisation Spatiale du Réseau Vélib Île-de-France

![Python](https://img.shields.io/badge/Python-3.9+-3776AB?style=for-the-badge&logo=python&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-Delaunay-8CAAE6?style=for-the-badge&logo=scipy&logoColor=white)
![Folium](https://img.shields.io/badge/Folium-Leaflet-77B829?style=for-the-badge&logo=leaflet&logoColor=white)
![Chart.js](https://img.shields.io/badge/Chart.js-Dashboard-FF6384?style=for-the-badge&logo=chart.js&logoColor=white)
![OpenData Paris](https://img.shields.io/badge/OpenData_Paris-API_Temps_Réel-0055A5?style=for-the-badge&logo=open-access&logoColor=white)
![Status](https://img.shields.io/badge/Statut-Opérationnel-10B981?style=for-the-badge)

---

> [!NOTE]
> **Cadre du projet** : Ce projet est issu d'un projet universitaire initial de **BUT Informatique**, continuellement enrichi, optimisé et transformé en un véritable **Cockpit Data Spatiale** sur mon temps libre par **Misbaou DIALLO**.

---

## Présentation du Projet

Le réseau **Vélib’ Métropole** en Île-de-France comprend **1 518 stations** et plus de **49 000 bornettes** réparties sur **69 communes**. 

Calculer toutes les connexions directes entre stations créerait plus de **1,15 million d'arêtes** ($O(V^2)$), ce qui engorgerait le réseau. Ce projet applique la théorie des graphes et l'analyse spatiale pour modéliser, optimiser et visualiser le réseau en temps réel.

### Fonctionnalités Clés du Cockpit Data
- **Connexion Live API OpenData Paris** : Synchronisation en temps réel des vélos disponibles (mécaniques vs électriques) et bornettes libres.
- **Triangulation de Delaunay (`SciPy`)** : Réduction du maillage spatial à **4 536 arêtes candidates** pondérées par la distance Haversine.
- **Coloration par Densité Spatiale** : Triangles colorés en échelle logarithmique (zones denses comme le centre de Paris en Midnight Navy sombre `#0F172A`).
- **Arbres Couvrants Minimums (MST)** : Algorithmes de **Kruskal** (*Union-Find par rang*) et **Prim** (*Min-Heap*) réduisant le réseau interconnecté global à **502,09 km**.
- **Calculateur d'Itinéraire Dijkstra Multi-Étapes** : Recherche du chemin le plus court avec étape optionnelle (*Waypoint*), temps estimé et **bilan d'économie de CO₂**.
- **Carte de Chaleur (HeatMap)** : Calque dynamique de densité de disponibilité des vélos.
- **Comparateur de Stations Côte-à-Côte** : Analyse comparative instantanée des capacités entre départ et arrivée.

---

## Performance des Algorithmes (Jeu de données : 1 518 stations)

| Algorithme | Usage & Objectif | Complexité | Temps d'Exécution |
|---|---|---|---|
| **Kruskal (Union-Find)** | Arbre Couvrant Minimum (MST) | $O(E \log E)$ | **~5.5 ms** |
| **Prim (Min-Heap)** | Validation Croisée MST | $O(E \log V)$ | **~49.0 ms** |
| **Dijkstra** | Plus Court Chemin Inter-Stations | $O((E + V) \log V)$ | **~4.0 ms** |

> [!TIP]
> **Gain de Calcul** : La triangulation de Delaunay permet d'économiser **99.6% des calculs de paires** tout en garantissant la présence de l'arbre couvrant minimum exact.

---

## Aperçu du Tableau de Bord

![Tableau de Bord Analytique Vélib](data/graphiques_velib.png)

---

## Installation et Exécution Local

```bash
# 1. Cloner le dépôt GitHub
git clone https://github.com/misbaou672/Velib-Optim.git
cd Velib-Optim

# 2. Installer les dépendances Python
pip install -r requirements.txt

# 3. Exécuter le générateur du Cockpit Data
python optimisation_velib.py
```

Ouvrez ensuite `carte_velib_optimisee.html` dans votre navigateur web pour accéder au Cockpit interactif.

---

## Architecture du Dépôt

```text
Velib-Optim/
├── optimisation_velib.py          # Moteur principal (Data API, Delaunay, Kruskal, Prim, Dijkstra, Folium)
├── carte_velib_optimisee.html      # Cockpit Data interactif web (Folium, Leaflet, Chart.js)
├── rapport_statistiques_velib.json # Rapport analytique au format JSON
├── plan.md                         # Roadmap & évolutions futures du projet
├── anomalies_et_limites.md         # Rapport technique d'architecture & vigilance
├── data/
│   ├── stations_velib_idf_complete.json # Dataset mis à jour en direct via l'API OpenData Paris
│   └── graphiques_velib.png        # Synthèse graphique Matplotlib
└── README.md
```

---

**Auteur** : **Misbaou DIALLO** *(BUT 3 Informatique)*  
Projet universitaire étendu et optimisé sur mon temps libre.
