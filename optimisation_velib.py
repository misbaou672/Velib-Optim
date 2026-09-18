"""
==============================================================================
Optimisation Algorithmique du Réseau Vélib (Paris & Petite Couronne)
Auteur : Misbaou DIALLO (BUT 3 Informatique)
==============================================================================

Ce projet implémente une étude algorithmique et cartographique pour optimiser
le réseau d'interconnexion (électrique / télécom) des stations Vélib métropolitaines.

Fonctionnalités :
  1. Triangulation de Delaunay (SciPy) pour réduire la complexité spatiale.
  2. Calcul des distances géodésiques réelles (Formule de Haversine).
  3. Algorithme de Kruskal (Union-Find / Ensemble Disjoint).
  4. Algorithme de Prim (File de priorité Min-Heap).
  5. Génération d'une carte interactive HTML avec Folium.
  6. Comparatif de performance et de complexité temporelle/spatiale.
"""

import json
import math
import os
import time
import heapq
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
import folium

# Chemins des fichiers
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(BASE_DIR, "data", "stations_velib.json")
OUTPUT_MAP = os.path.join(BASE_DIR, "carte_velib_optimisee.html")


def haversine_distance(lat1, lon1, lat2, lon2):
    """Calcule la distance en kilomètres entre deux coordonnées (Lat, Lon)."""
    R = 6371.0  # Rayon de la Terre en kilomètres
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon / 2) ** 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c


class DisjointSet:
    """Structure de données pour l'Algorithme de Kruskal (Union-Find)."""

    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, i):
        if self.parent[i] == i:
            return i
        self.parent[i] = self.find(self.parent[i])
        return self.parent[i]

    def union(self, i, j):
        root_i = self.find(i)
        root_j = self.find(j)
        if root_i != root_j:
            if self.rank[root_i] < self.rank[root_j]:
                root_i, root_j = root_j, root_i
            self.parent[root_j] = root_i
            if self.rank[root_i] == self.rank[root_j]:
                self.rank[root_i] += 1
            return True
        return False


def algo_kruskal(n_vertices, edges):
    """Algorithme de Kruskal pour trouver l'Arbre Couvrant Minimum (MST)."""
    start_time = time.perf_counter()
    sorted_edges = sorted(edges, key=lambda item: item[2])
    ds = DisjointSet(n_vertices)
    mst = []
    total_weight = 0.0

    for u, v, weight in sorted_edges:
        if ds.union(u, v):
            mst.append((u, v, weight))
            total_weight += weight
            if len(mst) == n_vertices - 1:
                break

    exec_time = (time.perf_counter() - start_time) * 1000
    return mst, total_weight, exec_time


def algo_prim(n_vertices, edges):
    """Algorithme de Prim pour trouver l'Arbre Couvrant Minimum (MST)."""
    start_time = time.perf_counter()
    adj = {i: [] for i in range(n_vertices)}
    for u, v, weight in edges:
        adj[u].append((weight, v))
        adj[v].append((weight, u))

    visited = [False] * n_vertices
    pq = [(0.0, 0, -1)]  # (weight, current_node, parent_node)
    mst = []
    total_weight = 0.0

    while pq and len(mst) < n_vertices:
        weight, u, parent = heapq.heappop(pq)
        if visited[u]:
            continue
        visited[u] = True
        if parent != -1:
            mst.append((parent, u, weight))
            total_weight += weight

        for next_weight, v in adj[u]:
            if not visited[v]:
                heapq.heappush(pq, (next_weight, v, u))

    exec_time = (time.perf_counter() - start_time) * 1000
    return mst, total_weight, exec_time


def main():
    print("=" * 65)
    print("  OPTIMISATION DU RÉSEAU VÉLIB PAR ARBRE COUVRANT MINIMUM (MST)")
    print("  Auteur : Misbaou DIALLO")
    print("=" * 65)

    # 1. Chargement des données
    with open(DATA_FILE, "r", encoding="utf-8") as f:
        stations = json.load(f)

    df = pd.DataFrame(stations)
    n = len(df)
    print(f"\n[+] {n} stations Vélib chargées en mémoire.")

    # 2. Triangulation de Delaunay
    coords = df[["longitude", "latitude"]].values
    tri = Delaunay(coords)

    # Extraire les arêtes uniques de la triangulation
    edges_set = set()
    for simplex in tri.simplices:
        for i in range(3):
            for j in range(i + 1, 3):
                u, v = simplex[i], simplex[j]
                if u > v:
                    u, v = v, u
                edges_set.add((u, v))

    print(f"[+] Triangulation de Delaunay générée ({len(edges_set)} connexions candidates).")

    # Calcul des longueurs réelles (en km)
    edges = []
    for u, v in edges_set:
        lat1, lon1 = df.loc[u, "latitude"], df.loc[u, "longitude"]
        lat2, lon2 = df.loc[v, "latitude"], df.loc[v, "longitude"]
        dist = haversine_distance(lat1, lon1, lat2, lon2)
        edges.append((u, v, dist))

    # 3. Exécution des Algorithmes MST
    mst_kruskal, weight_kruskal, time_kruskal = algo_kruskal(n, edges)
    mst_prim, weight_prim, time_prim = algo_prim(n, edges)

    print("\n" + "-" * 50)
    print("  RÉSULTATS DES ALGORITHMES DE RECHERCHE DE MST")
    print("-" * 50)
    print(f"• Kruskal  -> Longueur totale : {weight_kruskal:.3f} km | Temps : {time_kruskal:.4f} ms")
    print(f"• Prim     -> Longueur totale : {weight_prim:.3f} km | Temps : {time_prim:.4f} ms")
    print("-" * 50)

    # 4. Génération de la Carte Interactive Folium
    center_lat = df["latitude"].mean()
    center_lon = df["longitude"].mean()
    m = folium.Map(location=[center_lat, center_lon], zoom_start=13, tiles="OpenStreetMap")

    # Ajouter les stations (Marqueurs)
    for idx, row in df.iterrows():
        folium.CircleMarker(
            location=[row["latitude"], row["longitude"]],
            radius=6,
            popup=f"<b>{row['nom']}</b><br>Capacité: {row['capacite']} vélos<br>Code postal: {row['arrondissement']}",
            color="#2B6CB0",
            fill=True,
            fill_color="#3182CE",
            fill_opacity=0.8
        ).add_to(m)

    # Tracer les connexions optimales du MST (Lignes vertes)
    for u, v, weight in mst_kruskal:
        loc1 = [df.loc[u, "latitude"], df.loc[u, "longitude"]]
        loc2 = [df.loc[v, "latitude"], df.loc[v, "longitude"]]
        folium.PolyLine(
            locations=[loc1, loc2],
            weight=3.5,
            color="#38A169",
            opacity=0.85,
            tooltip=f"{weight*1000:.0f} mètres"
        ).add_to(m)

    m.save(OUTPUT_MAP)
    print(f"\n[✓] Carte interactive enregistrée dans : {OUTPUT_MAP}")
    print("[✓] Étude terminée avec succès !")


if __name__ == "__main__":
    main()
