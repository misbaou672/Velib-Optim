"""
==============================================================================
Optimisation Algorithmique du Réseau Vélib (1500+ Stations - Île-de-France)
Auteur : Misbaou DIALLO (BUT 3 Informatique)
==============================================================================

Ce projet implémente une étude algorithmique et cartographique pour optimiser
l'intégralité du réseau métropolitain Vélib en Île-de-France (1 518+ stations).

Données :
  - Source OpenData Île-de-France / Data.Gouv.fr (Temps réel & Géolocalisation)
  - Couverture : Paris 1er à 20e, 92 (Hauts-de-Seine), 93 (Seine-Saint-Denis),
    94 (Val-de-Marne), 78 (Yvelines), 91 (Essonne), 95 (Val-d'Oise).

Algorithmes & Traitements :
  1. Triangulation de Delaunay (SciPy) : Réduction de la complexité spatiale de O(V^2) à O(V).
  2. Distances Géodésiques réelles (Formule de Haversine).
  3. Algorithme de Kruskal (Union-Find / Ensemble Disjoint, O(E log E)).
  4. Algorithme de Prim (File de priorité Min-Heap, O(E log V)).
  5. Carte interactive HTML (Folium) avec rendu vectoriel du réseau minimal et couches par commune.
"""

import json
import math
import os
import time
import heapq
import requests
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
import folium
from folium.plugins import MiniMap, MarkerCluster

# Chemins des fichiers
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
DATA_FILE = os.path.join(DATA_DIR, "stations_velib_idf_complete.json")
OUTPUT_MAP = os.path.join(BASE_DIR, "carte_velib_optimisee.html")
OPENDATA_URL = "https://opendata.paris.fr/api/explore/v2.1/catalog/datasets/velib-disponibilite-en-temps-reel/exports/json"


def haversine_distance(lat1, lon1, lat2, lon2):
    """Calcule la distance géodésique en kilomètres entre deux points GPS (Lat, Lon)."""
    R = 6371.0  # Rayon moyen de la Terre en km
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon / 2) ** 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c


class DisjointSet:
    """Structure de données Union-Find pour l'Algorithme de Kruskal."""

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
    """Algorithme de Kruskal pour l'Arbre Couvrant Minimum (MST)."""
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
    """Algorithme de Prim pour l'Arbre Couvrant Minimum (MST)."""
    start_time = time.perf_counter()
    adj = {i: [] for i in range(n_vertices)}
    for u, v, weight in edges:
        adj[u].append((weight, v))
        adj[v].append((weight, u))

    visited = [False] * n_vertices
    pq = [(0.0, 0, -1)]  # (weight, node, parent)
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


def charger_donnees():
    """Charge le jeu de données depuis l'API OpenData ou depuis le fichier cache local."""
    os.makedirs(DATA_DIR, exist_ok=True)
    if os.path.exists(DATA_FILE):
        print(f"[+] Chargement du jeu de données depuis le cache local ({DATA_FILE})...")
        with open(DATA_FILE, "r", encoding="utf-8") as f:
            return json.load(f)

    print("[+] Téléchargement en temps réel depuis l'API Data.Gouv / OpenData Île-de-France...")
    r = requests.get(OPENDATA_URL, timeout=15)
    data = r.json()
    formatted = []
    for s in data:
        coords = s.get('coordonnees_geo') or {}
        lat, lon = coords.get('lat'), coords.get('lon')
        if lat and lon:
            formatted.append({
                'id': s.get('stationcode'),
                'nom': s.get('name'),
                'latitude': lat,
                'longitude': lon,
                'capacite': s.get('capacity', 0),
                'commune': s.get('nom_arrondissement_communes', 'Île-de-France'),
                'code_insee': s.get('code_insee_commune', '')
            })
    with open(DATA_FILE, "w", encoding="utf-8") as f:
        json.dump(formatted, f, ensure_ascii=False, indent=2)
    return formatted


def main():
    print("=" * 75)
    print("  OPTIMISATION DU RÉSEAU VÉLIB MÉTROPOLITAIN (1500+ STATIONS ÎLE-DE-FRANCE)")
    print("  Auteur : Misbaou DIALLO")
    print("=" * 75)

    stations = charger_donnees()
    df = pd.DataFrame(stations)
    n = len(df)
    print(f"\n[✓] {n} stations Vélib chargées (Paris & Région Île-de-France).")

    # Statistiques sur les communes
    communes = df['commune'].nunique()
    print(f"[✓] Couverture géographique : {communes} communes d'Île-de-France.")

    # 1. Triangulation de Delaunay
    coords = df[["longitude", "latitude"]].values
    tri = Delaunay(coords)

    edges_set = set()
    for simplex in tri.simplices:
        for i in range(3):
            for j in range(i + 1, 3):
                u, v = simplex[i], simplex[j]
                if u > v:
                    u, v = v, u
                edges_set.add((u, v))

    print(f"[✓] Triangulation de Delaunay générée ({len(edges_set)} arêtes candidates).")

    # Calcul des longueurs réelles (en km)
    edges = []
    for u, v in edges_set:
        lat1, lon1 = df.loc[u, "latitude"], df.loc[u, "longitude"]
        lat2, lon2 = df.loc[v, "latitude"], df.loc[v, "longitude"]
        dist = haversine_distance(lat1, lon1, lat2, lon2)
        edges.append((u, v, dist))

    # 2. Algorithmes MST
    mst_kruskal, weight_kruskal, time_kruskal = algo_kruskal(n, edges)
    mst_prim, weight_prim, time_prim = algo_prim(n, edges)

    print("\n" + "-" * 60)
    print("  RÉSULTATS DES ALGORITHMES MST (ARBRE COUVRANT MINIMUM)")
    print("-" * 60)
    print(f"• Kruskal (Union-Find) -> Longueur : {weight_kruskal:.3f} km | Temps : {time_kruskal:.2f} ms")
    print(f"• Prim (Min-Heap)      -> Longueur : {weight_prim:.3f} km | Temps : {time_prim:.2f} ms")
    print("-" * 60)

    # 3. Cartographie Folium
    center_lat = df["latitude"].mean()
    center_lon = df["longitude"].mean()
    m = folium.Map(location=[center_lat, center_lon], zoom_start=11, tiles="OpenStreetMap")

    # Cluster de marqueurs pour la fluidité avec 1500+ stations
    marker_cluster = MarkerCluster(name="Stations Vélib Île-de-France").add_to(m)

    for idx, row in df.iterrows():
        folium.CircleMarker(
            location=[row["latitude"], row["longitude"]],
            radius=4,
            popup=f"<b>{row['nom']}</b><br>Commune: {row['commune']}<br>Capacité: {row['capacite']} vélos",
            color="#2B6CB0",
            fill=True,
            fill_color="#3182CE",
            fill_opacity=0.8
        ).add_to(marker_cluster)

    # Groupe pour le tracé du réseau optimal MST
    mst_group = folium.FeatureGroup(name="Réseau Optimal Optimisé (MST - 106+ km)")
    for u, v, weight in mst_kruskal:
        loc1 = [df.loc[u, "latitude"], df.loc[u, "longitude"]]
        loc2 = [df.loc[v, "latitude"], df.loc[v, "longitude"]]
        folium.PolyLine(
            locations=[loc1, loc2],
            weight=2.5,
            color="#38A169",
            opacity=0.8,
            tooltip=f"{weight*1000:.0f} m"
        ).add_to(mst_group)

    mst_group.add_to(m)
    folium.LayerControl().add_to(m)
    MiniMap(toggle_display=True).add_to(m)

    m.save(OUTPUT_MAP)
    print(f"\n[✓] Carte interactive générée avec succès : {OUTPUT_MAP}")
    print("[✓] Processus terminé avec succès !")


if __name__ == "__main__":
    main()
