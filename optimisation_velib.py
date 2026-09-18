"""
Optimisation et analyse spatiale du réseau Vélib en Île-de-France (1 518 stations).

Projet de théorie des graphes et d'optimisation :
- Triangulation de Delaunay avec coloration selon la surface (densité spatiale)
- Algorithmes d'Arbre Couvrant Minimum (Kruskal & Prim)
- Recherche de plus court chemin (Dijkstra)
- Visualisation interactive avec Folium, Leaflet et Chart.js

Auteur : Misbaou DIALLO (BUT 3 Informatique)
"""

import argparse
import json
import math
import os
import sys
import time
import heapq
import requests
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import folium
from folium.plugins import MiniMap, MarkerCluster

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
DATA_FILE = os.path.join(DATA_DIR, "stations_velib_idf_complete.json")
OUTPUT_MAP = os.path.join(BASE_DIR, "carte_velib_optimisee.html")
REPORT_FILE = os.path.join(BASE_DIR, "rapport_statistiques_velib.json")
GRAPH_IMAGE = os.path.join(DATA_DIR, "graphiques_velib.png")
OPENDATA_URL = "https://opendata.paris.fr/api/explore/v2.1/catalog/datasets/velib-disponibilite-en-temps-reel/exports/json"


def haversine_distance(lat1, lon1, lat2, lon2):
    """Calcule la distance en kilomètres entre deux points GPS."""
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon / 2) ** 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c


def calculer_surface_triangle(p1, p2, p3):
    """Calcule la superficie approximative en km² d'un triangle (lon, lat)."""
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    deg_area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    # 1 deg² ≈ 8103 km² à la latitude de Paris (~48.85°)
    return deg_area * 8103.0


def obtenir_couleur_delaunay(area, min_area, max_area):
    """
    Retourne la couleur et l'opacité selon la superficie du triangle (échelle logarithmique) :
    - Petite surface (densité forte, centre-ville) -> Couleur sombre & opaque
    - Grande surface (densité faible, périphérie) -> Couleur claire & translucide
    """
    safe_area = max(area, 1e-7)
    safe_min = max(min_area, 1e-7)
    safe_max = max(max_area, 1e-7)

    log_area = math.log10(safe_area)
    log_min = math.log10(safe_min)
    log_max = math.log10(safe_max)

    if log_max == log_min:
        norm = 0.5
    else:
        norm = (log_area - log_min) / (log_max - log_min)
        norm = max(0.0, min(1.0, norm))

    if norm < 0.15:
        return "#0F172A", 0.80
    elif norm < 0.30:
        return "#1E1B4B", 0.70
    elif norm < 0.45:
        return "#312E81", 0.60
    elif norm < 0.60:
        return "#4338CA", 0.50
    elif norm < 0.75:
        return "#6D28D9", 0.40
    elif norm < 0.88:
        return "#A855F7", 0.28
    else:
        return "#E0E7FF", 0.18


class DisjointSet:
    """Structure Union-Find pour l'algorithme de Kruskal."""

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
    pq = [(0.0, 0, -1)]
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


def algo_dijkstra(n_vertices, edges, start_node, target_node):
    """Algorithme de Dijkstra pour trouver le plus court chemin."""
    start_time = time.perf_counter()
    adj = {i: [] for i in range(n_vertices)}
    for u, v, weight in edges:
        adj[u].append((weight, v))
        adj[v].append((weight, u))

    distances = {i: float('inf') for i in range(n_vertices)}
    distances[start_node] = 0.0
    predecessors = {i: None for i in range(n_vertices)}
    pq = [(0.0, start_node)]

    while pq:
        curr_dist, u = heapq.heappop(pq)
        if u == target_node:
            break
        if curr_dist > distances[u]:
            continue

        for weight, v in adj[u]:
            distance = curr_dist + weight
            if distance < distances[v]:
                distances[v] = distance
                predecessors[v] = u
                heapq.heappush(pq, (distance, v))

    path = []
    curr = target_node
    while curr is not None:
        path.append(curr)
        curr = predecessors[curr]
    path.reverse()

    exec_time = (time.perf_counter() - start_time) * 1000
    return path, distances[target_node], exec_time


def charger_donnees():
    """Charge les données réelles des stations depuis le cache local ou l'API OpenData."""
    os.makedirs(DATA_DIR, exist_ok=True)
    if os.path.exists(DATA_FILE):
        with open(DATA_FILE, "r", encoding="utf-8") as f:
            return json.load(f)

    print("Téléchargement des stations depuis l'API OpenData Paris...")
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


def generer_graphiques_matplotlib(df, edges):
    """Génère le tableau de bord analytique en PNG avec Matplotlib."""
    plt.style.use('dark_background')
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Analyse du Réseau Vélib Île-de-France (1 518 stations)",
                 fontsize=14, fontweight='bold', color='#6366F1')

    # Top 10 par capacité
    top10 = df.sort_values(by='capacite', ascending=False).head(10)
    axes[0, 0].barh(top10['nom'].str[:25], top10['capacite'], color='#38BDF8')
    axes[0, 0].set_title("Top 10 des stations par capacité", fontsize=11, fontweight='bold')
    axes[0, 0].set_xlabel("Capacité (vélos)")
    axes[0, 0].invert_yaxis()

    # Répartition par commune
    communes = df['commune'].value_counts().head(8)
    axes[0, 1].pie(communes.values, labels=communes.index, autopct='%1.1f%%',
                   colors=['#818CF8', '#34D399', '#FBBF24', '#F87171', '#A78BFA', '#F472B6', '#38BDF8', '#4ADE80'])
    axes[0, 1].set_title("Répartition des stations par commune", fontsize=11, fontweight='bold')

    # Distribution des distances
    distances_m = [e[2] * 1000 for e in edges]
    axes[1, 0].hist(distances_m, bins=30, color='#34D399', edgecolor='#111827')
    axes[1, 0].set_title("Distances inter-stations (Delaunay)", fontsize=11, fontweight='bold')
    axes[1, 0].set_xlabel("Distance (mètres)")
    axes[1, 0].set_ylabel("Nombre d'arêtes")

    # Capacités
    axes[1, 1].hist(df['capacite'], bins=20, color='#FBBF24', edgecolor='#111827')
    axes[1, 1].set_title("Distribution des capacités des stations", fontsize=11, fontweight='bold')
    axes[1, 1].set_xlabel("Nombre de bornettes")
    axes[1, 1].set_ylabel("Nombre de stations")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(GRAPH_IMAGE, dpi=200)
    plt.close()


def generer_statistiques(df, edges, weight_mst):
    """Exporte les métriques du réseau au format JSON."""
    total_stations = len(df)
    total_capacite = int(df['capacite'].sum())
    moyenne_capacite = float(df['capacite'].mean())

    top_capacites = df.sort_values(by='capacite', ascending=False).head(10)[
        ['nom', 'commune', 'capacite']
    ].to_dict(orient='records')

    par_commune = df['commune'].value_counts().head(15).to_dict()
    distances = [e[2] for e in edges]

    rapport = {
        "metriques_generales": {
            "total_stations": total_stations,
            "total_communes": df['commune'].nunique(),
            "capacite_totale_velos": total_capacite,
            "capacite_moyenne_station": round(moyenne_capacite, 2),
            "distance_mst_totale_km": round(weight_mst, 2),
            "nombre_connexions_delaunay": len(edges),
            "distance_inter_station_moyenne_km": round(sum(distances) / len(distances), 3)
        },
        "top_10_stations_capacite": top_capacites,
        "repartition_par_commune": par_commune
    }

    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        json.dump(rapport, f, ensure_ascii=False, indent=2)

    return rapport


def generer_carte_html_interactive(df, tri, mst_edges, edges, default_start_idx=0, default_target_idx=15):
    """Génère la carte web interactive Folium / Leaflet."""
    center_lat = df["latitude"].mean()
    center_lon = df["longitude"].mean()
    m = folium.Map(location=[center_lat, center_lon], zoom_start=11, tiles="OpenStreetMap")

    marker_cluster = MarkerCluster(name="Stations Vélib (1 518)").add_to(m)

    stations_js_data = []
    coords_list = []
    for idx, row in df.iterrows():
        coords_list.append((row["longitude"], row["latitude"]))
        stations_js_data.append({
            "idx": idx,
            "id": row["id"],
            "nom": row["nom"],
            "lat": row["latitude"],
            "lon": row["longitude"],
            "capacite": row["capacite"],
            "commune": row["commune"]
        })

        popup_html = f"""
        <div style="font-family: system-ui, -apple-system, sans-serif; min-width:180px;">
            <div style="font-weight:700; color:#1E40AF; font-size:14px; margin-bottom:4px;">{row['nom']}</div>
            <div style="color:#475569; font-size:12px; display:flex; align-items:center; gap:5px; margin-bottom:2px;">
                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M20 10c0 6-8 12-8 12s-8-6-8-12a8 8 0 0 1 16 0Z"/><circle cx="12" cy="10" r="3"/></svg>
                <span>Commune : <b>{row['commune']}</b></span>
            </div>
            <div style="color:#475569; font-size:12px; display:flex; align-items:center; gap:5px; margin-bottom:8px;">
                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="5.5" cy="17.5" r="3.5"/><circle cx="18.5" cy="17.5" r="3.5"/><path d="M15 6a1 1 0 1 0 0-2 1 1 0 0 0 0 2zm-3 11.5L9.5 10l-3 3.5M12 17.5V10l3.5-4H18"/></svg>
                <span>Capacité : <b>{row['capacite']}</b> vélos</span>
            </div>
            <button onclick="selectStationByClick({idx})" style="
                width:100%; background:#2563EB; color:white; border:none; padding:6px 8px; border-radius:4px; font-weight:600; font-size:12px; cursor:pointer; display:flex; align-items:center; justify-content:center; gap:5px;
            ">
                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="10"/><line x1="12" y1="8" x2="12" y2="16"/><line x1="8" y1="12" x2="16" y2="12"/></svg>
                <span>Sélectionner pour itinéraire</span>
            </button>
        </div>
        """
        folium.CircleMarker(
            location=[row["latitude"], row["longitude"]],
            radius=5,
            popup=folium.Popup(popup_html, max_width=250),
            color="#2B6CB0",
            fill=True,
            fill_color="#3182CE",
            fill_opacity=0.8
        ).add_to(marker_cluster)

    triangle_areas = []
    triangle_data = []
    for simplex in tri.simplices:
        p1 = coords_list[simplex[0]]
        p2 = coords_list[simplex[1]]
        p3 = coords_list[simplex[2]]
        area = calculer_surface_triangle(p1, p2, p3)
        triangle_areas.append(area)
        triangle_data.append((simplex, area))

    min_area = min(triangle_areas)
    max_area = max(triangle_areas)

    delaunay_group = folium.FeatureGroup(name="Maillage Delaunay (Coloration par superficie)")
    for simplex, area in triangle_data:
        p1 = [df.loc[simplex[0], "latitude"], df.loc[simplex[0], "longitude"]]
        p2 = [df.loc[simplex[1], "latitude"], df.loc[simplex[1], "longitude"]]
        p3 = [df.loc[simplex[2], "latitude"], df.loc[simplex[2], "longitude"]]

        color, opacity = obtenir_couleur_delaunay(area, min_area, max_area)
        area_str = f"{area * 100:.1f} ha" if area < 1.0 else f"{area:.2f} km²"

        folium.Polygon(
            locations=[p1, p2, p3],
            color="#4C1D95",
            weight=1.0,
            fill=True,
            fill_color=color,
            fill_opacity=opacity,
            tooltip=f"Triangle Delaunay<br>Surface : <b>{area_str}</b><br>Densité : {'Forte (couleur sombre)' if opacity > 0.5 else 'Faible (couleur claire)'}"
        ).add_to(delaunay_group)

    mst_group = folium.FeatureGroup(name="Réseau Optimal (MST - 502 km)")
    for u, v, weight in mst_edges:
        loc1 = [df.loc[u, "latitude"], df.loc[u, "longitude"]]
        loc2 = [df.loc[v, "latitude"], df.loc[v, "longitude"]]
        folium.PolyLine(
            locations=[loc1, loc2],
            weight=3.2,
            color="#10B981",
            opacity=0.95,
            tooltip=f"MST: {weight*1000:.0f} m"
        ).add_to(mst_group)

    delaunay_group.add_to(m)
    mst_group.add_to(m)

    folium.LayerControl(position='topleft', collapsed=False).add_to(m)
    MiniMap(toggle_display=True, position='bottomleft').add_to(m)

    m.save(OUTPUT_MAP)

    with open(OUTPUT_MAP, "r", encoding="utf-8") as f:
        html_content = f.read()

    edges_js = [{"u": int(u), "v": int(v), "w": round(float(w), 4)} for u, v, w in edges]

    top10_df = df.sort_values(by='capacite', ascending=False).head(8)
    chart_top_labels = top10_df['nom'].str[:20].tolist()
    chart_top_values = top10_df['capacite'].tolist()

    communes_top = df['commune'].value_counts().head(6)
    chart_commune_labels = communes_top.index.tolist()
    chart_commune_values = communes_top.values.tolist()

    dashboard_ui_html = f"""
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>

    <style>
    .leaflet-top.leaflet-left {{
        top: 80px !important;
        left: 15px !important;
    }}
    .kpi-card {{
        background: rgba(15, 23, 42, 0.92);
        color: white;
        padding: 8px 14px;
        border-radius: 8px;
        backdrop-filter: blur(8px);
        border: 1px solid #334155;
        display: flex;
        align-items: center;
        gap: 10px;
    }}
    .kpi-icon {{
        display: flex;
        align-items: center;
        justify-content: center;
        width: 32px;
        height: 32px;
        border-radius: 6px;
        background: rgba(255, 255, 255, 0.1);
    }}
    .ui-flex {{
        display: flex;
        align-items: center;
        gap: 6px;
    }}
    </style>

    <!-- Bandeau KPI avec icônes SVG -->
    <div id="kpi-banner" style="
        position: fixed; top: 15px; left: 50%; transform: translateX(-50%); display: flex; gap: 12px; z-index: 9999; font-family: system-ui, -apple-system, sans-serif;
    ">
        <div class="kpi-card">
            <div class="kpi-icon" style="color: #38BDF8;">
                <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M20 10c0 6-8 12-8 12s-8-6-8-12a8 8 0 0 1 16 0Z"/><circle cx="12" cy="10" r="3"/></svg>
            </div>
            <div>
                <div style="font-size: 10px; color: #94A3B8; text-transform: uppercase; font-weight: bold;">Stations</div>
                <div style="font-size: 17px; font-weight: bold; color: #38BDF8;">1 518</div>
            </div>
        </div>
        <div class="kpi-card">
            <div class="kpi-icon" style="color: #34D399;">
                <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="5.5" cy="17.5" r="3.5"/><circle cx="18.5" cy="17.5" r="3.5"/><path d="M15 6a1 1 0 1 0 0-2 1 1 0 0 0 0 2zm-3 11.5L9.5 10l-3 3.5M12 17.5V10l3.5-4H18"/></svg>
            </div>
            <div>
                <div style="font-size: 10px; color: #94A3B8; text-transform: uppercase; font-weight: bold;">Bornettes & Vélos</div>
                <div style="font-size: 17px; font-weight: bold; color: #34D399;">49 060</div>
            </div>
        </div>
        <div class="kpi-card">
            <div class="kpi-icon" style="color: #FBBF24;">
                <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="6" cy="19" r="3"/><circle cx="18" cy="5" r="3"/><circle cx="12" cy="6" r="3"/><path d="M8.5 17l2.5-8.5M15.5 17l-2.5-8.5"/></svg>
            </div>
            <div>
                <div style="font-size: 10px; color: #94A3B8; text-transform: uppercase; font-weight: bold;">Réseau MST</div>
                <div style="font-size: 17px; font-weight: bold; color: #FBBF24;">502.09 km</div>
            </div>
        </div>
        <div class="kpi-card">
            <div class="kpi-icon" style="color: #F472B6;">
                <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><polygon points="1 6 1 22 8 18 16 22 23 18 23 2 16 6 8 2 1 6"/><line x1="8" y1="2" x2="8" y2="18"/><line x1="16" y1="6" x2="16" y2="22"/></svg>
            </div>
            <div>
                <div style="font-size: 10px; color: #94A3B8; text-transform: uppercase; font-weight: bold;">Communes IDF</div>
                <div style="font-size: 17px; font-weight: bold; color: #F472B6;">69</div>
            </div>
        </div>
    </div>

    <!-- Panneau de contrôle latéral -->
    <div id="route-panel" style="
        position: fixed; top: 80px; right: 15px; width: 340px; max-height: calc(100vh - 100px); overflow-y: auto; background: rgba(255, 255, 255, 0.96); border-radius: 12px; padding: 16px; box-shadow: 0 10px 25px rgba(0,0,0,0.15); z-index: 9999; font-family: system-ui, -apple-system, sans-serif; font-size: 13px; backdrop-filter: blur(8px); border: 1px solid #E2E8F0;
    ">
        <h3 style="margin:0 0 8px 0; color:#1E293B; font-size:15px; display:flex; align-items:center; gap:6px;">
            <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#2563EB" stroke-width="2"><circle cx="12" cy="12" r="10"/><path d="m16.24 7.76-2.12 6.36-6.36 2.12 2.12-6.36z"/></svg>
            <span>Calculateur d'Itinéraire (Dijkstra)</span>
        </h3>
        
        <p style="margin:0 0 10px 0; color:#64748B; font-size:11px;">
            Sélectionnez deux stations dans les listes ci-dessous ou directement en cliquant sur la carte.
        </p>
        
        <label style="font-weight:600; color:#475569;">Station de départ :</label>
        <select id="start-station" style="width:100%; padding:7px; margin:4px 0 10px 0; border-radius:6px; border:1px solid #CBD5E0;"></select>
        
        <label style="font-weight:600; color:#475569;">Station d'arrivée :</label>
        <select id="target-station" style="width:100%; padding:7px; margin:4px 0 12px 0; border-radius:6px; border:1px solid #CBD5E0;"></select>
        
        <button onclick="calculateRoute()" style="
            width:100%; background:#2563EB; color:white; border:none; padding:9px; border-radius:6px; font-weight:bold; cursor:pointer; margin-bottom:8px; display:flex; align-items:center; justify-content:center; gap:6px;
        ">
            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/></svg>
            <span>Calculer le trajet le plus court</span>
        </button>
        
        <button id="toggle-delaunay-btn" onclick="toggleDelaunayLayer()" style="
            width:100%; background:#6B21A8; color:white; border:none; padding:8px; border-radius:6px; font-weight:600; cursor:pointer; display:flex; align-items:center; justify-content:center; gap:6px; transition: background 0.2s;
        ">
            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M12 2 2 22h20L12 2z"/></svg>
            <span>Afficher / Masquer Triangulation Delaunay</span>
        </button>

        <div id="route-results" style="margin-top:12px; display:none; padding:12px; background:#F8FAFC; border-radius:8px; border:1px solid #E2E8F0;">
            <div style="font-weight:bold; color:#1D4ED8; margin-bottom:6px; display:flex; align-items:center; gap:6px;">
                <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M22 11.08V12a10 10 0 1 1-5.93-9.14"/><polyline points="22 4 12 14.01 9 11.01"/></svg>
                <span>Résultats du trajet :</span>
            </div>
            <div style="display:flex; align-items:center; gap:6px; margin-bottom:3px;">
                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="#64748B" stroke-width="2"><path d="M21.3 15.3a2.4 2.4 0 0 1 0 3.4l-2.6 2.6a2.4 2.4 0 0 1-3.4 0L2.7 8.7a2.4 2.4 0 0 1 0-3.4l2.6-2.6a2.4 2.4 0 0 1 3.4 0l12.6 12.6z"/></svg>
                <span>Distance : <b id="route-dist" style="color:#0F172A;">-</b></span>
            </div>
            <div style="display:flex; align-items:center; gap:6px; margin-bottom:3px;">
                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="#64748B" stroke-width="2"><circle cx="12" cy="12" r="10"/><polyline points="12 6 12 12 16 14"/></svg>
                <span>Temps estimé (15 km/h) : <b id="route-time" style="color:#0F172A;">-</b></span>
            </div>
            <div style="display:flex; align-items:center; gap:6px;">
                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="#64748B" stroke-width="2"><path d="M20 10c0 6-8 12-8 12s-8-6-8-12a8 8 0 0 1 16 0Z"/><circle cx="12" cy="10" r="3"/></svg>
                <span>Stations traversées : <b id="route-hops" style="color:#0F172A;">-</b></span>
            </div>
        </div>

        <hr style="margin: 16px 0; border: 0; border-top: 1px solid #E2E8F0;">

        <h3 style="margin:0 0 12px 0; color:#1E293B; font-size:14px; display:flex; align-items:center; gap:6px;">
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="#1E293B" stroke-width="2"><line x1="18" y1="20" x2="18" y2="10"/><line x1="12" y1="20" x2="12" y2="4"/><line x1="6" y1="20" x2="6" y2="14"/></svg>
            <span>Statistiques du réseau</span>
        </h3>

        <div style="margin-bottom: 16px;">
            <div style="font-size:11px; font-weight:bold; color:#64748B; margin-bottom:6px;">Top stations par capacité</div>
            <canvas id="chartTopCapacity" height="160"></canvas>
        </div>

        <div>
            <div style="font-size:11px; font-weight:bold; color:#64748B; margin-bottom:6px;">Répartition par commune</div>
            <canvas id="chartCommunes" height="160"></canvas>
        </div>
    </div>

    <script>
    const STATIONS = {json.dumps(stations_js_data)};
    const EDGES = {json.dumps(edges_js)};
    let activeRouteLayer = null;
    let clickSelectionStep = 0;
    let delaunayLayerRef = null;

    document.addEventListener("DOMContentLoaded", function() {{
        const selectStart = document.getElementById("start-station");
        const selectTarget = document.getElementById("target-station");
        
        STATIONS.forEach(s => {{
            let opt1 = document.createElement("option");
            opt1.value = s.idx;
            opt1.textContent = s.nom + " (" + s.commune + ")";
            selectStart.appendChild(opt1);
            
            let opt2 = document.createElement("option");
            opt2.value = s.idx;
            opt2.textContent = s.nom + " (" + s.commune + ")";
            selectTarget.appendChild(opt2);
        }});
        
        selectStart.selectedIndex = {default_start_idx};
        selectTarget.selectedIndex = {default_target_idx};

        new Chart(document.getElementById('chartTopCapacity'), {{
            type: 'bar',
            data: {{
                labels: {json.dumps(chart_top_labels)},
                datasets: [{{
                    label: 'Capacité',
                    data: {json.dumps(chart_top_values)},
                    backgroundColor: '#3B82F6',
                    borderRadius: 4
                }}]
            }},
            options: {{
                responsive: true,
                plugins: {{ legend: {{ display: false }} }},
                scales: {{ y: {{ beginAtZero: true }} }}
            }}
        }});

        new Chart(document.getElementById('chartCommunes'), {{
            type: 'doughnut',
            data: {{
                labels: {json.dumps(chart_commune_labels)},
                datasets: [{{
                    data: {json.dumps(chart_commune_values)},
                    backgroundColor: ['#3B82F6', '#10B981', '#F59E0B', '#EF4444', '#8B5CF6', '#EC4899']
                }}]
            }},
            options: {{
                responsive: true,
                plugins: {{ legend: {{ position: 'right', labels: {{ boxWidth: 10, font: {{ size: 10 }} }} }} }}
            }}
        }});
    }});

    function toggleDelaunayLayer() {{
        const mapObj = Object.values(window).find(v => v && v.addLayer && v.eachLayer);
        if (!mapObj) return;

        if (!delaunayLayerRef) {{
            mapObj.eachLayer(layer => {{
                if (layer.options && (
                    (layer.options.name && layer.options.name.toLowerCase().includes("delaunay")) ||
                    (layer.options.overlayName && layer.options.overlayName.toLowerCase().includes("delaunay"))
                )) {{
                    delaunayLayerRef = layer;
                }}
            }});
        }}

        if (delaunayLayerRef) {{
            const btn = document.getElementById("toggle-delaunay-btn");
            if (mapObj.hasLayer(delaunayLayerRef)) {{
                mapObj.removeLayer(delaunayLayerRef);
                btn.style.background = "#64748B";
            }} else {{
                mapObj.addLayer(delaunayLayerRef);
                btn.style.background = "#6B21A8";
            }}
        }}
    }}

    function selectStationByClick(stationIdx) {{
        const selectStart = document.getElementById("start-station");
        const selectTarget = document.getElementById("target-station");

        if (clickSelectionStep === 0 || clickSelectionStep === 2) {{
            selectStart.value = stationIdx;
            clickSelectionStep = 1;
            console.log("Station de départ : " + STATIONS[stationIdx].nom);
        }} else if (clickSelectionStep === 1) {{
            selectTarget.value = stationIdx;
            clickSelectionStep = 2;
            calculateRoute();
        }}
    }}

    function calculateRoute() {{
        const uStart = parseInt(document.getElementById("start-station").value);
        const uTarget = parseInt(document.getElementById("target-station").value);
        
        if (uStart === uTarget) {{
            alert("Veuillez choisir deux stations différentes.");
            return;
        }}
        
        const n = STATIONS.length;
        const adj = Array.from({{ length: n }}, () => []);
        EDGES.forEach(e => {{
            adj[e.u].push({{ node: e.v, w: e.w }});
            adj[e.v].push({{ node: e.u, w: e.w }});
        }});
        
        const dist = Array(n).fill(Infinity);
        const parent = Array(n).fill(null);
        dist[uStart] = 0;
        const visited = Array(n).fill(false);
        
        for (let i = 0; i < n; i++) {{
            let u = -1;
            for (let j = 0; j < n; j++) {{
                if (!visited[j] && (u === -1 || dist[j] < dist[u])) u = j;
            }}
            if (dist[u] === Infinity || u === uTarget) break;
            visited[u] = true;
            
            adj[u].forEach(edge => {{
                if (dist[u] + edge.w < dist[edge.node]) {{
                    dist[edge.node] = dist[u] + edge.w;
                    parent[edge.node] = u;
                }}
            }});
        }}
        
        const path = [];
        let curr = uTarget;
        while (curr !== null) {{
            path.push(curr);
            curr = parent[curr];
        }}
        path.reverse();
        
        const totalDistKm = dist[uTarget];
        const minutes = Math.round((totalDistKm / 15) * 60);
        
        document.getElementById("route-results").style.display = "block";
        document.getElementById("route-dist").textContent = totalDistKm.toFixed(2) + " km";
        document.getElementById("route-time").textContent = minutes + " min";
        document.getElementById("route-hops").textContent = path.length + " stations";
        
        const mapObj = Object.values(window).find(v => v && v.addLayer && v.on);
        if (mapObj) {{
            if (activeRouteLayer) mapObj.removeLayer(activeRouteLayer);
            
            const routeCoords = path.map(idx => [STATIONS[idx].lat, STATIONS[idx].lon]);
            activeRouteLayer = L.polyline(routeCoords, {{
                color: '#EF4444',
                weight: 5,
                opacity: 0.9,
                dashArray: '6, 6'
            }}).addTo(mapObj);
            
            mapObj.fitBounds(activeRouteLayer.getBounds(), {{ padding: [50, 50] }});
        }}
    }}
    </script>
    """

    html_content = html_content.replace("</body>", f"{dashboard_ui_html}</body>")
    with open(OUTPUT_MAP, "w", encoding="utf-8") as f:
        f.write(html_content)


def main():
    parser = argparse.ArgumentParser(description="Optimisation et analyse du réseau Vélib Île-de-France")
    parser.add_argument("--depart", type=int, default=None, help="Index de la station de départ")
    parser.add_argument("--arrivee", type=int, default=None, help="Index de la station d'arrivée")
    args = parser.parse_args()

    print("Chargement des données Vélib...")
    stations = charger_donnees()
    df = pd.DataFrame(stations)
    n = len(df)

    print(f"-> {n} stations chargées sur {df['commune'].nunique()} communes d'Île-de-France.")

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

    edges = []
    for u, v in edges_set:
        lat1, lon1 = df.loc[u, "latitude"], df.loc[u, "longitude"]
        lat2, lon2 = df.loc[v, "latitude"], df.loc[v, "longitude"]
        dist = haversine_distance(lat1, lon1, lat2, lon2)
        edges.append((u, v, dist))

    print(f"-> Graphe Delaunay généré : {len(edges)} arêtes candidates.")

    # 2. Algorithmes MST
    mst_kruskal, weight_kruskal, time_kruskal = algo_kruskal(n, edges)
    mst_prim, weight_prim, time_prim = algo_prim(n, edges)

    print("\n--- Performances MST ---")
    print(f"Kruskal (Union-Find) : {weight_kruskal:.2f} km (exécuté en {time_kruskal:.2f} ms)")
    print(f"Prim (Min-Heap)      : {weight_prim:.2f} km (exécuté en {time_prim:.2f} ms)")

    # 3. Visualisation Matplotlib
    generer_graphiques_matplotlib(df, edges)

    # 4. Statistiques JSON
    generer_statistiques(df, edges, weight_kruskal)

    # 5. Carte HTML Folium
    generer_carte_html_interactive(df, tri, mst_kruskal, edges, args.depart or 0, args.arrivee or 15)
    print(f"\nCarte interactive générée avec succès : {OUTPUT_MAP}")


if __name__ == "__main__":
    main()
