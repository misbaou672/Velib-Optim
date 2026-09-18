"""
Optimisation et analyse spatiale du réseau Vélib en Île-de-France (1 518 stations).

Projet de théorie des graphes et d'optimisation :
- Connexion en direct à l'API OpenData Paris (Mise à jour en temps réel)
- Triangulation de Delaunay avec coloration selon la surface (densité spatiale)
- Algorithmes d'Arbre Couvrant Minimum (Kruskal & Prim)
- Recherche de plus court chemin (Dijkstra avec station étape / Waypoint)
- Cockpit Data & Analytics Spatiales complet (HeatMap, Filtrage Dynamique, Comparateur, Santé Réseau)

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
from folium.plugins import MiniMap, MarkerCluster, HeatMap

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
DATA_FILE = os.path.join(DATA_DIR, "stations_velib_idf_complete.json")
OUTPUT_MAP = os.path.join(BASE_DIR, "carte_velib_optimisee.html")
REPORT_FILE = os.path.join(BASE_DIR, "rapport_statistiques_velib.json")
GRAPH_IMAGE = os.path.join(DATA_DIR, "graphiques_velib.png")
OPENDATA_URL = "https://opendata.paris.fr/api/explore/v2.1/catalog/datasets/velib-disponibilite-en-temps-reel/exports/json?limit=-1"


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


def charger_donnees():
    """Charge les données réelles et temps réel des stations depuis l'API OpenData Paris."""
    os.makedirs(DATA_DIR, exist_ok=True)
    print("Connexion en direct à l'API OpenData Paris (Disponibilité temps réel)...")
    try:
        r = requests.get(OPENDATA_URL, timeout=15)
        if r.status_code == 200:
            data = r.json()
            formatted = []
            for s in data:
                coords = s.get('coordonnees_geo') or {}
                lat, lon = coords.get('lat'), coords.get('lon')
                if lat and lon:
                    formatted.append({
                        'id': str(s.get('stationcode', '')),
                        'nom': s.get('name', 'Station Vélib'),
                        'latitude': float(lat),
                        'longitude': float(lon),
                        'capacite': int(s.get('capacity', 0)),
                        'numbikesavailable': int(s.get('numbikesavailable', 0)),
                        'numdocksavailable': int(s.get('numdocksavailable', 0)),
                        'ebike': int(s.get('ebike', 0)),
                        'mechanical': int(s.get('mechanical', 0)),
                        'is_renting': s.get('is_renting', 'OUI'),
                        'is_returning': s.get('is_returning', 'OUI'),
                        'commune': s.get('nom_arrondissement_communes', 'Île-de-France'),
                        'code_insee': s.get('code_insee_commune', ''),
                        'duedate': s.get('duedate', '')
                    })
            if formatted:
                with open(DATA_FILE, "w", encoding="utf-8") as f:
                    json.dump(formatted, f, ensure_ascii=False, indent=2)
                print(f"-> API Réussie : {len(formatted)} stations synchronisées en direct.")
                return formatted
    except Exception as e:
        print(f"Avertissement API Live ({e}). Chargement depuis le cache local...")

    if os.path.exists(DATA_FILE):
        with open(DATA_FILE, "r", encoding="utf-8") as f:
            return json.load(f)

    raise RuntimeError("Impossible d'obtenir les données Vélib en direct ou depuis le cache.")


def generer_graphiques_matplotlib(df, edges):
    """Génère le tableau de bord analytique en PNG avec Matplotlib."""
    plt.style.use('dark_background')
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Analyse Temps Réel du Réseau Vélib Île-de-France",
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


def generer_statistiques(df, edges, weight_mst, time_kruskal, time_prim):
    """Exporte les métriques du réseau au format JSON."""
    total_stations = len(df)
    total_capacite = int(df['capacite'].sum()) if 'capacite' in df else 0
    total_velos_dispo = int(df['numbikesavailable'].sum()) if 'numbikesavailable' in df else 0
    total_bornettes_libres = int(df['numdocksavailable'].sum()) if 'numdocksavailable' in df else 0
    total_ebikes = int(df['ebike'].sum()) if 'ebike' in df else 0

    stations_vides = int((df['numbikesavailable'] == 0).sum()) if 'numbikesavailable' in df else 0
    stations_saturees = int((df['numdocksavailable'] == 0).sum()) if 'numdocksavailable' in df else 0

    top_capacites = df.sort_values(by='capacite', ascending=False).head(10)[
        ['nom', 'commune', 'capacite', 'numbikesavailable']
    ].to_dict(orient='records')

    par_commune = df['commune'].value_counts().head(15).to_dict()
    distances = [e[2] for e in edges]

    rapport = {
        "metriques_generales": {
            "total_stations": total_stations,
            "total_communes": df['commune'].nunique(),
            "total_velos_dispo_temps_reel": total_velos_dispo,
            "total_velos_electriques": total_ebikes,
            "total_bornettes_libres_temps_reel": total_bornettes_libres,
            "capacite_totale_velos": total_capacite,
            "stations_penurie_vides": stations_vides,
            "stations_saturees_pleines": stations_saturees,
            "distance_mst_totale_km": round(weight_mst, 2),
            "nombre_connexions_delaunay": len(edges),
            "distance_inter_station_moyenne_km": round(sum(distances) / len(distances), 3),
            "temps_kruskal_ms": round(time_kruskal, 2),
            "temps_prim_ms": round(time_prim, 2)
        },
        "top_10_stations_capacite": top_capacites,
        "repartition_par_commune": par_commune
    }

    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        json.dump(rapport, f, ensure_ascii=False, indent=2)

    return rapport


def generer_carte_html_interactive(df, tri, mst_edges, edges, time_kruskal=5.5, time_prim=50.2, default_start_idx=0, default_target_idx=15):
    """Génère le Cockpit Data spatial Vélib Île-de-France complet."""
    center_lat = df["latitude"].mean()
    center_lon = df["longitude"].mean()
    m = folium.Map(
        location=[center_lat, center_lon],
        zoom_start=11,
        tiles="OpenStreetMap",
        prefer_canvas=True
    )

    # Calcul dynamique des métriques Cockpit & Santé Réseau
    total_stations = len(df)
    total_capacite = int(df['capacite'].sum()) if 'capacite' in df else 0
    total_velos_dispo = int(df['numbikesavailable'].sum()) if 'numbikesavailable' in df else 0
    total_ebikes = int(df['ebike'].sum()) if 'ebike' in df else 0
    total_bornettes_libres = int(df['numdocksavailable'].sum()) if 'numdocksavailable' in df else 0
    total_communes = df['commune'].nunique()
    weight_mst_sum = round(sum(w for _, _, w in mst_edges), 1)

    stations_vides_count = int((df['numbikesavailable'] == 0).sum()) if 'numbikesavailable' in df else 0
    stations_saturees_count = int((df['numdocksavailable'] == 0).sum()) if 'numdocksavailable' in df else 0
    taux_remplissage_pct = round((total_velos_dispo / max(1, total_capacite)) * 100, 1)

    # Add Live Availability Heatmap Layer
    heat_data = [[row['latitude'], row['longitude'], max(1, int(row.get('numbikesavailable', 0)))] for _, row in df.iterrows()]
    heatmap_group = folium.FeatureGroup(name="Carte de Chaleur (Disponibilité Vélos)", show=False)
    HeatMap(heat_data, radius=12, blur=15, max_zoom=13).add_to(heatmap_group)
    heatmap_group.add_to(m)

    marker_cluster = MarkerCluster(name=f"Stations Vélib ({total_stations:,})").add_to(m)

    stations_js_data = []
    coords_list = []
    for idx, row in df.iterrows():
        coords_list.append((row["longitude"], row["latitude"]))
        bikes_dispo = int(row.get('numbikesavailable', 0))
        docks_dispo = int(row.get('numdocksavailable', 0))
        ebikes = int(row.get('ebike', 0))
        mech = int(row.get('mechanical', 0))
        capa = int(row.get('capacite', 0))

        stations_js_data.append({
            "idx": idx,
            "id": row["id"],
            "nom": row["nom"],
            "lat": row["latitude"],
            "lon": row["longitude"],
            "capacite": capa,
            "bikes": bikes_dispo,
            "docks": docks_dispo,
            "ebike": ebikes,
            "mech": mech,
            "commune": row["commune"]
        })

        popup_html = f"""
        <div style="font-family: system-ui, -apple-system, sans-serif; min-width:210px;">
            <div style="font-weight:700; color:#1E40AF; font-size:14px; margin-bottom:4px;">{row['nom']}</div>
            <div style="color:#475569; font-size:12px; display:flex; align-items:center; gap:5px; margin-bottom:6px;">
                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M20 10c0 6-8 12-8 12s-8-6-8-12a8 8 0 0 1 16 0Z"/><circle cx="12" cy="10" r="3"/></svg>
                <span>Commune : <b>{row['commune']}</b></span>
            </div>
            
            <div style="background:#F8FAFC; border:1px solid #E2E8F0; border-radius:8px; padding:8px; margin:6px 0; font-size:11px;">
                <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:4px;">
                    <span style="color:#334155; font-weight:600;">🚲 Vélos dispo :</span>
                    <b style="color:#059669; font-size:13px;">{bikes_dispo}</b>
                </div>
                <div style="display:flex; justify-content:space-between; font-size:10px; color:#64748B; margin-bottom:4px; padding-left:8px;">
                    <span>⚡ Elec: <b>{ebikes}</b></span>
                    <span>🚲 Méca: <b>{mech}</b></span>
                </div>
                <div style="display:flex; justify-content:space-between; align-items:center; border-top:1px dashed #CBD5E1; padding-top:4px;">
                    <span style="color:#334155; font-weight:600;">🔌 Bornettes libres :</span>
                    <b style="color:#2563EB; font-size:12px;">{docks_dispo} / {capa}</b>
                </div>
            </div>

            <button onclick="selectStationByClick({idx})" style="
                width:100%; background:#2563EB; color:white; border:none; padding:7px 8px; border-radius:6px; font-weight:700; font-size:12px; cursor:pointer; display:flex; align-items:center; justify-content:center; gap:5px; margin-top:6px;
            ">
                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="10"/><line x1="12" y1="8" x2="12" y2="16"/><line x1="8" y1="12" x2="16" y2="12"/></svg>
                <span>Sélectionner pour itinéraire</span>
            </button>
        </div>
        """

        # Color indicator based on availability status
        marker_color = "#3182CE"
        if bikes_dispo == 0:
            marker_color = "#E53E3E"  # Red for empty alert
        elif ebikes > 5:
            marker_color = "#805AD5"  # Purple for high e-bike availability
        elif bikes_dispo > 10:
            marker_color = "#38A169"  # Green for high bike availability

        folium.CircleMarker(
            location=[row["latitude"], row["longitude"]],
            radius=5,
            popup=folium.Popup(popup_html, max_width=260),
            color=marker_color,
            fill=True,
            fill_color=marker_color,
            fill_opacity=0.85
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
            weight=0.9,
            fill=True,
            fill_color=color,
            fill_opacity=opacity,
            tooltip=f"Triangle Delaunay<br>Surface : <b>{area_str}</b><br>Densité : {'Forte' if opacity > 0.5 else 'Faible'}"
        ).add_to(delaunay_group)

    mst_group = folium.FeatureGroup(name=f"Réseau Optimal (MST - {weight_mst_sum} km)")
    for u, v, weight in mst_edges:
        loc1 = [df.loc[u, "latitude"], df.loc[u, "longitude"]]
        loc2 = [df.loc[v, "latitude"], df.loc[v, "longitude"]]
        folium.PolyLine(
            locations=[loc1, loc2],
            weight=2.8,
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

    delaunay_var_name = str(delaunay_group.get_name())
    edges_js = [{"u": int(u), "v": int(v), "w": round(float(w), 4)} for u, v, w in edges]

    top10_df = df.sort_values(by='capacite', ascending=False).head(8)
    chart_top_labels = top10_df['nom'].str[:20].tolist()
    chart_top_values = top10_df['capacite'].tolist()

    communes_top = df['commune'].value_counts().head(6)
    chart_commune_labels = communes_top.index.tolist()
    chart_commune_values = communes_top.values.tolist()

    dashboard_ui_html = f"""
    <link href="https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@400;500;600;700;800&display=swap" rel="stylesheet">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>

    <style>
    * {{
        font-family: 'Plus Jakarta Sans', system-ui, -apple-system, sans-serif !important;
        box-sizing: border-box;
    }}
    html, body {{
        width: 100vw !important;
        height: 100vh !important;
        margin: 0 !important;
        padding: 0 !important;
        overflow: hidden !important;
    }}
    .leaflet-container {{
        width: 100% !important;
        height: 100% !important;
    }}
    .leaflet-top.leaflet-left {{
        top: 70px !important;
        left: 15px !important;
    }}
    
    /* Animations & Glassmorphism Cockpit Design */
    @keyframes floatIn {{
        from {{ opacity: 0; transform: translateY(-12px); }}
        to {{ opacity: 1; transform: translateY(0); }}
    }}
    @keyframes pulseGlow {{
        0%, 100% {{ box-shadow: 0 8px 32px rgba(15, 23, 42, 0.4), inset 0 0 0 1px rgba(255, 255, 255, 0.1); }}
        50% {{ box-shadow: 0 12px 36px rgba(99, 102, 241, 0.25), inset 0 0 0 1px rgba(139, 92, 246, 0.3); }}
    }}

    .glass-panel {{
        background: rgba(15, 23, 42, 0.88) !important;
        backdrop-filter: blur(18px) saturate(180%) !important;
        -webkit-backdrop-filter: blur(18px) saturate(180%) !important;
        border: 1px solid rgba(255, 255, 255, 0.12) !important;
        box-shadow: 0 20px 40px rgba(0, 0, 0, 0.35) !important;
        color: #F8FAFC !important;
    }}

    .kpi-card {{
        background: rgba(15, 23, 42, 0.90);
        color: white;
        padding: 6px 14px;
        border-radius: 12px;
        backdrop-filter: blur(16px);
        border: 1px solid rgba(255, 255, 255, 0.1);
        display: flex;
        align-items: center;
        gap: 10px;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        animation: floatIn 0.5s ease-out forwards;
    }}
    .kpi-card:hover {{
        transform: translateY(-2px);
        border-color: rgba(99, 102, 241, 0.4);
        box-shadow: 0 8px 25px rgba(99, 102, 241, 0.25);
    }}

    .kpi-icon {{
        display: flex;
        align-items: center;
        justify-content: center;
        width: 32px;
        height: 32px;
        border-radius: 8px;
        background: rgba(255, 255, 255, 0.08);
        transition: transform 0.3s ease;
    }}
    .kpi-card:hover .kpi-icon {{
        transform: scale(1.1);
    }}

    .btn-gradient {{
        background: linear-gradient(135deg, #3B82F6 0%, #2563EB 100%) !important;
        color: white !important;
        border: none !important;
        padding: 8px 12px !important;
        border-radius: 8px !important;
        font-weight: 700 !important;
        cursor: pointer !important;
        display: flex !important;
        align-items: center !important;
        justify-content: center !important;
        gap: 6px !important;
        font-size: 11px !important;
        transition: all 0.25s ease !important;
        box-shadow: 0 4px 14px rgba(37, 99, 235, 0.35) !important;
    }}
    .btn-gradient:hover {{
        transform: translateY(-1px) !important;
        box-shadow: 0 6px 20px rgba(37, 99, 235, 0.5) !important;
        filter: brightness(1.1) !important;
    }}

    .btn-purple {{
        background: linear-gradient(135deg, #6366F1 0%, #8B5CF6 100%) !important;
        box-shadow: 0 4px 14px rgba(99, 102, 241, 0.35) !important;
    }}

    .btn-emerald {{
        background: linear-gradient(135deg, #10B981 0%, #059669 100%) !important;
        box-shadow: 0 4px 14px rgba(16, 185, 129, 0.35) !important;
    }}

    .custom-select {{
        width: 100%;
        padding: 7px 10px;
        margin: 4px 0 10px 0;
        border-radius: 8px;
        background: rgba(15, 23, 42, 0.95);
        color: #F8FAFC;
        border: 1px solid rgba(255, 255, 255, 0.15);
        font-size: 12px;
        outline: none;
        transition: border 0.2s;
    }}
    .custom-select:focus {{
        border-color: #6366F1;
        box-shadow: 0 0 0 2px rgba(99, 102, 241, 0.25);
    }}

    .tab-btn {{
        background: transparent;
        border: none;
        color: #94A3B8;
        padding: 6px 10px;
        font-size: 10px;
        font-weight: 700;
        cursor: pointer;
        border-bottom: 2px solid transparent;
        transition: all 0.2s;
    }}

    /* Responsive & Zoom-Proof CSS */
    #route-panel {{
        max-height: calc(100vh - 105px) !important;
        overflow-y: auto !important;
        width: min(310px, 88vw) !important;
    }}

    #route-panel::-webkit-scrollbar {{
        width: 5px;
    }}
    #route-panel::-webkit-scrollbar-thumb {{
        background: rgba(255, 255, 255, 0.2);
        border-radius: 4px;
    }}

    #kpi-banner {{
        max-width: calc(100vw - 640px) !important;
        flex-wrap: wrap !important;
        justify-content: center !important;
        top: 14px !important;
        left: 50% !important;
        transform: translateX(-50%) !important;
    }}

    #bottom-legend-banner {{
        max-width: min(650px, 90vw) !important;
        flex-wrap: wrap !important;
        justify-content: center !important;
        gap: 10px !important;
        left: 50% !important;
        transform: translateX(-50%) !important;
        bottom: 16px !important;
    }}

    @media (max-width: 1250px) {{
        #kpi-banner {{
            display: none !important;
        }}
    }}

    @media (max-height: 750px) {{
        #kpi-banner {{
            display: none !important;
        }}
        #route-panel {{
            top: 65px !important;
            max-height: calc(100vh - 85px) !important;
        }}
        #bottom-legend-banner {{
            bottom: 10px !important;
            padding: 4px 12px !important;
        }}
    }}

    @media (max-width: 650px) {{
        #cockpit-search-bar {{
            width: calc(100vw - 30px) !important;
        }}
        #route-panel {{
            width: calc(100vw - 30px) !important;
            right: 15px !important;
            left: 15px !important;
            top: 65px !important;
            max-height: calc(100vh - 120px) !important;
        }}
        #bottom-legend-banner {{
            width: calc(100vw - 30px) !important;
            bottom: 10px !important;
            font-size: 10px !important;
            padding: 6px 12px !important;
            border-radius: 14px !important;
        }}
    }}
    </style>

    <!-- Cockpit Live Search Bar (Top Left Panel) -->
    <div id="cockpit-search-bar" class="glass-panel" style="
        position: fixed; top: 14px; left: 15px; width: min(270px, 80vw); z-index: 9999; padding: 6px 12px; border-radius: 12px; display: flex; align-items: center; gap: 8px;
    ">
        <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="#94A3B8" stroke-width="2.2"><circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/></svg>
        <input id="search-input" onkeyup="filterCockpitSearch(this.value)" placeholder="Chercher une station..." style="
            background: transparent; border: none; outline: none; color: #F8FAFC; font-size: 12px; width: 100%; font-weight: 600;
        " />
        <div id="search-results" style="
            position: absolute; top: 40px; left: 0; right: 0; background: rgba(15, 23, 42, 0.95); border: 1px solid rgba(255, 255, 255, 0.15); border-radius: 10px; max-height: 200px; overflow-y: auto; display: none; backdrop-filter: blur(12px);
        "></div>
    </div>

    <!-- Top KPI Banner (Métriques Dynamiques Cockpit Temps Réel) -->
    <div id="kpi-banner" style="
        position: fixed; top: 14px; left: 50%; transform: translateX(-50%); display: flex; gap: 10px; z-index: 9999;
    ">
        <div class="kpi-card">
            <div class="kpi-icon" style="color: #38BDF8; background: rgba(56, 189, 248, 0.15);">
                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2"><path d="M20 10c0 6-8 12-8 12s-8-6-8-12a8 8 0 0 1 16 0Z"/><circle cx="12" cy="10" r="3"/></svg>
            </div>
            <div>
                <div style="font-size: 9px; color: #94A3B8; text-transform: uppercase; font-weight: 800; letter-spacing: 0.5px;">Stations</div>
                <div style="font-size: 15px; font-weight: 800; color: #F8FAFC;">{total_stations:,}</div>
            </div>
        </div>
        <div class="kpi-card">
            <div class="kpi-icon" style="color: #34D399; background: rgba(52, 211, 153, 0.15);">
                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2"><circle cx="5.5" cy="17.5" r="3.5"/><circle cx="18.5" cy="17.5" r="3.5"/><path d="M15 6a1 1 0 1 0 0-2 1 1 0 0 0 0 2zm-3 11.5L9.5 10l-3 3.5M12 17.5V10l3.5-4H18"/></svg>
            </div>
            <div>
                <div style="font-size: 9px; color: #94A3B8; text-transform: uppercase; font-weight: 800; letter-spacing: 0.5px; display:flex; align-items:center; gap:4px;">
                    <span>Vélos Dispo (Live)</span>
                    <span style="width:6px; height:6px; border-radius:50%; background:#34D399; box-shadow:0 0 6px #34D399;"></span>
                </div>
                <div style="font-size: 15px; font-weight: 800; color: #F8FAFC;">{total_velos_dispo:,}</div>
            </div>
        </div>
        <div class="kpi-card">
            <div class="kpi-icon" style="color: #FBBF24; background: rgba(251, 191, 36, 0.15);">
                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2"><circle cx="6" cy="19" r="3"/><circle cx="18" cy="5" r="3"/><circle cx="12" cy="6" r="3"/><path d="M8.5 17l2.5-8.5M15.5 17l-2.5-8.5"/></svg>
            </div>
            <div>
                <div style="font-size: 9px; color: #94A3B8; text-transform: uppercase; font-weight: 800; letter-spacing: 0.5px;">Réseau MST</div>
                <div style="font-size: 15px; font-weight: 800; color: #F8FAFC;">{weight_mst_sum} km</div>
            </div>
        </div>
        <div class="kpi-card">
            <div class="kpi-icon" style="color: #A855F7; background: rgba(168, 85, 247, 0.15);">
                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2"><path d="M13 2L3 14h9l-1 8 10-12h-9l1-8z"/></svg>
            </div>
            <div>
                <div style="font-size: 9px; color: #94A3B8; text-transform: uppercase; font-weight: 800; letter-spacing: 0.5px;">Vélos Élec.</div>
                <div style="font-size: 15px; font-weight: 800; color: #F8FAFC;">{total_ebikes:,}</div>
            </div>
        </div>
    </div>

    <!-- Compact Route Calculator Panel (Top Right - Glassmorphism) -->
    <div id="route-panel" class="glass-panel" style="
        position: fixed; top: 75px; right: 15px; width: min(310px, 90vw); border-radius: 16px; padding: 16px; z-index: 9999; font-size: 12px; animation: floatIn 0.6s ease-out forwards;
    ">
        <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px;">
            <div style="font-weight:800; color:#F8FAFC; font-size:14px; display:flex; align-items:center; gap:8px;">
                <div style="width:24px; height:24px; border-radius:6px; background:rgba(99, 102, 241, 0.2); display:flex; align-items:center; justify-content:center; color:#818CF8;">
                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2"><circle cx="12" cy="12" r="10"/><path d="m16.24 7.76-2.12 6.36-6.36 2.12 2.12-6.36z"/></svg>
                </div>
                <span>Cockpit Itinéraire</span>
            </div>
            <button onclick="toggleStatsDrawer()" style="background:rgba(255,255,255,0.08); border:1px solid rgba(255,255,255,0.15); color:#94A3B8; padding:4px 10px; border-radius:8px; font-size:11px; font-weight:700; cursor:pointer; display:flex; align-items:center; gap:4px; transition:all 0.2s;">
                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><line x1="18" y1="20" x2="18" y2="10"/><line x1="12" y1="20" x2="12" y2="4"/><line x1="6" y1="20" x2="6" y2="14"/></svg>
                <span>Analytics</span>
            </button>
        </div>
        
        <label style="font-weight:700; color:#94A3B8; font-size:11px; text-transform:uppercase; letter-spacing:0.4px;">Départ</label>
        <select id="start-station" class="custom-select" onchange="updateStationComparison()"></select>

        <label style="font-weight:700; color:#94A3B8; font-size:11px; text-transform:uppercase; letter-spacing:0.4px;">Étape / Waypoint (Optionnel)</label>
        <select id="waypoint-station" class="custom-select" onchange="updateStationComparison()">
            <option value="-1">-- Aucune étape (Direct) --</option>
        </select>
        
        <label style="font-weight:700; color:#94A3B8; font-size:11px; text-transform:uppercase; letter-spacing:0.4px;">Arrivée</label>
        <select id="target-station" class="custom-select" onchange="updateStationComparison()"></select>
        
        <button onclick="calculateRoute()" class="btn-gradient" style="width:100%; margin-top:4px;">
            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2"><circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/></svg>
            <span>Calculer trajet (Dijkstra)</span>
        </button>

        <div id="route-results" style="margin-top:12px; display:none; padding:12px; background:rgba(0,0,0,0.3); border-radius:10px; border:1px solid rgba(255,255,255,0.08);">
            <div style="display:flex; justify-content:space-between; margin-bottom:5px;">
                <span style="color:#94A3B8;">Distance :</span>
                <b id="route-dist" style="color:#38BDF8;">-</b>
            </div>
            <div style="display:flex; justify-content:space-between; margin-bottom:5px;">
                <span style="color:#94A3B8;">Temps (15km/h) :</span>
                <b id="route-time" style="color:#34D399;">-</b>
            </div>
            <div style="display:flex; justify-content:space-between; margin-bottom:5px;">
                <span style="color:#94A3B8;">Stations traversées :</span>
                <b id="route-hops" style="color:#F8FAFC;">-</b>
            </div>
            <div style="display:flex; justify-content:space-between; border-top:1px dashed rgba(255,255,255,0.1); padding-top:4px;">
                <span style="color:#94A3B8;">Écon. CO₂ (vs Auto) :</span>
                <b id="route-co2" style="color:#FBBF24;">-</b>
            </div>
        </div>

        <!-- Floating Comparison Card Side-by-Side -->
        <div id="station-comparison-card" style="margin-top:12px; padding:10px; background:rgba(255,255,255,0.04); border-radius:10px; border:1px solid rgba(255,255,255,0.08); font-size:11px;">
            <div style="color:#38BDF8; font-weight:800; text-transform:uppercase; font-size:10px; margin-bottom:6px; display:flex; justify-content:space-between;">
                <span>Comparateur Départ / Arrivée</span>
                <span style="color:#94A3B8;">Direct</span>
            </div>
            <div style="display:grid; grid-template-columns: 1fr 1fr; gap:8px;">
                <div id="comp-start-box" style="background:rgba(0,0,0,0.25); padding:6px; border-radius:6px;">
                    <div id="comp-start-name" style="font-weight:700; color:#F8FAFC; font-size:10px; white-space:nowrap; overflow:hidden; text-overflow:ellipsis;">Départ</div>
                    <div style="color:#34D399; font-weight:800; margin-top:2px;">🚲 <span id="comp-start-bikes">-</span> dispo</div>
                    <div style="color:#94A3B8; font-size:9px;">⚡ <span id="comp-start-ebike">-</span> elec</div>
                </div>
                <div id="comp-target-box" style="background:rgba(0,0,0,0.25); padding:6px; border-radius:6px;">
                    <div id="comp-target-name" style="font-weight:700; color:#F8FAFC; font-size:10px; white-space:nowrap; overflow:hidden; text-overflow:ellipsis;">Arrivée</div>
                    <div style="color:#38BDF8; font-weight:800; margin-top:2px;">🔌 <span id="comp-target-docks">-</span> docks</div>
                    <div style="color:#94A3B8; font-size:9px;">⚡ <span id="comp-target-ebike">-</span> elec</div>
                </div>
            </div>
        </div>

        <!-- Collapsible Multi-Tab Cockpit Analytics Drawer inside panel -->
        <div id="stats-drawer" style="display:none; margin-top:14px; border-top:1px solid rgba(255,255,255,0.1); padding-top:12px;">
            <div style="display:flex; border-bottom:1px solid rgba(255,255,255,0.1); margin-bottom:10px;">
                <button class="tab-btn active" onclick="switchCockpitTab('tab-charts', this)">Graphiques</button>
                <button class="tab-btn" onclick="switchCockpitTab('tab-algos', this)">Performances</button>
                <button class="tab-btn" onclick="switchCockpitTab('tab-health', this)">Santé Réseau</button>
            </div>

            <div id="tab-charts">
                <div style="font-size:11px; font-weight:800; color:#94A3B8; text-transform:uppercase; letter-spacing:0.5px; margin-bottom:6px;">Top Capacités</div>
                <canvas id="chartTopCapacity" height="130"></canvas>

                <div style="font-size:11px; font-weight:800; color:#94A3B8; text-transform:uppercase; letter-spacing:0.5px; margin:10px 0 6px 0;">Répartition Communes</div>
                <canvas id="chartCommunes" height="130"></canvas>
            </div>

            <div id="tab-algos" style="display:none; font-size:11px;">
                <div style="color:#94A3B8; font-weight:800; text-transform:uppercase; margin-bottom:6px;">Comparatif Algorithmes MST</div>
                <div style="background:rgba(0,0,0,0.3); padding:8px; border-radius:8px; border:1px solid rgba(255,255,255,0.08); margin-bottom:8px;">
                    <div style="display:flex; justify-content:space-between; margin-bottom:4px;">
                        <span>Kruskal (Union-Find) :</span>
                        <b style="color:#34D399;">{round(time_kruskal, 2)} ms</b>
                    </div>
                    <div style="display:flex; justify-content:space-between; margin-bottom:4px;">
                        <span>Prim (Min-Heap) :</span>
                        <b style="color:#FBBF24;">{round(time_prim, 2)} ms</b>
                    </div>
                    <div style="display:flex; justify-content:space-between;">
                        <span>Distance MST Totale :</span>
                        <b style="color:#38BDF8;">{weight_mst_sum} km</b>
                    </div>
                </div>

                <div style="color:#94A3B8; font-weight:800; text-transform:uppercase; margin-bottom:6px;">Réduction Maillage Delaunay</div>
                <div style="background:rgba(0,0,0,0.3); padding:8px; border-radius:8px; border:1px solid rgba(255,255,255,0.08);">
                    <div style="display:flex; justify-content:space-between; margin-bottom:4px;">
                        <span>Graphe Complet V(V-1)/2 :</span>
                        <b style="color:#F87171;">1 151 403 arêtes</b>
                    </div>
                    <div style="display:flex; justify-content:space-between;">
                        <span>Graphe Delaunay :</span>
                        <b style="color:#34D399;">{len(edges):,} arêtes</b>
                    </div>
                </div>
            </div>

            <div id="tab-health" style="display:none; font-size:11px;">
                <div style="color:#94A3B8; font-weight:800; text-transform:uppercase; margin-bottom:6px;">Santé & Satiété Réseau</div>
                <div style="background:rgba(0,0,0,0.3); padding:8px; border-radius:8px; border:1px solid rgba(255,255,255,0.08); margin-bottom:8px;">
                    <div style="display:flex; justify-content:space-between; margin-bottom:4px;">
                        <span>Taux d'occupation global :</span>
                        <b style="color:#34D399;">{taux_remplissage_pct}%</b>
                    </div>
                    <div style="display:flex; justify-content:space-between; margin-bottom:4px;">
                        <span>Stations en alerte pénurie :</span>
                        <b style="color:#EF4444;">{stations_vides_count}</b>
                    </div>
                    <div style="display:flex; justify-content:space-between;">
                        <span>Stations saturées (0 dock) :</span>
                        <b style="color:#F59E0B;">{stations_saturees_count}</b>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <!-- Mini Bottom Banner Legend (Glossy Pill with Subtle Glow) -->
    <div id="bottom-legend-banner" class="glass-panel" style="
        position: fixed; bottom: 22px; left: 50%; transform: translateX(-50%); display: flex; align-items: center; gap: 16px; padding: 8px 20px; border-radius: 9999px; z-index: 9999; font-size: 11px; animation: pulseGlow 4s infinite ease-in-out;
    ">
        <div style="display:flex; align-items:center; gap:6px; font-weight:800; color:#94A3B8; border-right:1px solid rgba(255,255,255,0.15); padding-right:12px; letter-spacing:0.5px;">
            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="#A855F7" stroke-width="2.2"><path d="M12 2 2 22h20L12 2z"/></svg>
            <span>DENSITÉ DELAUNAY</span>
        </div>

        <div style="display:flex; align-items:center; gap:6px;">
            <span style="width:12px; height:12px; border-radius:4px; background:#0F172A; border:1px solid rgba(255,255,255,0.3); box-shadow:0 0 8px rgba(15,23,42,0.8);"></span>
            <span style="color:#F8FAFC; font-weight:700;">Haute (Sombre)</span>
        </div>

        <div style="display:flex; align-items:center; gap:6px;">
            <span style="width:12px; height:12px; border-radius:4px; background:#4338CA; box-shadow:0 0 8px rgba(67,56,202,0.6);"></span>
            <span style="color:#CBD5E1; font-weight:600;">Moyenne</span>
        </div>

        <div style="display:flex; align-items:center; gap:6px;">
            <span style="width:12px; height:12px; border-radius:4px; background:#E0E7FF; border:1px solid rgba(255,255,255,0.5); box-shadow:0 0 8px rgba(224,231,255,0.5);"></span>
            <span style="color:#CBD5E1; font-weight:600;">Faible (Claire)</span>
        </div>

        <button id="toggle-delaunay-btn" onclick="toggleDelaunayLayer()" class="btn-gradient btn-purple" style="
            margin-left:4px; padding:5px 14px !important; border-radius:9999px !important; font-size:11px !important;
        ">
            <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2"><path d="M12 2 2 22h20L12 2z"/></svg>
            <span>Delaunay ON/OFF</span>
        </button>
    </div>

    <script>
    const STATIONS = {json.dumps(stations_js_data)};
    const EDGES = {json.dumps(edges_js)};
    const DELAUNAY_VAR_NAME = "{delaunay_var_name}";

    let activeRouteLayer = null;
    let clickSelectionStep = 0;
    let chartTopCapacityInst = null;
    let chartCommunesInst = null;

    function getLeafletMap() {{
        return Object.values(window).find(v => v && v.fitBounds && v.addLayer && v.removeLayer && v._layers);
    }}

    function filterCockpitSearch(query) {{
        const box = document.getElementById("search-results");
        if (!box) return;
        if (!query || query.trim().length < 2) {{
            box.style.display = "none";
            return;
        }}
        const q = query.toLowerCase();
        const matches = STATIONS.filter(s => s.nom.toLowerCase().includes(q) || s.commune.toLowerCase().includes(q)).slice(0, 6);
        if (matches.length === 0) {{
            box.style.display = "none";
            return;
        }}
        box.innerHTML = "";
        matches.forEach(m => {{
            const div = document.createElement("div");
            div.style.padding = "6px 10px";
            div.style.cursor = "pointer";
            div.style.fontSize = "11px";
            div.style.borderBottom = "1px solid rgba(255,255,255,0.05)";
            div.innerHTML = "<b style='color:#38BDF8;'>" + m.nom + "</b> <span style='color:#94A3B8;'>(" + m.commune + ")</span>";
            div.onclick = function() {{
                box.style.display = "none";
                document.getElementById("search-input").value = m.nom;
                const mapObj = getLeafletMap();
                if (mapObj) mapObj.flyTo([m.lat, m.lon], 16);
            }};
            box.appendChild(div);
        }});
        box.style.display = "block";
    }}

    function switchCockpitTab(tabId, btn) {{
        document.getElementById("tab-charts").style.display = "none";
        document.getElementById("tab-algos").style.display = "none";
        document.getElementById("tab-health").style.display = "none";
        document.querySelectorAll(".tab-btn").forEach(b => b.classList.remove("active"));
        document.getElementById(tabId).style.display = "block";
        btn.classList.add("active");
        if (tabId === 'tab-charts') setTimeout(initCharts, 50);
    }}

    function updateStationComparison() {{
        const selectStart = document.getElementById("start-station");
        const selectTarget = document.getElementById("target-station");
        if (!selectStart || !selectTarget) return;

        const uStart = parseInt(selectStart.value);
        const uTarget = parseInt(selectTarget.value);

        if (!isNaN(uStart) && STATIONS[uStart]) {{
            const s = STATIONS[uStart];
            document.getElementById("comp-start-name").textContent = s.nom;
            document.getElementById("comp-start-bikes").textContent = s.bikes;
            document.getElementById("comp-start-ebike").textContent = s.ebike;
        }}

        if (!isNaN(uTarget) && STATIONS[uTarget]) {{
            const t = STATIONS[uTarget];
            document.getElementById("comp-target-name").textContent = t.nom;
            document.getElementById("comp-target-docks").textContent = t.docks;
            document.getElementById("comp-target-ebike").textContent = t.ebike;
        }}
    }}

    function initVelibApp() {{
        const selectStart = document.getElementById("start-station");
        const selectTarget = document.getElementById("target-station");
        const selectWaypoint = document.getElementById("waypoint-station");
        if (!selectStart || !selectTarget) return;

        selectStart.innerHTML = "";
        selectTarget.innerHTML = "";
        if (selectWaypoint) selectWaypoint.innerHTML = "<option value='-1'>-- Aucune étape (Direct) --</option>";
        
        STATIONS.forEach(s => {{
            let opt1 = document.createElement("option");
            opt1.value = s.idx;
            opt1.textContent = s.nom + " (" + s.commune + ")";
            selectStart.appendChild(opt1);
            
            let opt2 = document.createElement("option");
            opt2.value = s.idx;
            opt2.textContent = s.nom + " (" + s.commune + ")";
            selectTarget.appendChild(opt2);

            if (selectWaypoint) {{
                let opt3 = document.createElement("option");
                opt3.value = s.idx;
                opt3.textContent = s.nom + " (" + s.commune + ")";
                selectWaypoint.appendChild(opt3);
            }}
        }});
        
        selectStart.selectedIndex = {default_start_idx};
        selectTarget.selectedIndex = {default_target_idx};
        updateStationComparison();
    }}

    function initCharts() {{
        if (chartTopCapacityInst && chartCommunesInst) {{
            chartTopCapacityInst.resize();
            chartCommunesInst.resize();
            return;
        }}

        const canvas1 = document.getElementById('chartTopCapacity');
        if (canvas1 && !chartTopCapacityInst) {{
            chartTopCapacityInst = new Chart(canvas1, {{
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
                    maintainAspectRatio: true,
                    plugins: {{
                        legend: {{ display: false }}
                    }},
                    scales: {{
                        x: {{ ticks: {{ color: '#94A3B8', font: {{ size: 9 }} }} }},
                        y: {{ ticks: {{ color: '#94A3B8', font: {{ size: 9 }} }}, beginAtZero: true }}
                    }}
                }}
            }});
        }}

        const canvas2 = document.getElementById('chartCommunes');
        if (canvas2 && !chartCommunesInst) {{
            chartCommunesInst = new Chart(canvas2, {{
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
                    maintainAspectRatio: true,
                    plugins: {{
                        legend: {{
                            position: 'right',
                            labels: {{ color: '#CBD5E1', boxWidth: 10, font: {{ size: 10 }} }}
                        }}
                    }}
                }}
            }});
        }}
    }}

    if (document.readyState === "loading") {{
        document.addEventListener("DOMContentLoaded", initVelibApp);
    }} else {{
        initVelibApp();
    }}

    function toggleStatsDrawer() {{
        const drawer = document.getElementById("stats-drawer");
        if (!drawer) return;
        if (drawer.style.display === "none" || !drawer.style.display) {{
            drawer.style.display = "block";
            setTimeout(initCharts, 50);
        }} else {{
            drawer.style.display = "none";
        }}
    }}

    function toggleDelaunayLayer() {{
        const mapObj = getLeafletMap();
        const btn = document.getElementById("toggle-delaunay-btn");
        const delaunayLayer = window[DELAUNAY_VAR_NAME];

        if (!mapObj || !delaunayLayer) {{
            console.warn("Delaunay layer non prêt");
            return;
        }}

        if (mapObj.hasLayer(delaunayLayer)) {{
            mapObj.removeLayer(delaunayLayer);
            if (btn) btn.style.opacity = "0.4";
        }} else {{
            mapObj.addLayer(delaunayLayer);
            if (btn) btn.style.opacity = "1.0";
        }}
    }}

    function selectStationByClick(stationIdx) {{
        const selectStart = document.getElementById("start-station");
        const selectTarget = document.getElementById("target-station");
        if (!selectStart || !selectTarget) return;

        if (clickSelectionStep === 0 || clickSelectionStep === 2) {{
            selectStart.value = stationIdx;
            clickSelectionStep = 1;
            updateStationComparison();
            alert("Départ sélectionné : " + STATIONS[stationIdx].nom + "\\n\\nCliquez sur une 2ème station pour l'arrivée ou cliquez sur 'Calculer trajet'.");
        }} else if (clickSelectionStep === 1) {{
            selectTarget.value = stationIdx;
            clickSelectionStep = 2;
            updateStationComparison();
            calculateRoute();
        }}
    }}

    function dijkstraPathBetween(uStart, uTarget) {{
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
            if (u === -1 || dist[u] === Infinity || u === uTarget) break;
            visited[u] = true;
            
            adj[u].forEach(edge => {{
                if (dist[u] + edge.w < dist[edge.node]) {{
                    dist[edge.node] = dist[u] + edge.w;
                    parent[edge.node] = u;
                }}
            }});
        }}
        
        if (dist[uTarget] === Infinity) return null;

        const path = [];
        let curr = uTarget;
        while (curr !== null) {{
            path.push(curr);
            curr = parent[curr];
        }}
        path.reverse();
        return {{ path, dist: dist[uTarget] }};
    }}

    function calculateRoute() {{
        const selectStart = document.getElementById("start-station");
        const selectTarget = document.getElementById("target-station");
        const selectWaypoint = document.getElementById("waypoint-station");
        if (!selectStart || !selectTarget) return;

        const uStart = parseInt(selectStart.value);
        const uTarget = parseInt(selectTarget.value);
        const uWaypoint = selectWaypoint ? parseInt(selectWaypoint.value) : -1;
        
        if (isNaN(uStart) || isNaN(uTarget)) {{
            alert("Veuillez sélectionner une station de départ et d'arrivée.");
            return;
        }}

        if (uStart === uTarget) {{
            alert("Veuillez choisir deux stations différentes.");
            return;
        }}

        updateStationComparison();

        let fullPath = [];
        let totalDistKm = 0;

        if (uWaypoint !== -1 && uWaypoint !== uStart && uWaypoint !== uTarget) {{
            const leg1 = dijkstraPathBetween(uStart, uWaypoint);
            const leg2 = dijkstraPathBetween(uWaypoint, uTarget);
            if (!leg1 || !leg2) {{
                alert("Aucun itinéraire trouvé passant par cette étape.");
                return;
            }}
            fullPath = leg1.path.concat(leg2.path.slice(1));
            totalDistKm = leg1.dist + leg2.dist;
        }} else {{
            const res = dijkstraPathBetween(uStart, uTarget);
            if (!res) {{
                alert("Aucun itinéraire trouvé entre ces deux stations.");
                return;
            }}
            fullPath = res.path;
            totalDistKm = res.dist;
        }}

        const minutes = Math.max(1, Math.round((totalDistKm / 15) * 60));
        const co2Saved = Math.round(totalDistKm * 120);
        
        const resultsBox = document.getElementById("route-results");
        if (resultsBox) resultsBox.style.display = "block";

        document.getElementById("route-dist").textContent = totalDistKm.toFixed(2) + " km";
        document.getElementById("route-time").textContent = minutes + " min";
        document.getElementById("route-hops").textContent = fullPath.length + " stations";
        document.getElementById("route-co2").textContent = co2Saved + " g CO₂";
        
        const mapObj = getLeafletMap();
        if (mapObj) {{
            if (activeRouteLayer) mapObj.removeLayer(activeRouteLayer);
            
            const routeCoords = fullPath.map(idx => [STATIONS[idx].lat, STATIONS[idx].lon]);
            activeRouteLayer = L.polyline(routeCoords, {{
                color: '#EF4444',
                weight: 6,
                opacity: 0.95,
                dashArray: '8, 8'
            }}).addTo(mapObj);
            
            mapObj.fitBounds(activeRouteLayer.getBounds(), {{ padding: [60, 60] }});
        }}
    }}
    </script>
    """

    # Placer l'UI personnalisée tout à la fin pour s'assurer que Folium et Leaflet ont totalement initialisé leurs variables
    html_content = html_content + f"\n{dashboard_ui_html}\n"
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
    generer_statistiques(df, edges, weight_kruskal, time_kruskal, time_prim)

    # 5. Carte HTML Folium Cockpit
    generer_carte_html_interactive(df, tri, mst_kruskal, edges, time_kruskal, time_prim, args.depart or 0, args.arrivee or 15)
    print(f"\nCockpit interactif généré avec succès : {OUTPUT_MAP}")


if __name__ == "__main__":
    main()
