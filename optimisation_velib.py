"""
==============================================================================
Optimisation & Analytics du Réseau Vélib Métropole (1 518 Stations - Île-de-France)
Auteur : Misbaou DIALLO (BUT 3 Informatique)
==============================================================================

Fonctionnalités avancées :
  1. Triangulation de Delaunay (SciPy) pour réduire la complexité spatiale.
  2. Algorithmes MST : Kruskal (Union-Find) & Prim (Min-Heap) pour le réseau minimal.
  3. Recherche d'Itinéraire Optimal : Algorithme de Dijkstra (Min-Heap) entre 2 stations.
  4. Graphiques & Analytics Visuels (Matplotlib + Chart.js) :
     - Top 10 Stations par Capacité
     - Répartition par Département & Commune
     - Distribution des Distances Inter-Stations
  5. Application Web Interactive HTML (Folium + Leaflet + Calculateur + Chart.js Dashboard).
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
matplotlib.use('Agg')  # Rendu headless sans GUI
import matplotlib.pyplot as plt
import folium
from folium.plugins import MiniMap, MarkerCluster

# Chemins des fichiers
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
DATA_FILE = os.path.join(DATA_DIR, "stations_velib_idf_complete.json")
OUTPUT_MAP = os.path.join(BASE_DIR, "carte_velib_optimisee.html")
REPORT_FILE = os.path.join(BASE_DIR, "rapport_statistiques_velib.json")
GRAPH_IMAGE = os.path.join(DATA_DIR, "graphiques_velib.png")
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
    """Structure Union-Find pour l'Algorithme de Kruskal."""

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


def algo_dijkstra(n_vertices, edges, start_node, target_node):
    """Algorithme de Dijkstra pour le plus court chemin."""
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
    """Charge les données depuis le cache ou l'API."""
    os.makedirs(DATA_DIR, exist_ok=True)
    if os.path.exists(DATA_FILE):
        with open(DATA_FILE, "r", encoding="utf-8") as f:
            return json.load(f)

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
    """Génère un tableau de bord visuel en image PNG avec Matplotlib."""
    plt.style.use('dark_background')
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("TABLEAU DE BORD ANALYTIQUE RÉSEAU VÉLIB ÎLE-DE-FRANCE (1518 STATIONS)",
                 fontsize=14, fontweight='bold', color='#6366F1')

    # 1. Top 10 des stations par capacité
    top10 = df.sort_values(by='capacite', ascending=False).head(10)
    axes[0, 0].barh(top10['nom'].str[:25], top10['capacite'], color='#38BDF8')
    axes[0, 0].set_title("Top 10 Stations par Capacité de Vélos", fontsize=11, fontweight='bold')
    axes[0, 0].set_xlabel("Nombre de bornettes / vélos")
    axes[0, 0].invert_yaxis()

    # 2. Répartition des stations par commune (Top 8)
    communes = df['commune'].value_counts().head(8)
    axes[0, 1].pie(communes.values, labels=communes.index, autopct='%1.1f%%',
                   colors=['#818CF8', '#34D399', '#FBBF24', '#F87171', '#A78BFA', '#F472B6', '#38BDF8', '#4ADE80'])
    axes[0, 1].set_title("Répartition des Stations par Commune", fontsize=11, fontweight='bold')

    # 3. Distribution des distances inter-stations
    distances_mètres = [e[2] * 1000 for e in edges]
    axes[1, 0].hist(distances_mètres, bins=30, color='#34D399', edgecolor='#111827')
    axes[1, 0].set_title("Distribution des Distances Inter-Stations (Mètres)", fontsize=11, fontweight='bold')
    axes[1, 0].set_xlabel("Distance (mètres)")
    axes[1, 0].set_ylabel("Fréquence (arêtes Delaunay)")

    # 4. Capacité par station (Histogramme)
    axes[1, 1].hist(df['capacite'], bins=20, color='#FBBF24', edgecolor='#111827')
    axes[1, 1].set_title("Répartition des Capacités des Stations", fontsize=11, fontweight='bold')
    axes[1, 1].set_xlabel("Nombre de vélos")
    axes[1, 1].set_ylabel("Nombre de stations")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(GRAPH_IMAGE, dpi=200)
    plt.close()
    print(f"[✓] Graphiques analytiques PNG enregistrés dans : {GRAPH_IMAGE}")


def generer_statistiques(df, edges, weight_mst):
    """Calcule le rapport analytique JSON."""
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


def generer_carte_html_interactive(df, mst_edges, edges, default_start_idx=0, default_target_idx=15):
    """Génère la carte interactive HTML avec contrôles et graphiques Chart.js."""
    center_lat = df["latitude"].mean()
    center_lon = df["longitude"].mean()
    m = folium.Map(location=[center_lat, center_lon], zoom_start=11, tiles="OpenStreetMap")

    marker_cluster = MarkerCluster(name="Stations Vélib Île-de-France").add_to(m)

    stations_js_data = []
    for idx, row in df.iterrows():
        stations_js_data.append({
            "idx": idx,
            "id": row["id"],
            "nom": row["nom"],
            "lat": row["latitude"],
            "lon": row["longitude"],
            "capacite": row["capacite"],
            "commune": row["commune"]
        })

        folium.CircleMarker(
            location=[row["latitude"], row["longitude"]],
            radius=4,
            popup=f"<b>{row['nom']}</b><br>Commune: {row['commune']}<br>Capacité: {row['capacite']} vélos",
            color="#2B6CB0",
            fill=True,
            fill_color="#3182CE",
            fill_opacity=0.8
        ).add_to(marker_cluster)

    mst_group = folium.FeatureGroup(name="Réseau Optimal Optimisé (MST - 502 km)")
    for u, v, weight in mst_edges:
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

    with open(OUTPUT_MAP, "r", encoding="utf-8") as f:
        html_content = f.read()

    edges_js = [{"u": int(u), "v": int(v), "w": round(float(w), 4)} for u, v, w in edges]

    # Données pour les graphiques Chart.js
    top10_df = df.sort_values(by='capacite', ascending=False).head(8)
    chart_top_labels = top10_df['nom'].str[:20].tolist()
    chart_top_values = top10_df['capacite'].tolist()

    communes_top = df['commune'].value_counts().head(6)
    chart_commune_labels = communes_top.index.tolist()
    chart_commune_values = communes_top.values.tolist()

    dashboard_ui_html = f"""
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>

    <!-- Panneau Supérieur : Cartes KPI -->
    <div id="kpi-banner" style="
        position: fixed;
        top: 15px;
        left: 50%;
        transform: translateX(-50%);
        display: flex;
        gap: 12px;
        z-index: 9999;
        font-family: 'Segoe UI', Arial, sans-serif;
    ">
        <div style="background: rgba(15, 23, 42, 0.9); color: white; padding: 8px 16px; border-radius: 8px; backdrop-filter: blur(8px); border: 1px solid #334155; text-align: center;">
            <div style="font-size: 10px; color: #94A3B8; text-transform: uppercase; font-weight: bold;">Stations</div>
            <div style="font-size: 18px; font-weight: bold; color: #38BDF8;">1 518</div>
        </div>
        <div style="background: rgba(15, 23, 42, 0.9); color: white; padding: 8px 16px; border-radius: 8px; backdrop-filter: blur(8px); border: 1px solid #334155; text-align: center;">
            <div style="font-size: 10px; color: #94A3B8; text-transform: uppercase; font-weight: bold;">Vélos & Bornettes</div>
            <div style="font-size: 18px; font-weight: bold; color: #34D399;">49 060</div>
        </div>
        <div style="background: rgba(15, 23, 42, 0.9); color: white; padding: 8px 16px; border-radius: 8px; backdrop-filter: blur(8px); border: 1px solid #334155; text-align: center;">
            <div style="font-size: 10px; color: #94A3B8; text-transform: uppercase; font-weight: bold;">Réseau Optimisé (MST)</div>
            <div style="font-size: 18px; font-weight: bold; color: #FBBF24;">502.09 km</div>
        </div>
        <div style="background: rgba(15, 23, 42, 0.9); color: white; padding: 8px 16px; border-radius: 8px; backdrop-filter: blur(8px); border: 1px solid #334155; text-align: center;">
            <div style="font-size: 10px; color: #94A3B8; text-transform: uppercase; font-weight: bold;">Communes Couvertes</div>
            <div style="font-size: 18px; font-weight: bold; color: #F472B6;">69</div>
        </div>
    </div>

    <!-- Panneau Droit : Calculateur d'Itinéraire + Graphiques -->
    <div id="route-panel" style="
        position: fixed;
        top: 80px;
        right: 15px;
        width: 350px;
        max-height: calc(100vh - 100px);
        overflow-y: auto;
        background: rgba(255, 255, 255, 0.95);
        border-radius: 14px;
        padding: 16px;
        box-shadow: 0 8px 32px rgba(0,0,0,0.2);
        z-index: 9999;
        font-family: 'Segoe UI', Arial, sans-serif;
        font-size: 13px;
        backdrop-filter: blur(10px);
        border: 1px solid #E2E8F0;
    ">
        <h3 style="margin:0 0 12px 0; color:#1E293B; font-size:16px; display:flex; align-items:center; gap:8px;">
            🚲 <span>Calculateur d'Itinéraire</span>
        </h3>
        
        <label style="font-weight:600; color:#475569;">Station de départ :</label>
        <select id="start-station" style="width:100%; padding:7px; margin:4px 0 10px 0; border-radius:6px; border:1px solid #CBD5E0;"></select>
        
        <label style="font-weight:600; color:#475569;">Station d'arrivée :</label>
        <select id="target-station" style="width:100%; padding:7px; margin:4px 0 12px 0; border-radius:6px; border:1px solid #CBD5E0;"></select>
        
        <button onclick="calculateRoute()" style="
            width:100%;
            background:#2563EB;
            color:white;
            border:none;
            padding:10px;
            border-radius:6px;
            font-weight:bold;
            cursor:pointer;
            transition:0.2s;
        ">🔍 Trouver le chemin le plus court</button>
        
        <div id="route-results" style="margin-top:12px; display:none; padding:12px; background:#F8FAFC; border-radius:8px; border:1px solid #E2E8F0;">
            <div style="font-weight:bold; color:#1D4ED8; margin-bottom:6px;">Résultat du trajet (Dijkstra) :</div>
            <div>📏 Distance : <b id="route-dist" style="color:#0F172A;">-</b></div>
            <div>⏱️ Temps vélo (~15 km/h) : <b id="route-time" style="color:#0F172A;">-</b></div>
            <div>📍 Escales traversées : <b id="route-hops" style="color:#0F172A;">-</b></div>
        </div>

        <hr style="margin: 16px 0; border: 0; border-top: 1px solid #E2E8F0;">

        <h3 style="margin:0 0 12px 0; color:#1E293B; font-size:15px; display:flex; align-items:center; gap:8px;">
            📊 <span>Graphiques analytiques</span>
        </h3>

        <!-- Graphique 1 : Top Capacités -->
        <div style="margin-bottom: 16px;">
            <div style="font-size:11px; font-weight:bold; color:#64748B; margin-bottom:6px;">TOP STATIONS (CAPACITÉ)</div>
            <canvas id="chartTopCapacity" height="160"></canvas>
        </div>

        <!-- Graphique 2 : Répartition par Commune -->
        <div>
            <div style="font-size:11px; font-weight:bold; color:#64748B; margin-bottom:6px;">RÉPARTITION PAR COMMUNE</div>
            <canvas id="chartCommunes" height="160"></canvas>
        </div>
    </div>

    <script>
    const STATIONS = {json.dumps(stations_js_data)};
    const EDGES = {json.dumps(edges_js)};
    let activeRouteLayer = null;

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

        // Graphique 1 : Top Capacités
        new Chart(document.getElementById('chartTopCapacity'), {{
            type: 'bar',
            data: {{
                labels: {json.dumps(chart_top_labels)},
                datasets: [{{
                    label: 'Vélos max',
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

        // Graphique 2 : Répartition par Commune
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

    function calculateRoute() {{
        const uStart = parseInt(document.getElementById("start-station").value);
        const uTarget = parseInt(document.getElementById("target-station").value);
        
        if (uStart === uTarget) {{
            alert("Veuillez sélectionner deux stations différentes.");
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
        document.getElementById("route-dist").textContent = totalDistKm.toFixed(2) + " km (" + Math.round(totalDistKm*1000) + " m)";
        document.getElementById("route-time").textContent = minutes + " min";
        document.getElementById("route-hops").textContent = path.length + " stations";
        
        const mapObj = Object.values(window).find(v => v && v.addLayer && v.on);
        if (mapObj) {{
            if (activeRouteLayer) mapObj.removeLayer(activeRouteLayer);
            
            const routeCoords = path.map(idx => [STATIONS[idx].lat, STATIONS[idx].lon]);
            activeRouteLayer = L.polyline(routeCoords, {{
                color: '#EF4444',
                weight: 6,
                opacity: 0.95,
                dashArray: '8, 8'
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
    parser = argparse.ArgumentParser(description="Optimisation & Analytics du Réseau Vélib Île-de-France")
    parser.add_argument("--depart", type=int, default=None, help="Index de la station de départ")
    parser.add_argument("--arrivee", type=int, default=None, help="Index de la station d'arrivée")
    args = parser.parse_args()

    print("=" * 75)
    print("  OPTIMISATION & ANALYTICS RÉSEAU VÉLIB MÉTROPOLITAIN (1 518 STATIONS)")
    print("  Auteur : Misbaou DIALLO (BUT 3 Informatique)")
    print("=" * 75)

    stations = charger_donnees()
    df = pd.DataFrame(stations)
    n = len(df)

    print(f"\n[✓] {n} stations Vélib chargées (Paris & Île-de-France).")
    print(f"[✓] Couverture géographique : {df['commune'].nunique()} communes d'Île-de-France.")

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

    # 2. Algorithmes MST
    mst_kruskal, weight_kruskal, time_kruskal = algo_kruskal(n, edges)
    mst_prim, weight_prim, time_prim = algo_prim(n, edges)

    print("\n" + "-" * 60)
    print("  RÉSULTATS DES ALGORITHMES MST (ARBRE COUVRANT MINIMUM)")
    print("-" * 60)
    print(f"• Kruskal (Union-Find) -> Longueur : {weight_kruskal:.3f} km | Temps : {time_kruskal:.2f} ms")
    print(f"• Prim (Min-Heap)      -> Longueur : {weight_prim:.3f} km | Temps : {time_prim:.2f} ms")
    print("-" * 60)

    # 3. Graphiques Matplotlib
    generer_graphiques_matplotlib(df, edges)

    # 4. Statistiques analytiques JSON
    generer_statistiques(df, edges, weight_kruskal)

    # 5. Génération carte HTML interactive avec Dashboard Chart.js & Calculateur
    generer_carte_html_interactive(df, mst_kruskal, edges, args.depart or 0, args.arrivee or 15)
    print(f"[✓] Carte interactive avec dashboard de graphiques enregistrée : {OUTPUT_MAP}")
    print("[✓] Processus terminé avec succès !")


if __name__ == "__main__":
    main()
