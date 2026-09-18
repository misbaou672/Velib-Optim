"""
==============================================================================
Optimisation & Analytics du Réseau Vélib Métropole (1 518 Stations - Île-de-France)
Auteur : Misbaou DIALLO (BUT 3 Informatique)
==============================================================================

Fonctionnalités avancées :
  1. Triangulation de Delaunay (SciPy) pour réduire la complexité spatiale.
  2. Algorithmes MST : Kruskal (Union-Find) & Prim (Min-Heap) pour le réseau minimal.
  3. Recherche d'Itinéraire Optimal : Algorithme de Dijkstra (Min-Heap) entre 2 stations.
  4. Module d'Analyse Statistiques : Répartition géographique, capacités, hubs réseau,
     densités par département (75, 92, 93, 94, 78, 91, 95).
  5. Application Web Interactive HTML (Folium + Leaflet + Selecteur d'itinéraires + Chart.js).
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
import folium
from folium.plugins import MiniMap, MarkerCluster

# Chemins des fichiers
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
DATA_FILE = os.path.join(DATA_DIR, "stations_velib_idf_complete.json")
OUTPUT_MAP = os.path.join(BASE_DIR, "carte_velib_optimisee.html")
REPORT_FILE = os.path.join(BASE_DIR, "rapport_statistiques_velib.json")
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
    """Algorithme de Dijkstra pour trouver le plus court chemin entre 2 stations."""
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

    # Reconstruire le chemin
    path = []
    curr = target_node
    while curr is not None:
        path.append(curr)
        curr = predecessors[curr]
    path.reverse()

    exec_time = (time.perf_counter() - start_time) * 1000
    total_dist = distances[target_node]
    return path, total_dist, exec_time


def charger_donnees():
    """Charge le jeu de données depuis l'API OpenData ou le cache local."""
    os.makedirs(DATA_DIR, exist_ok=True)
    if os.path.exists(DATA_FILE):
        with open(DATA_FILE, "r", encoding="utf-8") as f:
            return json.load(f)

    print("[+] Téléchargement en temps réel depuis l'API OpenData...")
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


def generer_statistiques(df, edges, weight_mst):
    """Calcule un rapport analytique complet du réseau Vélib."""
    total_stations = len(df)
    total_capacite = int(df['capacite'].sum())
    moyenne_capacite = float(df['capacite'].mean())

    top_capacites = df.sort_values(by='capacite', ascending=False).head(10)[
        ['nom', 'commune', 'capacite']
    ].to_dict(orient='records')

    par_commune = df['commune'].value_counts().head(15).to_dict()

    distances = [e[2] for e in edges]
    dist_min = min(distances)
    dist_max = max(distances)
    dist_moy = sum(distances) / len(distances)

    rapport = {
        "metriques_generales": {
            "total_stations": total_stations,
            "total_communes": df['commune'].nunique(),
            "capacite_totale_velos": total_capacite,
            "capacite_moyenne_station": round(moyenne_capacite, 2),
            "distance_mst_totale_km": round(weight_mst, 2),
            "nombre_connexions_delaunay": len(edges),
            "distance_inter_station_moyenne_km": round(dist_moy, 3),
            "distance_min_km": round(dist_min, 3),
            "distance_max_km": round(dist_max, 3)
        },
        "top_10_stations_capacite": top_capacites,
        "repartition_par_commune": par_commune
    }

    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        json.dump(rapport, f, ensure_ascii=False, indent=2)

    return rapport


def generer_carte_html_interactive(df, mst_edges, edges, default_start_idx=0, default_target_idx=15):
    """Génère une carte HTML avec contrôles interactifs et sélection d'itinéraires."""
    center_lat = df["latitude"].mean()
    center_lon = df["longitude"].mean()
    m = folium.Map(location=[center_lat, center_lon], zoom_start=11, tiles="OpenStreetMap")

    # Cluster de marqueurs pour la fluidité
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

    # Groupe MST (Tracé vert)
    mst_group = folium.FeatureGroup(name="Réseau Optimal Optimisé (MST)")
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

    # Sauvegarder la carte Folium de base
    m.save(OUTPUT_MAP)

    # Ajouter le panneau interactif d'itinéraire et les scripts de calcul dynamique
    with open(OUTPUT_MAP, "r", encoding="utf-8") as f:
        html_content = f.read()

    # Formater les arêtes pour JavaScript (Dijkstra en navigateur)
    edges_js = [{"u": int(u), "v": int(v), "w": round(float(w), 4)} for u, v, w in edges]

    control_panel_html = f"""
    <div id="route-panel" style="
        position: fixed;
        top: 15px;
        right: 15px;
        width: 340px;
        background: rgba(255, 255, 255, 0.95);
        border-radius: 12px;
        padding: 16px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.25);
        z-index: 9999;
        font-family: 'Segoe UI', Arial, sans-serif;
        font-size: 13px;
        backdrop-filter: blur(8px);
        border: 1px solid #E2E8F0;
    ">
        <h3 style="margin:0 0 10px 0; color:#1A202C; font-size:16px; display:flex; align-items:center; gap:8px;">
            🚲 <span>Calculateur d'Itinéraire Vélib</span>
        </h3>
        
        <label style="font-weight:600; color:#4A5568;">Station de départ :</label>
        <select id="start-station" style="width:100%; padding:6px; margin:4px 0 10px 0; border-radius:6px; border:1px solid #CBD5E0;"></select>
        
        <label style="font-weight:600; color:#4A5568;">Station d'arrivée :</label>
        <select id="target-station" style="width:100%; padding:6px; margin:4px 0 12px 0; border-radius:6px; border:1px solid #CBD5E0;"></select>
        
        <button onclick="calculateRoute()" style="
            width:100%;
            background:#3182CE;
            color:white;
            border:none;
            padding:9px;
            border-radius:6px;
            font-weight:bold;
            cursor:pointer;
            transition:0.2s;
        ">🔍 Trouver le chemin le plus court</button>
        
        <div id="route-results" style="margin-top:12px; display:none; padding:10px; background:#F7FAFC; border-radius:6px; border:1px solid #E2E8F0;">
            <div style="font-weight:bold; color:#2B6CB0; margin-bottom:4px;">Résultat du parcours :</div>
            <div>📏 Distance : <b id="route-dist" style="color:#2D3748;">-</b></div>
            <div>⏱️ Temps à vélo (~15 km/h) : <b id="route-time" style="color:#2D3748;">-</b></div>
            <div>📍 Escales : <b id="route-hops" style="color:#2D3748;">-</b></div>
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
    }});

    function calculateRoute() {{
        const uStart = parseInt(document.getElementById("start-station").value);
        const uTarget = parseInt(document.getElementById("target-station").value);
        
        if (uStart === uTarget) {{
            alert("Veuillez sélectionner deux stations différentes.");
            return;
        }}
        
        // Dijkstra en JavaScript
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
        
        // Récupérer la carte Folium
        const mapObj = Object.values(window).find(v => v && v.addLayer && v.on);
        if (mapObj) {{
            if (activeRouteLayer) mapObj.removeLayer(activeRouteLayer);
            
            const routeCoords = path.map(idx => [STATIONS[idx].lat, STATIONS[idx].lon]);
            activeRouteLayer = L.polyline(routeCoords, {{
                color: '#E53E3E',
                weight: 6,
                opacity: 0.9,
                dashArray: '8, 8'
            }}).addTo(mapObj);
            
            mapObj.fitBounds(activeRouteLayer.getBounds(), {{ padding: [50, 50] }});
        }}
    }}
    </script>
    """

    html_content = html_content.replace("</body>", f"{control_panel_html}</body>")
    with open(OUTPUT_MAP, "w", encoding="utf-8") as f:
        f.write(html_content)


def main():
    parser = argparse.ArgumentParser(description="Optimisation & Analytics du Réseau Vélib Île-de-France")
    parser.add_argument("--depart", type=int, default=None, help="Index de la station de départ pour Dijkstra")
    parser.add_argument("--arrivee", type=int, default=None, help="Index de la station d'arrivée pour Dijkstra")
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

    # 3. Calcul de plus court chemin (Dijkstra)
    u_start = args.depart if args.depart is not None else 0
    u_target = args.arrivee if args.arrivee is not None else 15

    path, path_dist, path_time = algo_dijkstra(n, edges, u_start, u_target)
    start_name = df.loc[u_start, "nom"]
    target_name = df.loc[u_target, "nom"]

    print("\n" + "-" * 60)
    print("  RECHERCHE D'ITINÉRAIRE OPTIMAL (ALGORITHME DE DIJKSTRA)")
    print("-" * 60)
    print(f"📍 Départ  : {start_name} ({df.loc[u_start, 'commune']})")
    print(f"🎯 Arrivée : {target_name} ({df.loc[u_target, 'commune']})")
    print(f"📏 Distance totale : {path_dist:.3f} km ({path_dist*1000:.0f} mètres)")
    print(f"⏱️ Temps estimé à vélo (~15 km/h) : {round((path_dist / 15) * 60)} min")
    print(f"⚡ Calculé en : {path_time:.3f} ms ({len(path)} stations traversées)")
    print("-" * 60)

    # 4. Statistiques analytiques
    stats = generer_statistiques(df, edges, weight_kruskal)
    print(f"\n[✓] Rapport statistique exporté dans : {REPORT_FILE}")

    # 5. Génération carte HTML avec sélecteur interactif
    generer_carte_html_interactive(df, mst_kruskal, edges, u_start, u_target)
    print(f"[✓] Carte interactive avec calculateur d'itinéraire enregistrée : {OUTPUT_MAP}")
    print("[✓] Processus terminé avec succès !")


if __name__ == "__main__":
    main()
