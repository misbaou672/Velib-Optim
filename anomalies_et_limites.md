# Rapport d'Anomalies, Limites et Points d'Attention

**Projet** : Optimisation et Analyse Spatiale du Réseau Vélib Île-de-France  
**Auteur** : Misbaou DIALLO (BUT 3 Informatique)  
**Date** : Septembre 2026  

---

## Points de Vigilance & Anomalies Identifiées

### 1. Connexion API OpenData Paris & Résilience Réseau
- **Comportement** : L'API OpenData Paris (`opendata.paris.fr`) peut occasionnellement subir des latences ou des limites de débit (*rate limiting*) en cas de requêtes trop fréquentes.
- **Impact** : En cas d'indisponibilité du réseau ou de réponse lente (>15s), l'application bascule automatiquement sur le cache local (`stations_velib_idf_complete.json`).
- **Piste d'amélioration** : Mettre en place un mécanisme de *retry* exponentiel et un cache navigateur (LocalStorage/IndexedDB).

### 2. Rendu des Graphiques Chart.js dans le Tiroir Masqué
- **Comportement** : Lorsque Chart.js est initialisé dans une div avec `display: none`, le canvas calcule une taille nulle (`0px × 0px`).
- **Solution appliquée** : Le script déclenche l'initialisation et le redimensionnement (`initCharts()`) avec un léger délai de 50 ms lors du clic sur le bouton **Stats**.
- **Point à surveiller** : Sur certains navigateurs mobiles anciens, un second clic sur **Stats** peut être nécessaire pour forcer le recalcul du canvas.

### 3. Limites Géométriques de la Triangulation de Delaunay
- **Comportement** : L'enveloppe convexe (*Convex Hull*) de Delaunay relie parfois des stations périphériques éloignées au-dessus des forêts ou des boucles de la Seine.
- **Impact** : Certaines arêtes candidates du maillage dépassent 5 km de distance inter-stations.
- **Correction apportée** : L'algorithme de Kruskal / Prim élimine automatiquement les arêtes trop longues pour ne conserver que l'Arbre Couvrant Minimum (MST) de **502,09 km**.

### 4. Dépendance de Test Navigateur Headless (Google Chrome)
- **Comportement** : L'outil de test automatisé par sous-agent sous environnement Linux conteneurisé renvoie `google-chrome not found in PATH`.
- **Impact** : Impossible d'effectuer des captures d'écran automatisées sans installer Chrome headless dans le conteneur.
- **Solution** : Les validations unitaires de l'HTML et du JavaScript sont exécutées via le moteur de validation Python et testées en direct sur le serveur local `http://127.0.0.1:8080/carte_velib_optimisee.html`.
