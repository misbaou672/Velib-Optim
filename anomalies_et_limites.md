# Rapport Technique d'Anomalies, Limites & Architecture

![Architecture](https://img.shields.io/badge/Architecture-Optimisée-6366F1?style=for-the-badge&logo=cpu&logoColor=white)
![Qualité](https://img.shields.io/badge/Qualité-Production-10B981?style=for-the-badge&logo=check-circle&logoColor=white)
![Projet](https://img.shields.io/badge/Projet-Universitaire-38BDF8?style=for-the-badge&logo=graduation-cap&logoColor=white)

---

> [!IMPORTANT]
> **Contexte Technique** : Document de synthèse d'architecture rédigé par **Misbaou DIALLO** dans le cadre de l'extension et de l'amélioration continue sur mon temps libre du projet universitaire initial de **BUT Informatique**.

---

## Synthèse des Points de Vigilance & Solutions Appliquées

### 1. Résilience de la Connexion API OpenData Paris
- **Constat** : L'API OpenData Paris peut parfois présenter des latences ou des erreurs réseau temporaires.
- **Solution mise en place** : Mecanisme d'interrogation synchrone avec bascule automatique vers le cache local `stations_velib_idf_complete.json` en cas d'échec de la requête.

### 2. Rendu des Canvas Chart.js en Container Masqué
- **Constat** : Les éléments `<canvas>` de Chart.js initialisés dans un conteneur masqué (`display: none`) calculent une taille nulle (`0px × 0px`).
- **Solution mise en place** : Initialisation différée `initCharts()` et appel à `.resize()` déclenchés **50 ms après** l'ouverture effective du tiroir d'analytics.

### 3. Ordre d'Injection des Scripts JavaScript (Folium)
- **Constat** : L'injection de code JS personnalisé avant les déclarations des variables de calques Folium (`var feature_group_...`) provoquait un arrêt de l'exécution par `ReferenceError`.
- **Solution mise en place** : Placement du bloc UI personnalisé à la toute fin du document HTML (`</html>`) et référencement sécurisé par nom de chaîne littérale.

---

**Auteur** : **Misbaou DIALLO** *(BUT 3 Informatique)*
