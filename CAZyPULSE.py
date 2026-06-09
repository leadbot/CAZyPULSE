# -*- coding: utf-8 -*-
"""
Created on Thu May 21 10:32:27 2026

@author: dl923 / leadbot
"""

import requests
from bs4 import BeautifulSoup
import re
import networkx as nx
from collections import Counter
import itertools
from pyvis.network import Network

#%% Basic network

URL = "https://www.cazy.org/PULDB/index.php?cazy_mod=GH20&sp_name2=Bacteroides+thetaiotaomicron&sp_ncbi2="
URL = "https://www.cazy.org/PULDB/index.php?cazy_mod=GH20&complementAgain=on&sp_name2=&sp_ncbi2="
FILTER_UNK = True  # true skips the ubiquitous 'unk' (unknown) proteins

print("Fetching and parsing PULDB data...")
headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}
response = requests.get(URL, headers=headers)

if response.status_code != 200:
    print(f"Failed to retrieve data. Status code: {response.status_code}")
    exit()

soup = BeautifulSoup(response.content, 'html.parser')
tables = soup.find_all('table')

puls_domains = []

for table in tables:
    for row in table.find_all('tr'):
        cols = row.find_all('td')
        if len(cols) >= 3:
            species_text = cols[0].get_text(strip=True)
            if "Bacteroides" in species_text and "thetaiotaomicron" in species_text:
                modularity_text = cols[2].get_text()
                cleaned_text = re.sub(r'[◀▶\s\xa0]+', ' ', modularity_text).strip()

                if not cleaned_text or cleaned_text == 'Modularity':
                    continue

                domains = []
                for item in cleaned_text.split():
                    for sub_item in item.split('|'):
                        sub_item = sub_item.strip()
                        if sub_item:
                            if FILTER_UNK and sub_item.lower() == 'unk':
                                continue
                            domains.append(sub_item)
                if domains:
                    puls_domains.append(domains)

print(f"Extracted {len(puls_domains)} PULs.")

# =============================================================================
# 2. NetworkX Graph 
# =============================================================================
node_counts = Counter()
edge_counts = Counter()

for domains in puls_domains:
    unique_domains = sorted(list(set(domains)))
    for d in unique_domains:
        node_counts[d] += 1
    for d1, d2 in itertools.combinations(unique_domains, 2):
        edge_counts[(d1, d2)] += 1

G = nx.Graph()

# Map metrics to Pyvis-recognized attributes ('size', 'value', 'title')
for node, count in node_counts.items():
    G.add_node(
        node, 
        size=max(10, count/2),                     # Controls node radius
        title=f"Domain: {node}<br>Occurrences: {count}", # Interactive hover tooltip
        label=node                                   # Visible text label
    )

for (u, v), weight in edge_counts.items():
    G.add_edge(
        u, v, 
        value=weight,                                # Controls edge thickness
        title=f"Co-occurrence Frequency: {weight}",   # Edge hover tooltip
        color="rgba(55, 55, 55, 1)"             # Subtle edge color
    )

# =============================================================================
# 3. Interactive Pyvis Visualization
# =============================================================================
print("Generating interactive HTML network...")

# Initialize Pyvis network object
# notebook=False generates a native local HTML file 
net = Network(
    height="800px", 
    width="100%", 
    bgcolor="#f8f9fa", 
    font_color="#212529", 
    notebook=False
)

# Populate Pyvis object using the NetworkX graph architecture
net.from_nx(G)

# Configure elegant node aesthetics
net.set_options("""
var options = {
  "nodes": {
    "borderWidth": 2,
    "color": {
      "border": "#4a5568",
      "background": "#63b3ed",
      "highlight": {
        "border": "#2b6cb0",
        "background": "#3182ce"
      }
    },
    "font": {
      "size": 14,
      "face": "monospace"
    }
  },
  "edges": {
    "smooth": {
      "type": "continuous"
    }
  },
  "physics": {
    "barnesHut": {
      "gravitationalConstant": -10000,
      "centralGravity": 0.3,
      "springLength": 150,
      "springConstant": 0.05
    },
    "minVelocity": 0.75
  }
}
""")


# Save and automatically open the interactive file
output_html = "pul_domain_interactive_network.html"
net.show(output_html, notebook=False)
print(f"Success! Open '{output_html}' in your web browser to view your network.")

#%% Scraping the DUFs

import os
import json
import re
import time
import requests
from bs4 import BeautifulSoup
import networkx as nx
import itertools
from collections import Counter
from pyvis.network import Network
import math 
# =============================================================================
# 1. Configuration & Native Cache Setup
# =============================================================================
BASE_URL = "https://www.cazy.org/PULDB/"
#INDEX_URL = "https://www.cazy.org/PULDB/index.php?cazy_mod=GH20&sp_name2=Bacteroides+thetaiotaomicron&sp_ncbi2="
INDEX_URL = "https://www.cazy.org/PULDB/index.php?cazy_mod=GH20&complementAgain=on&sp_name2=&sp_ncbi2="

HEADERS = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'}

CACHE_FILE = "pfam_local_cache.json"

def load_local_cache():
    if os.path.exists(CACHE_FILE):
        try:
            with open(CACHE_FILE, 'r') as f:
                return json.load(f)
        except Exception:
            return {}
    return {}

def save_local_cache(cache_data):
    try:
        with open(CACHE_FILE, 'w') as f:
            json.dump(cache_data, f, indent=4)
    except Exception as e:
        print(f"  [!] Failed to write cache file: {e}")

def get_pfam_name_from_ebi(pfam_id, cache):
    if pfam_id in cache:
        return cache[pfam_id]

    api_url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}/"
    try:
        time.sleep(0.3)
        response = requests.get(api_url, headers=HEADERS, timeout=10)
        if response.status_code == 200:
            data = response.json()
            name = data.get('metadata', {}).get('name', {}).get('name', '')
            if name:
                cache[pfam_id] = name
                save_local_cache(cache)
                return name
    except Exception:
        print("Failed")
        pass
    return "Unknown Structure"

pfam_cache = load_local_cache()

# =============================================================================
# 2. Scrape Index Page for First 10 PUL Links
# =============================================================================
print("Fetching main page and locating target PUL sub-pages...")
try:
    index_res = requests.get(INDEX_URL, headers=HEADERS, timeout=15)
    index_res.raise_for_status()
except Exception as e:
    print(f"Failed to access the main PULDB page: {e}")
    exit()

soup = BeautifulSoup(index_res.content, 'html.parser')

pul_links = []
for table in soup.find_all('table'):
    for row in table.find_all('tr'):
        cols = row.find_all('td')
        if len(cols) >= 3:
            species_text = cols[0].get_text(strip=True)
            
            # Skip the table header row
            if species_text == "Species" or not species_text:
                continue
                
            n0_col = cols[1]
            link_tag = n0_col.find('a')
            if link_tag and 'href' in link_tag.attrs:
                href = link_tag['href']
                full_link = href if href.startswith('http') else BASE_URL + href
                pul_name = link_tag.get_text(strip=True)
                
                # We now store a 3-item tuple: (PUL Name, URL, Species Name)
                if full_link not in [l[1] for l in pul_links]:
                    pul_links.append((pul_name, full_link, species_text))

target_puls = pul_links[:]
print(f"Building complete co-occurrence functional map across {len(target_puls)} PULs...\n")

node_counts = Counter()
edge_counts = Counter()
node_types = {}       # tracks if a node is a standard module or a resolved unk function
node_to_proteins = {} # maps node annotations to list of source protein IDs

# =============================================================================
# 3. Parse and Unify Standard Modules + Resolved unk Functions
# =============================================================================
c = 1
for name, url, species in target_puls:
    # Print the species directly to the console during processing
    print(f"Processing n={c} | {name} | Species: {species}")
    try:
        pul_res = requests.get(url, headers=HEADERS, timeout=15)
        pul_soup = BeautifulSoup(pul_res.content, 'html.parser')
        
        all_rows = pul_soup.find_all('tr')
        pul_functional_units = []
        
        for row in all_rows:
            cells_text = [cell.get_text(strip=True) for cell in row.find_all('td') if cell.get_text(strip=True)]
            
            # Filter 1: Check if row has enough columns
            if len(cells_text) < 2:
                continue
                
            protein_id = cells_text[0]
            
            # Filter 2: The Fix! Remove 'Btheta' and accept ANY valid locus tag 
            # (letters, numbers, underscores) while ignoring headers like "Species"
            if not re.match(r'^[A-Za-z0-9_\.]+$', protein_id):
                continue
                
            modularity_raw = cells_text[1]
            
            # Clean up the modularity string
            cleaned_mod = re.sub(r'(unk)', r' \1 ', modularity_raw, flags=re.IGNORECASE)
            cleaned_mod = re.sub(r'[◀▶\s\xa0|]+', ' ', cleaned_mod).strip()
            
            has_unk = 'unk' in modularity_raw.lower()
            
            # FORMATTING FIX: Inject the species string into the protein ID so it saves to the dictionary
            formatted_protein_entry = f"{protein_id} <i>[{species}]</i>"
            
            if has_unk:
                # --- SCENARIO A: Extract the true identity behind the unk ---
                whole_row_text = row.get_text()
                pfam_ids = re.findall(r'PF\d{5}', whole_row_text)
                
                pure_descriptions = []
                for item in cells_text:
                    if item == protein_id or re.match(r'PF\d{5}', item) or 'unk' in item.lower():
                        continue
                    pure_descriptions.append(item)
                
                if len(pure_descriptions) >= 2:
                    primary_annotation = pure_descriptions[1]
                elif len(pure_descriptions) == 1:
                    primary_annotation = pure_descriptions[0]
                else:
                    primary_annotation = "Uncharacterized Functional Domain"
                
                # Double check: hypothetical classifications with missing/malformed Pfams
                if "hypothetical" in primary_annotation.lower() and not pfam_ids:
                    primary_annotation = "Uncharacterized Functional Domain"
                
                # Check for hypothetical proteins to run Pfam resolving
                if "hypothetical" in primary_annotation.lower() and pfam_ids:
                    resolved_items = []
                    for pfid in pfam_ids:
                        resolved_desc = get_pfam_name_from_ebi(pfid, pfam_cache)
                        resolved_items.append(f"{pfid} ({resolved_desc})")
                    final_identity = ", ".join(resolved_items)
                else:
                    final_identity = primary_annotation
                
                if final_identity:
                    pul_functional_units.append(final_identity)
                    node_counts[final_identity] += 1
                    node_types[final_identity] = "resolved_unk"
                    
                    # Save the formatted entry (with species) instead of just the raw protein_id
                    if final_identity not in node_to_proteins:
                        node_to_proteins[final_identity] = []
                    node_to_proteins[final_identity].append(formatted_protein_entry)
                    
            else:
                # --- SCENARIO B: Capture standard functional modules ---
                for sub_item in cleaned_mod.split():
                    sub_item = sub_item.strip()
                    if sub_item:
                        pul_functional_units.append(sub_item)
                        node_counts[sub_item] += 1
                        node_types[sub_item] = "standard_module"
                        
                        # Save the formatted entry (with species)
                        if sub_item not in node_to_proteins:
                            node_to_proteins[sub_item] = []
                        node_to_proteins[sub_item].append(formatted_protein_entry)
                        
        # Form connections between all unique functional labels inside this single PUL cluster
        unique_units = sorted(list(set(pul_functional_units)))
        if len(unique_units) > 1:
            for unit1, unit2 in itertools.combinations(unique_units, 2):
                edge_counts[(unit1, unit2)] += 1
                
        # --- TRACKING METRICS ---
        unique_nodes = len(node_counts)
        total_appended_counts = sum(node_counts.values()) # Sums up every single increment step
        
        print(f"  -> Unique Nodes Discovered: {unique_nodes}")
        print(f"  -> Total Cumulative Appended Counts: {total_appended_counts}")
        print("-" * 60) # Visual separator for clean terminal scrolling
        
        time.sleep(0.4)
        c += 1
    except Exception as e:
        print(f"  [X] Error map parsing: {e}")
# =============================================================================
# 4. Build and Color the NetworkX Graph Architecture
# =============================================================================
G = nx.Graph()

for node, count in node_counts.items():
    # Assign distinct colors based on where the functional element originated
    if node_types.get(node) == "resolved_unk":
        node_color = "#f6ad55"  # Soft orange for the newly discovered unk functions
        border_color = "#dd6b20"
        group_name = "Resolved Unknown Feature"
    else:
        node_color = "#63b3ed"  # Clean blue for standard CAZyme/Sus modules
        border_color = "#3182ce"
        group_name = "Known Modularity Domain"
        
    G.add_node(
        node,
        label=node,
        size=max(12, count * 1.5),
        title=f"Functional Identity: {node}<br>Classification: {group_name}<br>Batched Group Appearances: {count}",
        color={"background": node_color, "border": border_color}
    )

for (u, v), weight in edge_counts.items():
    G.add_edge(
        u, v,
        value=weight,
        title=f"Co-occurrence frequency across PULs: {weight}",
        color="rgba(160, 174, 192, 0.5)"
    )

# =============================================================================
# 5. Initialize Pyvis Plot Layout
# =============================================================================
print("\nPlotting interactive integrated functional map...")

net = Network(
    height="800px",
    width="100%",
    bgcolor="#f7fafc",
    font_color="#1a202c",
    notebook=False
)

net.from_nx(G)

net.set_options("""
var options = {
  "nodes": {
    "borderWidth": 2,
    "font": {
      "size": 13,
      "face": "monospace"
    }
  },
  "edges": {
    "smooth": {
      "type": "continuous"
    }
  },
  "physics": {
    "barnesHut": {
      "gravitationalConstant": -12000,
      "centralGravity": 0.25,
      "springLength": 160,
      "springConstant": 0.04
    },
    "minVelocity": 0.6
  }
}
""")

output_html = "resolved_functional_cooccurrence_network_ALL_BACTERIA.html"
net.show(output_html, notebook=False)
print(f"\nDone! Open '{output_html}' in your browser. Orange nodes show your previously unknown annotations seamlessly linked to the blue CAZyme families.")
#%% Plot only
import math 
G = nx.Graph()

# Setup explicit boundaries
MIN_COUNT = 8
MAX_COUNT = 500
MIN_CONNECTIONS = 5 
TOP_N_EDGES = 10  # Retain only the N strongest connections per node

# New: Define keywords to filter OUT (case-insensitive)
IGNORE_KEYWORDS = ["SusC", "SusD"]

# 1. First pass: Filter nodes by raw abundance AND keyword blacklist
filtered_nodes = {}
dropped_by_keyword = 0

for node, count in node_counts.items():
    # Check if the node contains any of our blacklisted strings
    if any(keyword.upper() in node.upper() for keyword in IGNORE_KEYWORDS):
        dropped_by_keyword += 1
        continue
        
    # Apply original abundance limits
    if MIN_COUNT <= count <= MAX_COUNT:
        filtered_nodes[node] = count

# 2. Populate the graph with nodes (including custom green rule for GH families)
for node, count in filtered_nodes.items():
    is_duf = "DUF" in node.upper()
    is_gh = "GH" in node.upper() # Target identifier group
    
    if is_gh:
        node_color = "#48bb78"  # Soft Emerald Green
        border_color = "#38a169"
        highlight_color = "#2f855a"
        classification = "Glycoside Hydrolase (GH) Family"
    elif is_duf:
        node_color = "#f56565"  # Bright Red
        border_color = "#e53e3e"
        highlight_color = "#c53030"
        classification = "Pfam-Resolved DUF Domain"
    elif node_types.get(node) == "resolved_unk":
        node_color = "#f6ad55"  # Soft Orange
        border_color = "#dd6b20"
        highlight_color = "#c05621"
        classification = "Resolved Unknown Feature"
    else:
        node_color = "#63b3ed"  # Denim Blue
        border_color = "#3182ce"
        highlight_color = "#2b6cb0"
        classification = "Known CAZyme/Sus Domain"
        
    G.add_node(
        node,
        label=node,
        size=max(12, count / 2),
        title=f"Functional Identity: {node}<br>Classification: {classification}<br>Batched Group Appearances: {count}",
        color={
            "background": node_color, 
            "border": border_color,
            "highlight": {
                "background": highlight_color,
                "border": highlight_color
            }
        }
    )

# 3. Populate the graph with edges (Modified to pre-filter top N connections per node)
# First, collect all potential candidate edges that match our valid nodes
valid_edges = []
for (u, v), weight in edge_counts.items():
    if u in filtered_nodes and v in filtered_nodes:
        valid_edges.append((u, v, weight))

# Group connections by node to easily identify top weights per domain
node_to_edges = {}
for u, v, weight in valid_edges:
    node_to_edges.setdefault(u, []).append((u, v, weight))
    node_to_edges.setdefault(v, []).append((u, v, weight))

# Determine which unique edges belong to the top N of either connected node
allowed_edges = set()
for node, edges in node_to_edges.items():
    # Sort edges descending by weight
    sorted_edges = sorted(edges, key=lambda x: x[2], reverse=True)
    for u, v, weight in sorted_edges[:TOP_N_EDGES]:
        # Store canonical representation to handle undirected nature safely
        allowed_edges.add(tuple(sorted((u, v))))

# Actually write the top-performing backbone edges to the graph
for u, v, weight in valid_edges:
    if tuple(sorted((u, v))) in allowed_edges:
        # Log-scaling compresses massive values, math.log(weight + 1) multiplied by a scaling factor
        calculated_width = max(2, math.log(weight + 1) * 10)
        
        G.add_edge(
            u, v,
            value=weight,
            width=calculated_width,
            title=f"Co-occurrence frequency across PULs: {weight}",
            color={
                "color": "rgba(113, 128, 150, 0.45)",
                "highlight": "#ff0000",             
                "inherit": False                     
            }
        )

# 4. Iterative Pruning Loop (Crucial now that major hubs are gone!)
print("\nStabilizing network architecture...")
total_dropped_in_cascade = 0
pass_count = 1

while True:
    low_connection_nodes = [node for node, degree in G.degree() if degree < MIN_CONNECTIONS]
    
    if not low_connection_nodes:
        break
        
    total_dropped_in_cascade += len(low_connection_nodes)
    G.remove_nodes_from(low_connection_nodes)
    print(f"  Pass {pass_count}: Pruned {len(low_connection_nodes)} nodes after losing neighbors.")
    pass_count += 1

print(f"\n--- Network Pruning Strategy ---")
print(f" Dropped by keyword filter {IGNORE_KEYWORDS}: {dropped_by_keyword} nodes")
print(f" Retained after abundance limits: {len(filtered_nodes)} nodes")
print(f" Total dropped during iterative cascade pruning: {total_dropped_in_cascade} nodes")
print(f" Final operational nodes in layout: {G.number_of_nodes()}")
print(f"--------------------------------")


# =============================================================================
# 5. Initialize Pyvis Plot Layout
# =============================================================================
print("Plotting interactive integrated functional map (Threshold: count >= 2)...")

net = Network(
    height="1200px",
    width="100%",
    bgcolor="#f7fafc",
    font_color="#1a202c",
    notebook=False
)

net.from_nx(G)

# Fix: Added 'interaction' options to make selection focus cleanly on edges
net.set_options("""
var options = {
  "nodes": {
    "borderWidth": 2,
    "font": {
      "size": 13,
      "face": "monospace"
    }
  },
  "edges": {
    "smooth": {
      "type": "continuous"
    }
  },
  "interaction": {
    "hover": true,
    "selectConnectedEdges": true
  },
  "physics": {
    "forceAtlas2Based": {
      "gravitationalConstant": -80,
      "centralGravity": 0.01,
      "springLength": 250,
      "springConstant": 0.04,
      "damping": 0.8
    },
    "solver": "forceAtlas2Based",
    "minVelocity": 0.75,
    "stabilization": {
      "enabled": true,
      "iterations": 2000
    }
  }
}
""")

output_html = "resolved_functional_cooccurrence_network_ALL_BACTERIA.html"
net.show(output_html, notebook=False)
print(f"\nDone! Open '{output_html}' in your browser.")

#%% Rank correlations
# =============================================================================
# 6. Generate Co-occurrence Profile Heatmaps from Source Data
# =============================================================================
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform

print("\nAssembling robust co-occurrence DataFrame directly from raw counts...")

# 1. Gather all unique nodes that survived your initial abundance and keyword filters
active_nodes = sorted(list(filtered_nodes.keys()))
n_nodes = len(active_nodes)

if n_nodes < 2:
    print("Warning: Too few nodes match your filtering criteria to plot a heatmap.")
else:
    # 2. Pre-allocate an empty square array filled with zeros
    matrix_array = np.zeros((n_nodes, n_nodes))
    node_to_idx = {node: i for i, node in enumerate(active_nodes)}
    
    # Populate the Diagonal: Individual domain baseline group counts
    for node in active_nodes:
        idx = node_to_idx[node]
        matrix_array[idx, idx] = filtered_nodes[node]
        
    # Populate the Off-Diagonals: True raw co-occurrence frequencies from edge_counts
    for (u, v), weight in edge_counts.items():
        if u in node_to_idx and v in node_to_idx:
            idx_u = node_to_idx[u]
            idx_v = node_to_idx[v]
            matrix_array[idx_u, idx_v] = weight
            matrix_array[idx_v, idx_u] = weight
            
    # 3. Pack the populated data matrix into a beautifully indexed DataFrame
    df_profile = pd.DataFrame(matrix_array, index=active_nodes, columns=active_nodes)
    
    # Establish default layout rules
    sns.set_theme(style="white")
    plt.rcParams['font.family'] = 'monospace'

    # -------------------------------------------------------------------------
    # PLOT 1: Pearson Correlation Profile Heatmap
    # -------------------------------------------------------------------------
    print("Building Plot 1: Pearson Correlation Matrix...")
    df_pearson = df_profile.corr(method='pearson')
    
    # Fill any isolated or NaN columns with 0 to prevent linkage failures
    df_pearson = df_pearson.fillna(0)
    
    cg_pearson = sns.clustermap(
        df_pearson,
        cmap='viridis',          
        metric='euclidean',      
        method='ward',           
        figsize=(28, 28),
        cbar_kws={'label': 'Pearson Correlation Coefficient ($r$)'},
       # linewidths=0.1,
       # linecolor='#ffffff'      
    )
    
    plt.setp(cg_pearson.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8, ha='center')
    plt.setp(cg_pearson.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)
    
    cg_pearson.figure.suptitle(
        "Functional Association Profile Matrix\n(Hierarchical Clustering via Pearson $r$)",
        fontsize=14, fontdict={'weight': 'bold'}, y=1.02
    )
    plt.show()
    pearson_out = "functional_cooccurrence_pearson_heatmap.png"
    cg_pearson.savefig(pearson_out, dpi=300, bbox_inches='tight')
    print(f" Successfully exported full Pearson dataset plot to: '{pearson_out}'")
    plt.close('all')

    # -------------------------------------------------------------------------
    # PLOT 2: Euclidean Distance Matrix Heatmap
    # -------------------------------------------------------------------------
    print("\nBuilding Plot 2: Euclidean Distance Matrix...")
    
    # Evaluate geometric spacing directly from our source co-occurrence table
    eucl_distances = squareform(pdist(df_profile, metric='euclidean'))
    df_euclidean = pd.DataFrame(eucl_distances, index=active_nodes, columns=active_nodes)
    
    cg_euclidean = sns.clustermap(
        df_euclidean,
        cmap='viridis_r',        # Flipped color map so strong/close spatial hits glow yellow
        metric='euclidean',
        method='ward',
        figsize=(14, 14),
        cbar_kws={'label': 'Euclidean Distance ($d$)'},
       # linewidths=0.1,
       # linecolor='#ffffff'      
    )
    
    plt.setp(cg_euclidean.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8, ha='center')
    plt.setp(cg_euclidean.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)
    plt.show()

    cg_euclidean.figure.suptitle(
        "Functional Profile Dissimilarity Matrix\n(Hierarchical Clustering via Euclidean Distance)",
        fontsize=14, fontdict={'weight': 'bold'}, y=1.02
    )
    
    euclidean_out = "functional_cooccurrence_euclidean_heatmap.png"
    cg_euclidean.savefig(euclidean_out, dpi=300, bbox_inches='tight')
    print(f" Successfully exported full Euclidean dataset plot to: '{euclidean_out}'")
    plt.close('all')
#%% Test
import networkx as nx
from pyvis.network import Network

G = nx.Graph()

# Setup explicit boundaries
MIN_COUNT = 5
MAX_COUNT = 200
MIN_CONNECTIONS = 5 
TOP_N_EDGES = 3  # Set to 0 for ALL edges, or e.g., 3 for most abundant 3 edges per node

# Define keywords to filter OUT (case-insensitive)
IGNORE_KEYWORDS = ["SusC", "SusD"]

# 1. First pass: Filter nodes by raw abundance AND keyword blacklist
filtered_nodes = {}
dropped_by_keyword = 0

for node, count in node_counts.items():
    if any(keyword.upper() in node.upper() for keyword in IGNORE_KEYWORDS):
        dropped_by_keyword += 1
        continue
        
    if MIN_COUNT <= count <= MAX_COUNT:
        filtered_nodes[node] = count

# 2. Populate the graph with nodes (with custom green rule for GH families)
for node, count in filtered_nodes.items():
    is_duf = "DUF" in node.upper()
    is_gh = "GH" in node.upper() # New target identifier group
    
    if is_gh:
        node_color = "#48bb78"  # Soft Emerald Green
        border_color = "#38a169"
        highlight_color = "#2f855a"
        classification = "Glycoside Hydrolase (GH) Family"
    elif is_duf:
        node_color = "#f56565"  # Bright Red
        border_color = "#e53e3e"
        highlight_color = "#c53030"
        classification = "Pfam-Resolved DUF Domain"
    elif node_types.get(node) == "resolved_unk":
        node_color = "#f6ad55"  # Soft Orange
        border_color = "#dd6b20"
        highlight_color = "#c05621"
        classification = "Resolved Unknown Feature"
    else:
        node_color = "#63b3ed"  # Denim Blue
        border_color = "#3182ce"
        highlight_color = "#2b6cb0"
        classification = "Known CAZyme/Sus Domain"
        
    G.add_node(
        node,
        label=node,
        size=max(12, count * 1.5),
        title=f"Functional Identity: {node}<br>Classification: {classification}<br>Batched Group Appearances: {count}",
        color={
            "background": node_color, 
            "border": border_color,
            "highlight": {
                "background": highlight_color,
                "border": highlight_color
            }
        }
    )

# 3. Populate the graph with edges
for (u, v), weight in edge_counts.items():
    if u in filtered_nodes and v in filtered_nodes:
        G.add_edge(u, v, weight=weight)

# 4. Iterative Pruning Loop
print("\nStabilizing network architecture...")
total_dropped_in_cascade = 0
pass_count = 1

while True:
    low_connection_nodes = [node for node, degree in G.degree() if degree < MIN_CONNECTIONS]
    
    if not low_connection_nodes:
        break
        
    total_dropped_in_cascade += len(low_connection_nodes)
    G.remove_nodes_from(low_connection_nodes)
    print(f"  Pass {pass_count}: Pruned {len(low_connection_nodes)} nodes after losing neighbors.")
    pass_count += 1

print(f"\n--- Network Pruning Strategy ---")
print(f" Dropped by keyword filter {IGNORE_KEYWORDS}: {dropped_by_keyword} nodes")
print(f" Retained after abundance limits: {len(filtered_nodes)} nodes")
print(f" Total dropped during iterative cascade pruning: {total_dropped_in_cascade} nodes")
print(f" Final operational nodes in layout: {G.number_of_nodes()}")
print(f"--------------------------------")

# 4.5 Compute the Top-N connections per node profile
top_edges_global = set()
if TOP_N_EDGES > 0:
    for node in G.nodes():
        edges = G.edges(node, data=True)
        sorted_edges = sorted(edges, key=lambda x: x[2]['weight'], reverse=True)
        
        for u, v, data in sorted_edges[:TOP_N_EDGES]:
            top_edges_global.add(tuple(sorted((u, v))))

# Re-inject filtered network configurations with custom JavaScript parameters
final_edges = list(G.edges(data=True))
G.remove_edges_from(list(G.edges()))

for u, v, data in final_edges:
    weight = data['weight']
    # If TOP_N_EDGES is 0, every edge passes automatically
    is_top_path = True if TOP_N_EDGES == 0 else tuple(sorted((u, v))) in top_edges_global
    
    G.add_edge(
        u, v,
        value=weight,
        width=max(1, weight * 0.8),
        title=f"Co-occurrence frequency across PULs: {weight}",
        is_top_limit=is_top_path, # Attribute picked up by JS below
        color={
            "color": "rgba(160, 174, 192, 0.4)",
            "highlight": "#ff0000",             
            "inherit": False                     
        }
    )

# =============================================================================
# 5. Initialize Pyvis Plot Layout
# =============================================================================
print("Plotting interactive integrated functional map...")

net = Network(
    height="800px",
    width="100%",
    bgcolor="#f7fafc",
    font_color="#1a202c",
    notebook=False
)

net.from_nx(G)

net.set_options("""
var options = {
  "nodes": {
    "borderWidth": 2,
    "font": { "size": 13, "face": "monospace" }
  },
  "edges": {
    "smooth": { "type": "continuous" }
  },
  "interaction": {
    "hover": true,
    "selectConnectedEdges": true
  },
  "physics": {
    "forceAtlas2Based": {
      "gravitationalConstant": -120,
      "centralGravity": 0.15,
      "springLength": 220,
      "springConstant": 0.06,
      "damping": 0.7
    },
    "solver": "forceAtlas2Based",
    "minVelocity": 0.1,
    "maxVelocity": 30,
    "stabilization": {
      "enabled": true,
      "iterations": 2000
    }
  }
}
""")

output_html = "resolved_functional_cooccurrence_network_ALL_BACTERIA.html"
net.write_html(output_html)

# Read generated file to patch HTML layout interface mechanics
with open(output_html, 'r') as file:
    html_content = file.read()

# Build interactive toggle checkbox UI banner string
label_text = f"Filter Top {TOP_N_EDGES} Connections Only" if TOP_N_EDGES > 0 else "Displaying All Edges (0 Active Limit)"
checkbox_state = "checked" if TOP_N_EDGES > 0 else "disabled"

toggle_ui_html = f"""
<div id="filter-control-panel" style="
    position: absolute; top: 10px; left: 10px; z-index: 999;
    background: white; padding: 12px; border-radius: 6px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.15); font-family: monospace;">
    <label style="font-weight: bold; cursor: pointer; display: flex; align-items: center; gap: 8px;">
        <input type="checkbox" id="top-edge-toggle" {checkbox_state} style="transform: scale(1.2); cursor: pointer;">
        {label_text}
    </label>
</div>
"""

# JavaScript processing routine mapping actions to structural DOM mutations
filter_js_script = """
<script type="text/javascript">
document.getElementById('top-edge-toggle').addEventListener('change', function(e) {
    var filterActive = e.target.checked;
    var allEdges = edges.get();
    
    var updatedEdges = allEdges.map(function(edge) {
        // If edge fails top verification flag while toggle switch state reads active -> collapse path visibility
        if (filterActive && edge.is_top_limit === false) {
            edge.hidden = true;
        } else {
            edge.hidden = false;
        }
        return edge;
    });
    
    edges.update(updatedEdges);
    network.stabilize();
});

window.addEventListener('load', function() {
    setTimeout(function() {
        if(document.getElementById('top-edge-toggle').disabled === false) {
            document.getElementById('top-edge-toggle').dispatchEvent(new Event('change'));
        }
    }, 500);
});
</script>
"""
import os
import webbrowser

# Read generated file to patch HTML layout interface mechanics
with open(output_html, 'r') as file:
    html_content = file.read()

# Build interactive toggle checkbox UI banner string
label_text = f"Filter Top {TOP_N_EDGES} Connections Only" if TOP_N_EDGES > 0 else "Displaying All Edges (0 Active Limit)"
checkbox_state = "checked" if TOP_N_EDGES > 0 else "disabled"

toggle_ui_html = f"""
<div id="filter-control-panel" style="
    position: absolute; top: 10px; left: 10px; z-index: 999;
    background: white; padding: 12px; border-radius: 6px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.15); font-family: monospace;">
    <label style="font-weight: bold; cursor: pointer; display: flex; align-items: center; gap: 8px;">
        <input type="checkbox" id="top-edge-toggle" {checkbox_state} style="transform: scale(1.2); cursor: pointer;">
        {label_text}
    </label>
</div>
"""

# JavaScript processing routine mapping actions to structural DOM mutations
filter_js_script = """
<script type="text/javascript">
document.getElementById('top-edge-toggle').addEventListener('change', function(e) {
    var filterActive = e.target.checked;
    var allEdges = edges.get();
    
    var updatedEdges = allEdges.map(function(edge) {
        if (filterActive && edge.is_top_limit === false) {
            edge.hidden = true;
        } else {
            edge.hidden = false;
        }
        return edge;
    });
    
    edges.update(updatedEdges);
    network.stabilize();
});

window.addEventListener('load', function() {
    setTimeout(function() {
        if(document.getElementById('top-edge-toggle').disabled === false) {
            document.getElementById('top-edge-toggle').dispatchEvent(new Event('change'));
        }
    }, 500);
});
</script>
"""

# Splice injections safely across body structures
html_content = html_content.replace('<body>', '<body>\n' + toggle_ui_html)
html_content = html_content.replace('</body>', filter_js_script + '\n</body>')

# Save the finalized, fully customized dashboard
with open(output_html, 'w') as file:
    file.write(html_content)

# The Fix: Use Python's built-in webbrowser to mimic net.show()'s opening trigger safely
webbrowser.open('file://' + os.path.realpath(output_html))
print(f"\nDone! Successfully launched '{output_html}' in your default web browser.")
#%% Interrogate nodes
# =============================================================================
# Print out the results specifically for your target node
# =============================================================================
target_node = "Uncharacterized Functional Domain"

if target_node in node_to_proteins:
    proteins = node_to_proteins[target_node]
    print(f"=== Proteins clustered inside '{target_node}' ({len(proteins)} total) ===")
    for pid in sorted(list(set(proteins))):  # set() deduplicates multiple domain occurrences on one protein
        print(f" - {pid}")
else:
    print(f"No proteins found for '{target_node}'. Check spelling or threshold conditions.")
    
#%% Save data
# Create the structural payload from your existing workspace variables
export_payload = {
    "node_counts": dict(node_counts),
    "node_types": node_types,
    "node_to_proteins": node_to_proteins,
    "edge_counts": {f"{u}|||{v}": weight for (u, v), weight in edge_counts.items()}
}

# Dump cleanly to a static local database file
with open("pul_network_data.json", "w") as f:
    json.dump(export_payload, f, indent=4)

print("Workspace dataset successfully exported to 'pul_network_data.json'.")

#%% LOAD and RUN
import json
import networkx as nx
from pyvis.network import Network

# Configurable parameter for your session
MIN_OCCURRENCE = 2  # Strip out nodes that occur less than twice

# 1. Load data from local storage
with open("pul_network_data.json", "r") as f:
    data = json.load(f)

node_counts = data["node_counts"]
node_types = data["node_types"]
node_to_proteins = data["node_to_proteins"]
edge_counts_raw = data["edge_counts"]

# 2. Filter nodes based on your threshold rule
filtered_nodes = {node: count for node, count in node_counts.items() if count >= MIN_OCCURRENCE}
print(f"Loaded {len(filtered_nodes)} nodes matching threshold rules (Count >= {MIN_OCCURRENCE}).")

# 3. Rebuild the NetworkX Graph Topology
G = nx.Graph()

for node, count in filtered_nodes.items():
    is_duf = "DUF" in node.upper()
    
    # Apply your visual coloring rule engine
    if is_duf:
        node_color, border_color, highlight_color = "#f56565", "#e53e3e", "#c53030"
        classification = "Pfam-Resolved DUF Domain"
    elif node_types.get(node) == "resolved_unk":
        node_color, border_color, highlight_color = "#f6ad55", "#dd6b20", "#c05621"
        classification = "Resolved Unknown Feature"
    else:
        node_color, border_color, highlight_color = "#63b3ed", "#3182ce", "#2b6cb0"
        classification = "Known CAZyme/Sus Domain"
    
    # Build out your custom interactive hover list
    associated_proteins = node_to_proteins.get(node, [])
    proteins_list_html = "".join([f"<li>{pid}</li>" for pid in sorted(list(set(associated_proteins)))])
    
    tooltip_content = (
        f"<b>Functional Identity:</b> {node}<br>"
        f"<b>Classification:</b> {classification}<br>"
        f"<b>Total Count:</b> {count}<br>"
        f"<b>Associated Loci ({len(associated_proteins)}):</b><br>"
        f"<ul style='margin:4px; padding-left:16px;'>{proteins_list_html}</ul>"
    )
        
    G.add_node(
        node,
        label=node,
        size=max(12, count * 1.5),
        title=tooltip_content,
        color={
            "background": node_color, 
            "border": border_color,
            "highlight": {"background": highlight_color, "border": highlight_color}
        }
    )

# Reconstruct edges safely
for edge_string, weight in edge_counts_raw.items():
    u, v = edge_string.split("|||")
    if u in filtered_nodes and v in filtered_nodes:
        G.add_edge(u, v, value=weight, title=f"Co-occurrence frequency: {weight}", color="rgba(160, 174, 192, 0.5)")

# 4. Spin up the Visual Interactive Viewport
net = Network(height="800px", width="100%", bgcolor="#f7fafc", font_color="#1a202c", notebook=False)
net.from_nx(G)
net.set_options("""
var options = {
  "nodes": { "borderWidth": 2, "font": { "size": 13, "face": "monospace" } },
  "edges": { "smooth": { "type": "continuous" } },
  "physics": {
    "barnesHut": { "gravitationalConstant": -15000, "centralGravity": 0.3, "springLength": 180, "springConstant": 0.04 },
    "minVelocity": 0.6
  }
}
""")

net.show("resolved_functional_cooccurrence_network.html", notebook=False)
print("Network fully loaded and generated without any remote web queries.")
