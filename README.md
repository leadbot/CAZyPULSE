# CAZyPULSE
<div align="center">
  <img src="https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png" width="100%" />
  
  <h1>🧬 CAZyPULSE</h1>
  <p><b>CAZy PUL Scraper and Evaluation Network</b></p>
  A Polysaccharide Utilization Loci Modularity & Association Mapper<br><p></p>
## Pipeline At A Glance

  <p>
    <a href="https://github.com/YOUR_GITHUB_USERNAME/PUL-Mapper/stargazers"><img src="https://img.shields.io/github/stars/YOUR_GITHUB_USERNAME/PUL-Mapper?color=blue&style=flat-square" alt="Stars"></a>
    <a href="https://github.com/YOUR_GITHUB_USERNAME/PUL-Mapper/network/members"><img src="https://img.shields.io/github/forks/YOUR_GITHUB_USERNAME/PUL-Mapper?color=violet&style=flat-square" alt="Forks"></a>
    <a href="https://github.com/YOUR_GITHUB_USERNAME/PUL-Mapper/blob/main/LICENSE"><img src="https://img.shields.io/github/license/YOUR_GITHUB_USERNAME/PUL-Mapper?color=green&style=flat-square" alt="License"></a>
    <img src="https://img.shields.io/badge/Python-3.8+-3776AB?style=flat-square&logo=python&logoColor=white" alt="Python Version">
  </p>
  
  <h4>Transforming raw CAZy PULDB modularity strings into stable, interactive functional networks.</h4>
</div>

<br />

## Pipeline At A Glance
**Input**<br>
CAZy PUL link e.g. https://www.cazy.org/PULDB/index.php?cazy_mod=GH1&sp_name2=&sp_ncbi2=
<table width="100%">
  <tr>
    <td width="33%" align="center"><b>1. Scrape & Extract</b></td>
    <td width="33%" align="center"><b>2. Resolve Unidentified Pfams</b></td>
    <td width="33%" align="center"><b>3. Network Optimization</b></td>
  </tr>
  <tr>
    <td>
      <ul>
        <li>Targeted crawling of CAZy PULDB indexes.</li>
        <li>Isolates clean genomic locus tags.</li>
        <li>Filters out ubiquitous background noise.</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Detects hidden <code>unk</code> functional units.</li>
        <li>Live EBI InterPro API integration.</li>
        <li>Native JSON caching layer to bypass throttling.</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Weighted Co-occurrence mapping via NetworkX.</li>
        <li>Iterative multi-pass cascade pruning.</li>
        <li>Interactive ForceAtlas2 Pyvis rendering.</li>
      </ul>
    </td>
  </tr>
</table>

<br />

## Interactive Visualizations & Charts

<details>
<summary><b>📈 Click to expand Output Previews</b></summary>
<br />

### 1. Interactive Force-Directed Network (`.html`)
The generated network cleanly separates classifications using customized color weights:
* 🟢 **Green Nodes:** Glycoside Hydrolase (GH) Families
* 🔴 **Red Nodes:** Pfam-Resolved DUF Domains
* 🟠 **Orange Nodes:** Resolved Unknown Genomic Features
* 🔵 **Blue Nodes:** Known CAZyme / Sus Structural Domains

### 2. Statistical Coregulation Profiles (`.png`)
The script processes matrix dependencies directly through Seaborn to construct two major clustering variants:

<table>
  <tr>
    <td align="center"><b>Pearson Correlation Profile</b></td>
    <td align="center"><b>Euclidean Distance Heatmap</b></td>
  </tr>
  <tr>
    <td><img src="functional_cooccurrence_pearson_heatmap.png" width="100%" alt="Pearson Heatmap Placeholder" /></td>
    <td><img src="functional_cooccurrence_euclidean_heatmap.png" width="100%" alt="Euclidean Heatmap Placeholder" /></td>
  </tr>
</table>

</details>

<br />

<img src="https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png" width="100%" />
