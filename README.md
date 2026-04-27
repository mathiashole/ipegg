# ipegg

[![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=R&logoColor=white&labelColor=101010)](https://www.r-project.org/about.html)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/mathiashole/ipegg?style=for-the-badge&labelColor=101010&color=white)
![GitHub last commit](https://img.shields.io/github/last-commit/mathiashole/ipegg?style=for-the-badge&labelColor=101010)

`ipegg` (**I**nter**P**roscan **E**nvironment & **G**enomic **G**rapher) is an R-based tool designed to visualize protein domain architectures from InterProScan TSV outputs. It creates publication-quality plots that respect the full protein length and allows for deep customization via YAML configuration files.

## 🚀 Features

### 📊 Advanced Visualization
- **True Backbone**: Represents the protein from amino acid 0 to its full length (V3 column), highlighting N- and C-terminal regions without domains.
- **Intelligent Scaling**: Automatically adjusts image height based on the number of sequences to maintain consistent block thickness.
- **Natural Sorting**: Handles alphanumeric IDs (e.g., Chr1, Chr2, Chr10) correctly.

### 🎨 Customization
- **Regex Normalization**: Automatically clean or rename complex domain names.
- **Domain Filtering**: Remove noisy or uninformative hits (e.g., non-cytoplasmic regions).
