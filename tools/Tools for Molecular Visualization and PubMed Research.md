## Overview

This repository contains two Python tools designed to assist with molecular visualization and literature research in the field of chemistry and drug development. The tools are:

1. **`2D_images.py`**: Generates 2D molecular structure images based on SMILES molecular formulas extracted from a CSV file.
2. **`PubMed.py`**: Retrieves research articles from PubMed, including titles, abstracts, and publication dates, based on a specified SMILES string or drug name.

---

## 1. 2D_images.py

### Purpose
The `2D_images.py` script processes a CSV file containing molecular data, specifically a column named `mol` with SMILES strings, and generates 2D molecular structure images for all SMILES strings found in the specified directory.

### Features
- Reads a CSV file with molecular data.
- Parses the `mol` column to extract SMILES strings.
- Generates high-quality 2D images for each SMILES formula.
- Saves the images in the specified output directory.

### Usage
1. Prepare a CSV file with a column named `mol` containing SMILES strings.
2. Run the script, specifying the input CSV file and output directory.
3. The script will generate and save 2D molecular structure images.


### Dependencies
- `rdkit` for molecular structure visualization.
- `pandas` for CSV handling.

---

## 2. PubMed.py

### Purpose
The `PubMed.py` script retrieves research articles from PubMed related to a specified SMILES string or drug name. It is designed to provide quick access to relevant scientific literature.

### Features
- Accepts a SMILES string or drug name as input.
- Queries PubMed for related articles.
- Retrieves the titles, abstracts, and publication dates of up to five articles.
- Outputs the information in a structured format.

### Usage
1. Specify a SMILES string or drug name as the input query.
2. Run the script to fetch the related articles.
3. Review the titles, abstracts, and publication dates of the retrieved articles.

### Dependencies
- `pandas` for data manipulation.
- `requests` for making API calls to PubMed.

---
