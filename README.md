# M<sup>3</sup>-20M: Multi-Modal Molecular Dataset

Welcome to the official repository for M<sup>3</sup>-20M, the first large-scale Multi-Modal Molecular dataset, containing over 20 million molecules! 🎉

## Overview
M<sup>3</sup>-20M (Multi-Modal Molecular dataset) is designed to support AI-driven drug design and discovery. It provides an unprecedented scale that highly benefits the training or fine-tuning of large models for superior performance in drug design and discovery tasks.

## Features
- **Scale**: Contains over 20 million molecules, 71 times more than the largest existing dataset.
- **Comprehensive Modalities**:
- One-dimensional SMILES strings
- Two-dimensional molecular graphs
- Three-dimensional molecular structures
- Physicochemical properties
- Text descriptions
- **Diverse Applications**: Supports various downstream tasks such as molecule generation, molecular property prediction, lead optimization, virtual screening, pharmacokinetics modeling, and drug-target interaction prediction.

## Dataset Details
M<sup>3</sup>-20M integrates data from multiple sources to provide a comprehensive view of each molecule. Here’s what you can find in the dataset:
- M^3_Original.csv: Descriptions from PubChem
- M^3_Physicochemical.csv: Physicochemical properties
- M^3_Description_Physicochemical.csv: Descriptions composed of physicochemical properties
- M^3_Multi.csv: Descriptions from PubChem, physicochemical properties, and those generated by GPT-3.5
- MPP folder: Contains multimodal datasets for molecular property prediction (BBBP-MM, BACE-MM, HIV-MM, ClinTox-MM, Tox21-MM)
- MOSES-Multi folder: Contains MOSES multimodal datasets for molecular generation
- QM9-Multi folder: Contains QM9 multimodal datasets

## Functions
We provide convenient functions that allow you to easily obtain the dataset, as well as the 2D and 3D representations of any molecule outside the dataset. The specific functions can be found in the Function folder.

## Download Links
The dataset is available for download from multiple sources:

- **Google Drive**: [Download Link](https://drive.google.com/drive/folders/1ai_HfcWfWoRdsfsDscfR1Dlg2OGEemWi?usp=sharing)
- **Baidu Cloud**:  [Download Link](https://pan.baidu.com/s/1kNL32Rj3r9PgdvMWQSDVhA?pwd=ADMS) password：ADMS
- **Hugging Face**:  [Download Link](https://huggingface.co/datasets/Alex99Gsy/M-3_Multi-Modal-Molecule)




## Example Usage
Here’s a simple example of how to load and explore the dataset:

```
import pandas as pd

# Load the dataset
df = pd.read_csv('path-to-dataset.csv')

# Display the first few rows
print(df.head())
```

## Contributing
We welcome contributions from the community! Feel free to submit issues or pull requests to help improve the dataset and its applications.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements
This dataset is a collaborative effort by researchers from Tongji University and Fudan University. We thank Siyuan Guo, Lexuan Wang, Chang Jin, Jinxian Wang, Han Peng, Huayang Shi, Wengen Li, Jihong Guan, and Shuigeng Zhou for their contributions and support.
## Contact
For any questions or inquiries, please reach out to gsy9901224@tongji.edu.cn.

Enjoy using M<sup>3</sup>-20M and happy researching! 🚀🔬
