# chemoinformatics
- RDKit - 강의 1 + Mol Mol2 포맷 | [youtube](https://www.youtube.com/watch?v=sxj56IQqhqM&list=PL30UV7ug7LwJYQgSp4THPjlb-9XAV4DCe&ab_channel=Prof.J.Lee)
# chemical file format
Ref : https://chem.libretexts.org/Courses/University_of_Arkansas_Little_Rock/ChemInformatics_(2017)%3A_Chem_4399_5399/2.2%3A_Chemical_Representations_on_Computer%3A_Part_II/2.2.2%3A_Anatomy_of_a_MOL_file  
- SMILES : 1D string
  e.g. C1(=C(C(=C(C(=C1C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4)C5=CC=CC=C5)C6=CC=CC=C6)C7=CC=CC=C7
- Mol file format  
  ![image](https://user-images.githubusercontent.com/48517782/142718373-2383e7db-8fdf-42df-8ada-8fec07e59e40.png)
  - hdeader block
  - atom block : the number of atoms and bonds, 3D coordinates of atoms
  - bond block : 
  - property block
- Mol2 format
  - atom types
  - partial charge
# database
- ZINC database
- PubChem | [link](https://pubchem.ncbi.nlm.nih.gov/) 
  - largest database 
  e.g. aspirin - download `SDF` format
    - Although 3D position of a simple molecule can calculated with VSEPR, it is not accurate for complex atom. Atomic position in 3D requires the sophisticate calculation with quantum chemistry.
- ChEMBL database  
  - chmicals with biological reacticity
