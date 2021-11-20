# chemoinformatics
- RDKit - 강의 1 + Mol Mol2 포맷 | [youtube](https://www.youtube.com/watch?v=sxj56IQqhqM&list=PL30UV7ug7LwJYQgSp4THPjlb-9XAV4DCe&ab_channel=Prof.J.Lee)
# chemical file format
Ref : https://chem.libretexts.org/Courses/University_of_Arkansas_Little_Rock/ChemInformatics_(2017)%3A_Chem_4399_5399/2.2%3A_Chemical_Representations_on_Computer%3A_Part_II/2.2.2%3A_Anatomy_of_a_MOL_file  
- SMILES : 1D string
  Ref : https://www.daylight.com/meetings/summerschool98/course/dave/smiles-isomers.html  
  e.g.   C1(=C(C(=C(C(=C1C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4)C5=CC=CC=C5)C6=CC=CC=C6)C7=CC=CC=C7
- Mol file format  
  ![image](https://user-images.githubusercontent.com/48517782/142718373-2383e7db-8fdf-42df-8ada-8fec07e59e40.png)
  - hdeader block
  - atom block : the number of atoms and bonds, 3D coordinates of atoms
  - bond block : 
  - property block
- Mol2 format
  - atom types
  - partial charge
e.g.
```
 OpenBabel10082008492D

 25 26  0  0  1  0  0  0  0  0999 V2000
    0.1340   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000   -0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    2.0000   -0.0000 OpenBabel10082008492D

 25 26  0  0  1  0  0  0  0  0999 V2000
    0.1340   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000   -0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    2.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1340    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.8660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0000    1.7321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    2.5981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    1.7321    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -3.0000    2.7321    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0000    1.7321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0000    1.7321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5000    2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5000    2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.7321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1340    0.2321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1340   -0.7679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000   -1.2679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8660   -0.7679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8660    0.2321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  3  5  1  0  0  0  0
  3  7  1  0  0  0  0
  5  6  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  1  0  0  0  0
 10 11  2  0  0  0  0
 10 12  1  0  0  0  0
 12 13  1  0  0  0  0
 12 14  1  0  0  0  0
 12 20  1  0  0  0  0
 14 19  1  0  0  0  0
 14 15  2  0  0  0  0
 15 16  1  0  0  0  0
 16 17  2  0  0  0  0
 17 18  1  0  0  0  0
 18 19  2  0  0  0  0
 20 25  1  0  0  0  0
 20 21  1  0  0  0  0
 21 22  1  0  0  0  0
 22 23  1  0  0  0  0
 23 24  1  0  0  0  0
 24 25  1  0  0  0  0
M  CHG  1   3   1
M  END
$$$$    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1340    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.8660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0000    1.7321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    2.5981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    1.7321    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -3.0000    2.7321    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0000    1.7321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0000    1.7321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5000    2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5000    2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.7321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1340    0.2321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1340   -0.7679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000   -1.2679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8660   -0.7679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8660    0.2321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  3  5  1  0  0  0  0
  3  7  1  0  0  0  0
  5  6  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  1  0  0  0  0
 10 11  2  0  0  0  0
 10 12  1  0  0  0  0
 12 13  1  0  0  0  0
 12 14  1  0  0  0  0
 12 20  1  0  0  0  0
 14 19  1  0  0  0  0
 14 15  2  0  0  0  0
 15 16  1  0  0  0  0
 16 17  2  0  0  0  0
 17 18  1  0  0  0  0
 18 19  2  0  0  0  0
 20 25  1  0  0  0  0
 20 21  1  0  0  0  0
 21 22  1  0  0  0  0
 22 23  1  0  0  0  0
 23 24  1  0  0  0  0
 24 25  1  0  0  0  0
M  CHG  1   3   1
M  END
$$$$
```
# database
- ZINC database
- PubChem | [link](https://pubchem.ncbi.nlm.nih.gov/) 
  - largest database 
  e.g. aspirin - download `SDF` format
    - Although 3D position of a simple molecule can calculated with VSEPR, it is not accurate for complex atom. Atomic position in 3D requires the sophisticate calculation with quantum chemistry.
- ChEMBL database  
  - chmicals with biological reacticity