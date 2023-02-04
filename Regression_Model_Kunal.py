#!/usr/bin/env python
# coding: utf-8

# # Importing Libraries

# In[73]:


import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, mean_absolute_percentage_error
from sklearn.metrics import accuracy_score


# # Reading CSV File

# In[74]:


df=pd.read_csv("Original_Dataset_csv.csv")
df.notna()
df.head()


# # RDKIT PROCEDURE (2D-Descriptor)

# In[75]:


logP_values = []
mol_weight = []
rot_bonds = []
aromatic_prop = []
tpsa = []
ringcount = []
hba = []
hbd = []
for i, smiles in enumerate(df['SMILES']):
    mol=Chem.MolFromSmiles(smiles)
    #logP Values
    logP=rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    logP_values.append(logP)
    #Molecular Weight
    weight = rdMolDescriptors.CalcExactMolWt(mol)
    mol_weight.append(weight)
    #No. of Rotate able Bonds
    n_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    rot_bonds.append(n_rot_bonds)
    #Aromatic Proposition
    prop = 1 - rdMolDescriptors.CalcFractionCSP3(mol)
    aromatic_prop.append(prop)
    #Topological Polar Surface
    tpsa_val = rdMolDescriptors.CalcTPSA(mol)
    tpsa.append(tpsa_val)
    #Ring Count of molecules
    ringcount_val = rdMolDescriptors.CalcNumRings(mol)
    ringcount.append(ringcount_val)
    #No. of Hydrogen Acceptor
    hba_val = rdMolDescriptors.CalcNumHBA(mol)
    hba.append(hba_val)
    #No. of Hydrogen Donor
    hbd_val = rdMolDescriptors.CalcNumHBD(mol)
    hbd.append(hbd_val)
    
df['logP']=logP_values
df['MolecularWeight'] = mol_weight
df['RotatableBonds'] = rot_bonds
df['AromaticProportion'] = aromatic_prop
df['TPSA'] = tpsa
df['RingCount'] = ringcount
df['HBA'] = hba
df['HBD'] = hbd


# In[76]:


#Edit the existing CSV value
df.to_csv('Original_Dataset_csv.csv', index=False)


# In[77]:


df


# # Random Forest Algorithm Implementation

# In[78]:


#Reading Train and Validation Dataset
train_data=pd.read_csv("Training_Dataset_csv.csv")
validate_data=pd.read_csv("Validation_Dataset_csv.csv")


# In[79]:


X_train = train_data[["logP", "MolecularWeight", "RotatableBonds", "AromaticProportion", "TPSA", "RingCount", "HBA", "HBD"]]
y_train = train_data["pIC50 (IC50 in microM)"]
X_valid= validate_data[["logP", "MolecularWeight", "RotatableBonds", "AromaticProportion", "TPSA", "RingCount", "HBA", "HBD"]]
y_valid = validate_data["pIC50 (IC50 in microM)"]


# In[80]:


#Random Forest for Validation Set
rf = RandomForestRegressor(n_estimators=90)
rf.fit(X_train, y_train)
y_pred_valid = rf.predict(X_valid)


# In[100]:


#Performance Checking using validation dataset"
mae = mean_absolute_error(y_valid, y_pred_valid)
mse = mean_squared_error(y_valid, y_pred_valid)
rmse = np.sqrt(mse)
r2 = r2_score(y_valid, y_pred_valid)
mape = mean_absolute_percentage_error(y_valid, y_pred_valid)
print("Mean Absolute Error            : ", mae)
print("Mean Squared Error             : ", mse)
print("R-squared scorer               :", r2)
print("Root mean squared Error        : ", rmse)
print("Mean absolute percentage error : ", mape)


# In[101]:


#Reading Test Dataset
test_data=pd.read_csv("Test_Dataset_csv.csv")
X_test= test_data[["logP", "MolecularWeight", "RotatableBonds", "AromaticProportion", "TPSA", "RingCount", "HBA", "HBD"]]
y_test = test_data["pIC50 (IC50 in microM)"]
#Random Forest for test dataset
rf = RandomForestRegressor(n_estimators=90)
rf.fit(X_train, y_train)
y_pred_test = rf.predict(X_test)


# In[102]:


y_pred_test


# Filling Values for Test_Dataset

# In[106]:


df_test=pd.read_csv("Test_Dataset_csv.csv")
df_test['pIC50 (IC50 in microM)']=y_pred_test


# In[107]:


df_test.head()


# In[108]:


#Save the test_dataset file
df_test.to_csv("C:\\Users\\HP\\Desktop\\Epidemio\\Epidemiology - Copy\\Test_Dataset_csv.csv", index=False)


# # Choosing Best Drug

# In[110]:


import urllib.request
import openbabel
import os
from openbabel import OBConversion, OBMol
url = "https://files.rcsb.org/download/6LU7.pdb"
urllib.request.urlretrieve(url, "6LU7.pdb")


# In[ ]:


conversion = OBConversion()
conversion.SetInAndOutFormats("smi", "pdb")

ligand = OBMol()
df1=pd.read_csv("Test_Dataset_csv.csv")
for i, smiles in enumerate(df['SMILES']):
    conversion.ReadString(ligand, smiles)
    
conversion.WriteFile(ligand, "ligand.pdb")


# In[ ]:


# Reading Protein and Ligand File
protein_file = "6LU7.pdb"
ligand_file = "ligand.pdb"
# Prepare the protein by removing the hydrogen atoms and writing to a new file
protein = openbabel.OBMol()
protein_conv = openbabel.OBConversion()
protein_conv.SetInFormat("pdb")
protein_conv.ReadFile(protein, protein_file)
protein.DeleteHydrogens()
protein_conv.SetOutFormat("pdb")
protein_conv.WriteFile(protein, "protein_noh.pdb")
# Prepare the ligand by adding hydrogens and writing to a new file
ligand = openbabel.OBMol()
ligand_conv = openbabel.OBConversion()
ligand_conv.SetInFormat("pdb")
ligand_conv.ReadFile(ligand, ligand_file)
ligand.AddHydrogens()
ligand_conv.SetOutFormat("pdb")
ligand_conv.WriteFile(ligand, "ligand_h.pdb")
# Running Autodock Vina using the prepared protein and ligand files
os.system("vina --config config.txt --ligand ligand_h.pdb --receptor protein_noh.pdb --out output.pdb")
#Binding Affinity predictions
binding_affinity = ad.get_binding_affinity()
print(binding_affinity)


# In[113]:


import openbabel

smiles = "CC(SC1=NC(C2=CC=CC=C2)=C(C#N)C(=O)N1)C(=O)NC1=CC=C(Cl)C=C1"
obConversion = openbabel.OBConversion()
obConversion.SetInFormat("smi")
obConversion.SetOutFormat("pdb")
mol = openbabel.OBMol()
obConversion.ReadString(mol, smiles)
pdb = obConversion.WriteString(mol)
print(pdb)


# In[ ]:




