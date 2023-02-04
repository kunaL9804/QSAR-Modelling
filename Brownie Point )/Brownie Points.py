#!/usr/bin/env python
# coding: utf-8

# # Brownie Points :)

# Importing Libraries

# In[33]:


import pandas as pd
import numpy as np
import csv
import random
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, mean_absolute_percentage_error
from sklearn.metrics import accuracy_score


# In[34]:


df_train=pd.read_csv("C:\\Users\\HP\\Desktop\\Epidemio\\Epidemiology - Copy\\Training_Dataset_csv.csv")
df_test=pd.read_csv("C:\\Users\\HP\\Desktop\\Epidemio\\Epidemiology - Copy\\Test_Dataset_csv.csv")


# In[35]:


df_train


# In[36]:


df_test


# In[37]:


#Merging Data
combined_df = pd.concat([df_train, df_test], axis=0)
combined_df.to_csv('merged_file.csv', index=False)


# In[38]:


combined_df


# In[39]:


# Load the original CSV file
df = pd.read_csv('merged_file.csv')

# Shuffle the rows of the dataframe
df = df.sample(frac=1)

# Split the dataframe into two new dataframes
split_index = int(df.shape[0] * 0.8) # 80% split
df1, df2 = df.iloc[:split_index], df.iloc[split_index:]

# Save the two new dataframes as CSV files
df1.to_csv('train_file.csv', index=False)
df2.to_csv('test_file.csv', index=False)


# # Internal Validation

# In[40]:


train_data=pd.read_csv("train_file.csv")
test_data=pd.read_csv("test_file.csv")


# In[41]:


X_train = train_data[["logP", "MolecularWeight", "RotatableBonds", "AromaticProportion", "TPSA", "RingCount", "HBA", "HBD"]]
y_train = train_data["pIC50 (IC50 in microM)"]
X_test= test_data[["logP", "MolecularWeight", "RotatableBonds", "AromaticProportion", "TPSA", "RingCount", "HBA", "HBD"]]
y_test = test_data["pIC50 (IC50 in microM)"]


# In[42]:


#Random Forest for Validation
rf = RandomForestRegressor(n_estimators=90)
rf.fit(X_train, y_train)
y_pred = rf.predict(X_test)


# In[43]:


#Performance Checking using validation dataset"
mae = mean_absolute_error(y_test, y_pred)
mse = mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)
r2 = r2_score(y_test, y_pred)
mape = mean_absolute_percentage_error(y_test, y_pred)
print("Mean Absolute Error            : ", mae)
print("Mean Squared Error             : ", mse)
print("R-squared scorer               :", r2)
print("Root mean squared Error        : ", rmse)
print("Mean absolute percentage error : ", mape)


# In[32]:


from sklearn.model_selection import LeaveOneOut, cross_val_score
from sklearn.metrics import make_scorer

# define the LOO-CV method
loo = LeaveOneOut()

# define the scoring metric, in this case R^2 
scoring = make_scorer(r2_score)

# calculate the LOO-CV score
loo_q2_score = cross_val_score(rf, X_train, y_train, cv=loo, scoring=scoring)
print(loo_q2_score)


# In[ ]:





# In[ ]:





# In[ ]:




