#!/usr/bin/env python
# coding: utf-8

# In[1]:

## Import necessary packages
import math
import pandas as pd 

## Read import file
file = open('chunk3_set3_modifiedmols_csearch.txt')
dat = file.readlines()


# In[2]:

## Create empty list to store descriptor names
descNames = []

## Add contents of data header line to descriptor names list
headerLine = dat[0].split(',')
for colName in range(len(headerLine)):
    descNames.append(headerLine[colName])
## Specify index of first MOE descriptor in names list
firstDesc = descNames.index('apol')


# In[3]:

## Create a dictionary for each molecule

## Dictionary assigns each molecule a database ID and saves an entry for each
## conformer containing corresponding descriptor values

d = {}
for line in dat[1:]:
    line = line.split(',')
    if '' in line:
        continue
    else:
        databaseID = line[11]
        dE = float(line[3])
        descriptorVals = [dE]
        for col in range(firstDesc, len(descNames)):
            descriptorVals.append(float(line[col]))
        if databaseID in d:
            d[databaseID].append(descriptorVals)
        if databaseID not in d:
            d[databaseID] = [descriptorVals]


# In[4]:

## Initialize list to store final Boltzmann averaged descriptors for each molecule
Boltzmann = []

## Iterate through each molecule and calculate Boltzmann average of descriptors
for key,item in d.items():
    ## Initialize list to store e^(-dE/(0.001986*298)) for each conformation
    individualList = []

    ## Calculate e^(-dE/(0.001986*298)) for each conformation
    for values in item:   
        ## save dE value of current conformation
        dE = values[0]

        ## add e^(-dE/(0.001986*298)) to list of e^(-dE/(0.001986*298)) and descriptors
        individualList.append((math.e)**(-dE/(0.001986*298)))
        values.append((math.e)**(-dE/(0.001986*298)))
    
    ## Sum e^(-dE/(0.001986*298)) results for all conformations
    individualSum = sum(individualList)

    ## Initialize value to calculate sum of weights of all conformations
    sumOfSum = 0 

    ## Create empty list to store weighted descriptor values for each conformation
    allDescVals = []

    ## Calculate weight and resulting weighted descriptor values for each conformation
    for values in item:
        ## Initialize list to store weighted descriptor values
        descriptorValSums = []
        ## Calculate the energy weight of conformation (= conformation e^(-dE/(0.001986*298)) value / sum of e^(-dE/(0.001986*298)) values for all conformations)
        weight = values[len(values) - 1]/individualSum
        ## Add energy weight of conformation to sum
        sumOfSum += weight

        ## Calculate conformation energy-weighted value for each descriptor
        for val in range(1, len(values)):
            descriptorValSums.append(values[val]*weight)

        ## Add weighted descriptor values to list of all molecule conformations
        allDescVals.append(descriptorValSums)
    
    ## Create dataframe to store all weighted descriptor values for molecule
    dataColumnNames = descNames[firstDesc:len(descNames)]
    dataColumnNames.append('Weight')      
    df = pd.DataFrame(allDescVals, columns = dataColumnNames)

    ## Initialize list to store each Boltzmanned descriptor value for molecule
    moleculesBoltzmannedDescriptors = []
    moleculesBoltzmannedDescriptors.append(key)

    ## Calculate the Boltzmann average of each descriptor (= sum of conformation weighted descriptor value / sum of conformation weights)
    for desc in range(0, len(dataColumnNames)-1):
        descBoltzmann = (df[dataColumnNames[desc]].sum())/sumOfSum
        moleculesBoltzmannedDescriptors.append(descBoltzmann)

    ## Add Boltzmanned descriptors for molecule to list of final Boltzmanned data
    Boltzmann.append(moleculesBoltzmannedDescriptors)


# In[5]:

## Create heading for output file with Ind and descriptor names
finalColNames = ['Ind']
finalColNames += descNames[firstDesc:len(descNames)]


# In[6]:

## Exports all Boltzmanned data to csv file
newDF = pd.DataFrame(Boltzmann, columns = finalColNames)
newDF.to_csv('chunk3_set3_modifiedmols_0506boltzmanned.csv')

