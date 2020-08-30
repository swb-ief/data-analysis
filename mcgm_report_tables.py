#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 18:44:16 2020

@author: wasilaq
"""

import tabula
import numpy as np
# import csv # for saving ward_identifiers dictionary

file = "http://stopcoronavirus.mcgm.gov.in/assets/docs/Dashboard.pdf"


# ward wise new case growth
new_case_growth = tabula.read_pdf(file, pages=24, multiple_tables=False)
new_cases = new_case_growth[0]

# check that correct table was scraped
nc_expected_header = ['Date of report', 'RC', 'HW', 'RS', 'RN', 'PS', 'A', 'C', 'D', 'KW', 'T', 'PN', 'N', 'FN', 'FS', 'MW', 'ME', 'B', 'E', 'GS', 'KE', 'GN', 'S', 'HE', 'L', 'Grand Total']

if (new_cases.columns == nc_expected_header).all():
    pass
else:
    print('Incorrect columns in ward cases table')
    
row_11_values = new_cases.loc[11, nc_expected_header[1:]]
if row_11_values.isnull().all(): # value for each ward should be null
    pass
else:
    print('Unexpected values in row 11')
    
for row in range(1,11):
    if np.array_equal(new_cases.loc[row, nc_expected_header[1:-1]], new_cases.loc[row, nc_expected_header[1:-1]].astype(str)): # value for each ward should be a str
        pass
    else:
        print('Unexpected values in rows 1-10')
        
        
# clean table
new_cases.drop(labels=11, inplace=True)


# save ward identifier with corresponding ward name in a dictionary
# wards = new_cases.iloc[0][1:] # series of ward names, index is identifiers
# identifiers = new_cases.iloc[0][1:].index
# ward_identifiers = {}
# wards.drop('Grand Total', inplace=True)

# for ward in wards:
#     text = str(ward)
#     name = text.replace('\r',' ')
#     index = wards[wards==text].index[0]
#     ward_identifiers[index] = name

# len(ward_identifiers) == len(wards) # check that dictionary has all wards

new_cases.drop(labels=0, inplace=True)


# elderly screening data
elderly_screening = tabula.read_pdf(file, pages=21, multiple_tables=False)
elderly = elderly_screening[0]

# check that correct table was scraped
e_expected_header = ['Wards', 'Total No. of Houses', 'Population', 'Total No. of', 'Sr Citizen', '* Sr Citizen', 'Unnamed: 6', 'Unnamed: 7']
if (elderly.columns == e_expected_header).all():
    pass
else:
    print('Incorrect columns in elderly table')
    
row_1_values = elderly.loc[1, e_expected_header[1:]]
if row_1_values.isnull().all(): # value for each ward should be null
    pass
else:
    print('Unexpected values in row 1')
    
for column in e_expected_header[:4]:
    if np.array_equal(elderly[column][2:], elderly[column][2:].astype(str)):
        pass
    else:
        print('Unexpected values in ' + str(column))
    
if np.array_equal(elderly['Unnamed: 6'][2:], elderly['Unnamed: 6'][2:].astype(float)):
    pass
else:
    print('Unexpected values in ' + str(column))
    
    
# clean table
elderly.loc[2, 'Wards'] = 'Daily Totals'
elderly.drop(labels=[0,1], inplace=True)
elderly.drop(columns=['* Sr Citizen','Unnamed: 7'], inplace=True)
elderly.columns = [
    'Wards', 'Total No. of Houses Surveyed', 'Population Covered', 'Total No. of Senior Citizens', 'Sr Citizen SPO2>95','Sr Citizen SPO2<95'
]
elderly.index = elderly['Wards']
elderly.drop(columns=['Wards'], inplace=True)



# save as csv files
saved_files_path = '/Users/wasilaq/SWB/data-analysis/'

new_cases.to_csv(saved_files_path + 'new_cases.csv')
elderly.to_csv(saved_files_path + 'elderly.csv')

# save ward identifiers dictionary
# with open('ward_identifiers.csv', 'w') as f:
#     w = csv.DictWriter(f, ward_identifiers.keys())
#     w.writeheader()
#     w.writerow(ward_identifiers)