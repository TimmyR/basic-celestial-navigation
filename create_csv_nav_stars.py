# -*- coding: utf-8 -*-
"""
Script to create a .csv file with the right ascension (ra) and declination (dec)
of all the common navigational stars.


Author: Timothe Rhein
13/12/2020
"""

import csv
from astropy.coordinates import SkyCoord

'''
####
FUNCTIONS
####
'''

def star_position(star) :
    '''
    Finds the star position from an online data base following SkyCoord.
    Converts to string
    '''
    star_position = SkyCoord.from_name(star)                     
    star_position_ra = star_position.ra.to_string(decimal = True)       
    star_position_dec = star_position.dec.to_string(decimal = True)
    
    return star_position_ra, star_position_dec


def nav_star_names(filename) :
    '''
    Import the star names from a file. Takes the first column excluding the 
    header to find the 58 navigational stars.
    '''
    file = []
    
    with open(filename) as csvfile :
        
        csvreader = csv.reader(csvfile, delimiter = ',')
        
        for row in csvreader :
            
            file.append(row[0])
            
    star_names = file[1:]
    
    return star_names
    


def write_csv(filename, star_names) :
    '''
    Writes csv file with ra and dec for the navigational stars from the 
    astropy function astropy.coordinates.SkyCoord.fromname()
    '''
    with open(filename, 'w', newline = '') as file :
        
        writer = csv.writer(file)
        writer.writerow(["Star","RA","Dec"])
        
        for i in range(len(star_names)) :
            
            star_position_ra, star_position_dec = star_position(star_names[i])
            
            writer.writerow([star_names[i] , star_position_ra , star_position_dec])

    
'''
####
MAIN CODE
####
'''


if __name__ == '__main__':
    write_csv('nav_stars.csv', nav_star_names('star_positions_2020.csv'))     



