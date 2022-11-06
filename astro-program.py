# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 14:55:55 2020

Code to calculate position based on Celestial Navigation.

Given the heights above the horizon of at least 2 stars, will output the 
user's possible latitudes and longitudes. The actual position can be inferred
from the user's most likely position (they know they are in the Atlantic
rather than the Indian Ocean, hopefully) or by verifying which position is 
repeated when the heights of more than 2 stars are used. 

@author: Timothe Rhein
"""
import sys
import csv
import numpy as np
from itertools import combinations
from scipy.optimize import fsolve
from astropy import units as u
from astropy.coordinates import Angle

from skyfield.api import Loader
load = Loader('.\skyfield_data')    #location of data downloaded using SkyField


### FUNCTIONS ###


#General


def check_exit(inp) :
    '''
    General function that exits the program.

    Parameters
    ----------
    inp : Input that is checked. If the input is 'quit' or 'exit', the program
            will be terminated.

    Returns
    -------
    None.

    '''
    if inp.lower() == 'quit' or inp.lower() == 'exit' :
        sys.exit('You have exited the program. Goodbye :)')
    else :
        pass


def value_check(value, minimum, maximum) :
    '''
    Checks that a value is within a certain range. Can be used for example to
    check that dates are correctly input.

    Parameters
    ----------
    value : number. The value that needs to be checked.
    minimum : number. The minimum value possible. Is considered an acceptable
                value.
    maximum : number. The maximum value possible. Is considered an acceptable
                value.

    Returns
    -------
    True if the value is within the correct range, False otherwise.

    '''
    if value < minimum or value > maximum : return False

    else : return True


#Load timescale

def load_timescale() :
    '''
    Tries to download an updated timescale file from online. Otherwise uses
    builtin file from the location '~/skyfield data'.

    Returns
    -------
    ts : SkyField timescale.

    '''
    try :
        ts = load.timescale(builtin = False)    #try downloading updated 
                                                #timescale file

    except :
        ts = load.timescale()                   #otherwise use builtin 
                                                #timescale file

    return ts


#Input

def input_sextant_altitude() :
    '''


    Returns
    -------
    None.

    '''
    while True :
        sextant_altitude_str = input(
            "    -Input the altitude reading for this star: ")  
        #TODO: make input criteria clear(though this is working pretty well)
        check_exit(sextant_altitude_str)    #check for input exit conditions

        try :
            #transform the string into an angle.
            sextant_altitude = Angle(sextant_altitude_str, unit = 'degree') 

        except ValueError :
            print('Invalid format, please try again.')
            continue

        if value_check(sextant_altitude.degree, 0, 90) :
            break

        else :
            print('Please enter an angle between 0 and 90 degrees.')

    #TODO: decide whether to output as a degree or in Angle() format
    return sextant_altitude.degree     


def input_datetime(ts) :
    '''
    User inputs date and time and the function outputs a Skyfield time value
    in the Terrestrial Time format.

    Parameters
    ----------
    ts : SkyField timescale. Should be loaded at the start of the program.

    Returns
    -------
    measurement_time : SkyField time value. Time at which the measurements were
                        taken. Terrestrial Time (TT) by default and displayed
                        as a Julian Date.

    '''
    print('Please input the date and time the star measurements were',
          ' taken at as integers below.')
    unit = ['year', 'month', 'day', 'hour', 'minute', 'second']

    while True :
        #datetime will contain values corresponding to 'unit' array after 
        #itteration
        datetime = []
        
        #Loop over inputs for each unit of time
        for i in unit :
            while True :
                try :
                    value = input('{}: '.format(i.title()))
                    check_exit(value)       #check for input exit conditions
                    value = int(value)
                    datetime.append(value)

                except ValueError:
                    print('ValueError: Please input an integer value!')
                    continue

                break

        #Check datetime has the right values ie. there is no 60th second, this 
        #adds another minute
        if (value_check(datetime[1],1,12) and value_check(datetime[2],1,31) and 
            value_check(datetime[3],0,23) and value_check(datetime[4],0,59) and 
            value_check(datetime[5],0,59) ):
            break

        else :
            print('\nPlease enter a valid date and time please. The format is',
                  ' 24 hour time.')

    #ts.utc(year,month,day,hour,minute,second)
    measurement_time = ts.utc(datetime[0],datetime[1],datetime[2],datetime[3],
                              datetime[4],datetime[5])

    return measurement_time


#Star data

def import_nav_stars(filename = 'nav_stars.csv') :
    '''
    Imports the file containing the right ascension (ra) and declination (dec)
    of the navigational stars. These are calculated using the astropy SkyCoord
    function which takes data (J200 position) from the Sesame name resolver:
        https://cds.u-strasbg.fr/cgi-bin/Sesame

    This program assumes the .csv file to contain 3 columns ordered
    Star, RA, Dec.

    #??? Might need to find a new source for stellar tables for increased
    precision since there is a slight variance in the positions of navigational
    stars. I beleive Sesame only gives the position on the day it is accessed.
    Try by using PyEphem.

    Parameters
    ----------
    filename : string. Must be a .csv file.

    Returns
    -------
    star_values : dictionary. Star names are the keys and the values are
                    [RA, Dec].

    '''
    file = []

    with open(filename) as csvfile :

        csvreader = csv.reader(csvfile, delimiter = ',')

        for row in csvreader :

            file.append(row)

    #Separate the star values into seperate lists, ignoring the first row.
    star_names = [file[i+1][0] for i in range(len(file)-1)]     
    star_ra_dec = [file[i+1][1:] for i in range(len(file)-1)]   # [RA, Dec]

    star_values = dict( zip( star_names , star_ra_dec ))

    return star_values


def number_of_stars() :
    '''
    Asks the user for the amount of stars that are being used to find your
    position.

    Returns
    -------
    number_of_stars : integer. The number of stars measurements have been
                        taken for.

    '''
    # TODO : make sure they input an integer greater than zero and less than
    # 10 cause that's just ridiculous...

    number_of_stars = input('How many stars have you taken measurements for? ')
    check_exit(number_of_stars)         #check for input exit conditions
    number_of_stars = int(number_of_stars)

    return number_of_stars


def finding_stars(number_of_stars) :
    '''
    Takes star ra and dec from the file nav_star.csv, for the stars that the
    observer has measured.

    Parameters
    -------
    number_of_stars : integer. The number of stars, found by input by the
                        function of the same name.

    Returns
    -------
    star_ra_dec : numpy array of size number_of_stars x 2. Each row contains
    [RA, Dec].

    altitudes : numpy array of size number_of_stars x 2. Contains sextant
    altitudes for each star.

    '''
    star_values = import_nav_stars('nav_stars.csv')

    #Values of the star will be held in a nx2 matrix. Each row will be [RA,Dec].
    star_ra_dec_list = []
    altitudes_list = []

    for star in range(number_of_stars) :
        while True :
            try :
                name = input('What is the star number {}? '.format(star+1))
                check_exit(name)
                values = star_values[name.lower()] 
                star_ra_dec_list.append(values)

            except KeyError :
                print('This star name was not found, please try again...')
                continue
            break

        #Input sextant altitude for each star
        sextant_altitude = input_sextant_altitude()
        altitudes_list.append(sextant_altitude)

    # star_ra_dec contains the values of the selected stars
    star_ra_dec = np.array(star_ra_dec_list)    
    altitudes = np.array(altitudes_list)

    return star_ra_dec, altitudes


#Conversion

def hrs_to_deg(hrs_value) :
    '''
    Converts a float value from hour-angle to degrees.

    Parameters
    ----------
    hrs_value : scalar. This will be turned into a float if it isn't already.

    Returns
    -------
    deg_value : float. Units are decimal degrees.

    '''
    deg_value = float(hrs_value) * (360/24)

    return deg_value


def deg_to_hrs(deg_value) :
    '''
    Converts a float value from degrees to hour-angles.

    Parameters
    ----------
    deg_value : scalar. This will be turned into a float if it isn't already.

    Returns
    -------
    hrs_value : float. Units are decimal hours.

    '''
    hrs_value = float(deg_value) * (24/360)

    return hrs_value


def cart_to_global(cartesian_coords) :
    '''
    Converts a position in cartesian coordinates into latitude and longitude,
    without taking into account elevation.

    Parameters
    ----------
    cartesian_coords : array. Contains the values in the order x, y, z.

    Returns
    -------
    latitude : float. Latitude of the position.
    longitude : float. Longitude of the position.

    '''
    # TODO: maybe using arctan2 instead of arctan will simplify this function.

    x, y, z = cartesian_coords

    #latitude = np.arctan(z / np.sqrt(x**2 + y**2))
    latitude_rad = np.arcsin(z)

    #Conditions to prevent ZeroDivisionError
    if x == 0 and y > 0 :
        longitude_rad = np.pi / 2
    elif x == 0 and y < 0 :
        longitude_rad = - np.pi / 2 #negative since its the opposite quadrant
    elif x == 0 and y == 0 :
        longitude_rad = 0
    #Conditions to make sure the position is in the right quadrant
    elif x < 0 and y <= 0:
        longitude_rad = np.arctan(y/x) - np.pi
    elif x < 0 and y >= 0:
        longitude_rad = np.arctan(y/x) + np.pi
    else:
        longitude_rad = np.arctan(y/x)

    latitude, longitude = np.degrees([latitude_rad, longitude_rad])
    return latitude, longitude


#Equations

def find_substellar_points(star_ra_dec, measurement_time) :
    '''
    Finds the substellar points of the previously selected stars whose RA and
    Dec are contained within the star_ra_dec array. The measurement_time is
    used to find the GHA of the vernal equinox (GHA_vernal) which is then used
    to find the longitude of the substellar point via:
        longitude = RA - GHA_vernal .
    The latitude is the same as the declination.

    Position of the vernal equinox is assumed to be most precise using the
    .gast attribute which outputs the Greenwich Apparent Sidereal Time.

    # ??? Need to check why there is a discrepancy between this subpoint and
    the one calculated using the skyfield subpoint() routine. Maybe it's the
    shape of the Earth? I'm assuming the Earth is spherical which isn't quite
    right but I need it that way for the method I'm using.

    Parameters
    ----------
    star_ra_dec : numpy array of size number_of_stars x 2. Each row contains
                    [RA, Dec].

    measurement_time : SkyField time value. Time at which the measurements were
                        taken. Terrestrial Time (TT) by default and displayed
                        as a Julian Date.

    Returns
    -------
    star_subpoints : numpy array of size number_of_stars x 2. Each row contains
                    [latitude, longitude].

    '''
    star_subpoints_list = []

    for i in range(len(star_ra_dec)) :
        latitude = float(star_ra_dec[i,1])   #Same as declination of the star

        #Longitude of the substellar point depends on the GHA of the vernal 
        #equinox (GHA_vernal)
        GHA_vernal = hrs_to_deg(measurement_time.gast)   
        longitude = float(star_ra_dec[i,0]) - GHA_vernal + 360 

        subpoint = [latitude, longitude]
        star_subpoints_list.append(subpoint)

    star_subpoints = np.array(star_subpoints_list)

    return star_subpoints


def plane_functions(variables, star_subpoints, altitudes, pair = [0,1]) :
    '''
    Returns an array of functions that each correspond to the equation of the
    plane of the circle of position of the observer. The circle of position is
    deduced based on the star's subpoint and the altitude of the star as
    measured by the sextant.

    Parameters
    ----------
    variables : list of floats of length 3. Corresponds to the variables
    x, y, and z.

    star_subpoints : numpy array of size number_of_stars x 2. Each row contains
    [latitude, longitude].

    altitude : array of length number_of_stars. Corresponds to the sextant
    altitude measurement for each star.

    pair : array of length 2. Defines for which two stars the position is
    being solved, as it contains the indices of the 2 stars to solve for.
    The default is [0,1] for the case of only 2 stars.

    Returns
    -------
    functions : list of length 2. Contains the function of the
    plane containing the circle of position deduced from each star of the pair.

    '''
    x, y, z = variables
    functions = []

    # Convert all angles to radians for trigonometry
    star_subpoints_rad = np.radians(star_subpoints)
    altitudes_rad = np.radians(altitudes)

    #Now outputs the plane_functions of only the stars defined in pair.
    for i in pair:
        latitude = star_subpoints_rad[i,0]
        longitude = star_subpoints_rad[i,1]

        a_i = np.cos(longitude) * np.cos(latitude)
        b_i = np.sin(longitude) * np.cos(latitude)
        c_i = np.sin(latitude)

        #zenith is the angle of the star to the vertical of the observer
        zenith = np.pi/2 - altitudes_rad[i] 
        p_i = np.cos(zenith)

        function = a_i*x + b_i*y + c_i*z - p_i      
        functions.append(function)

    return functions


def unit_sphere_function(variables):
    '''
    Function corresponding to the unit sphere.

    Parameters
    ----------
    variables : list of floats of length 3. Corresponds to the variables
    x, y, and z.

    Returns
    -------
    function : Single function corresponding to the unit sphere.

    '''
    x, y, z = variables
    function = x**2 + y**2 + z**2 - 1

    return function


def simultaneous_functions(variables, star_subpoints, altitudes, pair = [0,1]) :
    '''
    Joins the functions defined by plane_functions(...) and
    unit_sphere_function(...) to create a list of all functions to be solved
    simultaneously.

    Parameters
    ----------
    variables : list of floats of length 3. Corresponds to the variables
    x, y, and z.

    star_subpoints : numpy array of size number_of_stars x 2. Each row contains
    [latitude, longitude].

    altitude : array of length number_of_stars. Corresponds to the sextant
    altitude measurement for each star.

    pair : array of length 2. Defines for which two stars the position is
    being solved, as it contains the indices of the 2 stars to solve for.
    The default is [0,1] for the case of only 2 stars.

    Returns
    -------
    simultaneous_functions : list of length (number_of_stars + 1). Contains the
    function of the plane containing the circle of position deduced from each
    star, as well as the function of the unit sphere.

    '''
    planes = plane_functions(variables, star_subpoints, altitudes, pair)
    unit_sphere = unit_sphere_function(variables)

    #Concatenate lists - need to make unit_sphere a list before adding them.
    simultaneous_functions = planes + [unit_sphere]

    return simultaneous_functions


def simultaneous_equations_solver(star_subpoints, altitudes) :
    '''
    Solves the simultaneous equation
                    a_i*x + b_i*y + c_i*z - p_i = 0
    where i is the star number and goes from 1 to number_of_stars, as well as
    the equation
                    x**2 + y**2 +z**2 -1 = 0
    which is the condition that the position is on the surface of the globe.

    # TODO : would be nice to also be able to identify which stars are
    # involved in each pair

    # TODO: Choose the observer's position by comparing each pair of solutions 
    # and seeing which point is consistent throughout.


    Returns
    -------
    None.

    '''
    #Uses 2 points on opposite sides of the globe as INITIAL_GUESS in order to 
    #resolve for both points
    INITIAL_GUESS = np.array([ [1,0,0] , [-1,0,0] ])

    if len(star_subpoints) < 2 :
        print('Please choose more than 2 stars in order to resolve a position.')
        return None

    elif len(star_subpoints) == 2 :

        #TODO Need to check that this method actually ALWAYS gives 2 distinct 
        #solutions, otherwise will need to switch methods.
        sol1 = fsolve(simultaneous_functions, INITIAL_GUESS[0], 
                      (star_subpoints, altitudes))
        sol2 = fsolve(simultaneous_functions, INITIAL_GUESS[1], 
                      (star_subpoints, altitudes))

        return np.array([sol1, sol2])

    else :
        #Create list of all possible pairs of stars
        comb = combinations(range(len(star_subpoints)),2)
        comb_array = np.array([i for i in comb])

        solutions_list = []

        for pair in comb_array :

            sol1 = fsolve(simultaneous_functions, INITIAL_GUESS[0], 
                          (star_subpoints, altitudes, pair))
            sol2 = fsolve(simultaneous_functions, INITIAL_GUESS[1], 
                          (star_subpoints, altitudes, pair))

            solutions_list.append([sol1,sol2])

        solutions = np.array(solutions_list)

        return solutions


def lat_long_positions(positions_cart) :
    '''


    Parameters
    ----------
    cartesian_coords : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    positions_list = []

    #Convert each position into lat/long
    for i in range(len(positions_cart)) :
        temp_pair = []
        for j in range(len(positions_cart[i])) :
            position = cart_to_global(positions_cart[i,j])
            temp_pair.append(position)
        positions_list.append(temp_pair)

    positions = np.array(positions_list)

    return positions




'''
####
MAIN CODE
####
'''
run_query = input('Run code? [y/n]\n    ')

if run_query == 'y' :
    pass
else :
    sys.exit()

#Load timescale
ts = load_timescale()

#Define varaibles : Time and altitudes
measurement_time = input_datetime(ts)
#altitudes = 90 - np.array([70.45,61.50,26.87,48.02])
#[19.55, 28.5 , 63.13, 41.98]

#Define measured stars
n = number_of_stars()
star_ra_dec, altitudes = finding_stars(n)  
star_subpoints = find_substellar_points(star_ra_dec, measurement_time)

#Solve equations
positions_cart = simultaneous_equations_solver(star_subpoints, altitudes)
positions = lat_long_positions(positions_cart)

print(positions)
