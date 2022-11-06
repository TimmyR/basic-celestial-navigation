# Celestial Navigation

This program uses functional Python programming to implement a basic celestial navigation routine. Given the 
observed heights of at least two stars above the horizon, measured with a sextant, the program outputs pairs 
of possible positions. The user should easily be able to deduce the correct one given a very rough idea of
their position on the globe, for example that they are in the Atlantic rather than the Indian Ocean. 

The program demonstrates a single computational method for determining ones position from stars, and is a 
building block for the much more complete GUI astal nav program I am currently developping using 
object-oriented programming, which will include several more methods and take as inputs the heights of the 
Sun, Moon, and several planets. 

In order to run, you will need the following packages:

	1. numpy (https://numpy.org/)
	2. astropy (https://www.astropy.org/)
	3. skyfield (https://rhodesmill.org/skyfield/)
	
Please see their documentation regarding installation. 

If the 'nav_star.csv' file (which contains the positions of the navigational stars) is not present or you 
would like to update it, just run 'create_csv_nav_stars.py' to recreate it. 