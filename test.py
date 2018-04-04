#!/usr/bin/python3


'''
    This module calls the quartic_rxn function that found the
    optimized coefficient of the equation a*x**4 + b*x**3 + c*x**2.
    It requires the activation energy and the reaction energy of the
    chemical reaction.

    Example
    -------
        Consider a reaction which has an activation energy of 25 kcal/mol
        and a reaction energy of -40 kcal/mol.
        The function must be called as:
        a,b,c = quarticrxn(25,-40)
        
    On the other hand the quartic_plot function allows to obtain a plot
    of the energy profile from the optimized a, b and c.
    If you want to save a plot image, you need to specify when the function is called.
    A png file will be saved. The name of the .png file must be specified when
    function is called. This is the default
    Else, the plot will be showed in the screen. This is the default.

    Example
    -------
        If we want to save a file called reaction1.png we must call the function
        in that way:
        
        quartic_plot(a,b,c,'reaction1','save')
        
        If we only want to see the profile, we don not to specify anything,
        because the default is no_save.
        
        quartic_plot(a,b,c,'reaction1')
        
'''

from sympy import *
from quartic_rxn_opt import quarticrxn
from quartic_rxn_opt import quartic_plot


a, b, c = quarticrxn(25,-10)
quartic_plot(a,b,c,'example','save')


