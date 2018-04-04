#!/usr/bin/python3

# -*- coding: utf-8 -*-

"""Quartic reaction coordinate
    
    This module contains one function, which allows to perform an optimization of the coefficients
    in the polinomial expansion of quartic reaction coordinate model (a*x**4 + b*x**3 + c*x**2)
    and another function, which plot the data.
    
    Example
    -------
    Import this module to run a coefficient optimization for any chemical reaction. You need to
    give the activation and reaction energy.
    
    from quartic_rxn_opt import quarticrxn
    
    The result will be a list, which contains the three optimized coefficients (a, b, c) and it
    will be saved a png file with the energy profile.
    
    """

from sympy import *
import matplotlib.pyplot as plt
from matplotlib import rcParams

DEFAULT_NAME = 'test'
WANT_PLT = 'no_save'

def quarticrxn(eact, erxn):
    """Perform a coefficient optimization of the quartic reaction coordinate.
        
        Parameters
        ----------
        eact: float
        Activation energy of the chemical reaction (E_transition_state - E_reactant).
        erxn: float
        Reaction energy of the chemical reaction (E_product - E_reactant).
        
        Returns
        -------
        The optimized functions will be in a list (a, b, c).
        
        """
    a, b, c = symbols('a, b, c')
    xts = ((-6 * b) / (8 * a) - 1)
    func1 = a + b + c - erxn
    func2 = 4 * a + 3 * b + 2 * c
    func3 = a * xts ** 4 + b * xts ** 3 + c * xts ** 2 - eact
    # Matrix of functions
    functions_matrix = Matrix([func1, func2, func3])
    # Inverted Jacobian
    jacobian_inverted = (functions_matrix.jacobian([a, b, c])) ** (-1)
    # Matrix of coefficients
    coef_matrix = Matrix([100, -200, 200])
    while (functions_matrix.subs(
            [(a, coef_matrix[0]),
             (b, coef_matrix[1]),
             (c, coef_matrix[2])]).norm() > 1e-8):
        # Iterative equation to be solved
        coef_matrix = coef_matrix - (jacobian_inverted.subs([(a, coef_matrix[0]),
                                                             (b, coef_matrix[1]),
                                                             (c, coef_matrix[2])]
                                                            )
                                     * functions_matrix.subs([(a, coef_matrix[0]),
                                                              (b, coef_matrix[1]),
                                                              (c, coef_matrix[2])]
                                                             )
                                     )
    optimized_coefficients = [float(coef_matrix[0]), float(coef_matrix[1]), float(coef_matrix[2])]
    return optimized_coefficients



def quartic_plot(a,b,c,name=DEFAULT_NAME,want_plot=WANT_PLT):
    """Function that allows to plot the data.
        
        Parameters
        ----------
        a: callable
        Coefficient a obtained from the quarticrxn function.
        b: callable
        Coefficient b obtained from the quarticrxn function.
        c: callable
        Coefficient c obtained from the quarticrxn function.
        name: string
        Name of png file that contains the energy profile.
        The default is test.
        
        Returns
        -------
        An image (name.png) will be saved in the same folder where you run the code.
        
        """
    
    def frange(ini, fin, delta):
        j = ini
        while j < fin:
            yield j
            j += delta
    xaxis = list(frange(0, 1.01, 0.01))
    # Contruct the function Y (Energy profile)
    yaxis = []
    for n in frange(0, 1.01, 0.01):
        yaxis.append(a * n ** 4
                        + b * n ** 3
                        + c * n ** 2
                        )
    rcParams['font.family'] = 'serif'
    plt.plot(xaxis, yaxis)
    plt.ylabel('Potential Energy', fontsize=12)
    plt.xlabel('Reaction coordinate', fontsize=12)
    if want_plot == 'save':
        return plt.savefig(name+'.png', dpi=1080)
    else:
        plt.show()
    return


