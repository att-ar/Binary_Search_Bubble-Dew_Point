from numpy import linspace


def two_species_bubble_pt(t_min, t_max, p_total, a, b, x_a, x_b):
    '''
    (int, int, int or float, list, list, float, float) -> dict or -1 for failed case

    This function takes the minimum and maximum of the temperature range you would like to check
        t_min and t_max
    the total pressure the system will be
        p_total
    A list for each of the antoine's equation constants for species "a" and "b" in the order [A,B,C].
        a and b
    the liquid mole fraction of species "a" and "b"
        x_a and x_b

    >>>m = [8.08097, 1582.271, 239.726]
    >>>one_p = [7.74416, 1437.686, 198.463]
    >>>two_species_bubble_pt(0, 110, 760, m, one_p, 0.15, 0.85)
    {'Temp (Celsius)': 89.421,
     'Vapour pressure of species a (mmHg)': 1878.374,
     'Vapour Pressure of species b (mmHg)': 562.576,
     'Solution Pressure': 760}
    '''
    num_checkpoints = (t_max - t_min) * 900
    t_range = linspace(t_min, t_max, num_checkpoints)
    i = 0 #left side index
    j = len(t_range) - 1 #right side index

    while i != j + 1:
        #loop for binary search
        m = (i + j) // 2 # divided the range in two so that the search is log(N)/log(2)

        #Use Raoult's law to trial and error until I find the temperature that matches the vapour pressures
        #to the total pressure

        #t_range[m] is the temperature being checked for validity.
        p_sol = x_a * 10**( a[0] - ( a[1] / (t_range[m] + a[2]) ) )     +    x_b * 10**( b[0] - ( b[1] / (t_range[m] + b[2]) ) )
        #The pressure of the solution mixture
        if p_total - 0.1 > p_sol: #the temperature is below the needed temperature
            i = m + 1 #moving the left index to the middle point, becuase everything below middle point is unneeded
        elif p_total + 0.1 < p_sol:
            j = m - 1 #the opposite of whats happenig above for i

        else:# p_total - 0.1 < p_sol < p_total + 0.1: #if the pressure match up, the temperature has been found
            return {
            "Temp (Celsius)" : round(t_range[m], 3),
            "Vapour pressure of species a (mmHg)" : round(10**( a[0] - ( a[1] / (t_range[m] + a[2]) ) ), 3),
            "Vapour Pressure of species b (mmHg)" : round(10**( b[0] - ( b[1] / (t_range[m] + b[2]) ) ), 3),
            "Solution Pressure" : round(p_total, 3)
            }
    return -1

def two_species_dew_pt(t_min, t_max, p_total, a, b, y_a, y_b):
    '''
    (int, int, int or float, list, list, float, float) -> dict or -1 for failed case

    This function takes the minimum and maximum of the temperature range you would like to check
        t_min and t_max
    the total pressure the system will be
        p_total
    A list for each of the antoine's equation constants for species "a" and "b" in the order [A,B,C].
        a and b
    the vapour mole fraction of species "a" and "b"
        y_a and y_b

    >>>m = [8.08097, 1582.271, 239.726]
    >>>one_p = [7.74416, 1437.686, 198.463]
    >>>two_species_bubble_pt(0, 110, 760, m, one_p, 0.15, 0.85)
    {'Temp (Celsius)': 73.691,
     'Vapour pressure of species a (mmHg)': 1077.759,
     'Vapour Pressure of species b (mmHg)': 289.432,
     'Sum of mol fraction': 0.9993002}
    '''
    num_checkpoints = (t_max - t_min) * 9000
    t_range = linspace(t_min, t_max, num_checkpoints)
    i = 0 #left side index
    j = len(t_range) - 1 #right side index
    #Uses a modification of Raoult's Law

    while i != j + 1:
        #loop for binary search
        m = (i + j) // 2 # divided the range in two so that the search is log(N)/log(2)

        #Use Raoult's law to trial and error until I find the temperature that matches the vapour pressures
        #to the total pressure

        #t_range[m] is the temperature being checked for validity.
        test_goal = y_a / 10**( a[0] - ( a[1] / (t_range[m] + a[2]) ) )     +     y_b / 10**( b[0] - ( b[1] / (t_range[m] + b[2]) ) )

        #The pressure of the solution mixture
        if 0.998 > test_goal * p_total: #the temperature is above the needed temperature
            j = m - 1 #moving the right index to the middle point, becuase everything above middle point is unneeded
        elif 1.002 < test_goal * p_total:
            i = m + 1 #the opposite of whats happenig above for i

        else: #If the sum of the liquid mole fractions equal 1 then it worked and the results are returned
            #I can take this part out of the loop and get a more accurate result, but it is not needed
            return {
            "Temp (Celsius)" : round(t_range[m], 3),
            "Vapour pressure of species a (mmHg)" : round(10**( a[0] - ( a[1] / (t_range[m] + a[2]) ) ), 3),
            "Vapour Pressure of species b (mmHg)" : round(10**( b[0] - ( b[1] / (t_range[m] + b[2]) ) ), 3),
            "Sum of mol fraction" : round(test_goal*p_total, 7)
            }
    return -1

def temperature_pressure_finder(t_min, t_max, a, b, y_a, y_b, x_a, x_b):
    '''
    (int, int, list, list, float, float, float, float) -> dict or -1 for failed case
    Precondition: The more volatile species should be input as species a

    This function takes the minimum and maximum of the temperature range you would like to check
        t_min and t_max
    A list for each of the antoine's equation constants for species "a" and "b" in the order [A,B,C].
        a and b
    the vapour mole fraction of species "a" and "b"
        y_a and y_b
    the liquid mole fraction of species "a" and "b"
        x_a and x_b
    This function, given the parameters above for a 2-species vapour and liquid, both of known composition,
    in equilibrium with each other, will return the temperature at which that system can exist
    This temperature can then be used to find the pressure that the system must be at.

    Note: I combined both the Dew point and Bubble point formulas to obtain the equation needed.
        I substituted P_total from bubble point into dew point.

    >>>a = [7.11714, 1210.595, 229.664]
    >>>w_high = [7.96681, 1668.210, 228.000]
    >>>temperature_pressure_finder(0, 130, a, w_high, 0.7, 0.3, 0.3, 0.7)
    {'Temp (Celsius)': 65.0,
     'Vapour pressure of species a (mmHg)': 1020.347,
     'Vapour Pressure of species b (mmHg)': 187.611,
     'System Pressure of Mixture (mmHg)': 437.432}
    '''
    num_checkpoints = (t_max - t_min) * 9000
    t_range = linspace(t_min, t_max, num_checkpoints)
    i = 0 #left side index
    j = len(t_range) - 1 #right side index
    #Uses a modification of Raoult's Law

    while i != j + 1:
        #loop for binary search
        m = (i + j) // 2 # divided the range in two so that the search is log(N)/log(2)

        #Use Raoult's law to trial and error until I find the temperature that matches the vapour pressures
        #to the total pressure

        #t_range[m] is the temperature being checked for validity.
        test_goal =  x_a * y_b * 10**( a[0] - ( a[1] / (t_range[m] + a[2]) ) )  \
        / 10**( b[0] - ( b[1] / (t_range[m] + b[2]) ) ) \
        + x_b * y_a * 10**( b[0] - ( b[1] / (t_range[m] + b[2]) ) )  \
        / 10**( a[0] - ( a[1] / (t_range[m] + a[2]) ) )

        #The derived equation being checked, if it is valid then the results are returned
        if (1 - x_a*y_a - x_b*y_b) - 0.002 > test_goal: #the temperature is above the needed temperature
            j = m - 1 #moving the right index to the middle point, becuase everything above middle point is unneeded
        elif (1 - x_a*y_a - x_b*y_b) + 0.002 < test_goal:
            i = m + 1 #the opposite of whats happenig above for i

        else:
            return {
            "Temp (Celsius)" : round(t_range[m], 3),
            "Vapour pressure of species a (mmHg)" : round(10**( a[0] - ( a[1] / (t_range[m] + a[2]) ) ), 3),
            "Vapour Pressure of species b (mmHg)" : round(10**( b[0] - ( b[1] / (t_range[m] + b[2]) ) ), 3),
            "System Pressure of Mixture (mmHg)" :\
            round(x_a * 10**( a[0] - ( a[1] / (t_range[m] + a[2]) ) )  +  x_b * 10**( b[0] - ( b[1] / (t_range[m] + b[2]) ) ), 3)
            }
    return -1
