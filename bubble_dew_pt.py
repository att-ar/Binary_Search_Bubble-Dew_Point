from numpy import linspace

#Antoine's constants: [A, B, C]
a = [7.11714, 1210.595, 229.664] #acetone
e = [8.11220, 1592.864, 226.184] #ethanol
#liquid composition (mol fraction)
x_a = 0.7
x_e = 0.3

def two_species_bubble_pt(t_min, t_max, p_total, a, b, x_a, x_b):
    '''
    (int, int, int or float, list, list, float, float) -> tuple(int, float) or -1 for failed case

    This function takes the minimum and maximum of the temperature range you would like to check
        t_min and t_max
    the total pressure the system will be in mmHg (millimeters of Mercury) or the equivalent : torr
        p_total
    A list for each of the antoine's equation constants for species "a" and "b" in the order [A,B,C].
        a and b
    the liquid mole fraction of species "a" and "b"
        x_a and x_b
    
    It returns a tuple holding the Temperature that satisfies Raoult's Law
    and the total pressure that the temperature results in.
    
    #Antoine's constants: [A, B, C]
    >>>a = [7.11714, 1210.595, 229.664] #acetone
    >>>e = [8.11220, 1592.864, 226.184] #ethanol
    >>>two_species_bubble_pt(50, 100, 760, a, b, 0.7, 0.3)
    (61.87395832281134, 760.0728705051014)
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

        else:# p_total - 0.5 < p_sol < p_total + 0.5: #if the pressure match up, the temperature has been found
            return (t_range[m], p_sol)
    return -1

def two_species_dew_pt(t_min, t_max, p_total, a, b, y_a, y_b):
    '''
    (int, int, int or float, list, list, float, float) -> float or -1 for failed case

    This function takes the minimum and maximum of the temperature range you would like to check
        t_min and t_max
    the total pressure the system will be
        p_total
    A list for each of the antoine's equation constants for species "a" and "b" in the order [A,B,C].
        a and b
    the vapour mole fraction of species "a" and "b"
        y_a and y_b
    
    It returns a tuple holding the Temperature that satisfies the condition.
    and the sum of the liquid mole fractions which should be close to the integer 1.
    
        
        #Antoine's constants: [A, B, C]
    >>>a = [7.11714, 1210.595, 229.664] #acetone
    >>>e = [8.11220, 1592.864, 226.184] #ethanol
    >>>two_species_dew_pt(50, 100, 760, a, b, 0.7, 0.3)
    (66.06428895382723, 1.0006609870410603)
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

        else:# p_total - 0.5 < p_sol < p_total + 0.5: #if the pressure match up, the temperature has been found
            return (t_range[m], test_goal*p_total)
    return -1

