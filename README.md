# Binary_Search_Bubble-Dew_Point

This script was made because my **CHE101: Chemical Engineering Concepts 2** course at the University of Waterloo
introduced the concept of Bubble and Dew point for a two-component mixture.

In certain problems, a trial and error process (or using MS Excel Goal Seek/Solver) was required to find the Bubble or Dew Point temperature
of the two-component system.

I decided to make this python script since I found it easier to use and I wanted to make a binary search function.
The binary search loop was heavily inspired from _Practical Programming, Third Edition: An Introduction to Computer Science Using Python 3.6 by Paul Gries et al._
Note that it requires Antoine's Equation's constants for the species being analyzed which I get from _Elementary Principles of Chemical Processes by Richard M. Felder_.

Uses numpy.linspace
