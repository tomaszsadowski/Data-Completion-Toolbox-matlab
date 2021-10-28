Data Completion Toolbox:
------------
Toolbox allows to test and compare methods for Image Completion and Data Completion problems in Matlab.
Presented methods use various Nonnegative Matrix Factorization and Tensor decomposition algorithms. 
It was based on research performed during realization of PhD. List of publications: https://dona.pwr.edu.pl/szukaj/default.aspx?nrewid=802432 . 

Keywords: image processing, data completion, nonnegative matrix factorization, tensor decompositions, low-rank models.

Requierments:
------------
- MATLAB Tensor Toolbox, Sandia Corporation, (https://www.tensortoolbox.org), 
- some algorithms require Image Prcessing Toolbox and Parallel Computing Toolbox.

Installation:
------------
- unpack in desired folder,
- download and copy MATLAB Tensor Toolbox folder to 'methods' folder. 

Usage:
------------
To test, simply run 'main_demo.m'. 
To learn how to change parameters check test_toolbox.m or execute 'help test_toolbox'. Also check help or description in other functions.

Avaiable datasets: 
------------
- 'lena' - 256 x 256 x 3 image (color),
- 'barbara' - 512 x 512 x 3 image (color),
- 'monarch' - 512 x 512 x 3 image (color).

Avaiable methods:
------------  
1 - SmNMF-MC, Smooth Nonnegative Matrix Factorization for Matrix Completion, 2017,
2 - HALS, Hierarchical Alterning Least Squeares for Image Completion, 2017,  
3 - BSA-IC, B-Splines-based Algorithm for Image Completion, 2018,   
4 - ICSA, Image Completion under Separability Assumption, 2018,  
5 - ALS Tucker for Image Completion, 2019,
6 - PTT, Permutation Tensor Train, 2019,
7 - KA-TT, Ket Augmented Tensor Train, 2019,
8 - fALS, filtered ALS for Image Completion, 2019,
9 - TT-IC - Tensor Intrepolation for Image Completion, 2020,
10 - GDC(H), Geometric Data Completion (horizontal), 2021,
11 - GDC(V), Geometric Data Completion (vertical), 2021,
12 - GDC(HV), Geometric Data Completion (Horizontal-Vertical), 2021.

Feedback:
------------
Please send bug reports, comments, or questions to Tomasz Sadowski (mailto:tomasz.sadowski@pwr.edu.pl) or Rafal Zdunek (mailto:rafal.zdunek@pwr.edu.pl).


© Tomasz Sadowski 2018-2021