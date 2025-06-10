MATLAB CODES:
------------
rd2svds.m
trrsvds.m
driver_trrsvds_rd2svds.m 

AUTHORS:
-------
James Baglama
Department of Mathematics and Applied Mathematical Sciences
University of Rhode Island
E-mail: jbaglama@uri.edu

Vasilije Perovic
Department of Mathematics and Applied Mathematical Sciences
University of Rhode Island
E-mail: perovic@uri.edu

Jennifer Picucci
Department of Mathematics and Applied Mathematical Sciences
University of Rhode Island
E-mail: jenniferpicucci@uri.edu and,
U.S. Army Engineer Research and Development Center
E-mail: Jennifer.r.picucci@usace.army.mil

REFERENCE:
---------
Hybrid Iterative Refined Restarted Lanczos Bidiagonalization Method 
submitted for publication to Numerical Algoritms. Additional 
information can be found on the website
http://www.math.uri.edu/~jbaglama/

SOFTWARE REVISION DATE:
----------------------
v1.0, July 2021

SOFTWARE LANGUAGE:
-----------------
MATLAB R2021a

SOFTWARE:
--------
1. Introduction:
The routines rd2svds.m and trrsvds.m are MATLAB programs for computing a few singular 
values and singular vectors of a m x n matrix A. The MATLAB code trrsvds.m is a hybrid 
restarted Lanczos bidiagonalization method using thick-restarting and restarting with 
linear combination to compute the k largest (or smallest) singular values and associated 
vectors. Although trrsvds.m can be used to directly compute the 
smallest singular values it is not recommended. trrsvds.m is best suited 
for computing few of largest singular values while using limitied memory. 
The MATLAB code rd2svds.m is a restarted Lanczos bidiagonalization method 
using deflation to compute the k largest singular values and associated vectors. 
The driver routine, driver_trrsvds_rd2svds.m is a user interactive driver that calls
rd2svds.m and trrsvds.m for different numerical examples presented in the referenced paper. 

2. Installation and Setup:
Download rd2svds.m, trrsvds.m, and driver_trrsvds_rd2svds.m and make sure
the files are in the MATLAB path. Optional execute the following commands: 
>> websave('trrsvds.m', 'http://www.math.uri.edu/~jbaglama/software/trrsvds.m')
>> websave('rd2svds.m', 'http://www.math.uri.edu/~jbaglama/software/rd2svds.m')
>> websave('driver_trrsvds_rd2svds.m', 'http://www.math.uri.edu/~jbaglama/software/driver_trrsvds_rd2svds.m')

3. Usage for trrsvd.m:
Execute the following command for matrix A already in MATLAB:
>> [U,S,V,STATS] = trrsvd(A); 
If the m x n matrix A is a function, 'Afunc',  then the structure of Afun
must be y = Afunc(x,m,n,'transpose') where  transpose = 'F', 
then y = A*x and if transpose = 'T', then y = A'*x.
Execute the following command:
>> [U,S,V,STATS] = trrsvds('Afunc',m,n);

To change parameters, e.g. number of singular triplets (k) to 3 and basis size (m_b) 
to 5, use a struct variable, OPTS.k = 3 and OPTS.m_b = 5
>> [U,S,V,STATS] = trrsvd(A,OPTS);   
or 
>> [U,S,V,STATS] = trrsvds('Afunc',m,n,OPTS);

For other parameter choices and output options, execute the following command:
>> help trrsvds

4. Usage for rd2svds.m:
Execute the following command for matrix A already in MATLAB:
>> [U,S,V,FLAG] = rd2svds(A,m,n,P1,k,tol);
where [m,n] = size(A), P1 is n x 1 starting vector, k is the number of
desired singular triplets and tol is the tolerance.
If the m x n matrix A is a function, 'Afunc'(see 3 above) then 
execute the following command:
>> [U,S,V,FLAG] = rd2svds('Afunc',m,n,P1,k,tol);
For information on rd2svds execute the following command:
>> help rd2svds

5. Usage for driver program driver_trrsvds_rd2svds.m
execute the following command:
>> driver_trrsvds_rd2svds
The program will ask for user to select matrix (options are some of the matrices 
used in the paper), if the matrix is not in the MATLAB path, driver will use 
websave to download the matrix from the SuiteSparseMatrix Collection 
https://sparse.tamu.edu/, if fails to download the driver  will use the 
default A=diag(1:500). The user can download the matrix directly. The user can then 
select k and basis size values (values are restricted to values used in the reference 
paper). The driver will output results in table format.

Please report any bugs or send comments via email to the authors.

***************************************************************************
* THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
* WARRANTY. ALL CODES ARE ILLUSTRATIONS OF THE ALGORITHMS IN THE PAPER    *
* HYBRID ITERATIVE REFINED RESTARTED LANCZOS BIDIAGONALIZATION METHOD.    *
* THE CODES ARE NOT OPTIMIZED FOR PERFORMANCE OR SET UP FOR COMMERICAL    *
* USE. ANY USE BEYOND ILLUSTRATIVE PURPOSES OR PUBLIC USE REQURIES        *
* CONSENT OF THE AUTHORS.                                                 *
***************************************************************************

