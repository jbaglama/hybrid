function driver_trrsvds_rd2svds
fprintf( [
'*******************************************************************\n' ...
'Driver:                                                            \n' ...
' Hybrid: Thick-Restarted and Restarted SVDS - Algorithm 2: trrsvds \n' ...      
' Hybrid: Restarted Deflation (2×2) SVDS     - Algorithm 3: rd2svds \n' ...
' Paper:                                                            \n' ...
' Baglama, J, Perovic, V, and Picucci, J, "Hybrid Iterative Refined \n' ...
' Restarted Lanczos Bidiagonalization Methods",  2021 submitted      \n' ... 
' Numerical Algorithms, preprint:                                   \n' ...
' http://www.math.uri.edu/~jbaglama/paper34.pdf                     \n' ...
' Matlab Codes                                                      \n' ...
' 1. rd2svds                                                        \n' ...
' 2. trrsvds(NOR,RITZ) - Hybrid method - Iterative refined on       \n' ...
'    NORmal sys. A^TA & thick-restart with RITZ                     \n' ...
' 3. trrsvds(AUG,RITZ) - Hybrid method - Iterative refined on       \n' ...
'    AUGmented sys. [0 A; A^T 0] & Thick-restart with RITZ          \n' ...
' 4. trrsvds(NOR,HARM) - Hybrid method - Iterative refined on       \n' ...
'    NORmal sys. A^TA & Thick-restarted with HARMoinc Ritz          \n' ...
' 5. trrsvds(AUG,HARM) - Hybrid method - Iterative refined on       \n' ...
'    AUGmented sys. [0 A; A^T 0] & Thick-restart with HARMonic Ritz \n' ...
' 6. trrsvds(THK,RITZ) - Non-Hybrid - No iterative refined          \n' ...
'    Only THicK-restart with RITZ                                   \n' ...
' 7. trrsvds(THK,RITZ) - Non-Hybrid - No iterative refined          \n' ...
'    Only THicK-restart with HARMonic Ritz                          \n' ...
'*******************************************************************\n'])
disp(sprintf(' '));

% *************************************************************************
% Programs tested on MATLAB version R2021a
% May not run correctly on older versions or on student version
%  DATE MODIFIED: 12/21/21
%  VER:  1.2
%  
% AUTHORS: 
% James Baglama     email: jbaglama@uri.edu
% Vasilije Perovic  email: perovic@uri.edu
% Jennifer Picucci  email: jenniferpicucci@uri.edu
% *************************************************************************

% Clear all variables before starting
clear all; 

% *************************************************************************
% This driver program will look for largest singular triplets for the
% following  matrices some are typically used in literature 
% from the SuiteSparse Matrix Collection https://sparse.tamu.edu/ 
% Matrices from SuiteSparse Matrix must be loaded in MATLAB path or 
% will be donwload via MATLAB command websave. Internet connection required
% for download websave. No internet connection only diagonal example
% will run. 
% 1. diag(1:500)   500 x 500      
% 2. illc1033     1033 x 320      
% 3. dwt1242     1,242 x 1,242
% 4. well1850    1,850 x 712
% 5. pde2961      2961 x 2961
% 6. amazon03002 262,111 x 262,111  - May run slow dependent on computer speed

% *************************************************************************
disp(sprintf('1. diag(1:500)   500 x 500'));
disp(sprintf('2. illc1033     1033 x 320'));
disp(sprintf('3. dwt1242     1,242 x 1,242'));
disp(sprintf('4. well1850    1,850 x 712'));
disp(sprintf('5. pde2961      2961 x 2961'));
disp(sprintf('6. amazon03002 262,111 x 262,111 - may run slow - depends on computer'));
disp(sprintf('7. well1850    1,850 x 712 - (Smallest - SS)'));
disp(sprintf(' '));
prompt = 'Select matrix (1 - 6) for LS or (7) for SS: ';
example = input(prompt);
if isempty(example), example = 1; end
if ~isnumeric(example)
    error('ERROR: Select matrix from given values.'); 
end
if ~any(example == [1 2 3 4 5 6 7]) 
    error('ERROR: Select matrix from given values.'); 
end 
disp(sprintf(' '));
prompt = 'Select number of singular triplets (k) - 1, 2, 3, or 4: ';
k = input(prompt);
if isempty(k), k = 1; end
if ~isnumeric(k)
    error('ERROR: Select number of singular triplets (k) from given values.'); 
end
if ~any(k == [1 2 3 4]) 
    error('ERROR: Select number of singular triplets (k) from given values.'); 
end   
mb_rd2 = k+1;
disp(sprintf(' '));
if example < 7 
   m1 = num2str(k+1); m2 = num2str(k+2); m3=num2str(k+3); m4=num2str(k+4);
   marray=[k+1  k+2 k+3 k+4];
else
   m1 = num2str(12); m2 = num2str(13); m3=num2str(14); m4=num2str(15);
   marray=[12  13 14 15];
end
str1 = 'Select basis size (m)'; 
str2 = ': ';
prompt = strcat(str1," ",m1,', ',m2,', ',m3,', ',m4," ",str2," ");
mb = input(prompt);
if isempty(mb) && example < 7,  mb = k+1; end
if isempty(mb) && example == 7, mb = 12; end
if ~isnumeric(mb)
    error('ERROR: Select  basis size (m) from given values.'); 
end
if ~any(mb == marray) 
    error('ERROR: Select basis size (m) from given values.'); 
end
prompt = 'Output statistics after each call of trrsvds? Y/N: ';
disp(sprintf(' '));
str = input(prompt,'s'); 
if isempty(str) || ~ischar(str), str = 'N'; end
str = upper(str); if ~strcmp(str,'Y'), str = 'N'; end
disp(sprintf(' '));
      
% *************************************************************************
% MATLAB programs used in driver program - http://www.math.uri.edu/~jbaglama/
% 1. rd2svds 
% 2. trrsvds
% Public codes and internal MATLAB codes used in the paper for comparision
% are *NOT* included in driver program 
% ****** NOT PART OF DRIVER *************
% eigs(C)       -internal MATLAB C = [0 A; A' 0].
% eigs(A^TA)    -internal MATLAB
% svds          -internal MATLAB
% irlba         - http://www.math.uri.edu/~jbaglama/
% svdifp        - https://github.com/wildstone/SVDIFP
% primme_svds   - https://github.com/primme/primme
% ************************************************************************* 

% Check for exists of MATLAB codes for driver program
if exist('trrsvds') ~= 2,  error('trrsvds is needed'); end
if exist('rd2svds') ~= 2,  error('rd2svds is needed'); end

% Initial values.
A=[]; m=[]; n=[];

% Sigma set to be 'LS' largest.
% Can be set to 'SS' smallest - only runs one example well1850 or
% diag(1:500) if fail to download example. 
% Not recommeded for smallest singular triplets. 
if example <= 6
  sigma = 'LS'; 
  % Set tolerance for convergence.
  tol = 1d-6;    
else
  sigma = 'SS';
  % Set tolerance for convergence.
  tol = 1d-10;     
end

% *************************************************************************
% Get the matrix to use for the examples.
% *************************************************************************
   
% Set value to determine if matrix cannot be use.
% If fail to load matrix - program defaults to diag(1:500).
fail = 0;

% 1. Diagonal matrix - used for examples 1, 2, and 3 in referenced paper
if example == 1
   A=sparse(diag(1:500)); name = 'diag(1:500)';
    disp(sprintf('Computing singular triplets for diagonal matrix diag(1:500)'));
end
  
% 2. illc1033  1033 x 320 - used for example 4 in referenced paper
if example == 2  
   name = 'illc1033'; name_mat = strcat(name,'.mat'); 
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/HB/illc1033.mat';
end
    
% 3. dwt1242  1,242 x 1,242 
if example == 3
   name = 'dwt1242'; name_mat = strcat(name,'.mat');
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/HB/dwt_1242.mat';
end

% 4 and 7. well1850  1,850 x 712 
if example == 4 || example == 7
   name = 'well1850'; name_mat = strcat(name,'.mat'); 
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/HB/well1850.mat';
end
   
% 5. pde2961 2961 x 2961  
if example == 5
   name ='pde2961'; name_mat = strcat(name,'.mat');
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/Bai/pde2961.mat';
end
      
% 6. amazon03002 262,111 x 262,111 - used for examples 1, 2, and 4 in referenced paper 
if example == 6
   name ='amazon0302'; name_mat = strcat(name,'.mat');
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/SNAP/amazon0302.mat';
end

% Check to see if matrices are in path - if not try to download from web.
% If not in path and fail to download switch to using A=sparse(diag(1:500));
if example ~= 1
   if exist(name_mat) ~= 2
      try
         str = strcat(name_mat," ",'not in path'); disp(sprintf(str));
         str = strcat('Try getting matrix'," ",name," ", 'from web ...');
         disp(sprintf(str));
         websave(name_mat,address); 
      catch
         str = strcat('FAILED getting'," ",name," ", 'from web and does not exist in path');
         disp(sprintf(str));
         disp(sprintf('Changing to diagonal matrix diag(1:500)'));
         A=sparse(diag(1:500)); name = 'diag(1:500)'; fail = 1;
      end       
   end
end

disp(sprintf(' '));
    
% Load matrix from cell and clear from memory to save space.
if example ~= 1 && ~fail, load(name_mat); A = Problem.A; clear Problem; end

% Check again matrix exists is in path.
if isempty(A)
   str = strcat('FAILED getting'," ",name," ", 'from web and does not exist in path');
   disp(sprintf(str));
   disp(sprintf('Changing to diagonal matrix diag(1:500)'));
   A=sparse(diag(1:500)); name = 'diag(1:500)';
end  
    
% Get size of the matrix A.
[m,n] = size(A);
    
% Set the common starting random vector for ALL routines.
rng(0,'twister');  p0 = randn(n,1); 
    
% Set the fixed common parameters for trrsvds 
OPTS_trrsvds       = [];     % Clear to remove all previous settings.
OPTS_trrsvds.tol   = tol;    % Set tolerance
OPTS_trrsvds.sigma = sigma;  % Set location based on sigma.
          
% Initialize the Data 4 x 7 matrix for output 
%  7 -  the number of methods that will be run for each example. 
%  Data(1,1:7)  = # mvp
%  Data(2,1:7)  = mvp timings 
%  Data(3,1:7)  = total cpu in sec. timings
%  Data(4,1:7)  = max error for that k and m value 
%  Data (:, j)  =  0  -> initialization
%  Data (:, j)  = -10 -> method j did not conv. in the max. iter
Data(1:4+k,1:7) = zeros(4+k,7);
           
% ***********************************************************
% 1. Call rd2svds ** ONLY USED WHEN m == k+1 and sigma = 'LS'
% ***********************************************************
           
dj = 1; % count used for data
           
if mb == k+1 && ~strcmp(sigma,'SS')  
              
   disp(sprintf(' rd2svds: ( k = %d ,m = %d )',k,mb));
           
   % Call rd2svds
   tStart = tic; mcount = 0; MV_tim=0;
   [U,S,V,STATS] = rd2svds(A,m,n,p0,k,tol,1d-1*tol); 
   if strcmp(str,'Y'), disp(STATS), end
              
   % Collect data
   Data(1,dj) = STATS.numMatProds; 
   Data(2,dj) = STATS.timeMatProds; 
   Data(3,dj) = STATS.timeTotal;
   Data(4,dj) = STATS.maxnormres;
   for i=1:k, Data(4+i,dj) = S(i,i); end
   if strcmp(STATS.convergedKVals,'FALSE') % Check if method failed to converge. 
     disp(sprintf('rd2svds FAILED TO CONV (code:-10)...'))
     Data(1:4+k,dj) = -10*ones(4+k,1);
   end
               
   % Clear variables.
   V=[]; S=[]; U=[]; STATS=[]; 
   
end
        
% ************************** 
% 2. Call  trrsvds(NOR,RITZ)
% **************************
disp(sprintf(' trrsvds(NOR,RITZ): ( k = %d ,m = %d )',k,mb));
dj = dj+1;  % count used for data
           
% Set the parameters.
OPTS = OPTS_trrsvds; % Reset parameters
OPTS.p0 = p0;        % set starting vector 
OPTS.method='NOR';   % set refined method to use normal equations 
OPTS.roh='RITZ';     % set thick-restarted to use Ritz
OPTS.k = k;          % number of singular triplets
OPTS.m_b = mb;       % fix size of Lanczos basis
           
% Call trrsvds
[U,S,V,STATS] = trrsvds(A,OPTS);
% Can also be called as [U,S,V,STATS] = trrsvds('ATprod',m,n,OPTS);
if strcmp(str,'Y'), disp(STATS), end
          
% Collect data
Data(1,dj) = STATS.numMatProds; 
Data(2,dj) = STATS.timeMatProds; 
Data(3,dj) = STATS.timeTotal;
Data(4,dj) = STATS.maxnormres;
for i=1:k, Data(4+i,dj) = S(i,i); end
if strcmp(STATS.convergedKVals,'FALSE') % Check if method failed to converge. 
   disp(sprintf('trrsvds(NOR,RITZ) FAILED TO CONV (code:-10)...'))
   Data(1:4+k,dj) = -10*ones(4+k,1);
end
                   
% Clear variables.
V=[]; S=[]; U=[]; STATS=[]; OPTS=[];
           
% ************************** 
% 3. Call  trrsvds(AUG,RITZ)
% **************************
disp(sprintf(' trrsvds(AUG,RITZ): ( k = %d ,m = %d )',k,mb));
dj = dj+1;  % count used for data
           
% Set the parameters.
OPTS = OPTS_trrsvds;  % Reset parameters
OPTS.p0 = p0;         % set starting vector. 
OPTS.method='AUG';    % set refined method to use augmented equations
OPTS.roh='RITZ';      % set thick-restarted to use Ritz 
OPTS.k = k;           % number of singular triplets
OPTS.m_b = mb;        % fix size of Lanczos basis
           
% Call trrsvds
[U,S,V,STATS] = trrsvds(A,OPTS);
% Can also be called as [U,S,V,STATS] = trrsvds('ATprod',m,n,OPTS);
if strcmp(str,'Y'), disp(STATS), end
           
% Collect data
Data(1,dj) = STATS.numMatProds; 
Data(2,dj) = STATS.timeMatProds; 
Data(3,dj) = STATS.timeTotal;
Data(4,dj) = STATS.maxnormres;
for i=1:k, Data(4+i,dj) = S(i,i); end 
if strcmp(STATS.convergedKVals,'FALSE') % Check if method failed to converge. 
   disp(sprintf('trrsvds(AUG,RITZ) FAILED TO CONV (code:-10)...'))
   Data(1:4+k,dj) = -10*ones(4+k,1);
end 
               
% Clear variables to save space and avoid mis-uses later on.
V=[]; S=[]; U=[]; STATS=[]; OPTS=[]; 
           
% *************************** 
% 4. Call  trrsvds(NOR,HARM)
% ***************************
           
disp(sprintf(' trrsvds(NOR,HARM): ( k = %d ,m = %d )',k,mb));
dj = dj+1;  % count used for data
           
% Set the parameters.
OPTS = OPTS_trrsvds;  % Reset parameters
OPTS.p0 = p0;         % set starting vector. 
OPTS.method='NOR';    % set refined method to use normal equations
OPTS.roh='HARM';      % set thick-restarted to use harmonic Ritz 
OPTS.k = k;           % number of singular triplets 
OPTS.m_b = mb;        % fix size of Lanczos basis
           
% Call trrsvds 
[U,S,V,STATS] = trrsvds(A,OPTS);
% Can also be called as [U,S,V,STATS] = trrsvds('ATprod',m,n,OPTS);
if strcmp(str,'Y'), disp(STATS), end
           
% Collect data
Data(1,dj) = STATS.numMatProds; 
Data(2,dj) = STATS.timeMatProds; 
Data(3,dj) = STATS.timeTotal;
Data(4,dj) = STATS.maxnormres;
for i=1:k, Data(4+i,dj) = S(i,i); end
if strcmp(STATS.convergedKVals,'FALSE') % Check if method failed to converge.
   disp(sprintf('trrsvds(NOR,HARM) FAILED TO CONV (code:-10)...'))
   Data(1:4+k,dj) = -10*ones(4+k,1);
end 
           
% Clear variables to save space and avoid mis-uses later on.
V=[]; S=[]; U=[]; STATS=[]; OPTS=[];
           
% *************************** 
% 5. Call  trrsvds(AUG,HARM)
% ***************************
           
disp(sprintf(' trrsvds(AUG,HARM): ( k = %d ,m = %d )',k,mb));
dj = dj+1;  % count used for data
           
% Set the parameters.
OPTS = OPTS_trrsvds;   % Reset parameters
OPTS.p0 = p0;          % set starting vector. 
OPTS.method='AUG';     % set refined method to use augmented equations
OPTS.roh='HARM';       % set thick-restarted to use harmonic Ritz 
OPTS.k = k;            % number of singular triplets 
OPTS.m_b = mb;         % fix size of Lanczos basis
           
% Call method trrsvds
[U,S,V,STATS] = trrsvds(A,OPTS);
% Can also be called as [U,S,V,STATS] = trrsvds('ATprod',m,n,OPTS);
if strcmp(str,'Y'), disp(STATS), end
        
% Collect data
Data(1,dj) = STATS.numMatProds; 
Data(2,dj) = STATS.timeMatProds; 
Data(3,dj) = STATS.timeTotal;
Data(4,dj) = STATS.maxnormres;
for i=1:k, Data(4+i,dj) = S(i,i); end   
if strcmp(STATS.convergedKVals,'FALSE')
   disp(sprintf('trrsvds(AUG,HARM) FAILED TO CONV (code:-10)...'))
   Data(1:4+k,dj) = -10*ones(4+k,1);
end
           
% Clear variables to save space and avoid mis-useslater on.
V=[]; S=[]; U=[]; STATS=[]; OPTS=[];
           
% *************************** 
% 6. Call  trrsvds(THK,RITZ)
% ***************************
           
disp(sprintf(' trrsvds(THK,RITZ): ( k = %d ,m = %d )',k,mb));
dj = dj+1;  % count used for data
           
% Set the parameters.
OPTS = OPTS_trrsvds; 
OPTS.p0 = p0; % set starting vector. 
OPTS.method='THK';     % set method to NOT use refined - only thick-restart
OPTS.roh='RITZ';       % set thick-restarted to use Ritz 
OPTS.k = k;            % number of singular triplets 
OPTS.m_b = mb;         % fix size of Lanczos basis
           
% Call method trrsvds
[U,S,V,STATS] = trrsvds(A,OPTS);
% Can also be called as [U,S,V,STATS] = trrsvds('ATprod',m,n,OPTS);
if strcmp(str,'Y'), disp(STATS), end
        
% Collect data
Data(1,dj) = STATS.numMatProds; 
Data(2,dj) = STATS.timeMatProds; 
Data(3,dj) = STATS.timeTotal;
Data(4,dj) = STATS.maxnormres;
for i=1:k, Data(4+i,dj) = S(i,i); end  
if strcmp(STATS.convergedKVals,'FALSE')
   disp(sprintf('trrsvds(THK,RITZ) FAILED TO CONV (code:-10)...'))
   Data(1:4+k,dj) = -10*ones(4+k,1);
end
           
% Clear variables to save space and avoid mis-useslater on.
V=[]; S=[]; U=[]; STATS=[]; OPTS=[];
           
% *************************** 
% 7. Call  trrsvds(THK,HARM)
% ***************************
           
disp(sprintf(' trrsvds(THK,HARM): ( k = %d ,m = %d )',k,mb));
dj = dj+1;  % count used for data
           
% Set the parameters.
OPTS = OPTS_trrsvds; 
OPTS.p0 = p0; % set starting vector. 
OPTS.method='THK';     % set method to NOT use refined - only thick-restart
OPTS.roh='HARM';       % set thick-restarted to use harmonic Ritz 
OPTS.k = k;            % number of singular triplets 
OPTS.m_b = mb;         % fix size of Lanczos basis
           
% Call method trrsvds 
[U,S,V,STATS] = trrsvds(A,OPTS);
% Can also be called as [U,S,V,STATS] = trrsvds('ATprod',m,n,OPTS);
if strcmp(str,'Y'), disp(STATS), end
        
% Collect data
Data(1,dj) = STATS.numMatProds; 
Data(2,dj) = STATS.timeMatProds; 
Data(3,dj) = STATS.timeTotal;
Data(4,dj) = STATS.maxnormres;
for i=1:k, Data(4+i,dj) = S(i,i); end 
if strcmp(STATS.convergedKVals,'FALSE')
   disp(sprintf('trrsvds(THK,HARM) FAILED TO CONV (code:-10)...'))
   Data(1:4+k,dj) = -10*ones(4+k,1);
end 
           
% Clear variables to save space and avoid mis-useslater on.
V=[]; S=[]; U=[]; STATS=[]; OPTS=[];
           
% Display example information. 
stats = ["MVPs"; "MVPs-CPU";"TOT-CPU";"MAX-NORM"];
for i=1:k, str = num2str(i); stats = [stats;strcat('S(',str,')')]; end
if mb == k+1 && ~strcmp(sigma,'SS')  
   VarNames = {name, 'rd2svds','trrsvds(NOR,RITZ)','trrsvds(AUG,RITZ)','trrsvds(NOR,HARM)','trrsvds(AUG,HARM)','trrsvds(THK,RITZ)','trrsvds(THK,HARM)'};
   T = table(stats(:,1),Data(:,1),Data(:,2),Data(:,3), Data(:,4),Data(:,5),Data(:,6),Data(:,7),'VariableNames',VarNames);
else
   VarNames = {name, 'trrsvds(NOR,RITZ)','trrsvds(AUG,RITZ)','trrsvds(NOR,HARM)','trrsvds(AUG,HARM)','trrsvds(THK,RITZ)','trrsvds(THK,HARM)'};
   T = table(stats(:,1),Data(:,2),Data(:,3), Data(:,4),Data(:,5),Data(:,6),Data(:,7),'VariableNames',VarNames);
end
disp(sprintf(' Size of matrix A: %d x %d',m,n));
disp(sprintf(' Number of values: k = %d',k));
disp(sprintf(' Tolerance: TOL = %0.5g',tol));
disp(sprintf(' Size of space: m = %d',mb));
disp(sprintf(' Location: sigma = %s',sigma)); 
disp(sprintf(' '));
disp(T);
disp(sprintf(' '));
disp(sprintf('**********************************************************************'));
 
% ******************************************************************
% Warning on using sigma = 'SS'
% ******************************************************************
if strcmp(sigma,'SS')
   disp(sprintf('*****************************************************'));
   disp(sprintf('trrsvds not recommended for smallest singular values.'));
   disp(sprintf('*****************************************************'));
end

end % driver_trrsvds_rd2svds function