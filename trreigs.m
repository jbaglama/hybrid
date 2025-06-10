function varargout = trreigs(varargin)

% TRREIGS: Computes the K extreme eigenvalue and associated eigenvector 
% of a N x N symmetric matrix A.
%
% PROGRAM INFORMATION:
% -------------------
%
%  ...  = TRREIGS(A)
%  ...  = TRREIGS('AFUN',N)
%  ...  = TRREIGS(@(x)A(x),N)
% 
% The first input value into IRREIGS can be a numeric matrix or a
% function/function handle that returns Y = A*X. If the input matrix 
% is a function/function handle the second input value must be 
% the size of the matrix A. The the input/output structure of M-file 
% Y = AFUN(X,N). 
%
% OUTPUT OPTIONS:
% ---------------
%
% I.)   TRREIGS(A)
%       If convergence, displays the K desired eigenvalues.    
%
% II.)  D = TRREIGS(A) 
%       If convergence, returns K eigenvalues in the vector D. 
%
% III.) [V,D] = TRREIGS(A)
%       If convergence, returns a diagonal matrix D that contains the
%       K eigenvalues along the diagonal and the matrix V contains the 
%       corresponding orthogonal eigenvectors, such that A*V = V*D. 
%       If TRREIGS reaches the maximum number of iterations before convergence 
%       then V = [] and D = []. Use output option IV.).
%
% IV.)  [V,D,FLAG] = TRREIGS(A) 
%       This option returns the same structural output as (III) plus a two dimensional array FLAG 
%       that reports if the algorithm converges and the number of matrix-vector 
%       products. If FLAG(1) = 0 then this implies normal return: all eigenvalues have 
%       converged. If FLAG(1) = 1 then the maximum number of iterations have been 
%       reached before all desired eigenvalues have converged. FLAG(2) contains the 
%       number of matrix-vector products used by the code. If the maximum number of 
%       iterations are reached (FLAG(1) = 1), then the matrices V and D contain 
%       the last available approximations for the eigenpairs.
%
% INPUT OPTIONS:
% --------------
%                                   
%       ... = TRREIGS(A,OPTS) or TRREIGS('AFUN',N,OPTS) or TRREIGS(@(x)A(x),N,OPTS)
%       OPTS is a structure containing input parameters. The input parameters can
%       be given in any order and can greatly influence convergence rates. The structure OPTS 
%       may contain some or all of the following input parameters. If parameter OPTS is missing 
%       or an input parameter in OPTS is not set, default value(s) are used. The string for the input parameters 
%       can contain upper or lower case characters. 
%       
%  INPUT PARAMETER      DESCRIPTION 
% 
%  OPTS.COEFF        When K > 1, COEFF is used to determine which matrix is
%                    use to compute the coefficients for the linear combination
%                    of iterative refined vectors. See Section 6 in the
%                    reference.
%                    coeff = 1 -> matrix (6.5)
%                    coeff = 2 -> Add one more term - matrix (6.9) (Default) 
%                    coeff = 3 -> Add all terms - matrix (6.10)
%                    DEFAULT VALUE COEFF = 2
%
%  OPTS.DELTA1       Used to determine when thick restarting with Ritz vectors can 
%                    start switching to restarting with iterative refined vectors.
%                    The value should be related to the overall convergence tolerance. 
%                    Therefore, we restrict the input to be in decimal form, 0<= DELTA1 <= 1
%                    and check for switching if the maximum residual norm
%                    of the desired Ritz pairs is <= TOL^(DELTA1). 
%                    EXAMPLE: If TOL = 1d-8 and user would like to start switching when convergence 
%                    has reached 25% of TOL or 1d-2, user inputs .25. The routine computes TOL^(.25). 
%                    Too large and only thick-restarted will be done. Too small and stagnation 
%                    results from using a poor Krylov space. 
%                    || A*V_RITZ - V_RITZ*D_RITZ ||_2 <= TOL^(DELTA1)*||A||_2. 
%                    ||A||_2 is approximated by largest absolute Ritz value.
%                    DEFAULT VALUE  DELTA1 = 0.1   
%
%  OPTS.DELTA2       Used to determine when the Iterative refined vectors
%                    are "close" enough to the Ritz vectors and can be  
%                    considered for restarting. Too large and only thick-restarted will be done.
%                    Too small and stagnation may result from using iterative refined vectors
%                    closer to the non-corresponding Ritz vector. See Example 5.2 in reference. 
%                    (all 1 <= j <= K) |(Refine vectors_j)^T*(Ritz vectors_j)| > DELTA2
%                    Must be 0< DELTA2 < 1
%                    DEFAULT VALUE  DELTA2 = 0.9 
%
%  OPTS.FLT          Toggle to restrict using restarting with iterative refined Ritz vectors when the 
%                    current computed iterative refined Ritz values are better than the past iterations 
%                    best Ritz approximation. This is only used when K=1 and
%                    should only be used when M is small. FLT = 1 uses
%                    the restriction, FLT= 0 do not use restriction. 
%                    DEFAULT VALUE FLT = 1 if M > 2 and FLT = 0 when M = 2.
%
%  OPTS.K            Number of desired eigenvalues. 
%                    DEFAULT VALUE  K = 1
%
%  OPTS.M            Number of Lanczos vectors, i.e. size tridiagonal
%                    matrix. Currently, full reorthogonalization is
%                    used. Large M will increase non-matrix-vector product
%                    CPU times.
%                    DEFAULT VALUE  M = 2  
%
%  OPTS.MAXIT        Maximum number of iterations, i.e. maximum number of Lanczos restarts.                           
%                    DEFAULT VALUE  MAXIT = 1000
%
%  OPTS.MAXITREF     Maximum number of iterations used to find iterative
%                    refined Ritz values, see Section 5 Algorithm 5.1 in
%                    reference. Recommend 10% of MAXIT. Stagnation prevention 
%                    is in place. Require MAXITREF > =100.
%                    DEFAULT VALUE  MAXITREF = 100
%
%  OPTS.SIGMA        Two letter string specifying the location of the desired eigenvalues.
%                     'LM' or 'LA' Largest Magnitude or Algebraic
%                             'SA' Smallest Algebraic                   
%                    DEFAULT VALUE   SIGMA = 'LM' 
%                                                          
%  OPTS.TOL          Tolerance used for convergence. Convergence is determined when             
%                    || A*V - V*D ||_2 <= TOL*||A||_2. ||A||_2 is approximated by
%                    largest absolute Ritz value. V and D are either iterative 
%                    refined or Ritz. If iterative refined V matrix is made orthogonal.
%                    See sections 2 and 7 in the reference. 
%                    DEFAULT VALUE    TOL = SQRT(EPS) (roughly 1d-8) 
%
%  OPTS.V0           A matrix of starting vectors.       
%                    DEFAULT VALUE  V0 = randn
%
%  DATE MODIFIED: 12/03/2020
%  VER:  1.0
%  
% AUTHORS: J. Baglama, T, Bella, and J. Picucci
%
% REFERENCE:
%  1. Baglama, J, Bella, T, Picucci, J, "An Iterative Method for Computing a Few 
%     Eigenpairs of a Large Sparse Symmetric Matrix" preprint 2020.
%
% *************************************************************************
% * This MATLAB code is provided to illustrate Algorithm 6.1 and is NOT   *
% * optimized for performance or set up for commercial use.               *
% * Any use beyond illustrate purposes (e.g. comparison for publications) * 
% * requires consent of the authors.                                      * 
% *************************************************************************

% Too many output arguments requested.
if (nargout >= 4), error('ERROR: Too many output arguments.'); end

%----------------------------%
% BEGIN: PARSE INPUT VALUES. %
%----------------------------%

% No input arguments, return help.
if nargin == 0, help trreigs, return, end

% Get matrix A. Check type (numeric or character) and dimensions.
if ischar(varargin{1}) || isa(varargin{1}, 'function_handle')
   sst_data = 3; 
   if nargin == 1, error('ERROR: Need N (size of matrix A).'); end  
   n = varargin{2};
   if ~isnumeric(n) || length(n) > 2 
      error('ERROR: Second argument N must be a numeric value.'); 
   end
else
   sst_data = 2;
   n = size(varargin{1},1);
   if any(size(varargin{1}) ~= n), error('ERROR: Matrix A is not square.'); end
   if ~isnumeric(varargin{1}), error('ERROR: A must be a numeric matrix.'); end
end   

% Square root of machine tolerance used in convergence testing.
sqrteps = sqrt(eps);  

% Set all input options to default values.
k=1; hybrid = 0; maxit = 1000; m=2;
maxitref=100; projt = 0; sigma = 'LM'; tol = sqrteps; 
delta2 = 0.9; delta1=0.1; coeff = 2; flt = 0; input_flt=[];

% Preallocate memory for large matrices.
v = spalloc(n,m,n*m); f = spalloc(n,1,n);

% Get input options from the data structure.
if nargin > 1 + ischar(varargin{1})
   options = varargin{sst_data:nargin};
   names = fieldnames(options);
   I = strmatch('COEFF',upper(names),'exact');
   if  ~isempty(I), coeff = getfield(options,names{I}); end
   I = strmatch('DELTA1',upper(names),'exact');
   if ~isempty(I), delta1 = getfield(options,names{I}); end
   I = strmatch('DELTA2',upper(names),'exact');
   if ~isempty(I), delta2 = getfield(options,names{I}); end
   I = strmatch('FLT',upper(names),'exact');
   if  ~isempty(I), input_flt = getfield(options,names{I}); end
   I = strmatch('K',upper(names),'exact');
   if  ~isempty(I), k = getfield(options,names{I}); end
   I = strmatch('M',upper(names),'exact');
   if  ~isempty(I), m = getfield(options,names{I}); end   
   I = strmatch('MAXIT',upper(names),'exact');
   if ~isempty(I), maxit = getfield(options,names{I}); end
   I = strmatch('MAXITREF',upper(names),'exact');
   if ~isempty(I), maxitref = getfield(options,names{I}); end
   I = strmatch('SIGMA',upper(names),'exact');
   if  ~isempty(I), sigma = upper(getfield(options,names{I})); end
   I = strmatch('TOL',upper(names),'exact');
   if ~isempty(I), tol = getfield(options,names{I}); end
   I = strmatch('V0',upper(names),'exact');
   if ~isempty(I), v = getfield(options,names{I}); end
end 

%***************************************************
% Check for some input errors in the data structure. 
% **** This is not an exhaustive check list. ******
%***************************************************

% Check that input values are numerical values.
if (~isnumeric(coeff)    || ~isnumeric(delta1) || ~isnumeric(delta2)  || ...
    ~isnumeric(k)        || ~isnumeric(m)      || ~isnumeric(maxit)  || ...
    ~isnumeric(maxitref) ||  ~isnumeric(tol))
   error('ERROR: Incorrect type for input value(s) in the structure.');
end

% Check value of COEFF
if coeff ~=1 && coeff ~=2 && coeff ~=3 , error('ERROR: COEFF must 1, 2, or 3'); end

% Check value of DELTA1
if delta1<0 || delta1 > 1, error('ERROR: 0<= DELTA1 <= 1'); end

% Check value of DELTA2
if delta2<=0 || delta2 >= 1, error('ERROR: 0< DELTA2 < 1'); end

% Check value of FLT.
if ~isempty(input_flt) % User input. 
  if input_flt ~=0 && input_flt ~=1, error('ERROR: FLT must 0 or 1'); end
  flt = input_flt; % set to user input.
end

% Check value of K.
if k <= 0,  error('ERROR: K must be a positive value.'); end
if k >= m,  error('ERROR: K is too large compared to M.'); end

% Check value of M.
if m  <= 1,  error('ERROR: M must be greater than 1.');  end
% Resize Krylov subspace if M (i.e. number of Lanczos vectors) is larger 
% than N (i.e. the size of the matrix A) then resize to <= N.
if m > n, m = n; end

% Check FLT again.
if k > 1, flt = 1; end     % overide input, cannot be used when k > 1.
if isempty(input_flt)      % no user input, check if m was changed
   if m > 2, flt = 1; end  % set to default
end    

% Check value of MAXIT
if maxit  <= 0,  error('ERROR: MAXIT must be a positive value.'); end

% Check value of MAXITREF
if maxitref  <= 0,   error('ERROR: MAXITREF must be a positive value.'); end
if maxitref > maxit, error('ERROR: MAXITREF should be less than MAXIT'); end 
if maxitref < 100,   error('ERROR: MAXITREF >= 100'); end 

% Check value of SIGMA. 
if isnumeric(sigma),   error('ERROR: SIGMA must be SA, LM or LA.'); end
if length(sigma) ~= 2, error('ERROR: SIGMA must be SA, LM or LA.'); end
if ~strcmp(sigma,'LM') && ~strcmp(sigma,'LA') &&  ~strcmp(sigma,'SM') && ~strcmp(sigma,'SA')
   error('ERROR: SIGMA must be LM, LA or SA.');
end
if strcmp(sigma, 'SM') 
    error('SM currently not supported. Recommend user use A\X and LM - output is 1/lambda'); 
end

% Check value of TOL.
% Set tolerance to machine precision if tol < eps.
if tol < eps, tol = eps; warning('Changing TOL to EPS');  end

% If starting vector V0 is not given then set starting vector V0 to be a 
% (n x 1) matrix of normally distributed random numbers.
if nnz(v) == 0, v = randn(n,1); end 

% Check starting vector V0.
if ~isnumeric(v), error('ERROR: Incorrect starting matrix V0.'); end
if ((size(v,1) ~= n) || (size(v,2) ~= 1))
    error('ERROR: Incorrect size of starting matrix V0.'); 
end

%--------------------------%
% END: PARSE INPUT VALUES. %
%--------------------------%

%-----------------------------------------------------------%
% BEGIN: DESCRIPTION AND INITIALIZATION OF LOCAL VARIABLES. %
%-----------------------------------------------------------%
% <- DO NOT MODIFY ->
dmax=zeros(m,1);      % Initialize array for max. abs. value of all Ritz values 
if strcmp(sigma,'SA'), dmax = Inf(m,1); end
if strcmp(sigma,'LA'), dmax = -Inf(m,1); end
f = zeros(n,1);       % Initialize residual vector f.
iter = 1;             % Main loop iteration count.
r_rf_0=zeros(m,1);    % Used for iterative refined comparsion from last iteration.
rztol=tol^(delta1);   % Tolerance for checking on using iterative refined
Smax = 0;             % Holds the maximum absolute value of all computed Ritz values.
T = zeros(m,m);       % Initialize T matrix.
thick=0;              % Indicate if thick restarted is to be used.
thk = 2;              % Indicate starting point for Lanczos.

%---------------------------------------------------------%
% END: DESCRIPTION AND INITIALIZATION OF LOCAL VARIABLES. %
%---------------------------------------------------------%

%----------------------------------------------%
% BEGIN: INITIALIZATION - ONE STEP OF LANCZOS. %
%----------------------------------------------%
v(:,1) = v(:,1)/norm(v(:,1));  
v(:,2) = matrixprod(varargin{1},v(:,1),n); mprod=1;
T(1,1)  = v(:,2)'*v(:,1); 
v(:,2) = v(:,2) - v(:,1)*T(1,1);
dotv = v(:,2)'*v(:,1); % Reorth step
v(:,2) = v(:,2) - v(:,1)*dotv;
T(1,1) = T(1,1) + dotv; 
T(2,1) = norm(v(:,2)); T(1,2) = T(2,1);
if abs(T(2,1)) < eps*abs(T(1,1))
   error('ERROR: VO caused |T(2,1)| < eps|T(1,1))-> Change VO.');
end
v(:,2) = v(:,2)/T(2,1); 
%--------------------------------------------%
% END: INITIALIZATION - ONE STEP OF LANCZOS. %
%--------------------------------------------%

%---------------------------%
% BEGIN: ITERATION PROCESS. %
%---------------------------%
while (iter <= maxit)
     
  %-------------------------%
  % BEGIN: LANCZOS PROCESS. %
  %-------------------------% 
  for j=thk:m
    f = matrixprod(varargin{1},v(:,j),n); mprod=mprod+1;
    if thick == 1
       f = f - v(:,1:j-1)*T(j,1:j-1)';
       thick = 0;
    else
        f = f - v(:,j-1)*T(j,j-1);
    end  
    T(j,j) = f'*v(:,j);
    f = f - (v(:,j)*T(j,j));
    dotv = (f'*v(:,1:j))'; % Reorth step
    f = f - (v(:,1:j)*dotv);
    if norm(dotv)>sqrteps        % 2nd reorth step if needed
       dotv2 = (f'*v(:,1:j))';
       dotv=dotv+dotv2; 
       f = f - (v(:,1:j)*dotv2);
    end 
    for i=1:j, T(j,j) = T(j,j) + dotv(i); end
    fnorm = norm(f);
    if fnorm < max(Smax*tol,eps) && j>= k && iter > 1
        T = T(1:j,1:j); m=j; % resize T for exit;
        warning('Early termination in Lanczos - results may not be trusted');
        break; 
    end
    if fnorm < max(Smax*eps,eps) && (j < k || iter == 1)
        error('ERROR: Lanczos breakdown.'); 
    end
    if j < m 
       v(:,j+1) = f/fnorm;  
       T(j+1,j) = fnorm; 
       T(j,j+1) = T(j+1,j); 
    end
  end
  %-----------------------%
  % END: LANCZOS PROCESS. %
  %-----------------------%
  
  %------------------------------------------------------------------------%
  % BEGIN: COMPUTE/SORT EIGENVALS AND EIGENVECS OF T AND TEST CONVERGENCE. %
  %------------------------------------------------------------------------%
  [x_rz, d_rz] = eig(T); d_rz = diag(d_rz);
  Smax = max([Smax;abs(d_rz)]); % Use to app. ||A||_2
  if strcmp(sigma,'LM')
     [~,I_eig] = sort(abs(d_rz)); 
     I_eig = I_eig(length(I_eig):-1:1);
     x_rz = x_rz(:,I_eig); d_rz = d_rz(I_eig);
     for j=1:k
        dmax(j) = max(dmax(j),abs(d_rz(j)));
        dmax(j) = sign(d_rz(j))*dmax(j);
     end
  elseif strcmp(sigma,'LA')
     [~,I_eig] = sort(d_rz); 
     I_eig = I_eig(length(I_eig):-1:1);
     x_rz = x_rz(:,I_eig); d_rz = d_rz(I_eig);
     for j=1:k
        dmax(j) = max(dmax(j),d_rz(j));
     end
   elseif strcmp(sigma,'SA')
     [~,I_eig] = sort(d_rz); 
     x_rz = x_rz(:,I_eig); d_rz = d_rz(I_eig);
      for j=1:k
        dmax(j) = min(dmax(j),d_rz(j));
     end
  end 
  
  % Convergence check.
  nc = 0; % number of convergence desired Ritz values.
  res_r = ones(k,1); % reset residual norms. 
  for j=1:k 
      res_r(j) = abs(x_rz(m,j))*fnorm; 
      if res_r(j) < tol*Smax, nc = nc+1; end
  end
  if nc ==k || iter == maxit
     v=v*x_rz(:,1:k); d = diag(d_rz(1:k)); 
     break; 
  end
  
  % Compute size matrix T.
  Tsz = size(T,1); 
   
  %----------------------------------------------------------------------%
  % END: COMPUTE/SORT EIGENVALS AND EIGENVECS OF T AND TEST CONVERGENCE. %
  %----------------------------------------------------------------------% 
  
  %--------------------------------------------------------------------%
  % BEGIN: COMPUTE ITERATIVE REFINED AND DETERMINE RESTARTING VECTORS. %
  %--------------------------------------------------------------------%
  
  % Compute Iterative refined Ritz values.
  rconv=0; 
  if max(res_r(1)) < rztol*Smax
     [r_rf,x_rf,s_rf,u_rf,rconv] = refined_ritz(dmax(1:k),fnorm,T,maxitref,Tsz,k);
  end
  
  % rconv returns 0 or k. 0 indicates not all converged. rconv = k indicates all k converged 
  % Compute inner product between the k converged iterative refined Ritz vector(s) 
  % and Ritz vector(s). Also, check if iterative Ritz values are less than
  % the previous Ritz values. If so, do not use. 
  ang_rz_rf = 0; decr = 0;
  for i=1:rconv
      ang_rz_rf(i) = abs(x_rz(:,i)'*x_rf(:,i));  % computer the innner prod. to check angle
      if flt % if true - requires iter. refined to be better approx. prev. eig. of T.
         if strcmp(sigma,'LM') && abs(r_rf(i)) < abs(r_rf_0(i)), decr = 1; end
         if strcmp(sigma,'LA') && r_rf(i) < r_rf_0(i), decr = 1;  end
         if strcmp(sigma,'SA') && r_rf(i) > r_rf_0(i), decr = 1; end
      end
      r_rf_0(i) = dmax(i); 
  end
   
  % Find all inner products between iterative refined and Ritz > delta2
  I = find(ang_rz_rf > delta2);
  
  % Convergence check of iterative refined Ritz vector(s)
  if length(I) == k
      norm_res=[]; % reset residual norms. 
      [q_rf,~] = qr(x_rf,0); % QR to ensure orthogonal vectors
      nr = 0;                % number of iter. ref. converged.
      for j=1:k 
        Tq = T*q_rf(:,j) - r_rf(j)*q_rf(:,j);
        norm_res(j) = sqrt(Tq'*Tq + fnorm^2*q_rf(m,j)^2); 
        if norm_res(j) < tol*Smax, nr = nr+1; end
      end
      if nr ==k
         v=v*q_rf(:,1:k); d = diag(r_rf(1:k));
         break; 
      end
  end 
  
  % Determine if restart with linear combination of Iterative refined Ritz
  % or thick-restarted with Ritz vectors
   if (length(I) == k &&  max(res_r(1)) < rztol*Smax && decr == 0) 
      
      %------------------------------------------%
      % BEGIN: RESTART ITERATIVE REFINED SECTION %
      %------------------------------------------% 
      % Re-compute residual vector of iterative refined Ritz. Done without QR and is
      % different than convergence check since x_rf are not orthogonal.
      % Needed for restart. res_rf is m x k "small" vector.
      res_rf=[]; norm_res = []; % reset residual norms.
      for j=1:k
          res_rf(:,j) = s_rf(j)*(u_rf(1:m,j) - x_rf(:,j)*(x_rf(:,j)'*u_rf(1:m,j)));
          norm_res(j) = sqrt(res_rf(:,j)'*res_rf(:,j) + s_rf(j)^2*u_rf(m+1,j)^2);
      end
      
      % Initialize the k coefficients c_rf for the linear combination of
      % the k iterative refined Ritz vectors x_rf.
      c_rf = ones(length(r_rf),1);
      if k > 1
         Bc = zeros(k-1,k); Tem = T(m,1:m); em = zeros(m,1); em(m) = 1;
         
         %  coeff matrix - (6.5)
         if coeff == 1
            Dc = diag(x_rf(m,1:k));
            for i=1:k-1
               for j=1:k
                  Bc(i,j) = r_rf(j)^(i-1);
               end 
            end
           Bc = Bc*Dc;
         end
         
         % Added one more term coeff matrix - (6.9)
         if coeff == 2
            Bc(1,1:k) = x_rf(m,1:k);
            for i=2:k-1
               for j=1:k
                  Bc(i,j) = Tem*x_rf(:,j)*r_rf(j)^(i-2);
               end 
            end
         end
         
         % Use all terms coeff matrix - (6.10)
         if coeff == 3
             Bc(1,1:k) = x_rf(m,1:k);  % if k=2 => 1 x 2 matrix end
             if k >=3                  % if k=3 => 2 x 3 matrix end
               Bc(2,1:k) = Tem*x_rf(1:m,1:k);
             end
             if k >= 4
               for i=3:k-1
                for j=1:k
                  Bc(i,j) = Tem*x_rf(:,j)*r_rf(j)^(i-2);
                  r_sum = 0;
                  for jj = 3:i
                      r_sum =r_sum + r_rf(j)^(i-jj)*em'*T^(jj-2)*u_rf(1:m,j);
                  end
                    Bc(i,j) = Bc(i,j) + s_rf(j)*r_sum;    
                end
               end
            end
         end
         % compute the solution of the k-1 x k homogeneous system using
         % null. Following code is needed to avoid zero column(s) in Bc
         % for numerically converged eigenvectors. 
         Iin = []; Iout=[]; % set variables for which columns are used.
         Bc_max = max(abs(Bc),[],'all');
            for i=1:k  % search for columns numerically zero. 
              if norm(Bc(:,i),Inf) < sqrt(eps)*Bc_max
                Iout = [Iout i]; % references which columns to remove.
              else
                Iin = [Iin i];
              end
            end  
            if ~isempty(Iout) % remove numerically zero columns from Bc
               Bc = Bc(1:k-length(Iout)-1,Iin);
               c_rf(1:k) = zeros(k,1);
               c_rf_in = null(Bc); % call Matlab's null to solve system
               c_rf(Iin) = c_rf_in(:,1);
               c_rf(Iout) = norm_res(Iout); % Place zero entries back in sol.
            else
               c_rf = null(Bc); % call Matlab's null to solve full system
            end   
         c_rf=c_rf(:,1); % matlab's null can return more than one sol.
      end  
      
      %----------------------------------------------------------%
      % BEGIN: ONE STEP OF LANCZOS WITHOUT MATRIX VECTOR PRODUCT %
      % (See section 6 equations (6.2) & (6.12) in reference [1] %
      %----------------------------------------------------------%
      f_rf = 0;
      % Used for computing the coefficients the residual vectors  
      for j=1:k
         f_rf = f_rf + c_rf(j)*x_rf(m,j);
      end
      % Compute left and right vectors to avoid matrix-vector product
      % v_left & v_right are  m x 1 "small" vectors. 
      v_left = x_rf*c_rf; v_right = x_rf*(r_rf.*c_rf)+ res_rf*c_rf;
      v_left_norm = norm(v_left);
      v_left = v_left/v_left_norm; 
      v_right = v_right/v_left_norm;
      % Compute the matrix T
      T = zeros(2,2);
      % Update Lanczos vectors v(:,1) & v(:,2). f is reused and avoids extra
      % storage needs.
      f = v*v_right + (f_rf/v_left_norm)*f;
      v(:,1) = v*v_left; 
      T(1,1) = f'*v(:,1);
      f = f - (v(:,1)*T(1,1));
      dotv = (f'*v(:,1))';  % Reorth step
      f = f - (v(:,1)*dotv);
      T(1,1) = T(1,1) + dotv;
      fnorm = norm(f);
      v(:,2) = f/fnorm; 
      T(2,1) = fnorm; 
      T(1,2) = T(2,1);
      thk=2; thick=0; % set values to start Lanczos at step 2
      %--------------------------------------------------------%
      % END: ONE STEP OF LANCZOS WITHOUT MATRIX VECTOR PRODUCT %
      %--------------------------------------------------------%
      %----------------------------------------%
      % END: RESTART ITERATIVE REFINED SECTION %
      %----------------------------------------% 
   else
      %-----------------------------------%
      % BEGIN: THICK-RESTART RITZ SECTION %
      %-----------------------------------%
    
      % Simple strategy to improve convergence.
      k_thk = max(floor(nc + (Tsz-nc)/2),k);
      if Tsz - 1 - k_thk < 0, k_thk = Tsz-1; end
    
      % Set up matrices and vectors for thick-restarted 
      v = [v*x_rz(:,1:k_thk) f/fnorm];
      T = zeros(k_thk+1,k_thk+1);
      T(1:k_thk+1,1:k_thk) = [diag(d_rz(1:k_thk)); x_rz(m,1:k_thk)*fnorm];
      T = tril(T,-1) + tril(T)'; thk = k_thk+1; thick=1;
      %---------------------------------%
      % END: THICK-RESTART RITZ SECTION %
      %---------------------------------%
  end 
  %------------------------------------------------------------------%
  % END: COMPUTE ITERATIVE REFINED AND DETERMINE RESTARTING VECTORS. %
  %------------------------------------------------------------------%
 
 % Update the main iteration loop count. 
 iter = iter+1;

end % end main loop
%-------------------------%
% END: ITERATION PROCESS. %
%-------------------------%

%-----------------------%
% BEGIN: OUTPUT RESULTS %
%-----------------------%

% Output option I: Display eigenvalues only.
if (nargout == 0) 
    if iter < maxit
       eigenvalues = diag(d)
    else
       eigenvalues =[]    
    end
end
% Output option II: Set eigenvalues equal to output vector.  
if (nargout == 1)
    if iter < maxit
       varargout{1} = diag(d); 
    else
        varargout{1} =[];
    end
end

% Output option III: Output diagonal matrix of eigenvalues and
% corresponding matrix of eigenvectors.
if (nargout == 2)
    if iter < maxit
       varargout{1} = v;
       varargout{2} = d;
    else
      varargout{1} =[];
      varargout{2} =[];
    end  
end

% Output option IV: Output diagonal matrix of eigenvalues and
% corresponding matrix of eigenvectors and FLAG.
FLAG(1) = 0; FLAG(2) = mprod;
if iter >= maxit,FLAG(1) = 1; end 
if nargout == 3
    varargout{1} = v;
    varargout{2} = d;
    varargout{3} = FLAG; 
end

%---------------------%
% END: OUTPUT RESULTS %
%---------------------%

%---------------------------------------%
% BEGIN: COMPUTE REFINED RITZ ITERATION %
%---------------------------------------%

function [rho,v_min,s_min,u_min,rconv] = refined_ritz(D_ritz,R,T,maxitref,Tsz,k)
% Computes the Iterative refined Ritz values and vectors. Also, computes any needed vectos for computing residuals.  
%
% Input:
%   D_ritz  - (K x 1) vector of eigenvalues of T (aka Ritz values).
%        R  - real number - norm of residual vector F.     
%        T  - (TSZ x TSZ)  Lanczos tridiagonal matrix. 
% MAXITREF  - Integer indicating the maximum number of iterations for  the iterative Refined Ritz values. 
%     TSZ  -  Integer indicates the size of the  tridiagonal matrix.
%
% Output:
%     RHO  - (K x 1) vector of iterative refined Ritz values.
%   V_MIN  - (TSZ x K) Matrix of right singular values of [T; 0 R] and iterative refined Ritz vectors.
%   S_MIN  - (K x 1) vector of minumum singular values of [T; 0 R] associated with V_MIN and U_MIN. 
%            Values are needed in computing residuals.
%   U_MIN  - (TSZ+1 x K) Matrix of left singular values of [T; 0 R]
%            Values are needed in computing residuals.
%   RCONV  - Integer indicate if all K iterative refined Ritz converge.
%
%  Algorithm 5.1 in reference [1].
%
%  DATE MODIFIED: 06/03/2020
%  VER:  1.0

% Set up the augmented matrix [T; 0 R]
T_aug = zeros(Tsz+1,Tsz); T_aug(1:Tsz,1:Tsz) = T;
T_aug(Tsz+1:Tsz+1,Tsz) = R;

% Initialize values.
s_min = zeros(k,1);           % Intialize min. singular values of [T; 0 R]. 
v_min = zeros(Tsz,k);         % Intialize right singular values of [T; 0 R]. 
u_min = zeros(Tsz+1,k);       % Intialize left singular values of [T; 0 R].
rconv = zeros(k,1);           % Intialize convergence.
rho = D_ritz;                 % Initialize rho value.
rho_0 = rho;                  % Used to check for convergence during iterations.
sqrteps = sqrt(eps);          % square root of machine precision - eps 

% Compute k number of iterative refined Ritz values/vectors
for j = 1:k

    %  Set difference in rho values to test for stagnataion.
    diff_rho_0 = -1; v_min_0 = zeros(Tsz,1);
   
    % Iteratation to compute the iterative refined Ritz values/vectors
    for i=1:maxitref 

      % Compute the SVD of  [T; 0 R] - rho* I
      [U,S,V] = svd((T_aug-rho(j)*eye(Tsz+1,Tsz)),0);
    
      % Need the smallest singular triplet of [T; 0 R] - rho* I. Matlab
      % returns order of singular values largest to smallest.
      s_min(j) = S(Tsz,Tsz); v_min(:,j) = V(:,Tsz); u_min(:,j) = U(:,Tsz);
      
      % Compute the new rho = v_min'*T*v_min (aka refined Ritz value)
      rho(j) = v_min(:,j)'*T*v_min(:,j);
   
      % Compute the difference of previous rho to check for convergence 
      diff_rho = abs( (rho_0(j) - rho(j))/rho(j) );
  
      % Check for convergence
      if (diff_rho < eps && abs(v_min(:,j)'*u_min(1:Tsz,j)) < sqrteps)...
         || s_min(j) < eps || norm(v_min_0 - v_min) < eps 
            rconv(j) = 1;  break; 
      end
    
      % Check for stagnation to avoid too many unecessary iterations. Care 
      % must be taken to avoid the situation where rho(k) is still changing 
      % very slightly at first and then an increase in change later. An
      % early termination due to stagnation with no converegnce 
      % may avoid increase in change later that converges. 
      if abs(diff_rho_0 - diff_rho)/diff_rho < eps,  break; end  % stagnate
      if i>= 10 && mod(i,10)==0 % update every 10 iterations to avoid early termination. 
         diff_rho_0 = diff_rho; 
      end 
      
     % Update rho_0 and v_min_0 to test convergence. 
     rho_0(j) = rho(j); v_min_0 = v_min;
    end 
end

% After K iteration finish check to see if *all* K iterative Refined Ritz values
% have all converged. Reset rconv to return an integer value, rconv = k all
% converged and rconv = 0, not all converged.
if all(rconv == 1), rconv = k; else, rconv = 0; end

%-------------------------------------%
% END: COMPUTE REFINED RITZ ITERATION %
%-------------------------------------%

%-------------------------------%
% BEGIN: MATRIX-VECTOR PRODUCT. %
%-------------------------------%
function x = matrixprod(A,x,n)
% Computes the matrix vector products.
%
% Input:
%       A  - Matrix A.
%       x  - (N x 1) vector to multiply OP (operator).
%       N  - size of A.
% Output:
%       x  - (N x 1) Product of OP*X (operator).
%
if ischar(A)
    x=feval(A,x,n); % does not accept input parameters
elseif (isa(A, 'function_handle'))
    x=A(x);        % simple function handle 
else
    x = A*x;
end 
%-----------------------------%
% END: MATRIX-VECTOR PRODUCT. %
%-----------------------------%
