function varargout = trrsvds(varargin)

% TRRSVDS: Computes the k largest or smallest singular values and associated 
%          singular vectors of a m x n matrix A such that A*V = U*S and 
%          A'*U = V*S, V'*V=I and U'*U = I and S is a diagonal matrix.
%
% PROGRAM INFORMATION:
% -------------------
%
%  ...  = TRRSVDS(A)
%  ...  = TRRSVDS('AFUN',m,n)
% 
% The first input argument into TRRSVDS can be a numeric matrix A or 
% a function. If the m x n matrix A is a function, 'Afunc', 
% then the structure must be y = Afunc(x,m,n,'transpose'). If transpose = 'F', 
% then y = A*x. If transpose = 'T', then y = A'*x.
%
% OUTPUT OPTIONS:
% ---------------
%
% I.)   TRRSVDS(A)
%       If convergence, displays the k desired singular values.    
%
% II.)  S = TRRSVDS(A) 
%       If convergence, returns k singular values in the vector S. 
%
% III.) [U,S,V] = TRRSVDS(A)
%       If convergence, S is a diagonal matrix that contains the k desired 
%       singular values in descending order along the diagonal, the matrix 
%       V contains the corresponding "right" singular vectors, and U contains 
%       the corresponding "left" singular vectors such that that A*V = U*S, 
%       A'*U = V*S, V'*V = I, and U'*U = I. If TRRSVDS reaches the maximum 
%       number of iterations before convergence then U=[], S = [] and V = []. 
%
% IV.)  [U,S,V,STATS] = TRRSVDS(A) 
%       This option returns the same as (III) plus approximation of U,S,V if
%       no convergence and a structure STATS that reports statistical
%       information:
%
%       STATS =
%       numMatProds:  -> number of matrix-vector products with A and A^T
%       timeMatProds: -> total time computing products with A and A^T 
%       numIterRefRestart: -> number of restarts with Iterative Refined vectors
%       timeCompIterRef:   -> time computing iterative Refined vectors
%                               (independent of using them to restart)
%       timeReorth: -> time spent on full reorthogonalization 
%       timeTotal:  -> total time elasped
%       estimateSVmax: -> estimate of the maximum singular value over all iterations
%       estimateSVmin: -> estimate of the maximum singular value over all iterations
%       convergedKVals: -> true if all k singular triplets converged, otherwise false
%       outputestRitzorIterRef: -> if converged, names which values were used to exit
%                                  Ritz, Harmonic Ritz, Refined Ritz, or Iterative Refined Ritz 
%       maxnormres -> maximum residual of values for exit
%
%       If the maximum number of iterations are reached before convergence of all k desired 
%       singular values, convergedKVals = 'False' then the matrices U, V, and S contain any 
%       singular triplets that have converged plus the last singular triplets 
%       approximation for the pairs that have not converged.
%
% INPUT OPTIONS:
% --------------
%                                   
%       ... = TRRSVDS(A,OPTS) or TRRSVDS('AFUN',m,n,OPTS)
%       OPTS is a structure containing input parameters. The input parameters can
%       be given in any order and can greatly influence convergence rates. The structure OPTS 
%       may contain some or all of the following input parameters. If parameter OPTS is missing 
%       or an input parameter in OPTS is not set, default value(s) are used. The string for the input parameters 
%       can contain upper or lower case characters. 
%       
%  INPUT PARAMETER      DESCRIPTION 
%
%  OPTS.RoH          Four letter string ('RITZ' or 'HARM') specifying the use of either Ritz vectors or 
%                    harmonic Ritz vectors for thick-restarting.
%                    DEFAULT VALUE    RoH = 'HARM'  if SIGMA = 'SS'
%                                     RoH = 'RITZ'  if SIGMA = 'LS' or if cond(B) > 1/sqrt(eps)
%
%  OPTS.K            Number of desired singular values. 
%                    DEFAULT VALUE  K = 1
%
%  OPTS.M_B          Number of Lanczos vectors, i.e. size bidiagonal Lanczos
%                    matrix. Full reorthogonalization is used. Large M_B will 
%                    increase non-matrix-vector product CPU times.
%                    DEFAULT VALUE  M_B = 2   if SIGMA = 'LS' 
%                    DEFAULT VALUE  M_B = 15  if SIGMA = 'SS'
%
%  OPTS.MAXIT        Maximum number of iterations, i.e. maximum number of restarts.                           
%                    DEFAULT VALUE  MAXIT = 2000
%
%  OPTS.MAXITREF     Maximum number of iterations used to find iterative
%                    refined Ritz singular values. 
%                    DEFAULT VALUE  MAXITREF = 100
%
%  OPTS.METHOD       Three letter string ('NOR', 'AUG', or 'THK') to determine which method to use.
%                    Hybrid method - 'NOR' or 'AUG'. Determines how iterative refined Ritz values/vectors 
%                    are computed in the hybrid method. 
%                      NOR - (Hybrid) Lanczos Bidiagonal decomposition of size m
%                            AP_m = Q_mB_m  and A'Q_m = P_m B_m' + f_me_m'
%                            Equivalent eigenvalue system 
%                            A'AP_m = P_m B_m'B_m + B_m(m,m)f_me_m'
%                            Compute iterative refined Ritz on [B'*B; 0 B_m(m,m)*norm(f)]
%                      AUG - (Hybrid) Equivalent eigenvalue system                    
%                            [0 A; A' 0][Q_m 0; 0 P_m] = [Q_m 0; 0 P_m][0 B; B' 0]+[ 0 0; f_me_m' 0] 
%                            Compute iterative refined Ritz on [0 B; B' 0; norm(f) 0]
%                      THK - (non-Hybrid) no iterative refined Ritz values are computed. Thick- 
%                            restarted only with either Ritz vectors or harmonic Ritz vectors 
%                            depending on parameter RoH value.
%                    DEFAULT VALUE AUG = 'NOR'
%
%  OPTS.REORTH       Three letter string ('ONE' or 'TWO') specifying whether to use one-sided 
%                    ('ONE') full reorthogonalization or two-sided ('TWO'). One-sided is performed 
%                    only on the "short" vectors. Two-sided orthogonality is always used when 
%                    cond(A) estimated by cond(B) > 1/sqrt(eps).
%                    DEFAULT VALUE  REORTH = 'ONE' 
%
%  OPTS.SIGMA        Two letter string ('LS' or 'SS') specifying the location of the desired singular values.
%                    'LS'  Largest singular values and 'SS'  Smallest singular values.               
%                    DEFAULT VALUE   SIGMA = 'LS' 
%                                                          
%  OPTS.TOL          Tolerance used for convergence. Convergence is determined when             
%                    MAX (SQRT(|| AV - US||^2 + || A'U - VS ||^2))<= TOL*||A||. 
%                    Norm is the two norm and [U,S,V] are approximated either by Ritz, harmonic Ritz, Refined
%                    Ritz, or Iterative Refined Ritz singular triplets. ||A|| is approximated by largest singular 
%                    value of all projection matrices.
%                    DEFAULT VALUE    TOL = SQRT(EPS) (roughly 1d-8) 
%
%  OPTS.P0           An n x 1 starting vector if m >= n and sigma = 'LS' or 'SS' and an m x 1 starting vector 
%                    if m < n and sigma = 'SS'. P0 is not explicitly created and use first column of matrix V.  
%                    DEFAULT VALUE  P0 = randn(n,1)
%
%  DATE MODIFIED: 12/21/21
%  VER:  1.2
%  
% AUTHORS: 
% James Baglama     email: jbaglama@uri.edu
% Vasilije Perovic  email: perovic@uri.edu
% Jennifer Picucci  email: jenniferpicucci@uri.edu
%
% REFERENCES:
% 1. Baglama, J, Perovic, V, and Picucci, J, "Hybrid Iterative Refined Restarted 
%    Lanczos Bidiagonalization Method",  2021 submitted Numerical
%    Algorithms, preprint: http://www.math.uri.edu/~jbaglama/paper34.pdf
% 2. Baglama, J, Bella, T, and Picucci, J, "Hybrid Iterative Refined Method for 
%    Computing a Few Extreme Eigenpairs of a Symmetric Matrix", SIAM J. Sci. 
%    Comput. Special session on iterative methods, May 4, 2021. 
%    https://doi.org/10.1137/20M1344834 

% Start timing count using the MATLAB tic command.
tStart = tic; 

% Incorrect number of output arguments requested.
if (nargout > 4 || nargout==2 ), error('ERROR: Incorrect number of output arguments.'); end

%----------------------------%
% BEGIN: PARSE INPUT VALUES. %
%----------------------------%

% No input arguments, return help.
if nargin == 0, help trrsvds, return, end


% Matrix A is stored in varargin{1}. Check type (numeric or character) and dimensions.
if (isstruct(varargin{1})), error('A must be a matrix.'), end
if ischar(varargin{1})
   if nargin == 1, error('Need dimension M for matrix A).'), end  
   m = varargin{2};
   if ~isnumeric(m) || length(m) ~= 1 
      error('Second argument M must be a numeric value.'); 
   end
   if nargin == 2, error('Need dimension N for matrix A).'); end  
   n = varargin{3}; 
   if ~isnumeric(n) || length(n) ~= 1 
      error('Third argument N must be a numeric value.'); 
   end
else
   if ~isnumeric(varargin{1}), error('ERROR: A must be a numeric matrix.'); end
   [m,n] = size(varargin{1});
end

% Square root of machine tolerance used in convergence testing.
sqrteps = sqrt(eps);  

% Set all input options to default values.
k=1; maxit = 2000; m_b=[]; reorth='ONE'; method='NOR';
maxitref=100; sigma = 'LS'; tol = sqrteps; V=[]; roh = []; 

% Get input options from the data structure.
if nargin > 1 + 2*ischar(varargin{1})
   options = varargin{2+2*ischar(varargin{1}):nargin};
   names = fieldnames(options);
   I = strmatch('ROH',upper(names),'exact');
   if  ~isempty(I), roh = upper(getfield(options,names{I})); end
   I = strmatch('K',upper(names),'exact');
   if  ~isempty(I), k = getfield(options,names{I}); end
   I = strmatch('M_B',upper(names),'exact');
   if  ~isempty(I), m_b = getfield(options,names{I}); end   
   I = strmatch('MAXIT',upper(names),'exact');
   if ~isempty(I), maxit = getfield(options,names{I}); end
   I = strmatch('MAXITREF',upper(names),'exact');
   if ~isempty(I), maxitref = getfield(options,names{I}); end
   I = strmatch('METHOD',upper(names),'exact');
   if  ~isempty(I), method = upper(getfield(options,names{I})); end
   I = strmatch('REORTH',upper(names),'exact');
   if ~isempty(I), reorth = upper(getfield(options,names{I})); end
   I = strmatch('SIGMA',upper(names),'exact');
   if  ~isempty(I), sigma = upper(getfield(options,names{I})); end
   I = strmatch('TOL',upper(names),'exact');
   if ~isempty(I), tol = getfield(options,names{I}); end
   I = strmatch('P0',upper(names),'exact');
   if ~isempty(I), V = getfield(options,names{I}); end
end 

%***************************************************
% Check for some input errors in the data structure. 
% **** This is not an exhaustive check list. ******
%***************************************************

% Check that input values are numerical or char values.
if (~isnumeric(k) || ~isnumeric(m_b) || ~isnumeric(maxit) || ...
    ~isnumeric(maxitref) || ~ischar(method)|| ~ischar(reorth) || ...
    ~ischar(sigma)  || ~isnumeric(tol))
       error('ERROR: Incorrect type for input value(s) in the structure.');
end

% Check value of MAXIT
if maxit  <= 0,  error('ERROR: MAXIT must be a positive value.'); end

% Check value of MAXITREF
if maxitref  <= 0,       error('ERROR: MAXITREF must be a positive value.'); end
if maxitref > maxit,     error('ERROR: MAXITREF should be less than MAXIT'); end
if maxitref == 1
     warning(['WARNING: MAXITREF = 1 - computing refined Ritz - not iterative refined Ritz.']);
     warning(['WARNING: MAXITREF = 1 - not recommended with this routine.']);
end
if maxitref < 100
    warning(['WARNING: MAXITREF should be at least 100 to ensure Iterative Refined Ritz are computed.']);
end
    
% Check value of METHOD
if length(method)  ~= 3, error('ERROR: METHOD must be NOR, AUG, or THK'); end
if (~strcmp(method,'NOR') && ~strcmp(method,'AUG') && ~strcmp(method,'THK'))
   error('ERROR: METHOD must be NOR, AUG, or THK.'); 
end

% Check the values of REORTH 
if length(reorth) ~= 3, error('ERROR: REORTH must be ONE or TWO.'); end
if (~strcmp(reorth,'ONE') && ~strcmp(reorth,'TWO') )
   error('ERROR: REORTH must be ONE or TWO.'); 
end

% Check value of SIGMA. 
if length(sigma)  ~= 2, error('ERROR: SIGMA must be LS or SS'); end
if (~strcmp(sigma,'SS') && ~strcmp(sigma,'LS') )
   error('ERROR: SIGMA must be LS or SS.'); 
end

% Interchange m and n so that size(A'A) = min(m,n). Avoids
% finding zero values when searching for the smallest singular values.
interchange = 0; if n > m && strcmp(sigma,'SS'), t=m; m=n; n=t; interchange = 1; end

% Determine value of m_b to use
if isempty(m_b)
  if strcmp(sigma,'LS'), m_b = 2; else, m_b = 15; end
end

% Preallocate memory for W and F. These matrices are full and resizing will cause
% an increase in cpu time.
W = zeros(m,m_b); F = zeros(n,1);

% If starting p0 is not given then set starting vector p0 to be a 
% (n x 1) matrix of normally distributed random numbers.
% p0 is not explicitly created and use first column of matrix V.
if isempty(V)
  V = zeros(n,m_b); % Preallocate memory for for all V. 
  V(:,1) = randn(n,1); 
else
  V(:,2:m_b) = zeros(n,m_b-1); % Preallocate memory for other columns of V.
end

% Check for input errors in the data structure for K, M_B, and TOL
if k       <= 0,   error('ERROR: K must be a positive value.'),  end
if k  >  min(n,m), error('ERROR: K must be less than min(n,m)'), end   
if m_b     <= 1,   error('ERROR: M_B must be greater than 1.'),  end
if tol     <  0,   error('ERROR: TOL must be non-negative.'),    end
if m_b >= min(n,m)
   m_b = floor(min(n,m)-0.1);
   warning(['Changing M_B to ',num2str(m_b)]);
   if m_b - k - 1  < 0
      k = m_b - 1; 
      warning(['Changing K to ',num2str(k)]);
   end
end
if m_b - k - 1  < 0
   m_b = k+1; 
   warning(['WARNING: Changing M_B to ',num2str(m_b)]);
end
if m_b >= min(n,m)
   m_b = floor(min(n,m)-0.1);
   k = m_b - 1; 
   warning(['WARNING: Changing K to ',num2str(k)]);
   warning(['WARNING: Changing M_B to ',num2str(m_b)]);
end
if ~isnumeric(V), error('ERROR: Incorrect starting vector P0.'), end
if (size(V,1) ~= n), error('ERROR: Incorrect size of starting vector P0.'), end

% Check value of TOL.
% Set tolerance to machine precision if tol < eps.
if tol < eps, tol = eps; warning('WARNING: Changing TOL to EPS');  end

% Determine which vectors to use for  thick-restarting
if isempty(roh)
   if strcmp(sigma,'LS'), roh = 'RITZ'; else, roh = 'HARM'; end 
else
  if length(roh) ~= 4, error('ERROR: Unknown value for RoH. RoH must be RITZ or HARM'); end
  if ~ischar(roh), error('ERROR: Unknown value for RoH. RoH must be RITZ or HARM'); end
  if (~strcmp(roh,'RITZ') && ~strcmp(roh,'HARM') )
     error('ERROR: Unknown value for RoH. RoH must be RITZ or HARM'); 
  end
end

%--------------------------%
% END: PARSE INPUT VALUES. %
%--------------------------%

%-----------------------------------------------------------%
% BEGIN: DESCRIPTION AND INITIALIZATION OF LOCAL VARIABLES. %
%-----------------------------------------------------------%
% <- DO NOT MODIFY ->

% Set converegnce values for testing to 0 or empty. 
conv_rz_tol =0; conv_hm_tol =0; conv_it_tol=0; conv_rz_sqrt=0; conv_hm_sqrt=0;
conv_it_sqrt = 0; num_conv_val=0; V_B_cv=[]; U_B_cv=[]; S_B_cv=[]; out_val = -1;
res_norm_rz=0;

% Set all values empty for computing singular triplets for matrix B
U_B_rz=[]; S_B_rz=[]; V_B_rz=[];        % Holds the singular triplets of B (aka Ritz)
U_B_hm1=[]; S_B_hm1=[];V_B_hm1=[];      % Holds the singular triplets and intermediate
R_hm=[]; V_B_hm2=[];V_B_hm_last=[];     % steps for computing harmonic Ritz
s=[]; V_B_hm3 =[];V_B_hm=[];S_B_hm=[];  % singular values of B
U_rf=[]; S_rf=[]; V_rf=[]; s_min=[];    % Holds iterative refined Ritz singular values of B
W_rf=[]; rconv=0; U_rf1=[]; S_rf1=[];   % and values used for convergence tests.
V_rf1=[]; s_min1=[]; W_rf1=[]; rconv1=0;
U_B_it=[]; S_B_it=[]; V_B_it=[];

% Initialize array for max/min. singular value of B depending on sigma. 
if strcmp(sigma,'SS'), dmax = Inf(m_b,1); end
if strcmp(sigma,'LS'), dmax = -Inf(m_b,1); end

% Initialization and description of local variables.   
B = [];               % Bidiagonal matrix.
Bsz =[];              % Size of the bidiagonal matrix (will be <= m_B)
delta = 0.9;          % Used to determine when the Iter. refined Ritz singular vectors
                      % are "close" enough to the desired Ritz singular vectors
eps23 = eps^(2/3);    % Two thirds of eps, used for Smax. Avoids using zero.
I=[];                 % Used for indexing.
iter = 1;             % Main loop iteration count.  
k_org = k;            % Set k_org to the original value of k. k may be adjusted during iterations.
Fnorm =[];            % Two norm of the residual vector F.
Ritz = 1;             % Toggle for Lanczos bidiagonalization on which vectors are used for restarting
Smax = -Inf;          % Holds the maximum value of all computed singular values of B est. ||A||_2.
Smin = Inf;           % Holds the minimum value of all computed singular values of B est. cond(A).
SVTol = min(sqrteps,tol);  % Tolerance to determine when Lanczos encounters linearly dependent vectors. 
S_rf_0=zeros(m_b,1);    % Used for iterative refined comparsion from last iteration.

% Initialization of values used for STATS output.
mprod = 0;              % Number of matrix vector products with A and A^T.
IterRefRestart = 0;     % Number of restarts with Iter Refined Vectors.
timeMatProds=0;         % Total time computing matrix vector products with A and A^T.
timeIterRef = 0;        % Total time computing Iter. Ref. values.
timeReorth = 0;         % Total time doing full reorthogonalization

%---------------------------------------------------------%
% END: DESCRIPTION AND INITIALIZATION OF LOCAL VARIABLES. %
%---------------------------------------------------------%

%---------------------------%
% BEGIN: ITERATION PROCESS. %
%---------------------------%
while (iter <= maxit)
    
    % Compute the Lanczos bidiagonalization decomposition.
    [V,W,F,B,mprod,timeMatProds,timeReorth] = ablanzbd(varargin{1},V,...
     W,F,B,k,interchange,m_b,n,m,mprod,Ritz,SVTol*Smax,reorth,iter,...
     timeMatProds,timeReorth);
 
    % Reset k back to the original value k_org.
    k = k_org;
    
    % Determine the size of the bidiagonal matrix B.
    Bsz = size(B,1);
    
    % Compute the norm of the vector F.
    Fnorm = norm(F); 
   
    % Compute singular triplets of B. MATLAB's svds orders the singular values
    % largest to smallest.
    [U_B_rz,S_B_rz,V_B_rz] = svd(B); S_B_rz = diag(S_B_rz); 
    
    % Estimate ||A|| using the largest singular value over all iterations
    % and estimate the cond(A) using approximations to the largest and smallest 
    % singular values. If a small singular value is less than sqrteps use only Ritz
    % vectors for thick restarted and require two-sided reorthogonalization.
    if iter==1 
       Smax = S_B_rz(1); Smin = S_B_rz(Bsz);
    else 
       Smax = max(Smax,S_B_rz(1)); Smin = min(Smin,S_B_rz(Bsz));
    end
    Smax = max(eps23,Smax);
    if Smin/Smax < sqrteps, reorth = 'TWO'; roh='RITZ'; end
    
    % Re-order the singular values accordingly. MATLAB's SVD orders the 
    % singular values largest to smallest.
    if strcmp(sigma,'SS')
        [~,I] = sort(S_B_rz,'ascend');
        U_B_rz = U_B_rz(:,I); V_B_rz = V_B_rz(:,I); S_B_rz = S_B_rz(I);  
    end
    
    % If selected compute Harmonic Ritz values and vectors. 
    U_B_hm=[]; S_B_hm=[]; V_B_hm=[]; % Reset to empty.
    if strcmp(roh,'HARM') 
      R_hm = Fnorm*U_B_rz(Bsz,:); 
      % Update the SVD of B to be the SVD of [B ||F||E_m].
      [U_B_hm1,S_B_hm1,V_B_hm1] = svd([diag(S_B_rz) R_hm']); S_B_hm1 = diag(S_B_hm1);   
      if strcmp(sigma,'SS') % reorder if looking for smallest.
         [~,I] = sort(S_B_hm1,'ascend'); 
         U_B_hm1 = U_B_hm1(:,I); V_B_hm1 = V_B_hm1(:,I); S_B_hm1 = S_B_hm1(I);
      end
      U_B_hm1 = U_B_rz*U_B_hm1; U_B_hm = U_B_hm1; 
      V_B_hm2 = [[V_B_rz; zeros(1,Bsz)] flipud(eye(Bsz+1,1))]*V_B_hm1;
      V_B_hm_last = V_B_hm2(Bsz+1,:); % Set equal to the last row of V_B_hm2.
      s = Fnorm*(B\flipud(eye(Bsz,1))); 
      V_B_hm3 = V_B_hm2(1:Bsz,:) + s*V_B_hm2(Bsz+1,:);
      
      % Compute the orthogonal harmonic to get the Rayleigh Quotient values
      % to test for convergence. 
      [V_B_hm,R_Q] = qr(V_B_hm3,0); 
      D = diag(sign(diag(R_Q))); % Compute signs, so diag(R_Q) > 0 and signs of
      V_B_hm = V_B_hm*D;         % cols. V_B_hm  & V_B_hm3 are same. Note: D*D = I.
      S_B_hm = diag(U_B_hm'*B*V_B_hm);  % Rayleigh Quotient values
      
      % Need to reorder the Harmonic values based on Rayleigh Quotient.
      if strcmp(sigma,'SS') 
         [~,I] = sort(S_B_hm,'ascend');
      else
         [~,I] = sort(S_B_hm,'descend');
      end
      U_B_hm = U_B_hm(:,I); V_B_hm = V_B_hm(:,I); S_B_hm = S_B_hm(I);
    end

    % Used later for calling iterative refined function
    if strcmp(sigma,'LS')
     for j=1:k
        dmax(j) = max([dmax(j),S_B_rz(j)]);
     end
    else
      for j=1:k 
        dmax(j) = min([dmax(j),S_B_rz(j)]);
      end
    end
    
    % Convergence tests for Ritz and if computed Harmonic Ritz singular values/vectors.
    res_norm_rz_pre = res_norm_rz;
    [conv_rz_tol,conv_hm_tol,conv_it_tol,conv_rz_sqrt,conv_hm_sqrt,conv_it_sqrt,V_B_cv,U_B_cv,S_B_cv,...
    res_norm_rz,res_norm_hm,res_norm_it,out_val] = convtests(B,Bsz,Fnorm,tol,k_org,U_B_rz,...
    S_B_rz,V_B_rz,U_B_hm,S_B_hm,V_B_hm,[],[],[],Smax,sqrteps);
    num_conv_val = max([num_conv_val,conv_rz_sqrt,conv_hm_sqrt,conv_it_sqrt]);
    
    % If all desired singular values converged then exit main loop.
    if conv_rz_tol >= k_org || conv_hm_tol >= k_org || iter >= maxit, break,  end
  
    %--------------------------------------------------------------------%
    % BEGIN: COMPUTE ITERATIVE REFINED AND DETERMINE RESTARTING VECTORS. %
    %--------------------------------------------------------------------%
  
    % Compute Iterative refined Ritz values.
    rconv = 0; % set the number of convergence iterative refined Ritz to 0.
    if ~strcmp(method,'THK') && iter > 1 
       IRstart = tic;
       [U_rf,S_rf,V_rf,rconv] = refined_ritz(dmax(1:k),Fnorm,abs(B(Bsz,Bsz)),B,maxitref,Bsz,k,method);
       timeIterRef = timeIterRef + toc(IRstart); 
    end
    
    % rconv returns 0 or k. 0 indicates not all converged. rconv = k indicates all k converged 
    % Compute inner product between the k converged iterative refined Ritz vector(s) 
    % and Ritz vector(s). Also, check if iterative refined Ritz values are less than
    % the previous Ritz values. If so, do not use iterative refined Ritz
    % vectors to restart.
    ang_rz_rf = 0; decr = 0;
    if rconv > 0 
        ang_rz_rf = abs(diag(V_B_rz(:,1:k)'*V_rf(:,1:k))); % compute the innner prod. to check angle 
        for i=1:k_org
           if k_org > 1 || (k_org == 1 && Bsz > 2) % Requires iter. refined to be better approx. prev. vals.
             if strcmp(sigma,'LS') && (S_rf(i) + eps < S_rf_0(i)), decr = 1;  end
             if strcmp(sigma,'SS') && (S_rf(i) > S_rf_0(i) + eps),  decr = 1; end
           end
        end
    end
    
    % Reset S_rf_0 to be prev. dmax values for next iteration comparision
    S_rf_0 = dmax;
     
    % Find all inner products between iterative refined Ritz and Ritz >
    % delta = 0.9
    I = find(ang_rz_rf > delta);
 
    % Convergence check. 
    if length(I) == k && decr == 0 
       U_rf1 = U_rf; V_rf1 = V_rf; % if converged use iterative refined Ritz
       iterrefined = 1;
    else 
       % Use refined Ritz on augmented system - when all iterative refined Ritz
       % failed to converge.
       [U_rf1,S_rf1,V_rf1,rconv1] = refined_ritz(dmax(1:k),Fnorm,abs(B(Bsz,Bsz)),B,1,Bsz,k,'AUG');
       iterrefined = 0;
    end
    
    [V_B_it,V_R] = qr(V_rf1(:,1:k),0); % QR to ensure orthogonal vectors
    D = diag(sign(diag(V_R)));         % Compute signs and set diag(R) > 0
    V_B_it = V_B_it*D;                 % cols.V_B_it & V_rf1 are same. Note: D*D = I
    [U_B_it,U_R] = qr(U_rf1(:,1:k),0); % QR to ensure orthogonal vectors
    D = diag(sign(diag(U_R)));         % Compute signs and set diag(R) > 0
    U_B_it = U_B_it*D;                 % cols. U_B_it & U_rf1 are same. Note: D*D = I
    S_B_it = diag(U_B_it'*B*V_B_it);   % Rayleigh Quotient values - no need to reorder
 
    % Double check vectors are close to Ritz vectors.
    [~,R_ch] = qr(V_B_rz(:,1:k)'*V_B_it(:,1:k),0);
    if abs(prod(diag(R_ch))) > sqrteps
      % Convergence tests for Ritz singular value and to determine if to compute Iter. Refined.
      [conv_rz_tol,conv_hm_tol,conv_it_tol,conv_rz_sqrt,conv_hm_sqrt,conv_it_sqrt,V_B_cv,U_B_cv,S_B_cv,...
       res_norm_rz,res_norm_hm,res_norm_it,out_val] = convtests(B,Bsz,Fnorm,tol,k,U_B_rz,...
       S_B_rz,V_B_rz,[],[],[],U_B_it,S_B_it,V_B_it,Smax,sqrteps);
      num_conv_val = max([num_conv_val,conv_rz_sqrt,conv_hm_sqrt,conv_it_sqrt]); 
    end
   
    % If all desired singular values converged then exit main loop.
    if conv_it_tol >= k_org, break, end
    
    % Do not restart with iterative refined Ritz vectors if no  change in
    % values.
    if rconv > 0
      all_conv = 0;
      for i=1:k_org
        if abs(S_rf(i) - dmax(i)) < eps*10, all_conv = all_conv +1; end       
      end 
      if all_conv >= k_org, rconv = 0; end
    end
    
    % Do not restart consecutively with iterative refined Ritz vectors if last restart with 
    % iterative refined Ritz vectors caused norm of Ritz values/vectors to increase.
    if length(I) == k &&  rconv > 0 && decr == 0 && ~strcmp(method,'THK') && Ritz == 0
        if max(res_norm_rz) < max(res_norm_rz_pre), Ritz = 1; end
    end
    
    % Determine if restart with linear combination of iterative refined Ritz
    % or thick-restarted with Ritz or harmoinc Ritz vectors
    
    % Check to determine if restart with iterative refined Ritz or thick-restart
    if length(I) == k && rconv > 0 && decr == 0 && ~strcmp(method,'THK') && (Ritz ~= 0 || k == 1) 
        
      % Count the number of restarts with Iterative Refined. 
      IterRefRestart = IterRefRestart + 1;
       
      %------------------------------------------%
      % BEGIN: RESTART ITERATIVE REFINED SECTION %
      %------------------------------------------% 
      
      % Find the coefficients c_rf to setp linear combination.
      % See reference [1] for details. 
      %
      % Follows the equivalent eigenvalue problem.
      % A^T*AV = V*(B^T*B) + B(Bsz,Bsz)*f*e^T
      % Therefore, the approximate eigenpair is (S_rf^2,V_rf) 
      
      % Re-compute residual error of iterative refined vectors. 
      c_rf = ones(length(S_rf),1); 
      if k > 1
         
         % For restarting need the square of singular value - i.e. eigenvalue of T = B'*B.
         S_rf = S_rf.^2; T = B'*B; 
         
         % Compute norm of residual for A^T*AV = V*(B^T*B) + B(Bsz,Bsz)*f*e^T
         for i=1:k
            res_norm_it(i,1) = sqrt(norm((T - S_rf(i)*eye(Bsz))*V_rf(:,i))^2 + (abs(B(Bsz,Bsz))*Fnorm*abs(V_rf(Bsz,i)))^2);
         end

         % Initialize the k coefficients c_rf for the linear combination of
         % the k iterative refined Ritz singlar vectors V_rf.
         Bc = zeros(k-1,k); Tem = T(Bsz,1:Bsz); Bc(1,1:k) = V_rf(Bsz,1:k);
         for i=2:k-1
            for j=1:k
                Bc(i,j) = Tem*V_rf(:,j)*S_rf(j)^(i-2);
             end 
         end
        
         % Compute the solution of the k-1 x k homogeneous system using
         % null. Following code is needed to avoid zero column(s) in Bc
         % for numerically converged eigenvectors of B'*B.
         Iin = []; Iout=[]; % set variables for which columns are used.
         Bc_max = max(abs(Bc),[],'all');
         for i=1:k  % search for columns numerically zero.
            if norm(Bc(:,i),Inf) < sqrteps*Bc_max
              Iout = [Iout i]; % references which columns to remove.
            else
              Iin = [Iin i];
            end
         end  
         if ~isempty(Iout) % remove numerically zero columns from Bc
            Bc = Bc(1:k-length(Iout)-1,Iin);
            if ~isempty(Bc) 
               c_rf(1:k) = zeros(k,1);
               c_rf_in = null(Bc); % call Matlab's null to solve system
               if ~isempty(c_rf_in), c_rf(Iin) = c_rf_in(:,1); end
               c_rf(Iout) = res_norm_it(Iout); % Place zero entries back in sol.
            else
               c_rf = res_norm_it(1:k);
            end   
         else
            c_rf = null(Bc); % call Matlab's null to solve full system
         end
         c_rf = c_rf(:,1); % matlab's null can return more than one sol.  
      end  
      % The coefficents c_rf have been computed at this point. 

      % Compute the starting vector 
      % v_bar = c_rf(1)*V_rf(:,1)+ .... + c_rf(k)*V_rf(:,k) 
      v_bar = V_rf*c_rf; v_bar = v_bar/norm(v_bar);
      
      % Compute u_bar = B*v_bar and B(1,1).
      % Now have A*V*v_bar = W*u_bar*alpha (avoids computing A*V*v_bar)
      u_bar = B*v_bar; alpha = norm(u_bar); u_bar = u_bar/alpha;
      
      % Compute F = A^T*W*u_bar - alpha without using A^T*W*u_bar
      F = V*(B'*u_bar - alpha*v_bar) + F*u_bar(Bsz);
      
      % Reset B matrix.
      B=[]; B(1,1) = alpha;
      
      % Updated first column of V. Do not resize - slows down code.
      V(:,1) = V*v_bar;
      
      % Updated first column of W. Do not resize - slows down code.
      W(:,1) = W*u_bar;
      
      % Reorthogonalize the residual vector (indepdendent of REORTH)
      ORstart = tic;
      F = F - (V(:,1)'*F)*V(:,1);
      timeReorth = timeReorth + toc(ORstart);
      
      % Compute the normn of F, scale and update V and B matrices for restarting. 
      beta = norm(F); V(:,2) = F/beta; B(1,2) = beta;
      
      % Set the Ritz for which vectors restarting
      % Lanczos bidiagonalization decomposition
      % Ritz = 0 => iterative refined Ritz vectors
      % Ritz = 1 => thick-restarting with harmonic or Ritz
      Ritz = 0; 
      %----------------------------------------%
      % END: RESTART ITERATIVE REFINED SECTION %
      %----------------------------------------% 
    else
      %-----------------------------------%
      % BEGIN: THICK-RESTART RITZ SECTION %
      %-----------------------------------%
      
      % Simple strategy to improve convergence. Adjust k value.
      k = k_org + num_conv_val;
      if k < Bsz-1 && k > 1 % used for small values of Bsz
          if abs(S_B_rz(k+1) - S_B_rz(k)) < abs(S_B_rz(k-1) - S_B_rz(k))
              k = k+1;
          end
      end
      k = max(floor((Bsz+num_conv_val)/2),k); % used for large Bsz 
      if k >= Bsz, k = Bsz - 1; end
         
      % Use Harmonic Ritz vectors for thick-restarting
      if strcmp(roh,'HARM') 
          % Vectors are not orthogonal.
          [V_B_hm,R] = qr([ [V_B_hm3(:,1:k); zeros(1,k)]  [-s; 1] ],0);
          V(:,1:k+1) = [V F/Fnorm]*V_B_hm;
  
          % Update and compute the K x K+1 part of B.
          B = diag(S_B_hm1(1:k))*triu((R(1:k+1,1:k)+ R(:,k+1)*V_B_hm_last(1:k))');  
          
         % Update W 
         W(:,1:k) = W*U_B_hm1(:,1:k);
         
      else
         
         % Use Ritz vectors for thick-restarting  
         R = Fnorm*U_B_rz(Bsz,:); F = F/Fnorm;  
         V(:,1:k+1) = [V*V_B_rz(:,1:k) F];    
         B = [diag(S_B_rz(1:k)), R(1:k)'];
         
         % Update W 
         W(:,1:k) = W*U_B_rz(:,1:k);
      
      end
      
      % Set the Ritz for which vectors restarting
      % Lanczos bidiagonalization decomposition
      % Ritz = 0 => iterative refined Ritz vectors
      % Ritz = 1 => thick-restarting with harmonic or Ritz
      Ritz = 1; 
      
      %---------------------------------%
      % END: THICK-RESTART RITZ SECTION %
      %---------------------------------%
    end 
 
 % Update the main iteration loop count. 
 iter = iter+1;

end % end main loop
%-------------------------%
% END: ITERATION PROCESS. %
%-------------------------%

%-----------------------%
% BEGIN: OUTPUT RESULTS %
%-----------------------%

% Output option I: Display singular values only.
if (nargout == 0)
    if iter >= maxit
       disp('Maximum number of iterations exceeded.'); 
       disp('Use option IV for best estimate.');
       SingularValues=[]
    else
      SingularValues = S_B_cv(1:k_org)
    end
end

% Output option II: Set singular values equal to output vector.  
if (nargout == 1)
    if iter >= maxit
      disp('Maximum number of iterations exceeded.'); 
      disp('Use option IV for best estimate.');
      varargout{1}=[];
    else
      varargout{1} = S_B_cv; 
    end
end

% Output option III and IV: Output singular triplets (U,S,V)
if nargout > 1
    if iter >= maxit && nargout == 3
       disp('Maximum number of iterations exceeded.'); 
       disp('Use option IV for best estimate.');
       varargout{1}=[]; varargout{2}=[]; varargout{3}=[];
    else   
       if interchange
         varargout{1} = V*V_B_cv(:,1:k_org); 
         varargout{2} = diag(S_B_cv(1:k_org)); 
         varargout{3} = W*U_B_cv(:,1:k_org);
       else
         varargout{1} = W*U_B_cv(:,1:k_org); 
         varargout{2} = diag(S_B_cv(1:k_org)); 
         varargout{3} = V*V_B_cv(:,1:k_org);
       end
    end
    % Output options IV: Output singular triplets (U,S,V) and STATS
   if nargout == 4
      STATS.numMatProds = mprod;
      STATS.timeMatProds = timeMatProds;
      STATS.numIterRefRestart = IterRefRestart;
      STATS.timeCompIterRef = timeIterRef;
      STATS.timeReorth = timeReorth;
      STATS.timeTotal = toc(tStart);
      STATS.estimateSVmax = Smax;
      STATS.estimateSVmin = Smin;
      if iter >= maxit
         STATS.convergedKVals='FALSE';
      else
         STATS.convergedKVals='TRUE';
      end
      if out_val == 1
        STATS.outputestRitzorIterRef='RITZ SING.';   
        STATS.maxnormres = max(res_norm_rz);
      elseif out_val == 2
          STATS.outputestRitzorIterRef='HARM SING.'; 
        STATS.maxnormres = max(res_norm_hm);
      elseif out_val == 3
        if iterrefined  
           STATS.outputestRitzorIterRef='ITER. REF. SING';
        else
           STATS.outputestRitzorIterRef='REF. SING';
        end
        STATS.maxnormres = max(res_norm_it);
      elseif out_val == -1
          STATS.outputestRitzorIterRef='ERROR';
          STATS.maxnormres = 'ERROR';
      end   
      varargout{4}=STATS;
   end
end

%---------------------%
% END: OUTPUT RESULTS %
%---------------------%
end

%------------------------------------------------%
% BEGIN: LANCZOS BIDIAGONALIZATION DECOMPOSITION %
%------------------------------------------------%

function [V,W,F,B,mprod,timeMatProds,timeReorth] = ablanzbd(A,V,W,F,B,K,interchange,m_b,n,m,mprod,Ritz,SVTol,reorth,iter,timeMatProds,timeReorth)
% Computes the Lanczos bidiagonalization decomposition
%  
%  A*V  = W*B
%  A'*W = V*B' + F*E^T
%  
% with full reorthogonalization. 
% If the m x n matrix A is a function, 'Afunc', 
% then the structure must be y = Afunc(x,m,n,'transpose'). If transpose = 'F', 
% then y = A*x. If transpose = 'T', then y = A'*x.

% James Baglama
% DATE: 6/22/21

% Initialization of main loop count J.
J = 1; 

% Normalize starting vector.
if iter == 1 
   V(:,1) = V(:,1)/norm(V(:,1)); B=[];
else
   if Ritz, J = K+1; else, J = 2; end
end

% Matrix A product with vector(s) V, (W=A*V).
Astart = tic;
if interchange
   if ischar(A)
     W(:,J) = feval(A,V(:,J),m,n,'T');
  else
     W(:,J) = A'*V(:,J);
  end
else
  if ischar(A)
     W(:,J) = feval(A,V(:,J),m,n,'F');
  else
     W(:,J) = A*V(:,J);
  end
end
timeMatProds = timeMatProds + toc(Astart); 

% Count the number of matrix vector products.
mprod = mprod + 1;

% Input vectors are singular vectors and AV(:,J) which must be orthogonalized.
if iter ~= 1
    W(:,J) = orthog(W(:,J),W(:,1:J-1));
   
   % Second orthogonalization step required. Input vectors not always 
   % strongly orthogonal.
   W(:,J) = orthog(W(:,J),W(:,1:J-1)); 

end

% Compute the norm of W.
S = norm(W(:,J));

% Check for linearly dependent vectors.
if S <= SVTol
   W(:,J) = randn(size(W,1),1);
   W(:,J) = orthog(W(:,J),W(:,1:J-1));
   W(:,J) = W(:,J)/norm(W(:,J));
   S = 0;
else
   W(:,J) = W(:,J)/S;
end

% Begin of main iteration loop for the block Lanczos bidiagonalization decomposition.
while (J <= m_b)

   % Matrix A' product with vector(s), (F = A'*W).
   Astart = tic;
   if interchange
      if ischar(A)
         F = feval(A,W(:,J),m,n,'F');
      else
         F = A*W(:,J);
      end
   else
      if ischar(A)
        F = feval(A,W(:,J),m,n,'T');
      else
        F = A'*W(:,J);
      end
   end
   timeMatProds = timeMatProds + toc(Astart); 
   
   % Count the number of matrix vector products.
   mprod = mprod + 1;
   
   % One step of the Gram-Schmidt process.
   F = F - V(:,J)*S;
   
   % Second step to maintain strong local orthogonality 
   S2 = F'*V(:,J); F = F - V(:,J)*S2; S = S + S2;
   
   % Always perform full reorthogonalization step of "Short vectors".
   ORstart = tic;
   if J > 1, F = orthog(F,V(:,1:J-1)); end
   timeReorth = timeReorth + toc(ORstart); 
   
   if J+1 <= m_b

      % Compute the norm of F.
      R = norm(F); 
      
      % Check for linearly dependent vectors.
      if R <= SVTol
         F = randn(size(V,1),1); 
         F = orthog(F,V(:,1:J));
         V(:,J+1) = F/norm(F);
         R = 0;
      else
        V(:,J+1) = F/R;
      end 

      % Compute bidiagonal matrix B. 
      if isempty(B) 
         B = [S R];
      else
         B = [B zeros(J-1,1); zeros(1,J-1) S R];
      end 
    
      % Matrix A product with vector(s), (W=A*V).
      Astart = tic;
      if interchange
         if ischar(A)
           W(:,J+1) = feval(A,V(:,J+1),m,n,'T');
         else
           W(:,J+1) = A'*V(:,J+1);
         end
      else
         if ischar(A)
            W(:,J+1) = feval(A,V(:,J+1),m,n,'F');
         else
            W(:,J+1) = A*V(:,J+1);
         end
      end
      timeMatProds = timeMatProds + toc(Astart); 

      % Count the number of matrix vector products.
      mprod = mprod + 1;

      % One step of the Gram-Schmidt process.
      W(:,J+1) =  W(:,J+1) - W(:,J)*R;
      
      % Second step for local orthogonality 
      R2 = W(:,J+1)'*W(:,J); W(:,J+1) = W(:,J+1) - W(:,J)*R2; R = R + R2;
      
      % Full Reorthogonalization step. "Long vectors"
      ORstart = tic;
      if ( iter == 1 || strcmp(reorth,'TWO') )
         if J > 1
           W(:,J+1) = orthog(W(:,J+1),W(:,1:J-1)); 
         end
      end
      timeReorth = timeReorth + toc(ORstart);
   
      % Compute the norm of W.
      S = norm(W(:,J+1)); 
      
     % Check for linearly dependent vectors.
     if S <= SVTol
        W(:,J+1) = randn(size(W,1),1);
        W(:,J+1) = orthog(W(:,J+1),W(:,1:J));
        W(:,J+1) = W(:,J+1)/norm(W(:,J+1));
        S = 0;
     else
        W(:,J+1) = W(:,J+1)/S;
     end 

   else
    
    % Add last element to matrix B
     B = [B; zeros(1,J-1) S]; 

   end

    % Update iteration count.
    J = J + 1;

end
end

function y = orthog(y,X)
% Simple re-orthogonalization of vector y against columns of matrix X.
for i=1:size(X,2)
  dotY = y'*X(:,i);
  y = y - X(:,i)*dotY;
end

end
 
%----------------------------------------------%
% END: LANCZOS BIDIAGONALIZATION DECOMPOSITION %
%----------------------------------------------%

%--------------------------%
% BEGIN: CONVERGENCE TESTS %
%--------------------------%

function [conv_rz_tol,conv_hm_tol,conv_it_tol,conv_rz_sqrt,conv_hm_sqrt,...
          conv_it_sqrt,V_B_cv,U_B_cv,S_B_cv,res_norm_rz,res_norm_hm,res_norm_it,out_val] = ...
          convtests(B,Bsz,Fnorm,tol,k,U_B_rz,S_B_rz,V_B_rz,U_B_hm,...
          S_B_hm,V_B_hm,U_B_it,S_B_it,V_B_it,Smax,sqrteps)
% This function checks the convergence of singular triplets 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convergence check for all input values.                          %
%                                                                  %
% SQRT ( || A v - s u ||^2 + || A^T u - s v ||^2 )<= TOL*||A||_2   %
%                                                                  %
% Use the Lanczos bidiagonalization decomposition relationship     %
%  A V   = W B                                                     %
%  A^T W = V B^T + F E^T                                           %
%                                                                  %
% || A v - s u ||^2 = ||A V v - s W u||^2 = || W B v - s W u ||^2  %
%                                         = || B v - s u ||^2      %
%                                                                  %
% || A^T u - s v ||^2 = ||A^T W u - s V v||^2                      %
%                     = ||V B^T u + F E^T u - s V v||^2            %
% (since V^TF = 0)    = ||B^T u -s v||^2 + ||F E^T u||^2           %
%                     = ||B^T u -s v||^2 + (| E^T u| ||F||)^2      %
%                                                                  %
% RITZ => B v = s u and B^T u = s v and                            %
% || A v - s u ||   = 0                                            %
% || A^T u - s v || = | E^T u| ||F||                               %
% Therefore,                                                       %
%                                                                  %
% SQRT( || A v - s u ||^2 + || A^T u - s v ||^2 ) = | E^T u| ||F|| % 
%                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Baglama
% DATE: 3/30/21

% Initialize output values. 
conv_rz_tol = 0;  conv_hm_tol=0;  conv_it_tol=0; 
conv_rz_sqrt = 0; conv_hm_sqrt=0; conv_it_sqrt=0;
V_B_cv=[]; U_B_cv=[]; S_B_cv=[]; 
res_norm_hm = ones(k,1); res_norm_it = ones(k,1);

% Compute the residual for Ritz values
res_norm_rz = (abs(U_B_rz(Bsz,1:k))*Fnorm)';
for i=1:k
   if res_norm_rz(i,1) < tol*Smax,  conv_rz_tol  = conv_rz_tol+1;  end
   if res_norm_rz(i,1) < sqrteps*Smax, conv_rz_sqrt = conv_rz_sqrt+1; end
end

% Compute the residual for Harmonic Ritz values.
if ~isempty(S_B_hm)
  for i=1:k
    res_norm_1(i,1) = norm(B*V_B_hm(:,i) - S_B_hm(i)*U_B_hm(:,i))^2;
    res_norm_2(i,1) = norm(B'*U_B_hm(:,i) - S_B_hm(i)*V_B_hm(:,i))^2 + (U_B_hm(Bsz,i)*Fnorm)^2;
    res_norm_hm(i,1) = sqrt(res_norm_1(i,1)+res_norm_2(i,1));
    if res_norm_hm(i,1) < tol*Smax,  conv_hm_tol  = conv_hm_tol+1;  end
    if res_norm_hm(i,1) < sqrteps*Smax, conv_hm_sqrt = conv_hm_sqrt+1; end
  end
end

% Compute the residual for Itertive Refined Ritz values.
if ~isempty(S_B_it)
  for i=1:k
    res_norm_1(i,1) = norm(B*V_B_it(:,i) - S_B_it(i)*U_B_it(:,i))^2;
    res_norm_2(i,1) = norm(B'*U_B_it(:,i) - S_B_it(i)*V_B_it(:,i))^2 + (U_B_it(Bsz,i)*Fnorm)^2;
    res_norm_it(i,1) = sqrt(res_norm_1(i,1)+res_norm_2(i,1));
    if res_norm_it(i,1) < tol*Smax,  conv_it_tol  = conv_it_tol+1;  end
    if res_norm_it(i,1) < sqrteps*Smax, conv_it_sqrt = conv_it_sqrt+1; end
  end  
end

% Output values.
[~,I] = max([conv_rz_tol conv_hm_tol conv_it_tol]);
if I == 1, out_val = 1; U_B_cv = U_B_rz; V_B_cv= V_B_rz; S_B_cv = S_B_rz; end  % Output Ritz (Default)
if I == 2, out_val = 2; U_B_cv = U_B_hm; V_B_cv= V_B_hm; S_B_cv = S_B_hm; end  % Output Harmonic
if I == 3, out_val = 3; U_B_cv = U_B_it; V_B_cv= V_B_it; S_B_cv = S_B_it; end  % Output Iter. Ref. or Ref.

end
%------------------------%
% END: CONVERGENCE TESTS %
%------------------------%

%---------------------------------------%
% BEGIN: COMPUTE REFINED RITZ ITERATION %
%---------------------------------------%

function [u,rho,v_min,rconv] = refined_ritz(D_ritz,R,alpha,B,maxitref,Bsz,k,method)
% Computes the Iterative refined Ritz values and vectors.   
%
% INPUT:
%
%   D_RITZ  - (K x 1) vector of approx. singular values.
%        R  - real number - norm of residual vector F.
%   ALPHA   - real number - last diagonal element of B -> B(Bsz,Bsz).
%      B    - Bsz x Bsz bidiagonal matrix
% MAXITREF  - Integer indicating the maximum number of iterations for the iterative Refined Ritz values.
%             set to 1 and the refined Ritz are computed - used in  convergence testing.
%     BSZ  -  Integer indicates the size of the  tridiagonal matrix.
%     K    - number of Iterative refined Ritz values/vectros to compute
%   METHOD -   NOR -  Equivalent eigenvalue system
%                     AP_m = Q_mB_m  and A'Q_m = P_m B_m' + f_me_m'
%                     Equivalent eigenvalue system 
%                     A'AP_m = P_m B_m'B_m + B_m(m,m)f_me_m'
%                     Compute iterative refined Ritz on [B'*B; 0 B_m(m,m)*norm(f)]
%              AUG -  Equivalent eigenvalue system                    
%                     [0 A; A' 0][Q_m 0; 0 P_m] = [Q_m 0; 0 P_m][0 B; B' 0]+[ 0 0; f_me_m' 0] 
%                     Compute iterative refined Ritz on [0 B; B' 0; norm(f) 0]
%
% OUTPUT:
%      U   - (BSZ x K) matrix of left singular vectors:
%            METHOD 
%               NOR =>  U = B*V_min; U = U/norm(U); 
%               AUG =>  U = V_min(1:Bsz/2,Bsz); U = U/norm(U);           
%     RHO  - (K x 1) vector of iterative refined Ritz values.
%   V_MIN  - (BSZ x K) matrix of right singular values
%            METHOD 
%               AUG =>  V_MIN = V_MIN(Bsz/2+1:Bsz,Bsz); V_MIN = V_MIN/norm(V_MIN); 
%   RCONV  - Integer indicate if all K iterative refined Ritz converge.
%          - 0 not all converged - do not use
%          - K all converged
%
%
%  DATE MODIFIED: 5/10/2021
%  VER:  1.0

% Set up the augmented matrix [B'*B; 0 alpha*R]
if strcmp(method,'NOR')
   B_aug = zeros(Bsz+1,Bsz); B_aug(1:Bsz,1:Bsz) = B'*B;
   B_aug(Bsz+1,Bsz) = alpha*R;
   
   % Initialize rho value.
   rho = D_ritz.^2; 
   
else
    
   % Set up the augmented matrix [0 B; B' 0; R 0] 
   B_aug = [zeros(Bsz) B; B' zeros(Bsz); zeros(1,2*Bsz)];
   B_aug(2*Bsz+1,Bsz) = R; 
   
   % Initialize rho value.
   rho = D_ritz;  
   
end


% Initialize values.
s_min = zeros(k,1);           % Intialize min. singular values of B_aug. 
v_min = zeros(Bsz,k);         % Intialize right singular vectors of B_aug. 
u     = zeros(Bsz,k);         % Intialize left singular vectors of B_aug.
rconv = zeros(k,1);           % Intialize convergence.
rho_0 = rho;                  % Used to check for convergence during iterations.
sqrteps = sqrt(eps);          % square root of machine precision - eps
v_zero = zeros(Bsz,1);        % Set starting vector to zero. 

% Change size for augmented system.
if strcmp(method,'AUG'), Bsz = 2*Bsz; end

% Compute k number of iterative refined Ritz values/vectors.
% Compute in reverse order. Less likely k, k-1, .. converge - avoids
% uneeded iteration. Exit if any iterative refined Ritz fail to converge. 
for j = k:-1:1

    %  Set difference in rho values to test for stagnation.
    diff_rho_0 = -1; v_min_0 = v_zero;
   
    % Iteratation to compute the iterative refined Ritz values/vectors
    for i=1:maxitref 
        
      % Compute the SVD of  B_aug - rho* I
      [U,S,V] = svd((B_aug-rho(j)*eye(Bsz+1,Bsz)),0);
      
      % Need the smallest singular triplet of B_aug - rho* I. Matlab
      % returns order of singular values largest to smallest.
      if strcmp(method,'NOR')
         s_min(j) = S(Bsz,Bsz); v_min(:,j) = V(:,Bsz); 
      
         % Compute the new rho = v_min'*B'*B*v_min = ||B*v_min||^2 
         rho(j) = norm(B*v_min(:,j))^2; diff1 = 0;
        
      else
        
        s_min(j) = S(Bsz,Bsz); 
        v_min(:,j) = V(Bsz/2+1:Bsz,Bsz); v_min(:,j) = v_min(:,j)/norm(v_min(:,j));    
        
        % Compute the new rho = V(:,Bsz)'*[0 B; B' 0]*V(:,Bsz) 
        rho(j) = V(:,Bsz)'*B_aug(1:Bsz,:)*V(:,Bsz);  
        % check norm singular vectors are close to 1/sqrt(2)
        diff1 = abs(norm(V(Bsz/2+1:Bsz,Bsz)) - 1/sqrt(2));  
        
      end    
          
      % Compute the difference of previous rho to check for convergence 
      den = max([rho_0(j),rho(j), eps]); 
      diff_rho = abs( (rho_0(j) - rho(j))/den );
      
      % Check for convergence 
      if ((((diff_rho < eps && abs(V(:,Bsz)'*U(1:Bsz,Bsz)) < sqrteps)...
         || s_min(j) < eps || norm(abs(v_min_0) - abs(v_min(:,j))) < eps) && diff1 < sqrteps) || maxitref == 1) 
            rconv(j) = 1;
            if strcmp(method,'NOR')
               u(:,j) = B*v_min(:,j); 
               u(:,j) = u(:,j)/norm(u(:,j));
               rho(j) = sqrt(rho(j));
               
            else
                
               u(:,j) = V(1:Bsz/2,Bsz); 
               u(:,j) = u(:,j)/norm(u(:,j)); 
              
            end
            break; 
      end
    
      % Check for stagnation to avoid too many unecessary iterations. Care 
      % must be taken to avoid the situation where rho(k) is still changing 
      % very slightly at first and then an increase in change later. An
      % early termination due to stagnation with no converegnce 
      % may avoid increase in change later that converges. 
      if abs(diff_rho_0 - diff_rho)/max(diff_rho,eps) < eps,  break; end  % stagnation
      if i>= 10 && mod(i,10)==0 % update every 10 iterations to avoid early termination. 
         diff_rho_0 = diff_rho; 
      end 
      
     % Update rho_0 and v_min_0 to test convergence. 
     rho_0(j) = rho(j); v_min_0 = v_min(:,j); 
    end 
    
    if rconv(j) == 0, rconv = 0; return; end
    
end

% After K iteration finish check to see if *all* K iterative refined Ritz values
% have all converged. Reset rconv to return an integer value, rconv = k all
% converged and rconv = 0, not all converged.

if all(rconv == 1), rconv = k; else rconv = 0; end

end

%-------------------------------------%
% END: COMPUTE REFINED RITZ ITERATION %
%-------------------------------------%