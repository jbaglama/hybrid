function [U,S,V,STATS] = rd2svds(A,m,n,P,k,tol,tol_def)
% 
%  RD2SVDS: Computes the k largest singular value and associated singular vectors
%           of a m x n matrix A such that A*V = U*S and A'*U = V*S, V'*V=I
%           and U'*U = I and S is a diagonal matrix.
%
%  CONVERGENCE: sqrt(|| A*v - u*s||^2 + || A^T*u - v*s ||^2))<= tol*||A||
%               where norm is the 2-norm and ||A|| is approximated by largest 
%               singular value of the projected matrix B over all iterations. 
%
%  INPUT:
%   A        -  m x n numeric matrix A or an M-file ('Afunc'). If the m x n matrix A 
%               is a filename then y = Afunc(x,m,n,'transpose'). If transpose = 'F', 
%               then y = A*x. If transpose = 'T', then y = A'*x.
%   m,n      - size of the m x n matrix A
%   P        -  n x 1 "right" starting vector & first column of Lanczos matrix P - no vector P = []
%   k        - number of desird singular triplets
%  tol       - user specified tolerance - no input - tol = 1d-6
%  tol_def   - user specified deflation tolerance when k > 1 - must be < tol - no input tol_def=1d-1*tol
%
% OUTPUT:
%  U   - m x k  left approximate singular vectors
%  V   - n x k  right approximate singular vectors
%  S   - k x k  diagonal matrix of singular values 
%  If k failed to converged, output k-1 converged of [] 
%  STATS =
%   LockingProblem:    Determine if locking issue arises - cannot continue
%   numMatProds:       Number of matrix-vector products with A and A^T
%   timeMatProds:      Total time computing products with A and A^T 
%   numIterRefRestart: Number of restarts with Iterative Refined vectors
%   timeCompIterRef:   Time computing iterative Refined vectors
%   timeTotal:         Total time elasped
%   estimateSVmax:     Estimate of the maximum singular value over all iterations
%   convergedKVals:    True if ALL k singular triplets converged, otherwise false
%   maxnormres         Maximum residual of all converged values for exit
%
%  DATE MODIFIED: 12/27/21
%  VER:  2.0
%  
% AUTHORS: 
% James Baglama     email: jbaglama@uri.edu
% Vasilije Perovic  email: perovic@uri.edu
% Jennifer Picucci  email: jenniferpicucci@uri.edu
%
% REFERENCE:
% Baglama, J, Perovic, V, and Picucci, J, "Hybrid Iterative Refined Restarted 
% Lanczos Bidiagonalization Methods",  2021 submitted Numerical
% Algorithms, preprint: http://www.math.uri.edu/~jbaglama/paper34.pdf

% Start timing count using the MATLAB tic command.
tStart = tic; 

% Initialize values.
maxit_outer = 2000; maxit_inner = 100; sqrteps = eps^(1/2); f = zeros(n,1);  
V=zeros(n,k); U=zeros(m,k); S=zeros(k); B = zeros(2,2); Q = zeros(m,2); 
if isempty(tol),tol = 1d-6; end; tol = max(tol,eps); tol_k = tol; 
if isempty(tol_def), tol = 1d-1*tol; else, tol = tol_def; end  
if tol >= tol_k, error('ERROR: tol_def >= tol'); end; Smax_all = 0;
if isempty(P), P = randn(n,1); end; UTAP1 = 0; UTAP2 = 0; VTATQ1 = 0; VTATQ2 = 0;

% Initialization of  values used for STATS output.
mprod = 0;                % Number of matrix vector products with A and A^T.
numIterRef = 0;           % Number of restarts with Iter Refined Vectors.
timeMatProds = 0;         % Total time computing matrix vector products with A and A^T.
timeIterRef = 0;          % Total time computing Iter. Ref. values.
maxnormres = -1;          % Maximum norm of all converged singular triplets
convergedKVals = 'TRUE';  % Initialize to TRUE
LockingProblem = 'FALSE'; % Initialize to FALSE - no locking problems

% Begin iter_k for loop
for iter_k = 1:k
   
   % Set values for each k value
   iter = 1;  Smax = 0; if iter_k == k, tol = max(tol_k,eps); end  

   % Set up starting vector P(:,1)
   if iter_k > 1, P = f - V(:,1:iter_k-1)*(V(:,1:iter_k-1)'*f); end 
   P = [P/norm(P) zeros(n,1)]; 

   % Matrix-vector product, orthogonalization, and deflation to get Q(:,1)
   Astart = tic; if ischar(A), Q(:,1) = feval(A,P(:,1),m,n,'F'); else, ... 
   Q(:,1) = A*P(:,1); end; mprod = mprod + 1; timeMatProds = timeMatProds + toc(Astart);
   if iter_k > 1 % Deflation step
      UTAP1 = U(:,1:iter_k-1)'*Q(:,1); Q(:,1) = Q(:,1) - U(:,1:iter_k-1)*UTAP1; 
   end 
   B(1,1) = norm(Q(:,1)); Q(:,1) = Q(:,1)*(1/B(1,1)); 

   % Matrix-vector product, orthogonalization, and deflation to get P(:,2)
   Astart = tic; if ischar(A), P(:,2) = feval(A,Q(:,1),m,n,'T'); else,... 
   P(:,2) = A'*Q(:,1); end; mprod = mprod + 1; timeMatProds = timeMatProds + toc(Astart);
   if iter_k > 1 % Deflation step
      VTATQ1 = V(:,1:iter_k-1)'*P(:,2);P(:,2) = P(:,2) - V(:,1:iter_k-1)*VTATQ1; 
   end 
   P(:,2) = P(:,2) - P(:,1)*B(1,1); B(1,2) = norm(P(:,2)); P(:,2) = P(:,2)*(1/B(1,2));

   % Matrix-vector product, orthogonalization, and deflation to get Q(:,2)
   Astart = tic; if ischar(A), Q(:,2) = feval(A,P(:,2),m,n,'F'); else, ...
   Q(:,2) = A*P(:,2); end; mprod = mprod + 1; timeMatProds = timeMatProds + toc(Astart);
   if iter_k > 1 % Deflation step
      UTAP2 = U(:,1:iter_k-1)'*Q(:,2); Q(:,2) = Q(:,2) - U(:,1:iter_k-1)*UTAP2; 
   end
   Q(:,2) =  Q(:,2) - Q(:,1)*B(1,2); B(2,2) = norm(Q(:,2)); Q(:,2) = Q(:,2)*(1/B(2,2));

   % Matrix-vector product, orthogonalization, and deflation to get f
   Astart = tic; if ischar(A), f = feval(A,Q(:,2),m,n,'T'); else, ...
   f = A'*Q(:,2); end; mprod = mprod + 1; timeMatProds = timeMatProds + toc(Astart);
   if iter_k > 1 % Deflation step.
      VTATQ2 = V(:,1:iter_k-1)'*f; f = f - V(:,1:iter_k-1)*VTATQ2; 
   end 
   f = f - P(:,2)*B(2,2); f = f - P*(P'*f); fnorm = norm(f);

   % Begin main iter while loop for each singular triplet
   while (iter <= maxit_outer) 
    
       % Compute the largest singular value and associated vector of B 
       [u_b,s_b,v_b] = svd(B); v_b = v_b(:,1); u_b = u_b(:,1); s_b = s_b(1,1); Smax = max(Smax,s_b);
       
       % Compute iterative refined Ritz value/vectors
       IRstart = tic;
       iter_refined = 0; rho = Smax^2; rho_0 = rho; rconv = 0; v_rf_0 = zeros(2,1);
       T = B'*B; B_aug = [T; 0 B(2,2)*fnorm]; diff_rho_0 = -1;
       for i=1:maxit_inner
          B_aug(1,1) = T(1,1) - rho; B_aug(2,2) = T(2,2) - rho;
          [u_rf,s_rf,v_rf] = svd(B_aug,0);
          v_rf = v_rf(:,2); s_rf = s_rf(2,2); rho = norm(B*v_rf)^2;  
          diff_rho = abs((rho_0 - rho)/max(rho,rho_0)); 
          if (diff_rho < eps && abs(v_rf'*u_rf(1:2,2)) < sqrteps) ...
          || norm(abs(v_rf_0) - abs(v_rf)) < eps || s_rf < eps
             u_rf = B*v_rf; alpha = norm(u_rf); u_rf = u_rf/alpha; 
             rho = sqrt(rho); rconv = 1; break; % Converge iter. ref. exit
          end 
          if abs(diff_rho_0 - diff_rho)/max(diff_rho,eps) < eps, break; end  % Stagnate exit
          if i >= 10 && mod(i,10) == 0, diff_rho_0 = diff_rho; end 
          rho_0 = rho; v_rf_0 = v_rf; 
       end
       if abs(v_b'*v_rf) > 0.9 && rconv == 1 && iter>1, iter_refined = 1; end
       timeIterRef = timeIterRef + toc(IRstart);
      
       % Convergence check
       if ~iter_refined 
         B_aug = [-Smax*eye(2,2) B; B' -Smax*eye(2,2); zeros(1,4)]; B_aug(5,2) = fnorm; 
         [~,~,v] = svd(B_aug,0); v_c = v(3:4,4); v_c = v_c/norm(v_c);
         u_c = v(1:2,4); u_c = u_c/norm(u_c);  
       else
         v_c = v_rf; u_c = u_rf;
       end    
       rho_c = u_c'*(B*v_c);
       norm_def = norm([UTAP1 UTAP2]*v_c)^2+norm([VTATQ1 VTATQ2]*u_c)^2;
       norm_res_def = norm(B*v_c - rho_c*u_c)^2 + norm(B'*u_c - rho_c*v_c)^2 + (u_c(2)*fnorm)^2;
       norm_res = sqrt(norm_res_def+norm_def); Smax_all = max([Smax_all,Smax,S(1,1),rho_c]);
       if sqrt(norm_def) > tol*Smax_all && sqrt(norm_res_def) < tol*Smax_all
          LockingProblem = 'TRUE'; convergedKVals='FALSE'; break; % Locking problem 
       end
       if norm_res < tol*Smax_all &&  abs(v_b'*v_c) > 0.9 
          U(:,iter_k) = Q*u_c; V(:,iter_k)= P*v_c; S(iter_k,iter_k) = rho_c;
          maxnormres = max(maxnormres,norm_res); break; % Converge exit
       end 
      
       if iter_refined 
           
          % Count number of iter. refined restarts
          numIterRef = numIterRef +1;
         
          % Restart with iterative refined Ritz value/vectors
          f = P*(B'*u_rf - alpha*v_rf) + f*u_rf(2); B=[]; B(1,1) = alpha; 
          P(:,1) = P*v_rf; Q(:,1) = Q*u_rf; f = f - (P(:,1)'*f)*P(:,1); 
          B(1,2) = norm(f); P(:,2) = f*(1/B(1,2)); v_p = v_rf; u_p=u_rf;
         
       else
      
         % Restart with Ritz value/vectors 
         P = [P*v_b f*(1/fnorm)]; Q(:,1) = Q*u_b; 
         B(1,1) = s_b; B(1,2) = fnorm*u_b(2); v_p = v_b; u_p=u_b;
        
       end
  
       % Matrix-vector product, orthogonalization, and deflation to get Q(:,2)
       Astart = tic; if ischar(A), Q(:,2) = feval(A,P(:,2),m,n,'F'); else,... 
       Q(:,2) = A*P(:,2); end; mprod = mprod + 1; timeMatProds = timeMatProds + toc(Astart);
       if iter_k > 1 % Deflation step
          UTAP1 = [UTAP1 UTAP2]*v_p; UTAP2 = U(:,1:iter_k-1)'*Q(:,2);
          Q(:,2) = Q(:,2) - U(:,1:iter_k-1)*UTAP2; 
       end  
       Q(:,2) = Q(:,2) - Q(:,1)*(Q(:,2)'*Q(:,1));   
       B(2,2) = norm(Q(:,2)); Q(:,2) = Q(:,2)*(1/B(2,2));
  
      % Matrix-vector product, orthogonalization, and deflation to get f
      Astart = tic; if ischar(A), f = feval(A,Q(:,2),m,n,'T'); else, ...
      f = A'*Q(:,2); end; mprod = mprod + 1;timeMatProds = timeMatProds + toc(Astart);
      if iter_k > 1 % Deflation step
         VTATQ1 = [VTATQ1 VTATQ2]*u_p; VTATQ2 = V(:,1:iter_k-1)'*f;
         f = f - V(:,1:iter_k-1)*VTATQ2; 
      end 
      f = f - P(:,2)*B(2,2); fnorm = norm(f); 
  
      iter = iter+1;
 
   end % end: iter while loop
   
   if strcmp(LockingProblem,'TRUE'), STATS.LockingProblem = LockingProblem; break; end
   if iter >= maxit_outer, convergedKVals='FALSE'; iter_k = iter_k-1; break; end 
   
end % end: iter_k for loop

if iter_k >= 1 
    U = U(:,1:iter_k); V = V(:,1:iter_k); S = S(1:iter_k,1:iter_k);
else
    U=[]; V=[]; S=[]; 
end
           
%
STATS.numMatProds = mprod;
STATS.timeMatProds = timeMatProds;
STATS.numIterRefRestart = numIterRef;
STATS.timeCompIterRef = timeIterRef;
STATS.timeTotal = toc(tStart);
STATS.estimateSVmax = Smax_all;
STATS.convergedKVals = convergedKVals;
STATS.maxnormres = maxnormres;
end % end: function