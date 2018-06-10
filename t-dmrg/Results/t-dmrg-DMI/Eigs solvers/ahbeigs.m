function varargout = ahbeigs(varargin)
%
% AHBEIGS: will find a few eigenvalues and eigenvectors for either the 
% standard eigenvalue problem A*x = lambda*x or the generalized eigenvalue 
% problem A*x = lambda*B*x. 
%
% [V,D,FLAG] = AHBEIGS(A,OPTIONS)
% [V,D,FLAG] = AHBEIGS(A,B,OPTIONS)
% [V,D,FLAG] = AHBEIGS('Afunc',N,B,OPTIONS)
%
% The first input argument must be a matrix A which can be passed as a numeric 
% matrix or as a M-file ('Afunc') that  computes the product A*X or inv(B)*A*X, 
% where X is a (N x blsz) matrix. If A is passed as a M-file then the second 
% input argument N is the size of the matrix A. For the generalized eigenvalue 
% problem the matrix B is passed only as a numeric matrix. If a linear solver is
% required then use the M-file ('Afunc') and compute inv(B)*A*X. In all the 
% implementations AHBEIGS(A,...) can be replaced with AHBEIGS('Afunc',N,...).
%
% The function 'Afunc' must be a M-file. The (N x blsz) matrix X is the first 
% input argument, the second input argument is N, (the size of the matrix A), 
% and the third input argument is blsz, i.e. Afunc(X,N,blsz). The output is
% the (N x blsz) matrix Y where, Y = A*X or Y = inv(B)*A*X.
%
% OUTPUT OPTIONS:
% ---------------
%
% I.)   AHBEIGS(A) or AHBEIGS(A,B)
%       Displays the desired eigenvalues.    
%
% II.)  D = AHBEIGS(A) or D = AHBEIGS(A,B)
%       Returns the desired eigenvalues in the vector D. 
%
% III.) [V,D] = AHBEIGS(A) or [V,D] = AHBEIGS(A,B)
%       D is a diagonal matrix that contains the desired eigenvalues along the 
%       diagonal and the matrix V contains the corresponding eigenvectors, such 
%       that A*V = V*D or A*V = B*V*D. If AHBEIGS reaches the maximum number of
%       iterations before convergence then the matrices V and D contain any 
%       eigenpairs that have converged plus the last Ritz pair approximations 
%       for the eigenpairs that have not converged. 
%
% IV.)  [V,D,FLAG] = AHBEIGS(A) or [V,D,FLAG] = AHBEIGS(A,B)
%       This option returns the same as (III) plus a two dimensional array FLAG 
%       that reports if the algorithm converges and the number of matrix vector 
%       products. If FLAG(1) = 0 this implies normal return, all eigenvalues have 
%       converged. If FLAG(1) = 1 then the maximum number of iterations have been 
%       reached before all desired eigenvalues have converged. FLAG(2) contains the 
%       number of matrix vector products used by the code. If the maximum number of 
%       iterations are reached then the matrices V and D contain any eigenpairs 
%       that have converged plus the last Ritz pair approximations for the eigenpairs
%       that have not converged.
%
% INPUT OPTIONS:
% --------------
%                                   
%       ... = AHBEIGS(A,OPTS) or  ... = AHBEIGS(A,B,OPTS)
%       OPTS is a structure containing input parameters. The input parameters can
%       be given in any order. The structure OPTS may contain some or all of the 
%       following input parameters. The string for the input parameters can contain
%       upper or lower case characters. 
%       
%  INPUT PARAMETER      DESCRIPTION     
%
%  OPTS.ADJUST       Initial number of vectors to add to the K restart vectors. After
%                    vectors start to converge more vectors are added to help increase 
%                    convergence. Should be set so that K + ADJUST is a multiple of BLSZ.
%                    DEFAULT VALUE    OPTS.ADJUST = 3                
%   
%  OPTS.BLSZ         Block size of the Arnoldi Hessenberg matrix.        
%                    DEFAULT VALUE    BLSZ = 3
%
%  OPTS.CHOLB        Indicates if the Cholesky factorization of the matrix B is available.
%                    If the Cholesky factorization matrix R is available then set 
%                    CHOLB = 1 and replace the input matrix B with R where B = R'*R.
%                    DEFAULT VALUE   CHOLB = 0
%
%  OPTS.DISPR        Indicates if K Ritz values and residuals are to be displayed on each 
%                    iteration. Set positive to display the Ritz values and residuals on 
%                    each iteration.
%                    DEFAULT VALUE   DISPR = 0 
%
%  OPTS.PERMB        Permutation vector for the Cholesky factorization of B(PERMB,PERMB). 
%                    When the input matrix B is replaced with R where B(PERMB,PERMB)=R'*R
%                    then the vector PERMB is the permutation vector. 
%                    DEFAULT VALUE   PERMB=1:N                     
%                     
%  OPTS.K            Number of desired eigenvalues.             
%                    DEFAULT VALUE  K = 6   
%
%  OPTS.MAXIT        Maximum number of iterations, i.e. maximum number of restarts.                           
%                    DEFAULT VALUE   MAXIT = 100
%
%  OPTS.NBLS         Number of blocks in the Arnoldi Hessenberg matrix. If number of
%                    blocks is not sufficiently large enough then AHBEIGS will not
%                    converge or miss some desired eigenvalues.                            
%                    DEFAULT VALUE    NBLS = 10
%
%  OPTS.SIGMA        Two letter string or numeric value specifying the location 
%                    of the desired eigenvalues.
%                     'LM' or 'SM' Largest or Smallest Magnitude
%                     'LR' or 'SR' Largest or Smallest Real part
%                     'LI' or 'SI' Largest or Smallest Imaginary part
%                     'LA' or 'SA' Largest or Smallest Algebraic (Symmetric Problems only)
%                     NVAL  A numeric value. The program searches for the K closest
%                           eigenvalues to the numeric value NVAL. (Program will factor
%                           the matrix. This will increase storage requirements.)
%                    DEFAULT VALUE   SIGMA = 'LM'
%                                                                                    
%  OPTS.TOL          Tolerance used for convergence. Convergence is determined when             
%                    || Ax - lambda*x ||_2 <= TOL*max(abs(Ritz)). 
%                    DEFAULT VALUE    TOL = 1d-6 
%                                                              
%  OPTS.V0           A matrix of starting vectors.       
%                    DEFAULT VALUE  V0 = randn
% 
%  DATE:  6/15/2005
%  VER:  1.0

%  AUTHOR:
%  James Baglama     
%  Department of Mathematics 
%  University of Rhode Island
%  Kingston, Rhode Island 02881 
%  E-Mail: jbaglama@math.uri.edu 
%
%
% Research supported by NSF grant DMS-0311786. 
%

% Too many output arguments requested.
if (nargout >= 4) error('ERROR: Too many output arguments.'); end

%----------------------------%
% BEGIN: PARSE INPUT VALUES. %
%----------------------------%

% No input arguments, return help.
if nargin == 0, help ahbeigs, return, end

% Get matrix A. Check type (numeric or character) and dimensions.
if (isstruct(varargin{1})), error('A must be a matrix.'); end
A = varargin{1};
if ischar(varargin{1})
   if nargin == 1, error('Need N (size of matrix A).'); end  
   n = varargin{2};
   if ~isnumeric(n) | length(n) > 2 
      error('Second argument N must be a numeric value.'); 
   end
else
   [n,n] = size(varargin{1});
   if any(size(varargin{1}) ~= n), error('Matrix A is not square.');      end
   if ~isnumeric(varargin{1}),     error('A must be a numeric matrix.');  end
end

% Set all input options to default values or set to empty arrays.
blsz = 3; adjust = []; cholB = 0;  dispr = 0;  K = 6; B=[]; 
maxit = 100; nbls = 10; permB = []; sigma = 'LM'; tol = 1d-6; V=[]; 
L=[]; U=[]; P=[]; nval=[];

% If a generalized eigenvalue problem is to be solved get matrix B. 
% Check type (numeric or character) and dimensions.
if ((ischar(varargin{1}) & nargin > 2) | (isnumeric(varargin{1}) & nargin > 1)) 
   if ~isstruct(varargin{2+ischar(varargin{1})})
      B = varargin{2+ischar(varargin{1})};
      if ~isnumeric(B),     error('B must be a numeric matrix.');          end   
      if any(size(B) ~= n), error('Matrix B must be the same size as A.'); end
   end
end   

% Get input options from the structure array.
if nargin > 1 + ischar(varargin{1}) + ~isempty(B)
   options = varargin{2+ischar(varargin{1}) + ~isempty(B):nargin};
   names = fieldnames(options);
   I = strmatch('ADJUST',upper(names),'exact');
   if ~isempty(I), adjust = getfield(options,names{I}); end
   I = strmatch('BLSZ',upper(names),'exact');
   if ~isempty(I), blsz = getfield(options,names{I}); end
   I = strmatch('CHOLB',upper(names),'exact');
   if ~isempty(I), cholB = getfield(options,names{I}); end
   I = strmatch('DISPR',upper(names),'exact');
   if ~isempty(I), dispr = getfield(options,names{I}); end
   I = strmatch('K',upper(names),'exact');
   if  ~isempty(I), K = getfield(options,names{I}); end
   I = strmatch('MAXIT',upper(names),'exact');
   if ~isempty(I), maxit = getfield(options,names{I}); end
   I = strmatch('NBLS',upper(names),'exact');
   if ~isempty(I), nbls = getfield(options,names{I}); end
   I = strmatch('PERMB',upper(names),'exact');
   if ~isempty(I), permB = getfield(options,names{I}); end
   I = strmatch('SIGMA',upper(names),'exact');
   if  ~isempty(I), sigma = upper(getfield(options,names{I})); end
   I = strmatch('TOL',upper(names),'exact');
   if ~isempty(I), tol = getfield(options,names{I}); end
   I = strmatch('V0',upper(names),'exact');
   if ~isempty(I), V = getfield(options,names{I}); end
end 

% Check type of input values and output error message if needed.
if (~isnumeric(blsz)   | ~isnumeric(K)     | ~isnumeric(dispr) | ...
    ~isnumeric(maxit)  | ~isnumeric(nbls)  | ~isnumeric(tol)   | ...
    ~isnumeric(permB)  | ~isnumeric(cholB))
   error('Incorrect type for input value(s) in the structure.');
end

% If a numeric value is given for sigma then set nval=sigma and 
% sigma = 'LM' to denote that the code is searching eigenvalues near nval.
% Using the transformation inv(A-sigma*B)*B x = theta * x where 
% theta = 1/(lambda-sigma).
if isnumeric(sigma), nval = sigma; sigma = 'LM'; end

% Check the length and value of the character sigma.
if length(sigma) ~= 2, error('SIGMA must be SA, LA, SM, LM, SR, LR, LI, or SI.'), end
if (~strcmp(sigma,'SA') & ~strcmp(sigma,'LA') & ~strcmp(sigma,'LI') & ...
    ~strcmp(sigma,'SM') & ~strcmp(sigma,'LM') & ~strcmp(sigma,'SI') & ...
    ~strcmp(sigma,'SR') & ~strcmp(sigma,'LR'))
   error ('SIGMA must be SA, LA, SM, LM, SR, LR LI, or SI.'); 
end

% Resize Krylov subspace if blsz*nbls (i.e. number of Arnoldi vectors) is larger 
% than n (i.e. the size of the matrix A).
if blsz*nbls >= n 
   nbls = floor(n/blsz); 
   warning(['Changing NBLS to ',num2str(nbls)]);
end

% Increase the number of desired values to help increase convergence.
% Set K+adjust to be next multiple of BLSZ > K. 
if isempty(adjust)
   adjust = 1:blsz; I = find(mod(K + adjust-1,blsz) == 0); adjust = adjust(I(1)); 
end 
K_org = K; K = K + adjust; % K_org to be orginal value of K.
 
% Check for input errors in the structure array.
if K <= 0,    error('K must be a positive value.');        end
if K > n,     error('K is too large. Use eig(full(A))! '); end   
if blsz <= 0, error('BLSZ must be a positive value.');     end
if nbls <=1,  error('NBLS must be greater than 1.');       end
if tol < 0,   error('TOL must be non-negative.');          end
if maxit <=0, error('MAXIT must be positive.');            end
if ((cholB ~= 0) & (cholB ~= 1)), error('Unknown value for CHOLB.'); end
if ~isempty(permB)
   if ((size(permB,1) ~= n) | (size(permB,2) ~= 1)) & ...
      ((size(permB,1) ~= 1) | (size(permB,2) ~= n))
      error('Incorrect size for PERMB.'); 
   end
end

% Automatically adjust Krylov subspace to accommodate larger values of K.
if blsz*(nbls-1) - blsz - K - 1 < 0
   nbls = ceil((K+1)/blsz+2.1);
   warning(['Changing NBLS to ',num2str(nbls)]);
end
if blsz*nbls >= n, 
   nbls = floor(n/blsz-0.1);
   warning(['Changing NBLS to ',num2str(nbls)]); 
end
if blsz*(nbls-1) - blsz - K - 1 < 0,  error('K is too large.'); end

% If starting matrix V0 is not given then set starting matrix V0 to be a 
% (n x blsz) matrix of normally distributed random numbers. Also,
% preallocate memory size, so no resizing of the large matrix V occurs.
if nnz(V) == 0
   V = [randn(n,blsz) zeros(n,nbls*blsz)]; 
else
   if ~isnumeric(V),     error('Incorrect starting matrix V0.');         end
   if (size(V,1) ~= n),  error('Incorrect size of starting matrix V0.'); end
   if nnz(V) == 0,       error('Matrix V contains all zeros.');          end
   V = [V zeros(n,nbls*blsz)];
end 

% Set tolerance to machine precision if tol < eps.
if tol < eps, tol = eps; end

% Using shift and invert if sigma is given as a real number.
if ~isempty(nval)
   if ~isnumeric(varargin{1}) 
      error('A must be a numeric matrix in order to use the internal shift and invert.');          
   end
   if ~isempty(B) 
      [L,U,P] = lu(varargin{1}-nval*B);
   else
      [L,U,P] = lu(varargin{1}-nval*eye(n));
   end            
end   
    
% If a generalized eigenvalue problem is to be solved then try to get the Cholesky 
% factorization of the matrix B (i.e. B=R'*R) and set B=R the upper triangular matrix.
if ~isempty(B) & isempty(nval)
    
   % If B is sparse then compute the sparse cholesky factorization. Using symmmd or 
   % symamd to get a permutation vector permB such that B(permB,permB) will have a 
   % sparser cholesky factorization than B. The permutation is given by
   % B(permB,permB) = I(:,permB)'*B*I(:,permB), where I is the identity matrix. 
   % To compute the product A(permB,permB)*X use Y(permB,:) = X; Y=A*Y; Y=Y(permB,:);      
   if isempty(permB), permB = 1:n; end
 
   % If cholB = 1, then the input matrix B is given as R in the Cholesky factorization.
   if ~cholB
      if ~isreal(B), error('Matrix B must be real.'); end
      [L,cholerr] = chol(B(permB,permB));
      
      % Will try to solve problem using the inverse of B, (i.e. inv(B)*A*x = lambda*x)
      % by computing the LU factorization of B. This may cause numerical instabilities.
      if cholerr 
         [L,U,P] = lu(B); % Used only if B is not symmetric positive definite. 
      else                 
         B = L; L=[];  
      end  
   end
end   

%--------------------------%
% END: PARSE INPUT VALUES. %
%--------------------------%

%-----------------------------------------------------------%
% BEGIN: DESCRIPTION AND INITIALIZATION OF LOCAL VARIABLES. %
%-----------------------------------------------------------%
 
check_norm_balance=[]; % Accuray of eigenpairs of H. Determines if nobalance in eig is needed.
conv = 'F';            % Boolean to determine if all desired eigenvalues have converged.
H = [];                % Block Hessenberg (Hsz_n x Hsz_m) matrix. 
H_R = [];              % Stores the last block of H. Used in updating H.
H_eig = [];            % Eigenvalues of H.
H_eig_no = [];         % Eigenvalues of H using nobalance in eig. Not always computed.
H_eigv = [];           % Eigenvectors of H.
H_eigv_no = [];        % Eigenvectors of H using nobalance in eig. Not always computed.
Hsz_m = [];            % Second dimension of H.
Hsz_n = [];            % First dimension of H.
iter = 1;              % Main loop iteration count.
m = blsz;              % Current number of columns in the matrix V.
mprod = 0;             % The number of matrix vector products.
Q = [];                % Schur vectors for H.
R = [];                % Is a diagonal matrix of 1 or -1's used in the function setuparnv.
sqrteps = sqrt(eps);   % Square root of machine tolerance.
T = [];                % T matrix in the WY representation of the Householder products.
T_schur=[];            % Schur matrix for H.
Tsz_m = [];            % Size of T.

%---------------------------------------------------------%
% END: DESCRIPTION AND INITIALIZATION OF LOCAL VARIABLES. %
%---------------------------------------------------------%

%-----------------------------%
% BEGIN: MAIN ITERATION LOOP. %
%-----------------------------%

while (iter <= maxit)
 
    % Compute the block Householder Arnoldi decomposition.
    [V,H,T,mprod] = ArnHouBlk(varargin{1},B,L,U,P,nval,permB,V,H,T,n,m,nbls,blsz,mprod);
    
    % Determine the size of the block Hessenberg matrix H. Possible
    % truncation may occur in ArnHouBlk if an invariant subspace
    % has been found or K + ADJUST is not a multiple of BLSZ.
    Hsz_n = size(H,1); Hsz_m = size(H,2);

    % Compute the eigenvalue decomposition of the block Hessenberg H(1:Hsz_m,:).
    [H_eigv,H_eig] = eig(H(1:Hsz_m,:));
    
    % Check the accuracy of the computation of the eigenvalues of the
    % Hessenberg matrix. This is used to monitor balancing.
    check_norm_balance = norm(H(1:Hsz_m,:)*H_eigv - H_eigv*H_eig);
    if check_norm_balance >= eps*norm(H(1:Hsz_m,:))
       [H_eigv_no,H_eig_no] = eig(H(1:Hsz_m,:),'nobalance');
       check_norm_nobalance = norm(H(1:Hsz_m,:)*H_eigv_no - H_eigv_no*H_eig_no);
       if check_norm_nobalance < check_norm_balance
          H_eigv = H_eigv_no; H_eig = H_eig_no;
       end   
    end
    H_eig = diag(H_eig);
    
    % Sort the eigenvalues and check for convergence.
    [H_eig,H_eigv,K,conv] = eigsortconv(Hsz_m,Hsz_n,blsz,dispr,H,iter,K_org,...
                                        K,adjust,H_eig,H_eigv,sigma,tol);

    % If all desired Ritz values converged then exit main loop.
    if strcmp(conv,'T'), break; end
    
    % Update the main iteration loop count.
    iter = iter + 1; if iter >= maxit, break; end

    % Convert the complex eigenvectors of the eigenvalue decomposition of H
    % to real vectors and convert the complex diagonal matrix to block diagonal.
    % Use MATLAB's internal function cdf2rdf to accomplish this task.
    if ~isreal(H_eig)
        [Q,T_Schur] = cdf2rdf(H_eigv,diag(H_eig)); % Convert complex vectors to real vectors.
        [Q,R] = qr(Q(:,1:K),0);       % Compute the QR factorization of H_eigv(:,1:K).
    else
        [Q,R] = qr(H_eigv(:,1:K),0);  % Compute the QR factorization of H_eigv(:,1:K). 
    end
    
    % The Schur matrix for H.
    T_Schur = triu(Q'*H(1:Hsz_m,:)*Q,-1);
    
    % Compute the starting vectors and the residual vectors from the Householder
    % WY form. The starting vectors will be the first K Schur vectors and the 
    % residual vectors are stored as the last BLSZ vectors in the Householder WY form.  
    Tsz_m = size(T,1);
    V(:,Hsz_n-blsz+1:Hsz_n)= V(:,1:Tsz_m)*(T*V(Tsz_m-blsz+1:Tsz_m,1:Tsz_m)'); 
    V(Tsz_m-blsz+1:Tsz_m,Hsz_n-blsz+1:Hsz_n) = eye(blsz,blsz)+V(Tsz_m-blsz+1:Tsz_m,Hsz_n-blsz+1:Hsz_n);
    %
    V(:,1:K) = V(:,1:Hsz_m)*(T(1:Hsz_m,1:Hsz_m)*(V(1:Hsz_m,1:Hsz_m)'*Q(:,1:K)));
    V(1:Hsz_m,1:K) = Q(:,1:K) + V(1:Hsz_m,1:K);
    
    % Set the size of the large matrix V and move the residual vectors.
    m = K + 2*blsz; V(:,K+1:K+blsz) = V(:,Hsz_n-blsz+1:Hsz_n);
    
    % Set the new starting vector(s) to be the desired vectors V(:,1:K) with the 
    % residual vectors V(:,Hsz_n-blsz+1:Hsz_n). Place all vectors in the compact 
    % WY form of the householder product. Compute the next set of vectors by 
    % computing A*V(:,Hsz_n-blsz+1:Hsz_n) and store this in V(:,Hsz_n+1:Hsz_n+blsz).
    [V,T,R,mprod] = SetupArnVecs(V,T,varargin{1},B,L,U,P,nval,permB,n,blsz,m-blsz,mprod);

    % Compute the first K columns and K+BLSZ rows of the matrix H, used in augmenting.
    H_R = H(Hsz_n-blsz+1:Hsz_n,Hsz_m-blsz+1:Hsz_m); H=zeros(K+blsz,K);
    H(1:K,1:K) = R(1:K,1:K)*T_Schur(1:K,1:K)*R(1:K,1:K);
    H(K+1:K+blsz,1:K) = R(K+1:K+blsz,K+1:K+blsz)*H_R*Q(Hsz_m-(blsz-1):Hsz_m,1:K)*R(1:K,1:K);
    
end 

%---------------------------%
% END: MAIN ITERATION LOOP. %
%---------------------------%

%------------------------%
% BEGIN: OUTPUT RESULTS. %
%------------------------%

% Test to see if maximum number of iterations have been reached.
FLAG(1) = 0; if iter >= maxit FLAG(1) = 1; end; 
if FLAG(1)==1, warning('AHBEIGS did not find all requested eigenvalues!'); end

% Truncated eigenvalue and eigenvector arrays to include only desired eigenpairs.
H_eig = H_eig(1:K_org); H_eigv = H_eigv(:,1:K_org);

% Convert the shift and invert eigenvalues back to orginal state using 1/(lambda-sigma).
if ~isempty(nval), H_eig = 1./H_eig + nval; mprod = []; end

% Sort output.
[sortval,J] = sort(abs(H_eig(1:K_org))); J = J(length(J):-1:1);

% Output option I: Display eigenvalues only.
if (nargout == 0), eigenvalues = H_eig(J), end
 
% Output option II: Set eigenvalues equal to output vector.
if (nargout == 1), varargout{1} = H_eig(J); end
 
% Output option III. Output diagonal matrix of eigenvalues and
% corresponding matrix of eigenvectors.
if nargout == 2 | nargout == 3
   H_eigv = H_eigv(:,1:K_org); H_eigv = H_eigv(:,J);
   Tsz_m = size(T,1); H_sz1 = size(H_eigv,1); H_sz2 = size(H_eigv,2);
   V = V(:,1:Tsz_m)*(T*(V(1:H_sz1,1:Tsz_m)'*H_eigv));
   V(1:H_sz1,1:H_sz2) = H_eigv + V(1:H_sz1,1:H_sz2);   

   % Must solve a linear system to extract generalized eigenvectors.  
   if ~isempty(B) & isempty(L) & isempty(nval), V = B\V(permB,:); end  
   varargout{1} = V; varargout{2} = diag(H_eig(J));
end

% Output option IV: Output FLAG. 
if nargout == 3, FLAG(2)=mprod; varargout{3}=FLAG; end

%----------------------%
% END: OUTPUT RESULTS. %
%----------------------%

%---------------------------------------------%
% BEGIN: AUGMENTED BLOCK HOUSEHOLDER ARNOLDI. %
%---------------------------------------------%

function [V,H,T,mprod] = ArnHouBlk(A,B,L,U,P,nval,permB,V,H,T,n,m,nbls,blsz,mprod)
% Computes the  Block Arnoldi decomposition using Householder reflections.
% 
% The matrix A can be passed as a numeric matrix or as a filename. Note 
% that if the matrix A is a filename then the file must accept [X,N,BLSZ] 
% as input in that order and the file must return the matrix-vector product A*X.
%
% Input:
%       A      - Matrix A.
%       B      - Matrix B for generalized eigenvalue problems.
%       L      - LU factorization of (A-nval*B) or  B.
%       U      - LU factorization of (A-nval*B) or  B.
%       P      - Permutation matrix for the LU factorization.
%     NVAL     - Numeric value of sigma, if given. 
%    PERMB     - Permutation vector for  generalized eigenvalue problems.           
%       V      - (N x M) Initial matrix to start the Arnoldi reduction.
%       H      - Empty or contains Schur values and residuals values on restart.
%       T      - Upper triangular matrix used in Householder WY representation.
%       N      - Size of the matrix A.
%       M      - Current size of V.
%    NBLS      - Maximum number of Arnoldi vectors.
%    BLSZ      - Block size.
%    MPROD     - Number matrix-vector products.
%
% Output:
%       V      - N x BLSZ*NBLS+BLSZ matrix in WY-Householder format.
%       H      - NBLS*BLSZ+BLSZ x NBLS*BLSZ upper block Hessenberg matrix.
%       T      - Upper triangular matrix used in Householder WY representation.
%       MPROD  - Number of matrix-vector products.

% James Baglama
% DATE: 7/20/2005

% LOCAL VARIABLES
T_blk =zeros(blsz,blsz); % Small BLSZ x BLSZ upper triangular matrix used in Householder WY representation.

% Check size of the block Arnoldi matrix H. If size is greater than n
% reduce size.
if blsz*(nbls+1) > n, nbls = floor(n/blsz)-1; end
 
 % Begin of main iteration loop for the Augmented Block Householder Arnoldi decomposition.
 while (m <= blsz*(nbls+1))
     
     % Compute Householder vectors for BLSZ vectors V(m-blsz+1:n,m-blsz+1:m) using 
     % MATLAB's internal qr algorithm. 
      V(m-blsz+1:n,m-blsz+1:m) = qr(V(m-blsz+1:n,m-blsz+1:m));

     % Compute the columns of the upper block hessenberg matrix H.
      if m > blsz
         H(1:m-blsz,m-2*blsz+1:m-blsz) = V(1:m-blsz,m-blsz+1:m);
         H(m-blsz+1:m,m-2*blsz+1:m-blsz) = triu(V(m-blsz+1:m,m-blsz+1:m));
      end 

     % Set up householder vectors.
     m2 = 2*m; if m2 > n, m2=n; end
     V(1:m2,m-blsz+1:m) = tril(V(1:m2,m-blsz+1:m),-(m-blsz+1)); 
     for i=m-blsz+1:m,V(i,i)=1; end
          
     % Using the compact WY form of the householder product. See,
     % "A storage-efficient WY representation for products of
     % householder transformations", by Schreiber and Van Loan,
     % SIAM J. Sci. Stat. Comput. Vol. 10, No. 1 (1989) pp. 53-57.
     % P_{1} * P_{2} * ... * P_{n} = I + V *T *V';
    
     % Place BLSZ householder vectors into I + V *T_blk *V' form.
     T_blk(1,1) = -2/(V(:,m-blsz+1)'*V(:,m-blsz+1));
     for i=2:blsz
        T_blk(1:i,i) = (V(:,m-blsz+i)'*V(:,m-blsz+1:m-blsz+i))';
        T_blk(i,i) = -2/T_blk(i,i);
        T_blk(1:i-1,i) = T_blk(i,i)*T_blk(1:i-1,1:i-1)*T_blk(1:i-1,i);
     end   
     
     % Update I + V *T *V'.
     T = [ [T; zeros(blsz,m-blsz)] [(T*(V(:,m-blsz+1:m)'*V(:,1:m-blsz))')*T_blk; T_blk]];
 
     if m <= blsz*nbls 
       
       % V_{m+blsz} = ( I + V_m * T_m * V_m') * E_m.
       V(:,m+1:m+blsz) = V(:,1:m)*(T*V(m-blsz+1:m,1:m)'); 
       V(m-blsz+1:m,m+1:m+blsz) = eye(blsz,blsz)+V(m-blsz+1:m,m+1:m+blsz);    

       % Matrix-vector product, V_{m+blsz} = A * V_{m+blsz}.
       [V(:,m+1:m+blsz),mprod] = matrixprod(A,B,L,U,P,nval,permB,V(:,m+1:m+blsz),n,blsz,mprod);
               
       % V_{m+blsz} = (I + V_m * T_m' * V_m')*V_{m+blsz}.
       V(:,m+1:m+blsz) = V(:,m+1:m+blsz) + V(:,1:m)*((V(:,m+1:m+blsz)'*V(:,1:m))*T)';
    
     end   
        
     % Update iteration count.
     m = m + blsz;
  
end % End main loop.  

%-------------------------------------------%
% END: AUGMENTED BLOCK HOUSEHOLDER ARNOLDI. %
%-------------------------------------------%

%--------------------------------------------------------%
% BEGIN: EIGENVALUE REORDERING AND CONVERGENCE TESTING.  %
%--------------------------------------------------------%

function [H_eig,H_eigv,K,conv] = eigsortconv(Hsz_m,Hsz_n,blsz,dispr,H,iter,K_org,K,...
                                             adjust,H_eig,H_eigv,sigma,tol);
% This reorders the eigenvalues of H and determines convergence.
%
% Input:
%       HSZ_M  - Second dimension of H.
%       HSZ_N  - First dimension of H.
%        BLSZ  - Block size of H.
%       DISPR  - Boolean to determine if Ritz pairs should be displayed.
%           H  - Block Hessenberg (Hsz_n x Hsz_m) matrix.  
%        ITER  - Number of iterations of the main loop.
%       K_ORG  - Orginal number of desired eigenpairs.
%           K  - Adjusted number of eigenpairs.
%      ADJUST  - Orginal value added to K_ORG
%       H_EIG  - Eigenvalues of H.
%      H_EIGV  - Eigenvectors of H.
%       SIGMA  - Two letter string specifying the location of the desired eigenvalues.
%         TOL  - Tolerance used for convergence. Convergence is determined when             
%                || Ax - lambda_i*x ||_2 <= TOL*abs(Ritz_i).
%
% Output:
%       H_EIG  - Sorted eigenvalues of H. 
%      H_EIGV  - Sorted eigenvectors of H.   
%           K  - Adjusted number of eigenpairs. Increased if eigenvectors converge.
%        CONV  - Boolean to indicate all eigenpairs have converged.

% James Baglama
% DATE: 7/20/2005

% LOCAL VARIABLES
conv = 'F';  % Boolean to determine if all desired eigenpairs have converged.
eig_conj=[]; % Indicator for when K splits a conjugate pair. K is increased by 1.
I_eig = [];  % Integer vector of indexes for the ordering of the eigenvalues.
Len_res=[];  % Number of converged values.


% If sigma = LI or SI and all eigenvalues are real then sort by the real
% value.
if isreal(H_eig) & (strcmp(sigma,'SI') | strcmp(sigma,'LI'))
    if strcmp(sigma,'SI'), sigma = 'SR'; end
    if strcmp(sigma,'LI'), sigma = 'LR'; end
end    

% Determine the ordering of the eigenvalues.
switch sigma

    % Eigenvalues of largest magnitude.
    case 'LM', [Y,I_eig] = sort(abs(H_eig)); I_eig = I_eig(length(I_eig):-1:1);
    
    % Eigenvalues of smallest magnitude. 
    case 'SM', [Y,I_eig] = sort(abs(H_eig));
  
    % Eigenvalues of largest real part.
    case 'LR', [Y,I_eig] = sort(real(H_eig)); I_eig = I_eig(length(I_eig):-1:1);

    % Eigenvalues of smallest real part.
    case 'SR', [Y,I_eig] = sort(real(H_eig));

    % Eigenvalues of largest imaginary part.
    case 'LI', [Y,I_eig] = sort((abs(imag(H_eig)))); I_eig = I_eig(length(I_eig):-1:1);

    % Eigenvalues of smallest imaginary part.
    case 'SI', [Y,I_eig] = sort((abs(imag(H_eig))));

    % Largest algebraic eigenvalues (symmetric problems only).
    case 'LA', [Y,I_eig] = sort(H_eig); I_eig = I_eig(length(I_eig):-1:1);
  
    % Smallest algebraic eigenvalues (symmetric problems only).
    case 'SA', [Y,I_eig] = sort(H_eig);

end

% Sort the eigenvalue and eigenvector arrays.
H_eig = H_eig(I_eig); H_eigv = H_eigv(:,I_eig);

% Compute the residuals for the K_org Ritz values.
residuals = (sum(abs(H(Hsz_n-blsz+1:Hsz_n,Hsz_m-blsz+1:Hsz_m)*...
                          H_eigv(Hsz_m-(blsz-1):Hsz_m,1:K_org)).^2, 1).^(0.5));

% Output intermediate results.
 if dispr ~= 0 
   format short e
   % Sort output values.
   [sortval,J] = sort(abs(H_eig(1:K_org))); J = J(length(J):-1:1);
   S = [H_eig(J) residuals(J)'];
   if ~isreal(S)
      disp(sprintf('  Ritz values                   Residual       Iteration: %d',iter));
   else
      disp(sprintf('  Ritz values     Residual       Iteration: %d',iter));
   end
   disp(S); disp(' '); disp(' ');
 end

% Check for convergence.
Len_res = length(find(residuals(1:K_org)' < tol*abs(H_eig(1:K_org)))); 
if Len_res >= K_org, conv = 'T'; end % Set convergence to true.
 
% Adjust K to include more vectors as the number of vectors converge.
Len_res = length(find(residuals(1:K_org)' < eps*abs(H_eig(1:K_org)))); 
K =  K_org + adjust+ Len_res; if K > Hsz_m - 2*blsz-1, K = Hsz_m - 2*blsz-1; end
   
% Determine if K splits a conjugate pair. If so replace K with K + 1.
if ~isreal(H_eig(K))
   eig_conj = 1;
   if K < Hsz_m
      if abs(imag(H_eig(K)) + imag(H_eig(K+1))) < sqrt(eps)
         K = K + 1; eig_conj = 0;
      end
   end   
   if K > 1 & eig_conj
      if abs(imag(H_eig(K)) + imag(H_eig(K-1))) < sqrt(eps)
         eig_conj = 0;
      end
   end   
   if eig_conj, error('%d th conjugate pair split.',K(i)); end
end  
 
%-----------------------------------------------------%
% END: EIGENVALUE REORDERING AND CONVERGENCE TESTING. %
%-----------------------------------------------------% 
  
%----------------------------------------------------------%
% BEGIN: WY HOUSEHOLDER REPRESENATION OF STARTING VECTORS. %
%----------------------------------------------------------%

function [V,T,R,mprod] = SetupArnVecs(V,T,A,B,L,U,P,nval,permB,n,blsz,m,mprod)
% Given an orthogonal matrix V, this function SetupArnVecs computes the 
% compact WY form of the householder product. See, "A storage-efficient 
% WY representation for products of householder transformations", by Schreiber 
% and Van Loan, SIAM J. Sci. Stat. Comput. Vol. 10, No. 1 (1989) pp. 53-57.
% Also to save on computations the next vector in the block Arnoldi process
% is also calculated.
%
% Input:  
%       V - (N x M) orthogonal Matrix.
%       T - Upper triangular matrix for WY-compact form.
%       A - Matrix A.
%       B - Matrix B for generalized eigenvalue problems.
%       L - LU factorization of (A-nval*B) or  B.
%       U - LU factorization of (A-nval*B) or  B.
%       P - Permutation matrix for the LU factorization.
%    NVAL - Numeric value of sigma, if given.        
%   PERMB - Permutation vector for  generalized eigenvalue problems.   
%       N - Size of the matrix A.
%    BLSZ - Block size
%       m - Current size of V.
%   MPROD - Number of matrix-vector products.
%
% Output:
%       V - (N x M+BLSZ) matrix such that (N x 1:M) are the Householder vectors
%           and (N x M+1:M+BLSZ) are the Arnoldi vectors needed to restart the 
%           Arnoldi process
%       T - (M x M) Upper triangular matrix such that P_{1} * ... * P_{n} = I + V*T *V'.
%       R - R is a diagonal matrix of 1's and -1's. R*R = I where
%           V(input) * R = (I + V * T * V') * I(:,m).
%   MPROD - Number of matrix-vector products.

% James Baglama
% DATE: 7/20/2005

% LOCAL VARIABLES
V_sign = []; % Used to determine Householder vectors, the sign of the diagonal value.
Vdot = [];   % Dot product of Householder vectors.
D = [];      % Used for scaling.

% Initialize matrix T.
T = V(1:m,1:m);

% Matrix-vector product, V_{m+1:m+blsz} = A * V_{m-blsz+1:m}.
[V(:,m+1:m+blsz),mprod] = matrixprod(A,B,L,U,P,nval,permB,V(:,m-blsz+1:m),n,blsz,mprod);

% Compute the householder vectors V. The input matrix V(:,1:m) is orthogonal
% and the length of each vector is 1. Hence the householder vector
% u_i = v_i + e_i*||v_i||_2 *sign(v_ii) = v_i + e_i*sign(v_ii) and the
% dot product u_i'*u_i = 2 + 2*sign(v_ii)*v_ii. Also, since the matrix is
% orthognal only the first m columns and rows are updated at this stage.
% At this stage V=[V_1; V_2], V_1 is a m x m matrix which is now computed.
for i =1:m
   V_sign(i) = sign(V(i,i)); if V_sign(i) == 0, V_sign(i)=1; end
    Vdot = 1 + V_sign(i)*V(i,i);      % Dot product of householder vectors. 
    V(i,i) = V(i,i) + V_sign(i);      % Reflection to the ith axis.   
    D(i,i) = 1/V(i,i);                % Used for scaling. Note: V(i,i) >= 1.

    % Apply the householder vectors to update the i+1:m vectors in V. Only
    % the first m rows are needed at this stage.
    V(1:m,i+1:m) = V(1:m,i+1:m) - V(1:m,i)*((V_sign(i)/Vdot)*V(i,i+1:m));  
end

% Scale the matrix V_1, so that diagonal elements are 1.
V(1:m,1:m) = V(1:m,1:m)*D;

% Compute the upper triangular matrix R such that the following relationship
% holds V(input) * R = (I + V * T * V') * I(:,m).
R = -(diag(V_sign));

% Update the new block.
V(:,m+1:m+blsz) = V(:,m+1:m+blsz)*(R(m-blsz+1:m,m-blsz+1:m));

% Set up householder vectors. Put exact zeros in place of computed zeros.
V(1:m,1:m) = tril(V(1:m,1:m),0);  

% Compute T. Note that V(input) * R = (I + V * T * V') * E_m. Let
% V = [V_1; V_2]. V_1 is already computed. Hence we have
% T = inv(V_1)*(E_m'*V_(input)*R - E_m)*inv(V_1)'. 
% V_1 is upper triangular matrix with 1's on the diagonal.
T = T*R - eye(m); T = T/V(1:m,1:m)';  T = triu(V(1:m,1:m)\T);

% Compute the V_2 n-m x m matrix in V = [V_1; V_2]. 
% V_2 = V(input)*[0; I(n-m,m)]*R*inv(T*V_1).
V(m+1:n,1:m) = V(m+1:n,1:m)*(R/(T*V(1:m,1:m)'));

% V_{m+blsz} = (I + V_m * T_m' * V_m')*V_{m+blsz}. Update the vectors to
% be used in the block Arnoldi subroutine.
V(:,m+1:m+blsz) = V(:,m+1:m+blsz) + V(:,1:m)*((V(:,m+1:m+blsz)'*V(:,1:m))*T)';

%---------------------------------------------------------%
% END: WY HOUSEHOLDER REPRESENATION OF STARTING VECTORS.  %
%---------------------------------------------------------%

%-------------------------------%
% BEGIN: MATRIX-VECTOR PRODUCT. %
%-------------------------------%

function [X,mprod] = matrixprod(A,B,L,U,P,nval,permB,X,n,blsz,mprod)
% Computes the matrix vector products.
%
% Input:
%       A  - Matrix A.
%       B  - Matrix B for generalized eigenvalue problems.
%       L  - LU factorization of (A-nval*B) or  B.
%       U  - LU factorization of (A-nval*B) or  B.
%       P  - Permutation matrix for the LU factorization.
%     NVAL - Numeric value of sigma, if given. 
%    PERMB - Permutation vector for  generalized eigenvalue problems. 
%       X  - (N x BLSZ) Matrix to multiply OP (operator).
%       N  - Size of the matrix A.
%    BLSZ  - Block size.
%   MPROD  - Number of matrix-vector products.
%
% Output:
%       X  - (N x BLSZ) Product of OP*X (operator).
%   MPROD  - Number of matrix-vector products.

% James Baglama
% DATE: 7/20/2005

% Generalized eigenvalue problem (Cholesky factored).
if ~isempty(B) & isempty(L) & isempty(nval) 
   X(permB,:) = B \ X;
   if ischar(A), X = feval(A,X,n,blsz); else X = A*X; end
   X = (X(permB,:)'/B)';
end

% Generalized eigenvalue problem (LU factored).
if ~isempty(B) & ~isempty(L) & isempty(nval)
   if ischar(A), X = feval(A,X,n,blsz); else X = A*X; end 
   X = U \ (L \ (P * X));
end   

% Generalized or standard eigenvalue problem (Shift & Invert).
if ~isempty(nval)
   if isempty(B)
      X = U \ (L \ (P * X));
   else
      X = U \ (L \ (P * B* X));
  end
end

% Standard eigenvalue problem.
if isempty(B) & isempty(nval)  
   if ischar(A), X = feval(A,X,n,blsz); else X = A*X; end  
end

% Number of matrix-vector products with A.
mprod = mprod+blsz; 

%-----------------------------%
% END: MATRIX-VECTOR PRODUCT. %
%-----------------------------%

