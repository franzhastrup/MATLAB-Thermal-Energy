%% *****************************************************************************
% File: FiniteVolumeConvectionDiffusion2D.m
%   Script for simulation of steady convection-diffusion heat problem on two-
%   dimensional square domain using the Finite Volume Method on a cartesian,
%   structured mesh with square cells. Central difference fluxes applied for
%   the diffusive terms, and either central of upwinded difference fluxes
%   applied for the convective terms.
%   Implementation does not include source terms and is limited to the
%   following two sets of convective velocity fields and boundary conditions:
%     1) uniform flow [u,v] = [1,0], homogeneous Neumann BC. at north and south
%     walls (dTn/dy=dTs/dy=0), homogeneous Dirichlet BC at west wall (Tw=0),
%     inhomogeneous Dirichlet BC at east wall (Te=1).
%     2) corner flow [u,v] = [-x,y], homogeneous Neumann BC. at north wall
%     (dTn/dy=0), homogeneous Dirichlet BC at west wall (Tw=0), inhomogeneous
%     Dirichlet BC at east wall (Te=1), inhomogeneous Dirichlet BC
%     at south wall (Ts=x)
%   Linear system of equations solved directly using Matlab's backslash
%   operator
%
% Input         :
%   n           :  Number of cells along x,y-axis
%   L           :  Size of square in x,y-direction
%   Pe          :  Global Peclet number
%   problem     :  Problem #: 1 or 2, selects case of convective field and BCs
%   linsolver   :  Linear solver: 'direct','jacobi','gs','sor'
%   fvscheme    :  Finite volume scheme for convection-diffusion, 
%                  either 'cds-cds' or 'uds-cds' 
%   omegaSOR    :  Relaxation parameter omega for SOR iterative method
%   imax        :  Maximum number of iterations for iterative linear solver
%   tol         :  Relative tolerance on for iterative linear solver,
%                  should be selected such that tol > condest(A)*eps
%
% Output        :
%   T           :   Temperature at cell nodes, T(1:n,1:n)
%   A           :   Convection-diffusion system matrix, A(1:n^2,1:n^2)
%   s           :   Source array with BC contributions, s(1:n,1:n)
%   TT          :   Temperature field extrapolated to walls, TT(1:n+2,1:n+2)
%   CF,DF       :   Conv. and diff. fluxes through walls, CF=[CFw,CFe,CFs,CFn]
%   GHC         :   Global heat conservation, scalar (computed from wall fluxes)
%   solve_time  :   CPU-time for solving system of linear equations A*T(:)=s(:)
%   error       :   Vector of relative errors, in case of iterative lin. solver
%   iter        :   Number of iterations made, in case of iterative lin. solver
%   flag        :   Convergence flag (=0: converged) for iterative lin. solver
%   Thist       :   Array of iterates, in case of iterative lin. solver
%   Plots of the temperature field and convective velocity field
%
%   Author      : Franz Hastrup-Nielsen
%   Date        : 6/10/2015
%   Version     : 1.3
%*******************************************************************************

%% Clear
clear all
close all
clc

%% =================================== INPUT ===================================
n         = 1000;               % number of cells along x,y-axis
L         = 1;                  % size of square in x,y-direction
Pe        = 7;                  % global Peclet number
problem   = 1;                  % problem to solve (1 or 2 - see file header)
linsolver = 'direct';           % linear solver 
fvscheme  = 'uds-cds';          % finite volume scheme ('cds-cds' or 'uds-cds')
omegaSOR  = 1.93;               % relaxation parameter for SOR (0<=omegaSOR<2)
imax      = 10000;              % maximum number of iterations for linsolver
tol       = 1e-12;              % relative tolerance for linsolver
% ================================ INPUT DONE =================================



%% ============================== PRE-PROCESSING ===============================
% << assemble matrix A and source array s based on the arguments >>
% << (n,L,Pe,problem,fvscheme), following the 6-step approach below. >>
% << After implementation as script move pre-processing code to the function: >>
%
% <<    [A,s] = FVconvdiff2D_preprocess(n,L,Pe,problem,fvscheme) >>

%% 1) Coordinate arrays
% Input  : n, L
% Output : dx cell size
%          xf cell face coordinate vector (of length n+1)
%          xc cell center coordinate vector (of length n)
dx = L/n;                       % cell size in x,y-direction
xf = 0:dx:L;                    % cell face coordinate vector along x,y-axis
xc = dx/2:dx:L-dx/2;            % cell center coordinates along x-axis

%% 2) Generate convective velocity field at cell faces
% Input  : n, xf
% Output : UFw, UFe, VFs, VFn cell faces velocities (arrays of size n*n)
if problem == 1                     % problem 1 - uniform flow
    UFw = ones(n);                  % western face velocity, u = 1
    UFe = ones(n);                  % eastern face velocity, u = 1    
    VFs = zeros(n);                 % southern face velocity, v = 0
    VFn = zeros(n);                 % northern face velocity, v = 0    
elseif problem == 2                 % problem 2 - corner flow
    UFw = repmat(-xf(1:n),n,1);     % western face velocity, u = -x
    UFe = repmat(-xf(2:n+1),n,1);   % eastern face velocity, u = -x   
    VFs = repmat(xf(1:n)',1,n);     % southern face velocity, v = y
    VFn = repmat(xf(2:n+1)',1,n);   % northern face velocity, v = y    
else                                % problem unknown
    disp('problem not implemented')
end

%% 3) Generate convective face flux matrices
% << Input  : n, dx, Pe, UFw, UFe, VFs, VFn >>
% << Output : Fw, Fe, Fs, Fn convetive fluxes (arrays of size n*n), e.g. >>
% <<          is Fe the convective fluxes at the east-faces of the cells. >>
Fw = Pe*UFw*dx;
Fe = Pe*UFe*dx;
Fn = Pe*VFn*dx;
Fs = Pe*VFs*dx; %Page 121, table 5.3
D = -1; %Page 121, table 5.3

%% 4) Generate coefficient arrays
% << Input  : n, dx, Pe, fvscheme, Fw, Fe, Fs, Fn >>
% << Output : aW, aS, aE, aN, aP coefficient arrays (arrays of size n*n), e.g. >>
% <<          is aE the coefficients for the cells east of the one considered. >>
if strcmp(fvscheme,'cds-cds')     % CDS-CDS FV-scheme applied
    if dx*Pe>=2
       fprintf('Warning cell Peclet number to small') 
    end
    aW = D-Fw/2;
    aS = D-Fs/2;
    aE = D+Fe/2;
    aN = D+Fn/2; % Page 120, table 5.2
    % << compute coefficients for cds-cds scheme here >>
    % << check whether the cell Peclet number dx*Pe>=2 - then print warning >>
elseif strcmp(fvscheme,'uds-cds') % UDS-CDS FV-scheme applied
    aW = D+min(0,-Fw);
    aS = D+min(0,-Fs);
    aE = D+min(0,Fe);  
    aN = D+min(0,Fn);% Page 120, table 5.2
    % << compute coefficients for uds-cds scheme here >>
else                              % fvscheme unknown
    disp('fvscheme not implemented')
end
aP = -(aW+aE+aN+aS)+Fe-Fw+Fn-Fs; %Eq. 5.59

%% 5) Impose boundary conditions and compute source array from BC
% << Input  : n, dx, problem, aW, aS, aE, aN, aP >>
% << Output : s source term array (size n*n) >>
% <<          aW, aS, aE, aN, aP modified coefficient arrays (size n*n) >>
s = zeros(n);   % initialize source array

if problem == 1 % BC: inhomogeneous Dirichlet BC Tw=0, homogeneous Dirichlet
    %%BC Te=1, homogeneous Neumann BC. at north and south dTn/dy=dTs/dy=0
    Te = 1;
    Tw = 0;
    Dn = 0;
    Ds = 0;
    %Pg. 108-110
    aP(:,1) = aP(:,1)-aW(:,1); %Dirichlet in west
    s(:,1) = s(:,1)-2*aW(:,1)*Tw; %Source term in west
    aW(:,1) = 0;
    
    aP(:,end) = aP(:,end)-aE(:,end); %Dirichlet in east
    s(:,end) = s(:,end)-2*aE(:,end)*Te; %Source term in west
    aE(:,end) = 0;
    
    aP(end,:) = aP(end,:)+ aN(end,:);%Neumann north
    s(end,:) = s(end,:)+aN(end,:)*dx*Dn;
    aN(end,:) = 0;
    
    aP(1,:) = aP(1,:)+ aS(1,:);%Neumann south
    s(1,:) = s(1,:)-aS(1,:)*dx*Ds;
    aS(1,:) = 0;
    
elseif problem == 2 % homogeneous Dirichlet BC Tw=0, inhomogeneous
    %     Dirichlet BC Te=1 homogeneous Neumann BC (dTn/dy=0) inhomogeneous
    %     Dirichlet BC Ts=x
    Tw = 0;
    Te = 1;
    Dn = 0;
    Ts = xc;
    
    aP(:,1) = aP(:,1)-aW(:,1); %Dirichlet in west
    s(:,1) = s(:,1)-2*aW(:,1)*Tw; %Source term in west
    aW(:,1) = 0;
    
    aP(:,end) = aP(:,end)-aE(:,end); %Dirichlet in east
    s(:,end) = s(:,end)-2*aE(:,end)*Te; %Source term in west
    aE(:,end) = 0;
    
    aP(end,:) = aP(end,:)+ aN(end,:);%Neumann north
    s(end,:) = s(end,:)+aN(end,:)*dx*Dn;
    aN(end,:) = 0;
    
    aP(1,:) = aP(1,:)- aS(1,:); %Dirichlet south
    s(1,:) = s(1,:)-2*aS(1,:).*Ts(1,:);
    aS(1,:) = 0;
end



%% 6) Assemble system matrix
% << Input  : n, aW, aS, aE, aN, aP  >>
% << Output : A 5-diagonal system matrix (array of size n^2*n^2) >>
aW = aW(:);
aN = aN(:);
aS = aS(:);
aE = aE(:);
aP = aP(:);
s = s(:);
DiaPos = [-n -1 0 1 n];
Adiags = [[aW(n+1:end);zeros(n,1)],[aS(2:end);0],[aP],[0;aN(1:end-1)],[zeros(n,1);aE(1:end-n)]];     
A      = spdiags(Adiags,DiaPos,n^2,n^2);  

% << hints  : >>
% <<  - use function 'spdiags' to assemble A. The coefficient arrays >>
% <<    aW, aS, aE, aN, aP must in advance have been turned into column >>
% <<    vectors by setting aW=aW(:); etc. >>
% <<  - since the spdiags function places the first element of each >>
% <<    coefficient vector in the first column of the matrix, the coefficient >>
% <<    vectors aW,aS,aE,aN will need to be appropriately shifted and padded, >>
% <<    that is removing existing elements and adding zeros to either the start >>
% <<    or the end of the vectors corresponding to the diagonals in A, where >>
% <<    the vectors will be located.
% <<  - e.g. 1D-laplacian of report 1, ex. 2.4.1 may be assembled as follows >>
% <<     n      = 5;                                         >>
% <<     ovec   = ones(n-1,1);                               >>
% <<     Adiags = [[ovec;0;0],[1;-2*ovec;1],[0;0;ovec]];     >>
% <<     A      = spdiags(Adiags,-1:1,n+1,n+1);              >>
% ============================= PRE-PROCESSING DONE ============================



%% ===================== SOLVE SYSTEM OF LINEAR EQUATIONS ======================
% Solve equations using either direct or iterative method
% Input  : n, linsolver, T, s, >>
% Output : T temperatur field array (size n*n) >>
%          solve_time wall-clock time to solve linear system 
%          error, iter, flag, Thist (in case of iter. lin. solver) >>

if exist('A','var') && exist('s','var')
    T = zeros(n);                   % initialize temperature field
    solve_time = tic;               % start timing using Matlab's tic-toc function
    if strcmp(linsolver,'direct')   % direct solver using Matlab's backslash [\]
        T(:) = A\s(:);              % temperature field
    elseif any(strcmp(linsolver,{'jacobi','gs','sor'}))% stationary iterative solver
        [T(:),error,iter,flag,Thist] = ...% temperature field, error-vector, iterations,
            StationaryIterativeLinearSolver(A,s(:),linsolver,imax,tol,omegaSOR);
    else                            % linear solver unknown
        disp('linear solver not implemented')
    end
    solve_time = toc(solve_time);   % stop timing using Matlab's tic-toc function
end
% ==================== SYSTEM OF LINEAR EQUATIONS SOLVED ======================


%% ============================== Post-processing ==============================
% << compute extended temperature field TT, net convective and diffusive  >>
% << fluxes through the walls CF,DF and global heat conservation GHC.  >>
% << after implementation as script move entire post-processing to function: >>
%
% <<   [TT,GHC,CF,DF] = FVconvdiff2D_postprocess(T,n,L,Pe,problem,fvscheme) >>
%
% << Note: some operations repeated from pre-processing to retain portability >>

%% 1) Coordinate arrays
% Input  : n, L
% Output : dx cell size
%          xf cell face coordinate vector (of length n+1)
%          xc cell center coordinate vector (of length n)
dx = L/n;                       % cell size in x,y-direction
xf = 0:dx:L;                    % cell face coordinate vector along x,y-axis
xc = dx/2:dx:L-dx/2;            % cell center coordinates along x-axis

%% 2) Generate convective velocity field
% Input  : n, xf
% Output : UFw, UFe, VFs, VFn cell faces velocities (arrays of size n*n)
if problem == 1                 % problem 1 - uniform flow
    UFw = ones(n);              % western face velocity, u = 1
    UFe = ones(n);              % eastern face velocity, u = 1    
    VFs = zeros(n);             % southern face velocity, v = 0
    VFn = zeros(n);             % northern face velocity, v = 0    
elseif problem == 2             % problem 2 - corner flow
    UFw = repmat(-xf(1:n),n,1); % western face velocity, u = -x
    UFe = repmat(-xf(2:n+1),n,1); % eastern face velocity, u = -x   
    VFs = repmat(xf(1:n)',1,n); % southern face velocity, v = y
    VFn = repmat(xf(2:n+1)',1,n); % northern face velocity, v = y    
else                            % problem unknown
    disp('problem not implemented')
end

%% 3) Extrapolate temperature field to walls
% << Input  : n, dx, fvscheme, problem, T >>
% << Output : TT extended temperature field array (of size (n+2)*(n+2) ) >>
% 
% << hints  : >>
% <<  - the temperature in the corners is not uniquely defined from the >>
% <<    extrapolation expressions, however, as the corner temperatures >>
% <<    are not used, just leave the corners equal to zero (plotting the >>
% <<    temperatur field with 'surf' will not show this) >>
% <<  - see notes section 5.3 >>
TT = zeros(n+2);              % initialize extended temperature array


if strcmp(fvscheme,'cds-cds')     % CDS-CDS FV-scheme applied
    if problem==1
        TT(2:end-1,1) = Tw;
        TT(2:end-1,end) = Te;
        TT(1,2:end-1) = T(1,:)-1/2*dx*Ds; %???
        TT(end,2:end-1) = T(end,:)+1/2*dx*Dn; %???
        
    elseif problem==2
        TT(2:end-1,1) = Tw;
        TT(2:end-1,end) = Te;
        TT(1,2:end-1) = Ts; %???
        TT(end,2:end-1) = T(end,:)-1/2*dx*Dn; %???
    end
elseif strcmp(fvscheme,'uds-cds')
    if problem==1
        TT(2:end-1,1) = 2*Tw-T(:,1); %west
        TT(2:end-1,end) = T(1,end); %east
        TT(1,2:end-1) = T(1,:)-dx*Ds; 
        TT(end,2:end-1) = T(end,:); 
        
    elseif problem==2
        TT(2:end-1,1) = T(:,1); %west
        TT(2:end-1,end) = 2*Te-T(:,end); %east
        TT(1,2:end-1) = 2*Ts-T(1,:); 
        TT(end,2:end-1) = T(end,:); 
    end
end
TT(2:end-1,2:end-1) = T;

%% 4) Compute temperature gradients at walls
% << Input  : n, dx, problem, T
% << Output : DTw, DTe, DTs, DTn temperature wall gradients (vectors of length n) >>
if problem==1
    DTw = 2/dx*(T(:,1)-Tw); 
    DTe = 2/dx*(Te-T(:,end)); 
    DTs = Ds;
    DTn = Dn;
elseif problem==2
    DTw = 2/dx*(T(:,1)-Tw); 
    DTe = 2/dx*(Te-T(:,end)); 
    DTs = 2/dx*(T(1,:)-Ts);
    DTn = Dn;
end
% << hints: >>
% <<  - see notes section 5.3 >>

%% 5) Convective and diffusive wall fluxes, and global conservation of heat
% << Input  : n, dx, TT, DTw, DTe, DTs, DTn, Pe, UFw, UFe, VFs, VFn >>
% << Output : CF=[CFw,CFe,CFs,CFn] convective fluxes summed up over walls >>
% <<          DF=[DFw,DFe,DFs,DFn] diffusive fluxes summed up over walls >>
% <<          GHC global heat conservation (scalar) >>

CFw = sum(Pe*dx*UFw(:,1).*TT(2:end-1,1));
CFe = sum(Pe*dx*UFe(:,end).*TT(2:end-1,end));
CFs = sum(Pe*dx*VFs(1,:).*TT(1,2:end-1));
CFn = sum(Pe*dx*VFn(end,:).*TT(end,2:end-1));
DFw = dx*sum(DTw);
DFe = dx*sum(DTe);
DFn = dx*sum(DTn);
DFs = dx*sum(DTs);
CF=[CFw,CFe,CFs,CFn];
DF=[DFw,DFe,DFs,DFn];

GHC = CF(1)+CF(3)-DF(1)-DF(3)-(CF(2)+CF(4)-DF(2)-DF(4))

%
% << hints: >>
% <<  - see notes section 5.4 >>
% ============================ POST-PROCESSING DONE ============================



%% ================================= PLOT RESULTS ===============================
% << plot various results, e.g. the obtained temperature field and the >>
% << convective velocity field >>

%% Coordinate arrays
% Input  : n, L
% Output : dx cell size
%          xc cell center coordinate vector (of length n)
%          xt cell center coor. vector extended with walls (of length n+2)
%          Xc,Yc cell center arrays (match T, size n*n)
%          Xt,Yt cell center arrays extended with walls (match TT, size (n+2)*(n+2))
dx = L/n;                       % cell size in x,y-direction
xc = dx/2:dx:L-dx/2;            % cell center coordinate vector along x,y-axis
xt = [0 xc L];                  % extended cell center coor. vector, incl. walls
[Xc,Yc] = meshgrid(xc);         % cell center coordinate arrays
[Xt,Yt] = meshgrid(xt);         % extended cell center coor. arrays, incl. walls

%% Plots
figure 
surf(Xc,Yc,T,'EdgeColor','none'), view(0,90), colorbar
caxis([0 1])
xlim([0 1])
ylim([0 1])
grid on
axis equal
title(sprintf('Temperature field, Problem %i, FV %s, n = %i, Pe = %-4.0f',problem,fvscheme,n,Pe),'Interpreter','Latex','Fontsize',14)
ylabel('y','Interpreter','Latex','Fontsize',14)
xlabel('x','Interpreter','Latex','Fontsize',14)

print('-f1','-depsc2',sprintf('Problem%i_FV%s_n%i_Pe%-4.0f.eps',problem,fvscheme,n,Pe))