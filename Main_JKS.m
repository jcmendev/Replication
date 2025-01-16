
clc
close all; 
clear all;
% Code to replicate Johri, Khan, Sosa-Padilla

fprintf('\n  Johri, Khan, Sosa-Padilla \n');
tic 


set(0,'DefaultFigureWindowStyle','docked');

%% Algorithm parameters
maxit   = 10000               ; % maximum number of iterations
tol     = 1e-8              ; % tolerance for q-convergence
prntitq = 10                   ; % print q-iterations
prntitv = 0                   ; % print v-iterations

%% Model parameters 
% % Standard parameters
mu      = 2         ; % Risk aversion                                 
rf      = 0.01      ; % Riskfree (annual)
qf      = 1/(1+rf)  ; % Risk free bond price
Rf      = (1+rf)    ; % Gross risk free rate
beta    = 0.96      ; % Household Discount Factor
% % Parameters from literature
lambda  = .0385     ; % Probability of re-entry Chatterjee 
delta   = 0.0742    ; % Average debt duration
% delta   = 1    ; % Average debt duration
coupon = (rf + delta)*qf; % Macaulay Decaying Rate 
%% Parameters set to match model moments  
kappa   = 8.75 ;
% defc0   = .0;
defc0   = -0.47;
defc1   = 0.60 ;
%% Model implied parameter values
% A     = 1+r*Eb-Ec                            ; % rescaling parameter
A     = 0;
%% Endogenous VFariable State-space                   
% %  
bmax = 0;
bmin = -2.5;
nb   = 101;
b    = linspace(bmin,bmax,nb)';
%% Exogenous Processes
% % Output Tauchen
ny      = 15        ; % # of discrete points in y   
rho_y   = 0.933  ; 
var_y   = 0.027^2;
mean_y  = -0.5*var_y;

[logy,Py] = MC_Tauchen(ny,mean_y,rho_y,sqrt(var_y),4);
y         = exp(logy);

% % Interest rate volatility
nr_vol        = 7 ;
std_r_vol     = 0.2632;
rho_r_vol     = 0.8742;
mean_r_vol    = -6.2869;
[r_vol,Pr_vol,dr_vol] = MC_Rouwenhorst(nr_vol,std_r_vol,mean_r_vol,rho_r_vol);

if nr_vol == 1 
    Pr_vol = 1;
end

% % Interest rate Level
nr_lev        = 7;
rho_r_lev     = 0.908;
mean_r_lev    = 0.0;

for irvol = 1:nr_vol
    [r_lev(:,irvol),Pr_lev,dr_lev] = MC_Rouwenhorst(nr_lev,exp(r_vol(irvol)),mean_r_lev,rho_r_lev);
end

if nr_lev == 1 
    Pr_lev = 1;
end
r_lev = rf+r_lev;
%% State space and actions 
% Build the space of ACTIONS: SxB. That is, the next periods value of the
% endogenous state variables. "sp" and "bp" denote "s prime" and "b prime"
bp =  b       ;    m = length(bp) ;  

% Build the space of exogenous states
aux        = [vec(r_lev) kron(r_vol,ones(nr_lev,1))];
[yy,rrvol] = gridmake(y,aux(:,2));
rrlev      = vec(repmat(aux(:,1),1,ny)');

nz = length(rrlev);
% Build the space of states
  
[B,Y,RVOL] = gridmake(b,y,aux(:,2))     ;    n = length(B);
RLEV       = vec(repmat(aux(:,1),1,ny*nb)');

S = [B,Y,RLEV,RVOL];
% Build Matrix
Pz = kron(Py,kron(Pr_lev,Pr_vol));  

%% Pre-allocate for speed
cpay   = zeros(n,m)          ;       cdef   = zeros(n,m);
i0  = getindex(0,bp);   % action positions where b'=0

%% Punishment setup 

phi_Y = max(0,defc0*Y + defc1*Y.^2);
Ydef  = Y-phi_Y  ;
%%  Consumption under default state
for j=1:m                 
    cdef(:,j)   = Ydef - A;
end

% return

%% Initialize
q   = qf*ones(n,m)  ;         % guessed bond price function (risk-free)

vs = zeros(n,1);
vc = zeros(n,1);
V  = zeros(n,1); 
rho_q_new = 0.1;
rho_V_new = 0.2;

rho_q_new = 0.7;
rho_V_new = 0.7;  
        
try
load Policies_JKSP
catch
end
% q   = qf*zeros(n,m)  ;  
qrf = q;
%% Iterate on q and v on one loop
  
fprintf('\nSolving Social Planner Problem     \n') 
fprintf('\nIter     q-change    V-change\n');
for it = 1:maxit
  Vold = V; 
  v0   = getv0(i0,V,n,m);  
  qold = q;    
  qrfold = qrf;
  
 % 1. Construct the reward function (with and w/o default), for a given q
  
   for j=1:m
         %Consumption under repayment state      
         cpay(:,j)   = Y + B*coupon - q(:,j).*(bp(j)-(1-delta)*B) - A;
   end
    
  upay = (cpay.^(1-mu))./(1-mu); upay(cpay<=0)=-Inf; 
  udef = (cdef.^(1-mu))./(1-mu); udef(cdef<=0)=-Inf;
   
 % 2. Solve the sovereign problem, for a given q

%% 
    EV  = kron(Pz*(reshape(V,m,n/m)'),ones(m,1));
    Evc = kron(Pz*(reshape(vc,m,n/m)'),ones(m,1));
    Evs = kron(Pz*(reshape(vs,m,n/m)'),ones(m,1));
    EV0 = kron(Pz*(reshape(v0,m,n/m)'),ones(m,1));
     
    [vc,xc] = max(upay+beta.*EV,[],2);
    [vs,xs] = max(udef+lambda*(beta.*EV0)+(1-lambda)*(beta.*Evs),[],2);

%%
    [V,xv]  = max([vc vs],[],2); 
    d       = zeros(n,1);
    d(xv==2)= 1; 

%   % 3. Compute q 
% 
    qp_aux = zeros(n,1);
    qprf_aux = zeros(n,1);    
    for i=1:n
%         qp_aux(i,1)=(1-d(i))*q(i,xc(i)); 
        qp_aux(i,1)=q(i,xc(i)); 
        qprf_aux(i,1)= qrf(i,xc(i)); 
    end
  
%     kappa = 0;
 qp    = reshape(qp_aux,m,nz)';
 qprf    = reshape(qprf_aux,m,nz)';
 defp  = (reshape(d,m,nz)');
 epsyp = log(yy)-rho_y*log(yy')-(1-rho_y)*mean_y; 
%  argp  = ((1-qf)/qf) + kappa*(epsyp+0.5*kappa*sqrt(var_y));
 argp  = rrlev + kappa*(epsyp+0.5*kappa*sqrt(var_y));
 msp   = exp(-argp);
 
 
 ERHS   = zeros(nz,m);
 ERHSrf = zeros(nz,m);
 for iz = 1:nz 
    aux = msp(:,iz).*(1-defp)*coupon + msp(:,iz).*(1-defp).*(1-delta).*qp;
    ERHS(iz,:) = Pz(iz,:)*aux;
    
    auxrf = msp(:,iz).*coupon + msp(:,iz).*(1-delta).*qprf;
    ERHSrf(iz,:) = Pz(iz,:)*auxrf;
    
 end
%  ERHS = Pz*((1-defp).*(1+(1-delta)*qp));
%  ERHS = kron(ERHS,ones(m,1));     % add b(t) to rows to make same size as q
%  q = qf*ERHS;
  % 4. Update the price "q" (i.e. find a new "q")

   q = kron(ERHS,ones(m,1));     % add b(t) to rows to make same size as q
   qrf = kron(ERHSrf,ones(m,1));     % add b(t) to rows to make same size as q
%    q = min(max(q,0),(1+(1-delta)*qf)*qf);

  % 5. Stop if the old "q" is similar to new "q"
   
   dq = norm(q-qold);
   dqrf = norm(qrf-qrfold);
   dV = norm(V-Vold);
   
   error_it = [dq dV dqrf];
   if mod(it, prntitq) == 0; fprintf(' %4.0f   %8.6f    %8.6f\n',it,dq,dV); end
   if max(error_it)<tol,break,end

   % 6. Updating smoothing
   q = rho_q_new*q+(1-rho_q_new)*qold;
   V = rho_V_new*V+(1-rho_V_new)*Vold; 

end

fprintf('\nSocial Planner Problem solved in %2.2f minutes. \n',toc/60);

%% Optimal policy & transition prob matrix, for the equilibrium "q"

bprp = bp(xc);
bppol = (1-d).*bprp;

for i=1:n
    qstar(i,1)=(1-d(i))*q(i,xc(i)); 
    c(i,1) = (1-d(i))*cpay(i,xc(i))+d(i)*cdef(i,xs(i));

    
    qrp(i,1) = q(i,xc(i));
    crp(i,1) = cpay(i,xc(i));  

    cdefault(i,1) = cdef(i,xs(i));
    
end
% 
% Edef2 = Pz*(reshape(d,m,ny)');
% qstar2 = vec((1-Edef2')/(1-rf));
%%  Simulation
simulation_JKS
%% 
rq = reshape(qstar,nb,nz);
figure()
hold on
for iz = 1:nz
plot(b,rq(:,iz))
end
rd = reshape(d,nb,ny,nr_lev,nr_vol);

%% Default Sets
close all 

figure('Color','w')
ix = 1;
for il = [1 nr_lev]
    for iv = [1 nr_vol]
        subplot(2,2,ix)
        mesh(y,b,rd(:,:,il,iv))
        xlabel('Income','Interpreter','Latex')
        ylabel('Debt','Interpreter','Latex')
        title(['DSets r=',num2str(r_lev(il),2),' $\sigma_{r}=$',num2str(r_vol(iv),2),', pct def states=',num2str(100*sum(sum(rd(:,:,il,iv)))/(nb*ny),4)],'Interpreter','Latex')
        xlim([y(1) y(end)])
        ylim([b(1) b(end)])
        view(0,90)
        ix = ix+1;
    end
end
% %% Issuance and debt t+1 Policies 

ixy = ceil(ny);
ixy = 15;
ixlev = [2 6];
ixvol = [6 7];
bd = bppol - (1-delta).*B;
bd_y = bd./Y;

rbdy = reshape(-bd_y,nb,ny,nr_lev,nr_vol);
rbppol = reshape(bppol,nb,ny,nr_lev,nr_vol);
rvc    = reshape(vc,nb,ny,nr_lev,nr_vol);
rvs    = reshape(vs,nb,ny,nr_lev,nr_vol);

figure('Color','w')
subplot(1,2,1)
hold on
plot(-b,rbdy(:,ixy,ixlev(1),ixvol(1)),'r')
plot(-b,rbdy(:,ixy,ixlev(2),ixvol(1)),'b')
% plot(b,b,'--k')
% xlim([0.1,0.7])
ylim([0.0,0.6])
ylabel('Issuance','Interpreter','Latex')
xlabel('Debt at t ','Interpreter','Latex')
title('Low Volatility','Interpreter','Latex')
legend({'Low r','High r'},'Interpreter','Latex')
hold off
subplot(1,2,2)
hold on
plot(-b,rbdy(:,ixy,ixlev(1),ixvol(2)),'r')
plot(-b,rbdy(:,ixy,ixlev(2),ixvol(2)),'b')
% plot(b,b,'--k')
hold off
% xlim([0.1,0.7])
ylim([0.0,0.6])
ylabel('Issuance','Interpreter','Latex')
xlabel('Debt at t ','Interpreter','Latex')
title('High Volatility','Interpreter','Latex')
legend({'Low r','High r'},'Interpreter','Latex')
print(['Fig_Issuance'],'-depsc','-r0')

figure('Color','w')
subplot(1,2,1)
hold on
plot(-b,-rbppol(:,ixy,1,1),'b')
plot(-b,-rbppol(:,ixy,end,1),'r')
plot(-b,-b,'--k')
hold off
ylabel('Debt at t+1','Interpreter','Latex')
xlabel('Debt at t ','Interpreter','Latex')
title('Low Volatility','Interpreter','Latex')
legend({'Low r','High r'},'Interpreter','Latex')
subplot(1,2,2)
hold on
plot(-b,-rbppol(:,ixy,1,end),'b')
plot(-b,-rbppol(:,ixy,end,end),'r')
plot(-b,-b,'--k')
hold off
ylabel('Debt at t+1','Interpreter','Latex')
xlabel('Debt at t ','Interpreter','Latex')
title('High Volatility','Interpreter','Latex')
legend({'Low r','High r'},'Interpreter','Latex')
print(['Fig_DebtPol'],'-depsc','-r0')


figure('Color','w')
subplot(1,2,1)
hold on
plot(-b,rvc(:,ixy,ixlev(1),ixvol(1)),'r')
plot(-b,rvc(:,ixy,ixlev(2),ixvol(1)),'b')
plot(-b,rvs(:,ixy,ixlev(1),ixvol(1)),'--r')
plot(-b,rvs(:,ixy,ixlev(2),ixvol(1)),'--b')
% plot(b,b,'--k')
% xlim([0.1,0.7])
% ylim([0.0,0.5])
ylabel('V','Interpreter','Latex')
xlabel('Debt at t ','Interpreter','Latex')
title('Low Volatility','Interpreter','Latex')
legend({'Low r','High r'},'Interpreter','Latex')
hold off
subplot(1,2,2)
hold on
plot(-b,rvc(:,ixy,ixlev(1),ixvol(2)),'r')
plot(-b,rvc(:,ixy,ixlev(2),ixvol(2)),'b')
plot(-b,rvc(:,ixy,ixlev(1),ixvol(2)),'--r')
plot(-b,rvs(:,ixy,ixlev(2),ixvol(2)),'--b')
% plot(b,b,'--k')
hold off
% xlim([0.1,0.7])
% ylim([0.0,0.5])
ylabel('V','Interpreter','Latex')
xlabel('Debt at t ','Interpreter','Latex')
title('High Volatility','Interpreter','Latex')
legend({'Low r','High r'},'Interpreter','Latex')


