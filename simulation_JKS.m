% SIMUL_MODEL
%
% Simulates JKSP
%% Simulation parameters and switches
  L               = 1           ;   % Lags for autocorrelation    
  T               = 1000       ;   % Number of draws simulated 5000
  T_drop          = 500 ;   % Number of draws to drop
% Report the type of model data treatment 
      fprintf('\nModel simulated no-detrended \n')   

%% Simulation of exogenous states 
% Grid for P and Y
  ygrid = y;
% Initial exogeneous state
  iz0 = getindex(exp(mean_y),ygrid);
  
% Markov Chain simulation
  iz = simulmarkov(Pz,T,iz0);
  y_sim = yy(iz); 
  rlev_sim = rrlev(iz);
  rvol_sim = rrvol(iz);

% Redemption simulation      
  redemp = ones(floor(lambda*T),1);
  redemp = [redemp; zeros(T-length(redemp),1)];
  redemp = redemp(randperm(T)');

%% Simulation of endogenous states   

% Empty Tx1 matrices to store simulated variables 
  b_sim   = zeros(T,1); 
  bp_sim  = zeros(T,1);
  q_sim   = zeros(T,1); 
  r_sim   = zeros(T,1); 
  d_sim   = zeros(T,1);
  c_sim   = zeros(T,1);
  history = zeros(T,1);
  state   = zeros(T,1);

% Initial endogenous state   
%   s_sim(1) = mean(s);
  b_sim(1) = mean(b);
%   b_sim(1) = b(getindex(Eb,b));
       
% Model simulation  
 for t = 1:T   
     
  state(t) = getindex([y_sim(t) rlev_sim(t) rvol_sim(t) b_sim(t)],[Y RLEV RVOL B]);
  
   if history(t)==0        
     
     d_sim(t)     = d(state(t));    
       
      if d_sim(t) == 1
      y_sim(t)     = Ydef(state(t));      
      q_sim(t)     = 0;
      r_sim(t)     = NaN;
      bp_sim(t)    = 0;
      c_sim(t)     = y_sim(t) ;    
      history(t+1) = 1;

      elseif d_sim(t) == 0 
            
      q_sim(t)   = qstar(state(t)); 
      r_sim(t)   = (1/q_sim(t))-1;
      bp_sim(t)  = bppol(state(t));
      c_sim(t)     = y_sim(t) + b_sim(t)*coupon - q_sim(t)*(bp_sim(t)-(1-delta)*b_sim(t)) -A;  
      
      history(t+1) = 0;

      end
        
   elseif  history(t) == 1 & redemp(t) == 0
      
      y_sim(t)     = Ydef(state(t));     
      q_sim(t)     = 0;
      r_sim(t)     = NaN;
      bp_sim(t)    = 0;
      c_sim(t)     = y_sim(t) ;
      
      history(t+1) = 1;
      
   elseif history(t)==1 & redemp(t)==1
      y_sim(t)     = Ydef(state(t));
      d_sim(t)   = 0;
      q_sim(t)   = qstar(state(t));
      r_sim(t)   = (1/q_sim(t))-1;
%       bp_sim(t)  = 0;
      bp_sim(t)  = bppol(state(t));      
     c_sim(t)     = y_sim(t) + b_sim(t)*coupon - q_sim(t)*(bp_sim(t)-(1-delta)*b_sim(t)) -A;  
      history(t+1) = 0;
      
       
   end
      % Update state values   
      b_sim(t+1) = bp_sim(t);
        
 end

 % Drop first tranch of simulation
 
 b_sim   = b_sim(T_drop+1:end-1); 
 history = history(T_drop+1:end-1);
 bp_sim  = bp_sim(T_drop+1:end);
 q_sim   = q_sim(T_drop+1:end); 
 r_sim   = r_sim(T_drop+1:end); 
 d_sim   = d_sim(T_drop+1:end);
 c_sim   = c_sim(T_drop+1:end);
 y_sim   = y_sim(T_drop+1:end);
 redemp  = redemp(T_drop+1:end);

 % Computation of other simulated variables 

 GDP_sim = y_sim;
 by_sim  = b_sim./GDP_sim;
 TB_sim  = GDP_sim-c_sim-A; % esta le quite los costos de extraccion 
 NX_sim  = GDP_sim-c_sim; 
 TBY_sim = TB_sim./GDP_sim;
 NXY_sim = NX_sim./GDP_sim; 
 CA_sim  = q_sim.*bp_sim-b_sim;
 CAY_sim = CA_sim./GDP_sim;

% Count default and repayment states
 defrate = nanmean(d_sim);
 defrate = mean(d_sim(find(b_sim<0)));
 repay   = find((history==0)&(d_sim==0));
 default = find((history==1)&(redemp==0));

 spread = ((1/(1+rf))*(1+r_sim)-1) ;


  [ct,cc] = hpfilter(c_sim,1600);
 [yt,yc] = hpfilter(y_sim,1600);
 [st,sc] = hpfilter(spread,1600);
 %% 

fprintf('\n')
fprintf('\n                   Data to Match                  \n') 
fprintf('\n   Variable                        Paper       Model ')
% fprintf('\nDefault rate (pct)                 %5.2f      %5.2f' , NaN, defrate*100)
fprintf('\nDebt to GDP (pct)                  %5.2f      %5.2f' , 44, -mean(b_sim./y_sim)*100)
fprintf('\n sd Spread (pp)                    %5.2f      %5.2f' , 2.1, std(sc))
% fprintf('\n Mean Spread (pp)                    %5.2f      %5.2f' , 2.1, std(sc))
fprintf('\nSd(c)/sd(y)                        %5.2f      %5.2f' , 1.6, std(cc)./std(yc))
fprintf('\n corr(c,y)                         %5.2f      %5.2f' , 1.0, corr(cc,yc))
fprintf('\n corr(Spread,y)                    %5.2f      %5.2f' , -0.8, corr(sc,yc))

fprintf('\n')
