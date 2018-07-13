%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function : SES_dyanmics_vanilla 
% For Evaluating the SES system.
% Modified for readability, otherwise untouched
% from original authors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: R. Kett (2017)
% Original Authors: A. Richter, V. Dakos (2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = SES_dynamics_vanilla(~,y)

	global  a b gamma beta n p q r w1 ema ...
            K alpha mvp noiseRes ...
            noisePrice sd_resource

    % Initialization.
	dy = zeros(2,1);

	% Price calculation.
	p_act = p + noisePrice;
	
	% Cooperator effort calculation.
	ec1 = min(ema,a+b*y(1)); 
	ec = max(0,ec1);

	% Defector effort calculation.
	if y(1)^alpha>w1/(p_act*q)
    	ed = ema;
	else 
		ed = 0;
	end

	% Cooperator (c) and defector (d) profit calculations.
	profit_c = p_act*y(1)^alpha*q*ec + w1*ema - w1*ec;
	profit_d = p_act*y(1)^alpha*q*ed + w1*ema - w1*ed;
    
    coopHarvest = y(1)^alpha*q*ec; % Extra
    defectHarvest = y(1)^alpha*q*ed; % Extra
    
    
dy(1) = y(1)*(y(1)-mvp)*r*(1-(y(1)/K)) - q*y(1)^alpha*ed*(n-y(2)) - y(1)^alpha*q*(y(2))*ec +(sd_resource)*noiseRes; % update state variables
dy(2) = (n-y(2))*(y(2))*gamma/n - beta*(y(2))*(1-profit_c/profit_d);
