% TODO: Rid function of it's global scope.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function : SES_dyanmics_vanilla 
% For Evaluating the SES system.
% Modified for readability, otherwise untouched
% from original authors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: R. Kett (2017)
% Original Authors: A. Richter, V. Dakos (2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = SES_dynamics_ross(~,y)

	global  a b gamma beta n p q r w1 ema K alpha mvp noiseRes noisePrice sd_resource eopt ...
            noisePriceSON   nSON   pSON   w1SON   gammaSON   betaSON   emaSON ...
            noisePriceOMNRF nOMNRF pOMNRF w1OMNRF gammaOMNRF betaOMNRF emaOMNRF ...
            harvestOMNRFTEST harvestOMNRFTESTCoop harvestOMNRFTESTDef T v eOMNRF ecOMNRFTemp ...
            initYear

    % Initialization.
	dy = zeros(3,1);
    Gamma_a = 1.5;
    Gamma_b = 0.5;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SON Segment: Mostly similar to vanilla.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % SON Actual Individual Fish Market Price.
    pActSON = pSON + noisePriceSON;

    % SON Cooperator Effort Calculation.
    ec1SON = min(emaSON, a+b*y(1)); 
    ecSON = max(0,ec1SON);

    % SON Defector Effort Calculation.
    if y(1)^alpha > w1SON/(pActSON*q)
        edSON = emaSON;
    else 
        edSON = 0;
    end

    % SON Cooperator (C) and Defector (D) Profit Calculations.
    profitCSON = pActSON*y(1)^alpha*q*ecSON + w1SON*emaSON - w1SON*ecSON;
    profitDSON = pActSON*y(1)^alpha*q*edSON + w1SON*emaSON - w1SON*edSON;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OMNRF Segment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % OMNRF Actual Individual Fish Market Price.
    pActOMNRF = pOMNRF + noisePriceOMNRF;

    % OMNRF Cooperator Effort Calculation.
    % This is the mathematic scheme presented in the paper.
    % Since I have limited matlab knowledge, I opted to continue using
    % global variables, as such at the end of each yearly range
    % we check the total harvest of OMNRF and compare it to the Coop. harvest
    % (this is done in EWS_Community through the Temp. globals.) we then set eOMNRF
    % (another temp. global) to either 0,1, or 2 for 0::below, 1::within, or 2:: above
    % the set range. 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initYear is set in EWS_Community; 
    % initYear = 1 IFF (j) -> (j+1) in EWS_Community
    % Each subsequent if-statement below sets
    % initYear = 0, thus any increment after the
    % first while solving the DE system will keep
    % ecOMNRF equal to their starting year value.
    % IE. it will always be the same effort for that
    % years quota going forward.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if initYear == 1 
        fprintf('\n initYear == %i\n', initYear)
        if eOMNRF == -1
            fprintf('eOMNRF :: %i\n', eOMNRF);
            ec1OMNRF    = min(emaOMNRF, eopt);
            ecOMNRF     = max(0, ec1OMNRF);
            ecOMNRFTemp = ecOMNRF;
            initYear = 0;
        elseif eOMNRF == 0
            fprintf('eOMNRF :: %i\n', eOMNRF);
            ec1OMNRF     = min(Gamma_a * ecOMNRFTemp, 1.0);
            ecOMNRF      = max(0, ec1OMNRF);
            ecOMNRFTemp  = ecOMNRF;
            initYear = 0;
        elseif eOMNRF == 1
            fprintf('eOMNRF :: %i\n', eOMNRF);
            ecOMNRF     = ecOMNRFTemp;
            ecOMNRFTemp = ecOMNRF;
            initYear = 0;
        elseif eOMNRF == 2
            fprintf('eOMNRF :: %i\n', eOMNRF);
            ec1OMNRF     = min(Gamma_b * ecOMNRFTemp, 1.0);
            ecOMNRF      = max(0, ec1OMNRF);
            ecOMNRFTemp  = ecOMNRF;
            initYear = 0;
        end
        fprintf('ecOMNRFTemp = %i\n[', ecOMNRFTemp)
    else
        ecOMNRF = ecOMNRFTemp;
        fprintf('|', ecOMNRF)
    end
    

    % OMNRF Defector Effort Calculation.
    if y(1)^alpha>w1OMNRF/(pActOMNRF*q)
        edOMNRF = emaOMNRF;
    else 
        edOMNRF = 0;
    end

    % OMNRF Cooperator (C) and Defector (D) Profit Calculations.
    profitCOMNRF = pActOMNRF*y(1)^alpha*q*ecOMNRF + w1OMNRF*emaOMNRF - w1OMNRF*ecOMNRF;
    profitDOMNRF = pActOMNRF*y(1)^alpha*q*edOMNRF + w1OMNRF*emaOMNRF - w1OMNRF*edOMNRF;

    % OMNRF Cooperator (C) and Defector (D) Harvest Calculations.
    % The Temps are global and interact with EWS_Community.
    % These variables are accessed by two functions at differing times, this is 
    % very bad programming.
    harvestOMNRF         = q*y(1)^alpha*edOMNRF*(nOMNRF-y(3)) + q*y(1)^alpha*ecOMNRF*y(3);
    harvestOMNRFCoop     = q*y(1)^alpha.*ecOMNRF*y(3);
    harvestOMNRFDef      = q*y(1)^alpha*edOMNRF*(nOMNRF-y(3));
    harvestOMNRFTESTCoop = [harvestOMNRFTESTCoop; harvestOMNRFCoop];
    harvestOMNRFTESTDef  = [harvestOMNRFTESTDef; harvestOMNRFDef];
    harvestOMNRFTEST     = [harvestOMNRFTEST; harvestOMNRF];
 

% Solutions.
% Restrict solution to positive orthant.
if y(1) <= 0 
    dy(1) = 0;
else + (sd_resource)*noiseRes;
    dy(1) = y(1)*(y(1)-mvp)*r*(1-(y(1)/K)) - q*y(1)^alpha*edSON*(nSON-y(2)) - q*y(1)^alpha*ecSON*y(2) - q*y(1)^alpha*edOMNRF*(nOMNRF-y(3)) - q*y(1)^alpha*ecOMNRF*y(3) + (sd_resource)*noiseRes;
end
dy(2) = (gammaSON*y(2)*(nSON-y(2)))/nSON - betaSON*y(2)*(1-profitCSON/profitDSON);
dy(3) = (gammaOMNRF*y(3)*(nOMNRF-y(3)))/nOMNRF - betaOMNRF*y(3)*(1-profitCOMNRF/profitDOMNRF);

