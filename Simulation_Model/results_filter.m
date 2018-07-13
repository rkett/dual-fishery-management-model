% Results Filter for EWS_Community output.
% Gives min, max, mean of the output given by EWS_Community in conjunction
% with either SES_Dynamics_Ross or SES_Dynamics_Vanilla.

storeMax = [];
storeMin = [];
storeMean = [];
dirString = '/Users/woof-redux/Documents/SESY_V2/data/21-Dec-2017_ROSS_run_7/';
csvFiles = dir([dirString '/*.csv']);
csvFilesCount = size(csvFiles,1);
num2str(csvFilesCount);

for fileCount = 1:csvFilesCount
    num2str(fileCount)
    resultsFile = csvread([dirString 'Ross_MC_',num2str(fileCount),'.csv']);
    storeMean = vertcat(storeMean, mean(resultsFile));
    storeMax = vertcat(storeMax, max(resultsFile));
    storeMin = vertcat(storeMin, min(resultsFile));
end

totalSimulationMean = mean(storeMean)
totalSimulationMax = max(storeMax)
totalSimulationMin = max(storeMin)


