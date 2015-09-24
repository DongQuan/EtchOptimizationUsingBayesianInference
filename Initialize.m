function [] = Initialize()

%Initialize global variables

global PlasmaUnknowns;
global mean;
global sd;
global data;
global ExpParameters;
global NoExpParameters;
global kb;
global q_e;
global m;
global Ncount;
global TeSet;
global PlasmaSet;
global NoTestCases;
global NoTrainingCases;


NoExpParameters = 7;
Ncount = 0;
kb = 1.381*10^-23;
q_e = 1.6e-19;
m = 35.45*1.660468e-27; 
NoUnknowns = 15; %7 k's (A and B coeff) plus standard error

%Define Synthetic data
AllExpParameters = [0.5	303	20	700	0	1	0
0.5	303	20	700	0	0.9	0;
0.5	303	20	700	0	0.8	0;
0.5	303	20	700	0	0.7	0;
0.5	303	20	700	0	0.6	0;
0.5	303	20	700	0	0.5	0;
0.5	303	20	700	0	0.4	0;
0.5	303	20	700	0	0.3	0;
0.5	303	20	700	0	0.2	0;
0.5	303	20	700	0	0.1	0;
0.5	303	20	700	0	0	0;
0.5	303	20	700	0	0.8	0;
1	303	20	700	0	0.8	0;
1.5	303	20	700	0	0.8	0;
2	303	20	700	0	0.8	0;
2.5	303	20	700	0	0.8	0;
3	303	20	700	0	0.8	0;
2	303	20	400	0	1	0;
2	303	20	400	0	0.2	0;
2	303	20	400	0	0.8	0;
2	303	20	400	0	0	0;
2	303	20	700	0	1	0;
2	303	20	700	0	0.9	0;
2	303	20	700	0	0.8	0;
2	303	20	700	0	0.7	0;
2	303	20	700	0	0.6	0;
2	303	20	700	0	0.5	0;
2	303	20	700	0	0.4	0;
2	303	20	700	0	0.3	0;
2	303	20	700	0	0.2	0;
2	303	20	700	0	0.1	0;
2	303	20	700	0	0	0];
TeSet = [3.65 3.49 3.25 3.2 3.05 3 2.95 2.9 2.85 2.83 2.8 3.25 2.75 2.5 2.4 2.3 2.25 2.75 2.6 2.4 2.23 2.78 2.65 2.55 2.45 2.3 2.26 2.24 2.23 2.21 2.2 2.2];
%ExpParameters = [1.33 303 20 700 0 0 300;1.33 303 20 700 0 1 300; 1.33 303 20 800 0 .7 250;.33 303 20 800 0 .7 300; 1.33 303 20 700 0 .5 300; 1.33 303 20 700 0 .7 300];
SynData = zeros(size(ExpParameters,1),1);

% for ExpNo = 1:size(ExpParameters,1)
%     SynData(ExpNo) = SimEtch(mean,ExpNo);
% end

AllData = [140.230972761902;130.410857280232;115.212188610325;112.233750440762;102.616977408830;99.5317573706136;98.2774194559154;92.6430362142662;89.3688919593096;87.8189438757718;85.9200553557556;115.212188610325;129.948827258397;102.711297162993;91.4399518046163;79.9332744171053;74.0964936062844;129.984600731042;113.783612812632;91.4396836679203;71.7117595603738;133.171616588336;119.234844964471;108.280487923322;97.1100657415797;79.9217456486127;75.2457181848829;72.8925760918235;71.7121746856209;69.3436933825009;68.1555817666271;68.1555724539584];
TrainingDataIndex = randperm(size(AllExpParameters,1)+1);
TrainingDataIndex = TrainingDataIndex(1:size(AllExpParameters,1)/2)-1;
NoTrainingCases = length(TrainingDataIndex);
NoTestCases = size(AllExpParameters,1) - NoTrainingCases;
TrainingData = zeros(NoTrainingCases,1);
TrainingExp = zeros(NoTrainingCases,NoExpParameters);
TestData = zeros(NoTestCases,1);
TestExp = zeros(NoTestCases,NoExpParameters);
for i=1:NoTrainingCases
    TrainingData(i) = AllData(TrainingDataIndex(i));
    TrainingExp(i,:) = AllExpParameters(TrainingDataIndex(i),:);
end
for i=1:NoTestCases
    if i~=(TrainingDataIndex(i))
    TestData(i) = AllData(i);
    TestExp(i,:) = AllExpParameters(i,:);
    end
end
%Merge training and test data to one matrix
data = cat(1,TrainingData,TestData);
ExpParameters = cat(1,TrainingExp,TestExp);
end
