%% COBRA
initCobraToolbox
%% model I/O
% load model
model =readCbModel('/home/jihun/matlab/GEM/Uhlen2017-colon_adenocarcinoma/MODEL1707110792.xml');

% to export a model
%writeCbModel(model, fileName)

% Use of rBioNet to add RXN
%ReconstructionTool

%model.c : the objective coefficents
%model.s : stoichiometric matrix
%model.lb : lower bound on RXN rates
%model.ub : upper bound on RXN rates
%% simulation
% Setting of simulation constraints
%model = changeRxnBounds(model, rxnNameList, value, boundType);
model = changeRxnBounds(model, 'HMR_4521', 0, 'b');

%% Flux balance Analysis
%model = changeObjective(model, rxnNameList, objectiveCoeff);
model = changeObjective(model, 'HCC_biomass');
FBAsolution = optimizeCbModel(model);

%FBAsolution.stat : 1 == optimal solution is found ; 2 == the lower and uppter bounds are insufficent ;
%FBAsolution.v : a flux vector such that the optimal value of the objective function is attained
%FBAsolution.y : the vector of dual variables for the equality constraints
%FBAsolution.w : the vector of dual variables for the inequality constraints
