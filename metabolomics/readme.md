# Genome Scale Modeling
## COBRA
Basic
```
%% COBRA
initCobraToolbox

% load model
model =readCbModel('Model_Address');

% add a reaction to model
Ref : https://opencobra.github.io/cobratoolbox/latest/modules/reconstruction/refinement/index.html
model = addReaction(model,'[ReactionName]','reactionFormula', 'MetaboliteID_1 + 2 MetaboliteID_2 -> 3 MetaboliteID_3', 'subSystem', 'Artificial reactions')

```

running code
```
%% COBRA
initCobraToolbox

% load model
model =readCbModel('/home/jihun/matlab/GEM/colon.xml');

% To add a reaciton
model = addReaction(model,'bioMass','reactionFormula', 'M_m02392c[C_c] + M_m01721c[C_c] + M_m01722n[C_n] + M_m02847c[C_c] + M_m01307c[C_c] + M_m01365c[C_c] + M_m01369c[C_c] + M_m01370c[C_c] + M_m01451c[C_c] + M_m01450c[C_c] + M_m01589m[C_m] + M_m01628c[C_c] + M_m01975c[C_c] + M_m01974c[C_c] + M_m01986c[C_c] + M_m03161c[C_c] + M_m02125c[C_c] + M_m02184c[C_c] + M_m02360c[C_c] + M_m02426c[C_c] + M_m02471c[C_c] + M_m02733c[C_c] + M_m02724c[C_c] + M_m02750c[C_c] + M_m02770c[C_c] + M_m02896c[C_c] + M_m02908c[C_c] + M_m02993c[C_c] + M_m03089c[C_c] + M_m03101c[C_c] + M_m03135c[C_c] + M_m03052c[C_c] -> biomass[x]', 'subSystem', 'Artificial reactions')
model = changeObjective(model, 'bioMass');
FBAsolution = optimizeCbModel(model);
```


# GEM paper
- Tissue-based map of the human proteome 
  - Tissue-specific GEM (HMR2)
- Importance of the biomass formulation for cancer metabolic modeling and drug prediction
