clear all;
funspath=[pwd,filesep,'data',filesep];
addpath(funspath);
%% load data to construct the relation matrix R
load('RbpGene');
R13 = Re_13;
load('RbpAS');
R14 = Re_14;
load('MiDOs')
R25 = miDOs;
load('GeneAS')
R34 = Re_34;
load('miRNAGene')
R23 = MGasso;
load('GeneDisease')
R35 = GDasso;
[npro,ndi] = size(GDasso);
load('GeneDrug')
R63 = GDrgasso;
nDrug = size(R63,1);
fun_stat = sum(R14,1);
sel_do_idx = find(fun_stat>0);
R14 = R14(:,sel_do_idx);
R34 = R34(:,sel_do_idx);
ndi = length(sel_do_idx);
Rcell={R13,R14,R23,R25,R34,R35,R63};

k1=90;
k2=170;
k3=30;
k4=10;
k5=190;
k6=70;
k = {k1,k2,k3,k4,k5,k6};   
R1 = [R13];
[G1,z1]=learnVector2(R1,k1);
[G2,z2]=learnVector2(R23,k2);
[G3,z3]=learnVector2(R13',k3);
[G4,z4]=learnVector2(R34',k4);
[G5,z5]=learnVector2(R35',k5);
[G6,z6]=learnVector2(R63,k6);

Gcell={G1,G2,G3,G4,G5,G6};

%Constraint matrix theta PPI DDI

load('HumanPPI')

load('DrugInter')
thetaCell=cell(6,1);
thetaCell{3}=PPI.*-1;
thetaCell{6}=DDI.*-1;

%instanseIdx = {3,4,5,9,11,16,17,31};  
instanseIdx = {3,4,9,11,16,17,33};% instanse index for the corresponding position in the whole relation matrix R
fprintf('begin,time:%s\n',datestr(now));
t0=clock;
nTypes=6;
newF_RBP_AS_sparse = DFSMTF_Demo_RBP_AS(nTypes,instanseIdx,Gcell,Rcell,thetaCell);
save newF_RBP_AS_sparse newF_RBP_AS_sparse;


