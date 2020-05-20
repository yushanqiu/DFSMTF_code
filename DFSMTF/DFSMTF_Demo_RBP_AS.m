function [newF] = MFLDA_Demo_RBP_AS_sparse_G_shan(nTypes,instanseIdx,Gcell,Rcell,thetaCell)
% DFSMTF: predicting novel RBP-AS event annotations by matrix
% factorization
% Yushan Qiu, Xiaoqing Cheng, Hao Jiang, Guoxian Yu, Wai-Ki Ching, College
% of Mathematics and Statistics, Shenzhen University. Contact
% xiaoqing9054@xjtu.edu.cn
%%%%%%%%%%
%   nTypes: the number of object types
%   instanseIdx: instanse index for the corresponding position in the whole relation matrix R
%   Gcell: the list of basis matrices G
%   Rcell: the list of relation matrices R
%   thetaCell: the list of constraint matrices theta


%%%%%%%%%%




LD=Rcell{2};
[NL,ND]=size(LD);
max_iter = 1;
alpha = 1000000;
theta_p = cell(size(thetaCell));
theta_n = cell(size(thetaCell));
G_enum = cell(size(Gcell));
G_denom = cell(size(Gcell));
epsilon=0.0000000001;
threshold = 0.00001;
lamda_G=0.01;

for ii = 1:length(Gcell)
    G_enum{ii}=0;
    G_denom{ii} =0;
end
Scell = cell(size(Gcell));
for ii = 1:length(thetaCell)
    theta = thetaCell{ii};
    t = abs(theta);
    theta_p{ii}=(t+theta)/2;
    theta_n{ii}=(t-theta)/2;
end

 for ii = 1:length(Gcell) 
    D_temp=eye(size(Gcell{ii},1),size(Gcell{ii},1));
       DD{ii}=D_temp;
      end
  
for iter = 1:max_iter
   
    for rr = 1:length(instanseIdx)
        i = fix(instanseIdx{rr}/nTypes)+1;
        j = mod(instanseIdx{rr},nTypes);
        if j ==0
                i = i-1;
                j = 6;
        end
        Gmatii = Gcell{i};
        Gmatjj = Gcell{j};

        Rmat = Rcell{rr};

        Smat =(Gmatii'*Gmatii)\Gmatii'*Rmat*Gmatjj/(Gmatjj'*Gmatjj);
        Smat(isnan(Smat))=0;
        Scell{rr}=Smat;

        temp1 = Rcell{rr}*Gcell{j}*Scell{rr}';
        temp1(isnan(temp1))=0;
        t = abs(temp1);

        temp1p = (t+temp1)/2;
        temp1n = (t-temp1)/2;

        temp2 = Scell{rr}*Gcell{j}'*Gcell{j}*Scell{rr}';
        temp2(isnan(temp2))=0;
        t = abs(temp2);

        temp2p = (t+temp2)/2;
        temp2n = (t-temp2)/2;

        temp3 = Rcell{rr}'*Gcell{i}*Scell{rr};
        temp3(isnan(temp3))=0;
        t = abs(temp3);
        t(t<0)=0;
        temp3p = (t+temp3)/2;
        temp3n = (t-temp3)/2;

        temp4 = Scell{rr}'*Gcell{i}'*Gcell{i}*Scell{rr};
        temp4(isnan(temp4))=0;
        t = abs(temp4);
        t(t<0)=0;
        temp4p = (t+temp4)/2;
        temp4n = (t-temp4)/2;


        G_enum{i} = G_enum{i}+temp1p+Gcell{i}*temp2n;
        G_denom{i}= G_denom{i}+temp1n+Gcell{i}*temp2p;

        G_enum{j} = G_enum{j}+ temp3p+Gcell{j}*temp4n;
        G_denom{j}= G_denom{j}+temp3n+Gcell{j}*temp4p;

    end

    for ii = 1:length(thetaCell)
        if ~isempty(thetaCell{ii})
            G_enum{ii} = G_enum{ii}+theta_n{ii}*Gcell{ii};
            G_denom{ii} = G_denom{ii}+theta_p{ii}*Gcell{ii};
        end
    end

 for ii = 1:length(Gcell) 

   D=DD{ii};

for jj=1:size(Gcell{ii},1)
      a=Gcell{ii}(jj,:);
    
    D(jj,jj)=1/(norm(Gcell{ii}(jj,:),2)+epsilon);
    end
    DD{ii}=D;
    d=abs(D);
    D_p=(d+D)/2;
    D_n =(d-D)/2;

    G_enum{ii} = G_enum{ii}+(lamda_G*D_n*Gcell{ii});
            G_denom{ii} = G_denom{ii}+(lamda_G*D_p*Gcell{ii});
end
    for ii = 1:length(Gcell)
        G_denom{ii}=G_denom{ii}+eps;
        factor = sqrt(G_enum{ii}./G_denom{ii});
        Gcell{ii}=Gcell{ii}.*factor;
        Gcell{ii}(isnan(Gcell{ii}))=0;
        Gcell{ii}(isinf(Gcell{ii}))=0;
    end

    %% compare the target approximation (||R15-G1S15G5'||^2) with threshold 

    result = Rcell{2}-Gcell{1}*Scell{2}*Gcell{4}';
    R = sum(sum(result.^2));
    R_iter_AS(iter)=R

      for rr=1:length(instanseIdx)
          i = fix(instanseIdx{rr}/nTypes)+1;
        j = mod(instanseIdx{rr},nTypes);
        if j ==0
                i = i-1;
                j = 6;
        end
      Gmatii = Gcell{i};
      Gmatjj = Gcell{j};
      result1 = Rcell{rr}-Gmatii*Scell{rr}*Gmatjj';
        R1 = sum(sum(result1.^2));
        mus1(rr)=R1;
    end
    sum_R=sum(mus1);
    sum_R_iter(iter)=sum_R;

    if R<threshold
        break;
    end
    Gcell_iter{iter}=Gcell;
Scell_iter{iter}=Scell;

end
   
newF=Gcell{1}*Scell{2}*Gcell{4}';
