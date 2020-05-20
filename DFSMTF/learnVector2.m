function [newX,Z,dim] = learnVector2(LabelMat,dim)

Q=LabelMat;

[U,S,V]=svd(Q);

line_S = diag(S);
sum_S = sum(line_S);

U = U(:,1:dim);
S = diag(line_S(1:dim));
V = V(:,1:dim);

percentage = sum(sum(S))/sum_S;
fprintf('this dim:%d  own percent: %-10.4f\n',dim,percentage);
X = U*sqrt(S);
Z = V*sqrt(S);



X = mapminmax(X,0,1);
minx = min(X,[],2);
maxx=max(X,[],2);
newX = X-repmat(minx,1,dim);

Z = mapminmax(Z,0,1);
