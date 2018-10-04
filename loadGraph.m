function [A,V,E] = loadGraph(filename)
X = importdata(filename);
Z = X.data(3:end,1:2);
E = size(Z,1);
V = max(max(Z));
temp = ones(size(Z(:,1)));
Z =Z+1;
A = sparse(Z(:,2), Z(:,1) , temp);
end