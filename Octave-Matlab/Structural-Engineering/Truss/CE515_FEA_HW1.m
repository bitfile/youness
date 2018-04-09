%+-----------------------------------------------------------+
%| Youness Bougteb      CE 515      Homework #1     9/6/2017 |
%+-----------------------------------------------------------+
%|         Assignment: Truss Stiffness Matrix = DONE         |
%+-----------------------------------------------------------+
close all;clc;
%[num, txt, raw] = xlsread('CE515_FEA_HW1.xlsx');
[num, txt, raw] = xlsread('CE515_FEA_HW1.xlsx');
Para=num(:,1);  Para=Para(all(~isnan(Para),2),:);
node=num(:,2:3);node=node(all(~isnan(node),2),:);
data=num(:,4:9);data=data(all(~isnan(data),2),:);

A=Para(1);E=Para(2);EA=E*A;

nelems = size(data,1);
nnodes = size(node,1);

elem = data(:,1:2);
ndof = 2*nnodes;
dof  = zeros(1,4);
Le   = zeros(nelems,1);
Kg   = zeros(ndof); % Global Stiffness Matrix
k    = zeros(4);

for i = 1:nelems
    %Le = element length
    Le(i)=sqrt((data(i,4)- data(i,3))^2 + (data(i,6)- data(i,5))^2);
    c =(data(i,4)- data(i,3))/Le(i); % Direction cosine in x direction
    s =(data(i,6)- data(i,5))/Le(i); % Direction cosine in y direction
    % Equation (2.4-6) R. D. Cook - textbook (c = cos and s = sin)
    k = (EA/Le(i))*[c^2 c*s -c^2 -c*s;
                    c*s s^2 -c*s -s^2;
                    -c^2 -c*s c^2 c*s;
                    -c*s -s^2 c*s s^2];
    % Get the nodes connecting member
    elemnode  = elem( i, 1:2);
    % Get the staring DOFs of the ith member
    dof = 2*elemnode(1)-1:2*elemnode;
    % Get the ending DOFs of the ith member + concatenate them
    dofs = [dof 2*elemnode(2)-1:2*elemnode(2)];
    % Populate global stiffness matrix
    Kg(dofs,dofs) = Kg(dofs,dofs) + k;
end