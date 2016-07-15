%Clear all previous data
clear all;
clc;

%Read input files
tdfread('NodesIF.txt');
nodes=[x y];
clear x;
clear y;
tdfread('ElementsIF.txt');
elements=[Node1 Node2 E A];
clear Node1;
clear Node2;
clear E;
clear A;
tdfread('NodalLoadsIF.txt');
nodalLoads=[Node Dof Load];
clear Node;
clear Dof;
clear Load;
tdfread('BoundaryConditionsIF.txt');
bcs=[Nodebc DoF Specified_Displacement];
clear Nodebc;
clear DoF;
clear Specified_Displacement;
clc;

%Calculation of number of nodes and elements and initialising matrices
nNodes=size(nodes,1);
nElements=size(elements,1);
alldofs=1:2*nNodes;
K=zeros(2*nNodes);
u=zeros(2*nNodes,1);
f=zeros(2*nNodes,1);

%Boundary Conditions
doffix=[];
for i=1:size(bcs,1)
    thisdof = 2*(bcs(i,1)-1)+bcs(i,2);
    doffix = [doffix thisdof];
    u(thisdof) = bcs(i,3);   
end

doffree=alldofs;
doffree(doffix)=[];%delete doffix elements from doffree

%initialise the force matrix
f(nNodes*2,1);
%Nodal Loads
for i=1:size(nodalLoads,1)
    f(2*(nodalLoads(i,1)-1)+nodalLoads(i,2)) = nodalLoads(i,3);
end


for i=1:nElements
    E=elements(i,3);
    A=elements(i,4);
    n1=elements(i,1);
    n2=elements(i,2);
    d1=2*n1-1;
    d2=2*n1;
    d3=2*n2-1;
    d4=2*n2;
    x1=nodes(n1,1);
    y1=nodes(n1,2);
    x2=nodes(n2,1);
    y2=nodes(n2,2);
    l=sqrt((x2-x1)^2 + (y2-y1)^2);
    lambdax=(x2-x1)/l;
    lambday=(y2-y1)/l;
    k=(E*A/l)*[lambdax^2 lambdax*lambday -(lambdax^2) -lambdax*lambday
        lambdax*lambday lambday^2 -(lambdax*lambday) -(lambday^2)
        -(lambdax)^2 -lambdax*lambday lambdax^2 lambdax*lambday
        -(lambdax*lambday) -(lambday)^2 lambdax*lambday lambday^2];
    K(d1,d1)=K(d1,d1)+k(1,1);
    K(d1,d2)=K(d1,d2)+k(1,2);
    K(d1,d3)=K(d1,d3)+k(1,3);
    K(d1,d4)=K(d1,d4)+k(1,4);
    
    K(d2,d1)=K(d2,d1)+k(2,1);
    K(d2,d2)=K(d2,d2)+k(2,2);
    K(d2,d3)=K(d2,d3)+k(2,3);
    K(d2,d4)=K(d2,d4)+k(2,4);
    
    K(d3,d1)=K(d3,d1)+k(3,1);
    K(d3,d2)=K(d3,d2)+k(3,2);
    K(d3,d3)=K(d3,d3)+k(3,3);
    K(d3,d4)=K(d3,d4)+k(3,4);
    
    K(d4,d1)=K(d4,d1)+k(4,1);
    K(d4,d2)=K(d4,d2)+k(4,2);
    K(d4,d3)=K(d4,d3)+k(4,3);
    K(d4,d4)=K(d4,d4)+k(4,4);
end

%Solution
u(doffree)=K(doffree,doffree)\(f(doffree)-K(doffree,doffix)*u(doffix));
f(doffix)=K(doffix,:)*u;
format long
disp(['Displacement vector:']);u
disp(['Reaction (Force vector):']);f

%Plot old shape
figure(1);
hold on;
plot(nodes(:,1),nodes(:,2),'k.')
hold on; axis equal;
for i=1:nElements
    elnodes=elements(i,1:2);
    nodexy=nodes(elnodes,:);
    plot(nodexy(:,1),nodexy(:,2),'k--')
end

%Plot new shape
Magnification=20;
newnodes=nodes+Magnification*reshape(u,2,nNodes)';
plot(newnodes(:,1),newnodes(:,2),'o', ...
    'MarkerEdgeColor','k', 'MarkerFaceColor','r', 'MarkerSize',10)
hold on; axis equal;
for i=1:nElements
    elnodes=elements(i,1:2);
    nodexy=newnodes(elnodes,:);
    plot(nodexy(:,1),nodexy(:,2),'k-','Linewidth',2)
end
title('Displacement');

