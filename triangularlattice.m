% Honeycomb Lattice Unit Example
% Script by Andrew Phillips

clear all
close all

%
% Section and material definition
%

% Assuming units are Newtons (N) and millimetres (mm)
% Assuming a rectangular cross-section
t = 1; % mm
b = 1; % unit value (mm)
E = 2000; % MPa
% Calculating cross-sectional properties
A = b*t;
I = b*t^3/12;
% Calculating structural stifness / rigidity
EA = E*A; % axial
EI = E*I; % bending
% to reduce axial deformation
EA=EA*1000;%experiment with why this line is here - comment it out in terms of %. 

%
% Mesh definition
%

% Setting up a matrix of node points for the hexagon
% Setting the number of structural members - which is 3 for triangles
NumSM = 3; 
% Setting the length of each structural member
LengthSM = 10;
% Setting the number of elements making up each structural member
NumElemSM = 8;
% Calculating the number of elements required
NumElems = NumSM*NumElemSM; %which is 24
% Calculating the number of node
NumNodes = NumElems; 
% Setting the angle alpha from the figure
alpha = 60; %change it to 60 degrees for triangular lattice. 
alpha = alpha*pi/180;
% Calculating the position of nodes at the end of structural members  
NodePoints = zeros(NumNodes,2);
NodePoints(1,:) = [0 0]; %node 1 origin 
NodePoints(1*NumElemSM+1,:) = [LengthSM*cos(alpha) -LengthSM*sin(alpha)]; %node 9
NodePoints(2*NumElemSM+1,:) = [LengthSM*cos(alpha) -LengthSM-LengthSM*sin(alpha)]; %node 17

% Calculating the position of nodes along the lengh of structural members
% using the positions of the start and end nodes of structural members
for n=1:1:NumElemSM-1
    NodePoints(n+1,:) = NodePoints(1,:)+(NodePoints(1*NumElemSM+1,:)-NodePoints(1,:)).*n/NumElemSM;%nodes 2 to 8
    NodePoints(1*NumElemSM+n+1,:) = NodePoints(1*NumElemSM+1,:)+(NodePoints(2*NumElemSM+1,:)-NodePoints(1*NumElemSM+1,:)).*n/NumElemSM; %nodes 10 to 16
    NodePoints(2*NumElemSM+n+1,:) = NodePoints(2*NumElemSM+1,:)+(NodePoints(1,:)-NodePoints(2*NumElemSM+1,:)).*n/NumElemSM; %nodes 17 to 24
end

% Setting up a matrix of element definitions (start and end nodes)
ElemNodes = zeros(NumElems,2);
ElemNodes(:,1) = 1:1:NumNodes;
ElemNodes(:,2) = [2:1:NumNodes 1];

% Calculating element lengths - change for triangular lattice.
ElemLengths = zeros(NumElems,1);
for n=1:NumElems
    ElemLengths(n) = sqrt(sum(abs(NodePoints(ElemNodes(n,2),:)-NodePoints(ElemNodes(n,1),:)).^2));
end

% Setting up a matrix of element angles (transforming from local to global
% coordinate systems) for the hexagon
ElemAngle = zeros(NumElems,1);
ElemAngle(1:NumElemSM) = -alpha;%elements 1 to 8
ElemAngle(1+1*NumElemSM:2*NumElemSM) = -pi/2; %elements 9 to 16
ElemAngle(1+2*NumElemSM:3*NumElemSM) = -(pi-alpha); %elements 17 to 24

%Building structure stiffness matrix K

% Setting up the element stiffness matrix k for each element
k = zeros(6,6,NumElems);
% Setting it up this way allows for easy adaptation for elements with
% different lengths and different cross-sectional properties
for n=1:NumElems
    k(:,:,n) = [EA/ElemLengths(n) 0 0 -EA/ElemLengths(n) 0 0;
        0 12*EI/(ElemLengths(n)^3) 6*EI/(ElemLengths(n)^2) 0 -12*EI/(ElemLengths(n)^3) 6*EI/(ElemLengths(n)^2);
        0 6*EI/(ElemLengths(n)^2) 4*EI/ElemLengths(n) 0 -6*EI/(ElemLengths(n)^2) 2*EI/ElemLengths(n);
        -EA/ElemLengths(n) 0 0 EA/ElemLengths(n) 0 0;
        0 -12*EI/(ElemLengths(n)^3) -6*EI/(ElemLengths(n)^2) 0 12*EI/(ElemLengths(n)^3) -6*EI/(ElemLengths(n)^2);
        0 6*EI/(ElemLengths(n)^2) 2*EI/ElemLengths(n) 0 -6*EI/(ElemLengths(n)^2) 4*EI/ElemLengths(n)];
end
% Transforming the stiffness matrix k for each element
% from local to global to give kdash for each element
% global) 
N = zeros(6,6,NumElems);
NQ = zeros(3,3,NumElems);
kdash = zeros(6,6,NumElems);
for n=1:NumElems
    % setting up a transformation matrix for all four quadrants
    % top-left and bottom-right quadrants of N contain trigonometric terms
    % top-right and bottom-left quadrants of N contain zero terms
    NQ(:,:,n) = [cos(ElemAngle(n)) sin(ElemAngle(n)) 0;...
        -sin(ElemAngle(n)) cos(ElemAngle(n)) 0;...
        0 0 1];
    N(1:3,1:3,n)=NQ(:,:,n);
    N(4:6,4:6,n)=NQ(:,:,n);
    kdash(:,:,n)=N(:,:,n)'*k(:,:,n)*N(:,:,n);
end

% Constructing the structure stiffness matrix K
% Noting that there are three degrees of freedom at every node
K = zeros(3*NumNodes,3*NumNodes);
% placing the top-left quadrant
for n=1:NumNodes
    clear a
    a=find(ElemNodes(:,1)==n);
    if isempty(a)~=1
        for m=1:length(a)
            K(1+(n-1)*3:n*3,1+(n-1)*3:n*3)=...
                K(1+(n-1)*3:n*3,1+(n-1)*3:n*3)+kdash(1:3,1:3,a(m));
        end
    end
end
%placing the bottom-right quadrant
for n=1:NumNodes
    clear a
    a=find(ElemNodes(:,2)==n);
    if isempty(a)~=1
        for m=1:length(a)
            K(1+(n-1)*3:n*3,1+(n-1)*3:n*3)=...
                K(1+(n-1)*3:n*3,1+(n-1)*3:n*3)+kdash(4:6,4:6,a(m));
        end
    end
end
% placing the top-right and bottom-left quadrant
for n=1:NumNodes
    clear a b
    a=find(ElemNodes(:,1)==n);
    if isempty(a)~=1
        b=ElemNodes(a,2);
        for m=1:length(a)
            K(1+(n-1)*3:n*3,1+(b(m)-1)*3:b(m)*3)=...
                K(1+(n-1)*3:n*3,1+(b(m)-1)*3:b(m)*3)+kdash(1:3,4:6,a(m));
            K(1+(b(m)-1)*3:b(m)*3,1+(n-1)*3:n*3)=...
                K(1+(b(m)-1)*3:b(m)*3,1+(n-1)*3:n*3)+kdash(4:6,1:3,a(m));
        end
    end
end

%
% Point load definition
% Loading depends on whether we are looking at behaviour
% in the global x or y direction

% Noting that at each node three loads could be applied
% Points loads in the x and y directions,
% and a moment load about the y axis

P = 1;
NodeLoads = zeros(NumNodes,1);

% Loading in the x direction
LoadedNodeA = 1*NumElemSM+1;
LoadedNodeB = 3*NumElemSM+1;
NodeLoads(LoadedNodeA*3-2:LoadedNodeA*3) = [P/2 0 0]; %[x y z directions]
NodeLoads(LoadedNodeB*3-2:LoadedNodeB*3)=[P/2 0 0];

% % Loading in the y direction
%LoadedNode = 1;
%NodeLoads(LoadedNode*3-2:LoadedNode*3) = [0 -P 0];

%
% Boundary condition definition
% Boundary conditions depend on whether we are looking a behaviour
% in the global x or y direction

% Noting that at each node three restraints could be applied
% Preventing displacement in the x and y directions,
% and preventing rotation about the y axis

% We adjust the K matrix by setting rows and columns associated
% with the restrained degree of freedom to zero
% and the intersection to 1

% boundary conditons are set by listing the node and then the three DOFs
% indicating which ones are fixed at the node with a "1"

% Rotation restrictions to replicate the hexagon being constrained with a
% honeycomb lattice network

FixRotNodes = [1*NumElemSM+1 0 0 0; 2*NumElemSM+1 0 0 1;3*NumElemSM+1 0 0 1];

% Loading in the x direction
%RollerNode = [2*NumElemSM+1 1 0]; % fixed in the x direction
%PinnedNode = [1*NumElemSM+1 1 1]; % fixed in the x and y direction%5BCNodes = [RollerNode; PinnedNode; FixRotNodes];

% % Loading in the y direction
RollerNode = [2*NumElemSM+1 1 1 0]; % fixed in the x direction
PinnedNode = [3*NumElemSM+1 1 1 0]; % fixed in the x and y direction
BCNodes = [RollerNode; PinnedNode; FixRotNodes];

% Finding the fixed DOFs
NumBCNodes = size(BCNodes,1);
FixedDOFs = zeros(1,sum(sum(BCNodes(:,2:4)~=0)));
a=0;
for n=1:NumBCNodes
    for m=2:4
        if BCNodes(n,m)==1
            a=a+1;
            FixedDOFs(a) = (BCNodes(n,1)-1)*3+(m-1);
        end
    end
end
clear a

% Adjusting the K matrix
for n=1:size(FixedDOFs,2)
    K(FixedDOFs(n),:)=0;
    K(:,FixedDOFs(n))=0;
    K(FixedDOFs(n),FixedDOFs(n))=1;
end

%
% Analysis
%

%
% Finding the node displacements
%

RawNodeDisps=(K^-1).*(NodeLoads); %only analysis step within the whole code
% Reshape the resulting matrix
NodeDisps=reshape(RawNodeDisps,1,[]);
NodeDisps=NodeDisps';
% Calculating the displaced node positions
DispNodePoints = NodePoints+NodeDisps(:,:);

%
% Finding node forces and moments
%

% Calculating beam end force and moments (in global axis system)
BeamEndForcesGAS = zeros(NumElems,6);
for n=1:NumElems
    BeamEndForcesGAS(n,:) = (kdash(:,:,n)*[NodeDisps(ElemNodes(n,1),:)'; NodeDisps(ElemNodes(n,2),:)'])';  
end
% Calculating node external forces (in global axis system)
% This gives applied loading and reaction forces
NodeExternalForcesGAS = zeros(NumNodes,3);
for n=1:NumNodes
    NodeExternalForcesGAS(n,:) = sum(BeamEndForcesGAS(ElemNodes(:,2)==n,4:6),1) + sum(BeamEndForcesGAS(ElemNodes(:,1)==n,1:3),1);
end

% Calculating beam end force and moments (in the beam axis system)
% This gives axial forces, shear forces and bending moments
BeamEndForcesLAS = zeros(NumElems,6);
for n=1:NumElems
    BeamEndForcesLAS(n,:) = (k(:,:,n)*N(:,:,n)*[NodeDisps(ElemNodes(n,1),:)'; NodeDisps(ElemNodes(n,2),:)'])';
end

%
% Plotting the displaced shape
%

% Settting a magnification factor
MagFact = 0.25*LengthSM/max(max(abs(NodeDisps(:,1:2))));
PlotDispNodePoints = NodePoints+MagFact.*NodeDisps(:,1:2);

figure(1)
hold on
plot(NodePoints([1:end 1],1),NodePoints([1:end 1],2),'o - b','MarkerFaceColor','b');
plot(PlotDispNodePoints([1:end 1],1),PlotDispNodePoints([1:end 1],2),'o - r','MarkerFaceColor','r');
title({'Displaced Shape', ['Magnification factor: ', num2str(MagFact)]});
xlabel('x'); ylabel('y');
set(gca,'DataAspectRatio',[1 1 1]);

%
% Plotting bending moment, axial force, and shear force diagrams
%

% Calculating beam direction and perpendicular vectors
ElemBeamVectors = zeros(NumElems,3);
ElemBeamNormVectors = zeros(NumElems,3);
ElemPerpVectors = zeros(NumElems,3);
for n=1:NumElems
    ElemBeamVectors(n,1:2) = NodePoints(ElemNodes(n,2),:)-NodePoints(ElemNodes(n,1),:);
    ElemBeamNormVectors(n,1:2) = ElemBeamVectors(n,1:2)./norm(ElemBeamVectors(n,1:2));
    ElemPerpVectors(n,:) = cross([0 0 1],ElemBeamNormVectors(n,:));
end

% Adopting a convention that means bending moments are plotted on the
% tension side of the element

% Setting a scale factor
BMScaleFact = 0.25*LengthSM/max(max(abs(BeamEndForcesLAS(:,[3 6]))));

figure(2)
hold on
plot(NodePoints([1:end 1],1),NodePoints([1:end 1],2),'o - k','MarkerFaceColor','k');
for n=1:NumElems
    pointA = NodePoints(ElemNodes(n,1),:);
    pointB = pointA + BeamEndForcesLAS(n,3)*ElemPerpVectors(n,1:2)*BMScaleFact;
    pointD = NodePoints(ElemNodes(n,2),:);
    pointC = pointD - BeamEndForcesLAS(n,6)*ElemPerpVectors(n,1:2)*BMScaleFact;
    plot([pointA(1) pointB(1) pointC(1) pointD(1)],...
        [pointA(2) pointB(2) pointC(2) pointD(2)],'- r');
end
title({'Bending Moment', '(plotted on the tension side)', ['Absolute maximum: ', num2str(max(max(abs(BeamEndForcesLAS(:,[3 6]))))), ' Nmm']});
xlabel('x'); ylabel('y');
set(gca,'DataAspectRatio',[1 1 1]);

% Adopting a convention that means tension axial force is plotted in red
% and compression axial force is plotted in blue

% Setting a scale factor
AFScaleFact = 0.25*LengthSM/max(max(abs(BeamEndForcesLAS(:,[1 4]))));

figure(3)
hold on
plot(NodePoints([1:end 1],1),NodePoints([1:end 1],2),'o - k','MarkerFaceColor','k');
for n=1:NumElems
    pointA = NodePoints(ElemNodes(n,1),:);
    pointB = pointA - BeamEndForcesLAS(n,1)*ElemPerpVectors(n,1:2)*AFScaleFact;
    pointD = NodePoints(ElemNodes(n,2),:);
    pointC = pointD + BeamEndForcesLAS(n,4)*ElemPerpVectors(n,1:2)*AFScaleFact;
    if -BeamEndForcesLAS(n,1) > 0 
        plot([pointA(1) pointB(1) pointC(1) pointD(1)],...
            [pointA(2) pointB(2) pointC(2) pointD(2)],'- r');
    elseif -BeamEndForcesLAS(n,1) < 0
        plot([pointA(1) pointB(1) pointC(1) pointD(1)],...
            [pointA(2) pointB(2) pointC(2) pointD(2)],'- b');
    end 
end
title({'Axial Force', '(blue: compression, red: tension)', ['Absolute maximum: ', num2str(max(max(abs(BeamEndForcesLAS(:,[1 4]))))), ' N']});
xlabel('x'); ylabel('y');
set(gca,'DataAspectRatio',[1 1 1]);

% Adopting a convention that means that positive shear forces point in the
% positive y axis direction on the beam

% Setting a scale factor
SFScaleFact = 0.25*LengthSM/max(max(abs(BeamEndForcesLAS(:,[2 5]))));

figure(4)
hold on
plot(NodePoints([1:end 1],1),NodePoints([1:end 1],2),'o - k','MarkerFaceColor','k');
for n=1:NumElems
    pointA = NodePoints(ElemNodes(n,1),:);
    pointB = pointA - BeamEndForcesLAS(n,2)*ElemPerpVectors(n,1:2)*SFScaleFact;
    pointD = NodePoints(ElemNodes(n,2),:);
    pointC = pointD + BeamEndForcesLAS(n,5)*ElemPerpVectors(n,1:2)*SFScaleFact;
    plot([pointA(1) pointB(1) pointC(1) pointD(1)],...
        [pointA(2) pointB(2) pointC(2) pointD(2)],'- r');
end
title({'Shear Force', ['Absolute maximum: ', num2str(max(max(abs(BeamEndForcesLAS(:,[2 5]))))), ' N']});
xlabel('x'); ylabel('y');
set(gca,'DataAspectRatio',[1 1 1]);

%
% Plotting the beam direction and perpendicular vectors
%

% Setting a scale factor
LASScaleFact = 0.5*max(ElemLengths);

figure(5)
hold on
plot(NodePoints([1:end 1],1),NodePoints([1:end 1],2),'o - k','MarkerFaceColor','k');
for n=1:NumElems
    pointA = NodePoints(ElemNodes(n,1),:);
    pointB = pointA+ElemBeamNormVectors(n,1:2).*LASScaleFact;
    pointC = pointA+ElemPerpVectors(n,1:2).*LASScaleFact;
    plot([pointA(1) pointB(1)], [pointA(2) pointB(2)], '- r');
    plot([pointA(1) pointC(1)], [pointA(2) pointC(2)], '- g');
    plot(pointA(1), pointA(2),'o b','MarkerFaceColor','b');
end
title('Beam Axes');
set(gca,'DataAspectRatio',[1 1 1]);

%%%% THE END %%%%


    
