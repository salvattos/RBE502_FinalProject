%% RBE502 Final Project
%% Trajectory Generation
%Quintic Trajectory Generation
%Points
P0=[0;0;0];P1=[0;0;1];P2=[1;0;1];
P3=[1;1;1];P4=[0;1;1];P5=[0;0;1];

%Calc Coefficents - P1-5 have the same A mat, P01 has different A mat
t0 = 5;
tfP01 = 5;
tf = 15;

counter = 1;
for t = t0:.1:tf
    trajPoint = calcTraj(t0,tf,t,P1,P2);
    xLog(:,counter) = trajPoint(:,1);
    yLog(:,counter) = trajPoint(:,2);
    zLog(:,counter) = trajPoint(:,3);
    counter = counter + 1;
end



%% Simulation
global U;
U = [0;0;0;0]
global w;
w = [];
w(1) = 0;
T = [0, 10];
y0 = [P1.',0,0,0,0,0,0,0,0,0];

[t,y] = ode45(@(t,x) odeFunc(t,x,P1,P2,T(2)),T,y0);


plot4(t,[y(:,1),y(:,2),y(:,3),y(:,7),y(:,8),y(:,9)],{'X Pos','Y Pos','Z Pos','x Vel', 'y Vel', 'Z Vel'})
figure;
w = real(w)
plot(w)

function dX = odeFunc(t,X,P0,PF,runTime)
% Constants
m = .027;L = .046;Ix=16.571710E-6;Iy=16.571710E-6;Iz=29.261652E-6; %L & m might need to be converted back to mm/g as given
Ip=12.65625E-8;Kf=1.28192E-8;Km=5.964552E-3;Wmax=2618;Wmin=0;t0=0;g=9.81;
dX = zeros(12,1);
X =  num2cell(X);
global U;

[x, y, z, phi, theta, psi, xDot, yDot, zDot, phiDot, thetaDot, psiDot] = deal(X{:});

% Set Gains
kp = 100;
kd = 10;
K = [10;140;140;25];
lambda = [5;13;13;5];
boundary = [10;1.5;1.5;1];

desiredPts = zeros(3,6);
desiredPts(1:3,1:3) = calcTraj(t0,runTime,t,P0,PF); 

% Calculate Omega based off of previous U input
% Allocation Matrix
allocMat = [1/(4*Kf), -sqrt(2)/(4*Kf*L), -sqrt(2)/(4*Kf*L), -1/(4*Km*Kf);
            1/(4*Kf), -sqrt(2)/(4*Kf*L), sqrt(2)/(4*Kf*L), 1/(4*Km*Kf);
            1/(4*Kf), sqrt(2)/(4*Kf*L), sqrt(2)/(4*Kf*L), -1/(4*Km*Kf);
            1/(4*Kf), sqrt(2)/(4*Kf*L), -sqrt(2)/(4*Kf*L), 1/(4*Km*Kf)];

Wdesired = allocMat*U
Omega = Wdesired(1) - Wdesired(2) + Wdesired(3) - Wdesired(4);

global w;
w(end+1) = Wdesired(1);

if(Omega ~= 0)
    str = 'not zero!'
end
% Z Control Law
eZ = [desiredPts(1,3)-z;desiredPts(2,3)-zDot]; %might need to switch
sZ = eZ(2) + lambda(1)*eZ(1); 
satZ = sat(sZ,boundary(1));
UrZ = K(1) * satZ;
U(1) = (m/(cos(phi)*cos(theta)))*(desiredPts(3,3)+g+lambda(1)*eZ(2) + UrZ);

% Calculate X and Y forces
Fx = m*(-kp*(x-desiredPts(1,1)) - kd*(xDot-desiredPts(2,1)) + desiredPts(3,1));
Fy = m*(-kp*(y-desiredPts(1,2)) - kd*(yDot-desiredPts(2,2)) + desiredPts(3,2));

%Update desired points
thetaDesired = asin(Fx/U(1));
phiDesired = asin(-Fy/U(1));
psiDesired = 0;
desiredPts(1,4:6) = [phiDesired,thetaDesired,psiDesired];

%Phi control law
ePhi = [desiredPts(1,4)-phi;desiredPts(2,4)-phiDot];
sPhi = ePhi(2) + lambda(2)*ePhi(1); 
satPhi = sat(sPhi,boundary(2));
UrPhi = (K(2) + (Ip/Ix)*Omega)* satPhi; %try omega hat
U(2) = Ix*desiredPts(3,4) - thetaDot*psiDot*(Iy-Iz) + Ip*Omega*thetaDot + Ix*lambda(2)*ePhi(2) + Ix*UrPhi;


%Theta control law
eTheta = [desiredPts(1,5)-theta;desiredPts(2,5)-thetaDot];
sTheta = eTheta(2) + lambda(3)*eTheta(1); 
satTheta = sat(sTheta,boundary(3));
UrTheta = (K(3)+(Ip/Iy)*Omega)*satTheta;
U(3) = Iy*desiredPts(3,5) - phiDot*psiDot*(Iz-Ix) - Ip*Omega*phiDot+ Iy*lambda(3)*eTheta(2) + Iy*UrTheta;


%Psi control law
ePsi = [desiredPts(1,6)-psi;desiredPts(2,6)-psiDot];
sPsi = ePsi(2) + lambda(4)*ePsi(1); 
satPsi = sat(sPsi,boundary(4));
UrPsi = K(4) * satPsi;
U(4) = Iz*desiredPts(3,6) - phiDot*psiDot*(Ix-Iy) + Iz*lambda(4)*ePsi(2) + Iz*UrPsi;


% Acclerations
Xdd = [(1/m)*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*U(1);
        (1/m)*(cos(phi)*sin(theta)*cos(psi) - sin(phi)*cos(psi))*U(1);
        (1/m)*(cos(phi)*cos(theta))*U(1) - g;
        thetaDot*psiDot*((Iy-Iz)/Ix) - (Ip/Ix)*Omega*thetaDot + (1/Ix)*U(2);
        phiDot*psiDot*((Iz-Ix)/Iy) + (Ip/Iy)*Omega*phiDot + (1/Iy)*U(3);
        phiDot*thetaDot*((Ix-Iy)/Iz) + (1/Iz)*U(4)];
U;
dX(1) = xDot;
dX(2) = yDot;
dX(3) = zDot;
dX(4) = phiDot;
dX(5) = thetaDot;
dX(6) = psiDot;
dX(7:12) = Xdd;

end

function ret = sat(s,boundary)
    ret = min(max(s/boundary,-1),1);
end

function ret = calcTraj(t0,tf,currT,P0,PF)
    Amat = [1 t0 t0^2 t0^3 t0^4 t0^5;
        0 1 2*t0 3*t0^2 4*t0^3 5*t0^4;
        0 0 2 6*t0 12*t0^2 20*t0^3;
        1 tf tf^2 tf^3 tf^4 tf^5;
        0 1 2*tf 3*tf^2 4*tf^3 5*tf^4;
        0 0 2 6*tf 12*tf^2 20*tf^3];

    coEffs = [inv(Amat)*[P0(1);0;0;PF(1);0;0], inv(Amat)*[P0(2);0;0;PF(2);0;0], inv(Amat)*[P0(3);0;0;PF(3);0;0]];

    A = [1 currT currT^2 currT^3 currT^4 currT^5;
        0 1 2*currT 3*currT^2 4*currT^3 5*currT^4;
        0 0 2 6*currT 12*currT^2 20*currT^3];
    ret = [A*coEffs(:,1),A*coEffs(:,2),A*coEffs(:,3)];
end



function plot4(t,y,titles)
    for i = 1 : size(titles,2)
        subplot(3,3,i)
        plot(t,y(:,i))
        grid on
        title(titles(i) + " over Time")
        ylabel(titles(i))
        xlabel("Time (s)")
        
    end

end









