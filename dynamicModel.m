%% RBE502 Final Project
%% Trajectory Generation
%Quintic Trajectory Generation
%Points
P0=[0;0;0];P1=[0;0;1];P2=[1;0;1];
P3=[1;1;1];P4=[0;1;1];P5=[0;0;1];

%Calc Coefficents - P1-5 have the same A mat, P01 has different A mat
t0 = 0;
tfP01 = 5;
tf = 15;

counter = 1;
for t = 0:.1:5
    trajPoint = calcTraj(t0,tfP01,t,P1,P2);
    xLog(:,counter) = trajPoint(:,1);
    yLog(:,counter) = trajPoint(:,2);
    zLog(:,counter) = trajPoint(:,3);
    counter = counter + 1;
end


% plot4(0:.1:5,xLog.',{'X Pos','X Vel','X Accel'})
% figure;
% plot3([1,0,1],[1 0 1],[1 0 1])
% grid on

%% Simulation
global U;
U = [0;0;0;0]

T = [0, 5];
y0 = [P0.',0,0,0,0,0,0,0,0,0];

[t,y] = ode45(@(t,x) odeFunc(t,x,P0,P1,T(2)),T,y0);

plot4(t,[y(:,3),y(:,9)],{'Z Pos','Z Vel'})



function dX = odeFunc(t,X,P0,PF,runTime)
% Constants
m = .027;L=.046;Ix=16.571710E-6;Iy=16.571710E-6;Iz=29.261652E-6; %L & m might need to be converted back to mm/g as given
Ip=12.65625E-8;Kf=1.28192E-8;Km=5.964552E-3;Wmax=2618;Wmin=0;t0=0;g=9.81;
dX = zeros(12,1);
X =  num2cell(X);
global U;
boundary = .1;

[x, y, z, phi, theta, psi, xDot, yDot, zDot, phiDot, thetaDot, psiDot] = deal(X{:});

% Set Gains
kp = 100;
kd = 5;
K = [10;140;140;25];
lambda = [5;13;13;5];

desiredPts = zeros(3,6);
desiredPts(1:3,1:3) = calcTraj(t0,runTime,t,P0,PF); 

% Calculate Omega based off of previous U input
% Allocation Matrix
allocMat = [1/(4*Kf*L), -sqrt(2)/(4*Kf*L), -sqrt(2)/(4*Kf*L), -1/(4*Km*Kf);
            1/(4*Kf*L), -sqrt(2)/(4*Kf*L), sqrt(2)/(4*Kf*L), 1/(4*Km*Kf);
            1/(4*Kf*L), sqrt(2)/(4*Kf*L), sqrt(2)/(4*Kf*L), -1/(4*Km*Kf);
            1/(4*Kf*L), sqrt(2)/(4*Kf*L), -sqrt(2)/(4*Kf*L), 1/(4*Km*Kf)];

Wdesired = sqrt(allocMat*U);
Omega = Wdesired(1) - Wdesired(2) + Wdesired(3) - Wdesired(4);

% Z Control Law
eZ = [z-desiredPts(1,3);zDot-desiredPts(2,3)];
sZ = eZ(2) + lambda(1)*eZ(1); 
satZ = sat(sZ,boundary);
UrZ = -K(1) * satZ;
U(1) = UrZ;
U(1) = m * (-lambda(1) * eZ(2) + desiredPts(3,3) + g - K(1) * satZ)/(cos(phi) * cos(theta));

% Calculate X and Y forces
Fx = m*(-kp*(x-desiredPts(1,1)) - kd*(xDot-desiredPts(1,2)) + desiredPts(1,3));
Fy = m*(-kp*(y-desiredPts(2,1)) - kd*(yDot-desiredPts(2,2)) + desiredPts(2,3));

%Update desired points
thetaDesired = asin(Fx/U(1));
phiDesired = asin(Fy/U(1));
psiDesired = 0;
desiredPts(1,4:6) = [phiDesired,thetaDesired,psiDesired];

%Phi control law
ePhi = [phi-desiredPts(1,4);phi-desiredPts(2,4)];
sPhi = ePhi(2) + lambda(2)*ePhi(1); 
satPhi = sat(sPhi,boundary);
UrPhi = -K(2) * satPhi;
U(2) = Ix*UrPhi - (lambda(2)*phiDot*Ix) - (Ip*Omega*thetaDot) - (thetaDot*psiDot*(Iy - Iz));  


%Theta control law
eTheta = [theta-desiredPts(1,5);theta-desiredPts(2,5)];
sTheta = eTheta(2) + lambda(3)*eTheta(1); 
satTheta = sat(sTheta,boundary);
UrTheta = -K(3) * satTheta;
U(3) = UrTheta*Iy - (lambda(3)*thetaDot*Iy) - (Ip*Omega*phiDot) - (phiDot * psiDot * (Iz - Ix)); 


%Psi control law
ePsi = [psi-desiredPts(1,6);psi-desiredPts(2,6)];
sPsi = ePsi(2) + lambda(4)*ePsi(1); 
satPsi = sat(sPsi,boundary);
UrPsi = -K(4) * satPsi;
U(4) = -Iz * (UrPsi + lambda(4)*ePsi(2));


% Acclerations
Xdd = [(1/m)*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*U(1);
        (1/m)*(cos(phi)*sin(theta)*cos(psi) - sin(phi)*sin(psi))*U(1);
        (1/m)*(cos(phi)*cos(theta))*U(1) - g;
        thetaDot*psiDot*((Iy-Iz)/Ix) - (Ip/Ix)*Omega*thetaDot + (1/Ix)*U(2);
        phiDot*psiDot*((Iz-Ix)/Iy) + (Ip/Iy)*Omega*phiDot + (1/Iy)*U(3);
        phiDot*thetaDot*((Ix-Iy)/Iz) + (1/Iz)*U(4)];
U
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
        subplot(2,2,i)
        plot(t,y(:,i))
        grid on
        title(titles(i) + " over Time")
        ylabel(titles(i))
        xlabel("Time (s)")
        
    end

end












