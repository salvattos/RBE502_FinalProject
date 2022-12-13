%% RBE502 Final Project
% Constants
m = .027;L=.046;Ix=16.571710E-6;Iy=16.571710E-6;Iz=29.261652E-6; %L & m might need to be converted back to mm/g as given
Ip=12.65625E-8;Kf=1.28192E-8;Km=5.964552E-3;Wmax=2618;Wmin=0;

% Dynamic Model 
syms x xDot y yDot z zDot phi phiDot theta thetaDot psi psiDot u1 u2 u3 u4 g
X = [x; y; z; phi; theta; psi; xDot; yDot; zDot; phiDot; thetaDot; psiDot];
U = [u1;u2;u3;u4];


% Allocation Matrix
allocMat = [1/(4*Kf*L), -sqrt(2)/(4*Kf*L), -sqrt(2)/(4*Kf*L), -1/(4*Km*Kf);
            1/(4*Kf*L), -sqrt(2)/(4*Kf*L), sqrt(2)/(4*Kf*L), 1/(4*Km*Kf);
            1/(4*Kf*L), sqrt(2)/(4*Kf*L), sqrt(2)/(4*Kf*L), -1/(4*Km*Kf);
            1/(4*Kf*L), sqrt(2)/(4*Kf*L), -sqrt(2)/(4*Kf*L), 1/(4*Km*Kf)];

Wdesired = sqrt(allocMat*U);

Omega = Wdesired(1) - Wdesired(2) + Wdesired(3) - Wdesired(4);

% Acclerations
Xdd = [(1/m)*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*U(1);
        (1/m)*(cos(phi)*sin(theta)*cos(psi) - sin(phi)*sin(psi))*U(1);
        (1/m)*(cos(phi)*cos(theta))*U(1) - g;
        thetaDot*psiDot*((Iy-Iz)/Ix) - (Ip/Ix)*Omega*thetaDot + (1/Ix)*U(2);
        phiDot*psiDot*((Iz-Ix)/Iy) + (Ip/Iy)*Omega*phiDot + (1/Iy)*U(3);
        phiDot*thetaDot*((Ix-Iy)/Iz) + (1/Iz)*U(4)];

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
    trajPoint = calcTraj(t0,tfP01,t,P0,P3);
    xLog(:,counter) = trajPoint(:,1);
    yLog(:,counter) = trajPoint(:,2);
    zLog(:,counter) = trajPoint(:,3);
    counter = counter + 1;
end


plot4(0:.1:5,xLog.',{'X Pos','X Vel','X Accel'})
figure;
plot3(xLog,yLog,zLog)
grid on


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












