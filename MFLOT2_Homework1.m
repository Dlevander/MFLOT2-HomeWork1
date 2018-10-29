close all;
clc;
clear all;
tic
%Donnee:

b = 50e-3; % largeur du canal [m]
d_e = 15e-3; % hauteur de la sortie du canal [m]
h_e = 0.5*d_e; % [m]
d_t = 6e-3; % hauteur de la gorge [m]
h_t = 0.5*d_t; % [m]
L_c = 30e-3; % longueur partie convergente [m]
L_d = 90e-3; % longueur partie divergente [m]
L = 270e-3; % longueur du canal [m]
alpha = 3.25*pi/180; % angle de pente divergente [rad]
r1 = 254.3e-3; % rayon de courbure de la gorge [m] 
r2 = 153.7e-3; % rayon de courbure fin de divergente [m]
T0 = 300; % [K]
pa = 1.01325e5; % [Pa]
At = d_t * b; % [m^2]
Ae = b * d_e; %[m^2]
S = 111; %[K]
Tref = 273.15; %[K]
muref = 1.716e-5; %[Ns/m^2]
R = 287.1;
gamma = 1.4;
%Preallocation
X = linspace(0,0.39,3900);
H = zeros(1,3900);
Mx = zeros(3,3900);
px = zeros(3,3900);
p0x = zeros(3,3900);
Tx = zeros(3,3900);
Ux = zeros(3,3900);
Rhox = zeros(3,3900);
Me = zeros(3,1);
ce = zeros(3,1);
Te = zeros(3,1);
rhoe = zeros(3,1);
ReD = zeros(3,1);
mu = zeros(3,1);
Ue = zeros(3,1);
%%%%% Modelisation Geometry nozzle %%%%%

X_p1 = L_c + r1/( sqrt( 1/(tan(alpha))^2 -1 ) );
X_p2 = L_c + L_d - r2/( sqrt( 1/(tan(alpha))^2 +1 ) );
for i=1:length(X)
    if X(i)<=X_p1
        H(i) = h_t + r1 - sqrt( (r1)^2 - (X(i) - L_c)^2 );
    elseif X(i)<=X_p2 && X(i)>X_p1
        H(i) = 0.05748375*X(i) + 0.85568e-3;
    elseif X(i)<=0.12 && X(i)>X_p2
        H(i) = h_e - r2 + sqrt( (r2)^2 - (X(i) -L_c -L_d)^2);
    else
        H(i) = h_e;
    end
end

Ax = 2*H*b; % Vecteur des aires
%% Calcul %%

%cas sonic a gorge et subsonic apres
for i=1:1200
    Mx(1,i) = iterativeMachNumber(0.5,At,Ax(i),'subsonic');
end

%cas supersonic a la divergence et shock a x=0.07
Ash = 2*H(700)*b;
Ms1 = iterativeMachNumber(1.5,At,Ash,'supersonic');
Ms2 = sqrt( (1+ Ms1^(2) *(gamma-1)/2 )/(gamma * Ms1^(2) - (gamma-1)/2 ) );
A2star = Ash*Ms2*(((gamma+1)/2) / (1 + Ms2^2 * (gamma-1)/2))^((gamma+1)/(2*gamma-2));
for i=1:1200
    if i<300
        Mx(2,i) = iterativeMachNumber(0.5,At,Ax(i),'subsonic');
    elseif i == 300
        Mx(2,i) = 1;
    elseif i > 300 && i < 700
        Mx(2,i) = iterativeMachNumber(1.1,At,Ax(i),'supersonic');
    elseif i >= 700
        Mx(2,i) = iterativeMachNumber(0.5,A2star,Ax(i),'subsonic');
    end
end

%%% Cas sonic a gorge et supersonic partout avec shock a la sortie
% Ash = Ae
Ms1 = iterativeMachNumber(1.5,At,Ae,'supersonic');
Ms2 = sqrt( (1+ Ms1^(2) *(gamma-1)/2 )/(gamma * Ms1^(2) - (gamma-1)/2 ) );
A2star = Ash*Ms2*(((gamma+1)/2) / (1 + Ms2^2 * (gamma-1)/2))^((gamma+1)/(2*gamma-2));

for i=1:1200
    if i<300
        Mx(3,i) = iterativeMachNumber(0.5,At,Ax(i),'subsonic');
    elseif i == 300
        Mx(3,i) = 1;
    elseif i > 300 && i < 1200
        Mx(3,i) = iterativeMachNumber(1.5,At,Ax(i),'supersonic');
    elseif i == 1200
        Mx(3,i) = Ms2;
    end
end

%% Canal a section constante %%

% Temperature a l'entree
Tx(:,1200) = T0 ./ (1 + Mx(:,1200).^(2) .*(gamma - 1)/2 );
% Temperature sonic dans le canal 
Tstar = Tx(:,1200) .* (1 + Mx(:,1200).^2 .*(gamma - 1)./2 ) ./ ((gamma + 1)./2);

pe = pa; % hypothese

lambda1 = ones(3,1) * 0.02;
lambda2 = ones(3,1) * 0.03;

% Etat de sortie
for i=1:3
    while ( abs(lambda1(i) - lambda2(i)) > 1e-6 )
      lambda1(i) = lambda2(i);
      Me(i) = fsolve(@(Me) flambda1(L,Me,Mx(i,1200),lambda1(i)),0.5,optimoptions('fsolve','Display','off'));
      %Me(2) = fsolve(@(Me) flambda1(L,Me,Mx(2,1200),lambda1(2)),0.5,optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','off'));
      %Me(3) = fsolve(@(Me) flambda1(L,Me,Mx(3,1200),lambda1(3)),0.5,optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','off'));
      Te(i) = (Tstar(i) .* ((gamma + 1)/2))' ./ (1 + Me(i).^(2)*(gamma - 1)/2 );
      mu(i) = muref .* (Te(i)./Tref).^(3/2) * (Tref + S)./(Te(i)+S);  

      ce(i) = sqrt(gamma.*R.*Te(i));
      Ue(i) = Me(i).*ce(i);
      rhoe(i) = pe./(R.*Te(i));
      ReD(i) = rhoe(i).*Ue(i).*d_e./mu(i);

      lambda2(i) = fsolve(@(x) flambda2(x,ReD(i)),lambda1(i),optimoptions('fsolve','Display','off'));
      %lambda2(2) = fsolve(@(x) flambda2(x,ReD(2)),lambda1(2),optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','off'));
      %lambda2(3) = fsolve(@(x) flambda2(x,ReD(3)),lambda1(3),optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','off'));

    end
end

lambda = lambda2; %hypothese lambda constant
Ux(:,3900) = Ue;
Mx(:,3900) = Me;
Tx(:,3900) = Te;
px(:,3900) = pe;
Rhox(:,3900) = rhoe;
p0x(:,3900) = px(:,3900) .* (T0./Tx(:,3900)).^(gamma/(gamma-1));

% Etat star ((gamma + 1)./2)./(1 + Me.^2 .*(gamma - 1)./2 ).^(1./2)
Ustar = (Ux(:,3900) ./ Mx(:,3900)) ./ (Tx(:,3900)./Tstar).^(1/2);
Rhostar = Rhox(:,3900) .* Ux(:,3900)./Ustar;
pstar = px(:,3900).*Mx(:,3900).*(Tx(:,3900)./Tstar).^(1/2);
p0star = p0x(:,3900) .* Mx(:,3900) ./ (Tstar./Tx(:,3900)).^((gamma+1)/(2*gamma-2));

% Reste du canal
for i =1201:3899
    L = (3900-i)*1e-4;
    Mx(1,i) = fsolve(@(x) flambda1(L,Mx(1,3900),x,lambda),0.5,optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','off'));
    Mx(2,i) = fsolve(@(x) flambda1(L,Mx(2,3900),x,lambda),0.5,optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','off'));
    Mx(3,i) = fsolve(@(x) flambda1(L,Mx(3,3900),x,lambda),0.5,optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','off'));
end
Tx(:,1201:3899) = Tstar .* ((gamma + 1)./2)./(1 + Mx(:,1201:3899).^2 .*(gamma - 1)./2 );
Ux(:,1200:3899) = Ustar .* Mx(:,1200:3899) .* (Tx(:,1200:3899)./Tstar).^(1/2);
Rhox(:,1200:3899) = (Rhostar./Mx(:,1200:3899)) .* (Tx(:,1200:3899)./Tstar).^(-1/2);
p0x(:,1200:3899) = (p0star./Mx(:,1200:3899)) .* (Tx(:,1200:3899)./Tstar).^(-(gamma+1)/(2*gamma-2));
px(:,1200:3899) = p0x(:,1200:3899)./ (T0./Tx(:,1200:3899)).^(gamma/(gamma -1));


Qm = Rhox(:,3900) .* Ux(:,3900) .* Ae;

%% Tuyere
% On determine les p0 partout
p0x(:,1:1199) = p0x(:,1200)* ones(1,1199);
p0x(2,699) = p0x(2,700) / ( (((gamma+1)/2) / (gamma*Mx(2,699)^2 - ((gamma-1)/2)) )^(1/(gamma-1)) ...
    * ( (((gamma+1)/2)*Mx(2,699)^2) / (1 + Mx(2,699)^2 *(gamma-1)/2))^(gamma/(gamma-1)) );
p0x(2,1:699) = p0x(2,699)* ones(1,699);
p0x(3,1199) = p0x(3,1200) / ( (((gamma+1)/2) / (gamma*Mx(3,1199)^2 - ((gamma-1)/2)) )^(1/(gamma-1)) ...
    * ( (((gamma+1)/2)*Mx(3,1199)^2) / (1 + Mx(3,1199)^2 *(gamma-1)/2))^(gamma/(gamma-1)) );
p0x(3,1:1199) = p0x(3,1199)* ones(1,1199);
px(:,1:1200) = p0x(:,1:1200) ./ ((1 + (gamma-1)/2 .* Mx(:,1:1200).^2).^(gamma/(gamma-1)));
Tx(:,1:1200) = T0 ./ (1 + (gamma-1)/2 .* Mx(:,1:1200).^2);
Rho0 = p0x(:,1)./ (R*T0);
Rhox(:,1:1199) = Rho0 ./((T0./Tx(:,1:1199)).^(1/(gamma-1)));
Ux(:,1:1200) = Mx(:,1:1200) .* sqrt(gamma.*R.*Tx(:,1:1200));
p0x = p0x/101325;
px = px/101325;


%% Plot part

figure;
hold on
plot(X,H,'LineWidth',2)
grid on
title('Geometry');
ylim([0 8e-3])
xlabel('Position [m]')
ylabel('Hauteur [m]')
hold off

figure;
title('M(x)');
grid on
hold on 
plot(X,Mx(1,:),'r','LineWidth',2);
plot(X,Mx(2,:),'g','LineWidth',2);
plot(X,Mx(3,:),'b','LineWidth',2);
xlabel('Position [m]')
ylabel('Nombre de Mach []')
legend('Sans Shock','shock a 0.07','shock a 0.12');
hold off

figure;
title('p(x)');
grid on
hold on 
plot(X,px(1,:),'r','LineWidth',2);
plot(X,px(2,:),'g','LineWidth',2);
plot(X,px(3,:),'b','LineWidth',2);
xlabel('Position [m]')
ylabel('Pression [bar]')
legend('Sans Shock','shock a 0.07','shock a 0.12');
hold off

figure;
title('p0(x)');
grid on
hold on 
plot(X,p0x(1,:),'r','LineWidth',2);
plot(X,p0x(2,:),'g','LineWidth',2);
plot(X,p0x(3,:),'b','LineWidth',2);
xlabel('Position [m]')
ylabel('Pression totale [bar]')
legend('Sans Shock','shock a 0.07','shock a 0.12');
hold off

figure;
title('p(x)/p0(x)');
grid on
hold on 
plot(X,px(1,:)./p0x(1,:),'r','LineWidth',2);
plot(X,px(2,:)./p0x(2,:),'g','LineWidth',2);
plot(X,px(3,:)./p0x(3,:),'b','LineWidth',2);
xlabel('Position [m]')
ylabel('p/p0 []')
legend('Sans Shock','shock a 0.07','shock a 0.12','Location','southeast');
hold off

figure;
title('T(x)');
grid on
hold on 
plot(X,Tx(1,:),'r','LineWidth',2);
plot(X,Tx(2,:),'g','LineWidth',2);
plot(X,Tx(3,:),'b','LineWidth',2);
xlabel('Position [m]')
ylabel('Temperature [K]')
legend('Sans Shock','shock a 0.07','shock a 0.12');
hold off

figure;
title('Rho(x)');
grid on
hold on 
plot(X,Rhox(1,:),'r','LineWidth',2);
plot(X,Rhox(2,:),'g','LineWidth',2);
plot(X,Rhox(3,:),'b','LineWidth',2);
xlabel('Position [m]')
ylabel('Masse volumique [kg/m^3]')
legend('Sans Shock','shock a 0.07','shock a 0.12');
hold off

toc

function f = fanno(M) % Fonction f(M)
    gamma = 1.4;
    f = (1./gamma.*((1-M.^2)./(M.^2) + (gamma+1)./2 .* log(((gamma+1).*M.^(2) ./2) ./ (1 + M.^2 .* (gamma-1)./2))));
end
function f = flambda1(L,Me,Mi,lambda1) % Chercher le lambda en fonction de f(Me)
    d_e = 15e-3;
    f = (lambda1/2) * (L/d_e) - fanno(Mi) + fanno(Me);
end
function Lambda_Col = flambda2(lambda,ReD) % Colebrooke
    Lambda_Col = 1./sqrt(lambda) + 3.*log10((2.03./ReD).*(1./sqrt(lambda)));
end
function [Mx] = iterativeMachNumber(M_init , At , Ax , mode)
    gamma = 1.4;
    if strcmp(mode,'subsonic')
        f = @(M)( (At/Ax) *( ((gamma+1)/2) / (1 + (gamma-1)/2 * M^2) )^(-(gamma+1)/(2*gamma - 2)) -M );
    elseif strcmp(mode,'supersonic')
        f = @(M)( sqrt( ( ((gamma+1)/2) * (At/(M*Ax))^(-2*(gamma-1)/(gamma+1)) - 1 ) * 2/(gamma-1) ) -M );

    else
        error('ERROR : wrong parameter "mode" given to the function');
    end
    Mx = fsolve(@(M) f(M),M_init,optimoptions('fsolve','Display','off'));
end


