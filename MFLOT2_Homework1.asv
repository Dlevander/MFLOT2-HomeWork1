close all;
clc;
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
pa = 1.01325; % [bar]
At = d_t * b; % [m^2]
Ae = b * d_e; %[m^2]
S = 111; %[K]
Tref = 273.15; %[K]
muref = 1.716e-5; %[Ns/m^2]
R = 287.1;
gamma = 1.4;
Ax = 2*H*b; % Vecteur des aires

X = linspace(0,0.39,3900);
H = zeros(1,3900);
Mx = zeros(3,3900);
p = zeros(3,3900);
p0 = zeros(3,3900);
Tx = zeros(3,3900);
Ux = zeros(3,3900);
Rhox = zeros(3,3900);


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

% plot nozzle
figure;
title('Geometry');
plot(X,H)
grid on
%xlim([0 0.12])
ylim([0 8e-3])

%% Calcul %%

%cas sonic a gorge et subsonic apres
Me = iterativeMachNumber(0.5,At,Ae,'subsonic');
for i=1:length(Ax)
    Mx(1,i) = iterativeMachNumber(0.5,At,Ax(i),'subsonic');
end

%cas supersonic a la divergence et shock a x=0.07
Ash = 2*H(700)*b;
Ms1 = iterativeMachNumber(1.5,At,Ash,'supersonic');
Ms2 = sqrt( (1+ Ms1^(2) *(gamma-1)/2 )/(gamma * Ms1^(2) - (gamma-1)/2 ) );
A2star = Ash*Ms2*(((gamma+1)/2) / (1 + Ms2^2 * (gamma-1)/2))^((gamma+1)/(2*gamma-2));
for i=1:length(Ax)
    if i<300
        Mx(2,i) = iterativeMachNumber(0.5,At,Ax(i),'subsonic');
    elseif i == 300
        Mx(2,i) = 1;
    elseif i > 300 && i <= 700
        Mx(2,i) = iterativeMachNumber(1.1,At,Ax(i),'supersonic');
    elseif i > 700
        Mx(2,i) = iterativeMachNumber(0.5,A2star,Ax(i),'subsonic');
    end
end

%%% Cas sonic a gorge et supersonic partout avec shock a la sortie
% Ash = Ae
Ms1 = iterativeMachNumber(1.5,At,Ae,'supersonic');
Ms2 = sqrt( (1+ Ms1^(2) *(gamma-1)/2 )/(gamma * Ms1^(2) - (gamma-1)/2 ) );
A2star = Ash*Ms2*(((gamma+1)/2) / (1 + Ms2^2 * (gamma-1)/2))^((gamma+1)/(2*gamma-2));

for i=1:length(Ax)
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

Me = 1; % GROSSE HYPOTHESE DE DEPART

pe = pa;
lambda1 = (fanno(Mi) - fanno(Me))*2*d_e/(0.27);
Te = T0/(1 + Me^2 * (gamma-1)/2);
ce

% Temperature
Tx(1,1200) = T0 / (1 + Mx(1,1200)^2 *(gamma - 1)/2 );
Tx(2,1200) = T0 / (1 + Mx(2,1200)^2 *(gamma - 1)/2 );
Tx(3,1200) = T0 / (1 + Mx(3,1200)^2 *(gamma - 1)/2 );

% Vitesse
C = sqrt(gamma * R .* Tx);
Ustar = C;
Ux = Mx .* sqrt(Tx./Tstar) .*Ustar;

% Rho
Rho0 = p0./(R*T0);
Rhox = Rho0 .* ((T0./Tx).^(-1./(gamma-1)));

% figure;
% title('Rho(x)');
% grid on
% hold on 
% plot(X,Rhox(1,:),'r');
% plot(X,Rhox(2,:),'g');
% plot(X,Rhox(3,:),'b');
% hold off
% 
% figure;
% title('M(x)');
% grid on
% hold on 
% plot(X,Mx(1,:),'r');
% plot(X,Mx(2,:),'g');
% plot(X,Mx(3,:),'b');
% hold off
% 
% figure;
% title('T(x)');
% grid on
% hold on 
% plot(X,Tx(1,:),'r');
% plot(X,Tx(2,:),'g');
% plot(X,Tx(3,:),'b');
% hold off

%Red = rho*u*d_e/mu
%mu = muref .* (Tx./Tref)^(3/2) .* (Tref + S)/(Tx+S);

function [f] = fanno(M)
    f = (1/gamma*((1-M^2)/(M^2) + (gamma+1)/2 * log(((gamma+1)*M^2/2) / (1 + M^2 * (gamma-1)/2))));
end

function [lambda] = iterativeFriction(lambda_init , Re_d)
    x = lambda_init;
    f = @(lambda)(-3.0*log10(2.03*lambda/Re_d))
    epsilon = 1; % difference entre f(i+1) et f(i)
    
    while (epsilon > 0.001) % Precision 1e-3
        y = f(x)
        epsilon = abs(y - x)
        x = y
    end
    
    lambda = x;
end

function [Mx] = iterativeMachNumber(M_init , At , Ax , mode)
    % Constantes
    gamma = 1.4;
    epsilon = 1;
    x = M_init;
    
    % Choix d'equation selon supersonique ou subsonique
    if strcmp(mode,'subsonic')
        f = @(M)( At/Ax *( ((gamma+1)/2) / (1 + (gamma-1)/2 * M^2) )^((-gamma+1)/(2*gamma - 2)) );
    elseif strcmp(mode,'supersonic')
        %f = @(M)( sqrt( (2/(gamma-1)) * ( ((gamma+1)/2) * (At/(Ax*M))^((-2*gamma+2)/(gamma+1)) - 1) ) );       
        f = @(M)( sqrt( ( ((gamma+1)/2) * (At/(M*Ax))^(-2*(gamma-1)/(gamma+1)) - 1 ) * 2/(gamma-1) ) );
    else
        error('ERROR : wrong parameter "mode" given to the function');
    end
    
    % Boucle iterative
    while (epsilon > 1e-6) % Precision 1e-6
        y = f(x);
        epsilon = abs(y - x);
        x = y;
    end
    
    % Renvoi de valeur
    Mx = x;
end


