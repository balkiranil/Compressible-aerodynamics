clc;
clear;
disp('Anıl Taha Balkır');
disp('---------------------------------------');

%% Initial data 
gamma=1.4;
m=3.9;
theta1=19;
theta2=6; 
%% Any input
% gamma=1.4;
% m=input('mach1=');
% theta1=input('Turn Angle 1=');
% theta2=input('Turn Angle 2=');
% disp('---------------------------------------');

if(theta1>0)
    [m3,beta3,p3_p1,ro3_ro1,t3_t1,p03_p01,ro03_ro01,t03_t01] = obliqueshock(m,theta1,gamma,0); % find pressure and mach back of first oblique shock
else
    [m3,p3_p1,ro3_ro1,t3_t1,p03_p01,ro03_ro01,t03_t01] = expansionfan(m,-theta1,gamma);
end
if(theta2>0)
    [m2,beta2,p2_p1,ro2_ro1,t2_t1,p02_p01,ro02_ro01,t02_t01] = obliqueshock(m,theta2,gamma,0); % find pressure and mach back of second oblique shock
else
    [m2,p2_p1,ro2_ro1,t2_t1,p02_p01,ro02_ro01,t03_t01] = expansionfan(m,-theta2,gamma);
end
phi_max=max(theta1,theta2)/2; % guess a range for slip line angle
for i=1:2000
    phi(i)=-phi_max+phi_max*(i-1)/500; 
    if(theta2+phi(i)>0)
        [m4(i),beta4(i),p4_p3(i),ro4_ro3(i),t4_t3(i),p04_p03(i),ro04_ro03(i),t04_t03(i)] = obliqueshock(m2,theta2+phi(i),gamma,0); 
    else
        [m4(i),p4_p3(i),ro4_ro3(i),t4_t3(i),p04_p03(i),ro04_ro03(i),t04_t03(i)] = expansionfan(m2,-theta2-phi(i),gamma);
    end
    if(theta1-phi(i)>0)
        [m4p(i),beta4p(i),p4p_p2(i),ro4p_ro2(i),t4p_t2(i),p04p_p02(i),ro04p_ro02(i),t04p_t02(i)] = obliqueshock(m3,theta1-phi(i),gamma,0); 
    else
        [m4p(i),p4p_p2(i),ro4p_ro2(i),t4p_t2(i),p04p_p02(i),ro04p_ro02(i),t04p_t02(i)] = expansionfan(m3,-theta1+phi(i),gamma);
    end
    p4(i)=p4_p3(i)*p2_p1;
    p4p(i)=p4p_p2(i)*p3_p1;
    delta_p(i)=abs(p4(i)-p4p(i)); 
end

minimum=min(delta_p);
[x,y]=find(delta_p==minimum);
fprintf('Slip line Angle (phi) : %.3f deg \n',phi(x,y))
fprintf('DeltaP                : %.7f atm \n',delta_p(x,y))
fprintf('M4                    : %.3f\n',m4(x,y))
fprintf('M4 prime              : %.3f\n',m4p(x,y))

%% figure 1:DeltaP - Phi
figure(1);
subplot(1,3,1);
figure(1);
plot(phi,delta_p,'b'),xlabel('Phi'),ylabel('delta P');
title('Delta P');
grid on;
subplot(1,3,2);
plot(phi,p4,'r'),xlabel('Phi'),ylabel('P'),hold on;
plot(phi,p4p,'b');
title('Back Pressure');
grid on;
legend('P4','P4p'); 
subplot(1,3,3);
plot(phi,m4,'r'),xlabel('Phi'),ylabel('M'),hold on;
plot(phi,m4p,'b');
title('Mach');
grid on;
legend('M4','M4p'); 

%% Prandtl Meyer Function
function nu = prandtlmeyer(mach,gamma)
    nu=(sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)*(mach^2-1)/(gamma+1)))-atan(sqrt(mach^2-1)))*180/pi;
end
%% Oblique Shock Function

function [m3,beta,p3_p1,ro3_ro1,t3_t1,p03_p01,ro03_ro01,t03_t01] = obliqueshock(m1,theta,gamma,strong)
    if(strong==1)
        beta=betastrong(m1,theta,gamma);
    else
        beta=betaweak(m1,theta,gamma);
    end
    mn1=m1*sin(beta*pi/180);
    [p01_p1,ro01_ro1,t01_t1] = isentropic(m1,gamma);
    [mn2,p3_p1,ro3_ro1,t3_t1] = normalshock(mn1,gamma);
    m3=mn2/sin((beta-theta)*pi/180);
    [p02_p2,ro02_ro2,t02_t2] = isentropic(m3,gamma);
    p03_p01=(p02_p2/p01_p1)*(p3_p1);
    ro03_ro01=(ro02_ro2/ro01_ro1)*(ro3_ro1);
    t03_t01=(t02_t2/t01_t1)*(t3_t1);
    
end
%% Normal Shock Function

function [m3,p3_p1,ro3_ro1,t3_t1] = normalshock(m,gamma)
    p3_p1=(1+(2*gamma)*(m^2-1)/(gamma+1)); % pressure ratio
    ro3_ro1=(gamma+1)*m^2/(2+(gamma-1)*m^2); % density ratio
    t3_t1=p3_p1/ro3_ro1; % temprature ratio
    m3=sqrt((1+(gamma-1)*m^2/2)/(gamma*m^2-(gamma-1)/2)); % amount of m2
end

%% Isentropic function

function [p0_p,ro0_ro,t0_t] = isentropic(m,gamma)
	p0_p = (m^2*(gamma-1)/2+1)^(gamma/(gamma-1));
	ro0_ro = (m^2*(gamma-1)/2+1)^(1/(gamma-1));
	t0_t = m^2*(gamma-1)/2+1;
end

%% Inverse Prandtl Meyer Function

function m = inverseprandtlmeyer(nu,gamma)
    error=10;
    for i=1:1000
        m_temp=1+i/100;
        nu_temp=(sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)*(m_temp^2-1)/(gamma+1)))-atan(sqrt(m_temp^2-1)))*180/pi;
        if(abs(nu_temp-nu)<error)
            m=m_temp;
            error=abs(nu_temp-nu);
        end
    end
end

%% Expansion Fan Function

function [m3,p3_p1,ro3_ro1,t3_t1,p03_p01,ro03_ro01,t03_t01] = expansionfan(m1,theta,gamma)
    nu1=prandtlmeyer(m1,gamma);
    nu2=nu1+theta;
    m3=inverseprandtlmeyer(nu2,gamma);
    [p01_p1,ro01_ro1,t01_t1] = isentropic(m1,gamma);
    [p02_p2,ro02_ro2,t02_t2] = isentropic(m3,gamma);
    p3_p1=p01_p1/p02_p2;
    ro3_ro1=ro01_ro1/ro02_ro2;
    t3_t1=t01_t1/t02_t2;
    p03_p01=1;
    ro03_ro01=1;
    t03_t01=1;
end

%% Beta Weak Function

function beta = betaweak(m,theta,gamma)
    theta=theta*pi/180;
    lambda=sqrt((m^2-1)^2-3*(1+(gamma-1)*m^2/2)*(1+(gamma+1)*m^2/2)*tan(theta)^2);
    x=((m^2-1)^3-9*(1+(gamma-1)*m^2/2)*(1+(gamma-1)*m^2/2+(gamma+1)*m^4/4)*tan(theta)^2)/(lambda^3);
    beta=real(atan((m^2-1+2*lambda*cos((4*pi+acos(x))/3))/(3*(1+(gamma-1)*m^2/2)*tan(theta)))*180/pi);
end