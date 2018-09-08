%
close all
clear all
clc
%tic
% Dados iniciais
Ybus=[6.25-18.695i,-5.0+15.0i,-1.25+3.75i,0,0;
    -5+15.0i,10.83334-32.415i,-1.66667+5.0i,-1.66667+5.0i,-2.5+7.5i;
    -1.25+3.75i,-1.66667+5.0i,12.91667-38.695i,-10.0+30.0i,0;
    0,-1.66667+5.0i,-10.0+30.0i,12.91667-38.695i,-1.25+3.75i;
    0,-2.5+7.5i,0,-1.25+3.75i,3.75-11.21i;];
H_1=50.0;
H_2=50.0;
freq=60.0;
% Valores do fluxo de potência (dados)
Vreal=[1.06;1.04621;1.02032;1.01917;1.01209];
Vimag=[0;-0.05128;-0.0892;-0.09506;-0.10906];
Vbus=abs(Vreal+1i*Vimag)';
Ang=angle(Vreal+1i*Vimag)';
Pger=[1.29565,0.4,0,0,0;];
Qger=[-0.0748,0.3,0,0,0;];
% Calculando das tensoes internas e angulos de carga das maquinas
iger1=(Pger(1)-1i*Qger(1))/(Vreal(1)-1i*Vimag(1));
iger2=(Pger(2)-1i*Qger(2))/(Vreal(2)-1i*Vimag(2));
Eger1=(Vreal(1)+1i*Vimag(1))+iger1/(1/(1i*0.25));
Eger2=(Vreal(2)+1i*Vimag(2))+iger2/(1/(1i*1.50));
Delt_1=angle(Eger1);
Delt_2=angle(Eger2);
% Transformando cargas Pconstante em Zconstante
Ybus(2,2)=Ybus(2,2)+((0.20-1i*0.10)/(Vbus(2)^2));
Ybus(3,3)=Ybus(3,3)+((0.45-1i*0.15)/(Vbus(3)^2));
Ybus(4,4)=Ybus(4,4)+((0.40-1i*0.05)/(Vbus(4)^2));
Ybus(5,5)=Ybus(5,5)+((0.60-1i*0.10)/(Vbus(5)^2));
% Considerando impedancias dos geradores na Ybus
yger_1=1/(1i*0.25);
yger_2=1/(1i*1.50);
Ybus(1,1)=Ybus(1,1)+yger_1;
Ybus(2,2)=Ybus(2,2)+yger_2;
Gbus=real(Ybus);
Bbus=imag(Ybus);
% Atualizando variaveis de rede
i1=real(yger_1*Eger1);
i2=real(yger_2*Eger2);
i6=imag(yger_1*Eger1);
i7=imag(yger_2*Eger2);
Vatualiza=[Gbus -Bbus; Bbus Gbus]\([i1 i2 0 0 0 i6 i7 0 0 0]');
Vreal=Vatualiza(1:5);
Vimag=Vatualiza(6:10);
% Inicializando
iger1= yger_1*(Eger1-(Vreal(1)+1i*Vimag(1)));
iger2= yger_2*(Eger2-(Vreal(2)+1i*Vimag(2)));
Pele1=real(Eger1)*real(iger1)+imag(Eger1)*imag(iger1);
Pele2=real(Eger2)*real(iger2)+imag(Eger2)*imag(iger2);
Pele1_ant=Pele1; % Valores anteriores das potencias
Pele2_ant=Pele2;
Pmec1=Pele1;
Pmec2=Pele2;
veloc1=1.0; % Velocidades iniciais arbitradas
veloc2=1.0;
% Laço principal
t=0.00;
passo=10e-6;
t_total=1;
t_curto=0.0;
t_elim=0.1;
atual=0;
ydelt1=zeros(1,100001);
ydelt2=zeros(1,100001);
yveloc1=zeros(1,100001);
yveloc2=zeros(1,100001);
temp=zeros(1,100001);
%coder.varsize('ydelt1', 'ydelt2', 'yveloc1', 'yveloc2', 'temp' t_total/passo+1)
while ( t < t_total || t > t_total-passo/2 && t < t_total+passo/2 )
    %display(t); % Tempo atual de simulação
    % Atualizando tensoes internas
    Eger1=abs(Eger1)*cos(Delt_1)+1i*abs(Eger1)*sin(Delt_1);
    Eger2=abs(Eger2)*cos(Delt_2)+1i*abs(Eger2)*sin(Delt_2);
    % Atualizando variaveis de rede
    i1=real(yger_1*Eger1);
    i2=real(yger_2*Eger2);
    i6=imag(yger_1*Eger1);
    i7=imag(yger_2*Eger2);
    Vatualiza=[Gbus -Bbus; Bbus Gbus]\([i1 i2 0 0 0 i6 i7 0 0 0]');
    Vreal=Vatualiza(1:5);
    Vimag=Vatualiza(6:10);
    % Atualizando correntes
    iger1=yger_1*(Eger1-(Vreal(1)+1i*Vimag(1)));
    iger2=yger_2*(Eger2-(Vreal(2)+1i*Vimag(2)));
    % Atualizando potências
    Pele1_ant=Pele1;
    Pele2_ant=Pele2;
    Pele1=real(Eger1)*real(iger1)+imag(Eger1)*imag(iger1);
    Pele2=real(Eger2)*real(iger2)+imag(Eger2)*imag(iger2);
    % Atualizando variáveis de estado
    est1=[Delt_1 Delt_2 veloc1 veloc2]';
    est2=[Pmec1 Pmec2 Pele1 Pele2]';
    est3=[Pmec1 Pmec2 Pele1_ant Pele2_ant]';
    estados=newtonraphson(freq,H_1,H_2,passo,est1,est2,est3);
    Delt_1=estados(1);
    Delt_2=estados(2);
    veloc1=estados(3);
    veloc2=estados(4);
    % Evento (curto na barra 2)
    if t > t_curto-passo/2 && t < t_curto+passo/2
        Ybus(2,2)=Ybus(2,2)-1i*999999;
        Gbus=real(Ybus);
        Bbus=imag(Ybus);
    end
    if t > t_elim-passo/2 && t < t_elim+passo/2
        Ybus(2,2)=Ybus(2,2)+1i*999999;
        Gbus=real(Ybus);
        Bbus=imag(Ybus);
    end
    % Armazenando valores atuais para plotagem posterior
    atual=atual+1;
    ydelt1(atual)=Delt_1*180/pi;
    ydelt2(atual)=Delt_2*180/pi;
    yveloc1(atual)=veloc1;
    yveloc2(atual)=veloc2;
    temp(atual)=t;
    % Avançando 1 passo na simulação
    t=t+passo;
end
%toc
figure(1)
plot(temp(1:atual),ydelt1(1:atual),'LineWidth',2);
hold on;
plot(temp(1:atual),ydelt2(1:atual),'r','LineWidth',2);
ylabel('Delt (graus)');
xlabel('Tempo (s)');
%legend('G1','G2')
grid;
figure(2)
plot(temp(1:atual),yveloc1(1:atual),'LineWidth',2);
hold on;
plot(temp(1:atual),yveloc2(1:atual),'r','LineWidth',2);
ylabel('Velocidade (pu)');
xlabel('Tempo (s)');
%legend('G1','G2')
grid;