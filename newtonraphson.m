function [z]=newtonraphson(freq,H_1,H_2,passo,est1,est2,est3)
% Inicializndo variaveis de estado
z=est1;
% Matrizes constantes
ident4x4=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
Ma=[0 0 2*pi*freq 0; 0 0 0 2*pi*freq; 0 0 0 0; 0 0 0 0];
Mb=[0 0 0 0; 0 0 0 0; 1/(2*H_1) 0 -1/(2*H_1) 0; 0 1/(2*H_2) 0 -1/(2*H_2)];
jacobiano=2*ident4x4/passo-Ma;
% Laço
itera_max=100;
itera=1;
flag1=1;
tol=1e-10;
while ( flag1 == 1 && itera <= itera_max )
    Func1=[((z(3)-1)*2*pi*freq) ((z(4)-1)*2*pi*freq) 0 0]';
    Func2=[((est1(3)-1)*2*pi*freq) ((est1(4)-1)*2*pi*freq) 0 0]';
    result=Func1+Func2-(2*ident4x4/passo)*(z-est1)+Mb*(est2+est3);
    jac_res=jacobiano\result;
    % Atualizando variaveis de estado
    z=z+jac_res;
    % Verificando flag1
    flag1=0;
    for c=1:4
        if( abs(jac_res(c)) > tol )
            flag1=1;
        end
    end
    itera=itera+1;
end
end