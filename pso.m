clear;clc
rng default

fob = @(x) x(1).*sin(4*pi.*x(1)) - x(2).*sin(4*pi.*x(2) + pi) + 1;

lb = [-1 -1];       % limites inferiores
ub = [2 2];         % limites superiores

n = 20;             % numero de individuos
d = 2;              % numero de variaveis
w = 0.5;            % fator de inercia
c1 = 2;             % peso (pbest) 
c2 = 2;             % peso (gbest)
niter = 100;        % numero maximo de iteracoes

% Inicializa a população e avalia cada individuo

x = zeros(n,d);
fitness = zeros(n,1);

for i=1:n
   
    x(i,:) = lb + (ub-lb).*rand(1,d);
    fitness(i) = fob(x(i,:));
    
end

% Encontra gbest e pbest atuais

[fmax,index] = max(fitness);
gbest = x(index,:);
pbest = x;

t=0;                            % Inicia contador de iteracoes
conv = zeros(1,niter); 

v = zeros(n,d);                 % Inicializa velocidade das particulas

while t < niter
    
    t = t + 1;
    
    conv(t) = fmax;
    
    for i=1:n
    
        % Calcula velocidades
             
        v(i,:) = w*v(i,:) + c1*rand(1,d).*(pbest(i,:) - x(i,:)) + c2*rand(1,d).*(gbest - x(i,:));
        
        % Calcula novas posicoes
        
        x(i,:) = x(i,:) + v(i,:);
        
        % Traz as particulas para dentro do espaco de solucao
        
        for j = 1:d
           
            if x(i,j) < lb(j) 
                
                x(i,j) = lb(j);
                
            end
            
            if x(i,j) > ub(j) 
                
                x(i,j) = ub(j);
                
            end
                
        end
        
        % Atualiza a aptidao da particula
        
        fnew = fob(x(i,:));
        
        % Atualiza o pbest
        
        if fnew > fitness(i) 
           
            pbest(i,:) = x(i,:);
            fitness(i) = fnew;
                        
        end
        
        % Atualiza o gbest
        
        if fnew > fmax
        
            gbest = x(i,:);
            fmax = fnew;            
            
        end      

    end   
    
end

plot(conv,'LineWidth',2)
% grid 


