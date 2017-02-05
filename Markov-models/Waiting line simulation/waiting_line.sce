// TP2 Modale
// Cressot Loic - Dupré Loic
clear

// simule une loi de géométrique de paramètre mu
function [res] = loi_expo(mu)
    res = -log(rand())/mu
endfunction

// QUESTION 1 :

// génère l'évolution d'une file d'attente avec K=1 jusqu'au temps T
function [t,Xt] = file_attente(T,lambda,mu)

    t = zeros(1);
    t(1)=0;
    Xt = zeros(1);
    Xt(1) = 0;
    p=[ mu/(lambda+mu) ; lambda/(lambda+mu) ]
        
    while t($)<T;
    
        if Xt($)==0 then
           Xt($+1) = Xt($) + 1
           t($+1) = t($) + loi_expo(mu)
        else
            if rand() <= p(1) then
               Xt($+1) = Xt($) - 1
            else
               Xt($+1) = Xt($) + 1
            end
            t($+1) = t($) + loi_expo(lambda+mu)
        end
    end

endfunction

lambda = 0.8;
mu = 1;
rho = lambda/mu;

[t,Xt] = file_attente(1000,lambda, mu);

clf
plot2d2(t,Xt)
legend(sprintf('rho = %f',rho))



// QUESTION 2 :

// esperance et variance theoriques
esp_th = rho/(1-rho);
var_th = rho/(1-(rho^2));

// methode 1 : theoreme ergodique pour l'esperance et la variance empirique
dt = t(2:$) - t(1:$-1);
esp_emp_ergo = sum(dt.*Xt(2:$))/sum(dt)
var_emp_ergo = sum(dt.*(( Xt(2:$) - esp_emp_ergo ).^2) )/sum(dt) // marche pas  ??

// methode 2 : N simulation
// On lance N expériences à jusqu'à t = T et on calcule les variances et moyennes empiriques et théoriques
N = 1000;
X_T = zeros(N);
T= 100

for i=1:N;
    [ti,X_Ti]= file_attente(T,lambda, mu);
    X_T(i)=X_Ti($);
end

// esperance et variance empiriques
esp_emp_sim = sum(X_T)/N;
var_emp_sim = sum( (X_T-esp_emp_sim).^2 )/N; // marche pas ??



