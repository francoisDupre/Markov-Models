// TP2 Modale
// Cressot Loic - Dupré Loic
clear

// simule une loi de géométrique de paramètre mu
function [res] = loi_expo(mu)
    res = -log(rand())/mu;
endfunction

// simule une loi uniforme sur [0;2/mu]
function [res] = loi_unif(mu)
    res = 2*rand()/mu;
endfunction



// QUESTION 1 :

// génère l'évolution d'une file d'attente avec K services jusqu'au temps T pour une loi donnée (exponetielle ou uniforme)
function [t,Xt] = file_attente(T,lambda,mu,K,loi)

    t = zeros(1);
    Xt = zeros(1);

        
    while t($)<T;
        
        p= min(Xt($),K)*mu / ( lambda + min(Xt($),K)*mu )
    
        if Xt($)==0 then
           Xt($+1) = Xt($) + 1
           t($+1) = t($) + loi(mu)
        else
            if rand() <= p then
               Xt($+1) = Xt($) - 1
            else
               Xt($+1) = Xt($) + 1
            end
            t($+1) = t($) + loi( lambda + min(Xt($-1),K)*mu )
        end
    end

endfunction

for lambda = [0.9,1,1.2];
    mu = 1;
    rho = lambda/mu;
    
    [t,Xt] = file_attente(10000,lambda,mu,1,loi_expo);
    
    clf
    plot2d2(t,Xt)
    legend(sprintf('Evolution de la file pour K=1 et une loi exponentielle et rho = %f',rho))

end


// QUESTION 2 :
// ============

// methode 1 : theoreme ergodique pour l'esperance et la variance empirique
// 
// on va regarder toutes les valeurs de rho dans 0.1 .. 0.9

precision_esp_ergo = zeros(1,9);
precision_var_ergo = zeros(1,9);

for i=1:9
    lambda = i/10;
    mu = 1;
    rho = lambda/mu;
    
    // esperance et variance theoriques
    esp_th = rho/(1-rho);
    var_th = rho/((1-rho)^2);
    

    [t,Xt] = file_attente(10000,lambda,mu,1,loi_expo);
    
    dt = t(2:$) - t(1:$-1);
    esp_emp_ergo = sum(dt.*Xt(2:$))/sum(dt)
    var_emp_ergo = sum(dt.*(( Xt(2:$) - esp_emp_ergo ).^2) )/sum(dt)
    
    precision_esp_ergo(i) = abs(esp_th - esp_emp_ergo)/esp_th
    precision_var_ergo(i) = abs(var_th - var_emp_ergo)/var_th

end

plot(precision_esp_ergo,'blue')
plot(precision_var_ergo,'red')
legend(['precision de la simulation l esperance par ergodicite','precision de la simulation la variance par ergodicite'])

// methode 2 : N simulation

lambda = 0.8
mu = 1;
rho = lambda/mu;

// esperance et variance theoriques
esp_th = rho/(1-rho);
var_th = rho/((1-rho)^2);

// On lance N expériences à jusqu'à t = T et on calcule les variances et moyennes empiriques et théoriques
N = 100;
X_T = zeros(N);
T= 1000

for i=1:N;
    [ti,X_Ti]= file_attente(T,lambda, mu);
    X_T(i)=X_Ti($);
end

// esperance et variance empiriques
esp_emp_sim = sum(X_T)/N;
var_emp_sim = sum( (X_T-esp_emp_sim).^2 )/N; // marche pas ??

precision_esp_sim = abs(esp_th - esp_emp_sim)/esp_th
precision_var_sim = abs(var_th - var_emp_sim)/var_th





// QUESTION 3 :
// ============
// on utilise le théorème ergodique avec des fonctions indicatrices pour calculer la distribution de Xt
lambda = 0.9
mu = 1;
rho = lambda/mu;

[t,Xt] = file_attente(10000,lambda,mu,1,loi_expo);
    
Xt_sta = Xt(1000:$); // regime sattionnaire
t_sta = t(1000:$);
dt_sta = t_sta(2:$) - t_sta(1:$-1);

maximum = max(Xt_sta);
distrib_Xt = zeros(maximum+1);
sum_dt_sta = sum(dt_sta)

for i=0:maximum;
    distrib_Xt(i+1) = sum(dt_sta.*  (Xt_sta(2:$)==i)  )/sum_dt_sta;
end

// affichage de la distribution
plot(0:maximum,distrib_Xt,'red')

// comparaison avec Pi
Pi = (rho*(ones(1,maximum+1))).^(0:maximum) * (1-rho)
plot(0:maximum,Pi,'blue')

legend(['distribution de Xt','Pi'])



// QUESTION 4 :
// ============
// On sait que la proba invariant de X(n) est aussi Pi (cf poly 9.2.9)
// On s'attend donc à avoir la même courbe de question 3




// QUESTION subsidiaire 1 :
// ========================
// on reprend notre fonction de simulation avec la loi uniforme cette fois (on peut montrer facilement que la simulation est bonne)

lambda = 0.9;
mu = 1;
rho = lambda/mu;

[t,Xt] = file_attente(10000,lambda, mu,1, loi_unif);

clf
plot2d2(t,Xt)
legend(sprintf('Evolution de la file pour K=1 et une loi uniforme et rho = %f',rho))




// QUESTION subsidiaire 2 :
// ========================
// on reprend notre fonction de simulation avec K > 1 cette fois

lambda = 0.9;
mu = 0.1;
K = 10;   

rho = lambda/(K*mu);


[t,Xt] = file_attente(10000,lambda, mu,K, loi_expo);

clf
plot2d2(t,Xt)
legend(sprintf('Evolution de la file pour K=10 et une loi exponentielle et rho = %f',rho))



