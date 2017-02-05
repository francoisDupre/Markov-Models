// MODALE TP1 CRESSOT DUPRE
clear

// test de chi2
function[proba]=test_chi2(N,p0)
  n=sum(N);// taille de l'echantillon observe
  // calcul de zeta_ n
  zeta_n=n*sum(((N/n-p0).^2)./p0);
  // nombre de degres de liberte  (= nombre de classes dans N-1)
  d= length(N)-1;
  // on calcule la proba pour un chi 2 à d-1 degres d'etre superieur a zeta
  [p,q]=cdfchi("PQ",zeta_n,d);
  proba=q;
endfunction;


// Loi normale
function [x]=normale(y,m,s2)
  x=%e^(-(y-m).^2/2/s2)/sqrt(2*%pi*s2)
endfunction;

// Ouvrir le fichier de données (nombre de crabes par intervalle)
x=fscanfMat('crabe.txt');
x=x;

// intervalles
y=.580+.002+.004*[0:28];
yM=y+.002;
ym=y-.002;
Max=25;

// Dessiner la loi normale correspondante
// A FAIRE
// moyenne
x1=x.*y'
mu=sum(x1)/sum(x)
sigma2 = sum(((y-mu).^2).*x')/sum(x)
norma = normale(y,mu,sigma2)

plot2d(y,norma);

// Tracer l'histogramme
// A FAIRE
div = 0.004*sum(x)
bar(y,x/div);

// test du chi 2
X2 = test_chi2(x/div,norma')


// Données
pi0=[1; 3]/2/2;
pi=pi0;
mu=[.57; .67];
s2=[1 ;1]/10000;



// Algorithme EM pour les crabes 
//------------------------------
clf(); // clear figure

N=1000;
P=10;
rho=ones(2,N);
R=zeros(5,P);
R(:,1)=[mu(1);mu(2);pi(1);s2(1);s2(2)];
Likelihood=zeros(2,N)

// calculate all crab sizes in Y
Y=zeros(1,1000) // crabs sizes
counter=1
class=1
for k=1:1000;
    Y(k)=y(class)

    if counter*1.==x(class) then
        class = class+1
        counter=1
    end
    counter = counter+1
end


for k=1:P-1;
  // Iteration k
  // A FAIRE
  vec_normal=[ normale(Y,R(1,k),R(4,k)) ; normale(Y,R(2,k),R(5,k))] 
  vect_init=[R(3,k) ; 1-R(3,k)]'
  f_mu=zeros(2,N)
  f_mu(1,:)=vect_init(1)*vec_normal(1,:)
  f_mu(2,:)=vect_init(2)*vec_normal(2,:)
  f_sigma=f_mu(1,:) + f_mu(2,:)
  
 
  rho(1,:)=(vec_normal(1,:)*vect_init(1))./f_sigma
  rho(2,:)=(vec_normal(2,:)*vect_init(2))./f_sigma

  // mise à jour des paramètres
  R(3,k+1)= sum(rho(1,:))/N
  R(1,k+1)=(rho(1,:)*Y')/sum(rho(1,:))
  R(2,k+1)=(rho(2,:)*Y')/sum(rho(2,:))
  R(4,k+1)=( rho(1,:)*((Y-R(1,k)).^2)')/sum(rho(1,:))
  R(5,k+1)=( rho(2,:)*((Y-R(2,k)).^2)')/sum(rho(2,:))
end;

// Affichages
// A FAIRE
plot([1:P],R(1,:));
plot([1:P],R(2,:));
plot([1:P],R(3,:));
plot([1:P],R(4,:));
plot([1:P],R(5,:));
