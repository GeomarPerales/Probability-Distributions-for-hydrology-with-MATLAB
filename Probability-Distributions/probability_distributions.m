clear
close all
clc
tic
%TABLA F(Z)
%coeficientes
ao=2.490895;
a1=1.466003;
a2=-0.024343;
a3=0.178257;

%tabla de ordenadas
%horizontal
th=0.00:0.01:0.09;

%vertical
tv=0.0:0.1:3.9;

%matriz zeros
Fz=zeros(length(tv),length(th));

%funcion
%Fz=(ao+a1*t^2+a2*t^4+a3*t^6)^(-1);
for i=1:length(tv)
    for j=1:length(th)
        Fz(i,j)=(ao+a1*(tv(i)+th(j))^2+a2*(tv(i)+th(j))^4+a3*(tv(i)+th(j))^6)^(-1);
    end
end
%-------------------------------------------------------------------------
%ABRAMOWITSY
%CALCULO DE T
%coeficientes
as_a3=0.33267;

%matriz zeros
t_as=zeros(length(tv),length(th));

%funcion
for i=1:length(tv)
    for j=1:length(th)
        t_as(i,j)=1/(1+as_a3*(tv(i)+th(j)));
    end
end

%CALCULO DE FZ
%Coeficientes
as_ao=0.43618;
as_a1=0.12017;
as_a2=0.9373;

%nueva vertical
tv1=0.0:0.1:3.1;

Fz_as=zeros(length(tv1),length(th));

%Fz=(ao+a1*t^2+a2*t^4+a3*t^6)^(-1);
for i=1:length(tv1)
    for j=1:length(th)
        Fz_as(i,j)=1-Fz(i,j).*(as_ao.*t_as(i,j)-as_a1.*t_as(i,j)^2+as_a2.*t_as(i,j)^3);
    end
end
%--------------------------------------------------------------
%MASTING
%CALCULO DE T
%coeficientes
ma_bo=0.232164;

%matriz zeros
ma_as=zeros(length(tv),length(th));

%funcion
for i=1:length(tv)
    for j=1:length(th)
        t_ma(i,j)=1/(1+ma_bo*(tv(i)+th(j)));
    end
end

%CALCULO DE FZ
%Coeficientes
m_b1=0.31938;
m_b2=-0.35656;
m_b3=1.78148;
m_b4=-1.82126;
m_b5=1.33027;

Fz_ma=zeros(length(tv1),length(th));

%Fz=(ao+a1*t^2+a2*t^4+a3*t^6)^(-1);
for i=1:length(tv1)
    for j=1:length(th)
        P(i,j)=m_b1*t_ma(i,j)+m_b2*t_ma(i,j)^2+m_b3*t_ma(i,j)^3+m_b4*t_ma(i,j)^4+m_b5*t_ma(i,j)^5;
        Fz_ma(i,j)=1-Fz(i,j).*P(i,j);
    end
end
%-------------------------------------------------------------------------
%FACTOR DE FRECUENCIA usando WEIBULL
%------------------------------------------------------------------------
%muestra 1:59
N1=[1:59];

%coeficientes
co=2.515517;
c1=0.802853;
c2=0.010328;

d1=1.432788;
d2=0.189269;
d3=0.001308;

for i=1:length(N1)
    Px1(i)=i/(length(N1)+1);
    Px11(i)=i/(length(N1)+1);
    TM1(i)=1/(1-Px1(i));
    if Px1(i)>=.5
        Px1(i)=1-Px1(i);
    end
    W1(i)=sqrt(log(1./Px1(i)^2));
    k1(i)=W1(i)-(co+c1*W1(i)+c2*W1(i)^2)/(1+d1*W1(i)+d2*W1(i)^2+d3*W1(i)^3);
end
FF59=[TM1' N1' Px1' W1' k1'];
%------------------------------------------------------------------------
%muestra 1-99
N2=[1:99];

for i=1:length(N2)
    Px2(i)=i/(length(N2)+1);
    Px22(i)=i/(length(N2)+1);
    TM2(i)=1/(1-Px2(i));
    if Px2(i)>=.5
        Px2(i)=1-Px2(i);
    end
    W2(i)=sqrt(log(1./Px2(i)^2));
    k2(i)=W2(i)-(co+c1*W2(i)+c2*W2(i)^2)/(1+d1*W2(i)+d2*W2(i)^2+d3*W2(i)^3);
end
FF99=[TM2' N2' Px2' W2' k2'];

%--------------------------------------------------------------------------
%muestra 1-199
N3=[1:199];

for i=1:length(N3)
    Px3(i)=i/(length(N3)+1);
    Px33(i)=i/(length(N3)+1);
    TM3(i)=1/(1-Px3(i));
    if Px3(i)>=.5
        Px3(i)=1-Px3(i);
    end
    W3(i)=sqrt(log(1./Px3(i)^2));
    k3(i)=W3(i)-(co+c1*W3(i)+c2*W3(i)^2)/(1+d1*W3(i)+d2*W3(i)^2+d3*W3(i)^3);
end
FF199=[TM3' N3' Px3' W3' k3'];
%--------------------------------------------------------------------------
%muestra 1-499
N4=[1:499];

for i=1:length(N4)
    Px4(i)=i/(length(N4)+1);
    Px44(i)=i/(length(N4)+1);
    TM4(i)=1/(1-Px4(i));
    if Px4(i)>=.5
        Px4(i)=1-Px4(i);
    end
    W4(i)=sqrt(log(1./Px4(i)^2));
    k4(i)=W4(i)-(co+c1*W4(i)+c2*W4(i)^2)/(1+d1*W4(i)+d2*W4(i)^2+d3*W4(i)^3);
end
FF499=[TM4' N4' Px4' W4' k4'];
%--------------------------------------------------------------------------
%muestra 1-999
N5=[1:999];

for i=1:length(N5)
    Px5(i)=i/(length(N5)+1);
    Px55(i)=i/(length(N5)+1);
    TM5(i)=1/(1-Px5(i));
    if Px5(i)>=.5
        Px5(i)=1-Px5(i);
    end
    W5(i)=sqrt(log(1./Px5(i)^2));
    k5(i)=W5(i)-(co+c1*W5(i)+c2*W5(i)^2)/(1+d1*W5(i)+d2*W5(i)^2+d3*W5(i)^3);
end
FF999=[TM5' N5' Px5' W5' k5'];
%--------------------------------------------------------------------------
%DISTRIBUCIÓN NORMAL
%información de años y caudales
file='data.xlsx';

a=xlsread(file,'a1:a59');
Q=xlsread(file,'b1:b59');

xm=mean(Q); %media de caudales
d_st=std(Q); %desv. standar de caudales
cv=d_st/xm; %coef. de variación

V_N=[xm d_st cv];

%para 2,5,20,50, 100, 200, 500 y 1000 años
TR=[2 5 10 20 50 100 200 500 1000]; %tiempo de retorno
Px=[Px22(50) Px22(80) Px22(90) Px22(95) Px22(98) Px22(99) Px33(199) Px44(499) Px55(999)];
VAE=[FF99(50,5) FF99(80,5) FF99(90,5) FF99(95,5) FF99(98,5) FF99(99,5) FF199(199,5) FF499(499,5) FF999(999,5)]; %variable aleatoria standard

for i=1:length(VAE)
    Q_n(i)=xm+VAE(i)*d_st; %caudales segun Dist. NORMAL
end

%-------------------------------------------------------------------------
%TABLA DE FACTOR DE FRECUENCIA para prob. 50%, 80%, 90% y 100% y
%coeficiente de variación (CV) que varian de 0.05 hasta 1.00.
cvn=[0.05:0.05:1];

k_cv=zeros(length(cvn),length(VAE));

for i=1:length(cvn)
    for j=1:length(VAE)
        k_cv(i,j)=(exp((log(1+cvn(i)^2))^0.5*VAE(j)-0.5*(log(1+cvn(i)^2)))-1)/cvn(i);
    end
end

%--------------------------------------------------------------------------
%LOG NORMAL 2P

lnQ=log(Q); % log natural de caudales
cvQ=d_st/xm; %coef. de variación

uy1=0.5*log(xm^2/(1+cvQ^2));
dy1=sqrt(log(1+cvQ^2));

V_LN2P=[cvQ uy1 dy1];

%usando VAE t
for i=1:length(VAE)
    Q_ln2t(i)=exp(uy1+VAE(i)*dy1); %caudales segun Dist. LOG PEARSON 2P
end

%usando k
for i=1:length(VAE)
    cv2p(i)=(exp((log(1+cvQ^2))^0.5*VAE(i)-0.5*(log(1+cvQ^2)))-1)/cvQ;
end

for i=1:length(cv2p)
    Q_ln2k(i)=xm+cv2p(i)*d_st; %caudales segun Dist. LOG PEARSON 2P
end

%-----------------------------------------------------------------------
%LOG NORMAL 3P
sg=0;

%coeficiente asimetria
for i=1:length(Q)
    sg=sg+((Q(i)-xm).^3)/length(Q);  
end
g=length(Q)^2*sg/(length(Q)-1)/(length(Q)-2)/d_st^3;
cs=g;

W=(-g+sqrt(g^2+4))*0.5;
Z2=(1-W^(2/3))/W^(1/3);
dy2=(log(Z2^2+1))^(1/2);
uy2=log(d_st/Z2)-0.5*log(Z2^2+1);
xo=xm-d_st/Z2;

V_LN3P=[cs W Z2 dy2 uy2 xo];

%usando t - hallamos Q
for i=1:length(VAE)
    Q_ln3t(i)=xo+exp(uy2+VAE(i)*dy2);
end

%usando k - hallamos Q
for i=1:length(VAE)
    k_cv3(i)=(exp((log(1+Z2^2))^0.5*VAE(i)-0.5*(log(1+Z2^2)))-1)/Z2;
end

for i=1:length(VAE)
    Q_ln3k(i)=xm+k_cv3(i)*d_st;
end

%----------------------------------------------------------------
%TABLA DE FACTOR DE FRECUENCIA PARA GUMBEL
%PARA PROB. DE 50%,80%,90%, 95% Y 100 CON INC. DE 5 U

TM=[10:5:100];

A=zeros(length(TM),100);
n=length(TM);

%generacion de las series
for i=1:length(TM)
        n(i)=length(1:TM(i));
    for j=1:TM(i)
        A(i,j)=-log(-log((n(i)+1-j)/(n(i)+1)));
    end
end

%media y desviación standard de cada serie
for i=1:length(TM)
        M(i)=mean(A(i,1:TM(i)));
        DS(i)=std(A(i,1:TM(i)),1);
end

tabFF=[TM' M' DS']; %tabla media y desv. standard

for i=1:length(TR)
    f_TR(i)=-log(-log((TR(i)-1)/TR(i)));
end

Yt=zeros(length(TM),length(f_TR));

for i=1:length(TM)
    for j=1:length(f_TR)
        Yt(i,j)=-(M(i)-f_TR(j))/DS(i);
    end
end
%------------------------------------------------------------
%DISTRIBUCION GUMBEL
al=1.2825/d_st; %alfa
u=xm-0.45*d_st; %u

%mediante t
for i=1:length(Px)
    Ytg(i)=-log(-log(Px(i)));
    Q_gt(i)=Ytg(i)/al+u;
end

%mediante k
 for i=1:length(Q)
        A59(i)=-log(-log((length(Q)+1-i)/(length(Q)+1)));
 end
 
xmg=mean(A59);
d_stg=std(A59,1);

VG=[xmg d_stg];

for i=1:length(Ytg)
        k_g(i)=(Ytg(i)-xmg)/d_stg;
        Q_gk(i)=xm+k_g(i)*d_st;
end
%----------------------------------------------------------------------
%TABLA DE FRECUENCIA PEARSON
%prob .... y coe;f de sesgo de 0.0 a 2.0 con inc de 0.1
csp=[0:0.1:2];
k_p=zeros(length(csp),length(VAE));

for i=1:length(csp)
    for j=1:length(VAE)
        gcp(i)=csp(i)/6;
        P1(i,j)=VAE(j)+(VAE(j)^2-1)*gcp(i)+(VAE(j)^3-6*VAE(j))*gcp(i)^2/3;
        P2(i,j)=-(VAE(j)^2-1)*gcp(i)^3+VAE(j)*gcp(i)^4+gcp(i)^5/3;
        k_p(i,j)=P1(i,j)-P2(i,j);
    end
end
%---------------------------------------------------------------------
%DISTRIBUCION PEARSON
%variables a usar
%media xm
%desv. std d_std
%coef. variacion cv
%coef. simetria g
N=length(Q);
gc=cs/(sqrt(N*(N-1))/(N-2)*(1+8.5/N));
be=(2/gc)^2; %beta 
alf=d_st/sqrt(be); %alfa
y=xm-d_st*sqrt(be);

VP=[gc be alf y];

%mediante t
for i=1:length(Px)
    Q_pt(i)=alf*be*(1-1/9/be+VAE(i)*sqrt(1/9/be))^3+y;
end

%mediante k

for j=1:length(VAE)
    gc1=gc/6;
    P11(j)=VAE(j)+(VAE(j)^2-1)*gc1+(VAE(j)^3-6*VAE(j))*gc1^2/3;
    P22(j)=-(VAE(j)^2-1)*gc1^3+VAE(j)*gc1^4+gc1^5/3;
    kp1(j)=P11(j)+P22(j);
end

for i=1:length(VAE)
    Q_pk(i)=xm+kp1(i)*d_st;
end
%------------------------------------------------------------------------
%DISTRIBUCIÓN LOG PEARSON
xm_lp=mean(lnQ);
ds_lp=std(lnQ);
cv_lp=ds_lp/xm_lp;

V_LP=[xm ds_lp cv_lp];

%coeficiente asimetria
sg_lp=0;
for i=1:length(Q)
    sg_lp=sg_lp+((log(Q(i))-xm_lp).^3)/length(Q);  
    endg_lp=length(Q)^2*sg_lp/(length(Q)-1)/(length(Q)-2)/ds_lp^3;
end
cs_lp=sg_lp;

gc_lp=cs_lp/(sqrt(N*(N-1))/(N-2)*(1+8.5/N));
be_lp=(2/gc_lp)^2; %beta 
sc=ds_lp*sqrt(N/(N-1));
al_lp=sc/sqrt(be_lp); %alfa
y_lp=xm_lp-al_lp*be_lp;

V_LP2=[cs_lp gc_lp be_lp sc al_lp y_lp];

%variable t
for i=1:length(VAE)
    Q_lpt(i)=exp(al_lp*be_lp*(1-1/9/be_lp+VAE(i)*(1/9/be_lp)^0.5)^3+y_lp);
end

%variable k
for j=1:length(VAE)
    gc1=gc_lp/6;
    P11(j)=VAE(j)+(VAE(j)^2-1)*gc1+(VAE(j)^3-6*VAE(j))*gc1^2/3;
    P22(j)=-(VAE(j)^2-1)*gc1^3+VAE(j)*gc1^4+gc1^5/3;
    kp2(j)=P11(j)+P22(j);
end

for i=1:length(VAE)
    Q_lpk(i)=exp(xm_lp+kp2(i)*ds_lp);
end

%exportando archivos generados a un doc
new_file='D:/hidrologia_estadistica.xlsx'; %archivo de destino

%FZ - tabla
sheet=1;
xlrange='A3';
A1=xlswrite(new_file,a,sheet,xlrange);

xlrange='B3';
A1=xlswrite(new_file,Q,sheet,xlrange);

%tabla abramowitsy
sheet=2;
xlrange='B4';
A1=xlswrite(new_file,t_as,sheet,xlrange);

xlrange='M4';
A1=xlswrite(new_file,Fz_as,sheet,xlrange);

%tabla masting
sheet=3;
xlrange='B4';
A1=xlswrite(new_file,t_ma,sheet,xlrange);

xlrange='M4';
A1=xlswrite(new_file,Fz_ma,sheet,xlrange);

%tabla frecuencia de tablas de weibull
sheet=4;
xlrange='B4';
A1=xlswrite(new_file,FF59,sheet,xlrange);

xlrange='H4';
A1=xlswrite(new_file,FF99,sheet,xlrange);

xlrange='N4';
A1=xlswrite(new_file,FF199,sheet,xlrange);

xlrange='T4';
A1=xlswrite(new_file,FF499,sheet,xlrange);

xlrange='Z4';
A1=xlswrite(new_file,FF999,sheet,xlrange);

%distribucion normal
sheet=5;
xlrange='D4';
A1=xlswrite(new_file,V_N,sheet,xlrange);

xlrange='D6';
A1=xlswrite(new_file,TR,sheet,xlrange);

xlrange='D8';
A1=xlswrite(new_file,VAE,sheet,xlrange);

xlrange='D10';
A1=xlswrite(new_file,Q_n,sheet,xlrange);

%tabla de frecuencias para de log normal
sheet=6;
xlrange='D4';
A1=xlswrite(new_file,k_cv,sheet,xlrange);

%DISTRIBUCION LOG NORMAL 2P
sheet=7;
xlrange='D4';
A1=xlswrite(new_file,V_LN2P,sheet,xlrange);

xlrange='D8';
A1=xlswrite(new_file,TR,sheet,xlrange);

xlrange='D10';
A1=xlswrite(new_file,VAE,sheet,xlrange);

xlrange='D12';
A1=xlswrite(new_file,Q_ln2t,sheet,xlrange);

xlrange='D14';
A1=xlswrite(new_file,cv2p,sheet,xlrange);

xlrange='D16';
A1=xlswrite(new_file,Q_ln2k,sheet,xlrange);

%DISTRIBUCION LOG NORMAL 3P
sheet=8;
xlrange='D4';
A1=xlswrite(new_file,V_LN3P,sheet,xlrange);

xlrange='D8';
A1=xlswrite(new_file,TR,sheet,xlrange);

xlrange='D10';
A1=xlswrite(new_file,VAE,sheet,xlrange);

xlrange='D12';
A1=xlswrite(new_file,Q_ln3t,sheet,xlrange);

xlrange='D14';
A1=xlswrite(new_file,k_cv3,sheet,xlrange);

xlrange='D16';
A1=xlswrite(new_file,Q_ln3k,sheet,xlrange);

%tabla de factor de frecuencia GUMBEL
sheet=9;
xlrange='B4';
A1=xlswrite(new_file,Yt,sheet,xlrange);

xlrange='B25';
A1=xlswrite(new_file,A,sheet,xlrange);

xlrange='L4';
A1=xlswrite(new_file,tabFF,sheet,xlrange);

%DISTRIBUCION GUMBEL
sheet=10;
xlrange='D4';
A1=xlswrite(new_file,VG,sheet,xlrange);

xlrange='D6';
A1=xlswrite(new_file,Px,sheet,xlrange);

xlrange='D8';
A1=xlswrite(new_file,TR,sheet,xlrange);

xlrange='D10';
A1=xlswrite(new_file,Ytg,sheet,xlrange);

xlrange='D12';
A1=xlswrite(new_file,Q_gt,sheet,xlrange);

xlrange='D14';
A1=xlswrite(new_file,k_g,sheet,xlrange);

xlrange='D16';
A1=xlswrite(new_file,Q_gk,sheet,xlrange);

%tabla de frecuencia PEARSON
sheet=11;
xlrange='B4';
A1=xlswrite(new_file,k_p,sheet,xlrange);

%DISTRIBUCION PEARSON
sheet=12;
xlrange='D4';
A1=xlswrite(new_file,VP,sheet,xlrange);

xlrange='D6';
A1=xlswrite(new_file,TR,sheet,xlrange);

xlrange='D8';
A1=xlswrite(new_file,VAE,sheet,xlrange);

xlrange='D10';
A1=xlswrite(new_file,Q_pt,sheet,xlrange);

xlrange='D12';
A1=xlswrite(new_file,kp1,sheet,xlrange);

xlrange='D14';
A1=xlswrite(new_file,Q_pk,sheet,xlrange);

%DISTRIBUCION LOG PEARSON
sheet=13;
xlrange='A2';
A1=xlswrite(new_file,Q,sheet,xlrange);

xlrange='B2';
A1=xlswrite(new_file,lnQ,sheet,xlrange);

xlrange='D4';
A1=xlswrite(new_file,V_LP,sheet,xlrange);

xlrange='D6';
A1=xlswrite(new_file,V_LP2,sheet,xlrange);

xlrange='D8';
A1=xlswrite(new_file,TR,sheet,xlrange);

xlrange='D10';
A1=xlswrite(new_file,VAE,sheet,xlrange);

xlrange='D12';
A1=xlswrite(new_file,Q_lpt,sheet,xlrange);

xlrange='D14';
A1=xlswrite(new_file,Q_lpk,sheet,xlrange);

toc