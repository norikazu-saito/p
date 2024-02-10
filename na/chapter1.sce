//第１章のプログラム
//N. Saito 2015.09.01

//ニュートン法
function [res]=newton1(fun,dfun,kmax,tol,xig)
    xv=[]; fx=[]; inc=[]; it=[]; 
    x=xig; k=0; dif=tol+1; 
    while k<kmax & dif>tol
        it=[it;k]; xv=[xv;x]; f=fun(x); fx=[fx;f]; df=dfun(x);
        if abs(df)<1.0D-15 then
           error('df/dx=0');
        end
        xnew=x-f/df; dif=abs(xnew-x)/abs(x); 
        inc=[inc;dif]; x=xnew; k=k+1;
    end
//結果のまとめ
res=[it,xv,fx,inc];
endfunction

// 関数
// y=x^3-3x-1
function [y]=func11(x)
    y=x.^3-3*x-1
endfunction 
// y'=3x^2-3
function [y]=dfunc11(x)
    y=3*x.^2-3
endfunction 

//簡易ニュートン法
function [res]=simp_newton1(fun,dfun,kmax,tol,xig)
    xv=[]; fx=[]; inc=[]; it=[]; 
    x=xig; k=0; dif=tol+1; 
    df=dfun(x);
    if abs(df)<1.0D-15 then
           error('df/dx=0');
    end
    while k<kmax & dif>tol
        it=[it;k]; xv=[xv;x]; f=fun(x); fx=[fx;f]; 
        xnew=x-f/df; dif=abs(xnew-x)/abs(x); 
        inc=[inc;dif]; x=xnew; k=k+1;
    end
//結果のまとめ
res=[it,xv,fx,inc];
endfunction

//セカント法
function [res]=secant1(fun,kmax,tol,xig0,xig1)
    x=xig1; xx= xig0; ff=fun(xx);  k=1; dif=tol+1; 
    xv=[xx]; fx=[fun(xx)]; inc=[0]; it=[0]; 
    while k<kmax & dif>tol
        // x = x_k, xx = x_{k-1}, f=f(x_k), ff=f(x_{k-1})
        f=fun(x); df=(f-ff)/(x-xx);
        it=[it;k]; xv=[xv;x]; fx=[fx;f]; 
        xnew=x-f ./ df; dif=abs(xnew-x)/abs(x); 
        inc=[inc;dif]; xx=x; x=xnew; ff=f; k=k+1;
    end
//結果のまとめ
res=[it,xv,fx,inc];
endfunction

//不動点反復法
function [res]=fixed1(fun,phi,kmax,tol,xig)
    xv=[]; fx=[]; inc=[]; it=[];
    x=xig; k=0; dif=tol+1; 
    while k<kmax & dif>tol
        it=[it;k]; xv=[xv;x]; f=fun(x); fx=[fx;f];
        xnew=phi(x); dif=abs(xnew-x)/abs(x); 
        inc=[inc;dif]; x=xnew; k=k+1;
    end
//結果のまとめ
res=[it,xv,fx,inc];
endfunction
//不動点反復1
function [y]=phi1(x)
    y=x-(1/36)*(x^3-3*x-1)
endfunction 
//不動点反復2
function [y]=phi2(x)
    y=x-(x^3-3*x-1)
endfunction 

//反復法の比較
// x^3-3x-1=0 の正根1.879385241571816
function [err1,err2,err3]=comparison1(fun,dfun)
    xig =3; tol=1e-5; kmax=20; kmax0=3; res=[];
    k1=[]; k2=[]; k3=[]; k4=[];  
    err1=[]; err2=[]; err3=[]; err4=[]; 
    //ニュートン法 [res]=newton1(func11,dfunc11,20,1e-15,3)
    value=1.87938524157181663; 
    x=xig; k=0; dif=tol+1; 
    while k<kmax & dif>tol
        f=fun(x); df=dfun(x);
        if abs(df)<1.0D-15 then
           error('df/dx=0');
        end
        xnew=x-f/df; dif=abs(xnew-x)/abs(x); 
        x=xnew; k=k+1;  err1=[err1;abs(x-value)]; k1=[k1;k]; 
    end
    //簡易ニュートン法[res]=simp_newton1(func11,dfunc11,100,1e-15,3)
    value=1.87938524157182196;
    x=xig; k=0; dif=tol+1; df=dfun(x);
        if abs(df)<1.0D-15 then
           error('df/dx=0');
        end
    while k<kmax & dif>tol
        f=fun(x); xnew=x-f/df; dif=abs(xnew-x)/abs(x);
        x=xnew; k=k+1;  err2=[err2;abs(x-value)]; k2=[k2;k]; 
    end
     //不動点反復法1 [res]=fixed1(func11,phi1,500,1e-15,3)
     value=1.87938524157182463; 
     x=xig; k=0; dif=tol+1; 
    while k<kmax & dif>tol
        f=fun(x); xnew=phi1(x); dif=abs(xnew-x)/abs(x); 
        x=xnew; k=k+1; err3=[err3;abs(x-value)]; k3=[k3;k]; 
    end
//結果の表示
xset('thickness',2); 
plot(k1,log(err1),"-^c");  
plot(k2,log(err2),"-xc");  
plot(k3,log(err3),"-sc");  
//
xset('thickness',1);
xlabel("k"); ylabel("log E_k");
a=get("current_axes");
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
xs2pdf(0,"comp1.pdf"); 
endfunction

//収束の速さの検証
// 厳密な解 sol がわかっていることが前提
function [pv,err]=rate1(fun,dfun,sol,xig)
    kmax = 50; tol = 1.0D-8; //適当に設定する
    res=[]; xv=[]; it=[]; err=[]; pv=[];
    //調べたい方法のみを残し，あとの２つはコメントアウト（現状はセカント法で計算）
    //[res]=newton1(fun,dfun,kmax,tol,xig); //ニュートン法
    //[res]=simp_newton1(fun,dfun,kmax,tol,xig); //簡易ニュートン法
    xig0=xig; xig1=xig-0.1;[res]=secant1(fun,kmax,tol,xig0,xig1);//セカント法
    it=res(:,1); xv=res(:,2); // it: 反復番号,  xv: 反復列
    [N,NN]=size(it); //反復列の長さ
    err=log(abs(xv-sol*ones(N,1))); //誤差の対数
    for i=1:N-2
        pv(i)=(err(i+2)-err(i+1))/(err(i+1)-err(i));
    end    
endfunction

//多変数のニュートン法
function [xv,fx,inc,k]=newtonMulti(fun,dfun,kmax,tol,xig)
    xv=[]; fx=[]; inc=[];
    x=xig; k=0; dif=tol+1; 
    while k<kmax & dif>tol
        xv=[xv,x]; f=fun(x); fx=[fx,f]; df=dfun(x);
        if rcond(df)<1.0D-15 then
            error('df/dx is singular');
        end
        xnew=x-df\f; dif=norm(xnew-x,1); 
        inc=[inc,dif]; x=xnew; k=k+1;
    end
endfunction

//２変数のニュートン法の例
function [xv,fx,inc,k]=td_newton(xig)
//初期値の例 xig = [-0.9;0.1],[0.4;0.4],[0.9;0.1],[-0.45;-1.2]
kmax = 100;tol = 1.0D-15;
[xv1,fx1,inc1,k] = newtonMulti(Ffunc,Jfunc,kmax,tol,xig);
xv=xv1'; fx=fx1'; inc=inc1';
xset('thickness',2);
plot2d4(xv1(1,:),xv1(2,:),style=4); 
xlabel("x"); ylabel("y");
xset('thickness',1); 
//初期値
plot(xig(1,1),xig(2,1),"ok");
//近似解
plot(xv1(1,$), xv1(2,$),"xk"); 
//軸の文字の大きさ
a=get("current_axes");
a.isoview = "on"; 
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//ファイル出力
xs2pdf(0,"mnew1.pdf");
endfunction

// (1.29)の関数f(x,y), g(x,y)
function [F] = Ffunc(x)
//F(x,y) = [x^2+y^2-1 ; sin(\pix/2)+y^3]
//F(x,y) = 0 \iif (x,y) \simeq \pm (-0.476095758,0.879393409)
    F = x(1,:).^2+x(2,:).^2-1;
    F= [F;sin(%pi*x(1,:)/2)+x(2,:).^3];
endfunction
// Ffucのヤコビ行列
function [J]=Jfunc(x)
//J(x,y) = [2x , 2y ; \pi/2 sin(\pix/2),3y^2]
// detJ = 0 \iff y=0 or y = \pi/6 sin(\pix/2))/2
    J = zeros(2,2);
    pi2 = 0.5*%pi;
    J(1,1) = 2*x(1);
    J(1,2) = 2*x(2);
    J(2,1) = pi2*cos(pi2*x(1));
    J(2,2) = 3*x(2)^2;
endfunction

//連立法による代数方程式
//  例題の多項式の作り方：
//  --> zz=[1 2 -2+%i -2-%i 3+2*%i]; 
//  --> f=poly(zz,"z");
function [k,xx]=dk1(zz, beta, r, kmax)
f=poly(zz,"z"); // f0=0: 解くべき方程式
n=degree(f);   // 解の個数
//Aberthの方法による初期値の設定
// betaとrは別に計算して与えておく
gamma=1.7;
x0=beta+r*exp(%i*(2*%pi/n*(0:n-1)+gamma));
xx=[x0]; //反復列
inc=ones(1,n); //増分
//
//反復計算  
k=0; tol=1.0e-12;
while k <= kmax & norm(inc)>tol 
    k=k+1; 
    f0=horner(f,x0);
    s=[];
    for j=1:n
        q=1; 
        for i=1:n
                if i ~= j then
                    q = q*(x0(j)-x0(i));
                end
         end
         s=[s,q]
    end
    x1=x0-f0 ./ s; xx=[xx;x1];
    inc=x0-x1;
    x0=x1;
end
// 初期値の円の描画
t=linspace(0,2*%pi,200);
xi=beta+r*exp(%i*t);
 xset('thickness',2); 
plot(real(xi),imag(xi),"k--");
//近似解の表示
plot(real(xx),imag(xx),"-xc");
//枠の表示
r1=1.1*r;
square(beta-r1,beta-r1,beta+r1,beta+r1);
xset('thickness',1); 
xlabel("Real", "fontsize", 4); ylabel("Imaginary", "fontsize", 4);
//ファイル出力
xs2pdf(0,"dk1.pdf");
endfunction 

//問題2.2用
function [y]=func61(x)
    y=cosh(x).*cos(x)-1;
endfunction 
//問題2.2用
function [y]=dfunc61(x)
    y=sinh(x).*cos(x)-cosh(x).*sin(x);
endfunction 

//問題2.4用
function [y]=func71(x)
    L=160;  h=15; 
    y=(cosh(L*x/2)-1) -h*x;
endfunction 
//問題2.4用
function [y]=dfunc71(x)
    L=160;  h=15; 
    y=(L/2)*sinh(L*x/2) -h;
endfunction 

//関数のプロット　問題2.2, 2.4用 
function fplot(a,b,fun)
n=200; x=linspace(a,b,n+1)'; y=fun(x);
xset('thickness',3); 
plot(x,y,"-c");
//
xset('thickness',1); 
xx=[a;b];yy=[0;0];
plot(xx,yy,"-k");
//
xlabel("x"); ylabel("y");
//軸の文字の大きさ
a=get("current_axes");
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//ファイル出力
xs2pdf(0,"fplot1.pdf");
endfunction
////////////////////////////////////
