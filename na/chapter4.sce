//第4章のプログラム
//N. Saito 2016.07.24

//Euler method
function [t, u]=Euler(odefunc, t0, T, init, N)
//時間刻み幅の設定
// N <1のときはNを刻み幅とする（h=Nとする）
if N<1 then
    h=N;
    else 
    h = (T-t0)/N; 
end
//初期時刻と初期値の設定
t=[t0]; w = init;
//解ベクトルの定義
u=[]; u=[u;w];
//反復計算
time = t0; i=0;
while abs(time) < abs(T) // T<t<t0の場合にも対応できるように
    w = w + h*odefunc(time,w);    
    i=i+1; time = t0 + i*h
    u = [u; w]; t = [t;time];
end
endfunction

//右辺の関数1 f(t,y)
function [f] = func1(t,y)
    //f = cos(2*y);
    f = cos(y);
endfunction

//func1に対する厳密解(初期値 a=0)                                                 
function [f] = func1sol(t)                               
    f = asin((exp(2*t)-1) ./ (exp(2*t)+1))
endfunction

//Euler methodの誤差上界のチェック
// func1, func1solを使う (T=1)
//この場合誤差上界：(e^{lam T}-1)/lam * B_2 h < (e-1)0.5h < 0.86 h
function [res]=Euler2(kmax)
NN=4; hv=[]; error=[]; Uv=[]; rate=[];pv=[];T=1;t0=0;
//
for k=1:kmax
  [t, u] = Euler(func1, 0, 1, 0, NN);
  err = norm(func1sol(t)-u,%inf); error = [error; err];
  h=(T-t0)/NN; hv=[hv; h]; Uv=[Uv;0.86*h]; rate=[rate;0.86*h/err];
  NN = 2*NN;
  //pの計算
  if k>=2
            p = (log(error($-1))-log(error($)))/log(2);
        else 
            p=1;
        end
        pv=[pv;p];
end
//
res = [hv,error,Uv,rate,pv];
endfunction

//Heun method
function [t, u]=Heun(odefunc, t0, T, init, N)
//時間刻み幅の設定
// N <1のときはNを刻み幅とする（h=Nとする）
if N<1 then
    h=N;
    else 
    h = (T-t0)/N; 
end
//初期時刻と初期値の設定
t=[t0]; w = init;
//解ベクトルの定義
u=[]; u=[u;w];
//反復計算
time = t0; i=0;
while abs(time) < abs(T) // T<t<t0の場合にも対応できるように
    k1 = odefunc(time,w); k2 = odefunc(time+h,w+h*k1) 
    w = w + h*(k1+k2)/2;     
    i=i+1; time = t0 + i*h
    u = [u; w]; t = [t;time];
end
endfunction

//Runge-Kutta method
function [t, u]=RK(odefunc, t0, T, init, N)
//時間刻み幅の設定
// N <1のときはNを刻み幅とする（h=Nとする）
if N<1 then
    h=N;
    else 
    h = (T-t0)/N; 
end
//初期時刻と初期値の設定
t=[t0]; w = init;
//解ベクトルの定義
u=[]; u=[u;w];
//反復計算
time = t0; i=0;
while abs(time) < abs(T) // T<t<t0の場合にも対応できるように
    k1 = odefunc(time, w);
    k2 = odefunc(time+h/2, w+(h/2)*k1); 
    k3 = odefunc(time+h/2, w+(h/2)*k2);
    k4 = odefunc(time+h, w+h*k3);
    w = w + h*(k1 + 2*k2 + 2*k3 + k4)/6;   
    i=i+1; time = t0 + i*h
    u = [u; w]; t = [t;time];
end
endfunction


//精度(pの値)の検証；対数グラフの表示（Euler, Heun, Runge-Kutta） 
function [res] = OrderPlot2(odefunc, exact, t0, T, init)
NN=4; kmax=8; hv=[]; 
error1=[]; error2=[]; error3=[]; 
//
for k=1:kmax
  [t, u1] = Euler(odefunc, t0, T, init, NN);
  [t, u2] = Heun(odefunc, t0, T, init, NN);
  [t, u3] = RK(odefunc, t0, T, init, NN);
  err = exact(t)-u1; error1 = [error1; norm(err,%inf)];
  err = exact(t)-u2; error2 = [error2; norm(err,%inf)];
  err = exact(t)-u3; error3 = [error3; norm(err,%inf)];
  hv=[hv; (T-t0)/NN]; 
  NN = 2*NN;
end
//傾き確認用のデータ
//h1=[10^(-2),10^(-1)]; y=[10^(-6),10^(-4)];
//結果の両log目盛り表示
xset('thickness',2); 
//plot(hv, error1,"+-c"); 
//傾き確認用のデータ
h1=[10^(-2),10^(-1)]; y1=[2*10^(-4),2*10^(-3)];
h2=[10^(-2),10^(-1)]; y2=[10^(-6),10^(-4)];
h3=[10^(-2),10^(-1)]; y3=[10^(-12),10^(-8)];
//結果の両log目盛り表示
xset('thickness',2); 
plot(hv, error1,"+-c"); 
plot(hv, error2,"o-c"); 
plot(hv, error3,"v-c"); 
plot(h1, y1,"-k"); 
plot(h2, y2,"--k"); 
plot(h3, y3,"-.k"); 
//対数目盛り
a = gca();
a.log_flags = "ll";
xset('thickness',1); xgrid(color(128,128,128));  
xlabel('log(h)'); ylabel('log(error)');
legend('オイラー', 'ホイン','ルンゲ・クッタ','傾き1','傾き2','傾き4',4);
//
b=get("current_axes");
t=b.x_label; t.font_size=4;
t=b.y_label; t.font_size=4;
//結果のまとめ
res=[hv,error1,error2,error3];
//ファイル出力
xs2pdf(0,"ode_plot2.pdf");
endfunction

//Lotka-Volterraの競合系（ベクトル値２次元）
function [fn] = LotVol(t, y)
    c1=1; c2=1; d1=1; d2=1; b1=0; b2=0;
    [n,m] = size(y); fn = zeros(n, m);
    fn(1) = c1*y(1)*(1-b1*y(1)-d2*y(2));
    fn(2) = -c2*y(2)*(1-b2*y(2)-d1*y(1));
endfunction

//結果のプロット（ベクトル値２次元の場合のみ）
function OdePlot3(Solver, odefunc, t0, T, init, N)
[t, u] = Solver(odefunc, t0, T, init, N);
xset('thickness',2); 
plot(t,u(:,1),"-c"); plot(t,u(:,2),"-k");
xset('thickness',1);  
xlabel('t'); ylabel('u,v');
//軸の文字の大きさ
a=get("current_axes");
a.isoview = "on"; 
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//legend('u', 'v',"in_upper_left");
//ファイル出力
xs2pdf(0,"ode_plot3.pdf");
endfunction

//結果の相図プロット（ベクトル値２次元の場合のみ）
function OdePlot2(Solver, odefunc, t0, T, init, N)
[t, u] = Solver(odefunc, t0, T, init, N);
xset('thickness',2); 
plot(u(:,1),u(:,2),"-c"); xlabel('u'); ylabel('v');
xset('thickness',1); 
//軸の文字の大きさ
a=get("current_axes");
a.isoview = "on"; 
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//ファイル出力
xs2pdf(0,"phase0.pdf");
endfunction

//ファン・デル・ポールの微分方程式の相図
function PlotvdPol(Solver, T)
xset('thickness',2); h=0.01;
[t, u] = Solver(vdPol, 0, T, [0,0.2], h); plot(u(:,1),u(:,2),"-c"); 
[t, u] = Solver(vdPol, 0, T, [-4,-3], h); plot(u(:,1),u(:,2),"-c"); 
[t, u] = Solver(vdPol, 0, T, [4,3], h); plot(u(:,1),u(:,2),"-c"); 
[t, u] = Solver(vdPol, 0, T, [1,2], h); plot(u(:,1),u(:,2),"-c"); 
[t, u] = Solver(vdPol, 0, T, [-1,-2], h); plot(u(:,1),u(:,2),"-c"); 
xset('thickness',1);  
xlabel('u'); ylabel('v');
//軸の文字の大きさ
a=get("current_axes");
a.isoview = "on"; 
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//legend('u', 'v',"in_upper_left");
//ファイル出力
xs2pdf(0,"vdPol0.pdf");
endfunction

//ファン・デル・ポールの微分方程式（ベクトル値２次元）
function [fn] = vdPol(t, y)
    omega=1; mu=1; 
    [n,m] = size(y); fn = zeros(n, m);
    fn(1) = y(2);
    fn(2) = -omega^2*y(1) + mu*(1-y(1)^2)*y(2);
endfunction

//線形方程式
function [fn] = Linear(A,y)
    fn(1,1) = A(1,1)*y(1)+A(1,2)*y(2);
    fn(1,2) = A(2,1)*y(1)+A(2,2)*y(2);
endfunction

//Runge-Kutta method for 線形方程式
function [t, u]=RK2(A, t0, T, init, h, L)
//初期時刻と初期値の設定
t=[t0]; w = init;
//解ベクトルの定義
u=[]; u=[u;w];
//反復計算
time = t0; i=0;
while (abs(time) < abs(T)) & (norm(w,%inf)<L)// T<t<t0の場合にも対応できるように
    k1 = Linear(A, w);
    k2 = Linear(A, w+(h/2)*k1); 
    k3 = Linear(A, w+(h/2)*k2);
    k4 = Linear(A, w+h*k3);
    w = w + h*(k1 + 2*k2 + 2*k3 + k4)/6;   
    i=i+1; time = t0 + i*h
    u = [u; w]; t = [t;time];
end
endfunction

//結果の相図プロット（ベクトル値２次元の場合のみ）
function dynamics1(A, Num,T)
h=0.01; L=5;  xset('thickness',1); 
for k=1:Num
    init=2*L*rand(1,2)-[L,L];
    [t, u] = RK2(A, 0, T, init, h, L);
    [n m]=size(u); m=int(n/10); 
    xx=u(:,1); yy=u(:,2); x=xx(1:m:$); y=yy(1:m:$); 
    plot2d4(x,y,style=4);
    [t, u] = RK2(A, 0, -T, init, -h, L);
    [n m]=size(u); m=int(n/10); 
    xx=u(:,1); yy=u(:,2); x=xx(1:m:$); y=yy(1:m:$); 
    plot2d4(x($:-1:1),y($:-1:1),style=4);
end
xlabel('u'); ylabel('v');
xset('thickness',1);  
square(-L,-L,L,L); plot([-L,L],[0,0],"--k");  plot([0,0],[-L,L],"--k");  
//軸の文字の大きさなど
a=get("current_axes");
a.isoview = "on"; 
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//ファイル出力
xs2pdf(0,"dynamic1.pdf");
endfunction

//惑星運動の方程式
function [fn] = Planet(t,y)
    [n,m] = size(y); fn = zeros(n, m);
    R=(y(1)^2+y(2)^2)^(3/2);
    fn(1) = y(3);
    fn(2) = y(4);
    fn(3) = -y(1)/R;
    fn(4) = -y(2)/R;    
endfunction

//惑星運動の方程式の軌道
function PlotPlanet(Solver, T,k)
xset('thickness',2); h=0.01;
u0=[1-k, 0, 0, sqrt((1+k)/(1-k))];
[t, u] = Solver(Planet, 0, T, u0, h); plot(u(:,1),u(:,2),"-c"); 
xset('thickness',1);  xlabel('x'); ylabel('y');
x0=min(u(:,1));x1=max(u(:,1));y0=min(u(:,2));y1=max(u(:,2));
lx=0.05*(x1-x0); ly=0.05*(y1-y0);
square(x0-lx,y0-ly,x1+lx,y1+ly);
//軸の文字の大きさ
a=get("current_axes");
a.isoview = "on"; 
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//legend('u', 'v',"in_upper_left");
//ファイル出力
xs2pdf(0,"Planet0.pdf");
endfunction

//単振動
function [fn] = harmonic(t, y)
    c1=0; c2=-1; c3=1; c4=0; 
    [n,m] = size(y); fn = zeros(n, m);
    fn(1) = c1*y(1) + c2*y(2);
    fn(2) = c3*y(1) + c4*y(2);    
endfunction

//単振動の厳密解
function [y] = harmonic_exact(t)
    y(1) = cos(t); 
    y(2) = sin(t);    
endfunction

//結果の相図プロット
// Heun for dx/dt=-y, dy/dt=x
function OdePlot4(T, N)
[t, u] = Heun(harmonic, 0, T, [1,0], N);
plot(u(:,1),u(:,2),"-c"); xlabel('x'); ylabel('y');
square(-1.2,-1.2,1.2,1.2);
//軸の文字の大きさ
a=get("current_axes");
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//ファイル出力
xs2pdf(0,"heun1.pdf");
endfunction

//Crank-Nicolson法：dx/dt=-y, dy/dt=x
function CN1(T, h)
//時間刻み幅の設定
//初期時刻と初期値の設定
t0=0; t=[t0]; 
//解ベクトルの定義
w=[1;0];u=w;
//反復行列の定義
c=h/2; H=(1/(1+c^2))*[1-c^2,-2*c;2*c,1-c^2];
//反復計算
time = t0; i=0;
while time < T 
    w = H*w;     
    i=i+1; time = t0 + i*h
    u = [u,w]; t = [t;time];
end
//解のプロット
plot(u(1,:),u(2,:),"-c"); xlabel('x'); ylabel('y');
square(-1.2,-1.2,1.2,1.2);
a=get("current_axes");
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//ファイル出力
xs2pdf(0,"cn1.pdf");
endfunction

//Crank-Nicolson method for dx/dt=-y, dy/dt=x
//計算のみ（OrderPlot5用）
function [t,v]=CN2(T, h)
//時間刻み幅の設定
//初期時刻と初期値の設定
t0=0; t=[t0]; 
//解ベクトルの定義
w=[1;0];u=w;
//反復行列の定義
c=h/2; H=(1/(1+c^2))*[1-c^2,-2*c;2*c,1-c^2];
//反復計算
time = t0; i=0;
while time < T 
    w = H*w;     
    i=i+1; time = t0 + i*h
    u = [u,w]; t = [t;time];
end
//解としてはuを転置したものを返す（他との兼ね合いのため）
v=u'; 
endfunction

//誤差の観察
//精度(pの値)の検証；対数グラフの表示
// Heun vs. Crank-Nicolson for dx/dt=-y, dy/dt=x
function [res] = OrderPlot5(T)
kmax=8; hv=[]; error1=[]; error2=[];
t0 = 0; init=[1,0]; hh=(T-t0)/4;
//
for k=1:kmax
  [t1, u1] = Heun(harmonic, t0, T, init, hh);
  err =[cos(t1),sin(t1)]-u1; error1 = [error1; norm(err,%inf)];
  [t2, u2] = CN2(T,hh);
  err =[cos(t2),sin(t2)]-u2; error2 = [error2; norm(err,%inf)];
  hv=[hv; hh]; 
  hh = 0.5*hh;
end
//傾き確認用のデータ
h1=[10^(-2),10^(-1)]; y=[10^(-6),10^(-4)];
//結果の両log目盛り表示
//plot2d(hv, [error1, error2], style=[1,4], logflag="ll"); 
xset('thickness',2); 
plot(hv, error1,"+-c"); 
plot(hv, error2,"o-c"); 
plot(h1, y,"-k"); 
//対数目盛り
a = gca();
a.log_flags = "ll";
xset('thickness',1); xgrid(color(128,128,128));  
xlabel('log(h)'); ylabel('log(error)');
legend('Heun', 'Crank-Nicolson','傾き2',4);
//
b=get("current_axes");
t=b.x_label; t.font_size=4;
t=b.y_label; t.font_size=4;
//結果のまとめ
res=[hv,error1,error2];
//ファイル出力
xs2pdf(0,"heun_cn1.pdf");
endfunction

//ロディスティック方程式の右辺
function [f] = func8(t,y)
    p=1;q=3;
    f = (p-q*y)*y
endfunction

//ロディスティック方程式の厳密解
function [f] = func8sol(t)
    p=1;q=3;u0=3.8;
    f = p*u0*exp(p*t) ./ (p + q*u0*(exp(p*t)-1));  
endfunction

//ロジスティック方程式　Euler法の解の描画
function OdePlot8(T)
init=3.8; 
[t1, u1] = Euler(func8, 0, T, init, 0.01);
[t2, u2] = Euler(func8, 0, T, init, 0.1);
xset('thickness',2); plot(t1, u1,"xc");plot(t2, u2,"oc"); 
//厳密解の描画
xset('thickness',1); 
tt=linspace(0,T,100)'; yy=func8sol(tt); plot(tt, yy,"--k");
//軸の文字の大きさ
xlabel('t'); ylabel('u');
a=get("current_axes");
//a.isoview = "on"; 
t=a.x_label; t.font_size=3;
t=a.y_label; t.font_size=3;
legend('h=0.01', 'h=0.1','厳密解',1);
//ファイル出力
xs2pdf(0,"logistic0.pdf");
endfunction

//境界値問題の右辺関数1
function [f] = func10(x)
    a=3*%pi;
    f = 2*a*cos(a*x) + a^2*(1-x) .* sin(a*x);   
endfunction

//右辺関数1の厳密解
function [f] = func10sol(x)
    f = (1-x) .* sin(3*%pi*x);   
endfunction


//境界値問題1
function [xx,uu]=bvp1(N, L, rfunc)
a = 0.0; b = L; ua=0.0; ub=0.0;
h = (b - a)/(N + 1); x = [a + h: h: b - h]'; xx = [a; x; b];
// Aとfの定義
A = 2*eye(N, N)-diag(ones(N-1 , 1), -1)-diag(ones(N-1 , 1), 1); 
A = (1/(h^2))*A;
f=rfunc(x);
//Au=fの解法
u=A\f; 
//
uu = [ua; u; ub]; 
endfunction

//境界値問題のプロット用1
function bvpPlot1(N, L, rfunc)
[xx,uu]=bvp1(N, L, rfunc);
xset('thickness',2);  plot(xx, uu, "-c");
xset('thickness',1);  xlabel('x'); ylabel('u'); 
//軸の文字の大きさ
a=get("current_axes");
//a.isoview = "on"; 
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//ファイル出力
xs2pdf(0,"bvp1.pdf");
endfunction

//境界値問題の誤差の観察
function [res]=bvpOrder1(NN)
kmax=8; hv=[]; error1=[]; N=NN;
//
for k=1:kmax
  [xx, uu] = bvp1(N,1,func10);
  err = func10sol(xx)-uu; error1 = [error1; norm(err,%inf)];
  hh=1/(N+1); hv=[hv; hh]; N = (k+1)*NN;
end
//傾き確認用のデータ
h1=[5*10^(-2),10^(-1)]; y=[5*10^(-3),2*10^(-2)];
//結果の両log目盛り表示
xset('thickness',2); 
plot(hv, error1,"o-c"); 
plot(h1, y,"-k"); 
//対数目盛り
a = gca();
a.log_flags = "ll";
xset('thickness',1); xgrid(color(128,128,128));  
xlabel('log(h)'); ylabel('log(error)');
legend('誤差', '傾き2',4);
//
b=get("current_axes");
t=b.x_label; t.font_size=4;
t=b.y_label; t.font_size=4;
//結果のまとめ
res=[hv,error1];
//ファイル出力
xs2pdf(0,"bvporder1.pdf");
endfunction

//境界値問題2
function bvp5(N, ep)
a = 0.0; b = 1; ua=0.0; ub=0.0;
h = (b - a)/(N + 1); x = [a + h: h: b - h]'; xx = [a; x; b];
// Aとfの定義
A = 2*eye(N, N)-diag(ones(N-1 , 1), -1)-diag(ones(N-1 , 1), 1); 
A = (1/(h^2))*A;
f=func11(x,ep);
//Au=fの解法
u=A\f; 
//
uu = [ua; u; ub]; xset('thickness',2);  plot(xx, uu, "-c");
xset('thickness',1);  xlabel('x'); ylabel('u'); square(0,0,1,0.26);
//軸の文字の大きさ
a=get("current_axes");
a.isoview = "on"; 
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//ファイル出力
xs2pdf(0,"bvp5.pdf");
endfunction


//境界値問題の右辺関数2; bvp5用
function [f] = func11(x,ep)
    [n,m] = size(x); f = zeros(n, m);
     for i=1:n
           if abs(x(i)-0.5)<ep then
            f(i)=1/(2*ep);
           end
    end
endfunction

//問題4.30用
function [res]=rate(func,sol,Solver,t0,T,a,kmax)
NN=4; hv=[]; error=[]; pv=[];
//
for k=1:kmax
  [t, u] = Solver(func, t0, T, a, NN);
  err = norm(sol(t)-u,%inf); error = [error; err];
  h=(T-t0)/NN; hv=[hv; h]; NN = 2*NN;
  //pの計算
  if k>=2
            p = (log(error($-1))-log(error($)))/log(2);
        else 
            p=1;
        end
        pv=[pv;p];
end
//
res = [hv,error,pv];
endfunction

//問題4.30 右辺の関数 f(t,y)
function [f] = func5(t,y)
    g = -t ./ sqrt(1-t.^2) -(1 - t.^2);
    f = y.^2+g;
endfunction

//問題4.30 func5に対する厳密解                                               
function [f] = func5sol(t)                               
    f = sqrt(1-t.^2);
endfunction


//問題4.31用  結果の2次元プロット
function OdePlot1(Solver, odefunc, t0, T, init, N)
[t, u] = Solver(odefunc, t0, T, init, N);
xset('thickness',2); 
plot2d(t, u,style=4);xlabel('t'); ylabel('u');
xset('thickness',1); 
//軸の文字の大きさ
a=get("current_axes");
a.isoview = "on"; 
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//ファイル出力
xs2pdf(0,"ode_plot0.pdf");
endfunction

//問題4.31 右辺のf(t,y)
function [f] = func6(t,y)
    f = -10.*y+t;
endfunction

//問題4.31 func1に対する厳密解                                                 
function [f] = func6sol(t)                               
    f = 0.1 .* t - 0.01 + 10.01 .* exp(-10.*t);
endfunction


// 3rd Heun
function [t, u]=heun3(odefunc, t0, T, init, N)
//時間刻み幅の設定
// N <1のときはNを刻み幅とする（h=Nとする）
if N<1 then
    h=N;
    else 
    h = (T-t0)/N; 
end
//初期時刻と初期値の設定
t=[t0]; w = init;
//解ベクトルの定義
u=[]; u=[u;w];
//反復計算
time = t0; i=0;
while abs(time) < abs(T) // T<t<t0の場合にも対応できるように
    k1 = odefunc(time, w);
    k2 = odefunc(time+h/4, w+(h/4)*k1); 
    k3 = odefunc(time+2*h/3, w+(h/9)*(-2*k1+8*k2));
    w = w + h*(k1 + 3*k3)/4;   
    i=i+1; time = t0 + i*h
    u = [u; w]; t = [t;time];
end
endfunction

// 3rd Kutta
function [t, u]=kutta3(odefunc, t0, T, init, N)
//時間刻み幅の設定
// N <1のときはNを刻み幅とする（h=Nとする）
if N<1 then
    h=N;
    else 
    h = (T-t0)/N; 
end
//初期時刻と初期値の設定
t=[t0]; w = init;
//解ベクトルの定義
u=[]; u=[u;w];
//反復計算
time = t0; i=0;
while abs(time) < abs(T) // T<t<t0の場合にも対応できるように
    k1 = odefunc(time, w);
    k2 = odefunc(time+h/2, w+(h/2)*k1); 
    k3 = odefunc(time+h, w+h*(-k1+2*k2));
    w = w + h*(k1 + 4*k2 + k3)/6;   
    i=i+1; time = t0 + i*h
    u = [u; w]; t = [t;time];
end
endfunction


//問題4.32 用精度(pの値)の検証；対数グラフの表示
function [res] = OrderPlot7(odefunc, exact, t0, T, init)
NN=4; kmax=8; hv=[]; 
error1=[]; error2=[]; error3=[]; 
//
for k=1:kmax
  [t, u1] = heun3(odefunc, t0, T, init, NN);
  [t, u2] = kutta3(odefunc, t0, T, init, NN);
  err = exact(t)-u1; error1 = [error1; norm(err,%inf)];
  err = exact(t)-u2; error2 = [error2; norm(err,%inf)];
  hv=[hv; (T-t0)/NN]; 
  NN = 2*NN;
end
//傾き確認用のデータ
//h1=[10^(-2),10^(-1)]; y=[10^(-6),10^(-4)];
//結果の両log目盛り表示
xset('thickness',2); 
//plot(hv, error1,"+-c"); 
//傾き確認用のデータ
h1=[10^(-2),10^(-1)]; y1=[10^(-10),10^(-7)];
//結果の両log目盛り表示
xset('thickness',2); 
plot(hv, error1,"+-c"); 
plot(hv, error2,"o-c"); 
plot(h1, y1,"-k"); 
//対数目盛り
a = gca();
a.log_flags = "ll";
xset('thickness',1); xgrid(color(128,128,128));  
xlabel('log(h)'); ylabel('log(error)');
legend('3次ホイン', '3次クッタ','傾き3',4);
//
b=get("current_axes");
t=b.x_label; t.font_size=4;
t=b.y_label; t.font_size=4;
//結果のまとめ
res=[hv,error1,error2,error3];
//ファイル出力
xs2pdf(0,"ode_plot7.pdf");
endfunction
////////////////////////////////////////////////////
