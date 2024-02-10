//第2章のプログラム
//N. Saito 2016.02.02

//例題の関数
function y = func1(x)
    y = sin(x);
endfunction
//厳密な積分値
Q1ex=1.0-cos(2.0);

// 複合中点則
function integ = c_mid(a,b,m,func) 
    x = linspace(a,b,m+1);
    integ = 0;
    for i=1:m
        h = x(i+1)-x(i); c = (x(i) + x(i+1))/2;
        integ = integ + func(c)*h;
    end
endfunction

// 複合矩形則 R^l_n
function integ = c_rect(a,b,m,func) 
    x = linspace(a,b,m+1);
    integ = 0;
    for i=1:m
        h = x(i+1)-x(i); c = x(i);
        integ = integ + func(c)*h;
    end
endfunction

// 複合台形則
function integ = c_trape(a,b,m,func)
    x = linspace(a,b,m+1);
    integ = 0;
    for i=1:m
        h = x(i+1)-x(i);
        integ = integ + (func(x(i)) + func(x(i+1)))*h/2;
    end
endfunction

// 複合シンプソン則
function integ = c_simps(a,b,m,func)
    x = linspace(a,b,m+1);
    integ = 0;
    for i=1:m
        h = x(i+1)-x(i); 
        c = (x(i) + x(i+1))/2;
        integ = integ + (func(x(i)) + 4*func(c) + func(x(i+1)))*h/6;
    end
endfunction

//  複合公式の収束の速さ
function [res] = NI_order(a,b,func,rule,Qex,m0)
kmax=5; evec=[]; qvec=[]; hvec=[]; pvec=[]; m=m0;
//
for k=1:kmax
    hvec=[hvec;(b-a)/m]; 
    q=rule(a,b,m,func); qvec=[qvec;q]; evec=[evec;abs(Qex - q)];
    //
    if k>=2
            p = (log(evec($-1))-log(evec($)))/log(2.0);
    else 
            p=1;
    end
    pvec=[pvec;p];
    //
    m=2*m; 
end
//
res=[hvec,qvec,evec,pvec];
endfunction

// 複合公式の収束の速さ（両対数グラフ）1 //使っていない？
function order_plot(a,b,func,rule,Qex,mmin,mmax)
    h = [];
    err = [];
    for m = mmin:mmax
        tmp = abs(Qex - rule(a,b,m,func));
        err = [err; tmp];
        h = [h; (b-a)/m];
    end
xset('thickness',2);
    plot2d(h, err, style=5,logflag="ll");
xset('thickness',1);
xlabel("log h"); ylabel("log E_h");
//軸の文字の大きさ
a=get("current_axes");
//a.isoview = "on"; 
t=a.x_label; t.font_size=4;
t=a.y_label; t.font_size=4;
//ファイル出力
xs2pdf(0,"nc0.pdf");
endfunction

// 複合公式の収束の速さ（両対数グラフ）2
function order_plot1(a,b,func,rule,Qex)
h = []; err = []; kmax=7; m=4;
//
for k =1:kmax
     m=m+4*k;
     tmp = abs(Qex - rule(a,b,m,func));
      err = [err; tmp]; h = [h; (b-a)/m];
end
xset('thickness',3); plot(h, err, "o-c");
//
h1=[10^(-2),10^(-1)]; y1=[10^(-11),10^(-7)];
plot(h1, y1,"--k"); 
//対数目盛り
a = gca();
a.log_flags = "ll";
xset('thickness',1); xgrid(color(128,128,128));  
legend('(A)', '傾き4',4);
xlabel("log h", "fontsize", 4); ylabel("log E_h", "fontsize", 4);
//ファイル出力
xs2pdf(0,"nc1.pdf");
endfunction

// 複合公式の収束の速さ（両対数グラフ）3
function order_plot2(a,b,func,rule,Qex)
h = []; err = []; kmax=7; m=5; 
for k =1:kmax
     m=m+3*k;
     tmp = abs(Qex - rule(a,b,m,func));
      err = [err; tmp]; h = [h; (b-a)/m];
end
xset('thickness',3); 
plot(h, err, "o-c");
//
h1=[10^(-2),10^(-1)]; y1=[10^(-11),10^(-7)];
plot(h1, y1,"--k"); 
//対数目盛り
a = gca();
a.log_flags = "ll";
xset('thickness',1); xgrid(color(128,128,128));  
legend('(B)', '傾き4',4);
xlabel("log h", "fontsize", 4); ylabel("log E_h", "fontsize", 4);
//ファイル出力
xs2pdf(0,"nc2.pdf");
endfunction

// 複合公式の収束の速さ（両対数グラフ）4
function order_plot3(a,b,func,rule,Qex)
h = []; err = []; kmax=10; m=0;
//
mmin=4; mmax=80;
for m = mmin:mmax
    tmp = abs(Qex - rule(a,b,m,func));
    err = [err; tmp]; h = [h; (b-a)/m];
end
xset('thickness',2); plot(h, err, "o-c");
//対数目盛り
a = gca();
a.log_flags = "ll";
legend('(C)',4);
xset('thickness',1); xgrid(color(128,128,128));  
xlabel("log h", "fontsize", 4); ylabel("log E_h", "fontsize", 4);
//ファイル出力
xs2pdf(0,"nc3.pdf");
endfunction

//
function y = func3(x)
    y = abs(sin(%pi*x - %pi/4));
endfunction

//
Qex2 = 2/%pi;


//output Lagrange interpolation by NS
function result = output_lag2(a,b,n,func,m)
    h=(b-a)/20;
    aa=a-h; bb=b+h;
    x = linspace(aa,bb,m+1)';
    y = func(x); 
    c=min(y); d=max(y); hh=(d-c)/20;
    cc=c-hh; dd=d+hh;
    xi = linspace(a,b,n+1);
    eta = func(xi);
    yl = lag(xi,eta,x);
    //plot 
    clf();xset('thickness',2); 
    plot(x,y,'k--');  //関数と補間データと補間多項式
    xset('thickness',3); 
    plot(xi,eta,'co', x,yl,'c-');  //関数と補間データと補間多項式
   xset('thickness',1); 
    a=gca(); //アクティブな軸のオブジェクトを取得
    a.data_bounds(:,1)=[aa;bb]; //X軸の範囲を aa--bbに
    a.data_bounds(:,2)=[cc;dd]; //X軸の範囲を cc--ddに
    xlabel("x", "fontsize", 4); ylabel("y", "fontsize", 4);
    result = norm(y-yl,%inf); //l^¥infty norm
//ファイル出力
xs2pdf(0,"lag2.pdf");
endfunction

//Lagrange interpolation  with datas (xi_i,eta_i)
function y = lag(xi,eta,x)
    n = size(xi,'c');
    m = size(x,'r');
    L = ones(m,n);
    y = zeros(x);
    for i=1:n
        for j=1:n
            if j<>i then
                L(:,i) = L(:,i).*(x-xi(j))/(xi(i)-xi(j));
            end
        end
        y = y+eta(i)*L(:,i);
    end
endfunction

//
function y = func4(x)
    y = x + sin(3*x);
endfunction

//高次のニュートン・コーツ公式
function [value, err]=nc1(a,b,func,n,Qex)
//n=4,5,6,7,8のみ
//if (n>8)|(n<4) then
 //   error("n should be >8!");    
//end
//係数
select n
  case 4 then 
    w=[14/45, 64/45, 8/15, 64/45, 14/45];
  case 5 then 
    w=[95/288, 125/96, 125/144, 125/144, 125/96,95/288];
  case 6 then 
    w=[41/140, 54/35, 27/140, 68/35, 27/140, 54/35, 41/140];
  case 7 then 
    w=[5257/17280, 25039/17280, 343/640, 20923/17280, 20923/17280, 343/640, 25039/17280, 5257/17280];
  case 8 then 
    w=[3956/14175, 23552/14175, -3712/14175, 41984/14175, -18160/14175, 41984/14175, -3712/14175, 23552/14175, 3956/14175];
  else
    error("n should be 3> and < 9 !");
end
//
h = (b-a)/n;  x=linspace(a,b,n+1)'; val=func(x);
value=h*w*val; err=(Qex-value)/Qex;
endfunction

//
function y = func5(x)
    y = 1 ./ (x.^2+1);
endfunction
//
Qex5=2*atan(5.0);

// weights of Gauss-Legendre formula
function [a, b] = coeflege(n)
    if(n<=1) disp('n must be > 1'); return; end
    for k = 1:n 
        a(k) = 0; b(k) = 0; 
    end 
    b(1)=2;
    for k = 2:n 
        b(k) = 1/(4-1/(k-1)^2); 
    end
endfunction

// 
function [x,w] = gaulege(n)
    if(n<=1) disp('n must be > 1'); return; end
    [a, b] = coeflege(n); 
    JacM = diag(a) + diag(sqrt(b(2:n)),1) +  diag(sqrt(b(2:n)),-1);
    [w, x] = spec(JacM); x = diag(x); scal = 2; w = w(1,:)'.^2*scal;
    [x,ind] = gsort(x); w=w(ind); 
endfunction

// 複合シンプソン則（パラメータ付き）
function integ = c_simps1(a,b,m,func,t0)
    x = linspace(a,b,m+1);
    integ = 0;
    for i=1:m
        h = x(i+1)-x(i); 
        c = (x(i) + x(i+1))/2;
        integ = integ + (func(x(i),t0) + 4*func(c,t0) + func(x(i+1),t0))*h/6;
    end
endfunction

//
function y = func7(x,t0)
    y = 1 ./ sqrt(1-sin(t0/2).^2 .* sin(x).^2);
endfunction

// Gauss-Chebyshev formula
function integ = GC(func,n)
    t = [0:n];
    x = cos((t+0.5)*%pi/(n+1))
    w = %pi/(n+1)
    integ = 0;
    for i=1:n+1
        integ = integ + w*func(x(i));
    end
endfunction

//
function y = func51(x)
    y = (5*sqrt(1-x.^2)) ./ (25*x.^2+1);
endfunction

//
function [res] = error_GC1(wfunc,a,b,func,Qex)
    res=[];
    for n=4:2:40
    approx = GC(wfunc,n);
    approx1 =  c_simps(a,b,n/2,func);
    err = abs(approx - Qex)/abs(Qex);
    err1 = abs(approx1 - Qex)/abs(Qex);
    res=[res;n,approx,err,approx1,err1];
    end
endfunction
