//第6章のプログラム
//NS 2016.02.29

//LU分解の計算時間
function lutime(N)
    x=[]; y1=[]; y2=[]; y3=[]; z=[];
    //
    for n=1:20:N+1
        //1
        A=n*eye(n,n)+rand(n,n); 
        timer();[L,U,P]=lu(A);m=timer(); x=[x;n]; y1=[y1;m]; z=[z;n^3/3]
        //2
        A=n*eye(n,n)+rand(n,n); timer();[L,U,P]=lu(A);m=timer(); y2=[y2;m];
        //3
        A=n*eye(n,n)+rand(n,n); timer();[L,U,P]=lu(A);m=timer(); y3=[y3;m];
    end
    //
    z=(3/(n^3))*m*z;
   //
    xset('thickness',2); 
    plot(x,y1,'c-');plot(x,y2,'c-');plot(x,y3,'c-');plot(x,z,'k--');
    xset('thickness',1); 
    legend("LU factorization 1","LU factorization 2","LU factorization 3","(Const)/n^3",2);
    xlabel("size n", "fontsize", 4); ylabel("CPU time [sec]", "fontsize", 4);
    //ファイル出力
    xs2pdf(0,"lutime.pdf");
endfunction

//output Lagrange interpolation
function result = output_lag6(a,b,n,func,m)
    //
    xx=linspace(a,b,n+1)'; yy=func(xx);
    A=lag_mat(xx); cc =A \ yy; pp = poly(cc','z','coeff');
    //plot 
    x=linspace(a,b,m+1)'; y=horner(pp,x);
    clf(); xset('thickness',2); 
    y1=func4(x); plot(x,y1,'k--', xx,yy,'co');  //
    plot(x,y,'c-');  //関数と補間データと補間多項式
    //描画用
    h=(x($)-x(1))/30; b1=x(1)-h; b2=x($)+h; y2=[y;y1]
    d1=min(y2); d2=max(y2); hh=(d2-d1)/30; d1=d1-hh; d2=d2+hh;
    //
    xset('thickness',1); 
    a=gca(); //アクティブな軸のオブジェクトを取得
    a.data_bounds(:,1)=[b1;b2]; //X軸の範囲
    a.data_bounds(:,2)=[d1;d2]; //Y軸の範囲
    xlabel("x", "fontsize", 4); ylabel("y", "fontsize", 4);
    result = [natori(A),cond(A,1),cond(A,2)]; //l^¥infty norm
//ファイル出力
xs2pdf(0,"lag6.pdf");
endfunction

//Lagrange matrix
function [A]= lag_mat(x)
    [n m]=size(x); A=ones(n,n); 
    for i=1:n
        for j=2:n
            A(i,j)=x(i)^(j-1);
         end
    end
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

//Lagrange interpolation of a function func on an interval [a,b] 
function y = lag_f(a,b,n,func,x)
    xi = linspace(a,b,n+1);
    eta = func(xi);
    y = lag(xi,eta,x);
endfunction

//
function y = func4(x)
    y = x + sin(3*x);
endfunction

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

//名取塚本の方法
function [num]=natori(A)
////LU分解
    [L,U,P]=lu(A); [n,m]=size(A);
////上からの後退代入でyの成分を決めながらwを決める
    w=[]; w=[w;1/U(1,1)];
    for k=2:n
        t=0;
        for j=1:k-1
            t = t + U(j,k)*w(j);
         end
         w=[w; (-sign(t)-t)/U(k,k)]
    end
////後退代入
// //   z = P \ (L' \ w);
    z = P'  * (L' \ w);
////条件数
    num=norm(A,1)*norm(z,%inf);
endfunction

//Hilbert matrix
function [A]= hil_mat(n)
    A=ones(n,n);
    for i=1:n
        for j=1:n
            A(i,j)=1/(i+j-1);
         end
    end
endfunction

//
function [A]= fra_mat(n)
    A=ones(n,n);
    for i=1:n
        for j=1:n
            A(i,j)=n+1-max(i,j);
         end
    end
endfunction

//漸化式
function rec1(N, iv1, iv2)
x = linspace(1,N,N)';
x1=iv1; x2=iv2; y=[x1;x2];
for n = 3:N
    x3=func(x1,x2); y=[y;x3]; x1=x2; x2=x3;
end
//scf() a new window for drawing, 
scf(); z=log(y); 
xset('thickness',3);  plot(x, z, "-oc");
// label
xset('thickness',1); 
xlabel("k", "fontsize", 4); ylabel("log(x_k)", "fontsize", 4);
//ファイル出力
xs2pdf(0,"rec1.pdf");
endfunction

//
function [y] =func(x1,x2)
    y=(7/3)*x2 - (2/3)*x1;
    //y=2.25*x2-0.5*x1;
endfunction
