//第5章のプログラム
//NS 2016.02.29

//計算機イプシロン
function eps = mepsilon(eps0)
eps = eps0; 
while((1.0+eps)-1.0 > 0.0)
    eps = eps/2.0; 
end 
eps = 2.0*eps;
endfunction

//関数のグラフ
function plot_func(a,b,func,N)
    x=linspace(a,b,N)';
    y=func(x);
    xset('thickness',2); plot(x,y,"-c");
    a = gca(); xset('thickness',1); 
    xlabel("x", "fontsize", 4); ylabel("y", "fontsize", 4);
//ファイル出力
    xs2pdf(0,"p1.pdf");
endfunction

//
function [y]=func1(x)
    y=x .^7 - 7 * x .^ 6 + 21* x .^ 5 - 35 * x .^ 4 + 35* x.^ 3 - 21* x .^ 2 + 7*x-1+15;
endfunction

function [y]=func2(x)
y=x-7; y=y.*x+21; y=y.*x-35; y=y.*x+35; y=y.*x-21; y=y.*x+7;  y=y.*x-1+15;
//y=((((((x  - 7) .* x + 21) .* x  - 35) .* x  + 35) .* x - 21) .* x  + 7) .* x-1+15;
endfunction

function [y]=func3(x)
    y=(x-1) .^7+15;
endfunction

//
function [res]=dcos(n)
    res=[]; //空のベクトルを用意
     a=1.0; h=0.2; 
    for i=1:n
        h=h/2; val= (cos(a+h)-cos(a))/h; 
        res=[res;h,val];
    end
endfunction

function [res]=dcos2(n)
    res=[]; //空のベクトルを用意
     a=1.0; h=0.2; 
    for i=1:n
        h=h/2; hh=h/2; val= -(1/hh)*sin(hh)*sin(1+hh); 
        res=[res;h,val];
    end
endfunction
