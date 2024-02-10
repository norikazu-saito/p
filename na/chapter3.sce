//第3章のプログラム
//N. Saito 2016.02.02

// CG iteration
function [vx,it,vp,vr,vres,al,be] = cg(A,b,tol)
[n,n0]=size(A); vx=[]; vres=[]; it=[]; al=[]; be=[]; vp=[]; vr=[];
// initial setting
x = zeros(n,1);  vx=[vx,x];//初期値x0
r=b-A*x; res=norm(r); vr=[vr,r]; vres=[vres,res]; //残差ベクトル
p=r; vp=[vp,p]; //探索ベクトル
kit = 0;  it=[it,kit]; al=[al,0]; be=[be,0];
kmax=5*n; 
while kit < kmax & res > tol
    kit = kit + 1; 
    q=A*p; alpha=res^2/((p')*q); x=x+alpha*p; r=r-alpha*q;
    res_new=norm(r); beta=-(res_new/res)^2; p=r-beta*p;
    res=res_new; 
    vx=[vx,x]; vr=[vr,r]; vp=[vp,p]; vres=[vres,res]; al=[al,alpha]; be=[be,beta];it=[it,kit]; 
end
//plot2d(it,vres);
endfunction
