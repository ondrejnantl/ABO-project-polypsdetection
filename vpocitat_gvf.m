function [u,v] = vpocitat_gvf(energie,iterace_gvf,krok_gvf)

[dx,dy]=gradient(energie);
u = dx; v = dy;
gradient_na_2 = dx.^2 + dy.^2;

for i=1:iterace_gvf
    u = u + krok_gvf*4*del2(u) - gradient_na_2.*(u-dx);
    v = v + krok_gvf*4*del2(v) - gradient_na_2.*(v-dy);
    if (mod(i,10) == 0)
        title([num2str(i/iterace_gvf*100) '% vypoctu gvf'])
        pause(0.01)
    end
end

gradientt = sqrt(u.*u+v.*v);
u= u./(gradientt+1e-10); 
v = v./(gradientt+1e-10);

end

