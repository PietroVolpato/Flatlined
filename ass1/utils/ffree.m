function F=ffree(D,p)

if D<0.5-p.omega
    F=-sin(pi/(0.5-p.omega)*D);
elseif D>0.5+p.omega
    F=sin(pi/(0.5-p.omega)*(D-0.5-p.omega));
else
    F=-p.psi*sin(pi/p.omega*(D-0.5));
end