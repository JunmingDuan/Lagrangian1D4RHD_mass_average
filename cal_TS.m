clear all;
%syms ul zl rl gl el pl sl us
%syms ur zr rr gr er pr sr
%eqn = ((ul - zl/(rl*gl))*us - 1)*((er*zr+pr*ur)*us - sr*zr - pr) ...
%- ((ur - zr/(rr*gr))*us - 1)*((el*zl+pl*ul)*us - sl*zl - pl);
%eqn = simplify(expand(eqn))
%collect(eqn, us)

syms r u p g G
A = [-1/g/r^2,   -u*g/r.   0;
     -

