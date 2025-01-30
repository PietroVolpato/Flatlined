function Dnew=dynupdate(D,pp,params)

delta = params.chi*(params.phi*ffree(D,params) + (1-params.phi)*fbmi(pp)); 
Dnew = D + delta;