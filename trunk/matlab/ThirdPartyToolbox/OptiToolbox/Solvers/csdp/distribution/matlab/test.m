load control1.mat
whos

%%
clc
pars.objtol=1.0e-9;
[x,y,z,info]=csdp(At',b,c,K,pars)



%%
system('csdp.exe test.dat-s test.sol')