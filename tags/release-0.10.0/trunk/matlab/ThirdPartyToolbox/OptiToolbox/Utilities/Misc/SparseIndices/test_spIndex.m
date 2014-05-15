%% Compile spIndex
% spIndex simply prints or returns sparse indices as MEX sees them.
clc
cdir = cd;
try
    cd Utilities/Misc/SparseIndices
    mex -largeArrayDims spIndex.cpp
    cd(cdir);
catch
    cd(cdir);
end


%% Testing
clc
H = speye(3);
A = sparse([1 1 1;3 -2 -3; 1 -3 2]); 

disp('H');
spIndex(H)
disp('A');
spIndex(A)
disp('At');
spIndex(A')

[jc,ir,pr] = spIndex(A');