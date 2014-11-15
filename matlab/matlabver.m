function vstr = matlabver()
A=ver;
for i=1:length(A)
    if strcmp(A(i).Name,'MATLAB')
        vstr = A(i).Release;
        break
    end
end
    
