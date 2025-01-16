function v0 = getv0(i0,vc,n,m)

v0 = reshape(vc,m,n/m);  
v0 = v0(i0,:);      
v0 = kron(v0,ones(m/size(v0,1),1));
v0 = reshape(v0,n,1); 