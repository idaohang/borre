%Solution to Opgave 6.5
format rat
[U,S,V] = svd([1 1 1 1])
Aplus = V*[.5; 0; 0; 0]*U'
          
[U,S,V] = svd([0 1 0;1 0 0])
Bplus = V*S'*U'

[U,S,V] = svd([1 1;0 0])
Cplus = V*[1/S(1,1) 0; 0 0]*U'
%%%%%%%%%%%%%% end opg65.m  %%%%%%%%%%%%%%%%%%%%%%%