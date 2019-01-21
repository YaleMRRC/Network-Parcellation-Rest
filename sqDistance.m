function D = sqDistance(X)
D = bsxfun(@plus,dot(X,X,1)',dot(X,X,1))-2*(X'*X);