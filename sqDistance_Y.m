function D = sqDistance_Y(X,Y)
D = bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y);