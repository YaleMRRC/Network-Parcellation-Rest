function f = myfunc(D,Si,n)
    if (isempty(Si))
        f = 0;
    else
        f = loss(D,0,n) - loss(D,[Si,n+1],n);
    end
end