function L = loss(D,Si,n)
    if Si == 0
%         sprintf('hi')
        L = (1/n)*sum(D(1:n,n+1));
    else
%         sprintf('bye')
        L = (1/n)*(sum(min(D(Si,1:n),[],1)));
    end
end
