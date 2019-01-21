function S_opt_all = exemplar(sqD,sqd,t,n,Karg,l)

K=Karg;
S_opt=[];
sort_ground = 1:n;
for k=1:K
    display(k)
    if k==1           
        rho = inf;
        f = zeros(n,l);
        f_sum=0;
        for i=sort_ground  
            for subj=1:l
                D = sqD{subj};
                d0 = sqd{subj};
                D(n+1,1:n) = d0;
                D(1:n,n+1) = d0;


                f1(i,subj,k) = myfunc(D,S_opt,n);
                S = [S_opt,i];               
                f1(i,subj,k+1) = myfunc(D,S,n);
                f(i,subj) = f1(i,subj,k+1) - f1(i,subj,k);
            end
            f_sum(i) = sum(f(i,:),2);
        end
        [rho, sort_ground] = sort(f_sum,'descend');
        i_opt = sort_ground(1);
        sort_ground(1)=[];
        rho(1)=[];
        S_opt = [S_opt,i_opt];


    elseif k>1
        f = zeros(n,l);
        f_sum=0;
        i = sort_ground(1);
        for subj=1:l
            D = sqD{subj};
            d0 = sqd{subj};
            D(n+1,1:n) = d0;
            D(1:n,n+1) = d0;
            S = [S_opt,i];                

            f1(i,subj,k+1) = myfunc(D,S,n);
            f(i,subj) = f1(i,subj,k+1) - f1(i_opt,subj,k);
        end
        f_sum(i) = sum(f(i,:),2);
        if f_sum(i) > rho(2)
            i_opt = i;
            S_opt = [S_opt,i_opt];
            continue;    
        else
            f = zeros(n,l);
            for i=sort_ground
                for subj=1:l
                    D = sqD{subj};
                    d0 = sqd{subj};
                    D(n+1,1:n) = d0;
                    D(1:n,n+1) = d0;
                    S = [S_opt,i];                

                    f1(i,subj,k+1) = myfunc(D,S,n);
                    f(i,subj) = f1(i,subj,k+1) - f1(i_opt,subj,k);

                end           
                f_sum(i) = sum(f(i,:),2);
            end
            [rho, sort_ground] = sort(f_sum,'descend');
            i_opt = sort_ground(1);                
            sort_ground(1)=[];
            rho(1)=[];
            S_opt = [S_opt,i_opt];
        end                
    end
end
S_opt_all=S_opt;

