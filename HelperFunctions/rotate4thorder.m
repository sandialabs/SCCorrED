function A_prime = rotate4thorder(A,q)
    A_prime = zeros(3,3,3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    for m=1:3
                        for n=1:3
                            for o=1:3
                                for p=1:3
                                    A_prime(i,j,k,l) = A_prime(i,j,k,l) + A(m,n,o,p)*q(i,m)*q(j,n)*q(k,o)*q(l,p);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end