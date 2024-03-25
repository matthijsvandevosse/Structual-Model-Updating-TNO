function y = CSAC(G1, G2)

%% CSAC
for w = 1:size(G1,3)
    for j = 1:size(G1,2)
        num = abs(squeeze(G1(:,j,w))'*squeeze(G2(:,j,w))).^2;
        den =(squeeze(G1(:,j,w))'*squeeze(G1(:,j,w)))  * (squeeze(G2(:,j,w))'*squeeze(G2(:,j,w)))  ;
        y(w,j) = num/den;
    end
end



figure
semilogx(mean(y,2))
y = mean(mean((y)));


end