function y = FRAC(G1, G2)
%% FRAC
for i = 1:size(G1,1)
    for j = 1:size(G1,2)
        num = abs(squeeze(G1(i,j,:))'*squeeze(G2(i,j,:))).^2;
        den = norm(squeeze(G1(i,j,:))).^2 * norm(squeeze(G2(i,j,:))).^2  ;
        y(i,j) = num/den;
    end
end

y = mean(mean((y)));
