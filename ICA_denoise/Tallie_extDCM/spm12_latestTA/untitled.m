function interp_badTA(X)

x = X;
x(1,1) = NaN;
x(10,10) = NaN;
x(1,10) = NaN;
x([20 21 22],[20 21 22]) = NaN;
x([10 11 12],1) = NaN;

[xi,xj] = find(isnan(x));
xi = xi+1;
xj = xj+1;

x2 = nan(size(x,1)+2,size(x,2)+2);
x2(2:end-1,2:end-1) = x;

w = 1;
c1 = 1;
while w==1
    for k = 1:length(xi)
        x2(xi(k),xj(k)) = nanmean(nanmean(x2((xi(k)-1):(xi(k)+1),(xj(k)-1):(xj(k)+1))));
    end
    if any(any(isnan(x2(2:end-1,2:end-1))))
        disp(['running again (' num2str(c1) ')'])
        c1 = c1+1;
    else w = 0;
    end
end
x2 = x2(2:end,2:end);


end


%fig; subplot(221), imagesc(X), subplot(222), imagesc(x), subplot(224), imagesc(x2)

