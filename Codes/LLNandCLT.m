clear all
% Tan Bui Sep 2015
% law of large numbers and central limit theorem

Compute_Flag = 2; 
% 1: law of large numbers
% 2: Central limite theorem

if Compute_Flag == 1,
% law of large numbers for uniform random numbers
N = 1:1000:100001; n = length(N); S = zeros(n,1);
for i=1:length(N)
    S(i) = sum(rand(N(i),1))/N(i);
end
figure
axes('fontsize',12)
plot(1:n,S,'b--s','linewidth',2)
xlabel('N'); ylabel('S')
grid on
title('Law of large numbers')

else
%stop
% Central limit theorem for N(a,sigma^2)
a = 2; sigma = 2.5;
N = 100000;
nreplicate = 1000;

for i = 1:nreplicate
   Z(i) =  1/sigma/sqrt(N)*sum(a + sigma*randn(N,1)) - a/sigma*sqrt(N);
%   Z(i) =  1/sigma/sqrt(N)*sum(rand(N,1)) - a/sigma*sqrt(N);
end
histogram(Z)
fprintf('meanZ = %f and covZ = %f\n', mean(Z), cov(Z))

end
