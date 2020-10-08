function A=Al(l,N)
%  A=zeros(size(l));
%  for j=1:length(l)
%  if mod(l(j)*pi, 2*pi)-mod(l(j)*pi./N, 2*pi)==0
%      A(j)=1;
%  end
A=zeros(size(l));
A(mod(l*pi./N, pi)==0)=1;
end