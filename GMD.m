function C = GMD( H,x )
% 几何均值分解解码MIMO系统中的发射信号
% H -- NR*NT维瑞利信道
% x -- 接收信号
% c -- 解码信号
[NR,NT,L]=size(H);
c=zeros(NT,L);
    C=zeros(NT,L);
for j=1:L
    HH=H(:,:,j);
    K=NT;
    [Q,R,P]=gmdv(HH);
    y=Q'*x(:,j);
    S=zeros(K,1);
    S(K)=y(K)/R(K,K);
    for n=K-1:-1:1
        sum=0;
        for o=n+1:K
            sum=sum+R(n,o)*S(o);
        end
        S(n)=(y(n)-sum)/R(n,n);
    end
    C(:,j)=P*S;
    C(:,j)=(C(:,j)>=0)-(C(:,j)<0)+0;
end
C=(C+1)/2;
end
