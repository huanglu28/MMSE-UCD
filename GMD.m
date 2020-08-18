function C = GMD( H,x )
% ���ξ�ֵ�ֽ����MIMOϵͳ�еķ����ź�
% H -- NR*NTά�����ŵ�
% x -- �����ź�
% c -- �����ź�
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
