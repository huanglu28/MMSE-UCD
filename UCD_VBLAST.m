%ģ��MIMOϵͳ 
%����������NT,����������NR��������󳤶�L
%NR>NT
%x=H*c+v
NT=4;
NR=6;
L=1000;
SNR=[0:3:30];%����ȣ�dB��
x_real=randint(NT,L);%NT*L�����ź�
x=zeros(NT,L);%��UCD�㷨�����õ��ķ����ź�

%ʵ�ʷ����źŵ�0ת��Ϊ-1,1����1
X=(-1).^(x_real+1);  

%%%%%%%%%%%%%%MIMO�ŵ�����
%�ŵ�ģ�ͣ�y=HFx+z
y=zeros(NR,L);

%���� H:��˥����NR*NT*Lά�����ŵ�
H=randn(NR,NT,L)+1i*randn(NR,NT,L);

%���Ӿ�ֵΪ0,����Ϊ1����̬�ֲ���NR*1ά�ĸ�˹������v
v=randn(NR,L)+1i*randn(NR,L);

%%%%%%%%%%%%%%%%%%%%5%%%%%%% UCD�㷨 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('UCD�㷨');
erate=[];
for m=SNR
    snr=10^(m/10);
    z=1/snr*v;
    
    %%%%%����F����%%%%%
    F=zeros(NT,NT,L);
    for(l=1:L)
        
    %============== transmitter��ģ�ⷢ�Ͷ� ==============%
        %step1����H����svd�ֽ�
        [U,S,V]=svd(H(:,:,l));

        %step2�����ʷ������I(��λ�����)
        I=1;

        %step3��
        M=S*I;

        %step4:
        J=[U*M;sqrt(1/snr)*ones(NT,NT)];
        [~,M_J,~]=svd(J);
        M_J=M_J(1:NT,1:NT);
        
        %step5:��M_J���ξ�ֵ�ֽ�
        [Q,R,P]=gmdv(M_J);
        
        %step6:����F
        F(:,:,l)=V*P;
        
        %���������ź�y
        y(:,l)=H(:,:,l)*F(:,:,l)*X(:,l)+z(:,l);
    
     %=============== receiver:ģ����ն� ==============%
        %step7:����G=HF ��QR�ֽ�
        Q_GA=U*M*inv(M_J)*Q;
        %step8
        HH=Q_GA*R;
        G=inv(HH'*HH+(1/snr)*eye(NT))*HH';
        w=G*y(:,l);
        x(:,l)=(w>=0)-(w<0)+0;
    end
    x=(x+1)/2;
    %����UCD�㷨��������
    [errbit,err_ratio]=biterr(x_real,x);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'d-b'); %��ɫ����
hold on;

%%%%%%%%%%%%%%%%%%%%5%%%%%%% GMD�㷨 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('GMD�㷨');
x_GMD=zeros(NT,L);%��GMD�㷨�����õ��ķ����ź�
erate_gmd=[];
for m=SNR
    snr=10^(m/10);
    z=1/snr*v;
    for l=1:L
        HH=H(:,:,l)*F(:,:,l);
        y(:,l)=HH*X(:,l)+z(:,l);
        x_GMD(:,l)=GMD(HH,y(:,l));
    end
    [errbit,err_ratio]=biterr(x_real,x_GMD);
    erate_gmd=[erate_gmd,err_ratio];
end
semilogy(SNR,erate_gmd,'x-r'); %��ɫ����
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%% MMSE�㷨 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('MMSE�㷨');
%��ͬ������µ�������
erate_mmse=[];
x_MMSE=zeros(NT,L);
%��������
for m=SNR
    snr=10^(m/10);
    z=1/snr*v;
    %������õ����ź�
    for l=1:L
        y(:,l)=H(:,:,l)*X(:,l)+z(:,l);
    end
    x=MMSE(H,y,snr);
    %����V-blast�㷨��������
    [errbit,err_ratio]=biterr(x_real,x);
    erate_mmse=[erate_mmse,err_ratio];
end
semilogy(SNR,erate_mmse,'k'); %��ɫֱ��
hold on;

xlabel('SNR/dB');
ylabel('BER');
title('NT=4��NR=6ʱUCD�㷨�������ʺ�����ȹ�ϵ����');
legend('UCD-VBLAST','GMD-VBLAST','MMSE-VBLAST')
