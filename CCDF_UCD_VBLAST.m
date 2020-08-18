%ģ��MIMOϵͳ 
%����������NT,����������NR��������󳤶�L
%NR>NT
%x=H*c+v
NT=10;
NR=10;
L=2000;
SNR=0;%����ȣ�dB��
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
m=SNR;
snr=10^(m/10);
z=1/snr*v;
    
    %%%%%����F����%%%%%
    F=zeros(NT,NT,L);
    
    %�ŵ�����
    C_IT=zeros(1,L);
    C_it=zeros(1,L);
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
        
        HH=H(:,:,l);
        FF=F(:,:,l);
        C_IT(l)=log2(det(eye(NR)+snr*HH*FF*FF'*HH'));
        C_it(l)=log2(det(ones(NR)+snr*HH*HH'));
    end
    figure(1);
 
%     model1=cdfplot(C_IT);
%     set(model1, 'LineStyle', '-', 'Color', 'b');

% %     
     model2=cdfplot(C_it);
%     set(model2, 'LineStyle', '.-', 'Color', 'g');

    legend('UCD','GMD');
    xlabel('Capacity(bit/sec/Hz)');
    ylabel('CDF');
    title('NR=10,NT=10.SNR=0 dB'); 
    hold on;