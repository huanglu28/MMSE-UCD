%模拟MIMO系统 
%发射天线数NT,接收天线数NR，发射矩阵长度L
%NR>NT
%x=H*c+v
NT=10;
NR=10;
L=2000;
SNR=0;%信噪比（dB）
x_real=randint(NT,L);%NT*L发射信号
x=zeros(NT,L);%经UCD算法计算后得到的发射信号

%实际发射信号的0转化为-1,1保持1
X=(-1).^(x_real+1);  

%%%%%%%%%%%%%%MIMO信道传输
%信道模型：y=HFx+z
y=zeros(NR,L);

%生成 H:快衰弱的NR*NT*L维瑞利信道
H=randn(NR,NT,L)+1i*randn(NR,NT,L);

%服从均值为0,方差为1的正态分布的NR*1维的高斯白噪声v
v=randn(NR,L)+1i*randn(NR,L);

%%%%%%%%%%%%%%%%%%%%5%%%%%%% UCD算法 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('UCD算法');
m=SNR;
snr=10^(m/10);
z=1/snr*v;
    
    %%%%%生成F矩阵%%%%%
    F=zeros(NT,NT,L);
    
    %信道容量
    C_IT=zeros(1,L);
    C_it=zeros(1,L);
    for(l=1:L)
        
    %============== transmitter：模拟发送端 ==============%
        %step1：对H进行svd分解
        [U,S,V]=svd(H(:,:,l));

        %step2：功率分配矩阵I(单位阵代替)
        I=1;

        %step3：
        M=S*I;

        %step4:
        J=[U*M;sqrt(1/snr)*ones(NT,NT)];
        [~,M_J,~]=svd(J);
        M_J=M_J(1:NT,1:NT);
        
        %step5:对M_J几何均值分解
        [Q,R,P]=gmdv(M_J);
        
        %step6:产生F
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