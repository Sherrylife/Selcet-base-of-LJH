function [H,send,noise]=Generate_Datas(n_R,n_T,sigma)
%功能：产生QPSK调制下的发送端和接收端数据
%变量说明：
%   H---瑞利衰落信道
%   send ---QPSK调制信号
%   noise      ---高斯白噪声
%   n_R         ---接收天线数量
%   n_T         ---发送天线数量
X=randn(n_R,n_T);         Y=randn(n_R,n_T);
s1=2^(-0.5)*randsrc(n_T,1);    s2=2^(-0.5)*randsrc(n_T,1);
send=s1+1i*s2;      %QPSK调制信号
H=X+1i*Y;            %4*4瑞利衰落信道,方差为1
%noise=wgn(n_R,1,2*sigma^2,'linear');       %高斯白噪声,平均功率为1W    
% noise=zeros(n_R,1);
%  复高斯白噪声
% wgn(m,n,10*log10(p))%产生一个m行n列的高斯白噪声的矩阵,所有元素构成的样点的信号功率为p瓦
N=1000;
%方法1
% noise_temp=wgn(N,1,2*sigma^2,'linear')+1j*wgn(N,1,2*sigma^2,'linear');
noise_temp=wgn(N,1,10*log10(sigma^2/2) )+1j*wgn(N,1,10*log10(sigma^2/2) );
noise=noise_temp(1:n_R);
% %方法2
% noise_temp = sqrt(sigma)*randn(1,N)+1j*sqrt(sigma)*randn(1,N);
% noise=noise_temp(1:4)';
%验证所产生的噪声是否正确
% noise_Power=sum(abs(noise_temp).^2)/length(noise_temp)   % 验证
% mean=sum(noise)/N
%方法3
% noise=normrnd(0,2*sigma,[n_R 1])+1j*normrnd(0,2*sigma,[n_R 1]);
end
