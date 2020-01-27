function [H,send,noise]=BPSK_generate(n_R,n_T,sigma)
%功能：产生QPSK调制下的发送端和接收端数据
%变量说明：
%   H---瑞利衰落信道
%   send ---QPSK调制信号
%   noise      ---高斯白噪声
%   n_R         ---接收天线数量
%   n_T         ---发送天线数量
H=normrnd(0,1,[n_R,n_T]);
send=randsrc(n_T,1);
N=1000;
% noise_temp=wgn(N,1,10*log10(sigma^2/2));
% noise=noise_temp(1:n_R);
noise_temp=normrnd(0,sigma,[N 1]);
noise=noise_temp(1:n_R,1);
% noise=normrnd(0,sigma,[n_R 1]);
% noise=zeros(n_R,1);
end
