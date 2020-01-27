function [H,send,noise]=BPSK_generate(n_R,n_T,sigma)
%���ܣ�����QPSK�����µķ��Ͷ˺ͽ��ն�����
%����˵����
%   H---����˥���ŵ�
%   send ---QPSK�����ź�
%   noise      ---��˹������
%   n_R         ---������������
%   n_T         ---������������
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
