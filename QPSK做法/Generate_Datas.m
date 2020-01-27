function [H,send,noise]=Generate_Datas(n_R,n_T,sigma)
%���ܣ�����QPSK�����µķ��Ͷ˺ͽ��ն�����
%����˵����
%   H---����˥���ŵ�
%   send ---QPSK�����ź�
%   noise      ---��˹������
%   n_R         ---������������
%   n_T         ---������������
X=randn(n_R,n_T);         Y=randn(n_R,n_T);
s1=2^(-0.5)*randsrc(n_T,1);    s2=2^(-0.5)*randsrc(n_T,1);
send=s1+1i*s2;      %QPSK�����ź�
H=X+1i*Y;            %4*4����˥���ŵ�,����Ϊ1
%noise=wgn(n_R,1,2*sigma^2,'linear');       %��˹������,ƽ������Ϊ1W    
% noise=zeros(n_R,1);
%  ����˹������
% wgn(m,n,10*log10(p))%����һ��m��n�еĸ�˹�������ľ���,����Ԫ�ع��ɵ�������źŹ���Ϊp��
N=1000;
%����1
% noise_temp=wgn(N,1,2*sigma^2,'linear')+1j*wgn(N,1,2*sigma^2,'linear');
noise_temp=wgn(N,1,10*log10(sigma^2/2) )+1j*wgn(N,1,10*log10(sigma^2/2) );
noise=noise_temp(1:n_R);
% %����2
% noise_temp = sqrt(sigma)*randn(1,N)+1j*sqrt(sigma)*randn(1,N);
% noise=noise_temp(1:4)';
%��֤�������������Ƿ���ȷ
% noise_Power=sum(abs(noise_temp).^2)/length(noise_temp)   % ��֤
% mean=sum(noise)/N
%����3
% noise=normrnd(0,2*sigma,[n_R 1])+1j*normrnd(0,2*sigma,[n_R 1]);
end
