clc;clear;
rng('shuffle');
number=10000;   n_R=6;  n_T=4;       point_number=10;
E_b=1.5;  
% sigma=[0.7071    0.6302    0.5617    0.5006    0.4462    0.3976    0.3544    0.3159    0.2815    0.2509    0.2236]; 
sigma=sqrt(0.5*E_b*10.^((-3):3/(point_number-1):(0)));     
BER_MMSE=zeros(1,point_number);    BER_QRD=zeros(1,point_number);   BER_SQRD=zeros(1,point_number);
BER_PS=zeros(1,point_number);    BER_zfQRD=zeros(1,point_number);    BER_zfSQRD=zeros(1,point_number);
BER_VB=zeros(1,point_number);
for it=1:point_number
for ig=1:number
[H,s,n]=BPSK_generate(n_R,n_T,sigma(it));
x=H*s+n;
% result_MMSE=BPSK_MMSE(H,x,sigma(it));
% result_QRD=BPSK_QRD(H,x,sigma(it));
% result_SQRD=BPSK_SQRD(H,x,sigma(it));
% BER_MMSE(it)=BER_MMSE(it)+sum(result_MMSE==s);
% BER_QRD(it)=BER_QRD(it)+sum(result_QRD==s);
% BER_SQRD(it)=BER_SQRD(it)+sum(result_SQRD==s);
result_zfQRD=BPSK_ZF_QRD(H,x);
result_zfSQRD=BPSK_ZF_SQRD(H,x);
BER_zfQRD(it)=BER_zfQRD(it)+sum(result_zfQRD==s);
BER_zfSQRD(it)=BER_zfSQRD(it)+sum(result_zfSQRD==s);
end
end
%º∆À„ŒÛ¬Î¬ 
% err_MMSE=ones(1,point_number)-BER_MMSE/(n_T*number);   
% err_QRD=ones(1,point_number)-BER_QRD/(n_T*number);
% err_SQRD=ones(1,point_number)-BER_SQRD/(n_T*number);   
% err_PS=BER_PS/(n_T*number); 
% err_VB=BER_VB/(n_T*number); 
err_zfQRD=ones(1,point_number)-BER_zfQRD/(n_T*number); 
err_zfSQRD=ones(1,point_number)-BER_zfSQRD/(n_T*number); 


%Plot
xx=10*log10(E_b./(2*(sigma.^2)));
% semilogy(xx,err_MMSE,'r-^');  hold on;
% semilogy(xx,err_QRD,'gs-');   hold on;
% semilogy(xx,err_SQRD,'b-d');  hold on;
% semilogy(xx,err_PS,'y--o');  hold on;
% semilogy(xx,err_VB,'k--');  hold on;
semilogy(xx,err_zfQRD,'--*');  hold on;
semilogy(xx,err_zfSQRD,'o-');  hold on;
% legend('Pseudo inverse','V-BLAST','ZF-QRD','ZF-SQRD');
grid on;


















