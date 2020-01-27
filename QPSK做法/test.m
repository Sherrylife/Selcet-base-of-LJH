clc;clear;
rng('shuffle');
number=100000;   n_R=4;  n_T=4;       point_number=15;
E_b=1.9;  sigma=sqrt(0.5*E_b*10.^((-3):3/(point_number-1):(0))); 
BER_MMSE=zeros(1,point_number);    BER_QRD=zeros(1,point_number);   BER_SQRD=zeros(1,point_number);
BER_PS=zeros(1,point_number);    BER_zfQRD=zeros(1,point_number);    BER_zfSQRD=zeros(1,point_number);
BER_VB=zeros(1,point_number);   BER_PSA=zeros(1,point_number);  BER_MMSE_BLAST=zeros(1,point_number);
for it=1:point_number
for ig=1:number
[H,s,n]=Generate_Datas(n_R,n_T,sigma(it));
%H=sqrt(E_b/N0(it))*H;
x=H*s+n;
result_MMSE=MMSE_detector(H,x,sigma(it));
result_MMSE_BLAST=MMSE_BLAST(H,x,sigma(it));
result_QRD=MMSE_QRD(H,x,sigma(it));
result_SQRD=MMSE_SQRD(H,x,sigma(it));
result_PSA=MMSE_PSA(H,x,sigma(it));
% result_ps=Pseudo_inverse(H,x);
% result_VB=VBLAST(H,x);
% result_zfQRD=ZF_QRD(H,x);
% result_zfSQRD=ZF_SQRD(H,x);
BER_MMSE(it)=BER_MMSE(it)+Calculate_error(result_MMSE,s);
BER_MMSE_BLAST(it)=BER_MMSE_BLAST(it)+Calculate_error(result_MMSE_BLAST,s);
BER_QRD(it)=BER_QRD(it)+Calculate_error(result_QRD,s);
BER_SQRD(it)=BER_SQRD(it)+Calculate_error(result_SQRD,s);
BER_PSA(it)=BER_PSA(it)+Calculate_error(result_PSA,s);
% BER_PS(it)=BER_PS(it)+Calculate_error(result_ps,s);
% BER_VB(it)=BER_VB(it)+Calculate_error(result_VB,s);
% BER_zfQRD(it)=BER_zfQRD(it)+Calculate_error(result_zfQRD,s);
% BER_zfSQRD(it)=BER_zfSQRD(it)+Calculate_error(result_zfSQRD,s);
end
end
%º∆À„ŒÛ¬Î¬ 
err_MMSE=BER_MMSE/(n_T*2*number);   
err_MMSE_BLAST=BER_MMSE_BLAST/(n_T*2*number);   
err_QRD=BER_QRD/(n_T*2*number);
err_SQRD=BER_SQRD/(n_T*2*number); 
err_PSA=BER_PSA/(n_T*2*number); 
% err_PS=BER_PS/(n_T*2*number); 
% err_VB=BER_VB/(n_T*2*number); 
% err_zfQRD=BER_zfQRD/(n_T*2*number); 
% err_zfSQRD=BER_zfSQRD/(n_T*2*number); 


%Plot
xx=10*log10(E_b./(2*(sigma.^2)));
semilogy(xx,err_MMSE,'r-^');  hold on;
semilogy(xx,err_MMSE_BLAST,'g-^');  hold on;
semilogy(xx,err_QRD,'gs-');   hold on;
semilogy(xx,err_SQRD,'b-d');  hold on;
semilogy(xx,err_PSA,'o-');  hold on;
% semilogy(xx,err_PS,'b--o');  hold on;
% semilogy(xx,err_VB,'k--');  hold on;
% semilogy(xx,err_zfQRD,'--*');  hold on;
% semilogy(xx,err_zfSQRD,'c--.');  hold on;
% legend('linear-MMSE','MMSE-QRD','MMSE-SQRD','MMSE-PSA','Pseudo inverse','V-BLAST','ZF-QRD','ZF-SQRD');
% legend('Pseudo inverse','V-BLAST','ZF-QRD','ZF-SQRD');
legend('linear-MMSE','MMSE-BLAST','MMSE-QRD','MMSE-SQRD','MMSE-PSA');
xlabel('\bf{$E_b/N_0$ in dB}');
ylabel('BER');
title('MMSEºÏ≤‚À„∑®');
grid on;


















