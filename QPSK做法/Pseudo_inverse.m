function result=Pseudo_inverse(H,x)
%���ܣ�α���㷨���QPSK
%����˵����
%result ---������
%H      ---����˥���ŵ�����
%x ---���ն�����
%estimate_s ---�����ź�
%result     ---�о���Ľ��

%���H��α�����
G=(H'*H)\H';
%�õ������ź�
estimate_s=G*x;
%�о�
[~,n_T]=size(H); 
result=zeros(n_T,1);
for ii=1:n_T
    real_s=real(estimate_s(ii));    imag_s=imag(estimate_s(ii));
    if (real_s>0&&imag_s>0) %��һ����   
            result(ii)=2^(-0.5)*(1+1i);   
    end
    if( real_s<0&&imag_s>0 )%�ڶ�����
             result(ii)=2^(-0.5)*(-1+1i);   
    end
    if (real_s<0&&imag_s<0) %��������
            result(ii)=2^(-0.5)*(-1-1i);   
    end
    if  (real_s>0&&imag_s<0) %��������
            result(ii)=2^(-0.5)*(1-1i);
    end
end
% disp('α��:�����ź�');disp(estimate_s);
end