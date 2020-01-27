function result=MMSE_detector(H,receiver_x,sigma)
%���ܣ�����MMSE���������QPSK�����ź�
%����˵����
%result ---������
%H      ---����˥���ŵ�����
%H_     ---��������������
%send   ---���Ͷ�����
%noise  ---��˹������
%receiver_x ---���ն�����
%receiver_x_---�����Ľ��ն˾���
%sigma      ---��˹�������ķ��ȡֵΪ1
%estimate_s ---�����ź�
%result     ---�о���Ľ��
%n_T        ---4������

%���ݾ��������С��ԭ������H��reveiver_x
[n_R,n_T]=size(H);       H_=H;
H_((1+n_R):(n_T+n_R),:)=sigma*eye(n_T);
receiver_x_=receiver_x;
receiver_x_((1+n_R):(n_R+n_T),1)=zeros(n_T,1);
%�˳������ź�
estimate_s=(H_'*H_)\(H_'*receiver_x_);
%�о�
result=zeros(4,1);
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
end