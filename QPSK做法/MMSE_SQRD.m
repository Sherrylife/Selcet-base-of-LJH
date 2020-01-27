function result2=MMSE_SQRD(H,receiver_x,sigma)
%���ܣ����������QR�ֽ��MMSE׼����(���PSA����㷨)
%����˵����
%   result2  ---������
%   H       ---����˥���ŵ�����
%   H_      ---��������������
%   reveiver_x  ---���յ����ź�
%   receiver_x_ ---�����Ľ����ź�
%   Q_          ---H_��QR�ֽ��Q����
%   R_          ---H_��QR�ֽ��R����
%   sigma       ---��˹�������ķ���
%   norm           ---����Q�������Ԫ�ص�ƽ����
%   estimate_s  ---�����ź�
%   Q1          ---Q_((1):(n_R),1:n_T)
%   Q2          ---Q_((n_R+1):(n_R+n_T),1:n_T)
%   error       ---����Q2�����������Ԫ�ص�ƽ����
%   k_min       ---error����Сֵ���ڵ��±�
%   Rotate           ---ĳ�������ĺ�˹�����¾���
%   p           ---��¼������˳��

%��ʼ��
[n_R,n_T]=size(H);      receiver_x_=receiver_x;   
result2=zeros(n_T,1);            H_=H;   
H_((1+n_R):(n_R+n_T),:)=sigma*eye(n_T);
Q_=H_;      R_=zeros(n_T);   receiver_x_((1+n_R):(n_R+n_T),1)=zeros(n_T,1);
p=1:n_T;
%����Q_������еĶ�������ƽ��
norm=zeros(n_T,1);
for jj=1:n_T
   norm(jj)=Q_(:,jj)'*Q_(:,jj); 
end
%���H_��QR�ֽ����
for ii=1:n_T
    %����
    temp=norm(ii:n_T);   [~,min_index]=min(temp);
    min_index=min_index+ii-1;
    %����
    R_(:,[ii,min_index])=R_(:,[min_index,ii]);
    p([ii,min_index])=p([min_index,ii]);   
    norm([ii,min_index])=norm([min_index,ii]);   
    Q_(1:(n_R+ii-1),[ii,min_index])=Q_(1:(n_R+ii-1),[min_index,ii]);
    %��������
    R_(ii,ii)=sqrt(norm(ii));       Q_(:,ii)=Q_(:,ii)/R_(ii,ii);
    for k=(ii+1):n_T
        R_(ii,k)=(Q_(:,ii))'*Q_(:,k);
        Q_(:,k)=Q_(:,k)-R_(ii,k)*Q_(:,ii);
        norm(k)=norm(k)-R_(ii,k)*R_(ii,k)';
    end
end

%�������źŲ���ԭ����
y=Q_'*receiver_x_;  estimate_s=zeros(n_T,1);
for ii=n_T:-1:1
    if ((ii+1)<=n_T)
        d_ii=R_(ii,(ii+1):n_T)*estimate_s((ii+1):n_T,1);
    else
        d_ii=0;
    end
    estimate_s(ii,1)=(y(ii)-d_ii)/R_(ii,ii);
    %�о�
    real_s=real(estimate_s(ii));    imag_s=imag(estimate_s(ii));
    if (real_s>0&&imag_s>0) %��һ����   
            result2(ii)=2^(-0.5)*(1+1i);   
            estimate_s(ii)=result2(ii);
    end
    if( real_s<0&&imag_s>0 )%�ڶ�����
             result2(ii)=2^(-0.5)*(-1+1i);  
             estimate_s(ii)=result2(ii);
    end
    if (real_s<0&&imag_s<0) %��������
            result2(ii)=2^(-0.5)*(-1-1i);   
            estimate_s(ii)=result2(ii);
    end
    if  (real_s>0&&imag_s<0) %��������
            result2(ii)=2^(-0.5)*(1-1i);
            estimate_s(ii)=result2(ii);
    end
end
for ii=1:n_T        %��ԭ˳��
    index=find(p==ii);  p([ii,index])=p([index,ii]);
    result2([ii,index])=result2([index,ii]);
end

end



















































































