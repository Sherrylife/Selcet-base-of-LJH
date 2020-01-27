function result2=BPSK_SQRD(H,receiver_x,sigma)
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
% %PSA���
k_min=n_T;    Q1=Q_((1):(n_R),1:n_T); Q2=Q_((1+n_R):(n_R+n_T),1:n_T);  
for ii=n_T:-1:2
    error=zeros(1,ii);  
    for r=1:ii
        error(r)=Q2(r,1:ii)*Q2(r,1:ii)';
    end
    [~,k_i]=min(error);     k_min=min([k_min,k_i]);
    if (k_i<ii)     %����
        Q2([ii,k_i],:)=Q2([k_i,ii],:);  p([ii,k_i])=p([k_i,ii]);
    end
    if (k_min<ii)  
        %�����˹�����¾���Rotate
        a=Q2(ii,k_min:(ii));      e=[zeros(1,ii-k_min),1];
        u=(a-sqrt(a*a')*e)/sqrt((a-sqrt(a*a')*e)*(a-sqrt(a*a')*e)');
        w=(u*a')/(a*u');
        Rotate=eye(ii-k_min+1)-(1+w)*(u'*u);
        %����Q1��Q2����
        Q2(1:ii,k_min:ii)=Q2(1:ii,k_min:ii)*Rotate;
        Q1(:,k_min:ii)=Q1(:,k_min:ii)*Rotate;
    end
end
R_=sigma*eye(n_T)/(Q2);
Q_(1:n_R,:)=Q1; Q_((1+n_R):(n_R+n_T),:)=Q2;
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
    if(estimate_s(ii,1)>0)
        result2(ii)=1;
    else
        result2(ii)=-1;
    end
end
for ii=1:n_T        %��ԭ˳��
    index=find(p==ii);  p([ii,index])=p([index,ii]);
    result2([ii,index])=result2([index,ii]);
end

end



















































































