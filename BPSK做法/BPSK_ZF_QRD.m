function result=BPSK_ZF_QRD(H,x)
%���ܣ�����������QR�ֽ�������㷨
%����˵����
%result ---������
%H      ---����˥���ŵ�����
%x      ---���ն�����
%estimate_s ---�����ź�
%result     ---�о���Ľ��
%Q          ---H��QR�ֽ�
%R          ---H��QR�ֽ�
%p          ---��¼����
%d_k          ---��k·�ź����ܵ�����������ȫ���źŵĸ���

%��ʼ��
 [~,n_T]=size(H); R=zeros(n_T);   Q=H;    %p=1:n_T;
 for ii=1:n_T
    R(ii,ii)=sqrt(Q(:,ii)'*Q(:,ii));
    Q(:,ii)=Q(:,ii)/R(ii,ii);
    for jj=(ii+1):n_T
        R(ii,jj)=Q(:,ii)'*Q(:,jj);
        Q(:,jj)=Q(:,jj)-R(ii,jj)*Q(:,ii);
    end
 end
 %SIC
 y=Q'*x;    estimate_s=zeros(n_T,1);    result=zeros(n_T,1);   
 for k=n_T:(-1):1
    if (k+1<=n_T)
        d_k=R(k,(k+1):n_T)*estimate_s((k+1):n_T,1);
    else
        d_k=0;
    end  
    estimate_s(k,1)=(y(k)-d_k)/R(k,k);
    %�о�
    if (estimate_s(k)>0)
        result(k)=1;        estimate_s(k)=1;
    else
        result(k)=-1;       estimate_s(k)=-1;
    end
 end
 %disp('QRD:�����ź�');disp(estimate_s);
 end

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 