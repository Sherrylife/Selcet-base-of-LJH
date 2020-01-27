function result=BPSK_ZF_QRD(H,x)
%功能：基于无排序QR分解的迫零算法
%参数说明：
%result ---解调结果
%H      ---瑞利衰落信道矩阵
%x      ---接收端数据
%estimate_s ---估计信号
%result     ---判决后的结果
%Q          ---H的QR分解
%R          ---H的QR分解
%p          ---记录排序
%d_k          ---第k路信号所受到的来自其他全部信号的干扰

%初始化
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
    %判决
    if (estimate_s(k)>0)
        result(k)=1;        estimate_s(k)=1;
    else
        result(k)=-1;       estimate_s(k)=-1;
    end
 end
 %disp('QRD:估计信号');disp(estimate_s);
 end

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 