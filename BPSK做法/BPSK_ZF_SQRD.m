function result=BPSK_ZF_SQRD(H,x)
%功能：基于排序QR分解的迫零算法
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
 [~,n_T]=size(H); R=zeros(n_T);   Q=H;    p=1:n_T;
 for ii=1:n_T
    %排序和交换
    temp=Q(:,ii:n_T);       [~,k_ii]=min(sum(temp.*conj(temp)));
    k_ii=k_ii+ii-1;         Q(:,[ii,k_ii])=Q(:,[k_ii,ii]);
    R(:,[ii,k_ii])=R(:,[k_ii,ii]);  p([k_ii,ii])=p([ii,k_ii]);
    %
    R(ii,ii)=sqrt(Q(:,ii)'*Q(:,ii));     Q(:,ii)=Q(:,ii)/R(ii,ii);
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
    if(estimate_s(k,1)>0)
        result(k)=1;    estimate_s(k)=1;
    else
        result(k)=-1;   estimate_s(k)=-1;
    end
 end
%  还原顺序
 for ii=1:n_T        %还原顺序
    index=find(p==ii);  p([ii,index])=p([index,ii]);
    result([ii,index])=result([index,ii]);
 end
 %disp('SQRD:估计信号');disp(estimate_s);
end

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 