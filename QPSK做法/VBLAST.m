function estimate_s=VBLAST(copy_H,copy_x)
%功能：V-BLAST算法解调QPSK
%参数说明：
%result ---解调结果
%H      ---瑞利衰落信道矩阵
%x      ---接收端数据
%estimate_s ---估计信号
%result     ---判决后的结果
%G          ---H的伪逆矩阵
%k          ---G中具有最小列二范数的下标

%初始化
x=copy_x;        H=copy_H;      [~,n_T]=size(H);  
estimate_s=zeros(n_T,1);    p=1:n_T;
for ii=1:n_T
   G=(H'*H)\H';     
   [~,k]=min(sum(G.*conj(G),2));
   estimate_s(p(k))=G(k,:)*x;     
   %判决
    real_s=real(estimate_s(p(k)));    imag_s=imag(estimate_s(p(k)));
    if (real_s>0&&imag_s>0) %第一象限   
            estimate_s(p(k))=2^(-0.5)*(1+1i);   
    end
    if( real_s<0&&imag_s>0 )%第二象限
             estimate_s(p(k))=2^(-0.5)*(-1+1i);   
    end
    if (real_s<0&&imag_s<0) %第三象限
            estimate_s(p(k))=2^(-0.5)*(-1-1i);   
    end
    if  (real_s>0&&imag_s<0) %第四象限
            estimate_s(p(k))=2^(-0.5)*(1-1i);
    end
    %更新x和矩阵H和顺序p
    x=x-H(:,k)*estimate_s(p(k));
    H(:,k)=[];  p(k)=[];
end
%还原顺序
end






















