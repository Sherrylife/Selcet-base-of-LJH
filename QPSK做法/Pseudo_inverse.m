function result=Pseudo_inverse(H,x)
%功能：伪逆算法解调QPSK
%参数说明：
%result ---解调结果
%H      ---瑞利衰落信道矩阵
%x ---接收端数据
%estimate_s ---估计信号
%result     ---判决后的结果

%求解H的伪逆矩阵
G=(H'*H)\H';
%得到估计信号
estimate_s=G*x;
%判决
[~,n_T]=size(H); 
result=zeros(n_T,1);
for ii=1:n_T
    real_s=real(estimate_s(ii));    imag_s=imag(estimate_s(ii));
    if (real_s>0&&imag_s>0) %第一象限   
            result(ii)=2^(-0.5)*(1+1i);   
    end
    if( real_s<0&&imag_s>0 )%第二象限
             result(ii)=2^(-0.5)*(-1+1i);   
    end
    if (real_s<0&&imag_s<0) %第三象限
            result(ii)=2^(-0.5)*(-1-1i);   
    end
    if  (real_s>0&&imag_s<0) %第四象限
            result(ii)=2^(-0.5)*(1-1i);
    end
end
% disp('伪逆:估计信号');disp(estimate_s);
end