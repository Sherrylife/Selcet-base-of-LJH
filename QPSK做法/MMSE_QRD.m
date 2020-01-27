function result1=MMSE_QRD(H,receiver_x,sigma)
%功能：基于无排序的QR分解的MMSE准则解调
%参数说明：
%   result1  ---解码结果
%   H       ---瑞利衰落信道矩阵
%   H_      ---扩充后的瑞利矩阵
%   reveiver_x  ---接收到的信号
%   receiver_x_ ---扩充后的接收信号
%   Q_          ---H_的QR分解的Q矩阵
%   R_          ---H_的QR分解的R矩阵
%   sigma       ---高斯白噪声的方差
%   norm           ---储存Q矩阵的列元素的平方和
%   estimate_s  ---估计信号

%初始化
[n_R,n_T]=size(H);      receiver_x_=receiver_x;   
result1=zeros(n_T,1);            H_=H;   
H_((1+n_R):(n_R+n_T),:)=sigma*eye(n_T);
Q_=H_;      R_=zeros(n_T);   receiver_x_((1+n_R):(n_R+n_T),1)=zeros(n_T,1);

%计算Q_矩阵的列的二范数的平方
norm=zeros(n_T,1);
for jj=1:n_T
   norm(jj)=Q_(:,jj)'*Q_(:,jj); 
end
%求解H_的QR分解矩阵
for ii=1:n_T
    R_(ii,ii)=sqrt(norm(ii));
    Q_(:,ii)=Q_(:,ii)/R_(ii,ii);
    for k=(ii+1):n_T
        R_(ii,k)=(Q_(:,ii))'*Q_(:,k);
        Q_(:,k)=Q_(:,k)-R_(ii,k)*Q_(:,ii);
        norm(k)=norm(k)-R_(ii,k)*R_(ii,k)';
    end
end

%求解估计信号
y=Q_'*receiver_x_;  estimate_s=zeros(n_T,1);
for ii=n_T:-1:1
    if ((ii+1)<=n_T)
        d_ii=R_(ii,(ii+1):n_T)*estimate_s((ii+1):n_T,1);
    else
        d_ii=0;
    end
    estimate_s(ii,1)=(y(ii)-d_ii)/R_(ii,ii);
    %判决
    real_s=real(estimate_s(ii));    imag_s=imag(estimate_s(ii));
    if (real_s>0&&imag_s>0) %第一象限   
            result1(ii)=2^(-0.5)*(1+1i);   estimate_s(ii)=result1(ii);
    end
    if( real_s<0&&imag_s>0 )%第二象限
             result1(ii)=2^(-0.5)*(-1+1i);  estimate_s(ii)=result1(ii); 
    end
    if (real_s<0&&imag_s<0) %第三象限
            result1(ii)=2^(-0.5)*(-1-1i);   estimate_s(ii)=result1(ii);
    end
    if  (real_s>0&&imag_s<0) %第四象限
            result1(ii)=2^(-0.5)*(1-1i);    estimate_s(ii)=result1(ii);
    end
end
end