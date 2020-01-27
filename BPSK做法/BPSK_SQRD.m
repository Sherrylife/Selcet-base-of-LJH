function result2=BPSK_SQRD(H,receiver_x,sigma)
%功能：基于排序的QR分解的MMSE准则解调(外加PSA检测算法)
%参数说明：
%   result2  ---解码结果
%   H       ---瑞利衰落信道矩阵
%   H_      ---扩充后的瑞利矩阵
%   reveiver_x  ---接收到的信号
%   receiver_x_ ---扩充后的接收信号
%   Q_          ---H_的QR分解的Q矩阵
%   R_          ---H_的QR分解的R矩阵
%   sigma       ---高斯白噪声的方差
%   norm           ---储存Q矩阵的列元素的平方和
%   estimate_s  ---估计信号
%   Q1          ---Q_((1):(n_R),1:n_T)
%   Q2          ---Q_((n_R+1):(n_R+n_T),1:n_T)
%   error       ---储存Q2矩阵的行向量元素的平方和
%   k_min       ---error的最小值所在的下标
%   Rotate           ---某行向量的豪斯霍尔德矩阵
%   p           ---记录交换的顺序

%初始化
[n_R,n_T]=size(H);      receiver_x_=receiver_x;   
result2=zeros(n_T,1);            H_=H;   
H_((1+n_R):(n_R+n_T),:)=sigma*eye(n_T);
Q_=H_;      R_=zeros(n_T);   receiver_x_((1+n_R):(n_R+n_T),1)=zeros(n_T,1);
p=1:n_T;
%计算Q_矩阵的列的二范数的平方
norm=zeros(n_T,1);
for jj=1:n_T
   norm(jj)=Q_(:,jj)'*Q_(:,jj); 
end
%求解H_的QR分解矩阵
for ii=1:n_T
    %排序
    temp=norm(ii:n_T);   [~,min_index]=min(temp);
    min_index=min_index+ii-1;
    %交换
    R_(:,[ii,min_index])=R_(:,[min_index,ii]);
    p([ii,min_index])=p([min_index,ii]);   
    norm([ii,min_index])=norm([min_index,ii]);   
    Q_(1:(n_R+ii-1),[ii,min_index])=Q_(1:(n_R+ii-1),[min_index,ii]);
    %交换结束
    R_(ii,ii)=sqrt(norm(ii));       Q_(:,ii)=Q_(:,ii)/R_(ii,ii);
    for k=(ii+1):n_T
        R_(ii,k)=(Q_(:,ii))'*Q_(:,k);
        Q_(:,k)=Q_(:,k)-R_(ii,k)*Q_(:,ii);
        norm(k)=norm(k)-R_(ii,k)*R_(ii,k)';
    end
end
% %PSA检测
k_min=n_T;    Q1=Q_((1):(n_R),1:n_T); Q2=Q_((1+n_R):(n_R+n_T),1:n_T);  
for ii=n_T:-1:2
    error=zeros(1,ii);  
    for r=1:ii
        error(r)=Q2(r,1:ii)*Q2(r,1:ii)';
    end
    [~,k_i]=min(error);     k_min=min([k_min,k_i]);
    if (k_i<ii)     %交换
        Q2([ii,k_i],:)=Q2([k_i,ii],:);  p([ii,k_i])=p([k_i,ii]);
    end
    if (k_min<ii)  
        %计算豪斯霍尔德矩阵Rotate
        a=Q2(ii,k_min:(ii));      e=[zeros(1,ii-k_min),1];
        u=(a-sqrt(a*a')*e)/sqrt((a-sqrt(a*a')*e)*(a-sqrt(a*a')*e)');
        w=(u*a')/(a*u');
        Rotate=eye(ii-k_min+1)-(1+w)*(u'*u);
        %更新Q1和Q2矩阵
        Q2(1:ii,k_min:ii)=Q2(1:ii,k_min:ii)*Rotate;
        Q1(:,k_min:ii)=Q1(:,k_min:ii)*Rotate;
    end
end
R_=sigma*eye(n_T)/(Q2);
Q_(1:n_R,:)=Q1; Q_((1+n_R):(n_R+n_T),:)=Q2;
%求解估计信号并还原排序
y=Q_'*receiver_x_;  estimate_s=zeros(n_T,1);
for ii=n_T:-1:1
    if ((ii+1)<=n_T)
        d_ii=R_(ii,(ii+1):n_T)*estimate_s((ii+1):n_T,1);
    else
        d_ii=0;
    end
    estimate_s(ii,1)=(y(ii)-d_ii)/R_(ii,ii);
    %判决
    if(estimate_s(ii,1)>0)
        result2(ii)=1;
    else
        result2(ii)=-1;
    end
end
for ii=1:n_T        %还原顺序
    index=find(p==ii);  p([ii,index])=p([index,ii]);
    result2([ii,index])=result2([index,ii]);
end

end



















































































