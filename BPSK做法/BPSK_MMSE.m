function result=BPSK_MMSE(H,receiver_x,sigma)
%功能：利用MMSE估计器解调QPSK调制信号
%变量说明：
%result ---解调结果
%H      ---瑞利衰落信道矩阵
%H_     ---扩充后的瑞利矩阵
%send   ---发送端数据
%noise  ---高斯白噪声
%receiver_x ---接收端数据
%receiver_x_---扩充后的接收端矩阵
%sigma      ---高斯白噪声的方差，取值为1
%estimate_s ---估计信号
%result     ---判决后的结果
%n_T        ---4根天线

%根据均方误差最小的原则扩充H和reveiver_x
[n_R,n_T]=size(H);       H_=H;
H_((1+n_R):(n_T+n_R),:)=sigma*eye(n_T);
receiver_x_=receiver_x;
receiver_x_((1+n_R):(n_R+n_T),1)=zeros(n_T,1);
%滤出估计信号
estimate_s=(H_'*H_)\(H_'*receiver_x_);
%判决
result=zeros(n_T,1);
for ii=1:n_T
   if (estimate_s(ii)>0)
       result(ii)=1;
   else
       result(ii)=-1;
   end
end
end