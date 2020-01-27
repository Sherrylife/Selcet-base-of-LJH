function count=Calculate_error(result,send)
%功能说明：计算result对应的二进制编码相对于send的二进制编码的错误个数
%参数说明:
%   result      ---解调出来的QPSK信号
%   send        ---实际发送的QPSK信号
%   count       ---解调信号的错误个数
count=0;    total=length(result);
for ii=1:total
   if (result(ii)~=send(ii))
        if (-result(ii)==send(ii))  %错两个二进制码
            count=count+2;
        else %只错一个编码
            count=count+1;
        end
   end
end
end
% count=0;    total=length(result);
% for ii=1:total
%     if (result(ii)~=send(ii))
%         count=count+1;
%     end
% end
% end