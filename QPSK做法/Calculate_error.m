function count=Calculate_error(result,send)
%����˵��������result��Ӧ�Ķ����Ʊ��������send�Ķ����Ʊ���Ĵ������
%����˵��:
%   result      ---���������QPSK�ź�
%   send        ---ʵ�ʷ��͵�QPSK�ź�
%   count       ---����źŵĴ������
count=0;    total=length(result);
for ii=1:total
   if (result(ii)~=send(ii))
        if (-result(ii)==send(ii))  %��������������
            count=count+2;
        else %ֻ��һ������
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