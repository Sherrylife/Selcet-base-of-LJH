function estimate_s=VBLAST(copy_H,copy_x)
%���ܣ�V-BLAST�㷨���QPSK
%����˵����
%result ---������
%H      ---����˥���ŵ�����
%x      ---���ն�����
%estimate_s ---�����ź�
%result     ---�о���Ľ��
%G          ---H��α�����
%k          ---G�о�����С�ж��������±�

%��ʼ��
x=copy_x;        H=copy_H;      [~,n_T]=size(H);  
estimate_s=zeros(n_T,1);    p=1:n_T;
for ii=1:n_T
   G=(H'*H)\H';     
   [~,k]=min(sum(G.*conj(G),2));
   estimate_s(p(k))=G(k,:)*x;     
   %�о�
    real_s=real(estimate_s(p(k)));    imag_s=imag(estimate_s(p(k)));
    if (real_s>0&&imag_s>0) %��һ����   
            estimate_s(p(k))=2^(-0.5)*(1+1i);   
    end
    if( real_s<0&&imag_s>0 )%�ڶ�����
             estimate_s(p(k))=2^(-0.5)*(-1+1i);   
    end
    if (real_s<0&&imag_s<0) %��������
            estimate_s(p(k))=2^(-0.5)*(-1-1i);   
    end
    if  (real_s>0&&imag_s<0) %��������
            estimate_s(p(k))=2^(-0.5)*(1-1i);
    end
    %����x�;���H��˳��p
    x=x-H(:,k)*estimate_s(p(k));
    H(:,k)=[];  p(k)=[];
end
%��ԭ˳��
end






















