function overallauc=positiontooverallauc1()
%��֪����������disease-microbe��������������֤ʱ������������������Ӧ����disease-microbe����������֤������AUCֵ
%��ѡ����Ϊ���м�����֪����ص�microbe
%��ȡ����disease-microbe��������������֤ʱ�������
load new_position.mat;
%��ȡdisease-microbe���������ڽӾ���,����disease������microbe
load new_interaction;
%n��ʾdisease��Ŀ��m��ʾmicrobe��Ŀ  ��disbiome�����ݿ��ά��
n=218;
m=1052;
%[n,m]=size(new_interaction);
%��¼��������֪disease-microbe������disease ID��microbe ID��ÿһ�еĵ�һ������Ϊmicrobe��ID���ڶ�������Ϊdisease��ID
sID=textread('microbe-disease association.txt');
%ppΪ������Ŀ
[pp,qq]=size(sID);


for i=1:pp
if new_sorted_data(i)>m*n-pp+1
new_sorted_data(i)=m*n-pp+1;
end
end
%���㿴ǰkλ��ʱ���TPR��FPR
for k=1:m*n-pp+1
    tp=0;
    for t=1:pp
        if new_sorted_data(1,t)<=k
            tp=tp+1;
        end
    end
    tpr(1,k)=tp/pp;
    fp=k*pp-tp;
    fpr(1,k)=fp/(pp*(m*n-pp));
end
plot(fpr,tpr,'b')
%����ÿ��С���ε����
clear area;
area(1,1)=tpr(1,1)*fpr(1,1)/2;
for k=2:m*n-pp+1
    area(1,k)=[tpr(1,k-1)+tpr(1,k)]*[fpr(1,k)-fpr(1,k-1)]/2;
end
%����AUCֵ
overallauc=sum(area);
end
          
