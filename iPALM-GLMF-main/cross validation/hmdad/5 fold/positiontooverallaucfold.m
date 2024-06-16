function overallauc=positiontooverallaucfold(position)
%��֪����������disease-microbe��������������֤ʱ������������������Ӧ����disease-microbe����������֤������AUCֵ
%��ѡ����Ϊ���м�����֪����ص�microbe
%��ȡ����disease-microbe��������������֤ʱ�������
load position
position1=position;
for i=1:100
    position=position1(i,:);
%��ȡdisease-microbe���������ڽӾ���,����disease������microbe
load interaction;
%n��ʾdisease��Ŀ��m��ʾmicrobe��Ŀ
[n,m]=size(interaction);
%��¼��������֪disease-microbe������disease ID��microbeID��ÿһ�еĵ�һ������Ϊmicrobe��ID���ڶ�������Ϊdisease��ID
sID=textread('microbe-disease association.txt');
%ppΪdisease-microbe������Ŀ
[pp,qq]=size(sID);

%���㿴ǰkλ��ʱ���TPR��FPR
for k=1:m*n-floor(pp/5)*4
    tp=0;
    for t=1:pp
        if position(1,t)<=k
            tp=tp+1;
        end
    end
    tpr(1,k)=tp/pp;
fp=k*pp-tp;
  
     fpr(1,k)=fp/(floor(pp/5)*4*(m*n-pp+floor(pp/5)-1)+(pp-floor(pp/5)*4)*(m*n-floor(pp/5)*4-1));
end

plot(fpr, tpr);

%����ÿ��С���ε����
clear area;
area(1,1)=tpr(1,1)*fpr(1,1)/2;
for k=2:m*n-floor(pp/5)*4
    area(1,k)=[tpr(1,k-1)+tpr(1,k)]*[fpr(1,k)-fpr(1,k-1)]/2;
end
plot(fpr,tpr,'b')
overallauc(i)=sum(area);
end
save overallauc overallauc
end
          


