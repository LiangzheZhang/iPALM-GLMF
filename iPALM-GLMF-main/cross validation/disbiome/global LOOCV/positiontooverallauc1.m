function overallauc=positiontooverallauc1()
%已知数据是所有disease-microbe关联在做交叉验证时候的排序，现在想算出对应所有disease-microbe关联交叉验证的整体AUC值
%候选集合为所有疾病已知不相关的microbe
%读取所有disease-microbe关联在做交叉验证时候的排序
load new_position.mat;
%读取disease-microbe关联网络邻接矩阵,行是disease，列是microbe
load new_interaction;
%n表示disease数目，m表示microbe数目  是disbiome的数据库的维度
n=218;
m=1052;
%[n,m]=size(new_interaction);
%记录下所有已知disease-microbe关联的disease ID和microbe ID，每一行的第一个数字为microbe的ID，第二个数字为disease的ID
sID=textread('microbe-disease association.txt');
%pp为关联数目
[pp,qq]=size(sID);


for i=1:pp
if new_sorted_data(i)>m*n-pp+1
new_sorted_data(i)=m*n-pp+1;
end
end
%计算看前k位的时候的TPR和FPR
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
%计算每个小梯形的面积
clear area;
area(1,1)=tpr(1,1)*fpr(1,1)/2;
for k=2:m*n-pp+1
    area(1,k)=[tpr(1,k-1)+tpr(1,k)]*[fpr(1,k)-fpr(1,k-1)]/2;
end
%计算AUC值
overallauc=sum(area);
end
          
