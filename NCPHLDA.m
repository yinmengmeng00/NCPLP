function [SM,SD]=NCPHLDA(interaction,sd,sl)
% [SM,SD]=NCPHLDA(knownre,AD,AM);

[nd,nl]=size(interaction);
SM = zeros(nd,nl)+10^-30;
% NCP= zeros(nd,nl);
SD = zeros(nd,nl)+10^-30;

for i = 1:nd
    for j = 1:nl
       SD(i,j)=sd(i,:)*interaction(:,j)/norm(interaction(:,j));
    end
end

for i = 1:nd
      for j = 1:nl
         SM(i,j)=interaction(i,:)*sl(:,j)/norm(interaction(i,:));
     end
end
save('SD.mat','SD');
save('SM.mat','SM');

% for i = 1:nd
%       for j = 1:nl
%          NCP(i,j)=NCP_d(i,j)/norm(sd(i,:));
%      end
% end
% for i = 1:nd
%       for j = 1:nl
%          NCP(i,j)=NCP_l(i,j)/norm(sl(:,j));
%      end
% end
% for i= 1:nd
%     for j= 1:nl
%         NCP(i,j)=NCP_d(i,j)+NCP_l(i,j)/(norm(sd(i,:))+norm(sl(:,j)));
%     end
% end
end