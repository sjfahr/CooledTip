%This parses the power further

function [P]=power_parser(power_array);

power_temp(:,1)=power_array(:,4); %record time
power_temp(:,2)=power_array(:,5); % Off = 0; On = 1
power_temp(:,3)=power_array(:,6); % Power %
power_temp(:,4)=power_temp(:,3)*15/100;  %Convert % into W

%Find where the power changes at all;
P=(diff(power_array(:,5))~=0)+(diff(power_array(:,6))~=0);
P(1+1:end+1,:) = P(1:end,:);  %add row back in because diff function eliminates the first one.
P(1,:)=P(2,:);



(diff(power_log(:,5))~=0)+(diff(power_log(:,6))~=0)~=0
% j=1;
% for i=1:size(power_temp,1)
% while power_temp(i,3)==power(i+1,3) && power_temp(i,2)==power_temp(i+1,2)
%     P(j+1)=P(j);
% end
%     
    