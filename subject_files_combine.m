function subject_files_combine
% 
% load Subject Subject
% count1=1;
% count2=1;
% count3=1;
% count4=1;
% count5=1;
% count6=1;
% for cs=1:length(Subject)
%     if(cs<=6)
%         Subject1(cs)=Subject(cs);
%     elseif(cs>6)&&(cs<=12)
%         Subject2(count1)=Subject(cs);
%         count1=count1+1;
%     elseif(cs>12)&&(cs<=18)
%         Subject3(count2)=Subject(cs);
%         count2=count2+1;
%     elseif(cs>18)&&(cs<=24)
%         Subject4(count3)=Subject(cs);
%         count3=count3+1;
%     elseif(cs>24)&&(cs<=30)
%         Subject5(count4)=Subject(cs);
%         count4=count4+1;
%     elseif(cs>30)&&(cs<=36)
%         Subject6(count5)=Subject(cs);
%         count5=count5+1;
%     elseif(cs>36)&&(cs<=42)
%         Subject7(count6)=Subject(cs);
%         count6=count6+1;
%     end
% end
% 
% save Subject1 Subject1
% save Subject2 Subject2
% save Subject3 Subject3
% save Subject4 Subject4
% save Subject5 Subject5
% save Subject6 Subject6
% save Subject7 Subject7

load Subject1 Subject1
load Subject2 Subject2
load Subject3 Subject3
load Subject4 Subject4
load Subject5 Subject5
load Subject6 Subject6
load Subject7 Subject7
count=1;
for cs=1:length(Subject1)
    Subject(count)=Subject1(cs);
    count=count+1;
end
for cs=1:length(Subject2)
    Subject(count)=Subject2(cs);
    count=count+1;
end
for cs=1:length(Subject3)
    Subject(count)=Subject3(cs);
    count=count+1;
end
for cs=1:length(Subject4)
    Subject(count)=Subject4(cs);
    count=count+1;
end
for cs=1:length(Subject5)
    Subject(count)=Subject5(cs);
    count=count+1;
end
for cs=1:length(Subject6)
    Subject(count)=Subject6(cs);
    count=count+1;
end
for cs=1:length(Subject7)
    Subject(count)=Subject7(cs);
    count=count+1;
end

save Subject Subject
save Subject Subject




% load Subject_RI Subject_RI
% count1=1;
% count2=1;
% count3=1;
% count4=1;
% count5=1;
% for cs=1:length(Subject_RI)
%     if(cs<=6)
%         Subject_RI1(cs)=Subject_RI(cs);
%     elseif(cs>6)&&(cs<=12)
%         Subject_RI2(count1)=Subject_RI(cs);
%         count1=count1+1;
%     elseif(cs>12)&&(cs<=16)
%         Subject_RI3(count2)=Subject_RI(cs);
%         count2=count2+1;
%     elseif(cs>16)&&(cs<=22)
%         Subject_RI4(count3)=Subject_RI(cs);
%         count3=count3+1;
%     elseif(cs>22)&&(cs<=28)
%         Subject_RI5(count4)=Subject_RI(cs);
%         count4=count4+1;
%     elseif(cs>28)&&(cs<=34)
%         Subject_RI6(count5)=Subject_RI(cs);
%         count5=count5+1;
%     end
% end
% 
% save Subject_RI1 Subject_RI1
% save Subject_RI2 Subject_RI2
% save Subject_RI3 Subject_RI3
% save Subject_RI4 Subject_RI4
% save Subject_RI5 Subject_RI5
% save Subject_RI6 Subject_RI6


load Subject_RI1 Subject_RI1
load Subject_RI2 Subject_RI2
load Subject_RI3 Subject_RI3
load Subject_RI4 Subject_RI4
load Subject_RI5 Subject_RI5
load Subject_RI6 Subject_RI6

count=1;
for cs=1:length(Subject_RI1)
    Subject_RI(count)=Subject_RI1(cs);
    count=count+1;
end
for cs=1:length(Subject_RI2)
    Subject_RI(count)=Subject_RI2(cs);
    count=count+1;
end
for cs=1:length(Subject_RI3)
    Subject_RI(count)=Subject_RI3(cs);
    count=count+1;
end
for cs=1:length(Subject_RI4)
    Subject_RI(count)=Subject_RI4(cs);
    count=count+1;
end
for cs=1:length(Subject_RI5)
    Subject_RI(count)=Subject_RI5(cs);
    count=count+1;
end
for cs=1:length(Subject_RI6)
    Subject_RI(count)=Subject_RI6(cs);
    count=count+1;
end


save Subject_RI Subject_RI
save Subject_RI Subject_RI
        