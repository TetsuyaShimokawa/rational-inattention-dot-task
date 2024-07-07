function DOTtask_some_additional

 load Subject Subject
 
 % make list for model selection in Section 4.2 and Supplemental Materials
 make_list(Subject); 
 
 % number of Dot and accuracy of answers in Section4.5
 analysis_treatment_accuracy(Subject);
end

function make_list(Subject)
for cs=1:length(Subject)
    % id
    A(cs,1)=cs;
    
    %NIAC,NIAS
    A(cs,2)=Subject(cs).NIAS;
    A(cs,3)=Subject(cs).NIAC;
    
    % Shannon entropy model rambda
    A(cs,4)=Subject(cs).lambda;    
    % Shannon entropy model  AIC
    A(cs,5)=Subject(cs).aic_MM;
    % Shannon entropy model  llh
    A(cs,6)=Subject(cs).llh_MM;
    
    % Tallis entropy model parameters
    A(cs,7)=Subject(cs).T.lambda;
    A(cs,8)=Subject(cs).T.sigma;       
    % Tallis entropy model  AIC
    A(cs,9)=Subject(cs).T.aic_MM;
    % Tallis entropy model  llh
    A(cs,10)=Subject(cs).T.llh_MM;
end

%subject average in Section 4.2 
nanmean(A(:,[4:end]))

% list of all subjects in Supplemental Materials
A

% write it down to excel file
%xlswrite('DOT_result.xlsx',A,'all')
%xlswrite('DOT_result.xlsx',nanmean(A),'average')
end

function analysis_treatment_accuracy(Subject)

% Make Treatment groups
count(1)=0;count(2)=0;count(3)=0;
for cs=1:length(Subject)
    for ct=1:length(Subject(cs).init.Trials)
        if(Subject(cs).init.Trials(ct).numb_dot==60)
            ctg=1;
        elseif(Subject(cs).init.Trials(ct).numb_dot==80)
            ctg=2;
        elseif(Subject(cs).init.Trials(ct).numb_dot==120)
            ctg=3;
        end
        count(ctg)=count(ctg)+1;
        if([Subject(cs).ans(ct).red,1-Subject(cs).ans(ct).red]*Subject(cs).trials(ct).rtn_red_blue'>0)
            TRT{ctg}.accuracy(count(ctg),1)=1;
        else
            TRT{ctg}.accuracy(count(ctg),1)=0;
        end
    end
end
% Accuracy Rate
disp('Dot and Accuracy')
[nanmean(TRT{1}.accuracy),nanmean(TRT{2}.accuracy),nanmean(TRT{3}.accuracy)]
 end