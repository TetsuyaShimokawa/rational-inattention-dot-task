function DOTtask_main(section_flg)

%You have to run subject_files_combine.m to reorganize Subject.mat at the very beginning.
subject_files_combine

switch section_flg
    case 'RI subject group & Shannon entropy'  %Section4.3, 4.4, 4.5, 4.6
        subject_flg='ri';entropy_flg='S';
    case 'RI subject group & Tsallis entropy'  %Section4.7
        subject_flg='ri';entropy_flg='T';
    case 'ALL subject group & Shannon entropy'  %Section4.8
        subject_flg='all';entropy_flg='S';
end

switch subject_flg
    case 'all'
         load Subject Subject
         
        % % Normalize NIRS and Gaze data.
        [Subject]=normalize_data(Subject);

        switch entropy_flg
            case 'S'
                % Parameter estimation for each individual.
                 Subject=parameter_estimation(Subject);
                
                %Section 4.8
                %Execute gaze's anova1 grouped by lambda. lambda estimation and gaze are done per individual.
                anova_gaze_Subject = anova_gaze_lambda(Subject);

                %Section 4.8
                %Execute anova1 of nirs, grouped by lambda. lambda estimation and nirs are done by individual.
                anova_nirs_Subject = anova_nirs_lambda(Subject);

                %Test for differences in lambda, gaze, and nirs grouped by Treatment (number of Dots).
                [anova_gaze_TG,anova_nirs_TG ]=analysis_treatment(Subject);

                %Treatment（Return数）ごとにグループ化して、lambda, gaze, nirsに差があるか
                [anova_gaze_RTG,anova_nirs_RTG]=analysis_return(Subject);
        
            case 'T'
                % Parameter estimation for each individual.
                 Subject_T=parameter_T_estimation(Subject);
                
                 %Execute gaze's anova1 grouped by lambda. lambda estimation and gaze are done per individual.
                anova_T_gaze_Subject = anova_gaze_lambda(Subject_T);

                %Execute anova1 of nirs, grouped by lambda. lambda estimation and nirs are done by individual.
                anova_T_nirs_Subject = anova_nirs_lambda(Subject_T);

                %Test for differences in lambda, gaze, and nirs grouped by Treatment (number of Dots).
                [anova_T_gaze_TG,anova_T_nirs_TG ]=analysis_treatment(Subject_T);

                %Test for differences in lambda, gaze, and nirs grouped by Treatment (number of Returns).
                [anova_T_gaze_RTG,anova_T_nirs_RTG]=analysis_return(Subject_T);

                S_aic=[];T_aic=[];diff_aic=[];
                for cs=1:length(Subject)
                    S_aic(cs)=Subject(cs).aic_MM;
                    T_aic(cs)=Subject_T(cs).aic_MM;
                    diff_aic(cs)=Subject(cs).aic_MM - Subject_T(cs).aic_MM;
                end
                [nanmean(S_aic),nanmean(T_aic),nanmean(diff_aic)]
                length(find(diff_aic<=0))%The number of people for whom the Shannon model is better.
        end
        

    case 'ri'
         
        % % Normalize NIRS and Gaze data.
        load Subject_RI Subject_RI
        [Subject_RI]=normalize_data(Subject_RI);
         

        switch entropy_flg
            case 'S'
                % Parameter estimation for each individual
                Subject_RI=parameter_estimation(Subject_RI);
                 
                %Section 4.4
                 %Execute gaze's anova1 grouped by lambda. lambda estimation and gaze are done per individual.
                anova_RI_gaze_Subject = anova_gaze_lambda(Subject_RI);

                %Section 4.3
                %Execute anova1 of nirs, grouped by lambda. lambda estimation and nirs are done by individual.
                anova_RI_nirs_Subject = anova_nirs_lambda(Subject_RI);

                %Section 4.5
                %Test for differences in lambda, gaze, and nirs grouped by Treatment (number of Dots).
                [anova_RI_gaze_TG,anova_RI_nirs_TG ]=analysis_treatment(Subject_RI);

                %Section 4.6
                %Test for differences in lambda, gaze, and nirs grouped by Treatment (number of Returns).
                [anova_RI_gaze_RTG,anova_RI_nirs_RTG]=analysis_return(Subject_RI);
                
            case 'T'
                % Parameter estimation for each individual
                 Subject_RI_T=parameter_T_estimation(Subject_RI);
                
                 %Section 4.7
                 %Execute gaze's anova1 grouped by lambda. lambda estimation and gaze are done per individual.
                anova_RI_T_gaze_Subject = anova_gaze_lambda(Subject_RI_T);

                %Section 4.7
                %Execute anova1 of nirs, grouped by lambda. lambda estimation and nirs are done by individual.
                anova_RI_T_nirs_Subject = anova_nirs_lambda(Subject_RI_T);
                
                %Test for differences in lambda, gaze, and nirs grouped by Treatment (number of Dots).
                [anova_RI_T_gaze_TG,anova_RI_T_nirs_TG ]=analysis_treatment(Subject_RI_T);

                %Test for differences in lambda, gaze, and nirs grouped by Treatment (number of Returns).
                [anova_RI_T_gaze_RTG,anova_RI_T_nirs_RTG]=analysis_return(Subject_RI_T);

                S_aic=[];T_aic=[];diff_aic=[];
                for cs=1:length(Subject_RI)
                    S_aic(cs)=Subject_RI(cs).aic_MM;
                    T_aic(cs)=Subject_RI_T(cs).aic_MM;
                    diff_aic(cs)=Subject_RI(cs).aic_MM - Subject_RI_T(cs).aic_MM;
                end
                [nanmean(S_aic),nanmean(T_aic),nanmean(diff_aic)]
                length(find(diff_aic<=0))%The number of people for whom the Shannon model is better.
        end
        
end
end

function [Subject]=normalize_data(Subject)

for cs=1:length(Subject)
    trials=Subject(cs).init.Trials;
    Nmat=[];
    for ct=1:length(Subject(cs).ans)
        
        for cch=1:3
            if(cch==1)
                ch=19;
            elseif(cch==2)
                ch=21;
            elseif(cch==3)
                ch=22;
            end
            clear dnirs adnirs
            dnirs=(Subject(cs).fnirs_data(ct).mat(ch,:) -Subject(cs).fnirs_data(ct).mat(ch,1))';
            adnirs=nanmean(dnirs);
            Nmat(ct,cch)=adnirs;
            
            Nmat_sum(ct,cch)=sum(dnirs);
            Nmat_max(ct,cch)=max(dnirs);
            
            cspan=[250,500,1000,1500,2000,3000];
            maxlg=length(Subject(cs).fnirs_data(ct).mat(ch,:));
            for csp=1:3%5sec,10sec,20sec
                cand_spn=cspan(csp);                
                span=min(cand_spn,maxlg);
                dnirs_spn{csp}=(Subject(cs).fnirs_data(ct).mat(ch,[1:span]) -Subject(cs).fnirs_data(ct).mat(ch,1))';
                adnirs=nanmean(dnirs_spn{csp});
                Nmat_spn{csp}(ct,cch)=adnirs;
            end                
        end
    end

    Subject(cs).average_nirs=Nmat;
    Subject(cs).normalized_average_nirs=(Nmat - nanmean(Nmat))./nanstd(Nmat);    
    for csp=1:3
        Subject(cs).average_nirs_spn{csp}=Nmat_spn{csp};
        Subject(cs).normalized_average_nirs_spn{csp}=(Nmat_spn{csp} - nanmean(Nmat_spn{csp}))./nanstd(Nmat_spn{csp});
    end
    
    Subject(cs).sum_nirs=Nmat_sum;
    Subject(cs).normalized_sum_nirs=(Nmat_sum - nanmean(Nmat_sum))./nanstd(Nmat_sum);   
    
    Subject(cs).max_nirs=Nmat_max;
    Subject(cs).normalized_max_nirs=(Nmat_max - nanmean(Nmat_max))./nanstd(Nmat_max);    
end
end

function Subject=parameter_estimation(Subject)
for cs=1:length(Subject)
    clear rtn ans_red
    for ct=1:length(Subject(cs).ans)
        ans_red(ct,1)=Subject(cs).ans(ct).red;
        rtn(ct,:)=Subject(cs).trials(ct).rtn_red_blue;
    end
    optnew = optimset('Display','notify','MaxFunEvals',2000,'TolFun',1e-6);
    params0=0;
    func=@(params)minus_log_likelyhood(rtn,ans_red,params);
    [params_ml, min_mllh, exitsig] = fminunc(func,params0,optnew);
    Subject(cs).lambda=params_ml;%this shows log(lambda)
    Subject(cs).llh_MM=-min_mllh;
    numb_param=length(params_ml);
    numb_d=length(ans_red);
    Subject(cs).aic_MM=-2*(-min_mllh)+2*numb_param;
    Subject(cs).R2_MM=-((-min_mllh)-(numb_d*log(1/2)))/(numb_d*log(1/2));
end
end
        
function [mllh]=minus_log_likelyhood(rtn,inv,params)
lambda=params(1);

t=[inv,1-inv];
[A]=rtn/exp(lambda);
[sigma_mat]=calc_softmax_prob(A);

numb_O=length(t(1,:));
numb_N=length(t(:,1));
llh=0;
for cn=1:numb_N
    for co=1:numb_O
        llh=llh + t(cn,co)*log(sigma_mat(cn,co));
    end
end
mllh=-1*llh;

end

function [sigma_mat]=calc_softmax_prob(A)
numb_O=length(A(1,:));
numb_N=length(A(:,1));

for cn=1:numb_N
    sumsigma=0;
    for co=1:numb_O        
        sumsigma=sumsigma + exp(A(cn,co));        
    end
    sumA(cn)=sumsigma;
end

for cn=1:numb_N
    for co=1:numb_O        
        sigma_mat(cn,co)=exp(A(cn,co))/sumA(cn);       
    end
end
end

function AnovaResult_gaze = anova_gaze_lambda(Subject)

for cr=1:1
    for cc=1:3
        G{cr,cc}=[];
    end
end
    

L=[];
for cs=1:length(Subject)
    L(cs)=Subject(cs).lambda;
end
[v,lambda_idx]=sort(L);
 spn=floor(length(L)/3);
 s1=spn;
 s2=2*spn;
 
 for cs=1:length(Subject)
     infoptn=1;
     cidx=find(lambda_idx==cs);
     if(cidx<=s1)
         G{infoptn,1}=[G{infoptn,1}, Subject(cs).ave_gaze_sum];
     elseif(cidx>s1)&&(cidx<=s2)
         G{infoptn,2}=[G{infoptn,2}, Subject(cs).ave_gaze_sum];
     elseif(cidx>s2)
         G{infoptn,3}=[G{infoptn,3}, Subject(cs).ave_gaze_sum];
     end
end

min_spl=min([length(G{1,1}),length(G{1,2}),length(G{1,3})]);
Gmat(:,1)=[G{1,1}([1:min_spl])]';
Gmat(:,2)=[G{1,2}([1:min_spl])]';
Gmat(:,3)=[G{1,3}([1:min_spl])]';
[AnovaResult_gaze.p,AnovaResult_gaze.tbl] = anova1(Gmat);   
disp('gaze_lambda')
nanmean(Gmat)
end

function AnovaResult_nirs = anova_nirs_lambda(Subject)

for cch=1:3
% ch=[19,21,22]
    for cr=1:1
        for cc=1:3
            G{cr,cc}=[];
        end
    end

    L=[];
    for cs=1:length(Subject)
        L(cs)=Subject(cs).lambda;
    end
    [v,lambda_idx]=sort(L);
     spn=floor(length(L)/3);
     s1=spn;
     s2=2*spn;

     for cs=1:length(Subject)
         infoptn=1;
         cidx=find(lambda_idx==cs);
         nirsd=nanmean(Subject(cs).max_nirs(:,cch)); 
         
         if(cidx<=s1)
             G{infoptn,1}=[G{infoptn,1}, nirsd];%全Treatmentの個人平均
         elseif(cidx>s1)&&(cidx<=s2)
             G{infoptn,2}=[G{infoptn,2}, nirsd];
         elseif(cidx>s2)
             G{infoptn,3}=[G{infoptn,3}, nirsd];
         end
     end

    min_spl=min([length(G{1,1}),length(G{1,2}),length(G{1,3})]);
    Gmat(:,1)=[G{1,1}([1:min_spl])]';
    Gmat(:,2)=[G{1,2}([1:min_spl])]';
    Gmat(:,3)=[G{1,3}([1:min_spl])]';  
    [AnovaResult_nirs.ch(cch).p,AnovaResult_nirs.ch(cch).tbl] = anova1(Gmat);   
    disp(['nirs_lambda ch',num2str(cch)])
    nanmean(Gmat)
end
end

function [AnovaResult_gaze_TG,AnovaResult_nirs_TG ]=analysis_treatment(Subject)

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
        TG{ctg}.ans_red(count(ctg),1)=Subject(cs).ans(ct).red;
        TG{ctg}.rtn(count(ctg),:)=Subject(cs).trials(ct).rtn_red_blue;
        TG{ctg}.Gaze(count(ctg),1)=Subject(cs).normalized_gaze_sum(ct);

        TG{ctg}.NIRS(count(ctg),:)=Subject(cs).max_nirs(ct,:);
    end
end

% Run parameter estimation
TG=parameter_estimation_TG(TG);
[TG{1}.lambda,TG{2}.lambda,TG{3}.lambda]

% ANOVA Gaze
AnovaResult_gaze_TG = anova_gaze_TG(TG);

% ANOVA NIRS
AnovaResult_nirs_TG = anova_nirs_TG(TG);
end

function TG=parameter_estimation_TG(TG)
for ctg=1:length(TG)
    clear rtn ans_red
    optnew = optimset('Display','notify','MaxFunEvals',2000,'TolFun',1e-6);
    params0=0;
    func=@(params)minus_log_likelyhood(TG{ctg}.rtn,TG{ctg}.ans_red,params);
    [params_ml, min_mllh, exitsig] = fminunc(func,params0,optnew);
    TG{ctg}.lambda=params_ml;%this shows log(lambda)
    TG{ctg}.llh_MM=-min_mllh;
    numb_param=length(params_ml);
    numb_d=length(TG{ctg}.ans_red);
    TG{ctg}.aic_MM=-2*(-min_mllh)+2*numb_param;
    TG{ctg}.R2_MM=-((-min_mllh)-(numb_d*log(1/2)))/(numb_d*log(1/2));
end
end

function AnovaResult_gaze_TG = anova_gaze_TG(TG)

for cc=1:3
    G{1,cc}=[];
end  
for ctg=1:length(TG)
    G{1,ctg}=TG{ctg}.Gaze(:,1)';
end
min_spl=min([length(G{1,1}),length(G{1,2}),length(G{1,3})]);
Gmat(:,1)=[G{1,1}([1:min_spl])]';
Gmat(:,2)=[G{1,2}([1:min_spl])]';
Gmat(:,3)=[G{1,3}([1:min_spl])]';
disp(['gaze'])
[AnovaResult_gaze_TG.p,AnovaResult_gaze_TG.tbl] = anova1(Gmat);    
nanmean(Gmat)
end

function AnovaResult_nirs_TG = anova_nirs_TG(TG)

for cch=1:3
    for cc=1:3
        G{1,cc}=[];
    end  
    for ctg=1:length(TG)
        G{1,ctg}=TG{ctg}.NIRS(:,cch)';
    end
    
    Gmat(:,1)=[G{1,1}'];
    Gmat(:,2)=[G{1,2}'];
    Gmat(:,3)=[G{1,3}'];
    [AnovaResult_nirs_TG.ch(cch).p,AnovaResult_nirs_TG.ch(cch).tbl] = anova1(Gmat);    
    disp(['nirs_ch',num2str(cch)])
    nanmean(Gmat)
end
end

function [AnovaResult_gaze_RTN,AnovaResult_nirs_RTN]=analysis_return(Subject)

count(1)=0;count(2)=0;count(3)=0;
for cs=1:length(Subject)
    for ct=1:length(Subject(cs).init.Trials)
        if(Subject(cs).init.Trials(ct).rtn==5)
            ctg=1;
        elseif(Subject(cs).init.Trials(ct).rtn==100)
            ctg=2;
        elseif(Subject(cs).init.Trials(ct).rtn==5000)
            ctg=3;
        end
        count(ctg)=count(ctg)+1;
        RTN{ctg}.ans_red(count(ctg),1)=Subject(cs).ans(ct).red;
        RTN{ctg}.rtn(count(ctg),:)=Subject(cs).trials(ct).rtn_red_blue;
        if([Subject(cs).ans(ct).red,1-Subject(cs).ans(ct).red]*Subject(cs).trials(ct).rtn_red_blue'>0)
            RTN{ctg}.accuracy(count(ctg),1)=1;
        else
            RTN{ctg}.accuracy(count(ctg),1)=0;
        end
        RTN{ctg}.Gaze(count(ctg),1)=Subject(cs).normalized_gaze_sum(ct);

        RTN{ctg}.NIRS(count(ctg),:)=Subject(cs).max_nirs(ct,:);
    end
end

% parameter estimation
RTN=parameter_estimation_TG(RTN);
[RTN{1}.lambda,RTN{2}.lambda,RTN{3}.lambda]

% Accuracy Rate
[nanmean(RTN{1}.accuracy),nanmean(RTN{2}.accuracy),nanmean(RTN{3}.accuracy)]

% ANOVA Gaze
AnovaResult_gaze_RTN = anova_gaze_TG(RTN);

% ANOVA NIRS
AnovaResult_nirs_RTN = anova_nirs_TG(RTN);
end


function Subject=parameter_T_estimation(Subject)
for cs=1:length(Subject)
    clear rtn ans_red
    for ct=1:length(Subject(cs).ans)
        ans_red(ct,1)=Subject(cs).ans(ct).red;
        rtn(ct,:)=Subject(cs).trials(ct).rtn_red_blue;
    end
    optnew = optimset('Display','notify','MaxFunEvals',500,'TolFun',1e-6);
    params0=[0,1];
    func=@(params)minus_log_likelyhood_T(rtn,ans_red,params);
    [params_ml, min_mllh, exitsig] = fminunc(func,params0,optnew);
    Subject(cs).lambda=params_ml(1);%this shows log(lambda)
    Subject(cs).sigma=params_ml(2);%this shows log(lambda)
    Subject(cs).llh_MM=-min_mllh;
    numb_param=length(params_ml);
    numb_d=length(ans_red);
    Subject(cs).aic_MM=-2*(-min_mllh)+2*numb_param;
    Subject(cs).R2_MM=-((-min_mllh)-(numb_d*log(1/2)))/(numb_d*log(1/2));
end
end

function [mllh]=minus_log_likelyhood_T(rtn,inv,params)
lambda=params(1);
sigma=params(2);

t=[inv,1-inv];
[sigma_mat]=calc_TEntrop_prob(rtn,lambda,sigma);

numb_O=length(t(1,:));
numb_N=length(t(:,1));
llh=0;
for cn=1:numb_N
    for co=1:numb_O
        llh=llh + t(cn,co)*max(log(sigma_mat(cn,co)),-100);
    end
end
mllh=-1*llh;
end

function [sigma_mat]=calc_TEntrop_prob(rtn,lambda,sigma)
numb_N=length(rtn(:,1));
for cn=1:numb_N   
    R=rtn(cn,1)-rtn(cn,2);

    if(sigma==1)
        P=exp(rtn(cn,1)/exp(lambda))/(exp(rtn(cn,2)/exp(lambda))+exp(rtn(cn,1)/exp(lambda)));
        sigma_mat(cn,:)=[P,1-P];
    else
        P = optimvar('P',1,"LowerBound",0,"UpperBound",1);%P=probabilty_red
        expr1 = R-(exp(lambda)*(-sigma/(sigma-1))*(-((1-P)^(sigma-1))+(P^(sigma-1))));
        eqn1 = expr1 == 0;
        prob = eqnproblem;
        prob.Equations.eqn1 = eqn1;
        P0.P = [0.5];
        [sol,fval,exitflag] = solve(prob,P0);
        sigma_mat(cn,:)=[sol.P,1-sol.P];
   
    end
end
end


