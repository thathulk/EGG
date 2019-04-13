clc;clear;
runs=30;
j=1; %DTLZ1
% j=2 %DTLZ2
% j=3 %DTLZ3
% j=4 %DTLZ4
% j=36 %WFG1
% j=37 %WFG2
% j=38 %WFG3
% j=39 %WFG4
% j=40 %WFG5
% j=41 %WFG6
% j=42 %WFG7
% j=43 %WFG8
% j=44 %WFG9
idx=1;
for i=j:7:35 %This is for DTLZ problems
% for i=j:9:80 %This is for WFG problems
    R=1.1;
    if(i==5 || i==12 || i==19 || i==26 || i==33)
        continue;
    end
    if(i==6 || i==13 || i==20 || i==27 || i==34)
        continue;
    end
    if(i==7 || i==14 || i==21 || i==28 || i==35)
        continue;
    end
    if(i>=1 && i<=35)%DTLZ1-DTLZ7
        if (i>=1 && i<=7)
            problem = sprintf('DTLZ%d',i);
            m=3;
            if(i==1)
                N=91*200;
            elseif(i==3)
                N=91*500;
            else
                N=91*100;
            end
        end
        if (i>=8 && i<=14)
            problem = sprintf('DTLZ%d',i-7);
            m=5;
            if(i==1+7)
                N=210*200;
            elseif(i==3+7)
                N=210*500;
            else
                N=210*80;
            end
        end
        if (i>=15 && i<=21)
            problem = sprintf('DTLZ%d',i-14);
            m=8;
            if(i==1+7+7)
                N=156*250;
            elseif(i==3+7+7)
                N=156*750;
            else
                N=156*150;
            end
        end
        if (i>=22 && i<=28)
            problem = sprintf('DTLZ%d',i-21);
            m=10;
            if(i==1+7+7+7)
                N=275*200;
            elseif(i==3+7+7+7)
                N=275*750;
            else
                N=275*100;
            end
        end
        if (i>=29 && i<=35)
            problem = sprintf('DTLZ%d',i-28);
            m=15;
            if(i==1+7+7+7+7)
                N=135*500;
            elseif(i==3+7+7+7+7)
                N=135*2000;
            else
                N=135*250;
            end
        end
        bounds=R*ones(1,m);
        if(i==1 || i==8 || i==15 || i==22 || i==29)%DTLZ1
            bounds=R*0.5*ones(1,m);
        end
        if(i==7 || i==14 || i==21 || i==28 || i==35)
            bounds=[R*ones(1,m-1) 1.1*2.0*m];
        end
    end
    if(i>=36 && i<=80)%WFG1-WFG9
        if (i>=36 && i<=44)
            problem = sprintf('WFG%d',i-35);
            m=3;N=91*500;
        end
        if (i>=45 && i<=53)
            problem = sprintf('WFG%d',i-44);
            m=5;N=210*500;
        end
        if (i>=54 && i<=62)
            problem = sprintf('WFG%d',i-53);
            m=8;N=156*750;
        end
        if (i>=63 && i<=71)
            problem = sprintf('WFG%d',i-62);
            m=10;N=275*750;
        end
        if (i>=72 && i<=80)
            problem = sprintf('WFG%d',i-71);
            m=15;N=135*2000;
        end
        bounds = [2:2:2*m];
        bounds=R*bounds;
    end
    value=zeros(runs,6);
    for k=1:runs
        PF = load(sprintf('MOEAD_EGG/MOEAD_%s_%d_T%d',problem,m,k));
        PF1 = load(sprintf('MOEAD_SBX/MOEAD_%s_%d_T%d',problem,m,k));
        [value(k,1)] = Coverage(PF1,PF);[value(k,2)] = Coverage(PF,PF1);
        PF2 = load(sprintf('MOEAD_DE/MOEAD_%s_%d_T%d',problem,m,j));
        [value(k,3)] = Coverage(PF2,PF);[value(k,4)] = Coverage(PF,PF2);
        PF3 = load(sprintf('MOEAD_EP/%s_%dM_%dfes/FUN%d.tsv',problem,m,N,j));
        [value(k,5)] = Coverage(PF3,PF);[value(k,6)] = Coverage(PF,PF3);
    end
    h=subplot(1,5,idx);
    if(idx==1)
        h1=h;
    end
    if(idx==2)
        h2=h;
    end
    if(idx==3)
        h3=h;
    end
    if(idx==4)
        h4=h;
    end
    if(idx==5)
        h5=h;
    end
    boxplot(value);title(sprintf('%s-%d',problem,m));ylim([-0.1 1.1]);
    idx=idx+1;
end
for idx=1:5
    if(idx==1)
        set(h1,'Position',[0.035 0.13 0.16 0.75]);
    end
    if(idx==2)
        set(h2,'Position',[0.235 0.13 0.16 0.75]);
    end
    if(idx==3)
        set(h3,'Position',[0.435 0.13 0.16 0.75]);
    end
    if(idx==4)
        set(h4,'Position',[0.635 0.13 0.16 0.75]);
    end
    if(idx==5)
        set(h5,'Position',[0.83 0.13 0.16 0.75]);
    end
end
set(gcf,'Position',[300 300 800 250]);