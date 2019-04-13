clc;clear;
runs=30;
i=1;
% i = 1:4 %DTLZ1-DTLZ4 with 3 objectives
% i = 8:11 %DTLZ1-DTLZ4 with 5 objectives
% i = 15:18 %DTLZ1-DTLZ4 with 8 objectives
% i = 22:25 %DTLZ1-DTLZ4 with 10 objectives
% i = 29:32 %DTLZ1-DTLZ4 with 15 objectives
% i = 36:44 %WFG1-WFG9 with 3 objectives
% i = 45:53 %WFG1-WFG9 with 5 objectives
% i = 54:62 %WFG1-WFG9 with 8 objectives
% i = 63:71 %WFG1-WFG9 with 10 objectives
% i = 72:80 %WFG1-WFG9 with 15 objectives
R=1.1;
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
value=zeros(runs,3);
for j=1:runs
    PF = load(sprintf('thetaDEA_SBX/thetaDEA_%s_%d_T%d',problem,m,j));
    if(size(unique(PF,'rows'),1)>1)
        if m<6
            [value(j,1), nrW] = hypeIndicatorExact8(PF, bounds, m);
        else
            ftemp =  hypeIndicatorSampled( PF, bounds, 10000);
            value(j,1) = sum(ftemp);%prod(bounds);
        end
    end

    PF = load(sprintf('thetaDEA_DE/thetaDEA_%s_%d_T%d',problem,m,j));
    if(size(unique(PF,'rows'),1)>1)
        if m<6
            [value(j,2), nrW] = hypeIndicatorExact8(PF, bounds, m);
        else
            ftemp =  hypeIndicatorSampled( PF, bounds, 10000);
            value(j,2) = sum(ftemp);
        end
    end

    PF = load(sprintf('thetaDEA_EGG/thetaDEA_%s_%d_T%d',problem,m,j));
    if(size(unique(PF,'rows'),1)>1)
        if m<6
            [value(j,3), nrW] = hypeIndicatorExact8(PF, bounds, m);
        else
            ftemp =  hypeIndicatorSampled( PF, bounds, 10000);
            value(j,3) = sum(ftemp);
        end
    end
end
value = value/prod(bounds);
%     s = value(:,1);
%     save(sprintf('thetaDEA_SBX/AA_%s_%d',problem,m),'s');
%     s = value(:,2);
%     save(sprintf('thetaDEA_DE/AA_%s_%d',problem,m),'s');
%     s = value(:,3);
%     save(sprintf('thetaDEA_EGG/AA_%s_%d',problem,m),'s');

%     s = load(sprintf('thetaDEA_SBX/AA_%s_%d.mat',problem,m));
%     value(:,1) = s.s;
%     s = load(sprintf('thetaDEA_DE/AA_%s_%d.mat',problem,m));
%     value(:,2) = s.s;
%     s = load(sprintf('thetaDEA_EGG/AA_%s_%d.mat',problem,m));
%     value(:,3) = s.s;

[p1,h(1)]=ranksum(value(:,1),value(:,3));
[p2,h(2)]=ranksum(value(:,2),value(:,3));
for a=1:2
    if(h(a)==1)
        if(median(value(:,a)) < median(value(:,3)))
            sb(a)='-';
        else
            sb(a)='+';
        end
    else
        sb(a) ='~';
    end
end
med(1)=median(value(:,1));
med(2)=median(value(:,2));
med(3)=median(value(:,3));
[Val, Idx]=sort(med,'descend');
disp(sprintf('%6.5f (%3.2e) %s %d\t%6.5f (%3.2e) %s %d\t%6.5f (%3.2e) %d',median(value(:,1)),iqr(value(:,1)),sb(1),find(Idx==1),median(value(:,2)),iqr(value(:,2)),sb(2),find(Idx==2),median(value(:,3)),iqr(value(:,3)),find(Idx==3)));