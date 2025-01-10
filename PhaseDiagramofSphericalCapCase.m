%%% Phase diagram of spherical cap case
clear;clc


%% Parameter assignment
% Original parameters
Lambda=3e-9;  % Reference experimental value: 3e−9~3e−6 m/s
Lambdaw=1e-13;  % Reference experimental value: 1e−14~1e−11 m^3/(N∙s)
Jp=6e16;  % Reference experimental value: 1e17~1e18 m^-2∙s^-1
Jw=5e-10;  % Reference experimental value: 3e-10~20e-10 m∙s^-1
jp=1e6;
jpw=1.5e-19;
Sigmar=-5e-11;
Gamma=0.5e-4;  % Reference experimental value: 10^−4 N/m
Gammaj=0.4e-4;  % Reference experimental value: 10^−4 N/m
KBT=4.1e-21;
L=1.2e-5;  % cell diameter
C_cell=7.3e25;  % cell ion concentration

% Dimensionless variables
R0=L;
T0=L/Lambda;
C0=Lambda/KBT/Lambdaw;

% Dimensionless parameters
ji=Jp*KBT*Lambdaw/Lambda^2;
jw=Jw/Lambda;
jpHat=jp*KBT*Lambdaw/L^2/Lambda^2;
jpwHat=jpw/L^2/Lambda;
gamma=Gamma*Lambdaw/L/Lambda;
gammajHat=Gammaj/Gamma;
SigmarHat=Sigmar/L/Gamma;
c_cell=C_cell/C0;
dimensionless_parameters=[ji jw jpHat jpwHat gamma gammajHat SigmarHat c_cell];
% sprintf('%e  %e  %e  %e  %e  %e  %e  %e',dimensionless_parameters)



%% Phase diagram -- consider leak (directly calculate a dense grid, accelerate the convergence rate by giving a precise value range of variables)
%%% Calculate all the possible solutions
for alternative_situation=1:8  % Phase diagram的横纵坐标类型(1为j和jw，2为j和GammajHat，3为j和SigmarHat，4为j和jpHat，5为j和jpwHat，6为jpHat和jpwHat)
    gridsize=100;  % 无量纲数变化的网格数量
    Lumen_phase=zeros(gridsize,gridsize,3);  % 储存二维参数平面中管腔所处的相(对应所有可能的稳定解)
    Lumen_phase_MS=zeros(gridsize,gridsize);  % 储存二维参数平面中的最稳定的稳定解所对应的管腔所处的相(MS:most stable)
    eigenvalue_allsolution=zeros(gridsize,gridsize,3,2);  % 储存所有稳定解的特征值
    Rs_NS=zeros(gridsize,gridsize,3);
    thetas_NS=zeros(gridsize,gridsize,3);
    deltaCs_NS=zeros(gridsize,gridsize,3);

    % Phase diagram横纵坐标范围，同时寻找所选参数范围中最接近原始数值的位置
    if alternative_situation==1  % j和jw
        ji_alter=linspace(-0.5*ji,2*ji,gridsize);
        jw_alter=linspace(0,10*jw,gridsize);
        [~,I]=min(abs(ji_alter(:)-ji));x_base=I;  % x轴范围中最接近原始数值的位置的选取
        [~,I]=min(abs(jw_alter(:)-jw));y_base=I;  % y轴范围中最接近原始数值的位置的选取
    elseif alternative_situation==2  % j和GammajHat
        ji_alter=linspace(-0.5*ji,2*ji,gridsize);
        gammajHat_alter=linspace(0,2*gammajHat,gridsize);
        [~,I]=min(abs(ji_alter(:)-ji));x_base=I;  % x轴范围中最接近原始数值的位置的选取
        [~,I]=min(abs(gammajHat_alter(:)-gammajHat));y_base=I;  % y轴范围中最接近原始数值的位置的选取
    elseif alternative_situation==3  % j和SigmarHat
        ji_alter=linspace(-0.5*ji,2*ji,gridsize);
        SigmarHat_alter=linspace(0,3*SigmarHat,gridsize);
        [~,I]=min(abs(ji_alter(:)-ji));x_base=I;  % x轴范围中最接近原始数值的位置的选取
        [~,I]=min(abs(SigmarHat_alter(:)-SigmarHat));y_base=I;  % y轴范围中最接近原始数值的位置的选取
    elseif alternative_situation==4  % j和jpHat
        ji_alter=linspace(0,8*ji,gridsize);
        jpHat_alter=linspace(0,5*jpHat,gridsize);
        [~,I]=min(abs(ji_alter(:)-ji));x_base=I;  % x轴范围中最接近原始数值的位置的选取
        [~,I]=min(abs(jpHat_alter(:)-jpHat));y_base=I;  % y轴范围中最接近原始数值的位置的选取
    elseif alternative_situation==5  % j和jpwHat
        ji_alter=linspace(0,8*ji,gridsize);
        jpwHat_alter=linspace(0,5*jpwHat,gridsize);
        [~,I]=min(abs(ji_alter(:)-ji));x_base=I;  % x轴范围中最接近原始数值的位置的选取
        [~,I]=min(abs(jpwHat_alter(:)-jpwHat));y_base=I;  % y轴范围中最接近原始数值的位置的选取
    elseif alternative_situation==6  % jpHat和jpwHat
        jpHat_alter=linspace(0,5*jpHat,gridsize);
        jpwHat_alter=linspace(0,5*jpwHat,gridsize);
        [~,I]=min(abs(jpHat_alter(:)-jpHat));x_base=I;  % x轴范围中最接近原始数值的位置的选取
        [~,I]=min(abs(jpwHat_alter(:)-jpwHat));y_base=I;  % y轴范围中最接近原始数值的位置的选取
    elseif alternative_situation==7  % gammajHat和SigmarHat
        gammajHat_alter=linspace(0,2*gammajHat,gridsize);
        SigmarHat_alter=linspace(0,3*SigmarHat,gridsize);
        [~,I]=min(abs(gammajHat_alter(:)-gammajHat));x_base=I;  % x轴范围中最接近原始数值的位置的选取
        [~,I]=min(abs(SigmarHat_alter(:)-SigmarHat));y_base=I;  % y轴范围中最接近原始数值的位置的选取
    elseif alternative_situation==8  % SigmarHat和jpwHat
        SigmarHat_alter=linspace(0,3*SigmarHat,gridsize);
        jpwHat_alter=linspace(0,5*jpwHat,gridsize);
        [~,I]=min(abs(SigmarHat_alter(:)-SigmarHat));x_base=I;  % x轴范围中最接近原始数值的位置的选取
        [~,I]=min(abs(jpwHat_alter(:)-jpwHat));y_base=I;  % y轴范围中最接近原始数值的位置的选取
    end
    
    cycle_count=0;  % 记录循环位置
    for search_region=1:4  % 搜索范围以绝对稳定的原始参数数值点为中心，将相平面分成四个部分计算
        x_search_start=x_base;y_search_start=y_base;
        if search_region==1
            x_search_end=gridsize;y_search_end=gridsize;
        elseif search_region==2
            x_search_end=gridsize;y_search_end=1;
        elseif search_region==3
            x_search_end=1;y_search_end=gridsize;
        elseif search_region==4
            x_search_end=1;y_search_end=1;
        end
        
        if x_search_end>=x_search_start
            x_search_interval=1;
        else
            x_search_interval=-1;
        end
        if y_search_end>=y_search_start
            y_search_interval=1;
        else
            y_search_interval=-1;
        end
        
        Novalue_count_before=zeros(1,3);  % 记录上一列寻找结果中的没有非零稳定解的连续网格数
        for i_phase1_alter=x_search_start:x_search_interval:x_search_end
            Novalue_count=zeros(1,3);  % 记录管腔没有非零稳定解的连续网格数
            
            for i_phase2_alter=y_search_start:y_search_interval:y_search_end
                %%% alternative dimensionless parameters assignment
                if alternative_situation==1
                    ji_alter_temp=ji_alter(i_phase1_alter);  jw_alter_temp=jw_alter(i_phase2_alter);
                    gammajHat_alter_temp=gammajHat;  SigmarHat_alter_temp=SigmarHat;
                    jpHat_alter_temp=jpHat;  jpwHat_alter_temp=jpwHat;  gamma_alter_temp=gamma;
                elseif alternative_situation==2
                    ji_alter_temp=ji_alter(i_phase1_alter);  jw_alter_temp=jw;
                    gammajHat_alter_temp=gammajHat_alter(i_phase2_alter);  SigmarHat_alter_temp=SigmarHat;
                    jpHat_alter_temp=jpHat;  jpwHat_alter_temp=jpwHat;  gamma_alter_temp=gamma;
                elseif alternative_situation==3
                    ji_alter_temp=ji_alter(i_phase1_alter);  jw_alter_temp=jw;
                    gammajHat_alter_temp=gammajHat;  SigmarHat_alter_temp=SigmarHat_alter(i_phase2_alter);
                    jpHat_alter_temp=jpHat;  jpwHat_alter_temp=jpwHat;  gamma_alter_temp=gamma;
                elseif alternative_situation==4
                    ji_alter_temp=ji_alter(i_phase1_alter);  jw_alter_temp=jw;
                    gammajHat_alter_temp=gammajHat;  SigmarHat_alter_temp=SigmarHat;
                    jpHat_alter_temp=jpHat_alter(i_phase2_alter);  jpwHat_alter_temp=jpwHat;  gamma_alter_temp=gamma;
                elseif alternative_situation==5
                    ji_alter_temp=ji_alter(i_phase1_alter);  jw_alter_temp=jw;
                    gammajHat_alter_temp=gammajHat;  SigmarHat_alter_temp=SigmarHat;
                    jpHat_alter_temp=jpHat;  jpwHat_alter_temp=jpwHat_alter(i_phase2_alter);  gamma_alter_temp=gamma;
                elseif alternative_situation==6
                    ji_alter_temp=ji;  jw_alter_temp=jw;
                    gammajHat_alter_temp=gammajHat;  SigmarHat_alter_temp=SigmarHat;
                    jpHat_alter_temp=jpHat_alter(i_phase1_alter);  jpwHat_alter_temp=jpwHat_alter(i_phase2_alter);
                    gamma_alter_temp=gamma;
                elseif alternative_situation==7
                    ji_alter_temp=ji;  jw_alter_temp=jw;
                    gammajHat_alter_temp=gammajHat_alter(i_phase1_alter);  SigmarHat_alter_temp=SigmarHat_alter(i_phase2_alter);
                    jpHat_alter_temp=jpHat;  jpwHat_alter_temp=jpwHat;
                    gamma_alter_temp=gamma;
                elseif alternative_situation==8
                    ji_alter_temp=ji;  jw_alter_temp=jw;
                    gammajHat_alter_temp=gammajHat;  SigmarHat_alter_temp=SigmarHat_alter(i_phase1_alter);
                    jpHat_alter_temp=jpHat;  jpwHat_alter_temp=jpwHat_alter(i_phase2_alter);
                    gamma_alter_temp=gamma;
                end
                
                %%% Searching for all the possible steady-state solution
                syms Rs thetas deltaCs
                eqns=[2*Rs*(1-cos(thetas))*(ji_alter_temp-deltaCs)==jpHat_alter_temp*sin(thetas)/(1-Rs*sin(thetas)),...
                    2*Rs*(1-cos(thetas))*(jw_alter_temp-2*gamma_alter_temp/Rs+deltaCs)==jpwHat_alter_temp*sin(thetas)/(1-Rs*sin(thetas)),...
                    -gammajHat_alter_temp+cos(thetas)-SigmarHat_alter_temp/Rs/sin(thetas)==0];
                vars=[Rs thetas deltaCs];
  
                calculation_times=0;
                for i_solution=1:3  % 依次寻找三个可能的解
                    if (Novalue_count(i_solution)<=3 || sum(thetas_NS(i_phase1_alter,y_search_start:y_search_interval:i_phase2_alter,i_solution))<=20) && Novalue_count_before(i_solution)~=abs(y_search_start-y_search_end)+1
                        if i_phase1_alter==x_search_start && i_phase2_alter==y_search_start  % 第一次初始寻找数值解时，对搜索范围进行分类，对于容易收敛的解给予较少循环次数，对于不容易收敛的解给予较多循环次数
                            for ii=1:500
                                temp=vpasolve(eqns,vars,[0 10;0 pi;0 ji_alter_temp],'Random',true);
                                if i_solution==1 && ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                    thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                    Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                    deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                    calculation_times=calculation_times+ii;
                                    break
                                elseif i_solution==2 && ~isempty(double(temp.thetas)) && ~ismember(double(temp.thetas),thetas_NS(i_phase1_alter,i_phase2_alter,:)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                    thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                    Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                    deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                    calculation_times=calculation_times+ii;
                                    break
                                elseif i_solution==3 && ~isempty(double(temp.thetas)) && ~ismember(double(temp.thetas),thetas_NS(i_phase1_alter,i_phase2_alter,:)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                    thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                    Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                    deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                    % 对三个解按照thetas大小排序
                                    [thetas_NS(i_phase1_alter,i_phase2_alter,:),I]=sort(thetas_NS(i_phase1_alter,i_phase2_alter,:));
                                    Rs_NS(i_phase1_alter,i_phase2_alter,:)=Rs_NS(i_phase1_alter,i_phase2_alter,I);
                                    deltaCs_NS(i_phase1_alter,i_phase2_alter,:)=deltaCs_NS(i_phase1_alter,i_phase2_alter,I);
                                    calculation_times=calculation_times+ii;
                                    break
                                end
                            end
                        elseif i_phase1_alter==x_search_start && i_phase2_alter~=y_search_start  % 后续循环寻找数值解时，以上一个循环的结果为基准，在周围小范围内搜索
                            for ii=1:5
                                calculation_times=calculation_times+1;
                                temp=vpasolve(eqns,vars,[Rs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)*0.9 Rs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)*1.1; ...
                                    max(0,thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)-pi/100) min(thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)+pi/100,pi); ...
                                    deltaCs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)*0.9 deltaCs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)*1.1],'Random',true);
                                if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                    thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                    Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                    deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                    break
                                end
                            end
                            if thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)==0  % 若上一个循环的结果为零，则以其为基础的大范围中寻找该循环下的数值解
                                for ii=1:5
                                    calculation_times=calculation_times+1;
                                    temp=vpasolve(eqns,vars,[Rs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)*0.75 Rs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)*1.5; ...
                                        max(0,thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)-pi/50) min(thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)+pi/50,pi); ...
                                        deltaCs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)*0.5 deltaCs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)*2],'Random',true);
                                    if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                        thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                        Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                        deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                        break
                                    end
                                end
                            end
                            if thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)==0  % 若解仍未找到，则扩大搜索范围
                                for ii=1:5*(1-logical(Novalue_count(i_solution)))+5
                                    calculation_times=calculation_times+1;
                                    temp=vpasolve(eqns,vars,[0 1000;0 pi;0 ji_alter_temp],'Random',true);
                                    if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && ~ismember(double(temp.thetas),thetas_NS(i_phase1_alter,i_phase2_alter,:))
                                        thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                        Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                        deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                        break
                                    end
                                end
                            end
                            if thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)==0  % 若上一个解的deltaCs接近0，则将deltaCs的搜索范围扩大到负值
                                if abs(deltaCs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution))<5e-2
                                    for ii=1:5
                                        calculation_times=calculation_times+1;
                                        temp=vpasolve(eqns,vars,[Rs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)*0.75 Rs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)*1.5; ...
                                            max(0,thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)-pi/50) min(thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)+pi/50,pi); ...
                                            -2 abs(deltaCs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution))*2],'Random',true);
                                        if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                            thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                            Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                            deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                            break
                                        end
                                    end
                                end
                            end
                        elseif i_phase1_alter~=x_search_start
                            logic_before=0;  % 存储是否依据前一个网格结果作为基准的逻辑值
                            if i_phase2_alter-y_search_interval==0 || i_phase2_alter-y_search_interval>gridsize
                                if thetas_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution)~=0
                                    Rs_NS_before=Rs_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution);
                                    thetas_NS_before=thetas_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution);
                                    deltaCs_NS_before=deltaCs_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution);
                                    logic_before=1;
                                end
                            else
                                if thetas_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution)<=thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution) && thetas_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution)~=0
                                    Rs_NS_before=Rs_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution);
                                    thetas_NS_before=thetas_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution);
                                    deltaCs_NS_before=deltaCs_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution);
                                    logic_before=1;
                                elseif thetas_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution)<=thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution) && thetas_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution)==0 && thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)~=0
                                    Rs_NS_before=Rs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution);
                                    thetas_NS_before=thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution);
                                    deltaCs_NS_before=deltaCs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution);
                                    logic_before=1;
                                elseif thetas_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution)>thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution) && thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)~=0
                                    Rs_NS_before=Rs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution);
                                    thetas_NS_before=thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution);
                                    deltaCs_NS_before=deltaCs_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution);
                                    logic_before=1;
                                elseif thetas_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution)>thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution) && thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,i_solution)==0 && thetas_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution)~=0
                                    Rs_NS_before=Rs_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution);
                                    thetas_NS_before=thetas_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution);
                                    deltaCs_NS_before=deltaCs_NS(i_phase1_alter-x_search_interval,i_phase2_alter,i_solution);
                                    logic_before=1;
                                end
                            end
                            
                            if logic_before==1
                                if thetas_NS_before<pi*3/4 && thetas_NS_before>1e-5  % 若上一个循环的结果不为零，且解相对较小(倾向于是稳定解)，则以其为基础的小范围中寻找该循环下的数值解
                                    for ii=1:5
                                        calculation_times=calculation_times+1;
                                        temp=vpasolve(eqns,vars,[Rs_NS_before*0.9 Rs_NS_before*1.1; ...
                                            max(0,thetas_NS_before-pi/100) min(thetas_NS_before+pi/100,pi); ...
                                            deltaCs_NS_before*0.9 deltaCs_NS_before*1.1],'Random',true);
                                        if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && ~ismember(double(temp.thetas),thetas_NS(i_phase1_alter,i_phase2_alter,:))
                                            thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                            Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                            deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                            break
                                        end
                                    end
                                elseif thetas_NS_before>=pi*3/4 || thetas_NS_before<1e-5  % 若上一个循环的结果不为零，且解相对较大(倾向于是不稳定解)，则以其为基础的大范围中寻找该循环下的数值解
                                    for ii=1:5
                                        calculation_times=calculation_times+1;
                                        temp=vpasolve(eqns,vars,[Rs_NS_before*0.75 Rs_NS_before*1.5; ...
                                            max(0,thetas_NS_before-pi/50) min(thetas_NS_before+pi/50,pi); ...
                                            deltaCs_NS_before*0.5 deltaCs_NS_before*2],'Random',true);
                                        if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && ~ismember(double(temp.thetas),thetas_NS(i_phase1_alter,i_phase2_alter,:))
                                            thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                            Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                            deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                            break
                                        end
                                    end
                                end
                                if thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)==0 && thetas_NS_before<pi*3/4 && thetas_NS_before>1e-5  % 若上一个循环的结果为零，则以其为基础的小范围中寻找该循环下的数值解
                                    for ii=1:5
                                        calculation_times=calculation_times+1;
                                        temp=vpasolve(eqns,vars,[Rs_NS_before*0.75 Rs_NS_before*1.5; ...
                                            max(0,thetas_NS_before-pi/50) min(thetas_NS_before+pi/50,pi); ...
                                            deltaCs_NS_before*0.5 deltaCs_NS_before*2],'Random',true);
                                        if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && ~ismember(double(temp.thetas),thetas_NS(i_phase1_alter,i_phase2_alter,:))
                                            thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                            Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                            deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                            break
                                        end
                                    end
                                elseif thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)==0 && thetas_NS_before>=pi*3/4 || thetas_NS_before<1e-5  % 若上一个循环的结果不为零，且解相对较大(倾向于是不稳定解)，则以其为基础的大范围中寻找该循环下的数值解
                                    for ii=1:5
                                        calculation_times=calculation_times+1;
                                        temp=vpasolve(eqns,vars,[Rs_NS_before*0.5 Rs_NS_before*2; ...
                                            max(0,thetas_NS_before-pi/50) min(thetas_NS_before+pi/50,pi); ...
                                            deltaCs_NS_before*0.3 deltaCs_NS_before*3],'Random',true);
                                        if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && ~ismember(double(temp.thetas),thetas_NS(i_phase1_alter,i_phase2_alter,:))
                                            thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                            Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                            deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                            break
                                        end
                                    end
                                end
                            end
                            if thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)==0  % 若容易收敛的解未找到，则扩大搜索范围
                                for ii=1:5*(1-logical(Novalue_count(i_solution)))+5
                                    calculation_times=calculation_times+1;
                                    temp=vpasolve(eqns,vars,[0 1000;0 pi;0 ji_alter_temp],'Random',true);
                                    if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && ~ismember(double(temp.thetas),thetas_NS(i_phase1_alter,i_phase2_alter,:))
                                        thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                        Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                        deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                        break
                                    end
                                end
                            end
                            if thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)==0  % 若上一个解的deltaCs接近0，则将deltaCs的搜索范围扩大到负值
                                if abs(deltaCs_NS_before)<5e-2
                                    for ii=1:5
                                        calculation_times=calculation_times+1;
                                        temp=vpasolve(eqns,vars,[Rs_NS_before*0.5 Rs_NS_before*2; ...
                                            max(0,thetas_NS_before-pi/50) min(thetas_NS_before+pi/50,pi); ...
                                            -2 abs(deltaCs_NS_before)*3],'Random',true);
                                        if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && ~ismember(double(temp.thetas),thetas_NS(i_phase1_alter,i_phase2_alter,:))
                                            thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.thetas);
                                            Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.Rs);
                                            deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution)=double(temp.deltaCs);
                                            break
                                        end
                                    end
                                end
                            end
                        end
                        
                        %%% Linear stability analysis
                        if Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)~=0
                            Rs_best=Rs_NS(i_phase1_alter,i_phase2_alter,i_solution);
                            thetas_best=thetas_NS(i_phase1_alter,i_phase2_alter,i_solution);
                            deltaCs_best=deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution);

                            Mstab11=(4*jw_alter_temp*Rs_best-4*gamma_alter_temp+4*Rs_best*deltaCs_best-4*(jw_alter_temp*Rs_best-gamma_alter_temp+Rs_best* ...
                                deltaCs_best)*cos(thetas_best)-(jpwHat_alter_temp*sin(thetas_best))/(-1+Rs_best*sin(thetas_best))^2+(SigmarHat_alter_temp* ...
                                csc(thetas_best)*(2*jw_alter_temp*Rs_best-4*gamma_alter_temp+2*Rs_best*deltaCs_best-(jpwHat_alter_temp*cot(thetas_best))/ ...
                                (-1+Rs_best*sin(thetas_best))^2))/(Rs_best-SigmarHat_alter_temp*cot(thetas_best)*csc(thetas_best)^2))/(Rs_best^2*(2-3* ...
                                cos(thetas_best)+cos(thetas_best)^3+(SigmarHat_alter_temp*sin(thetas_best))/(Rs_best-SigmarHat_alter_temp*cot(thetas_best)* ...
                                csc(thetas_best)^2)));
                            Mstab12=(4*sin(thetas_best/2)^2)/(2-3*cos(thetas_best)+cos(thetas_best)^3+(SigmarHat_alter_temp*sin(thetas_best))/(Rs_best- ...
                                SigmarHat_alter_temp*cot(thetas_best)*csc(thetas_best)^2));
                            Mstab21=(3*csc(thetas_best/2)^5*sec(thetas_best/2)^8*((-5*Rs_best^2+64*SigmarHat_alter_temp^2)*cos(thetas_best)+Rs_best*(-40* ...
                                Rs_best+60*Rs_best*cos(2*thetas_best)+9*Rs_best*cos(3*thetas_best)-24*Rs_best*cos(4*thetas_best)-5*Rs_best*cos(5*thetas_best)+ ...
                                4*Rs_best*cos(6*thetas_best)+Rs_best*cos(7*thetas_best)-40*SigmarHat_alter_temp*sin(thetas_best)+32*SigmarHat_alter_temp* ...
                                sin(2*thetas_best)+20*SigmarHat_alter_temp*sin(3*thetas_best)-16*SigmarHat_alter_temp*sin(4*thetas_best)-4*SigmarHat_alter_temp* ...
                                sin(5*thetas_best)))*(-2*(jpHat_alter_temp*(39*Rs_best^2-48*SigmarHat_alter_temp)+4*Rs_best*(ji_alter_temp*Rs_best* ...
                                (6+5*Rs_best^2-10*SigmarHat_alter_temp)+2*c_cell*(jw_alter_temp*Rs_best*(6+5*Rs_best^2-7*SigmarHat_alter_temp)+Rs_best* ...
                                deltaCs_best*(6+5*Rs_best^2-7*SigmarHat_alter_temp)+gamma_alter_temp*(-6-5*Rs_best^2+16*SigmarHat_alter_temp))+deltaCs_best*(-10* ...
                                Rs_best^2*gamma_alter_temp+5*Rs_best^3*(-1+2*jw_alter_temp+2*deltaCs_best)+4*gamma_alter_temp*(-3+8*SigmarHat_alter_temp)-2*Rs_best* ...
                                (3-6*deltaCs_best-5*SigmarHat_alter_temp+7*deltaCs_best*SigmarHat_alter_temp+jw_alter_temp*(-6+7*SigmarHat_alter_temp)))))* ...
                                cos(thetas_best/2)+Rs_best*(42*jpHat_alter_temp*Rs_best+80*c_cell*jw_alter_temp*Rs_best+70*c_cell*jw_alter_temp*Rs_best^3- ...
                                80*c_cell*gamma_alter_temp-70*c_cell*Rs_best^2*gamma_alter_temp-40*Rs_best*deltaCs_best+80*c_cell*Rs_best*deltaCs_best+80*jw_alter_temp* ...
                                Rs_best*deltaCs_best-35*Rs_best^3*deltaCs_best+70*c_cell*Rs_best^3*deltaCs_best+70*jw_alter_temp*Rs_best^3*deltaCs_best-80* ...
                                gamma_alter_temp*deltaCs_best-70*Rs_best^2*gamma_alter_temp*deltaCs_best+80*Rs_best*deltaCs_best^2+70*Rs_best^3*deltaCs_best^2+ ...
                                ji_alter_temp*Rs_best*(40+35*Rs_best^2-48*SigmarHat_alter_temp)-144*c_cell*jw_alter_temp*Rs_best*SigmarHat_alter_temp+192*c_cell* ...
                                gamma_alter_temp*SigmarHat_alter_temp+48*Rs_best*deltaCs_best*SigmarHat_alter_temp-144*c_cell*Rs_best*deltaCs_best*SigmarHat_alter_temp- ...
                                144*jw_alter_temp*Rs_best*deltaCs_best*SigmarHat_alter_temp+192*gamma_alter_temp*deltaCs_best*SigmarHat_alter_temp-144*Rs_best* ...
                                deltaCs_best^2*SigmarHat_alter_temp)*cos((3*thetas_best)/2)+57*jpHat_alter_temp*Rs_best^2*cos((5*thetas_best)/2)+ ...
                                24*ji_alter_temp*Rs_best^2*cos((5*thetas_best)/2)+48*c_cell*jw_alter_temp*Rs_best^2*cos((5*thetas_best)/2)+25*ji_alter_temp* ...
                                Rs_best^4*cos((5*thetas_best)/2)+50*c_cell*jw_alter_temp*Rs_best^4*cos((5*thetas_best)/2)-48*c_cell*Rs_best*gamma_alter_temp* ...
                                cos((5*thetas_best)/2)-50*c_cell*Rs_best^3*gamma_alter_temp*cos((5*thetas_best)/2)-24*Rs_best^2*deltaCs_best*cos((5*thetas_best)/2)+ ...
                                48*c_cell*Rs_best^2*deltaCs_best*cos((5*thetas_best)/2)+48*jw_alter_temp*Rs_best^2*deltaCs_best*cos((5*thetas_best)/2)- ...
                                25*Rs_best^4*deltaCs_best*cos((5*thetas_best)/2)+50*c_cell*Rs_best^4*deltaCs_best*cos((5*thetas_best)/2)+50*jw_alter_temp* ...
                                Rs_best^4*deltaCs_best*cos((5*thetas_best)/2)-48*Rs_best*gamma_alter_temp*deltaCs_best*cos((5*thetas_best)/2)-50*Rs_best^3* ...
                                gamma_alter_temp*deltaCs_best*cos((5*thetas_best)/2)+48*Rs_best^2*deltaCs_best^2*cos((5*thetas_best)/2)+50*Rs_best^4*deltaCs_best^2* ...
                                cos((5*thetas_best)/2)-16*ji_alter_temp*Rs_best^2*SigmarHat_alter_temp*cos((5*thetas_best)/2)+16*c_cell*jw_alter_temp* ...
                                Rs_best^2*SigmarHat_alter_temp*cos((5*thetas_best)/2)+64*c_cell*Rs_best*gamma_alter_temp*SigmarHat_alter_temp*cos((5*thetas_best)/2)+ ...
                                16*Rs_best^2*deltaCs_best*SigmarHat_alter_temp*cos((5*thetas_best)/2)+16*c_cell*Rs_best^2*deltaCs_best*SigmarHat_alter_temp* ...
                                cos((5*thetas_best)/2)+16*jw_alter_temp*Rs_best^2*deltaCs_best*SigmarHat_alter_temp*cos((5*thetas_best)/2)+64*Rs_best* ...
                                gamma_alter_temp*deltaCs_best*SigmarHat_alter_temp*cos((5*thetas_best)/2)+16*Rs_best^2*deltaCs_best^2*SigmarHat_alter_temp* ...
                                cos((5*thetas_best)/2)-3*jpHat_alter_temp*Rs_best^2*cos((7*thetas_best)/2)-12*ji_alter_temp*Rs_best^2*cos((7*thetas_best)/2)- ...
                                24*c_cell*jw_alter_temp*Rs_best^2*cos((7*thetas_best)/2)-16*ji_alter_temp*Rs_best^4*cos((7*thetas_best)/2)- ...
                                32*c_cell*jw_alter_temp*Rs_best^4*cos((7*thetas_best)/2)+24*c_cell*Rs_best*gamma_alter_temp*cos((7*thetas_best)/2)+ ...
                                32*c_cell*Rs_best^3*gamma_alter_temp*cos((7*thetas_best)/2)+12*Rs_best^2*deltaCs_best*cos((7*thetas_best)/2)- ...
                                24*c_cell*Rs_best^2*deltaCs_best*cos((7*thetas_best)/2)-24*jw_alter_temp*Rs_best^2*deltaCs_best*cos((7*thetas_best)/2)+ ...
                                16*Rs_best^4*deltaCs_best*cos((7*thetas_best)/2)-32*c_cell*Rs_best^4*deltaCs_best*cos((7*thetas_best)/2)-32*jw_alter_temp* ...
                                Rs_best^4*deltaCs_best*cos((7*thetas_best)/2)+24*Rs_best*gamma_alter_temp*deltaCs_best*cos((7*thetas_best)/2)+32*Rs_best^3* ...
                                gamma_alter_temp*deltaCs_best*cos((7*thetas_best)/2)-24*Rs_best^2*deltaCs_best^2*cos((7*thetas_best)/2)-32*Rs_best^4*deltaCs_best^2* ...
                                cos((7*thetas_best)/2)-16*ji_alter_temp*Rs_best^2*SigmarHat_alter_temp*cos((7*thetas_best)/2)+16*c_cell*jw_alter_temp* ...
                                Rs_best^2*SigmarHat_alter_temp*cos((7*thetas_best)/2)+16*Rs_best^2*deltaCs_best*SigmarHat_alter_temp*cos((7*thetas_best)/2)+ ...
                                16*c_cell*Rs_best^2*deltaCs_best*SigmarHat_alter_temp*cos((7*thetas_best)/2)+16*jw_alter_temp*Rs_best^2*deltaCs_best* ...
                                SigmarHat_alter_temp*cos((7*thetas_best)/2)+16*Rs_best^2*deltaCs_best^2*SigmarHat_alter_temp*cos((7*thetas_best)/2)- ...
                                15*jpHat_alter_temp*Rs_best^2*cos((9*thetas_best)/2)-4*ji_alter_temp*Rs_best^2*cos((9*thetas_best)/2)-8*c_cell*jw_alter_temp* ...
                                Rs_best^2*cos((9*thetas_best)/2)-8*ji_alter_temp*Rs_best^4*cos((9*thetas_best)/2)-16*c_cell*jw_alter_temp*Rs_best^4* ...
                                cos((9*thetas_best)/2)+8*c_cell*Rs_best*gamma_alter_temp*cos((9*thetas_best)/2)+16*c_cell*Rs_best^3*gamma_alter_temp* ...
                                cos((9*thetas_best)/2)+4*Rs_best^2*deltaCs_best*cos((9*thetas_best)/2)-8*c_cell*Rs_best^2*deltaCs_best* ...
                                cos((9*thetas_best)/2)-8*jw_alter_temp*Rs_best^2*deltaCs_best*cos((9*thetas_best)/2)+8*Rs_best^4*deltaCs_best* ...
                                cos((9*thetas_best)/2)-16*c_cell*Rs_best^4*deltaCs_best*cos((9*thetas_best)/2)-16*jw_alter_temp*Rs_best^4* ...
                                deltaCs_best*cos((9*thetas_best)/2)+8*Rs_best*gamma_alter_temp*deltaCs_best*cos((9*thetas_best)/2)+16*Rs_best^3*gamma_alter_temp* ...
                                deltaCs_best*cos((9*thetas_best)/2)-8*Rs_best^2*deltaCs_best^2*cos((9*thetas_best)/2)-16*Rs_best^4*deltaCs_best^2* ...
                                cos((9*thetas_best)/2)-3*jpHat_alter_temp*Rs_best^2*cos((11*thetas_best)/2)+3*ji_alter_temp*Rs_best^4* ...
                                cos((11*thetas_best)/2)+6*c_cell*jw_alter_temp*Rs_best^4*cos((11*thetas_best)/2)-6*c_cell*Rs_best^3*gamma_alter_temp* ...
                                cos((11*thetas_best)/2)-3*Rs_best^4*deltaCs_best*cos((11*thetas_best)/2)+6*c_cell*Rs_best^4*deltaCs_best* ...
                                cos((11*thetas_best)/2)+6*jw_alter_temp*Rs_best^4*deltaCs_best*cos((11*thetas_best)/2)-6*Rs_best^3*gamma_alter_temp* ...
                                deltaCs_best*cos((11*thetas_best)/2)+6*Rs_best^4*deltaCs_best^2*cos((11*thetas_best)/2)+ji_alter_temp*Rs_best^4* ...
                                cos((13*thetas_best)/2)+2*c_cell*jw_alter_temp*Rs_best^4*cos((13*thetas_best)/2)-2*c_cell*Rs_best^3*gamma_alter_temp* ...
                                cos((13*thetas_best)/2)-Rs_best^4*deltaCs_best*cos((13*thetas_best)/2)+2*c_cell*Rs_best^4*deltaCs_best* ...
                                cos((13*thetas_best)/2)+2*jw_alter_temp*Rs_best^4*deltaCs_best*cos((13*thetas_best)/2)-2*Rs_best^3*gamma_alter_temp* ...
                                deltaCs_best*cos((13*thetas_best)/2)+2*Rs_best^4*deltaCs_best^2*cos((13*thetas_best)/2)+48*jpHat_alter_temp* ...
                                Rs_best*sin(thetas_best/2)+24*c_cell*jpwHat_alter_temp*Rs_best*sin(thetas_best/2)+88*ji_alter_temp*Rs_best^3* ...
                                sin(thetas_best/2)+176*c_cell*jw_alter_temp*Rs_best^3*sin(thetas_best/2)-176*c_cell*Rs_best^2*gamma_alter_temp* ...
                                sin(thetas_best/2)+24*jpwHat_alter_temp*Rs_best*deltaCs_best*sin(thetas_best/2)-88*Rs_best^3*deltaCs_best* ...
                                sin(thetas_best/2)+176*c_cell*Rs_best^3*deltaCs_best*sin(thetas_best/2)+176*jw_alter_temp*Rs_best^3*deltaCs_best* ...
                                sin(thetas_best/2)-176*Rs_best^2*gamma_alter_temp*deltaCs_best*sin(thetas_best/2)+176*Rs_best^3*deltaCs_best^2* ...
                                sin(thetas_best/2)-48*jpHat_alter_temp*Rs_best*SigmarHat_alter_temp*sin(thetas_best/2)-64*ji_alter_temp*Rs_best* ...
                                SigmarHat_alter_temp*sin(thetas_best/2)-128*c_cell*jw_alter_temp*Rs_best*SigmarHat_alter_temp*sin(thetas_best/2)- ...
                                32*ji_alter_temp*Rs_best^3*SigmarHat_alter_temp*sin(thetas_best/2)-64*c_cell*jw_alter_temp*Rs_best^3*SigmarHat_alter_temp* ...
                                sin(thetas_best/2)+192*c_cell*gamma_alter_temp*SigmarHat_alter_temp*sin(thetas_best/2)+112*c_cell*Rs_best^2* ...
                                gamma_alter_temp*SigmarHat_alter_temp*sin(thetas_best/2)+64*Rs_best*deltaCs_best*SigmarHat_alter_temp*sin(thetas_best/2)- ...
                                128*c_cell*Rs_best*deltaCs_best*SigmarHat_alter_temp*sin(thetas_best/2)-128*jw_alter_temp*Rs_best*deltaCs_best* ...
                                SigmarHat_alter_temp*sin(thetas_best/2)+32*Rs_best^3*deltaCs_best*SigmarHat_alter_temp*sin(thetas_best/2)-64*c_cell* ...
                                Rs_best^3*deltaCs_best*SigmarHat_alter_temp*sin(thetas_best/2)-64*jw_alter_temp*Rs_best^3*deltaCs_best*SigmarHat_alter_temp* ...
                                sin(thetas_best/2)+192*gamma_alter_temp*deltaCs_best*SigmarHat_alter_temp*sin(thetas_best/2)+112*Rs_best^2* ...
                                gamma_alter_temp*deltaCs_best*SigmarHat_alter_temp*sin(thetas_best/2)-128*Rs_best*deltaCs_best^2*SigmarHat_alter_temp* ...
                                sin(thetas_best/2)-64*Rs_best^3*deltaCs_best^2*SigmarHat_alter_temp*sin(thetas_best/2)+56*jpHat_alter_temp* ...
                                Rs_best*sin((3*thetas_best)/2)+28*c_cell*jpwHat_alter_temp*Rs_best*sin((3*thetas_best)/2)+72*ji_alter_temp* ...
                                Rs_best^3*sin((3*thetas_best)/2)+144*c_cell*jw_alter_temp*Rs_best^3*sin((3*thetas_best)/2)-144*c_cell*Rs_best^2* ...
                                gamma_alter_temp*sin((3*thetas_best)/2)+28*jpwHat_alter_temp*Rs_best*deltaCs_best*sin((3*thetas_best)/2)-72*Rs_best^3* ...
                                deltaCs_best*sin((3*thetas_best)/2)+144*c_cell*Rs_best^3*deltaCs_best*sin((3*thetas_best)/2)+144*jw_alter_temp*Rs_best^3* ...
                                deltaCs_best*sin((3*thetas_best)/2)-144*Rs_best^2*gamma_alter_temp*deltaCs_best*sin((3*thetas_best)/2)+144*Rs_best^3*deltaCs_best^2* ...
                                sin((3*thetas_best)/2)-48*jpHat_alter_temp*Rs_best*SigmarHat_alter_temp*sin((3*thetas_best)/2)-16*ji_alter_temp* ...
                                Rs_best*SigmarHat_alter_temp*sin((3*thetas_best)/2)+16*c_cell*jw_alter_temp*Rs_best*SigmarHat_alter_temp*sin((3*thetas_best)/2)- ...
                                24*ji_alter_temp*Rs_best^3*SigmarHat_alter_temp*sin((3*thetas_best)/2)-24*c_cell*jw_alter_temp*Rs_best^3*SigmarHat_alter_temp* ...
                                sin((3*thetas_best)/2)+64*c_cell*gamma_alter_temp*SigmarHat_alter_temp*sin((3*thetas_best)/2)+80*c_cell*Rs_best^2*gamma_alter_temp* ...
                                SigmarHat_alter_temp*sin((3*thetas_best)/2)+16*Rs_best*deltaCs_best*SigmarHat_alter_temp*sin((3*thetas_best)/2)+16*c_cell* ...
                                Rs_best*deltaCs_best*SigmarHat_alter_temp*sin((3*thetas_best)/2)+16*jw_alter_temp*Rs_best*deltaCs_best*SigmarHat_alter_temp* ...
                                sin((3*thetas_best)/2)+24*Rs_best^3*deltaCs_best*SigmarHat_alter_temp*sin((3*thetas_best)/2)-24*c_cell*Rs_best^3* ...
                                deltaCs_best*SigmarHat_alter_temp*sin((3*thetas_best)/2)-24*jw_alter_temp*Rs_best^3*deltaCs_best*SigmarHat_alter_temp* ...
                                sin((3*thetas_best)/2)+64*gamma_alter_temp*deltaCs_best*SigmarHat_alter_temp*sin((3*thetas_best)/2)+80*Rs_best^2*gamma_alter_temp* ...
                                deltaCs_best*SigmarHat_alter_temp*sin((3*thetas_best)/2)+16*Rs_best*deltaCs_best^2*SigmarHat_alter_temp*sin((3*thetas_best)/2)- ...
                                24*Rs_best^3*deltaCs_best^2*SigmarHat_alter_temp*sin((3*thetas_best)/2)-8*jpHat_alter_temp*Rs_best* ...
                                sin((5*thetas_best)/2)-4*c_cell*jpwHat_alter_temp*Rs_best*sin((5*thetas_best)/2)-52*ji_alter_temp* ...
                                Rs_best^3*sin((5*thetas_best)/2)-104*c_cell*jw_alter_temp*Rs_best^3*sin((5*thetas_best)/2)+104*c_cell*Rs_best^2* ...
                                gamma_alter_temp*sin((5*thetas_best)/2)-4*jpwHat_alter_temp*Rs_best*deltaCs_best*sin((5*thetas_best)/2)+52*Rs_best^3* ...
                                deltaCs_best*sin((5*thetas_best)/2)-104*c_cell*Rs_best^3*deltaCs_best*sin((5*thetas_best)/2)-104*jw_alter_temp*Rs_best^3* ...
                                deltaCs_best*sin((5*thetas_best)/2)+104*Rs_best^2*gamma_alter_temp*deltaCs_best*sin((5*thetas_best)/2)-104*Rs_best^3*deltaCs_best^2* ...
                                sin((5*thetas_best)/2)-16*ji_alter_temp*Rs_best*SigmarHat_alter_temp*sin((5*thetas_best)/2)+16*c_cell*jw_alter_temp* ...
                                Rs_best*SigmarHat_alter_temp*sin((5*thetas_best)/2)+8*ji_alter_temp*Rs_best^3*SigmarHat_alter_temp*sin((5*thetas_best)/2)+ ...
                                40*c_cell*jw_alter_temp*Rs_best^3*SigmarHat_alter_temp*sin((5*thetas_best)/2)-48*c_cell*Rs_best^2*gamma_alter_temp*SigmarHat_alter_temp* ...
                                sin((5*thetas_best)/2)+16*Rs_best*deltaCs_best*SigmarHat_alter_temp*sin((5*thetas_best)/2)+16*c_cell*Rs_best*deltaCs_best* ...
                                SigmarHat_alter_temp*sin((5*thetas_best)/2)+16*jw_alter_temp*Rs_best*deltaCs_best*SigmarHat_alter_temp*sin((5*thetas_best)/2)- ...
                                8*Rs_best^3*deltaCs_best*SigmarHat_alter_temp*sin((5*thetas_best)/2)+40*c_cell*Rs_best^3*deltaCs_best*SigmarHat_alter_temp* ...
                                sin((5*thetas_best)/2)+40*jw_alter_temp*Rs_best^3*deltaCs_best*SigmarHat_alter_temp*sin((5*thetas_best)/2)-48*Rs_best^2* ...
                                gamma_alter_temp*deltaCs_best*SigmarHat_alter_temp*sin((5*thetas_best)/2)+16*Rs_best*deltaCs_best^2*SigmarHat_alter_temp* ...
                                sin((5*thetas_best)/2)+40*Rs_best^3*deltaCs_best^2*SigmarHat_alter_temp*sin((5*thetas_best)/2)-20*jpHat_alter_temp* ...
                                Rs_best*sin((7*thetas_best)/2)-10*c_cell*jpwHat_alter_temp*Rs_best*sin((7*thetas_best)/2)-28*ji_alter_temp* ...
                                Rs_best^3*sin((7*thetas_best)/2)-56*c_cell*jw_alter_temp*Rs_best^3*sin((7*thetas_best)/2)+56*c_cell*Rs_best^2*gamma_alter_temp* ...
                                sin((7*thetas_best)/2)-10*jpwHat_alter_temp*Rs_best*deltaCs_best*sin((7*thetas_best)/2)+28*Rs_best^3* ...
                                deltaCs_best*sin((7*thetas_best)/2)-56*c_cell*Rs_best^3*deltaCs_best*sin((7*thetas_best)/2)-56*jw_alter_temp*Rs_best^3* ...
                                deltaCs_best*sin((7*thetas_best)/2)+56*Rs_best^2*gamma_alter_temp*deltaCs_best*sin((7*thetas_best)/2)-56*Rs_best^3*deltaCs_best^2* ...
                                sin((7*thetas_best)/2)+4*ji_alter_temp*Rs_best^3*SigmarHat_alter_temp*sin((7*thetas_best)/2)-4*c_cell*jw_alter_temp* ...
                                Rs_best^3*SigmarHat_alter_temp*sin((7*thetas_best)/2)-16*c_cell*Rs_best^2*gamma_alter_temp*SigmarHat_alter_temp*sin((7*thetas_best)/2)- ...
                                4*Rs_best^3*deltaCs_best*SigmarHat_alter_temp*sin((7*thetas_best)/2)-4*c_cell*Rs_best^3*deltaCs_best*SigmarHat_alter_temp* ...
                                sin((7*thetas_best)/2)-4*jw_alter_temp*Rs_best^3*deltaCs_best*SigmarHat_alter_temp*sin((7*thetas_best)/2)-16*Rs_best^2* ...
                                gamma_alter_temp*deltaCs_best*SigmarHat_alter_temp*sin((7*thetas_best)/2)-4*Rs_best^3*deltaCs_best^2*SigmarHat_alter_temp* ...
                                sin((7*thetas_best)/2)-4*jpHat_alter_temp*Rs_best*sin((9*thetas_best)/2)-2*c_cell*jpwHat_alter_temp* ...
                                Rs_best*sin((9*thetas_best)/2)+12*ji_alter_temp*Rs_best^3*sin((9*thetas_best)/2)+24*c_cell*jw_alter_temp*Rs_best^3* ...
                                sin((9*thetas_best)/2)-24*c_cell*Rs_best^2*gamma_alter_temp*sin((9*thetas_best)/2)-2*jpwHat_alter_temp*Rs_best*deltaCs_best* ...
                                sin((9*thetas_best)/2)-12*Rs_best^3*deltaCs_best*sin((9*thetas_best)/2)+24*c_cell*Rs_best^3*deltaCs_best* ...
                                sin((9*thetas_best)/2)+24*jw_alter_temp*Rs_best^3*deltaCs_best*sin((9*thetas_best)/2)-24*Rs_best^2*gamma_alter_temp*deltaCs_best* ...
                                sin((9*thetas_best)/2)+24*Rs_best^3*deltaCs_best^2*sin((9*thetas_best)/2)+4*ji_alter_temp*Rs_best^3*SigmarHat_alter_temp* ...
                                sin((9*thetas_best)/2)-4*c_cell*jw_alter_temp*Rs_best^3*SigmarHat_alter_temp*sin((9*thetas_best)/2)-4*Rs_best^3* ...
                                deltaCs_best*SigmarHat_alter_temp*sin((9*thetas_best)/2)-4*c_cell*Rs_best^3*deltaCs_best*SigmarHat_alter_temp* ...
                                sin((9*thetas_best)/2)-4*jw_alter_temp*Rs_best^3*deltaCs_best*SigmarHat_alter_temp*sin((9*thetas_best)/2)-4*Rs_best^3* ...
                                deltaCs_best^2*SigmarHat_alter_temp*sin((9*thetas_best)/2)+4*ji_alter_temp*Rs_best^3*sin((11*thetas_best)/2)+8*c_cell* ...
                                jw_alter_temp*Rs_best^3*sin((11*thetas_best)/2)-8*c_cell*Rs_best^2*gamma_alter_temp*sin((11*thetas_best)/2)-4*Rs_best^3* ...
                                deltaCs_best*sin((11*thetas_best)/2)+8*c_cell*Rs_best^3*deltaCs_best*sin((11*thetas_best)/2)+8*jw_alter_temp*Rs_best^3* ...
                                deltaCs_best*sin((11*thetas_best)/2)-8*Rs_best^2*gamma_alter_temp*deltaCs_best*sin((11*thetas_best)/2)+8*Rs_best^3*deltaCs_best^2* ...
                                sin((11*thetas_best)/2)))/(16384*Rs_best^3*(-1+cos(thetas_best))^3*(2+cos(thetas_best))^2*(Rs_best-csc(thetas_best))^2* ...
                                (Rs_best-SigmarHat_alter_temp*cot(thetas_best)*csc(thetas_best)^2)^2*(8*SigmarHat_alter_temp+12*Rs_best*sin(thetas_best)+ ...
                                2*Rs_best*sin(2*thetas_best)-4*Rs_best*sin(3*thetas_best)-Rs_best*sin(4*thetas_best)));
                            Mstab22=-((3*(1+c_cell+deltaCs_best)*csc(thetas_best/2)^2)/(Rs_best*(2+cos(thetas_best))));
                            Mstab=[Mstab11 Mstab12;Mstab21 Mstab22];
                            eigenvalue_allsolution(i_phase1_alter,i_phase2_alter,i_solution,:)=eig(Mstab);
                        end
                        
                        %%% 判断lumen所处的相
                        if Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)~=0
                            if isreal(eigenvalue_allsolution(i_phase1_alter,i_phase2_alter,i_solution,:))
                                if eigenvalue_allsolution(i_phase1_alter,i_phase2_alter,i_solution,1)<0 && eigenvalue_allsolution(i_phase1_alter,i_phase2_alter,i_solution,2)<0
                                    Lumen_phase(i_phase1_alter,i_phase2_alter,i_solution)=1;  % Lumen is stable
                                else
                                    Lumen_phase(i_phase1_alter,i_phase2_alter,i_solution)=2;  % Lumen is unstable
                                end
                            elseif ~isreal(eigenvalue_allsolution(i_phase1_alter,i_phase2_alter,i_solution,:)) && sum(real(eigenvalue_allsolution(i_phase1_alter,i_phase2_alter,i_solution,:))<0)==2  % 两个特征值的实部均小于零，且存在虚部
                                Lumen_phase(i_phase1_alter,i_phase2_alter,i_solution)=3;  % Lumen is oscillated
                            end
                            cycle_count=cycle_count+logical(i_solution==3);
                            Novalue_count(i_solution)=0;
                        elseif Rs_NS(i_phase1_alter,i_phase2_alter,i_solution)==0
                            Lumen_phase(i_phase1_alter,i_phase2_alter,i_solution)=4;  % No steady-state lumen
                            Novalue_count(i_solution)=Novalue_count(i_solution)+1;
                            cycle_count=cycle_count+logical(i_solution==3);
                        end
                    else
                        Lumen_phase(i_phase1_alter,i_phase2_alter,i_solution)=4;  % No steady-state lumen
                        cycle_count=cycle_count+logical(i_solution==3);
                    end
                    
                    % 若上一列寻找结果中的没有非零稳定解的连续网格数等于上一列网格数，则将该列的非零稳定解的连续网格数设为该列网格数
                    if Novalue_count_before(i_solution)==abs(y_search_start-y_search_end)+1
                        Novalue_count(i_solution)=Novalue_count_before(i_solution);
                    end
                    
                    % 输出当前循环计算结果
                    if i_solution==3
                        sprintf('calculation_times=%d, alternative_situation=%d, grid_localization=%d (total: %d), Lumen_phase=%d & %d & %d, thetas=%d & %d & %d (%d,%d)', ...
                            calculation_times,alternative_situation,cycle_count,gridsize^2,Lumen_phase(i_phase1_alter,i_phase2_alter,:),thetas_NS(i_phase1_alter,i_phase2_alter,:),i_phase1_alter,i_phase2_alter)
                    end
                end
            end
            Novalue_count_before=Novalue_count;
        end
    end


    %%% 对二维参数平面的每个点对应的三个解进行排序
    for search_region=1:4  % 搜索范围以绝对稳定的原始参数数值点为中心，将相平面分成四个部分计算
        x_search_start=x_base;y_search_start=y_base;
        if search_region==1
            x_search_end=gridsize;y_search_end=gridsize;
        elseif search_region==2
            x_search_end=gridsize;y_search_end=1;
        elseif search_region==3
            x_search_end=1;y_search_end=gridsize;
        elseif search_region==4
            x_search_end=1;y_search_end=1;
        end

        if x_search_end>=x_search_start
            x_search_interval=1;
        else
            x_search_interval=-1;
        end
        if y_search_end>=y_search_start
            y_search_interval=1;
        else
            y_search_interval=-1;
        end

        for i_phase1_alter=x_search_start:x_search_interval:x_search_end
            for i_phase2_alter=y_search_start:y_search_interval:y_search_end
                % 首先删除重复解
                if thetas_NS(i_phase1_alter,i_phase2_alter,1)~=0 && thetas_NS(i_phase1_alter,i_phase2_alter,1)==thetas_NS(i_phase1_alter,i_phase2_alter,2)
                    Rs_NS(i_phase1_alter,i_phase2_alter,2)=0;
                    thetas_NS(i_phase1_alter,i_phase2_alter,2)=0;
                    deltaCs_NS(i_phase1_alter,i_phase2_alter,2)=0;
                    Lumen_phase(i_phase1_alter,i_phase2_alter,2)=4;
                end
                if thetas_NS(i_phase1_alter,i_phase2_alter,1)~=0 && thetas_NS(i_phase1_alter,i_phase2_alter,1)==thetas_NS(i_phase1_alter,i_phase2_alter,3)
                    Rs_NS(i_phase1_alter,i_phase2_alter,3)=0;
                    thetas_NS(i_phase1_alter,i_phase2_alter,3)=0;
                    deltaCs_NS(i_phase1_alter,i_phase2_alter,3)=0;
                    Lumen_phase(i_phase1_alter,i_phase2_alter,3)=4;
                end
                if thetas_NS(i_phase1_alter,i_phase2_alter,2)~=0 && thetas_NS(i_phase1_alter,i_phase2_alter,2)==thetas_NS(i_phase1_alter,i_phase2_alter,3)
                    Rs_NS(i_phase1_alter,i_phase2_alter,3)=0;
                    thetas_NS(i_phase1_alter,i_phase2_alter,3)=0;
                    deltaCs_NS(i_phase1_alter,i_phase2_alter,3)=0;
                    Lumen_phase(i_phase1_alter,i_phase2_alter,3)=4;               
                end

                % 排序
                if sum(thetas_NS(i_phase1_alter,i_phase2_alter,:)==0)==0  % 三个解均不为零，按thetas从小到大排序
                    [thetas_NS(i_phase1_alter,i_phase2_alter,:),I_sol]=sort(thetas_NS(i_phase1_alter,i_phase2_alter,:));
                    Rs_NS(i_phase1_alter,i_phase2_alter,:)=Rs_NS(i_phase1_alter,i_phase2_alter,I_sol);
                    deltaCs_NS(i_phase1_alter,i_phase2_alter,:)=deltaCs_NS(i_phase1_alter,i_phase2_alter,I_sol);
                    Lumen_phase(i_phase1_alter,i_phase2_alter,:)=Lumen_phase(i_phase1_alter,i_phase2_alter,I_sol);
                elseif sum(thetas_NS(i_phase1_alter,i_phase2_alter,:)==0)==1  % 存在一个零解，将非零解按thetas从小到大排序
                    [thetas_NS(i_phase1_alter,i_phase2_alter,:),I_sol]=sort(thetas_NS(i_phase1_alter,i_phase2_alter,:));
                    I_sol_temp=zeros(size(I_sol));
                    I_sol_temp(1:2)=I_sol(2:3);
                    I_sol_temp(3)=I_sol(1);
                    I_sol=I_sol_temp;
                    Rs_NS(i_phase1_alter,i_phase2_alter,:)=Rs_NS(i_phase1_alter,i_phase2_alter,I_sol);
                    deltaCs_NS(i_phase1_alter,i_phase2_alter,:)=deltaCs_NS(i_phase1_alter,i_phase2_alter,I_sol);
                    Lumen_phase(i_phase1_alter,i_phase2_alter,:)=Lumen_phase(i_phase1_alter,i_phase2_alter,I_sol);
                else  % 存在两个零解，寻找与存在三个解的参数坐标点排序后的解中最接近的解位置作为该坐标点的解位置
                    if i_phase1_alter==x_search_start && i_phase2_alter~=y_search_start
                        thetas_NS_before=thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval,:);

                        Rs_NS_temp=zeros(1,3);  % 储存正确解顺序的临时变量
                        thetas_NS_temp=zeros(1,3);
                        deltaCs_NS_temp=zeros(1,3);
                        Lumen_phase_temp=ones(1,3).*4;

                        for i_solution=1:3
                            if thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)~=0
                                if min(abs(thetas_NS_before-thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)))<=0.2  % 若最接近的角度的差异小于0.2，说明该解计算准确，直接将对应位置赋给该解
                                    [~,I_sol]=min(abs(thetas_NS_before-thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)));
                                    Rs_NS_temp(I_sol)=Rs_NS(i_phase1_alter,i_phase2_alter,i_solution);
                                    thetas_NS_temp(I_sol)=thetas_NS(i_phase1_alter,i_phase2_alter,i_solution);
                                    deltaCs_NS_temp(I_sol)=deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution);
                                    Lumen_phase_temp(I_sol)=Lumen_phase(i_phase1_alter,i_phase2_alter,i_solution);
                                else  % 若最接近的角度的差异大于等于0.2，说明该解计算不准确，往上一个网格搜索，直到找出正确解为止
                                    i_temp=1;
                                    while min(abs(thetas_NS_before-thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)))>0.2
                                        i_temp=i_temp+1;
                                        thetas_NS_before=thetas_NS(i_phase1_alter,i_phase2_alter-y_search_interval*i_temp,:);
                                    end

                                    [~,I_sol]=min(abs(thetas_NS_before-thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)));
                                    Rs_NS_temp(I_sol)=Rs_NS(i_phase1_alter,i_phase2_alter,i_solution);
                                    thetas_NS_temp(I_sol)=thetas_NS(i_phase1_alter,i_phase2_alter,i_solution);
                                    deltaCs_NS_temp(I_sol)=deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution);
                                    Lumen_phase_temp(I_sol)=Lumen_phase(i_phase1_alter,i_phase2_alter,i_solution);
                                end
                            end
                        end

                        Rs_NS(i_phase1_alter,i_phase2_alter,:)=Rs_NS_temp;
                        thetas_NS(i_phase1_alter,i_phase2_alter,:)=thetas_NS_temp;
                        deltaCs_NS(i_phase1_alter,i_phase2_alter,:)=deltaCs_NS_temp;
                        Lumen_phase(i_phase1_alter,i_phase2_alter,:)=Lumen_phase_temp;
                    elseif i_phase1_alter~=x_search_start
                        thetas_NS_before=thetas_NS(i_phase1_alter-x_search_interval,i_phase2_alter,:);

                        Rs_NS_temp=zeros(1,3);  % 储存正确解顺序的临时变量
                        thetas_NS_temp=zeros(1,3);
                        deltaCs_NS_temp=zeros(1,3);
                        Lumen_phase_temp=ones(1,3).*4;
                        for i_solution=1:3
                            if thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)~=0
                                if min(abs(thetas_NS_before-thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)))<=0.2  % 若最接近的角度的差异小于0.2，说明该解计算准确，直接将对应位置赋给该解
                                    [~,I_sol]=min(abs(thetas_NS_before-thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)));
                                    Rs_NS_temp(I_sol)=Rs_NS(i_phase1_alter,i_phase2_alter,i_solution);
                                    thetas_NS_temp(I_sol)=thetas_NS(i_phase1_alter,i_phase2_alter,i_solution);
                                    deltaCs_NS_temp(I_sol)=deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution);
                                    Lumen_phase_temp(I_sol)=Lumen_phase(i_phase1_alter,i_phase2_alter,i_solution);
                                else  % 若最接近的角度的差异大于等于0.2，说明该解计算不准确，往上一个网格搜索，直到找出正确解为止
                                    i_temp=1;
                                    while min(abs(thetas_NS_before-thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)))>0.2
                                        i_temp=i_temp+1;
                                        thetas_NS_before=thetas_NS(i_phase1_alter-x_search_interval*i_temp,i_phase2_alter,:);
                                    end

                                    [~,I_sol]=min(abs(thetas_NS_before-thetas_NS(i_phase1_alter,i_phase2_alter,i_solution)));
                                    Rs_NS_temp(I_sol)=Rs_NS(i_phase1_alter,i_phase2_alter,i_solution);
                                    thetas_NS_temp(I_sol)=thetas_NS(i_phase1_alter,i_phase2_alter,i_solution);
                                    deltaCs_NS_temp(I_sol)=deltaCs_NS(i_phase1_alter,i_phase2_alter,i_solution);
                                    Lumen_phase_temp(I_sol)=Lumen_phase(i_phase1_alter,i_phase2_alter,i_solution);
                                end
                            end
                        end

                        Rs_NS(i_phase1_alter,i_phase2_alter,:)=Rs_NS_temp;
                        thetas_NS(i_phase1_alter,i_phase2_alter,:)=thetas_NS_temp;
                        deltaCs_NS(i_phase1_alter,i_phase2_alter,:)=deltaCs_NS_temp;
                        Lumen_phase(i_phase1_alter,i_phase2_alter,:)=Lumen_phase_temp;
                    end
                end
            end
        end
    end

    %%% 计算二维参数平面中的最稳定的稳定解所对应的管腔所处的相
    for i_phase1_alter=1:gridsize
        for i_phase2_alter=1:gridsize
            % 根据稳定解稳定性的优先级判断管腔所处的相
            if ismember(1,Lumen_phase(i_phase1_alter,i_phase2_alter,:))
                Lumen_phase_MS(i_phase1_alter,i_phase2_alter)=1;  % Lumen is stable
            elseif ismember(3,Lumen_phase(i_phase1_alter,i_phase2_alter,:))
                Lumen_phase_MS(i_phase1_alter,i_phase2_alter)=3;  % Lumen is oscillated
            elseif ismember(2,Lumen_phase(i_phase1_alter,i_phase2_alter,:))
                Lumen_phase_MS(i_phase1_alter,i_phase2_alter)=2;  % Lumen is unstable
            else
                Lumen_phase_MS(i_phase1_alter,i_phase2_alter)=4;  % No steady-state lumen
            end

            % 若在该网格存在两个稳定的稳定解(稳定或振荡)，则优先选择半径较大的解的稳定性作为该网格的管腔所处的相
            if ismember(1,Lumen_phase(i_phase1_alter,i_phase2_alter,:)) && ismember(3,Lumen_phase(i_phase1_alter,i_phase2_alter,:))
                i_position_stable=find(Lumen_phase(i_phase1_alter,i_phase2_alter,:)==1);
                i_position_oscillation=find(Lumen_phase(i_phase1_alter,i_phase2_alter,:)==3);
                if Rs_NS(i_phase1_alter,i_phase2_alter,i_position_stable)>Rs_NS(i_phase1_alter,i_phase2_alter,i_position_oscillation)
                    Lumen_phase_MS(i_phase1_alter,i_phase2_alter)=1;  % Lumen is stable
                else
                    Lumen_phase_MS(i_phase1_alter,i_phase2_alter)=3;  % Lumen is oscillated
                end
            end
        end
    end


    %%% 画图
    for aa=1
        if alternative_situation==1  % j和jw
            x_plot=ji_alter;  y_plot=jw_alter;
        elseif alternative_situation==2  % j和GammajHat
            x_plot=ji_alter;  y_plot=gammajHat_alter;
        elseif alternative_situation==3  % j和SigmarHat
            x_plot=ji_alter;  y_plot=SigmarHat_alter;
        elseif alternative_situation==4  % j和jpHat
            x_plot=ji_alter;  y_plot=jpHat_alter;
        elseif alternative_situation==5  % j和jpwHat
            x_plot=ji_alter;  y_plot=jpwHat_alter;
        elseif alternative_situation==6  % jpHat和jpwHat
            x_plot=jpHat_alter;  y_plot=jpwHat_alter;
        elseif alternative_situation==7  % gammajHat和SigmarHat
            x_plot=gammajHat_alter;  y_plot=SigmarHat_alter;
        elseif alternative_situation==8  % SigmarHat和jpwHat
            x_plot=SigmarHat_alter;  y_plot=jpwHat_alter;
        end

        % 修正若干计算错误的网格点
        if alternative_situation==1
            Lumen_phase_MS(351,680)=2;Lumen_phase_MS(350,681:685)=2;Lumen_phase_MS(357,706)=1;Lumen_phase_MS(361:415,732:968)=3;
        elseif alternative_situation==3
            for i_phase1_alter=437:461
                for i_phase2_alter=618:799
                    if Lumen_phase_MS(i_phase1_alter,i_phase2_alter)==4
                        Lumen_phase_MS(i_phase1_alter,i_phase2_alter)=1;
                    end
                end
            end
            Lumen_phase_MS(490,361:473)=1;Lumen_phase_MS(219:222,2:18)=2;Lumen_phase_MS(526:599,1)=1;
        elseif alternative_situation==4
            Lumen_phase_MS(50,22)=1;Lumen_phase_MS(49,21)=2;Lumen_phase_MS(49,1:14)=1;
        elseif alternative_situation==5
        end

        % 画阴影区域
        for i_phase1_alter=1:gridsize
            i_phase2_start=1;
            for y_cycle=1:10
                if i_phase2_start<=gridsize-1
                    for i_phase2_alter=i_phase2_start:gridsize-1
                        if Lumen_phase_MS(i_phase1_alter,i_phase2_alter)~=Lumen_phase_MS(i_phase1_alter,i_phase2_alter+1)
                            break
                        end
                    end
                    if i_phase2_alter==gridsize-1 && Lumen_phase_MS(i_phase1_alter,i_phase2_alter)==Lumen_phase_MS(i_phase1_alter,i_phase2_alter+1)
                        i_phase2_alter=gridsize;
                    end
                else
                    i_phase2_start=gridsize;
                    i_phase2_alter=gridsize;
                end

                x_interval=x_plot(2)-x_plot(1);
                y_interval=y_plot(2)-y_plot(1);

                if Lumen_phase_MS(i_phase1_alter,i_phase2_start)==1
                    fill([x_plot(i_phase1_alter)-x_interval x_plot(i_phase1_alter) x_plot(i_phase1_alter) x_plot(i_phase1_alter)-x_interval], ...
                        [y_plot(i_phase2_start)-y_interval y_plot(i_phase2_start)-y_interval y_plot(i_phase2_alter) y_plot(i_phase2_alter)],[0.95 0.4 0.4],'EdgeColor','none')  % red
                elseif Lumen_phase_MS(i_phase1_alter,i_phase2_start)==2
                    fill([x_plot(i_phase1_alter)-x_interval x_plot(i_phase1_alter) x_plot(i_phase1_alter) x_plot(i_phase1_alter)-x_interval], ...
                        [y_plot(i_phase2_start)-y_interval y_plot(i_phase2_start)-y_interval y_plot(i_phase2_alter) y_plot(i_phase2_alter)],[0.6 0.8 1],'EdgeColor','none')  % blue
                elseif Lumen_phase_MS(i_phase1_alter,i_phase2_start)==3
                    fill([x_plot(i_phase1_alter)-x_interval x_plot(i_phase1_alter) x_plot(i_phase1_alter) x_plot(i_phase1_alter)-x_interval], ...
                        [y_plot(i_phase2_start)-y_interval y_plot(i_phase2_start)-y_interval y_plot(i_phase2_alter) y_plot(i_phase2_alter)],[0.2 1 0.4],'EdgeColor','none')  % green
                elseif Lumen_phase_MS(i_phase1_alter,i_phase2_start)==4
                    fill([x_plot(i_phase1_alter)-x_interval x_plot(i_phase1_alter) x_plot(i_phase1_alter) x_plot(i_phase1_alter)-x_interval], ...
                        [y_plot(i_phase2_start)-y_interval y_plot(i_phase2_start)-y_interval y_plot(i_phase2_alter) y_plot(i_phase2_alter)],[0.8 0.8 0.8],'EdgeColor','none')  % gray
                end
                hold on

                i_phase2_start=i_phase2_alter+1;
                if i_phase2_start==gridsize+1
                    break
                end
            end
        end

        % 画边缘交界线
        boundary_x_summarize=ones(gridsize^2,2)*1e5;  % 相边界x坐标汇总矩阵
        boundary_y_summarize=ones(gridsize^2,2)*1e5;  % 相边界y坐标汇总矩阵
        boundary_colorsummarize=zeros(gridsize^2,2);  % 相边界绘制颜色汇总矩阵
        boundary_count=1;  % 记录相边界坐标汇总矩阵赋值的位置
        for i_phase1_alter=1:gridsize-1
            for i_phase2_alter=1:gridsize-1
                if Lumen_phase_MS(i_phase1_alter,i_phase2_alter)~=Lumen_phase_MS(i_phase1_alter,i_phase2_alter+1)  % 画与x轴方向平行的边界线
                    I_x=[find(abs(boundary_x_summarize-x_plot(i_phase1_alter)+x_interval)<1e-10)' find(abs(boundary_x_summarize-x_plot(i_phase1_alter))<1e-10)'];
                    I_y=find(abs(boundary_y_summarize(I_x)-y_plot(i_phase2_alter))<1e-10);

                    if ~isempty(I_x(I_y)) && length(I_x(I_y))==1 && (boundary_colorsummarize(I_x(I_y))==Lumen_phase_MS(i_phase1_alter,i_phase2_alter) || boundary_colorsummarize(I_x(I_y))==Lumen_phase_MS(i_phase1_alter,i_phase2_alter+1))
                        if boundary_colorsummarize(I_x(I_y))==1
                            plot([x_plot(i_phase1_alter)-x_interval x_plot(i_phase1_alter)],[y_plot(i_phase2_alter) y_plot(i_phase2_alter)],'LineWidth',3,'color',[0.5 0 0])
                        elseif boundary_colorsummarize(I_x(I_y))==2
                            plot([x_plot(i_phase1_alter)-x_interval x_plot(i_phase1_alter)],[y_plot(i_phase2_alter) y_plot(i_phase2_alter)],'LineWidth',3,'color','b')
                        elseif boundary_colorsummarize(I_x(I_y))==3
                            plot([x_plot(i_phase1_alter)-x_interval x_plot(i_phase1_alter)],[y_plot(i_phase2_alter) y_plot(i_phase2_alter)],'LineWidth',3,'color',[0 0.5 0])
                        elseif boundary_colorsummarize(I_x(I_y))==4
                            plot([x_plot(i_phase1_alter)-x_interval x_plot(i_phase1_alter)],[y_plot(i_phase2_alter) y_plot(i_phase2_alter)],'LineWidth',3,'color',[0.2 0.2 0.2])
                        end
                        boundary_colorsummarize(boundary_count,:)=boundary_colorsummarize(I_x(I_y));
                    else
                        if Lumen_phase_MS(i_phase1_alter,i_phase2_alter)==1
                            plot([x_plot(i_phase1_alter)-x_interval x_plot(i_phase1_alter)],[y_plot(i_phase2_alter) y_plot(i_phase2_alter)],'LineWidth',3,'color',[0.5 0 0])
                        elseif Lumen_phase_MS(i_phase1_alter,i_phase2_alter)==2
                            plot([x_plot(i_phase1_alter)-x_interval x_plot(i_phase1_alter)],[y_plot(i_phase2_alter) y_plot(i_phase2_alter)],'LineWidth',3,'color','b')
                        elseif Lumen_phase_MS(i_phase1_alter,i_phase2_alter)==3
                            plot([x_plot(i_phase1_alter)-x_interval x_plot(i_phase1_alter)],[y_plot(i_phase2_alter) y_plot(i_phase2_alter)],'LineWidth',3,'color',[0 0.5 0])
                        elseif Lumen_phase_MS(i_phase1_alter,i_phase2_alter)==4
                            plot([x_plot(i_phase1_alter)-x_interval x_plot(i_phase1_alter)],[y_plot(i_phase2_alter) y_plot(i_phase2_alter)],'LineWidth',3,'color',[0.2 0.2 0.2])
                        end
                        boundary_colorsummarize(boundary_count,:)=Lumen_phase_MS(i_phase1_alter,i_phase2_alter);
                    end

                    boundary_x_summarize(boundary_count,:)=x_plot(i_phase1_alter:i_phase1_alter+1)-x_interval;
                    boundary_y_summarize(boundary_count,:)=y_plot(i_phase2_alter);
                    boundary_count=boundary_count+1;
                end
                hold on

                if Lumen_phase_MS(i_phase1_alter,i_phase2_alter)~=Lumen_phase_MS(i_phase1_alter+1,i_phase2_alter)  % 画与y轴方向平行的边界线
                    I_x=find(abs(boundary_x_summarize-x_plot(i_phase1_alter)-x_interval)<1e-10);
                    I_y=[find(abs(boundary_y_summarize(I_x)-y_plot(i_phase2_alter)+y_interval)<1e-10) find(abs(boundary_y_summarize(I_x)-y_plot(i_phase2_alter))<1e-10)];

                    if ~isempty(I_x(I_y)) && length(I_x(I_y))==1 && (boundary_colorsummarize(I_x(I_y))==Lumen_phase_MS(i_phase1_alter,i_phase2_alter) || boundary_colorsummarize(I_x(I_y))==Lumen_phase_MS(i_phase1_alter+1,i_phase2_alter))
                        if boundary_colorsummarize(I_x(I_y))==1
                            plot([x_plot(i_phase1_alter) x_plot(i_phase1_alter)],[y_plot(i_phase2_alter)-y_interval y_plot(i_phase2_alter)],'LineWidth',3,'color',[0.5 0 0])
                        elseif boundary_colorsummarize(I_x(I_y))==2
                            plot([x_plot(i_phase1_alter) x_plot(i_phase1_alter)],[y_plot(i_phase2_alter)-y_interval y_plot(i_phase2_alter)],'LineWidth',3,'color','b')
                        elseif boundary_colorsummarize(I_x(I_y))==3
                            plot([x_plot(i_phase1_alter) x_plot(i_phase1_alter)],[y_plot(i_phase2_alter)-y_interval y_plot(i_phase2_alter)],'LineWidth',3,'color',[0 0.5 0])
                        elseif boundary_colorsummarize(I_x(I_y))==4
                            plot([x_plot(i_phase1_alter) x_plot(i_phase1_alter)],[y_plot(i_phase2_alter)-y_interval y_plot(i_phase2_alter)],'LineWidth',3,'color',[0.2 0.2 0.2])
                        end
                        boundary_colorsummarize(boundary_count,:)=boundary_colorsummarize(I_x(I_y));
                    else
                        if Lumen_phase_MS(i_phase1_alter,i_phase2_alter)==1
                            plot([x_plot(i_phase1_alter) x_plot(i_phase1_alter)],[y_plot(i_phase2_alter)-y_interval y_plot(i_phase2_alter)],'LineWidth',3,'color',[0.5 0 0])
                        elseif Lumen_phase_MS(i_phase1_alter,i_phase2_alter)==2
                            plot([x_plot(i_phase1_alter) x_plot(i_phase1_alter)],[y_plot(i_phase2_alter)-y_interval y_plot(i_phase2_alter)],'LineWidth',3,'color','b')
                        elseif Lumen_phase_MS(i_phase1_alter,i_phase2_alter)==3
                            plot([x_plot(i_phase1_alter) x_plot(i_phase1_alter)],[y_plot(i_phase2_alter)-y_interval y_plot(i_phase2_alter)],'LineWidth',3,'color',[0 0.5 0])
                        elseif Lumen_phase_MS(i_phase1_alter,i_phase2_alter)==4
                            plot([x_plot(i_phase1_alter) x_plot(i_phase1_alter)],[y_plot(i_phase2_alter)-y_interval y_plot(i_phase2_alter)],'LineWidth',3,'color',[0.2 0.2 0.2])
                        end
                        boundary_colorsummarize(boundary_count,:)=Lumen_phase_MS(i_phase1_alter,i_phase2_alter);
                    end

                    boundary_x_summarize(boundary_count,:)=x_plot(i_phase1_alter);
                    boundary_y_summarize(boundary_count,:)=y_plot(i_phase2_alter:i_phase2_alter+1)-y_interval;
                    boundary_count=boundary_count+1;
                end
                hold on
            end
        end

        % 画零坐标虚线
        if x_interval>0
            plot([0 max(x_plot)],[0 0],'--','LineWidth',1,'Color','k')  % 画y=0横线
        else
            plot([0 min(x_plot)],[0 0],'--','LineWidth',1,'Color','k')  % 画y=0横线
        end
        if y_interval>0
            plot([0 0],[0 max(y_plot)],'--','LineWidth',1,'Color','k')  % 画x=0横线
        else
            plot([0 0],[0 min(y_plot)],'--','LineWidth',1,'Color','k')  % 画x=0横线
        end

        % x_label和y_label内容赋值
        if alternative_situation==1  % j和jw
            x_label='j';y_label='jw';
            x_label_latex='$J^P\kappa_BT\Lambda^w/\Lambda^2$';y_label_latex='$J^W/\Lambda$';
            x_label_min=min(ji_alter);x_label_max=max(ji_alter);
            y_label_min=min(jw_alter);y_label_max=max(jw_alter);
        elseif alternative_situation==2  % j和GammajHat
            x_label='j';y_label='GammajHat';
            x_label_latex='$J^P\kappa_BT\Lambda^w/\Lambda^2$';y_label_latex='$\gamma_j/\Gamma$';
            x_label_min=min(ji_alter);x_label_max=max(ji_alter);
            y_label_min=min(gammajHat_alter);y_label_max=max(gammajHat_alter);
        elseif alternative_situation==3  % j和SigmarHat
            x_label='j';y_label='SigmarHat';
            x_label_latex='$J^P\kappa_BT\Lambda^w/\Lambda^2$';y_label_latex='$\sigma_r/L\Gamma$';
            x_label_min=min(ji_alter);x_label_max=max(ji_alter);
            y_label_min=min(SigmarHat_alter);y_label_max=max(SigmarHat_alter);
        elseif alternative_situation==4  % j和jpHat
            x_label='j';y_label='jpHat';
            x_label_latex='$J^P\kappa_BT\Lambda^w/\Lambda^2$';y_label_latex='$j_p\kappa_BT\Lambda^w/L^2\Lambda^2$';
            x_label_min=min(ji_alter);x_label_max=max(ji_alter);
            y_label_min=min(jpHat_alter);y_label_max=max(jpHat_alter);
        elseif alternative_situation==5  % j和jpwHat
            x_label='j';y_label='jpwHat';
            x_label_latex='$J^P\kappa_BT\Lambda^w/\Lambda^2$';y_label_latex='$j_p^w/L^2\Lambda$';
            x_label_min=min(ji_alter);x_label_max=max(ji_alter);
            y_label_min=min(jpwHat_alter);y_label_max=max(jpwHat_alter);
        elseif alternative_situation==6  % jpHat和jpwHat
            x_label='jpHat';y_label='jpwHat';
            x_label_latex='$j_p\kappa_BT\Lambda^w/L^2\Lambda^2$';y_label_latex='$j_p^w/L^2\Lambda$';
            x_label_min=min(jpHat_alter);x_label_max=max(jpHat_alter);
            y_label_min=min(jpwHat_alter);y_label_max=max(jpwHat_alter);
        elseif alternative_situation==7  % gammajHat和SigmarHat
            x_label='gammajHat';y_label='jpwHat';
            x_label_latex='$\gamma_j/\Gamma$';y_label_latex='$j_p^w/L^2\Lambda$';
            x_label_min=min(gammajHat_alter);x_label_max=max(gammajHat_alter);
            y_label_min=min(SigmarHat_alter);y_label_max=max(SigmarHat_alter);
        elseif alternative_situation==8  % SigmarHat和jpwHat
            x_label='SigmarHat';y_label='jpwHat';
            x_label_latex='$\sigma_r/L\Gamma$';y_label_latex='$j_p^w/L^2\Lambda$';
            x_label_min=min(SigmarHat_alter);x_label_max=max(SigmarHat_alter);
            y_label_min=min(jpwHat_alter);y_label_max=max(jpwHat_alter);
        end

        set(gca,'FontName','Arial','FontSize',20);
        filename_plot=strcat("Phase diagram of consider leak (",x_label," and ",y_label,")");
        % title(filename_plot)
        hold off

        axis ([x_label_min x_label_max y_label_min y_label_max]);
        xlabel(x_label_latex,'interpreter','latex')
        ylabel(y_label_latex,'interpreter','latex')
        label_interval_temp=[1e-8 2e-8 5e-8 1e-7 2e-7 5e-7 1e-6 2e-6 5e-6 1e-5 2e-5 5e-5 1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2 2e-2 5e-2 1e-1 2e-1 5e-1];  % 待选取的坐标轴标注数值间隔
        label_interval=[label_interval_temp label_interval_temp.*1e8];
        [~,I]=min(abs(label_interval(:)-(x_label_max-x_label_min)/5));x_label_interval=label_interval(I);  % x轴坐标标注间隔的选取
        [~,I]=min(abs(label_interval(:)-(y_label_max-y_label_min)/5));y_label_interval=label_interval(I);  % y轴坐标标注间隔的选取
        if max(x_plot)>0
            set(gca,'XTick',0:x_label_interval:x_label_max);  % 修改x坐标刻度
        else
            set(gca,'XTick',x_label_min:x_label_interval:0);  % 修改x坐标刻度
        end
        if max(y_plot)>0
            set(gca,'YTick',0:y_label_interval:y_label_max);  % 修改y坐标刻度
        else
            set(gca,'YTick',y_label_min:y_label_interval:0);  % 修改y坐标刻度
        end
        set(gca,'FontName','Arial','FontSize',18)
        set(gcf,'unit','centimeters','position',[0,0,20,18])  % 调整照片大小
    end
end


