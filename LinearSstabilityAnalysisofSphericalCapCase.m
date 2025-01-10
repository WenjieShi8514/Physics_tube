%%% Linear stability analysis of spherical cap case
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



%% Linear Stability Analysis
for alternative_situation=1:7
    step=1000;
    thetas_NS=zeros(step,3);
    Rs_NS=zeros(step,3);
    deltaCs_NS=zeros(step,3);
    eigenvalue_allsolution=zeros(step,3,2);  % Stores eigenvalues for each parameter value. The second dimension represents two stable solutions, and the third dimension represents two eigenvalues for each stable solution.

    %%% Types of x-axis in Linear Stability Analysis (1 for j, 2 for jw, 3 for GammajHat, 4 for SigmarHat, 5 for jpHat, 6 for jpwHat)
    if alternative_situation==1
        ji_alter=linspace(0.5*ji,1.5*ji,step);
        [~,I]=min(abs(ji_alter(:)-ji));x_base=I;  % Selection of the closest position to the original value in the x-axis range
    elseif alternative_situation==2
        jw_alter=linspace(0,10*jw,step);
        [~,I]=min(abs(jw_alter(:)-jw));x_base=I;  % Selection of the closest position to the original value in the x-axis range
    elseif alternative_situation==3
        gammajHat_alter=linspace(0,2*gammajHat,step);
        [~,I]=min(abs(gammajHat_alter(:)-gammajHat));x_base=I;  % Selection of the closest position to the original value in the x-axis range
    elseif alternative_situation==4
        SigmarHat_alter=linspace(0,2*SigmarHat,step);
        [~,I]=min(abs(SigmarHat_alter(:)-SigmarHat));x_base=I;  % Selection of the closest position to the original value in the x-axis range
    elseif alternative_situation==5
        jpHat_alter=linspace(0,2*jpHat,step);
        [~,I]=min(abs(jpHat_alter(:)-jpHat));x_base=I;  % Selection of the closest position to the original value in the x-axis range
    elseif alternative_situation==6
        jpwHat_alter=linspace(0,2*jpwHat,step);
        [~,I]=min(abs(jpwHat_alter(:)-jpwHat));x_base=I;  % Selection of the closest position to the original value in the x-axis range
    elseif alternative_situation==7
        gamma_alter=linspace(0,200*gamma,step);
        [~,I]=min(abs(gamma_alter(:)-gamma));x_base=I;  % Selection of the closest position to the original value in the x-axis range
    end

    % Stable solution calculation and stability analysis for subsequent baseline parameter values
    Novalue_count1=0;  % Stores the number of grids without stable solutions
    Novalue_count2=0;  % Stores the number of grids without stable solutions
    Novalue_count3=0;  % Stores the number of grids without stable solutions
    for i=x_base:step
        %%% alternative dimensionless parameters assignment
        ji_alter_temp=ji;jw_alter_temp=jw;gammajHat_alter_temp=gammajHat;SigmarHat_alter_temp=SigmarHat;jpHat_alter_temp=jpHat;jpwHat_alter_temp=jpwHat;gamma_alter_temp=gamma;
        if alternative_situation==1
            ji_alter_temp=ji_alter(i);
        elseif alternative_situation==2
            jw_alter_temp=jw_alter(i);
        elseif alternative_situation==3
            gammajHat_alter_temp=gammajHat_alter(i);
        elseif alternative_situation==4
            SigmarHat_alter_temp=SigmarHat_alter(i);
        elseif alternative_situation==5
            jpHat_alter_temp=jpHat_alter(i);
        elseif alternative_situation==6
            jpwHat_alter_temp=jpwHat_alter(i);
        elseif alternative_situation==7
            gamma_alter_temp=gamma_alter(i);
        end

        %%% Search for numerical solutions of all possible stable solutions, stored in thetas_NS, Rs_NS, and deltaCs_NS
        syms Rs thetas deltaCs
        eqns=[2*Rs*(1-cos(thetas))*(ji_alter_temp-deltaCs)==jpHat_alter_temp*sin(thetas)/(1-Rs*sin(thetas)),...
            2*Rs*(1-cos(thetas))*(jw_alter_temp-2*gamma_alter_temp/Rs+deltaCs)==jpwHat_alter_temp*sin(thetas)/(1-Rs*sin(thetas)),...
            -gammajHat_alter_temp+cos(thetas)-SigmarHat_alter_temp/Rs/sin(thetas)==0];
        vars=[Rs thetas deltaCs];

        if i==x_base  % During the first initial search for numerical solutions, the search range is categorized, with fewer iterations for easily converging solutions and more iterations for solutions that are difficult to converge.
            for ii=1:500
                temp=vpasolve(eqns,vars,[0 10;0 pi;0 ji_alter_temp],'Random',true);
                if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                    thetas_NS(i,1)=double(temp.thetas);
                    Rs_NS(i,1)=double(temp.Rs);
                    deltaCs_NS(i,1)=double(temp.deltaCs);
                    break
                end
            end
            for ii=1:500
                temp=vpasolve(eqns,vars,[0 10;0 pi;0 ji_alter_temp],'Random',true);
                if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && double(temp.thetas)~=thetas_NS(i,1)
                    thetas_NS(i,2)=double(temp.thetas);
                    Rs_NS(i,2)=double(temp.Rs);
                    deltaCs_NS(i,2)=double(temp.deltaCs);
                    break
                end
            end
            for ii=1:500
                temp=vpasolve(eqns,vars,[0 10;0 pi;0 ji_alter_temp],'Random',true);
                if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && (double(temp.thetas)~=thetas_NS(i,1) && double(temp.thetas)~=thetas_NS(i,2))
                    thetas_NS(i,3)=double(temp.thetas);
                    Rs_NS(i,3)=double(temp.Rs);
                    deltaCs_NS(i,3)=double(temp.deltaCs);
                    break
                end
            end
        else  % In subsequent iterations, the search for numerical solutions is based on the results of the previous iteration, searching within a small surrounding range.
            % 第一个解
            for aa=1
                % If the parameter values deviate too far from reality, leading to no stable solution, terminate the iteration early.
                if thetas_NS(i-1,1)==0
                    Novalue_count1=Novalue_count1+1;
                else
                    Novalue_count1=0;
                end
                if Novalue_count1>=3 && sum(thetas_NS(:,1))>20
                    break
                end

                for ii=1:10
                    temp=vpasolve(eqns,vars,[Rs_NS(i-1,1)*0.9 Rs_NS(i-1,1)*1.1;max(0,thetas_NS(i-1,1)-pi/100) thetas_NS(i-1,1)+pi/100;deltaCs_NS(i-1,1)*0.9 deltaCs_NS(i-1,1)*1.1],'Random',true);
                    if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                        thetas_NS(i,1)=double(temp.thetas);
                        Rs_NS(i,1)=double(temp.Rs);
                        deltaCs_NS(i,1)=double(temp.deltaCs);
                        break
                    end
                end
                if thetas_NS(i,1)==0  % If the previous deltaCs is close to 0, expand the search range of deltaCs to negative values; otherwise, only slightly expand the search range.
                    if abs(deltaCs_NS(i-1,1))<5e-2
                        for ii=1:100
                            temp=vpasolve(eqns,vars,[Rs_NS(i-1,1)*0.99 Rs_NS(i-1,1)*1.01;max(0,thetas_NS(i-1,1)-pi/500) thetas_NS(i-1,1)+pi/500;deltaCs_NS(i-1,1)-1 abs(deltaCs_NS(i-1,1))*2],'Random',true);
                            if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                thetas_NS(i,1)=double(temp.thetas);
                                Rs_NS(i,1)=double(temp.Rs);
                                deltaCs_NS(i,1)=double(temp.deltaCs);
                                break
                            end
                        end
                    else
                        for ii=1:20
                            temp=vpasolve(eqns,vars,[Rs_NS(i-1,1)*0.8 Rs_NS(i-1,1)*1.2;max(0,thetas_NS(i-1,1)-pi/50) thetas_NS(i-1,1)+pi/50;deltaCs_NS(i-1,1)*0.5 deltaCs_NS(i-1,1)*2],'Random',true);
                            if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                thetas_NS(i,1)=double(temp.thetas);
                                Rs_NS(i,1)=double(temp.Rs);
                                deltaCs_NS(i,1)=double(temp.deltaCs);
                                break
                            end
                        end
                    end
                end
                if thetas_NS(i,1)==0  % 若容易收敛的解未找到，则扩大搜索范围
                    for ii=1:20
                        temp=vpasolve(eqns,vars,[0 1000;0 pi;0 ji_alter_temp],'Random',true);
                        if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                            thetas_NS(i,1)=double(temp.thetas);
                            Rs_NS(i,1)=double(temp.Rs);
                            deltaCs_NS(i,1)=double(temp.deltaCs);
                            break
                        end
                    end
                end
            end
            % 第二个解
            for aa=1
                % If the parameter values deviate too far from reality, leading to no stable solution, terminate the iteration early.
                if thetas_NS(i-1,2)==0
                    Novalue_count2=Novalue_count2+1;
                else
                    Novalue_count2=0;
                end
                if Novalue_count2>=3 && sum(thetas_NS(:,2))>20
                    break
                end

                for ii=1:20
                    temp=vpasolve(eqns,vars,[Rs_NS(i-1,2)*0.9 Rs_NS(i-1,2)*1.1;thetas_NS(i-1,2)-pi/100 min(pi,thetas_NS(i-1,2)+pi/100);deltaCs_NS(i-1,2)*0.9 deltaCs_NS(i-1,2)*1.1],'Random',true);
                    if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && double(temp.thetas)~=thetas_NS(i,1)
                        thetas_NS(i,2)=double(temp.thetas);
                        Rs_NS(i,2)=double(temp.Rs);
                        deltaCs_NS(i,2)=double(temp.deltaCs);
                        break
                    end
                end
                if thetas_NS(i,2)==0  % If the previous deltaCs is close to 0, expand the search range of deltaCs to negative values; otherwise, only slightly expand the search range.
                    if abs(deltaCs_NS(i-1,2))<5e-2
                        for ii=1:100
                            temp=vpasolve(eqns,vars,[Rs_NS(i-1,2)*0.9 Rs_NS(i-1,2)*1.1;max(0,thetas_NS(i-1,2)-pi/100) thetas_NS(i-1,2)+pi/100;deltaCs_NS(i-1,2)-1 abs(deltaCs_NS(i-1,2))*2],'Random',true);
                            if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                thetas_NS(i,2)=double(temp.thetas);
                                Rs_NS(i,2)=double(temp.Rs);
                                deltaCs_NS(i,2)=double(temp.deltaCs);
                                break
                            end
                        end
                    else
                        for ii=1:20
                            temp=vpasolve(eqns,vars,[Rs_NS(i-1,2)*0.8 Rs_NS(i-1,2)*1.2;max(0,thetas_NS(i-1,2)-pi/50) thetas_NS(i-1,2)+pi/50;deltaCs_NS(i-1,2)*0.5 deltaCs_NS(i-1,2)*2],'Random',true);
                            if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                thetas_NS(i,2)=double(temp.thetas);
                                Rs_NS(i,2)=double(temp.Rs);
                                deltaCs_NS(i,2)=double(temp.deltaCs);
                                break
                            end
                        end
                    end
                end
                if thetas_NS(i,2)==0  % 若不容易收敛的解仍然未找到，则扩大搜索范围
                    for ii=1:20
                        temp=vpasolve(eqns,vars,[0 1000;0 pi;0 ji_alter_temp],'Random',true);
                        if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && double(temp.thetas)~=thetas_NS(i,1)
                            thetas_NS(i,2)=double(temp.thetas);
                            Rs_NS(i,2)=double(temp.Rs);
                            deltaCs_NS(i,2)=double(temp.deltaCs);
                            break
                        end
                    end
                end
            end
            % 第三个解
            for aa=1
                % If the parameter values deviate too far from reality, leading to no stable solution, terminate the iteration early.
                if thetas_NS(i-1,3)==0
                    Novalue_count3=Novalue_count3+1;
                else
                    Novalue_count3=0;
                end
                if Novalue_count3>=3 && sum(thetas_NS(:,3))>20
                    break
                end

                for ii=1:20
                    temp=vpasolve(eqns,vars,[Rs_NS(i-1,3)*0.9 Rs_NS(i-1,3)*1.1;thetas_NS(i-1,3)-pi/100 min(pi,thetas_NS(i-1,3)+pi/100);deltaCs_NS(i-1,3)*0.9 deltaCs_NS(i-1,3)*1.1],'Random',true);
                    if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && (double(temp.thetas)~=thetas_NS(i,1) && double(temp.thetas)~=thetas_NS(i,2))
                        thetas_NS(i,3)=double(temp.thetas);
                        Rs_NS(i,3)=double(temp.Rs);
                        deltaCs_NS(i,3)=double(temp.deltaCs);
                        break
                    end
                end
                if thetas_NS(i,3)==0  % If the previous deltaCs is close to 0, expand the search range of deltaCs to negative values; otherwise, only slightly expand the search range.
                    if abs(deltaCs_NS(i-1,3))<5e-2
                        for ii=1:100
                            temp=vpasolve(eqns,vars,[Rs_NS(i-1,3)*0.9 Rs_NS(i-1,3)*1.1;max(0,thetas_NS(i-1,3)-pi/100) thetas_NS(i-1,3)+pi/100;deltaCs_NS(i-1,3)-1 abs(deltaCs_NS(i-1,3))*2],'Random',true);
                            if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                thetas_NS(i,3)=double(temp.thetas);
                                Rs_NS(i,3)=double(temp.Rs);
                                deltaCs_NS(i,3)=double(temp.deltaCs);
                                break
                            end
                        end
                    else
                        for ii=1:20
                            temp=vpasolve(eqns,vars,[Rs_NS(i-1,3)*0.8 Rs_NS(i-1,3)*1.2;max(0,thetas_NS(i-1,3)-pi/50) thetas_NS(i-1,3)+pi/50;deltaCs_NS(i-1,3)*0.5 deltaCs_NS(i-1,3)*2],'Random',true);
                            if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                thetas_NS(i,3)=double(temp.thetas);
                                Rs_NS(i,3)=double(temp.Rs);
                                deltaCs_NS(i,3)=double(temp.deltaCs);
                                break
                            end
                        end
                    end
                end
                if thetas_NS(i,3)==0  % 若不容易收敛的解仍然未找到，则扩大搜索范围
                    for ii=1:20
                        temp=vpasolve(eqns,vars,[0 1000;0 pi;0 ji_alter_temp],'Random',true);
                        if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && (double(temp.thetas)~=thetas_NS(i,1) && double(temp.thetas)~=thetas_NS(i,2))
                            thetas_NS(i,3)=double(temp.thetas);
                            Rs_NS(i,3)=double(temp.Rs);
                            deltaCs_NS(i,3)=double(temp.deltaCs);
                            break
                        end
                    end
                end
            end
        end

        sprintf('alternative_situation=%d, step=%d (total: %d), thetas=%d & %d & %d',alternative_situation,i,step,thetas_NS(i,1),thetas_NS(i,2),thetas_NS(i,3))

        %%% Stability analysis
        for ii=1:3
            if thetas_NS(i,ii)~=0
                Rs_best=Rs_NS(i,ii);
                thetas_best=thetas_NS(i,ii);
                deltaCs_best=deltaCs_NS(i,ii);

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
                eigenvalue_allsolution(i,ii,:)=eig(Mstab);
            end
        end
    end

    % 基底参数值往前的稳定解计算和稳定性分析
    Novalue_count1=0;  % Stores the number of grids without stable solutions
    Novalue_count2=0;  % Stores the number of grids without stable solutions
    Novalue_count3=0;  % Stores the number of grids without stable solutions
    for i=x_base:-1:1
        %%% alternative dimensionless parameters assignment
        ji_alter_temp=ji;jw_alter_temp=jw;gammajHat_alter_temp=gammajHat;SigmarHat_alter_temp=SigmarHat;jpHat_alter_temp=jpHat;jpwHat_alter_temp=jpwHat;
        if alternative_situation==1
            ji_alter_temp=ji_alter(i);
        elseif alternative_situation==2
            jw_alter_temp=jw_alter(i);
        elseif alternative_situation==3
            gammajHat_alter_temp=gammajHat_alter(i);
        elseif alternative_situation==4
            SigmarHat_alter_temp=SigmarHat_alter(i);
        elseif alternative_situation==5
            jpHat_alter_temp=jpHat_alter(i);
        elseif alternative_situation==6
            jpwHat_alter_temp=jpwHat_alter(i);
        end

        %%% Search for numerical solutions of all possible stable solutions, stored in thetas_NS, Rs_NS, and deltaCs_NS
        syms Rs thetas deltaCs
        eqns=[2*Rs*(1-cos(thetas))*(ji_alter_temp-deltaCs)==jpHat_alter_temp*sin(thetas)/(1-Rs*sin(thetas)),...
            2*Rs*(1-cos(thetas))*(jw_alter_temp-2*gamma_alter_temp/Rs+deltaCs)==jpwHat_alter_temp*sin(thetas)/(1-Rs*sin(thetas)),...
            -gammajHat_alter_temp+cos(thetas)-SigmarHat_alter_temp/Rs/sin(thetas)==0];
        vars=[Rs thetas deltaCs];

        if i==x_base  % During the first initial search for numerical solutions, the search range is categorized, with fewer iterations for easily converging solutions and more iterations for solutions that are difficult to converge.
            for ii=1:500
                temp=vpasolve(eqns,vars,[0 1000;0 pi;0 ji_alter_temp],'Random',true);
                if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                    thetas_NS(i,1)=double(temp.thetas);
                    Rs_NS(i,1)=double(temp.Rs);
                    deltaCs_NS(i,1)=double(temp.deltaCs);
                    break
                end
            end
            for ii=1:500
                temp=vpasolve(eqns,vars,[0 1000;0 pi;0 ji_alter_temp],'Random',true);
                if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && double(temp.thetas)~=thetas_NS(i,1)
                    thetas_NS(i,2)=double(temp.thetas);
                    Rs_NS(i,2)=double(temp.Rs);
                    deltaCs_NS(i,2)=double(temp.deltaCs);
                    break
                end
            end
            for ii=1:500
                temp=vpasolve(eqns,vars,[0 1000;0 pi;0 ji_alter_temp],'Random',true);
                if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && (double(temp.thetas)~=thetas_NS(i,1) && double(temp.thetas)~=thetas_NS(i,2))
                    thetas_NS(i,3)=double(temp.thetas);
                    Rs_NS(i,3)=double(temp.Rs);
                    deltaCs_NS(i,3)=double(temp.deltaCs);
                    break
                end
            end
        else  % In subsequent iterations, the search for numerical solutions is based on the results of the previous iteration, searching within a small surrounding range.
            % 第一个解
            for aa=1
                % If the parameter values deviate too far from reality, leading to no stable solution, terminate the iteration early.
                if thetas_NS(i+1,1)==0
                    Novalue_count1=Novalue_count1+1;
                else
                    Novalue_count1=0;
                end
                if Novalue_count1>=3 && sum(thetas_NS(:,1))>20
                    break
                end

                for ii=1:10
                    temp=vpasolve(eqns,vars,[Rs_NS(i+1,1)*0.9 Rs_NS(i+1,1)*1.1;max(0,thetas_NS(i+1,1)-pi/100) thetas_NS(i+1,1)+pi/100;deltaCs_NS(i+1,1)*0.9 deltaCs_NS(i+1,1)*1.1],'Random',true);
                    if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                        thetas_NS(i,1)=double(temp.thetas);
                        Rs_NS(i,1)=double(temp.Rs);
                        deltaCs_NS(i,1)=double(temp.deltaCs);
                        break
                    end
                end
                if thetas_NS(i,1)==0  % If the previous deltaCs is close to 0, expand the search range of deltaCs to negative values; otherwise, only slightly expand the search range.
                    if abs(deltaCs_NS(i+1,1))<5e-2
                        for ii=1:100
                            temp=vpasolve(eqns,vars,[Rs_NS(i+1,1)*0.99 Rs_NS(i+1,1)*1.01;max(0,thetas_NS(i+1,1)-pi/500) thetas_NS(i+1,1)+pi/500;deltaCs_NS(i+1,1)+1 abs(deltaCs_NS(i+1,1))*2],'Random',true);
                            if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                thetas_NS(i,1)=double(temp.thetas);
                                Rs_NS(i,1)=double(temp.Rs);
                                deltaCs_NS(i,1)=double(temp.deltaCs);
                                break
                            end
                        end
                    else
                        for ii=1:20
                            temp=vpasolve(eqns,vars,[Rs_NS(i+1,1)*0.8 Rs_NS(i+1,1)*1.2;max(0,thetas_NS(i+1,1)-pi/50) thetas_NS(i+1,1)+pi/50;deltaCs_NS(i+1,1)*0.5 deltaCs_NS(i+1,1)*2],'Random',true);
                            if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                thetas_NS(i,1)=double(temp.thetas);
                                Rs_NS(i,1)=double(temp.Rs);
                                deltaCs_NS(i,1)=double(temp.deltaCs);
                                break
                            end
                        end
                    end
                end
                if thetas_NS(i,1)==0  % 若容易收敛的解未找到，则扩大搜索范围
                    for ii=1:20
                        temp=vpasolve(eqns,vars,[0 1000;0 pi;0 ji_alter_temp],'Random',true);
                        if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                            thetas_NS(i,1)=double(temp.thetas);
                            Rs_NS(i,1)=double(temp.Rs);
                            deltaCs_NS(i,1)=double(temp.deltaCs);
                            break
                        end
                    end
                end
            end
            % 第二个解
            for aa=1
                % If the parameter values deviate too far from reality, leading to no stable solution, terminate the iteration early.
                if thetas_NS(i+1,2)==0
                    Novalue_count2=Novalue_count2+1;
                else
                    Novalue_count2=0;
                end
                if Novalue_count2>=3 && sum(thetas_NS(:,2))>20
                    break
                end

                for ii=1:20
                    temp=vpasolve(eqns,vars,[Rs_NS(i+1,2)*0.9 Rs_NS(i+1,2)*1.1;thetas_NS(i+1,2)-pi/100 min(pi,thetas_NS(i+1,2)+pi/100);deltaCs_NS(i+1,2)*0.9 deltaCs_NS(i+1,2)*1.1],'Random',true);
                    if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && double(temp.thetas)~=thetas_NS(i,1)
                        thetas_NS(i,2)=double(temp.thetas);
                        Rs_NS(i,2)=double(temp.Rs);
                        deltaCs_NS(i,2)=double(temp.deltaCs);
                        break
                    end
                end
                if thetas_NS(i,2)==0  % If the previous deltaCs is close to 0, expand the search range of deltaCs to negative values; otherwise, only slightly expand the search range.
                    if abs(deltaCs_NS(i+1,2))<5e-2
                        for ii=1:100
                            temp=vpasolve(eqns,vars,[Rs_NS(i+1,2)*0.9 Rs_NS(i+1,2)*1.1;max(0,thetas_NS(i+1,2)-pi/100) thetas_NS(i+1,2)+pi/100;deltaCs_NS(i+1,2)+1 abs(deltaCs_NS(i+1,2))*2],'Random',true);
                            if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                thetas_NS(i,2)=double(temp.thetas);
                                Rs_NS(i,2)=double(temp.Rs);
                                deltaCs_NS(i,2)=double(temp.deltaCs);
                                break
                            end
                        end
                    else
                        for ii=1:20
                            temp=vpasolve(eqns,vars,[Rs_NS(i+1,2)*0.8 Rs_NS(i+1,2)*1.2;max(0,thetas_NS(i+1,2)-pi/50) thetas_NS(i+1,2)+pi/50;deltaCs_NS(i+1,2)*0.5 deltaCs_NS(i+1,2)*2],'Random',true);
                            if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                thetas_NS(i,2)=double(temp.thetas);
                                Rs_NS(i,2)=double(temp.Rs);
                                deltaCs_NS(i,2)=double(temp.deltaCs);
                                break
                            end
                        end
                    end
                end
                if thetas_NS(i,2)==0  % 若不容易收敛的解仍然未找到，则扩大搜索范围
                    for ii=1:20
                        temp=vpasolve(eqns,vars,[0 1000;0 pi;0 ji_alter_temp],'Random',true);
                        if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && double(temp.thetas)~=thetas_NS(i,1)
                            thetas_NS(i,2)=double(temp.thetas);
                            Rs_NS(i,2)=double(temp.Rs);
                            deltaCs_NS(i,2)=double(temp.deltaCs);
                            break
                        end
                    end
                end
            end
            % 第三个解
            for aa=1
                % If the parameter values deviate too far from reality, leading to no stable solution, terminate the iteration early.
                if thetas_NS(i+1,3)==0
                    Novalue_count3=Novalue_count3+1;
                else
                    Novalue_count3=0;
                end
                if Novalue_count3>=3 && sum(thetas_NS(:,3))>20
                    break
                end

                for ii=1:20
                    temp=vpasolve(eqns,vars,[Rs_NS(i+1,3)*0.9 Rs_NS(i+1,3)*1.1;thetas_NS(i+1,3)-pi/100 min(pi,thetas_NS(i+1,3)+pi/100);deltaCs_NS(i+1,3)*0.9 deltaCs_NS(i+1,3)*1.1],'Random',true);
                    if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && (double(temp.thetas)~=thetas_NS(i,1) && double(temp.thetas)~=thetas_NS(i,2))
                        thetas_NS(i,3)=double(temp.thetas);
                        Rs_NS(i,3)=double(temp.Rs);
                        deltaCs_NS(i,3)=double(temp.deltaCs);
                        break
                    end
                end
                if thetas_NS(i,3)==0  % If the previous deltaCs is close to 0, expand the search range of deltaCs to negative values; otherwise, only slightly expand the search range.
                    if abs(deltaCs_NS(i+1,3))<5e-2
                        for ii=1:100
                            temp=vpasolve(eqns,vars,[Rs_NS(i+1,3)*0.9 Rs_NS(i+1,3)*1.1;max(0,thetas_NS(i+1,3)-pi/100) thetas_NS(i+1,3)+pi/100;deltaCs_NS(i+1,3)+1 abs(deltaCs_NS(i+1,3))*2],'Random',true);
                            if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                thetas_NS(i,3)=double(temp.thetas);
                                Rs_NS(i,3)=double(temp.Rs);
                                deltaCs_NS(i,3)=double(temp.deltaCs);
                                break
                            end
                        end
                    else
                        for ii=1:20
                            temp=vpasolve(eqns,vars,[Rs_NS(i+1,3)*0.8 Rs_NS(i+1,3)*1.2;max(0,thetas_NS(i+1,3)-pi/50) thetas_NS(i+1,3)+pi/50;deltaCs_NS(i+1,3)*0.5 deltaCs_NS(i+1,3)*2],'Random',true);
                            if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi)
                                thetas_NS(i,3)=double(temp.thetas);
                                Rs_NS(i,3)=double(temp.Rs);
                                deltaCs_NS(i,3)=double(temp.deltaCs);
                                break
                            end
                        end
                    end
                end
                if thetas_NS(i,3)==0  % 若不容易收敛的解仍然未找到，则扩大搜索范围
                    for ii=1:20
                        temp=vpasolve(eqns,vars,[0 1000;0 pi;0 ji_alter_temp],'Random',true);
                        if ~isempty(double(temp.thetas)) && (double(temp.thetas)>1e-8 && double(temp.thetas)<pi) && (double(temp.thetas)~=thetas_NS(i,1) && double(temp.thetas)~=thetas_NS(i,2))
                            thetas_NS(i,3)=double(temp.thetas);
                            Rs_NS(i,3)=double(temp.Rs);
                            deltaCs_NS(i,3)=double(temp.deltaCs);
                            break
                        end
                    end
                end
            end
        end

        sprintf('alternative_situation=%d, step=%d (total: %d), thetas=%d & %d & %d',alternative_situation,i,step,thetas_NS(i,1),thetas_NS(i,2),thetas_NS(i,3))

        %%% Stability analysis
        for ii=1:3
            if thetas_NS(i,ii)~=0
                Rs_best=Rs_NS(i,ii);
                thetas_best=thetas_NS(i,ii);
                deltaCs_best=deltaCs_NS(i,ii);

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
                eigenvalue_allsolution(i,ii,:)=eig(Mstab);
            end
        end
    end

    %%% Plotting
    % subplot version
    for aa=1
        if alternative_situation==1
            x_plot=ji_alter;
            x_label='j';
            x_label_latex='$J^P\kappa_BT\Lambda^w/\Lambda^2$';
        elseif alternative_situation==2
            x_plot=jw_alter;
            x_label='jw';
            x_label_latex='$J^W/\Lambda$';
        elseif alternative_situation==3
            x_plot=gammajHat_alter;
            x_label='gammajHat';
            x_label_latex='$\gamma_j/\Gamma$';
        elseif alternative_situation==4
            x_plot=-SigmarHat_alter;
            x_label='SigmarHat';
            x_label_latex='$\sigma_r/L\Gamma$';
        elseif alternative_situation==5
            x_plot=jpHat_alter;
            x_label='jpHat';
            x_label_latex='$j_p\kappa_BT\Lambda^w/L^2\Lambda^2$';
        elseif alternative_situation==6
            x_plot=jpwHat_alter;
            x_label='jpwHat';
            x_label_latex='$j_p^w/L^2\Lambda$';
        elseif alternative_situation==7
            x_plot=gamma_alter;
            x_label='gamma';
            x_label_latex='$\Gamma\Lambda^w/L\Lambda$';
        end

        % 管腔半径
        subplot(3,1,1)
        Rs_NS_plot=Rs_NS.*R0/L;  % Temporary variable for plotting lumen radius
        for i=2:step-1
            for ii=1:3
                if thetas_NS(i,ii)~=0
                    if eigenvalue_allsolution(i,ii,1)<0 && eigenvalue_allsolution(i,ii,2)<0
                        if isreal(eigenvalue_allsolution(i,ii,:))
                            scatter(x_plot(i),Rs_NS_plot(i,ii),5,'filled','MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',0.75);
                            hold on
                        else
                            scatter(x_plot(i),Rs_NS_plot(i,ii),5,'filled','MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'LineWidth',0.75);
                            hold on
                        end
                    else
                        scatter(x_plot(i),Rs_NS_plot(i,ii),5,'filled','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',0.75);
                        hold on
                    end
                end
            end
        end

        set(gca,'FontName','Arial','FontSize',12);
        filename_plot=strcat("Tube geometry change with ",x_label," (consider leak)");
        %         title(filename_plot)
        hold off
        % 计算y坐标最优最大值
        for ii=1:1
            [~,I]=max(Rs_NS_plot);Rs_NS_plot(I)=mean(Rs_NS_plot);
        end
        axis ([min(x_plot) max(x_plot) 0 max(max(Rs_NS_plot))*1.1]);
        xlabel(x_label_latex,'interpreter','latex')
        ylabel('$R_s/L$','interpreter','latex')
        label_interval_temp=[1e-8 2e-8 5e-8 1e-7 2e-7 5e-7 1e-6 2e-6 5e-6 1e-5 2e-5 5e-5 1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2 2e-2 5e-2 1e-1 2e-1 5e-1];  % Candidate intervals for axis label values
        label_interval=[label_interval_temp label_interval_temp.*1e8];
        if alternative_situation==1
            [~,I]=min(abs(label_interval(:)-max(abs(min(x_plot)),abs(max(x_plot)))/10));x_interval=label_interval(I);  % x轴坐标标注间隔的选取
        elseif alternative_situation==3
            x_interval=0.4;
        else
            [~,I]=min(abs(label_interval(:)-max(abs(min(x_plot)),abs(max(x_plot)))/5));x_interval=label_interval(I);  % x轴坐标标注间隔的选取
        end
        [~,I]=min(abs(label_interval(:)-max(max(Rs_NS_plot))/5));y_interval=label_interval(I);  % y轴坐标标注间隔的选取
        if max(x_plot)>0
            set(gca,'XTick',0:x_interval:max(x_plot));  % Modify x-axis ticks
        else
            set(gca,'XTick',min(x_plot):x_interval:0);  % Modify x-axis ticks
        end
        set(gca,'YTick',0:y_interval:max(max(Rs_NS_plot))*1.1);  % 修改y坐标刻度
        set(gca,'FontName','Arial','FontSize',10)

        % 管腔角度
        subplot(3,1,2)
        for i=2:step-1
            for ii=1:3
                if thetas_NS(i,ii)~=0
                    if eigenvalue_allsolution(i,ii,1)<0 && eigenvalue_allsolution(i,ii,2)<0
                        if isreal(eigenvalue_allsolution(i,ii,:))
                            scatter(x_plot(i),thetas_NS(i,ii),5,'filled','MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',0.75);
                            hold on
                        else
                            scatter(x_plot(i),thetas_NS(i,ii),5,'filled','MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'LineWidth',0.75);
                            hold on
                        end
                    else
                        scatter(x_plot(i),thetas_NS(i,ii),5,'filled','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',0.75);
                        hold on
                    end
                end
            end
        end
        plot([min(x_plot) max(x_plot)],[pi/2 pi/2],'--','LineWidth',0.5,'Color','k')  % 画pi/2横线

        set(gca,'FontName','Arial','FontSize',12);
        hold off
        axis ([min(x_plot) max(x_plot) 0 pi]);
        xlabel(x_label_latex,'interpreter','latex')
        ylabel('$\theta_s$','interpreter','latex')
        label_interval_temp=[1e-8 2e-8 5e-8 1e-7 2e-7 5e-7 1e-6 2e-6 5e-6 1e-5 2e-5 5e-5 1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2 2e-2 5e-2 1e-1 2e-1 5e-1];  % Candidate intervals for axis label values
        label_interval=[label_interval_temp label_interval_temp.*1e8];
        if alternative_situation==1
            [~,I]=min(abs(label_interval(:)-max(abs(min(x_plot)),abs(max(x_plot)))/10));x_interval=label_interval(I);  % x轴坐标标注间隔的选取
        elseif alternative_situation==3
            x_interval=0.4;
        else
            [~,I]=min(abs(label_interval(:)-max(abs(min(x_plot)),abs(max(x_plot)))/5));x_interval=label_interval(I);  % x轴坐标标注间隔的选取
        end
        [~,I]=min(abs(label_interval(:)-pi/5));y_interval=label_interval(I);  % y轴坐标标注间隔的选取
        if max(x_plot)>0
            set(gca,'XTick',0:x_interval:max(x_plot));  % Modify x-axis ticks
        else
            set(gca,'XTick',min(x_plot):x_interval:0);  % Modify x-axis ticks
        end
        set(gca,'YTick',0:y_interval:pi);  % 修改y坐标刻度
        set(gca,'FontName','Arial','FontSize',10)

        % 管腔体积
        subplot(3,1,3)
        volume=2*pi/3.*Rs_NS_plot.^3.*(2+cos(thetas_NS)).*(1-cos(thetas_NS)).^2;
        volume_log=log10(volume);
        for i=2:step-1
            for ii=1:3
                % 若干计算错误的点排除
                if alternative_situation==1 && i==501
                    continue
                elseif alternative_situation==1 && i==744
                    volume_log(i,1)=volume_log(i-1,1);
                    continue
                elseif alternative_situation==1 && i==999
                    continue
                elseif alternative_situation==2 && i==101
                    continue
                elseif alternative_situation==2 && i==102
                    continue
                elseif alternative_situation==3 && i==383
                    continue
                elseif alternative_situation==3 && i==705
                    continue
                elseif alternative_situation==3 && i==772
                    continue
                elseif alternative_situation==4 && i==501
                    continue
                elseif alternative_situation==4 && i==635
                    continue
                elseif alternative_situation==5 && i==305
                    continue
                elseif alternative_situation==6 && i==29
                    continue
                elseif alternative_situation==6 && i==37
                    continue
                elseif alternative_situation==6 && i==322
                    continue
                end

                if thetas_NS(i,ii)~=0
                    if eigenvalue_allsolution(i,ii,1)<0 && eigenvalue_allsolution(i,ii,2)<0
                        if isreal(eigenvalue_allsolution(i,ii,:))
                            scatter(x_plot(i),volume_log(i,ii),5,'filled','MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',0.75);
                            hold on
                        else
                            scatter(x_plot(i),volume_log(i,ii),5,'filled','MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'LineWidth',0.75);
                            hold on
                        end
                    else
                        scatter(x_plot(i),volume_log(i,ii),5,'filled','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',0.75);
                        hold on
                    end
                end
            end
        end

        set(gca,'FontName','Arial','FontSize',12);
        hold off
        % 计算y坐标最优最大值
        for ii=1:1
            [~,I]=max(volume_log);volume_log(I)=mean(volume_log);
        end
        axis ([min(x_plot) max(x_plot) min(volume_log(~isinf(volume_log)))-1 max(volume_log(~isinf(volume_log)))+1]);
        xlabel(x_label_latex,'interpreter','latex')
        ylabel('volume (a.u.)','interpreter','latex')
        label_interval_temp=[1e-8 2e-8 5e-8 1e-7 2e-7 5e-7 1e-6 2e-6 5e-6 1e-5 2e-5 5e-5 1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2 2e-2 5e-2 1e-1 2e-1 5e-1];  % Candidate intervals for axis label values
        label_interval=[label_interval_temp label_interval_temp.*1e8];
        if alternative_situation==1
            [~,I]=min(abs(label_interval(:)-max(abs(min(x_plot)),abs(max(x_plot)))/10));x_interval=label_interval(I);  % x轴坐标标注间隔的选取
        elseif alternative_situation==3
            x_interval=0.4;
        else
            [~,I]=min(abs(label_interval(:)-max(abs(min(x_plot)),abs(max(x_plot)))/5));x_interval=label_interval(I);  % x轴坐标标注间隔的选取
        end
        [~,I]=min(abs(label_interval(:)-max(max(volume_log))/10));y_interval=label_interval(I);  % y轴坐标标注间隔的选取
        if max(x_plot)>0
            set(gca,'XTick',0:x_interval:max(x_plot));  % Modify x-axis ticks
        else
            set(gca,'XTick',min(x_plot):x_interval:0);  % Modify x-axis ticks
        end
        Y_labelpoint=log10([1e-15 1e-14 1e-13 1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1]);
        set(gca,'Ytick',Y_labelpoint)  %实际的值
        set(gca,'YtickLabel',{'10^{-15}','10^{-14}','10^{-13}','10^{-12}','10^{-11}','10^{-10}','10^{-9}','10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'});  %希望显示的值
        set(gca,'FontName','Arial','FontSize',10)
        set(gcf,'unit','centimeters','position',[0,0,8,22.5])  % Adjust image size
    end
end

