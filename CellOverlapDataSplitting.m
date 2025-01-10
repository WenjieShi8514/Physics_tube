function [data_1,data_2,data_3] = CellOverlapDataSplitting(data_1,data_2,data_3,C_AP,protein_number,data_changdubi,number_changdubi,j)

%%% Splitting Original Cell Data
%%% If the third file of the previous cell is the same as the first file of the current cell, the two cells are adjacent; process data_1.
%%% If the first file of the next cell is the same as the third file of the current cell, the two cells are adjacent; process data_3.
%%% Data Splitting Method:
%%% 1. Determine the adjacency of the current cell with the previous and next cells.
%%% 2. Calculate the fluorescence intensity ratio of the previous and next cells based on the average fluorescence intensity 
%%%    of the basal domain (data_2).
%%% 3. Use different strategies to split the original data based on various conditions:
%%%    (3.1) If there is no adjacent cell on the anterior side but there is one on the posterior side, it indicates this cell
%%%          is the first cell in the chain. Fit the anterior contractile ring distribution and determine the stage of the 
%%%          cell based on the lumen length/cell length. Use the functional relationship and confidence interval between 
%%%          rho_anterior and rho_posterior to provide the intensity of the posterior contractile ring.
%%%         (3.1.1) If the adjacent cell on the posterior side is the last cell in the chain, fit the posterior contractile 
%%%                 ring distribution and determine the stage of the posterior cell based on the lumen length/cell length.
%%%                 Use the functional relationship and confidence interval between rho_anterior and rho_posterior to give
%%%                 the intensity of the anterior contractile ring. Compare the rho_posterior of this cell and rho_anterior
%%%                 of the posterior cell with the observation results of the two contractile ring overlap distributions. 
%%%                 Assign the difference to the two rings so that the average rho_posterior/rho_anterior of these two cells
%%%                 is closest to C_AP. Traverse the posterior contractile ring range of this cell and the anterior contractile
%%%                 ring range of the next cell to find the combination that best matches the observed results.
%%%         (3.1.2) If the adjacent cell on the posterior side is not the last cell in the chain, it is not possible to fit 
%%%                 the posterior contractile ring distribution. Directly use rho_anterior of this cell to provide rho_posterior.
%%%                 Based on the observation results of the two overlapping contractile rings, determine the rho_anterior of 
%%%                 the next cell. Traverse the posterior contractile ring width range of this cell and the anterior contractile
%%%                 ring width range of the next cell to find the best match.
%%%    (3.2) If there are adjacent cells on both the anterior and posterior sides, this cell is a middle cell in the chain. Use 
%%%          the posterior contractile ring data of the previous cell to determine rho_anterior and w_anterior of this cell.
%%%         (3.2.1) If the adjacent cell on the posterior side is the last cell in the chain, fit the posterior contractile ring
%%%                 distribution and determine the stage of the posterior cell based on the lumen length/cell length. Directly 
%%%                 use rho_posterior of this cell to provide rho_anterior. Based on the rho_anterior of the next cell and the
%%%                 observed overlap results, determine rho_posterior of this cell. Traverse the posterior contractile ring 
%%%                 range of this cell and the anterior contractile ring range of the next cell to find the best match.
%%%         (3.2.2) If the adjacent cell on the posterior side is not the last cell in the chain, it is not possible to fit the
%%%                 posterior contractile ring distribution. Determine the stage of the next cell based on the 
%%%                 lumen length/cell length. Directly use rho_anterior of this cell to provide rho_posterior. Based on the 
%%%                 observed overlap results, determine rho_anterior of the next cell. Traverse the posterior contractile ring
%%%                 width range of this cell and the anterior contractile ring width range of the next cell to find the best match.
%%%    (3.3) If there is an adjacent cell on the anterior side but not on the posterior side, this cell is the last cell in the 
%%%          chain. Fit the posterior contractile ring distribution and determine the stage of the cell based on the 
%%%          lumen length/cell length. Use the functional relationship and sliding standard deviation between rho_anterior and
%%%          rho_posterior to provide the intensity of the anterior contractile ring.
%%%         (3.3.1) If the adjacent cell on the anterior side is the first cell in the chain, fit the anterior contractile ring
%%%                 distribution and determine the stage of the anterior cell based on the lumen length/cell length. Use the 
%%%                 functional relationship and confidence interval between rho_anterior and rho_posterior to give the intensity
%%%                 of the posterior contractile ring. Compare the rho_anterior of this cell and rho_posterior of the anterior 
%%%                 cell with the observation results of the two contractile ring overlap distributions. Assign the difference 
%%%                 to the two rings so that the average rho_posterior/rho_anterior of these two cells is closest to C_AP. 
%%%                 Traverse the anterior contractile ring range of this cell and the posterior contractile ring range of the 
%%%                 previous cell to find the best match.
%%%         (3.3.2) If the adjacent cell on the anterior side is not the first cell in the chain, it is not possible to fit the 
%%%                 anterior contractile ring distribution. Directly use rho_posterior of this cell to provide rho_anterior. 
%%%                 Based on the observed overlap results, determine rho_posterior of the previous cell. Traverse the anterior
%%%                 contractile ring width range of this cell and the posterior contractile ring width range of the previous 
%%%                 cell to find the best match.
%%%    (3.4) If the cells are not adjacent, skip the splitting step.

global framepath2 filename_0 filename_xiahuaxian filename_fanxiegang filename_1 filename_2 filename_3 filename_xlsx filename_csv filename_month filename_day filename1 filename2 filename3
global data_3_split1 data_3_split2 data_3_AP_split1 data_3_AP_split2

for bb=1
    %%% step 1
    % 拼接上一个细胞和下一个细胞的数据文件路径，和上上一个细胞和下下一个细胞的数据文件路径
    filename2_previous=num2str(j-1);
    filename_previous=[filename_month filename_day filename_xiahuaxian filename1 filename2_previous filename3];
    framepath_filename_previous1=[framepath2 filename_previous filename_xiahuaxian filename_1 filename_csv]; % 上一个细胞的第一个文件
    framepath_filename_previous2=[framepath2 filename_previous filename_xiahuaxian filename_2 filename_csv]; % 上一个细胞的第二个文件
    framepath_filename_previous3=[framepath2 filename_previous filename_xiahuaxian filename_3 filename_csv]; % 上一个细胞的第三个文件
    filename2_previousprevious=num2str(j-2);
    filename_previousprevious=[filename_month filename_day filename_xiahuaxian filename1 filename2_previousprevious filename3];
    framepath_filename_previousprevious1=[framepath2 filename_previousprevious filename_xiahuaxian filename_1 filename_csv]; % 上上一个细胞的第一个文件
    framepath_filename_previousprevious2=[framepath2 filename_previousprevious filename_xiahuaxian filename_2 filename_csv]; % 上上一个细胞的第二个文件
    framepath_filename_previousprevious3=[framepath2 filename_previousprevious filename_xiahuaxian filename_3 filename_csv]; % 上上一个细胞的第三个文件
    filename2_next=num2str(j+1);
    filename_next=[filename_month filename_day filename_xiahuaxian filename1 filename2_next filename3];
    framepath_filename_next1=[framepath2 filename_next filename_xiahuaxian filename_1 filename_csv]; % 下一个细胞的第一个文件
    framepath_filename_next2=[framepath2 filename_next filename_xiahuaxian filename_2 filename_csv]; % 下一个细胞的第二个文件
    framepath_filename_next3=[framepath2 filename_next filename_xiahuaxian filename_3 filename_csv]; % 下一个细胞的第三个文件
    filename2_nextnext=num2str(j+2);
    filename_nextnext=[filename_month filename_day filename_xiahuaxian filename1 filename2_nextnext filename3];
    framepath_filename_nextnext1=[framepath2 filename_nextnext filename_xiahuaxian filename_1 filename_csv]; % 下下一个细胞的第一个文件
    framepath_filename_nextnext2=[framepath2 filename_nextnext filename_xiahuaxian filename_2 filename_csv]; % 下下一个细胞的第二个文件
    framepath_filename_nextnext3=[framepath2 filename_nextnext filename_xiahuaxian filename_3 filename_csv]; % 下下一个细胞的第三个文件
    
    % 判断相邻情况
    neighbour_situation=0;  % 该细胞与前后细胞相邻情况的判定变量，值为1表示仅与后侧相邻，值为2表示与前后都相邻，值为3表示仅与前侧相邻，值为4表示该细胞不为相邻细胞
    if exist(framepath_filename_previous3,'file')  % 存在上一个属于相邻类型的细胞
        data_previous3=xlsread(framepath_filename_previous3);
        if length(data_previous3)==length(data_1) && ~sum(sum(abs(data_previous3-data_1)))  % 上一个细胞为相邻细胞
            if exist(framepath_filename_next1,'file')  % 存在下一个属于相邻类型的细胞
                data_next1=xlsread(framepath_filename_next1);
                if length(data_next1)==length(data_3) && ~sum(sum(abs(data_next1-data_3)))  % 下一个细胞为相邻细胞
                    neighbour_situation=2;  % 前后细胞均为相邻细胞，则属于情况2
                else
                    neighbour_situation=3;  % 上一个细胞为相邻细胞，下一个细胞为相邻类型的细胞但不相邻，则属于情况3
                end
            else
                neighbour_situation=3;  % 上一个细胞为相邻细胞，下一个细胞非相邻类型的细胞，则属于情况3
            end
        else
            if exist(framepath_filename_next1,'file')  % 存在下一个属于相邻类型的细胞
                data_next1=xlsread(framepath_filename_next1);
                if length(data_next1)==length(data_3) && ~sum(sum(abs(data_next1-data_3)))  % 下一个细胞为相邻细胞
                    neighbour_situation=1;  % 上一个细胞为相邻类型的细胞但不相邻，下一个细胞为相邻细胞，则属于情况1
                else
                    neighbour_situation=4;  % 上一个细胞为相邻类型的细胞但不相邻，下一个细胞为相邻类型的细胞但不相邻，则属于情况4
                end
            else
                neighbour_situation=4;  % 上一个细胞为相邻类型的细胞但不相邻，下一个细胞非相邻类型的细胞，则属于情况4
            end
        end
    else
        if exist(framepath_filename_next1,'file')  % 存在下一个属于相邻类型的细胞
            data_next1=xlsread(framepath_filename_next1);
            if length(data_next1)==length(data_3) && ~sum(sum(abs(data_next1-data_3)))  % 下一个细胞为相邻细胞
                neighbour_situation=1;  % 上一个细胞为非相邻类型的细胞，下一个细胞为相邻细胞，则属于情况1
            else
                neighbour_situation=4;  % 上一个细胞为非相邻类型的细胞，下一个细胞为相邻类型的细胞但不相邻，则属于情况4
            end
        else
            neighbour_situation=4;  % 上一个细胞为非相邻类型的细胞，下一个细胞非相邻类型的细胞，则属于情况4
        end
    end
    
    %%% step 2
    average_data_previous2=0;  % 前一个细胞的中间部分(basal domain)平均荧光强度
    average_data_2=0;  % 该细胞的中间部分(basal domain)平均荧光强度
    average_data_next2=0;  % 后一个细胞的中间部分(basal domain)平均荧光强度
    
    if neighbour_situation==1  % 仅与后侧相邻
        data_next2=xlsread(framepath_filename_next2);
        average_data_2=sum(data_2)/length(data_2);  % 计算该细胞的中间部分(basal domain)平均荧光强度
        average_data_next2=sum(data_next2)/length(data_next2);  % 计算后一个细胞的中间部分(basal domain)平均荧光强度
    elseif neighbour_situation==2  % 与前后都相邻
        data_previous2=xlsread(framepath_filename_previous2);
        data_next2=xlsread(framepath_filename_next2);
        average_data_previous2=sum(data_previous2)/length(data_previous2);  % 计算前一个细胞的中间部分(basal domain)平均荧光强度
        average_data_2=sum(data_2)/length(data_2);  % 计算该细胞的中间部分(basal domain)平均荧光强度
        average_data_next2=sum(data_next2)/length(data_next2);  % 计算后一个细胞的中间部分(basal domain)平均荧光强度
    elseif neighbour_situation==3  % 仅与前侧相邻
        data_previous2=xlsread(framepath_filename_previous2);
        average_data_previous2=sum(data_previous2)/length(data_previous2);  % 计算前一个细胞的中间部分(basal domain)平均荧光强度
        average_data_2=sum(data_2)/length(data_2);  % 计算该细胞的中间部分(basal domain)平均荧光强度
    elseif neighbour_situation==4  % 该细胞不为相邻细胞，则不做计算
    end
    
    
    
    %%% step 3
    % 先将该细胞的data_1, data_2和data_3的x轴归一化
    data_x_1=zeros(length(data_1),1);
    data_x_2=zeros(length(data_2),1);
    data_x_3=zeros(length(data_3),1);
    for ii=1:length(data_1)
        data_x_1(ii)=(data_1(ii,1)-(data_1(end,1)+data_2(end,1)+data_3(end,1))/2)/(data_1(end,1)+data_2(end,1)+data_3(end,1));
    end
    for ii=1:length(data_2)
        data_x_2(ii)=(data_1(end,1)+data_2(ii,1)-(data_1(end,1)+data_2(end,1)+data_3(end,1))/2)/(data_1(end,1)+data_2(end,1)+data_3(end,1));
    end
    for ii=1:length(data_3)
        data_x_3(ii)=(data_1(end,1)+data_2(end,1)+data_3(ii,1)-(data_1(end,1)+data_2(end,1)+data_3(end,1))/2)/(data_1(end,1)+data_2(end,1)+data_3(end,1));
    end
    data_1(:,1)=data_x_1;
    data_2(:,1)=data_x_2;
    data_3(:,1)=data_x_3;
    
    if neighbour_situation==1  % 仅与后侧相邻
        %%%  step (3.1)
        % 拟合前端收缩环分布，得到rho_interior
        for ab=1
            % 若该细胞所处时期为lumen形成早期，则在lateral domain全域拟合收缩环数学期望，若该细胞所处时期为lumen形成中后期，则在lateral domain后20%范围内拟合收缩环数学期望
            if data_changdubi(number_changdubi)<=0.2
                E_start=-0.2;
                E_end=0.2;
            else
                E_start=-0.1;
                E_end=0.1;
            end
            if data_changdubi(number_changdubi)<=0.2
                E_lateral_start=abs(data_1(end,1));
            else
                E_lateral_start=abs(data_1(end,1)+(data_1(1,1)-data_1(end,1))*4/5);
            end
            
            data_1and2=[data_1;data_2];
            x=data_1and2(:,1);  % 拟合函数的x轴
            y=data_1and2(:,2);  % 拟合函数的y轴
            wucha_best=100000000000;

            % 使用Levenberg-Marquardt method拟合
            for i_fit=1:10  % 循环10次，寻找最优拟合解
                ft=fittype('rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x+E_lateral)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_lateral)/w_posterior).^2)','independent','x','dependent','y');
                opts=fitoptions('Method','NonlinearLeastSquares');
                opts.Algorithm='Levenberg-Marquardt';
                opts.Display='Off';
                opts.MaxFunEvals=600;
                opts.MaxIter=400;
                opts.Lower=[E_start E_lateral_start 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
                opts.Upper=[E_end 0.5 max(y) max(y) max(y) max(y) 1 1 1];
                opts.StartPoint=[0 0.5 max(y)/5 max(y)/2 max(y)/2 max(y)/2 0.1 0.1 0.1];
                [xData,yData]=prepareCurveData(x,y);
                % Fit parameter result to variables
                [fitresult,gof]=fit(xData,yData,ft,opts);
                rho0=fitresult.rho0;
                rho_inf=fitresult.rho_inf;
                rho_inf_anterior=fitresult.rho_inf_anterior;
                rho_inf_posterior=fitresult.rho_inf_posterior;
                w=fitresult.w;
                w_anterior=fitresult.w_anterior;
                w_posterior=fitresult.w_posterior;
                E=fitresult.E;
                E_lateral=fitresult.E_lateral;
                
                wucha=0;
                for ii=1:length(x)
                    wucha=wucha+(rho0+rho_inf*exp(-1/2*((x(ii)-E)/w)^2)+rho_inf_anterior*exp(-1/2*((x(ii)+E_lateral)/w_anterior)^2)+...
                        rho_inf_posterior*exp(-1/2*((x(ii)-E_lateral)/w_posterior)^2)-data_1and2(ii,2))^2;
                end
                if wucha<wucha_best
                    wucha_best=wucha;
                    rho0_best=rho0;
                    rho_inf_best=rho_inf;
                    w_best=w;
                    E_best=E;
                    rho_inf_anterior_best=rho_inf_anterior;
                    w_anterior_best=w_anterior;
                    rho_inf_posterior_best=rho_inf_posterior;
                    w_posterior_best=w_posterior;
                    E_lateral_best=E_lateral;
                end
            end
        end
        
        % 根据lumen length/cell length判断该细胞所处时期，根据rho_anterior和rho_posterior的函数关系和置信区间给出后侧收缩环的强度
        rho_posterior_fit=rho_inf_anterior_best/C_AP;  % rho_posterior的回归均值
        
        % 判断下一个细胞的相邻情况
        neighbour_situation_next=0;  % 下一个细胞的相邻情况的判定变量，值为1表示不与下下一个细胞相邻(即为此细胞链的最后一个细胞)，值为2表示与下下一个细胞相邻
        data_next3=xlsread(framepath_filename_next3);
        if exist(framepath_filename_nextnext1,'file')  % 存在下下一个属于相邻类型的细胞
            data_nextnext1=xlsread(framepath_filename_nextnext1);
            if length(data_next3)==length(data_nextnext1) && ~sum(sum(abs(data_next3-data_nextnext1)))  % 下一个细胞与下下一个细胞为相邻细胞
                neighbour_situation_next=2;
            else
                neighbour_situation_next=1;
            end
        else
            neighbour_situation_next=1;
        end
        
        % 将下一个细胞的data_1和data_3的x轴归一化
        data_x_1=zeros(length(data_next1),1);
        data_x_2=zeros(length(data_next2),1);
        data_x_3=zeros(length(data_next3),1);
        for ii=1:length(data_next1)
            data_x_1(ii)=(data_next1(ii,1)-(data_next1(end,1)+data_next2(end,1)+data_next3(end,1))/2)/(data_next1(end,1)+data_next2(end,1)+data_next3(end,1));
        end
        for ii=1:length(data_next2)
            data_x_2(ii)=(data_next1(end,1)+data_next2(ii,1)-(data_next1(end,1)+data_next2(end,1)+data_next3(end,1))/2)/(data_next1(end,1)+data_next2(end,1)+data_next3(end,1));
        end
        for ii=1:length(data_next3)
            data_x_3(ii)=(data_next1(end,1)+data_next2(end,1)+data_next3(ii,1)-(data_next1(end,1)+data_next2(end,1)+data_next3(end,1))/2)/(data_next1(end,1)+data_next2(end,1)+data_next3(end,1));
        end
        data_next1(:,1)=data_x_1;
        data_next2(:,1)=data_x_2;
        data_next3(:,1)=data_x_3;
        
        if neighbour_situation_next==1
            %%% step (3.1.1)
            
            % 拟合后一个细胞后端收缩环分布，得到rho_posterior
            for ab=1
                % 若后一个细胞所处时期为lumen形成早期，则在lateral domain全域拟合收缩环数学期望，若后一个细胞所处时期为lumen形成中后期，则在lateral domain后20%范围内拟合收缩环数学期望
                if data_changdubi(number_changdubi+1)<=0.2
                    E_next_start=-0.2;
                    E_next_end=0.2;
                else
                    E_next_start=-0.1;
                    E_next_end=0.1;
                end
                if data_changdubi(number_changdubi+1)<=0.2
                    E_lateral_next_start=abs(data_next3(1,1));
                else
                    E_lateral_next_start=abs(data_next3(end,1)+(data_next3(1,1)-data_next3(end,1))*4.5/5);
                end

                data_next2and3=[data_next2;data_next3];
                x=data_next2and3(:,1);  % 拟合函数的x轴
                y=data_next2and3(:,2);  % 拟合函数的y轴
                wucha_best=100000000000;

                % 使用Levenberg-Marquardt method拟合
                for i_fit=1:10  % 循环10次，寻找最优拟合解
                    ft=fittype('rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x+E_lateral)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_lateral)/w_posterior).^2)','independent','x','dependent','y');
                    opts=fitoptions('Method','NonlinearLeastSquares');
                    opts.Algorithm='Levenberg-Marquardt';
                    opts.Display='Off';
                    opts.MaxFunEvals=600;
                    opts.MaxIter=400;
                    opts.Lower=[E_next_start E_lateral_next_start 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
                    opts.Upper=[E_next_end 0.5 max(y) max(y) max(y) max(y) 1 1 1];
                    opts.StartPoint=[0 0.5 max(y)/5 max(y)/2 max(y)/2 max(y)/2 0.1 0.1 0.1];
                    [xData,yData]=prepareCurveData(x,y);
                    % Fit parameter result to variables
                    [fitresult,gof]=fit(xData,yData,ft,opts);
                    rho0=fitresult.rho0;
                    rho_inf=fitresult.rho_inf;
                    rho_inf_anterior=fitresult.rho_inf_anterior;
                    rho_inf_posterior=fitresult.rho_inf_posterior;
                    w=fitresult.w;
                    w_anterior=fitresult.w_anterior;
                    w_posterior=fitresult.w_posterior;
                    E=fitresult.E;
                    E_lateral=fitresult.E_lateral;

                    wucha=0;
                    for ii=1:length(x)
                        wucha=wucha+(rho0+rho_inf*exp(-1/2*((x(ii)-E)/w)^2)+rho_inf_anterior*exp(-1/2*((x(ii)+E_lateral)/w_anterior)^2)+...
                            rho_inf_posterior*exp(-1/2*((x(ii)-E_lateral)/w_posterior)^2)-data_next2and3(ii,2))^2;
                    end
                    if wucha<wucha_best
                        wucha_best=wucha;
                        rho0_next_best=rho0;
                        rho_inf_best=rho_inf;
                        w_best=w;
                        E_best=E;
                        rho_inf_anterior_best=rho_inf_anterior;
                        w_anterior_best=w_anterior;
                        rho_inf_posterior_best=rho_inf_posterior;
                        w_posterior_best=w_posterior;
                        E_lateral_best=E_lateral;
                    end
                end
            end

            % 根据lumen length/cell length判断后一个细胞所处时期，根据rho_anterior和rho_posterior的函数关系和置信区间给出前侧收缩环的强度
            rho_anterior_next_fit=rho_inf_posterior_best*C_AP;  % rho_anterior的回归均值

            % 拟合该细胞posterior-lateral domain和后一个细胞anterior-lateral domain叠加观测的区域拟合高斯分布
            for ab=1
                E_lateral_start=data_3(end,1)+(data_3(1,1)-data_3(end,1))*1/5;
                rho0_start=max((rho0_best+rho0_next_best)*0.8,0.001);
                rho0_end=(rho0_best+rho0_next_best)*1.2;

                x=data_3(:,1);  % 拟合函数的x轴
                y=data_3(:,2);  % 拟合函数的y轴
                wucha_best=100000000000;

                % 使用Levenberg-Marquardt method拟合
                for i_fit=1:10  % 循环10次，寻找最优拟合解
                    ft=fittype('rho0_overlap+rho_inf_overlap*exp(-1/2*((x-E_lateral)/w_overlap).^2)','independent','x','dependent','y');
                    opts=fitoptions('Method','NonlinearLeastSquares');
                    opts.Algorithm='Levenberg-Marquardt';
                    opts.Display='Off';
                    opts.MaxFunEvals=600;
                    opts.MaxIter=400;
                    opts.Lower=[E_lateral_start rho0_start 0.01 0.01];
                    opts.Upper=[0.5 rho0_end max(y) 1];
                    opts.StartPoint=[0.5 (rho0_start+rho0_end)/2 max(y)/2 0.1];
                    [xData,yData]=prepareCurveData(x,y);
                    % Fit parameter result to variables
                    [fitresult,gof]=fit(xData,yData,ft,opts);
                    rho0_overlap=fitresult.rho0_overlap;
                    rho_inf_overlap=fitresult.rho_inf_overlap;
                    w_overlap=fitresult.w_overlap;
                    E_lateral=fitresult.E_lateral;

                    wucha=0;
                    for ii=1:length(x)
                        wucha=wucha+(rho0_overlap+rho_inf_overlap*exp(-1/2*((x(ii)-E_lateral)/w_overlap)^2)+-data_1and2(ii,2))^2;
                    end
                    if wucha<wucha_best
                        wucha_best=wucha;
                        rho0_overlap_best=rho0_overlap;
                        rho_inf_overlap_best=rho_inf_overlap;
                        w_overlap_best=w_overlap;
                        E_lateral_best=E_lateral;
                    end
                end
            end
            y=rho0_overlap_best+rho_inf_overlap_best*exp(-1/2*((x-E_lateral_best)./w_overlap_best).^2);
            
            % 根据该细胞的拟合rho0和rho_inf_posterior，和后一个细胞的rho0和rho_inf_anterior，与该细胞posterior-lateral domain和后一个细胞anterior-lateral domain叠加
            % 观测的区域数据的拟合函数(去噪后)极大值对比,将差值部分分配给两者，使得这两个细胞的平均rho_anterior/rho_posterior最接近C_AP
            rho_inf_chazhi=rho0_overlap_best+rho_inf_overlap_best-rho0_next_best-rho_anterior_next_fit-rho0_best-rho_posterior_fit;  % 叠加观测结果和前后环拟合结果和的差值
            wucha_best=1000;
            for rho_inf_chazhi_part1=rho_inf_chazhi/1000:rho_inf_chazhi/1000:rho_inf_chazhi  % rho_inf_chazhi_part1指将插值部分分给该细胞rho_posterior的部分
                rho_inf_chazhi_part2=rho_inf_chazhi-rho_inf_chazhi_part1;  % rho_inf_chazhi_part1指将插值部分分给后一个细胞rho_anterior的部分
                wucha=rho_inf_anterior_best/(rho_posterior_fit+rho_inf_chazhi_part1)+(rho_anterior_next_fit+rho_inf_chazhi_part2)/rho_inf_posterior_best...
                    -rho_anterior_next_fit/rho_inf_posterior_best-rho_inf_anterior_best/rho_posterior_fit;
                if wucha<wucha_best
                    wucha_best=wucha;
                    rho_inf_chazhi_part1_best=rho_inf_chazhi_part1;
                end
            end
            rho_posterior_fit=rho_posterior_fit+rho_inf_chazhi_part1_best;  % 将最优分配方法应用到前后两个高斯分布上
            rho_anterior_next_fit=rho_anterior_next_fit+rho_inf_chazhi-rho_inf_chazhi_part1_best;
            
            % 遍历该细胞的后侧收缩环范围，和后一个细胞的前侧收缩环范围，找出最符合这两个收缩环叠加分布观测结果的组合
            for ab=1
                % 若该细胞所处时期为lumen形成早期，则在lateral domain全域拟合收缩环数学期望，若该细胞所处时期为lumen形成中后期，则在lateral domain后20%范围内拟合收缩环数学期望
                if data_changdubi(number_changdubi)<=0.2
                    E_lateral_start=abs(data_3(1,1));
                else
                    E_lateral_start=abs(data_3(end,1)+(data_3(1,1)-data_3(end,1))*4/5);
                end
                % 若后一个细胞所处时期为lumen形成早期，则在lateral domain全域拟合收缩环数学期望，若后一个细胞所处时期为lumen形成中后期，则在lateral domain后20%范围内拟合收缩环数学期望
                if data_changdubi(number_changdubi+1)<=0.2
                    E_lateral_next_start=abs(data_3(1,1));
                else
                    E_lateral_next_start=abs(data_3(end,1)+(data_3(1,1)-data_3(end,1))*4/5);
                end

                rho0_start=max((rho0_best+rho0_next_best)*0.8,0.001);
                rho0_end=(rho0_best+rho0_next_best)*1.2;

                x=data_3(:,1);  % 拟合函数的x轴
                y=data_3(:,2);  % 拟合函数的y轴
                wucha_best=100000000000;

                % 使用Levenberg-Marquardt method拟合
                for i_fit=1:10  % 循环10次，寻找最优拟合解
                    ft=fittype('rho0+rho_inf_anterior_next*exp(-1/2*((x-E_lateral_next)/w_anterior_next).^2)+rho_inf_posterior*exp(-1/2*((x-E_lateral)/w_posterior).^2)','independent','x','dependent','y');
                    opts=fitoptions('Method','NonlinearLeastSquares');
                    opts.Algorithm='Levenberg-Marquardt';
                    opts.Display='Off';
                    opts.MaxFunEvals=600;
                    opts.MaxIter=400;
                    opts.Lower=[E_lateral_next_start E_lateral_start rho0_start rho_anterior_next_fit rho_posterior_fit 0.01 0.01];
                    opts.Upper=[0.5 0.5 rho0_end rho_anterior_next_fit rho_posterior_fit 1 1];
                    opts.StartPoint=[0.5 0.5 (rho0_start+rho0_end)/2 rho_anterior_next_fit rho_posterior_fit 0.1 0.1];
                    [xData,yData]=prepareCurveData(x,y);
                    % Fit parameter result to variables
                    [fitresult,gof]=fit(xData,yData,ft,opts);
                    rho0=fitresult.rho0;
                    w_anterior_next=fitresult.w_anterior_next;
                    w_posterior=fitresult.w_posterior;
                    E_lateral=fitresult.E_lateral;
                    E_lateral_next=fitresult.E_lateral_next;

                    wucha=0;
                    for ii=1:length(x)
                        wucha=wucha+(rho0+rho_posterior_fit*exp(-1/2*((x(ii)-E_lateral)/w_posterior)^2)+rho_anterior_next_fit*exp(-1/2*((x(ii)-E_lateral_next)/w_anterior_next)^2)-data_1and2(ii,2))^2;
                    end
                    if wucha<wucha_best
                        wucha_best=wucha;
                        rho0_overlap_best=rho0;
                        w_anterior_next_best=w_anterior_next;
                        w_posterior_best=w_posterior;
                        E_lateral_best=E_lateral;
                        E_lateral_next_best=E_lateral_next;
                    end
                end
            end
            
            % 将叠加分布观测结果进行拆分
            if protein_number==1
                data_3_split1=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给该细胞的部分
                data_3_split2=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给后一个细胞的部分
                rho0_split1=rho0_best/(rho0_best+rho0_next_best)*rho0_overlap_best;  % 该细胞的rho0和下一个细胞的rho0对比，取rho0_overlap的加权数值
                for ii=1:length(data_3)
                    residual_error=data_3(ii,2)-(rho0_overlap_best+rho_posterior_fit*exp(-1/2*((x(ii)-E_lateral_best)/w_posterior_best)^2)+rho_anterior_next_fit*exp(-1/2*((x(ii)-E_lateral_next_best)/w_anterior_next_best)^2));
                    data_3_split1(ii,1)=data_3(ii,1);
                    data_3_split1(ii,2)=rho0_split1+rho_posterior_fit*exp(-1/2*((data_3(ii,1)-E_lateral_best)/w_posterior_best)^2)+residual_error*average_data_2/(average_data_2+average_data_next2);
                    if data_3_split1(ii,2)<max(data_3(:,2))/200
                        data_3_split1(ii,2)=min(data_3(ii,2)/2,1);
                    end
                    if data_3_split1(ii,2)>data_3(ii,2)
                        data_3_split1(ii,2)=data_3(ii,2)*9/10;
                    end
                    data_3_split2(ii,1)=data_3(ii,1);
                    data_3_split2(ii,2)=data_3(ii,2)-data_3_split1(ii,2);
                end
            elseif protein_number==2
                data_3_AP_split1=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给该细胞的部分
                data_3_AP_split2=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给后一个细胞的部分
                rho0_split1=rho0_best/(rho0_best+rho0_next_best)*rho0_overlap_best;  % 该细胞的rho0和下一个细胞的rho0对比，取rho0_overlap的加权数值
                for ii=1:length(data_3)
                    residual_error=data_3(ii,2)-(rho0_overlap_best+rho_posterior_fit*exp(-1/2*((x(ii)-E_lateral_best)/w_posterior_best)^2)+rho_anterior_next_fit*exp(-1/2*((x(ii)-E_lateral_next_best)/w_anterior_next_best)^2));
                    data_3_AP_split1(ii,1)=data_3(ii,1);
                    data_3_AP_split1(ii,2)=rho0_split1+rho_posterior_fit*exp(-1/2*((data_3(ii,1)-E_lateral_best)/w_posterior_best)^2)+residual_error*average_data_2/(average_data_2+average_data_next2);
                    if data_3_AP_split1(ii,2)<max(data_3(:,2))/200
                        data_3_AP_split1(ii,2)=min(data_3(ii,2)/2,1);
                    end
                    if data_3_AP_split2(ii,2)>data_3(ii,2)
                        data_3_AP_split2(ii,2)=data_3(ii,2)*9/10;
                    end
                    data_3_AP_split2(ii,1)=data_3(ii,1);
                    data_3_AP_split2(ii,2)=data_3(ii,2)-data_3_AP_split1(ii,2);
                end
            end

        elseif neighbour_situation_next==2
            %%% step (3.1.2)
            % 直接根据该细胞的rho_anterior，以及rho_anterior和rho_posterior的函数关系和置信区间给出后侧收缩环的强度
            rho_posterior_fit=rho_inf_anterior_best/C_AP;  % rho_posterior的回归均值
            
            % 根据这两个收缩环叠加分布观测结果确定后一个细胞的rho_anterior
            % 拟合该细胞posterior-lateral domain和后一个细胞anterior-lateral domain叠加观测的区域拟合高斯分布
            for ab=1
                E_lateral_start=data_3(end,1)+(data_3(1,1)-data_3(end,1))*1/5;
                rho0_start=max(rho0_best,0.001);
                rho0_end=max(data_3(:,2));

                x=data_3(:,1);  % 拟合函数的x轴
                y=data_3(:,2);  % 拟合函数的y轴
                wucha_best=100000000000;

                % 使用Levenberg-Marquardt method拟合
                for i_fit=1:10  % 循环10次，寻找最优拟合解
                    ft=fittype('rho0_overlap+rho_inf_overlap*exp(-1/2*((x-E_lateral)/w_overlap).^2)','independent','x','dependent','y');
                    opts=fitoptions('Method','NonlinearLeastSquares');
                    opts.Algorithm='Levenberg-Marquardt';
                    opts.Display='Off';
                    opts.MaxFunEvals=600;
                    opts.MaxIter=400;
                    opts.Lower=[E_lateral_start rho0_start 0.01 0.01];
                    opts.Upper=[0.5 rho0_end max(y) 1];
                    opts.StartPoint=[0.5 (rho0_start+rho0_end)/2 max(y)/2 0.1];
                    [xData,yData]=prepareCurveData(x,y);
                    % Fit parameter result to variables
                    [fitresult,gof]=fit(xData,yData,ft,opts);
                    rho0_overlap=fitresult.rho0_overlap;
                    rho_inf_overlap=fitresult.rho_inf_overlap;
                    w_overlap=fitresult.w_overlap;
                    E_lateral=fitresult.E_lateral;

                    wucha=0;
                    for ii=1:length(x)
                        wucha=wucha+(rho0_overlap+rho_inf_overlap*exp(-1/2*((x(ii)-E_lateral)/w_overlap)^2)+-data_1and2(ii,2))^2;
                    end
                    if wucha<wucha_best
                        wucha_best=wucha;
                        rho0_overlap_best=rho0_overlap;
                        rho_inf_overlap_best=rho_inf_overlap;
                        w_overlap_best=w_overlap;
                        E_lateral_best=E_lateral;
                    end
                end
            end
            rho_anterior_next_fit=rho0_overlap_best+rho_inf_overlap_best-rho0_best-rho_posterior_fit;  % 后一个细胞的rho_anterior
            
            % 遍历该细胞的后侧收缩环宽度范围，和后一个细胞的前侧收缩环宽度范围，找出最符合这两个收缩环叠加分布观测结果的组合
            for ab=1
                % 若该细胞所处时期为lumen形成早期，则在lateral domain全域拟合收缩环数学期望，若该细胞所处时期为lumen形成中后期，则在lateral domain后20%范围内拟合收缩环数学期望
                if data_changdubi(number_changdubi)<=0.2
                    E_lateral_start=abs(data_3(1,1));
                else
                    E_lateral_start=abs(data_3(end,1)+(data_3(1,1)-data_3(end,1))*4/5);
                end
                % 若后一个细胞所处时期为lumen形成早期，则在lateral domain全域拟合收缩环数学期望，若后一个细胞所处时期为lumen形成中后期，则在lateral domain后20%范围内拟合收缩环数学期望
                if data_changdubi(number_changdubi+1)<=0.2
                    E_lateral_next_start=abs(data_3(1,1));
                else
                    E_lateral_next_start=abs(data_3(end,1)+(data_3(1,1)-data_3(end,1))*4/5);
                end

                rho0_start=max(rho0_best,0.001);
                rho0_end=max(data_3(:,2));

                x=data_3(:,1);  % 拟合函数的x轴
                y=data_3(:,2);  % 拟合函数的y轴
                wucha_best=100000000000;

                % 使用Levenberg-Marquardt method拟合
                for i_fit=1:10  % 循环10次，寻找最优拟合解
                    ft=fittype('rho0+rho_inf_anterior_next*exp(-1/2*((x-E_lateral_next)/w_anterior_next).^2)+rho_inf_posterior*exp(-1/2*((x-E_lateral)/w_posterior).^2)','independent','x','dependent','y');
                    opts=fitoptions('Method','NonlinearLeastSquares');
                    opts.Algorithm='Levenberg-Marquardt';
                    opts.Display='Off';
                    opts.MaxFunEvals=600;
                    opts.MaxIter=400;
                    opts.Lower=[E_lateral_next_start E_lateral_start rho0_start rho_anterior_next_fit rho_posterior_fit 0.01 0.01];
                    opts.Upper=[0.5 0.5 rho0_end rho_anterior_next_fit rho_posterior_fit 1 1];
                    opts.StartPoint=[0.5 0.5 (rho0_start+rho0_end)/2 rho_anterior_next_fit rho_posterior_fit 0.1 0.1];
                    [xData,yData]=prepareCurveData(x,y);
                    % Fit parameter result to variables
                    [fitresult,gof]=fit(xData,yData,ft,opts);
                    rho0=fitresult.rho0;
                    w_anterior_next=fitresult.w_anterior_next;
                    w_posterior=fitresult.w_posterior;
                    E_lateral=fitresult.E_lateral;
                    E_lateral_next=fitresult.E_lateral_next;

                    wucha=0;
                    for ii=1:length(x)
                        wucha=wucha+(rho0+rho_posterior_fit*exp(-1/2*((x(ii)-E_lateral)/w_posterior)^2)+rho_anterior_next_fit*exp(-1/2*((x(ii)-E_lateral_next)/w_anterior_next)^2)-data_1and2(ii,2))^2;
                    end
                    if wucha<wucha_best
                        wucha_best=wucha;
                        rho0_overlap_best=rho0;
                        w_anterior_next_best=w_anterior_next;
                        w_posterior_best=w_posterior;
                        E_lateral_best=E_lateral;
                        E_lateral_next_best=E_lateral_next;
                    end
                end
            end
            
            % 将叠加分布观测结果进行拆分
            if protein_number==1
                data_3_split1=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给该细胞的部分
                data_3_split2=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给后一个细胞的部分
                rho0_split1=rho0_best;  % 该细胞的rho0和下一个细胞的rho0对比，取rho0_overlap的加权数值
                for ii=1:length(data_3)
                    residual_error=data_3(ii,2)-(rho0_overlap_best+rho_posterior_fit*exp(-1/2*((x(ii)-E_lateral_best)/w_posterior_best)^2)+rho_anterior_next_fit*exp(-1/2*((x(ii)-E_lateral_next_best)/w_anterior_next_best)^2));
                    data_3_split1(ii,1)=data_3(ii,1);
                    data_3_split1(ii,2)=rho0_split1+rho_posterior_fit*exp(-1/2*((data_3(ii,1)-E_lateral_best)/w_posterior_best)^2)+residual_error*average_data_2/(average_data_2+average_data_next2);
                    if data_3_split1(ii,2)<max(data_3(:,2))/200
                        data_3_split1(ii,2)=min(data_3(ii,2)/2,1);
                    end
                    if data_3_split1(ii,2)>data_3(ii,2)
                        data_3_split1(ii,2)=data_3(ii,2)*9/10;
                    end
                    data_3_split2(ii,1)=data_3(ii,1);
                    data_3_split2(ii,2)=data_3(ii,2)-data_3_split1(ii,2);
                end
            elseif protein_number==2
                data_3_AP_split1=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给该细胞的部分
                data_3_AP_split2=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给后一个细胞的部分
                rho0_split1=rho0_best;  % 该细胞的rho0和下一个细胞的rho0对比，取rho0_overlap的加权数值
                for ii=1:length(data_3)
                    residual_error=data_3(ii,2)-(rho0_overlap_best+rho_posterior_fit*exp(-1/2*((x(ii)-E_lateral_best)/w_posterior_best)^2)+rho_anterior_next_fit*exp(-1/2*((x(ii)-E_lateral_next_best)/w_anterior_next_best)^2));
                    data_3_AP_split1(ii,1)=data_3(ii,1);
                    data_3_AP_split1(ii,2)=rho0_split1+rho_posterior_fit*exp(-1/2*((data_3(ii,1)-E_lateral_best)/w_posterior_best)^2)+residual_error*average_data_2/(average_data_2+average_data_next2);
                    if data_3_AP_split1(ii,2)<max(data_3(:,2))/200
                        data_3_AP_split1(ii,2)=min(data_3(ii,2)/2,1);
                    end
                    if data_3_AP_split1(ii,2)>data_3(ii,2)
                        data_3_AP_split1(ii,2)=data_3(ii,2)*9/10;
                    end
                    data_3_AP_split2(ii,1)=data_3(ii,1);
                    data_3_AP_split2(ii,2)=data_3(ii,2)-data_3_AP_split1(ii,2);
                end
            end
        end
        
        % 将拆分后的属于该细胞的部分重新赋值给该细胞的posterior-lateral domain部分
        if protein_number==1
            data_3=data_3_split1;
        elseif protein_number==2
            data_3=data_3_AP_split1;
        end
    elseif neighbour_situation==3  % 仅与前侧相邻
        %%% 由于上一个细胞已经拆分完成，因此无需根据原始方法再进行拟合、回归和拆分，查找前一个细胞拟合后侧收缩环的数据，直接使用前一个细胞的posterior-lateral
        %%% domain数据的拆分结果
        if protein_number==1
            data_1=data_3_split2; % 上一个细胞的posterior-lateral domain拆分后属于后一个细胞的部分赋值给该细胞的anterior-lateral domain
        elseif protein_number==2
            data_1=data_3_AP_split2; % 上一个细胞的posterior-lateral domain拆分后属于后一个细胞的部分赋值给该细胞的anterior-lateral domain
        end

        % 将该细胞第一个数据文件的数据进行颠倒
        data_1_fuben=data_1;
        for ii=1:length(data_1)
            data_1(ii,2)=data_1_fuben(length(data_1)-ii+1,2);
        end
    elseif neighbour_situation==2  % 与前后侧都相邻
        %%% step (3.2)
        if protein_number==1
            data_1=data_3_split2; % 上一个细胞的posterior-lateral domain拆分后属于后一个细胞的部分赋值给该细胞的anterior-lateral domain
        elseif protein_number==2
            data_1=data_3_AP_split2; % 上一个细胞的posterior-lateral domain拆分后属于后一个细胞的部分赋值给该细胞的anterior-lateral domain
        end
        
        % 将该细胞第一个数据文件的数据进行颠倒
        data_1_fuben=data_1;
        for ii=1:length(data_1)
            data_1(ii,2)=data_1_fuben(length(data_1)-ii+1,2);
        end
        
        % 将该细胞anterior-lateral domain拆分完成后的数据拟合rho0_best
        for ab=1
            % 若该细胞所处时期为lumen形成早期，则在lateral domain全域拟合收缩环数学期望，若该细胞所处时期为lumen形成中后期，则在lateral domain后20%范围内拟合收缩环数学期望
            if data_changdubi(number_changdubi)<=0.2
                E_start=-0.2;
                E_end=0.2;
            else
                E_start=-0.1;
                E_end=0.1;
            end
            if data_changdubi(number_changdubi)<=0.2
                E_lateral_start=abs(data_1(end,1));
            else
                E_lateral_start=abs(data_1(end,1)+(data_1(1,1)-data_1(end,1))*4/5);
            end
            
            data_1and2=[data_1;data_2];
            x=data_1and2(:,1);  % 拟合函数的x轴
            y=data_1and2(:,2);  % 拟合函数的y轴
            wucha_best=100000000000;

            % 使用Levenberg-Marquardt method拟合
            for i_fit=1:10  % 循环10次，寻找最优拟合解
                ft=fittype('rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x+E_lateral)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_lateral)/w_posterior).^2)','independent','x','dependent','y');
                opts=fitoptions('Method','NonlinearLeastSquares');
                opts.Algorithm='Levenberg-Marquardt';
                opts.Display='Off';
                opts.MaxFunEvals=600;
                opts.MaxIter=400;
                opts.Lower=[E_start E_lateral_start 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
                opts.Upper=[E_end 0.5 max(y) max(y) max(y) max(y) 1 1 1];
                opts.StartPoint=[0 0.5 max(y)/5 max(y)/2 max(y)/2 max(y)/2 0.1 0.1 0.1];
                [xData,yData]=prepareCurveData(x,y);
                % Fit parameter result to variables
                [fitresult,gof]=fit(xData,yData,ft,opts);
                rho0=fitresult.rho0;
                rho_inf=fitresult.rho_inf;
                rho_inf_anterior=fitresult.rho_inf_anterior;
                rho_inf_posterior=fitresult.rho_inf_posterior;
                w=fitresult.w;
                w_anterior=fitresult.w_anterior;
                w_posterior=fitresult.w_posterior;
                E=fitresult.E;
                E_lateral=fitresult.E_lateral;
                
                wucha=0;
                for ii=1:length(x)
                    wucha=wucha+(rho0+rho_inf*exp(-1/2*((x(ii)-E)/w)^2)+rho_inf_anterior*exp(-1/2*((x(ii)+E_lateral)/w_anterior)^2)+...
                        rho_inf_posterior*exp(-1/2*((x(ii)-E_lateral)/w_posterior)^2)-data_1and2(ii,2))^2;
                end
                if wucha<wucha_best
                    wucha_best=wucha;
                    rho0_best=rho0;
                    rho_inf_best=rho_inf;
                    w_best=w;
                    E_best=E;
                    rho_inf_anterior_best=rho_inf_anterior;
                    w_anterior_best=w_anterior;
                    rho_inf_posterior_best=rho_inf_posterior;
                    w_posterior_best=w_posterior;
                    E_lateral_best=E_lateral;
                end
            end
        end

        % 判断下一个细胞的相邻情况
        neighbour_situation_next=0;  % 下一个细胞的相邻情况的判定变量，值为1表示不与下下一个细胞相邻(即为此细胞链的最后一个细胞)，值为2表示与下下一个细胞相邻
        data_next3=xlsread(framepath_filename_next3);
        if exist(framepath_filename_nextnext1,'file')  % 存在下下一个属于相邻类型的细胞
            data_nextnext1=xlsread(framepath_filename_nextnext1);
            if length(data_next3)==length(data_nextnext1) && ~sum(sum(abs(data_next3-data_nextnext1)))  % 下一个细胞与下下一个细胞为相邻细胞
                neighbour_situation_next=2;
            else
                neighbour_situation_next=1;
            end
        else
            neighbour_situation_next=1;
        end

        % 将下一个细胞的data_1和data_3的x轴归一化
        data_x_1=zeros(length(data_next1),1);
        data_x_2=zeros(length(data_next2),1);
        data_x_3=zeros(length(data_next3),1);
        for ii=1:length(data_next1)
            data_x_1(ii)=(data_next1(ii,1)-(data_next1(end,1)+data_next2(end,1)+data_next3(end,1))/2)/(data_next1(end,1)+data_next2(end,1)+data_next3(end,1));
        end
        for ii=1:length(data_next2)
            data_x_2(ii)=(data_next1(end,1)+data_next2(ii,1)-(data_next1(end,1)+data_next2(end,1)+data_next3(end,1))/2)/(data_next1(end,1)+data_next2(end,1)+data_next3(end,1));
        end
        for ii=1:length(data_next3)
            data_x_3(ii)=(data_next1(end,1)+data_next2(end,1)+data_next3(ii,1)-(data_next1(end,1)+data_next2(end,1)+data_next3(end,1))/2)/(data_next1(end,1)+data_next2(end,1)+data_next3(end,1));
        end
        data_next1(:,1)=data_x_1;
        data_next2(:,1)=data_x_2;
        data_next3(:,1)=data_x_3;

        if neighbour_situation_next==1
            %%% step (3.2.1)
            % 拟合后一个细胞后端收缩环分布，得到rho_posterior
            for ab=1
                % 若后一个细胞所处时期为lumen形成早期，则在lateral domain全域拟合收缩环数学期望，若后一个细胞所处时期为lumen形成中后期，则在lateral domain后20%范围内拟合收缩环数学期望
                if data_changdubi(number_changdubi+1)<=0.2
                    E_next_start=-0.2;
                    E_next_end=0.2;
                else
                    E_next_start=-0.1;
                    E_next_end=0.1;
                end
                if data_changdubi(number_changdubi+1)<=0.2
                    E_lateral_next_start=abs(data_next3(1,1));
                else
                    E_lateral_next_start=abs(data_next3(end,1)+(data_next3(1,1)-data_next3(end,1))*4.5/5);
                end

                data_next2and3=[data_next2;data_next3];
                x=data_next2and3(:,1);  % 拟合函数的x轴
                y=data_next2and3(:,2);  % 拟合函数的y轴
                wucha_best=100000000000;

                % 使用Levenberg-Marquardt method拟合
                for i_fit=1:10  % 循环10次，寻找最优拟合解
                    ft=fittype('rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x+E_lateral)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_lateral)/w_posterior).^2)','independent','x','dependent','y');
                    opts=fitoptions('Method','NonlinearLeastSquares');
                    opts.Algorithm='Levenberg-Marquardt';
                    opts.Display='Off';
                    opts.MaxFunEvals=600;
                    opts.MaxIter=400;
                    opts.Lower=[E_next_start E_lateral_next_start 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
                    opts.Upper=[E_next_end 0.5 max(y) max(y) max(y) max(y) 1 1 1];
                    opts.StartPoint=[0 0.5 max(y)/5 max(y)/2 max(y)/2 max(y)/2 0.1 0.1 0.1];
                    [xData,yData]=prepareCurveData(x,y);
                    % Fit parameter result to variables
                    [fitresult,gof]=fit(xData,yData,ft,opts);
                    rho0=fitresult.rho0;
                    rho_inf=fitresult.rho_inf;
                    rho_inf_anterior=fitresult.rho_inf_anterior;
                    rho_inf_posterior=fitresult.rho_inf_posterior;
                    w=fitresult.w;
                    w_anterior=fitresult.w_anterior;
                    w_posterior=fitresult.w_posterior;
                    E=fitresult.E;
                    E_lateral=fitresult.E_lateral;

                    wucha=0;
                    for ii=1:length(x)
                        wucha=wucha+(rho0+rho_inf*exp(-1/2*((x(ii)-E)/w)^2)+rho_inf_anterior*exp(-1/2*((x(ii)+E_lateral_best)/w_anterior)^2)+...
                            rho_inf_posterior*exp(-1/2*((x(ii)-E_lateral_best)/w_posterior)^2)-data_next2and3(ii,2))^2;
                    end
                    if wucha<wucha_best
                        wucha_best=wucha;
                        rho0_next_best=rho0;
                        rho_inf_best=rho_inf;
                        w_best=w;
                        E_best=E;
                        rho_inf_anterior_best=rho_inf_anterior;
                        w_anterior_best=w_anterior;
                        rho_inf_posterior_best=rho_inf_posterior;
                        w_posterior_best=w_posterior;
                        E_lateral_best=E_lateral;
                    end
                end
            end

            % 根据lumen length/cell length判断后一个细胞所处时期，直接根据该细胞的rho_posterior给出该细胞的rho_anterior
            rho_anterior_next_fit=rho_inf_posterior_best*C_AP;

            % 根据后一个细胞的rho_anterior和这两个收缩环叠加分布观测结果确定该细胞的rho_posterior
            % 拟合该细胞posterior-lateral domain和后一个细胞anterior-lateral domain叠加观测的区域拟合高斯分布
            for ab=1
                E_lateral_start=data_3(end,1)+(data_3(1,1)-data_3(end,1))*1/5;
                rho0_start=max((rho0_best+rho0_next_best)*0.8,0.001);
                rho0_end=(rho0_best+rho0_next_best)*1.2;

                x=data_3(:,1);  % 拟合函数的x轴
                y=data_3(:,2);  % 拟合函数的y轴
                wucha_best=100000000000;

                % 使用Levenberg-Marquardt method拟合
                for i_fit=1:10  % 循环10次，寻找最优拟合解
                    ft=fittype('rho0_overlap+rho_inf_overlap*exp(-1/2*((x-E_lateral)/w_overlap).^2)','independent','x','dependent','y');
                    opts=fitoptions('Method','NonlinearLeastSquares');
                    opts.Algorithm='Levenberg-Marquardt';
                    opts.Display='Off';
                    opts.MaxFunEvals=600;
                    opts.MaxIter=400;
                    opts.Lower=[E_lateral_start rho0_start 0.01 0.01];
                    opts.Upper=[0.5 rho0_end max(y) 1];
                    opts.StartPoint=[0.5 (rho0_start+rho0_end)/2 max(y)/2 0.1];
                    [xData,yData]=prepareCurveData(x,y);
                    % Fit parameter result to variables
                    [fitresult,gof]=fit(xData,yData,ft,opts);
                    rho0_overlap=fitresult.rho0_overlap;
                    rho_inf_overlap=fitresult.rho_inf_overlap;
                    w_overlap=fitresult.w_overlap;
                    E_lateral=fitresult.E_lateral;

                    wucha=0;
                    for ii=1:length(x)
                        wucha=wucha+(rho0_overlap+rho_inf_overlap*exp(-1/2*((x(ii)-E_lateral)/w_overlap)^2)+-data_1and2(ii,2))^2;
                    end
                    if wucha<wucha_best
                        wucha_best=wucha;
                        rho0_overlap_best=rho0_overlap;
                        rho_inf_overlap_best=rho_inf_overlap;
                        w_overlap_best=w_overlap;
                        E_lateral_best=E_lateral;
                    end
                end
            end
            rho_posterior_fit=rho0_overlap_best+rho_inf_overlap_best-rho0_best-rho_anterior_next_fit-rho0_next_best;  % 该细胞的rho_posterior
            
            % 遍历该细胞的后侧收缩环宽度范围，和后一个细胞的前侧收缩环宽度范围，找出最符合这两个收缩环叠加分布观测结果的组合
            for ab=1
                % 若该细胞所处时期为lumen形成早期，则在lateral domain全域拟合收缩环数学期望，若该细胞所处时期为lumen形成中后期，则在lateral domain后20%范围内拟合收缩环数学期望
                if data_changdubi(number_changdubi)<=0.2
                    E_lateral_start=abs(data_3(1,1));
                else
                    E_lateral_start=abs(data_3(end,1)+(data_3(1,1)-data_3(end,1))*4/5);
                end
                % 若后一个细胞所处时期为lumen形成早期，则在lateral domain全域拟合收缩环数学期望，若后一个细胞所处时期为lumen形成中后期，则在lateral domain后20%范围内拟合收缩环数学期望
                if data_changdubi(number_changdubi+1)<=0.2
                    E_lateral_next_start=abs(data_3(1,1));
                else
                    E_lateral_next_start=abs(data_3(end,1)+(data_3(1,1)-data_3(end,1))*4/5);
                end

                rho0_start=max((rho0_best+rho0_next_best)*0.8,0.001);
                rho0_end=(rho0_best+rho0_next_best)*1.2;

                x=data_3(:,1);  % 拟合函数的x轴
                y=data_3(:,2);  % 拟合函数的y轴
                wucha_best=100000000000;

                % 使用Levenberg-Marquardt method拟合
                for i_fit=1:10  % 循环10次，寻找最优拟合解
                    ft=fittype('rho0+rho_inf_anterior_next*exp(-1/2*((x-E_lateral_next)/w_anterior_next).^2)+rho_inf_posterior*exp(-1/2*((x-E_lateral)/w_posterior).^2)','independent','x','dependent','y');
                    opts=fitoptions('Method','NonlinearLeastSquares');
                    opts.Algorithm='Levenberg-Marquardt';
                    opts.Display='Off';
                    opts.MaxFunEvals=600;
                    opts.MaxIter=400;
                    opts.Lower=[E_lateral_next_start E_lateral_start rho0_start rho_anterior_next_fit rho_posterior_fit 0.01 0.01];
                    opts.Upper=[0.5 0.5 rho0_end rho_anterior_next_fit rho_posterior_fit 1 1];
                    opts.StartPoint=[0.5 0.5 (rho0_start+rho0_end)/2 rho_anterior_next_fit rho_posterior_fit 0.1 0.1];
                    [xData,yData]=prepareCurveData(x,y);
                    % Fit parameter result to variables
                    [fitresult,gof]=fit(xData,yData,ft,opts);
                    rho0=fitresult.rho0;
                    w_anterior_next=fitresult.w_anterior_next;
                    w_posterior=fitresult.w_posterior;
                    E_lateral=fitresult.E_lateral;
                    E_lateral_next=fitresult.E_lateral_next;

                    wucha=0;
                    for ii=1:length(x)
                        wucha=wucha+(rho0+rho_posterior_fit*exp(-1/2*((x(ii)-E_lateral)/w_posterior)^2)+rho_anterior_next_fit*exp(-1/2*((x(ii)-E_lateral_next)/w_anterior_next)^2)-data_1and2(ii,2))^2;
                    end
                    if wucha<wucha_best
                        wucha_best=wucha;
                        rho0_overlap_best=rho0;
                        w_anterior_next_best=w_anterior_next;
                        w_posterior_best=w_posterior;
                        E_lateral_best=E_lateral;
                        E_lateral_next_best=E_lateral_next;
                    end
                end
            end
            
            % 将叠加分布观测结果进行拆分
            if protein_number==1
                data_3_split1=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给该细胞的部分
                data_3_split2=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给后一个细胞的部分
                rho0_split1=rho0_best/(rho0_best+rho0_next_best)*rho0_overlap_best;  % 该细胞的rho0和下一个细胞的rho0对比，取rho0_overlap的加权数值
                for ii=1:length(data_3)
                    residual_error=data_3(ii,2)-(rho0_overlap_best+rho_posterior_fit*exp(-1/2*((x(ii)-E_lateral_best)/w_posterior_best)^2)+rho_anterior_next_fit*exp(-1/2*((x(ii)-E_lateral_next_best)/w_anterior_next_best)^2));
                    data_3_split1(ii,1)=data_3(ii,1);
                    data_3_split1(ii,2)=rho0_split1+rho_posterior_fit*exp(-1/2*((data_3(ii,1)-E_lateral_best)/w_posterior_best)^2)+residual_error*average_data_2/(average_data_2+average_data_next2);
                    if data_3_split1(ii,2)<max(data_3(:,2))/200
                        data_3_split1(ii,2)=min(data_3(ii,2)/2,1);
                    end
                    if data_3_split1(ii,2)>data_3(ii,2)
                        data_3_split1(ii,2)=data_3(ii,2)*9/10;
                    end
                    data_3_split2(ii,1)=data_3(ii,1);
                    data_3_split2(ii,2)=data_3(ii,2)-data_3_split1(ii,2);
                end
            elseif protein_number==2
                data_3_AP_split1=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给该细胞的部分
                data_3_AP_split2=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给后一个细胞的部分
                rho0_split1=rho0_best/(rho0_best+rho0_next_best)*rho0_overlap_best;  % 该细胞的rho0和下一个细胞的rho0对比，取rho0_overlap的加权数值
                for ii=1:length(data_3)
                    residual_error=data_3(ii,2)-(rho0_overlap_best+rho_posterior_fit*exp(-1/2*((x(ii)-E_lateral_best)/w_posterior_best)^2)+rho_anterior_next_fit*exp(-1/2*((x(ii)-E_lateral_next_best)/w_anterior_next_best)^2));
                    data_3_AP_split1(ii,1)=data_3(ii,1);
                    data_3_AP_split1(ii,2)=rho0_split1+rho_posterior_fit*exp(-1/2*((data_3(ii,1)-E_lateral_best)/w_posterior_best)^2)+residual_error*average_data_2/(average_data_2+average_data_next2);
                    if data_3_AP_split1(ii,2)<max(data_3(:,2))/200
                        data_3_AP_split1(ii,2)=min(data_3(ii,2)/2,1);
                    end
                    if data_3_AP_split1(ii,2)>data_3(ii,2)
                        data_3_AP_split1(ii,2)=data_3(ii,2)*9/10;
                    end
                    data_3_AP_split2(ii,1)=data_3(ii,1);
                    data_3_AP_split2(ii,2)=data_3(ii,2)-data_3_AP_split1(ii,2);
                end
            end
        elseif neighbour_situation_next==2
            %%% step (3.2.2)
            % 根据lumen length/cell length判断后一个细胞所处时期，直接根据该细胞的前侧收缩环rho_anterior给出该细胞的rho_posterior
            rho_posterior_fit=rho_inf_anterior_best/C_AP;
            
            % 根据这两个收缩环叠加分布观测结果确定后一个细胞的rho_anterior
            % 拟合该细胞posterior-lateral domain和后一个细胞anterior-lateral domain叠加观测的区域拟合高斯分布
            for ab=1
                E_lateral_start=data_3(end,1)+(data_3(1,1)-data_3(end,1))*1/5;
                rho0_start=max(rho0_best,0.001);
                rho0_end=max(data_3(:,2));

                x=data_3(:,1);  % 拟合函数的x轴
                y=data_3(:,2);  % 拟合函数的y轴
                wucha_best=100000000000;

                % 使用Levenberg-Marquardt method拟合
                for i_fit=1:10  % 循环10次，寻找最优拟合解
                    ft=fittype('rho0_overlap+rho_inf_overlap*exp(-1/2*((x-E_lateral)/w_overlap).^2)','independent','x','dependent','y');
                    opts=fitoptions('Method','NonlinearLeastSquares');
                    opts.Algorithm='Levenberg-Marquardt';
                    opts.Display='Off';
                    opts.MaxFunEvals=600;
                    opts.MaxIter=400;
                    opts.Lower=[E_lateral_start rho0_start 0.01 0.01];
                    opts.Upper=[0.5 rho0_end max(y) 1];
                    opts.StartPoint=[0.5 (rho0_start+rho0_end)/2 max(y)/2 0.1];
                    [xData,yData]=prepareCurveData(x,y);
                    % Fit parameter result to variables
                    [fitresult,gof]=fit(xData,yData,ft,opts);
                    rho0_overlap=fitresult.rho0_overlap;
                    rho_inf_overlap=fitresult.rho_inf_overlap;
                    w_overlap=fitresult.w_overlap;
                    E_lateral=fitresult.E_lateral;

                    wucha=0;
                    for ii=1:length(x)
                        wucha=wucha+(rho0_overlap+rho_inf_overlap*exp(-1/2*((x(ii)-E_lateral)/w_overlap)^2)+-data_1and2(ii,2))^2;
                    end
                    if wucha<wucha_best
                        wucha_best=wucha;
                        rho0_overlap_best=rho0_overlap;
                        rho_inf_overlap_best=rho_inf_overlap;
                        w_overlap_best=w_overlap;
                        E_lateral_best=E_lateral;
                    end
                end
            end
            rho_anterior_next_fit=rho0_overlap_best+rho_inf_overlap_best-rho0_best-rho_posterior_fit;  % 后一个细胞的rho_anterior
            
            % 遍历该细胞的后侧收缩环宽度范围，和后一个细胞的前侧收缩环宽度范围，找出最符合这两个收缩环叠加分布观测结果的组合
            for ab=1
                % 若该细胞所处时期为lumen形成早期，则在lateral domain全域拟合收缩环数学期望，若该细胞所处时期为lumen形成中后期，则在lateral domain后20%范围内拟合收缩环数学期望
                if data_changdubi(number_changdubi)<=0.2
                    E_lateral_start=abs(data_3(1,1));
                else
                    E_lateral_start=abs(data_3(end,1)+(data_3(1,1)-data_3(end,1))*4/5);
                end
                % 若后一个细胞所处时期为lumen形成早期，则在lateral domain全域拟合收缩环数学期望，若后一个细胞所处时期为lumen形成中后期，则在lateral domain后20%范围内拟合收缩环数学期望
                if data_changdubi(number_changdubi+1)<=0.2
                    E_lateral_next_start=abs(data_3(1,1));
                else
                    E_lateral_next_start=abs(data_3(end,1)+(data_3(1,1)-data_3(end,1))*4/5);
                end

                rho0_start=max(rho0_best,0.001);
                rho0_end=max(data_3(:,2));

                x=data_3(:,1);  % 拟合函数的x轴
                y=data_3(:,2);  % 拟合函数的y轴
                wucha_best=100000000000;

                % 使用Levenberg-Marquardt method拟合
                for i_fit=1:10  % 循环10次，寻找最优拟合解
                    ft=fittype('rho0+rho_inf_anterior_next*exp(-1/2*((x-E_lateral_next)/w_anterior_next).^2)+rho_inf_posterior*exp(-1/2*((x-E_lateral)/w_posterior).^2)','independent','x','dependent','y');
                    opts=fitoptions('Method','NonlinearLeastSquares');
                    opts.Algorithm='Levenberg-Marquardt';
                    opts.Display='Off';
                    opts.MaxFunEvals=600;
                    opts.MaxIter=400;
                    opts.Lower=[E_lateral_next_start E_lateral_start rho0_start rho_anterior_next_fit rho_posterior_fit 0.01 0.01];
                    opts.Upper=[0.5 0.5 rho0_end rho_anterior_next_fit rho_posterior_fit 1 1];
                    opts.StartPoint=[0.5 0.5 (rho0_start+rho0_end)/2 rho_anterior_next_fit rho_posterior_fit 0.1 0.1];
                    [xData,yData]=prepareCurveData(x,y);
                    % Fit parameter result to variables
                    [fitresult,gof]=fit(xData,yData,ft,opts);
                    rho0=fitresult.rho0;
                    w_anterior_next=fitresult.w_anterior_next;
                    w_posterior=fitresult.w_posterior;
                    E_lateral=fitresult.E_lateral;
                    E_lateral_next=fitresult.E_lateral_next;

                    wucha=0;
                    for ii=1:length(x)
                        wucha=wucha+(rho0+rho_posterior_fit*exp(-1/2*((x(ii)-E_lateral)/w_posterior)^2)+rho_anterior_next_fit*exp(-1/2*((x(ii)-E_lateral_next)/w_anterior_next)^2)-data_1and2(ii,2))^2;
                    end
                    if wucha<wucha_best
                        wucha_best=wucha;
                        rho0_overlap_best=rho0;
                        w_anterior_next_best=w_anterior_next;
                        w_posterior_best=w_posterior;
                        E_lateral_best=E_lateral;
                        E_lateral_next_best=E_lateral_next;
                    end
                end
            end
            
            % 将叠加分布观测结果进行拆分
            if protein_number==1
                data_3_split1=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给该细胞的部分
                data_3_split2=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给后一个细胞的部分
                rho0_split1=rho0_best;  % 该细胞的rho0和下一个细胞的rho0对比，取rho0_overlap的加权数值
                for ii=1:length(data_3)
                    residual_error=data_3(ii,2)-(rho0_overlap_best+rho_posterior_fit*exp(-1/2*((x(ii)-E_lateral_best)/w_posterior_best)^2)+rho_anterior_next_fit*exp(-1/2*((x(ii)-E_lateral_next_best)/w_anterior_next_best)^2));
                    data_3_split1(ii,1)=data_3(ii,1);
                    data_3_split1(ii,2)=rho0_split1+rho_posterior_fit*exp(-1/2*((data_3(ii,1)-E_lateral_best)/w_posterior_best)^2)+residual_error*average_data_2/(average_data_2+average_data_next2);
                    if data_3_split1(ii,2)<max(data_3(:,2))/200
                        data_3_split1(ii,2)=min(data_3(ii,2)/2,1);
                    end
                    if data_3_split1(ii,2)>data_3(ii,2)
                        data_3_split1(ii,2)=data_3(ii,2)*9/10;
                    end
                    data_3_split2(ii,1)=data_3(ii,1);
                    data_3_split2(ii,2)=data_3(ii,2)-data_3_split1(ii,2);
                end
            elseif protein_number==2
                data_3_AP_split1=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给该细胞的部分
                data_3_AP_split2=zeros(size(data_3));  % 该细胞posterior-lateral domain数据 (data_3)拆分给后一个细胞的部分
                rho0_split1=rho0_best;  % 该细胞的rho0和下一个细胞的rho0对比，取rho0_overlap的加权数值
                for ii=1:length(data_3)
                    residual_error=data_3(ii,2)-(rho0_overlap_best+rho_posterior_fit*exp(-1/2*((x(ii)-E_lateral_best)/w_posterior_best)^2)+rho_anterior_next_fit*exp(-1/2*((x(ii)-E_lateral_next_best)/w_anterior_next_best)^2));
                    data_3_AP_split1(ii,1)=data_3(ii,1);
                    data_3_AP_split1(ii,2)=rho0_split1+rho_posterior_fit*exp(-1/2*((data_3(ii,1)-E_lateral_best)/w_posterior_best)^2)+residual_error*average_data_2/(average_data_2+average_data_next2);
                    if data_3_AP_split1(ii,2)<max(data_3(:,2))/200
                        data_3_AP_split1(ii,2)=min(data_3(ii,2)/2,1);
                    end
                    if data_3_AP_split2(ii,2)>data_3(ii,2)
                        data_3_AP_split2(ii,2)=data_3(ii,2)*9/10;
                    end
                    data_3_AP_split2(ii,1)=data_3(ii,1);
                    data_3_AP_split2(ii,2)=data_3(ii,2)-data_3_AP_split1(ii,2);
                end
            end
        end

        % 将拆分后的属于该细胞的部分重新赋值给该细胞的posterior-lateral domain部分
        if protein_number==1
            data_3=data_3_split1;
        elseif protein_number==2
            data_3=data_3_AP_split1;
        end
    else  % 与前后侧都不相邻，则跳过拆分步骤
    end
end


