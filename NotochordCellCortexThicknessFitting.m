%%% Notochord cell cortex thickness fitting

%%% Function:
%%% 1. Pre-assessment: Compare the number of cell data files in each folder and the number of cells in the length ratio file to determine if there are any missing cells.
%%% 2. Fitting characteristic parameters of the contractile ring:
%%%    (Fitting parameters: rho_inf, rho_inf_anterior/posterior represent amplitude; w, w2 are standard deviations; rho0 is the cortical basal thickness; 
%%%                         E is the mathematical expectation of the equatorial axis contractile ring, E_lateral is the mathematical expectation of the lateral contractile rings.)
%%%       Step 1: Fit the characteristic parameters of cells with no adjacent cells on either side and calculate the cell A-P axis planar cell polarity (PCP) intensity ratio.
%%%          (1) Search for all cell data files (csv files) in different folders under a given directory, distinguishing between cell data with and without adjacent cells.
%%%          (2) Identify cell data with no adjacent cells on either side, and remove excess data on both sides from the raw contractile ring grayscale data (except when lumen length/cell length < 0.2).
%%%          (3) Perform Tri-Gaussian distribution fitting on the processed cell data to obtain characteristic parameter values.
%%%          (4) Calculate the mean A-P axis (PCP) intensity ratio (A_P_ratio = rho_inf_anterior / rho_inf_posterior) using the fitted parameters.
%%%       Step 2: Fit the characteristic parameters of cells with adjacent cells on both sides.
%%%          (1) Determine the A-P polarity direction of each cell: 
%%%              If the current image contains data of cells without adjacent cells on either side, use those cells to determine the A-P axis direction.
%%%              If no data of cells without adjacent cells exist, test both A-P directions. If the processed cells (most) consistently exhibit the same A-P polarity, use that A-P polarity direction.
%%%          (2) Preprocess the three data files of each cell and store them in a matrix.
%%%          (3) Remove excess data on both sides from the raw contractile ring grayscale data (except when lumen length/cell length < 0.2).
%%%          (4) Perform Tri-Gaussian distribution fitting on the processed cell data to obtain characteristic parameter values.
%%% 3. Output of fitting results:
%%%       (1) Save the contractile ring characteristic parameters (9 parameters), error values, file names, and lumen length/cell length ratio in an xlsx file located in the root directory.
%%%       (2) Generate and save plots of each cell's measured values and fitted curves in subfolders for each date.

%%%       Classify and output data based on the length ratio measurement data corresponding to each csv file name, and save the fitted images to the same path as the csv file.
%%% Example of file name meaning:
%%% 0909_1058_1
%%% 09: September
%%% 09: 9th day
%%% 10: 10th image taken that day
%%% 5: The fifth cell in the image
%%% 8: Represents a specific protein
%%% 1: The first measurement data file for this cell in the left-to-right measurement sequence

warning('off','curvefit:fit:usingTrustRegion')
warning('off','curvefit:fit:noStartPoint')

clear all;close all;clc

% 定义文件名char型变量
global framepath2 filename_0 filename_xiahuaxian filename_fanxiegang filename_1 filename_2 filename_3 filename_xlsx filename_csv filename_month filename_day filename1 filename2 filename3
global data_3_split1 data_3_split2 data_3_AP_split1 data_3_AP_split2

framepath='C:\Users\Desktop\';
filename_0='0';
filename_xiahuaxian='_';
filename_fanxiegang='\';
filename_1='1';
filename_2='2';
filename_3='3';
filename_lifeact='lifeact';
filename_actin='actin';
filename_MLC='MLC';
filename_RhoA='RhoA';
filename_RBD='RBD';
filename_VD1='VD1';
filename_talinA='talinA';
filename_Cdc42WT='Cdc42WT';
filename_Cdc42D118A='Cdc42D118A';
filename_xlsx='.xlsx';
filename_csv='.csv';
filename_changdubi='长度比';
filename_Radius='半径 (Radius, 测量结果为lumen面积)';
filename_ContactAngle='接触角 (Contact Angle, 测量结果为圆周角)';

% 定义不同蛋白5个时期的文件名数组(用于写入到excel文件时使用)
filenamearray_all=cell(1000,5,9); % 所有文件名的分类汇总，其中第一维表示某一蛋白在某一时期的总数，第二维表示时期，第三维表示蛋白类型

% 定义不同蛋白5个时期的拟合数据数组(用于写入到excel文件时使用)
EquivalentTriGaussian_data_all=zeros(1000,10,5,9); % 等效三高斯拟合参数的分类汇总，其中第一维表示某一蛋白在某一时期的总数，第二维表示拟合参数，第三维表示时期，第四维表示蛋白类型
TriGaussian_data_all=zeros(1000,10,5,9);  % 直接三高斯拟合参数的分类汇总
MultiGaussian_data_all=zeros(1000,22,5,9);  % 多高斯拟合参数的分类汇总
data_changdubi_all=zeros(1000,5,9);  % 所有细胞lumen/cell长度比数据汇总
wucha_best_all=zeros(1000,5,9);
k_position1=ones(1,9);k_position2=ones(1,9);k_position3=ones(1,9);k_position4=ones(1,9);k_position5=ones(1,9);  % 用于统计每个时期参数输出到data_all矩阵中每个时期对应的位置

data_FittingParameterName=cell(1,10);
data_FittingParameterName{1,1}='基础厚度';
data_FittingParameterName{1,2}='赤道轴收缩环强度';
data_FittingParameterName{1,3}='赤道轴收缩环聚集程度';
data_FittingParameterName{1,4}='赤道轴收缩环数学期望';
data_FittingParameterName{1,5}='前侧收缩环强度';
data_FittingParameterName{1,6}='前侧收缩环聚集程度';
data_FittingParameterName{1,7}='前侧收缩环数学期望';
data_FittingParameterName{1,8}='后侧收缩环强度';
data_FittingParameterName{1,9}='后侧收缩环聚集程度';
data_FittingParameterName{1,10}='后侧收缩环数学期望';


%%% 拟合数据输出文件的路径创建和预写入
warning off MATLAB:xlswrite:AddSheet
for aa=1
    framepath_filename_lifeact=[framepath filename_lifeact filename_xlsx];
    framepath_filename_actin=[framepath filename_actin filename_xlsx];
    framepath_filename_MLC=[framepath filename_MLC filename_xlsx];
    framepath_filename_RhoA=[framepath filename_RhoA filename_xlsx];
    framepath_filename_RBD=[framepath filename_RBD filename_xlsx];
    framepath_filename_VD1=[framepath filename_VD1 filename_xlsx];
    framepath_filename_talinA=[framepath filename_talinA filename_xlsx];
    framepath_filename_Cdc42WT=[framepath filename_Cdc42WT filename_xlsx];
    framepath_filename_Cdc42D118A=[framepath filename_Cdc42D118A filename_xlsx];
    % lifeact xlsx文件预写入
    xlswrite(framepath_filename_lifeact,data_FittingParameterName,'sheet1','B1');
    % writematrix(M,'M.xls','Sheet',2,'Range','A3:E8')
    xlswrite(framepath_filename_lifeact,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_lifeact,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_lifeact,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_lifeact,data_FittingParameterName,'sheet5','B1');
    % actin xlsx文件预写入
    xlswrite(framepath_filename_actin,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_actin,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_actin,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_actin,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_actin,data_FittingParameterName,'sheet5','B1');
    % MLC xlsx文件预写入
    xlswrite(framepath_filename_MLC,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_MLC,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_MLC,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_MLC,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_MLC,data_FittingParameterName,'sheet5','B1');
    % RhoA xlsx文件预写入
    xlswrite(framepath_filename_RhoA,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_RhoA,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_RhoA,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_RhoA,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_RhoA,data_FittingParameterName,'sheet5','B1');
    % RBD xlsx文件预写入
    xlswrite(framepath_filename_RBD,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_RBD,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_RBD,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_RBD,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_RBD,data_FittingParameterName,'sheet5','B1');
    % VD1 xlsx文件预写入
    xlswrite(framepath_filename_VD1,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_VD1,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_VD1,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_VD1,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_VD1,data_FittingParameterName,'sheet5','B1');
    % talinA xlsx文件预写入
    xlswrite(framepath_filename_talinA,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_talinA,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_talinA,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_talinA,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_talinA,data_FittingParameterName,'sheet5','B1');
    % Cdc42WT xlsx文件预写入
    xlswrite(framepath_filename_Cdc42WT,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_Cdc42WT,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_Cdc42WT,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_Cdc42WT,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_Cdc42WT,data_FittingParameterName,'sheet5','B1');
    % Cdc42D118A xlsx文件预写入
    xlswrite(framepath_filename_Cdc42D118A,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_Cdc42D118A,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_Cdc42D118A,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_Cdc42D118A,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_Cdc42D118A,data_FittingParameterName,'sheet5','B1');
end

%%% 拟合数据前首先判断文件数和细胞数是否一致
for k=1:9
    filename3=num2str(k);
    
    for month=1:12
        if month<10
            filename_month=[filename_0 num2str(month)];
        else
            filename_month=num2str(month);
        end
        for day=1:31
            if day<10
                filename_day=[filename_0 num2str(day)];
            else
                filename_day=num2str(day);
            end
            number_changdubi=0;  % 用于输出文件时指示长度比文件中的数据位置
            
            filename_changdubi='长度比';
            framepath2=[framepath num2str(k) filename_fanxiegang filename_month filename_day filename_fanxiegang];  % 子文件夹路径，如 'C:\Users\Desktop\1\0928'
            framepath_filename_changdubi=[framepath2 filename_changdubi filename_xlsx];
            
            if exist(framepath_filename_changdubi,'file')   % 长度比文件是否存在可用于快速筛选是否存在该天拍摄该蛋白的数据，缩短查找时间
                file_number=0;
                data_lumenlength=xlsread(framepath_filename_changdubi,'sheet1','A1:A100');  % 读入该文件下的lumen长度
                data_celllength=xlsread(framepath_filename_changdubi,'sheet1','B1:B100');  % 读入该文件下的cell长度
                data_changdubi=data_lumenlength./data_celllength;
                cell_number=length(data_changdubi);
                
                for i=1:30  % 外侧循环用于查找目标文件是否存在(i表示第几张photo，j表示photo中的第几个细胞，k表示荧光通道)
                    filename1=num2str(i);
                    for j=1:10
                        filename2=num2str(j);
                        filename=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3];  % 如 0928_111
                        framepath_filename=[framepath2 filename];
                        
                        % 拼接完整的文件路径，用于下一步exist函数的查找
                        framepath_filename_1=[framepath_filename filename_xiahuaxian filename_1 filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111_1.csv
                        framepath_filename_2=[framepath_filename filename_xiahuaxian filename_2 filename_csv];
                        framepath_filename_3=[framepath_filename filename_xiahuaxian filename_3 filename_csv];
                        framepath_filename_danyi=[framepath_filename filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111.csv
                        
                        if exist(framepath_filename_danyi,'file') || exist(framepath_filename_1,'file')
                            file_number=file_number+1;
                        end
                    end
                end
                if cell_number~=file_number
                    fprintf('缺少数据文件');
                    framepath2
                    abc
                end
            end
        end
    end
end


%% Step 1
for k=1:9  % 1表示lifeact，2表示actin，3表示MLC，4表示RhoA，5表示Anillin RBD，6表示VD1,7表示talinA，8表示Cdc42WT，9表示Cdc42D118A
    filename3=num2str(k);
    
    if k==1
        filename_protein='lifeact';
    elseif k==2
        filename_protein='actin';
    elseif k==3
        filename_protein='MLC';
    elseif k==4
        filename_protein='RhoA';
    elseif k==5
        filename_protein='Anillin RBD';
    elseif k==6
        filename_protein='VD1';
    elseif k==7
        filename_protein='talinA';
    elseif k==8
        filename_protein='Cdc42WT';
    elseif k==9
        filename_protein='Cdc42D118A';
    end
    filename_photo='photo';
    filename_cell='cell';
    filename_kongge=' ';
    
    for month=1:12
        if month<10
            filename_month=[filename_0 num2str(month)];
        else
            filename_month=num2str(month);
        end
        for day=1:31
            if day<10
                filename_day=[filename_0 num2str(day)];
            else
                filename_day=num2str(day);
            end
            number_changdubi=0;  % 用于输出文件时指示长度比文件中的数据位置
            
            framepath2=[framepath num2str(k) filename_fanxiegang filename_month filename_day filename_fanxiegang];  % 子文件夹路径，如 'C:\Users\Desktop\1\0928'
            framepath_filename_changdubi=[framepath2 filename_changdubi filename_xlsx];
            framepath_filename_Radius=[framepath2 filename_Radius filename_xlsx];
            framepath_filename_ContactAngle=[framepath2 filename_ContactAngle filename_xlsx];
            
            if exist(framepath_filename_changdubi,'file')   % 长度比文件是否存在可用于快速筛选是否存在该天拍摄该蛋白的数据，缩短查找时间
                data_lumenradius=xlsread(framepath_filename_Radius,'sheet1','A1:A100');  % 读入lumen半径
                data_lumenradius=sqrt(data_lumenradius./pi);
                data_lumen_ContactAngle=xlsread(framepath_filename_ContactAngle,'sheet1','A1:A100');  % 读入lumen接触角
                data_lumenlength=xlsread(framepath_filename_changdubi,'sheet1','A1:A100');  % 读入该文件下的lumen长度
                data_celllength=xlsread(framepath_filename_changdubi,'sheet1','B1:B100');  % 读入该文件下的cell长度
                % 若Lumen length大于cell length，则cell length=Lumen length+0.5
                for i=1:length(data_lumenradius)
                    if 2*data_lumenradius(i)*sin(data_lumen_ContactAngle(i)*pi./180)>data_celllength(i)
                        data_celllength(i)=2*data_lumenradius(i)*sin(data_lumen_ContactAngle(i)*pi./180)+0.5;
                    end
                end
                data_changdubi=2.*data_lumenradius.*sin(data_lumen_ContactAngle.*pi./180)./data_celllength;

                for i=1:30  % 外侧循环用于查找目标文件是否存在(i表示第几张photo，j表示photo中的第几个细胞，k表示荧光通道)
                    filename1=num2str(i);
                    for j=1:10
                        filename2=num2str(j);
                        filename=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3];  % 如 0928_111
                        framepath_filename=[framepath2 filename];
                        
                        % 拼接完整的文件路径，用于下一步exist函数的查找
                        framepath_filename_1=[framepath_filename filename_xiahuaxian filename_1 filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111_1.csv
                        framepath_filename_2=[framepath_filename filename_xiahuaxian filename_2 filename_csv];
                        framepath_filename_3=[framepath_filename filename_xiahuaxian filename_3 filename_csv];
                        framepath_filename_danyi=[framepath_filename filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111.csv
                        
                        if exist(framepath_filename_danyi,'file')   % 一个细胞的数据有单个文件
                            number_changdubi=number_changdubi+1;
                            for aa=1
                                %%% 对原始细胞收缩环灰度值数据进行两侧多余数据的删减处理
                                % 根据测量数据拟合basal-lateral domain双高斯分布
                                % 数据点拉伸
                                data=xlsread(framepath_filename_danyi); % 导入数据
                                [X,Y]=size(data); % 丈量数据矩阵尺寸
                                if mod(X,2)==0   % 若数据为偶数个，则进行线性插值，方便后续的平移工作中的对称性
                                    data2=zeros(X+1,Y);
                                    data2(:,1)=imresize(data(:,1),[X+1,1], 'bilinear'); % 将X轴进行线性插值
                                    data2(:,2)=imresize(data(:,2),[X+1,1], 'bilinear'); % 将Y轴进行线性插值
                                else
                                    data2=data;
                                end
                                
                                % 计算细胞的anterior-lateral basal posterior-lateral domain的长度
                                basal_lateral_length=data(end,1);  % 细胞的basal-lateral domain总长度
                                anterior_lateral_length=data_celllength(number_changdubi)-data_lumenlength(number_changdubi);
                                posterior_lateral_length=anterior_lateral_length;
                                basal_length=basal_lateral_length-anterior_lateral_length-posterior_lateral_length;
                                
                                % X轴平移至关于原点对称的区间，并归一化
                                bias=0.5; % 平移量
                                x=data2(:,1)/max(data2(:,1)); % 对x轴进行归一化
                                x=x-bias; % x轴平移至关于坐标原点对称
                                % Y轴数据归一化
                                y=data2(:,2)/max(data2(:,2));
                                
                                
                                Y_anterior=zeros(1,round((X-1)*0.2));  % 定义细胞两侧的收缩环灰度值数据，用于单高斯分布拟合
                                Y_posterior=zeros(1,round((X-1)*0.2));
                                Y_anterior_AP=zeros(1,round((X-1)*0.2));  % 定义细胞另一蛋白两侧的收缩环灰度值数据，用于单高斯分布拟合
                                Y_posterior_AP=zeros(1,round((X-1)*0.2));

                                if data_changdubi(number_changdubi)<=0.2  % lumen形成初期两侧收缩环极性未建立，不删除数据
                                    for k_AP=1:10
                                        if k_AP~=k
                                            filename3_AP=num2str(k_AP);  % AP: anothor protein
                                            filename_AP=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3_AP];  % 如 0928_111
                                            framepath_filename_AP=[framepath2 filename_AP];
                                            framepath_filename_danyi_AP=[framepath_filename_AP filename_csv];
                                            if exist(framepath_filename_danyi_AP,'file')
                                                break
                                            end
                                        end
                                    end
                                    if k_AP~=10  % 表示该细胞有另一个蛋白数据
                                        % 读取该蛋白的数据
                                        % 数据点拉伸
                                        data_AP=xlsread(framepath_filename_danyi_AP); % 导入数据
                                        [X_AP,Y_AP]=size(data_AP); % 丈量数据矩阵尺寸
                                        if mod(X_AP,2)==0   % 若数据为偶数个，则进行线性插值，方便后续的平移工作中的对称性
                                            data2_AP=zeros(X_AP+1,Y_AP);
                                            data2_AP(:,1)=imresize(data_AP(:,1),[X_AP+1,1], 'bilinear'); % 将X轴进行线性插值
                                            data2_AP(:,2)=imresize(data_AP(:,2),[X_AP+1,1], 'bilinear'); % 将Y轴进行线性插值
                                        else
                                            data2_AP=data_AP;
                                        end

                                        % X轴平移至关于原点对称的区间，并归一化
                                        bias=0.5; % 平移量
                                        x_AP=data2_AP(:,1)/max(data2_AP(:,1)); % 对x轴进行归一化
                                        x_AP=x_AP-bias; % x轴平移至关于坐标原点对称
                                        % Y轴数据归一化
                                        y_AP=data2_AP(:,2)/max(data2_AP(:,2));
                                    end
                                    x_min=-0.5;x_max=0.5;
                                else  % lumen形成中后期两侧收缩环极性已建立，需删除数据
                                    for k_AP=1:10
                                        if k_AP~=k
                                            filename3_AP=num2str(k_AP);  % AP: anothor protein
                                            filename_AP=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3_AP];  % 如 0928_111
                                            framepath_filename_AP=[framepath2 filename_AP];
                                            framepath_filename_danyi_AP=[framepath_filename_AP filename_csv];
                                            if exist(framepath_filename_danyi_AP,'file')
                                                break
                                            end
                                        end
                                    end
                                    if k_AP~=10  % 表示该细胞有另一个蛋白数据
                                        % 读取该蛋白的数据
                                        % 数据点拉伸
                                        data_AP=xlsread(framepath_filename_danyi_AP); % 导入数据
                                        [X_AP,Y_AP]=size(data_AP); % 丈量数据矩阵尺寸
                                        if mod(X_AP,2)==0   % 若数据为偶数个，则进行线性插值，方便后续的平移工作中的对称性
                                            data2_AP=zeros(X_AP+1,Y_AP);
                                            data2_AP(:,1)=imresize(data_AP(:,1),[X_AP+1,1], 'bilinear'); % 将X轴进行线性插值
                                            data2_AP(:,2)=imresize(data_AP(:,2),[X_AP+1,1], 'bilinear'); % 将Y轴进行线性插值
                                        else
                                            data2_AP=data_AP;
                                        end

                                        % X轴平移至关于原点对称的区间，并归一化
                                        bias=0.5; % 平移量
                                        x_AP=data2_AP(:,1)/max(data2_AP(:,1)); % 对x轴进行归一化
                                        x_AP=x_AP-bias; % x轴平移至关于坐标原点对称
                                        % Y轴数据归一化
                                        y_AP=data2_AP(:,2)/max(data2_AP(:,2));

                                        % Levenberg-Marquardt method拟合两侧高斯分布，确定要删去的部分
                                        for ab=1
                                            if data_changdubi(number_changdubi)<=0.2  % lumen形成早期，lateral domain收缩环未建立，共有9个特征参数拟合
                                                E_lateral_start=0.3;
                                                E_lateral_end=0.5;
                                            else  % lumen形成中后期，lateral domain收缩环建立，共有8个特征参数拟合
                                                E_lateral_start=0.5;
                                                E_lateral_end=0.5;
                                            end

                                            wucha_best=100000000;
                                            % 使用Levenberg-Marquardt method拟合
                                            for i_fit=1:5  % 循环10次，寻找最优拟合解
                                                ft=fittype('rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)','independent','x','dependent','y');
                                                opts=fitoptions('Method','NonlinearLeastSquares');
                                                opts.Algorithm='Levenberg-Marquardt';
                                                opts.Display='Off';
                                                opts.MaxFunEvals=600;
                                                opts.MaxIter=400;
                                                opts.Lower=[-0.2 -0.5 0.45 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
                                                opts.Upper=[0.2 -0.45 0.5 1 5 5 5 1 1 1];
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
                                                E_anterior=fitresult.E_anterior;
                                                E_posterior=fitresult.E_posterior;

                                                wucha=0;
                                                for ii=1:length(x)
                                                    wucha=wucha+(rho0+rho_inf*exp(-1/2*((x(ii)-E)/w)^2)+rho_inf_anterior*exp(-1/2*((x(ii)-E_anterior)/w_anterior)^2)+...
                                                        rho_inf_posterior*exp(-1/2*((x(ii)-E_posterior)/w_posterior)^2)-y(ii))^2;
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
                                                    E_anterior_best=E_anterior;
                                                    E_posterior_best=E_posterior;
                                                end
                                            end
                                        end
                                        % 另一个蛋白数据拟合
                                        for ab=1
                                            if data_changdubi(number_changdubi)<=0.2  % lumen形成早期，lateral domain收缩环未建立，共有9个特征参数拟合
                                                E_lateral_start=0.3;
                                                E_lateral_end=0.5;
                                            else  % lumen形成中后期，lateral domain收缩环建立，共有8个特征参数拟合
                                                E_lateral_start=0.5;
                                                E_lateral_end=0.5;
                                            end

                                            wucha_best=100000000;
                                            % 使用Levenberg-Marquardt method拟合
                                            for i_fit=1:5  % 循环10次，寻找最优拟合解
                                                ft=fittype('rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)','independent','x','dependent','y');
                                                opts=fitoptions('Method','NonlinearLeastSquares');
                                                opts.Algorithm='Levenberg-Marquardt';
                                                opts.Display='Off';
                                                opts.MaxFunEvals=600;
                                                opts.MaxIter=400;
                                                opts.Lower=[-0.2 -0.5 0.45 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
                                                opts.Upper=[0.2 -0.45 0.5 1 5 5 5 1 1 1];
                                                [xData,yData]=prepareCurveData(x_AP,y_AP);
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
                                                E_anterior=fitresult.E_anterior;
                                                E_posterior=fitresult.E_posterior;

                                                wucha=0;
                                                for ii=1:length(x_AP)
                                                    wucha=wucha+(rho0+rho_inf*exp(-1/2*((x_AP(ii)-E)/w)^2)+rho_inf_anterior*exp(-1/2*((x_AP(ii)-E_anterior)/w_anterior)^2)+...
                                                        rho_inf_posterior*exp(-1/2*((x_AP(ii)-E_posterior)/w_posterior)^2)-y_AP(ii))^2;
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
                                                    E_anterior_AP_best=E_anterior;
                                                    E_posterior_AP_best=E_posterior;
                                                end
                                            end
                                        end
                                        x_min=min(E_anterior_best,E_anterior_AP_best);  % 前侧边界取两个蛋白的拟合前边界中较小的
                                        x_max=max(E_posterior_best,E_posterior_AP_best);  % 后侧边界取两个蛋白的拟合后边界中较大的
                                    else  % kk=9表示该细胞仅有一种蛋白的数据文件
                                        % 读取该蛋白的数据
                                        % 数据点拉伸

                                        % Levenberg-Marquardt method拟合两侧高斯分布，确定要删去的部分
                                        for ab=1
                                            if data_changdubi(number_changdubi)<=0.2  % lumen形成早期，lateral domain收缩环未建立，共有9个特征参数拟合
                                                E_lateral_start=0.3;
                                                E_lateral_end=0.5;
                                            else  % lumen形成中后期，lateral domain收缩环建立，共有8个特征参数拟合
                                                E_lateral_start=0.5;
                                                E_lateral_end=0.5;
                                            end

                                            wucha_best=100000000;
                                            % 使用Levenberg-Marquardt method拟合
                                            for i_fit=1:5  % 循环10次，寻找最优拟合解
                                                ft=fittype('rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)','independent','x','dependent','y');
                                                opts=fitoptions('Method','NonlinearLeastSquares');
                                                opts.Algorithm='Levenberg-Marquardt';
                                                opts.Display='Off';
                                                opts.MaxFunEvals=600;
                                                opts.MaxIter=400;
                                                opts.Lower=[-0.2 -0.5 0.45 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
                                                opts.Upper=[0.2 -0.45 0.5 1 5 5 5 1 1 1];
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
                                                E_anterior=fitresult.E_anterior;
                                                E_posterior=fitresult.E_posterior;

                                                wucha=0;
                                                for ii=1:length(x)
                                                    wucha=wucha+(rho0+rho_inf*exp(-1/2*((x(ii)-E)/w)^2)+rho_inf_anterior*exp(-1/2*((x(ii)-E_anterior)/w_anterior)^2)+...
                                                        rho_inf_posterior*exp(-1/2*((x(ii)-E_posterior)/w_posterior)^2)-y(ii))^2;
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
                                                    E_anterior_best=E_anterior;
                                                    E_posterior_best=E_posterior;
                                                end
                                            end
                                        end
                                        x_min=E_anterior_best;
                                        x_max=E_posterior_best;
                                    end
                                end
                                % 若x_min或x_max偏离边界过多(超过总细胞长度的1/10)，说明数据拟合效果不好，仍取原始边界值
                                if (x_min+0.5)/0.5>0.1
                                    x_min=-0.5;
                                elseif (0.5-x_max)/0.5>0.1
                                    x_max=0.5;
                                end
                                
                                %%% 对删减处理完成的细胞灰度值数据进行三高斯分布拟合
                                % 找到原X轴上和x_min x_max最接近的点
                                wucha_best=100;
                                for ii=1:X
                                    wucha=abs(x(ii)-x_min);
                                    if wucha<wucha_best
                                        wucha_best=wucha;
                                        x_min_best=x(ii);  % 将最接近x_min的原X坐标值暂时保存在x_min_best变量中
                                    end
                                end
                                x_min=x_min_best;  % x_min改成X轴上最接近x_min点的值
                                wucha_best=100;
                                for ii=1:X
                                    wucha=abs(x(ii)-x_max);
                                    if wucha<wucha_best
                                        wucha_best=wucha;
                                        x_max_best=x(ii);  % 将最接近x_max的原X坐标值暂时保存在x_max_best变量中
                                    end
                                end
                                x_max=x_max_best;  % x_max改成X轴上最接近x_max点的值
                                
                                % 删去X轴和对应灰度值(y)上小于x_min和大于x_max的部分
                                x_delete=zeros(round((X-1)*(x_max-x_min)+1),1);
                                y_delete=zeros(round((X-1)*(x_max-x_min)+1),1);
                                for ii=1:round((X-1)*(x_max-x_min)+1)
                                    x_delete(ii)=x(ii+round((x_min+0.5)*(X-1)));
                                    y_delete(ii)=y(ii+round((x_min+0.5)*(X-1)));
                                end
                                x=x_delete;
                                y=y_delete;
                                % X轴重新平移至关于原点对称的区间，并拉伸至(-0.5,0.5)
                                bias=x_max+x_min; % 平移量
                                stretch=1/(x_max-x_min); % 拉伸量
                                x=x-bias/2; % x轴平移至关于坐标原点对称
                                x=x*stretch;
                                
                                %%% 使用Levenberg-Marquardt method拟合
                                max_y=0;min_y=100;
                                for iy=floor(length(y)*0.25):floor(length(y)*0.75)
                                    if y(iy)>max_y
                                        max_y=y(iy);
                                    elseif y(iy)<min_y
                                        min_y=y(iy);
                                    end
                                end
                                if max_y/min_y<1.5  % 表示赤道轴收缩环已经几乎消失，lumen时期处于后期，为防止拟合出现误差，将收缩环位置固定在中央
                                    E_bottom=-0.05;E_top=0.05;
                                else  % 收缩环较强，其位置不固定
                                    E_bottom=-0.1;E_top=0.1;
                                end

                                % 多高斯拟合
                                % 调用MultiGaussian_fit子函数计算Multi-Gaussian distribution
                                cell_type=1;
                                [rho0_best,rho_inf_array,E_array,w_array,break_point,wucha_best,TriGaussian_parameter]=MultiGaussian_fit(x,y,data_changdubi(number_changdubi),anterior_lateral_length,basal_length,posterior_lateral_length,cell_type);
                                framepath_filename
                                % 调用MultiGaussianTransfertoEquivalentTriGaussian子函数计算等效Tri-Gaussian distribution
                                [EquivalentTriGaussian_parameter,area_sumsum]=MultiGaussianTransfertoEquivalentTriGaussian(x,y,rho0_best,rho_inf_array,E_array,w_array,break_point);
                            end
                            
                            % 判断细胞是否是左侧为anterior part，右侧为posterior part。若不是，将拟合参数、拟合曲线、原始数据关于原点做对称处理
                            if rho_inf_array(end-1)<rho_inf_array(end)
                                % MultiGaussian参数做对称处理
                                rho_inf_array_reverse=zeros(size(rho_inf_array));
                                E_array_reverse=zeros(size(rho_inf_array));
                                w_array_reverse=zeros(size(rho_inf_array));
                                for i_Gaussian=1:length(rho_inf_array)-2
                                    rho_inf_array_reverse(i_Gaussian)=rho_inf_array(end-i_Gaussian-1);
                                    E_array_reverse(i_Gaussian)=E_array(end-i_Gaussian-1);
                                    w_array_reverse(i_Gaussian)=w_array(end-i_Gaussian-1);
                                end
                                rho_inf_array_reverse(end-1)=rho_inf_array(end);rho_inf_array_reverse(end)=rho_inf_array(end-1);
                                E_array_reverse(end-1)=E_array(end);E_array_reverse(end)=E_array(end-1);
                                w_array_reverse(end-1)=w_array(end);w_array_reverse(end)=w_array(end-1);
                                
                                rho_inf_array=rho_inf_array_reverse;
                                E_array=-E_array_reverse;
                                w_array=w_array_reverse;

                                % EquivalentTriGaussian参数做对称处理
                                anterior_posterior_ring_reverse=EquivalentTriGaussian_parameter(5:6);
                                EquivalentTriGaussian_parameter(5:6)=EquivalentTriGaussian_parameter(8:9);
                                EquivalentTriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                EquivalentTriGaussian_parameter(4)=-EquivalentTriGaussian_parameter(4);

                                % TriGaussian参数做对称处理
                                anterior_posterior_ring_reverse=TriGaussian_parameter(5:6);
                                TriGaussian_parameter(5:6)=TriGaussian_parameter(8:9);
                                TriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                TriGaussian_parameter(4)=-TriGaussian_parameter(4);
                                
                                y_reverse=zeros(length(y),1);  % 原始数据翻转
                                for ii=1:length(y)
                                    y_reverse(length(y)-ii+1)=y(ii);
                                end
                                y=y_reverse;
                            end

                            % 荧光强度值(纵坐标)归一化
                            % 将拟合完成的多高斯分布在定义域上的数学期望变换为1，同时将原始测量数据和等效三高斯分布按等比例拉伸
                            if length(rho_inf_array)-2==1
                                f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                    +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2);
                            elseif length(rho_inf_array)-2==2
                                f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                    +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2);
                            elseif length(rho_inf_array)-2==3
                                f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                    +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2)...
                                    +rho_inf_array(5)*exp(-1/2*((y-E_array(5))/w_array(5)).^2);
                            elseif length(rho_inf_array)-2==4
                                f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                    +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2)+...
                                    rho_inf_array(5)*exp(-1/2*((y-E_array(5))/w_array(5)).^2)+rho_inf_array(6)*exp(-1/2*((y-E_array(6))/w_array(6)).^2);
                            elseif length(rho_inf_array)-2>=5
                                f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                    +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2)...
                                    +rho_inf_array(5)*exp(-1/2*((y-E_array(5))/w_array(5)).^2)+rho_inf_array(6)*exp(-1/2*((y-E_array(6))/w_array(6)).^2)...
                                    +rho_inf_array(7)*exp(-1/2*((y-E_array(7))/w_array(7)).^2);
                            end

                            a0=integral(f,-1/2,1/2);
                            rho0_best=rho0_best/a0;
                            rho_inf_array=rho_inf_array./a0;
                            EquivalentTriGaussian_parameter(2:3:8)=EquivalentTriGaussian_parameter(2:3:8)./a0;EquivalentTriGaussian_parameter(1)=EquivalentTriGaussian_parameter(1)./a0;
                            TriGaussian_parameter(2:3:8)=TriGaussian_parameter(2:3:8)./a0;TriGaussian_parameter(1)=TriGaussian_parameter(1)./a0;
                            y_MultiGaussian_fit=ones(size(x))*rho0_best;  % 多高斯拟合曲线
                            for i_Gaussian=1:length(rho_inf_array)
                                y_MultiGaussian_fit=y_MultiGaussian_fit+rho_inf_array(i_Gaussian)*exp(-1/2*((x-E_array(i_Gaussian))/w_array(i_Gaussian)).^2);
                            end
                            y_EquivalentTriGaussian_fit=ones(size(x))*EquivalentTriGaussian_parameter(1);  % 等效三高斯拟合曲线
                            for i_Gaussian=1:3
                                y_EquivalentTriGaussian_fit=y_EquivalentTriGaussian_fit+EquivalentTriGaussian_parameter(3*i_Gaussian-1)*exp(-1/2*((x-EquivalentTriGaussian_parameter(3*i_Gaussian+1))/EquivalentTriGaussian_parameter(3*i_Gaussian)).^2);
                            end
                            y=y./a0;
                            

                            %%% 绘图
                            filename_plot=[filename_photo filename_kongge filename1 filename_kongge filename_cell filename_kongge filename2 filename_kongge filename_protein];
                            
                            figure
                            plot(x,y,'.','MarkerSize',15,'Color','k');  % 绘制数据点
                            hold on;
                            plot(x,y_MultiGaussian_fit,'b','linewidth',3);  % 绘制拟合的高斯曲线
                            plot(x,y_EquivalentTriGaussian_fit,'r','linewidth',3);  % 绘制拟合的高斯曲线
                            hold off;
                            axis ([-0.5 0.5 0 max(y)+0.1]);
                            title(filename_plot)
                            xlabel('x')
                            ylabel('cortex thickness')
                            set(gca,'YTick',0:0.2:max(y)+0.1);  % 修改y坐标刻度
                            set(gca,'XTick',-0.5:0.1:0.5);  % 修改y坐标刻度
                            legend('Experiments','Multi-Gaussian Curve','Tri-Gaussian Curve')
                            set(gca,'FontName','Arial','FontSize',12)
                            saveas(gcf,[framepath2,filename_plot,'.jpg']);
                            close

                            %%% 拟合参数汇总
                            for aa=1
                                if data_changdubi(number_changdubi)<=0.2
                                    data_changdubi_all(k_position1(k),1,k)=data_changdubi(number_changdubi);
                                    filenamearray_all{k_position1(k),1,k}=filename;
                                    MultiGaussian_data_all(k_position1(k),1,1,k)=rho0_best;
                                    MultiGaussian_data_all(k_position1(k),2:length(rho_inf_array)+1,1,k)=rho_inf_array;
                                    MultiGaussian_data_all(k_position1(k),length(rho_inf_array)+2:2*length(rho_inf_array)+1,1,k)=w_array;
                                    MultiGaussian_data_all(k_position1(k),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,1,k)=E_array;
                                    EquivalentTriGaussian_data_all(k_position1(k),:,1,k)=EquivalentTriGaussian_parameter;
                                    TriGaussian_data_all(k_position1(k),:,1,k)=TriGaussian_parameter;
                                    wucha_best_all(k_position1(k),1,k)=wucha_best;
                                    k_position1(k)=k_position1(k)+1;
                                elseif data_changdubi(number_changdubi)>0.2 && data_changdubi(number_changdubi)<=0.4
                                    data_changdubi_all(k_position2(k),2,k)=data_changdubi(number_changdubi);
                                    filenamearray_all{k_position2(k),2,k}=filename;
                                    MultiGaussian_data_all(k_position2(k),1,2,k)=rho0_best;
                                    MultiGaussian_data_all(k_position2(k),2:length(rho_inf_array)+1,2,k)=rho_inf_array;
                                    MultiGaussian_data_all(k_position2(k),length(rho_inf_array)+2:2*length(rho_inf_array)+1,2,k)=w_array;
                                    MultiGaussian_data_all(k_position2(k),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,2,k)=E_array;
                                    EquivalentTriGaussian_data_all(k_position2(k),:,2,k)=EquivalentTriGaussian_parameter;
                                    TriGaussian_data_all(k_position2(k),:,2,k)=TriGaussian_parameter;
                                    wucha_best_all(k_position2(k),2,k)=wucha_best;
                                    k_position2(k)=k_position2(k)+1;
                                elseif data_changdubi(number_changdubi)>0.4 && data_changdubi(number_changdubi)<=0.6
                                    data_changdubi_all(k_position3(k),3,k)=data_changdubi(number_changdubi);
                                    filenamearray_all{k_position3(k),3,k}=filename;
                                    MultiGaussian_data_all(k_position3(k),1,3,k)=rho0_best;
                                    MultiGaussian_data_all(k_position3(k),2:length(rho_inf_array)+1,3,k)=rho_inf_array;
                                    MultiGaussian_data_all(k_position3(k),length(rho_inf_array)+2:2*length(rho_inf_array)+1,3,k)=w_array;
                                    MultiGaussian_data_all(k_position3(k),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,3,k)=E_array;
                                    EquivalentTriGaussian_data_all(k_position3(k),:,3,k)=EquivalentTriGaussian_parameter;
                                    TriGaussian_data_all(k_position3(k),:,3,k)=TriGaussian_parameter;
                                    wucha_best_all(k_position3(k),3,k)=wucha_best;
                                    k_position3(k)=k_position3(k)+1;
                                elseif data_changdubi(number_changdubi)>0.6 && data_changdubi(number_changdubi)<=0.8
                                    data_changdubi_all(k_position4(k),4,k)=data_changdubi(number_changdubi);
                                    filenamearray_all{k_position4(k),4,k}=filename;
                                    MultiGaussian_data_all(k_position4(k),1,4,k)=rho0_best;
                                    MultiGaussian_data_all(k_position4(k),2:length(rho_inf_array)+1,4,k)=rho_inf_array;
                                    MultiGaussian_data_all(k_position4(k),length(rho_inf_array)+2:2*length(rho_inf_array)+1,4,k)=w_array;
                                    MultiGaussian_data_all(k_position4(k),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,4,k)=E_array;
                                    EquivalentTriGaussian_data_all(k_position4(k),:,4,k)=EquivalentTriGaussian_parameter;
                                    TriGaussian_data_all(k_position4(k),:,4,k)=TriGaussian_parameter;
                                    wucha_best_all(k_position4(k),4,k)=wucha_best;
                                    k_position4(k)=k_position4(k)+1;
                                else
                                    data_changdubi_all(k_position5(k),5,k)=data_changdubi(number_changdubi);
                                    filenamearray_all{k_position5(k),5,k}=filename;
                                    MultiGaussian_data_all(k_position5(k),1,5,k)=rho0_best;
                                    MultiGaussian_data_all(k_position5(k),2:length(rho_inf_array)+1,5,k)=rho_inf_array;
                                    MultiGaussian_data_all(k_position5(k),length(rho_inf_array)+2:2*length(rho_inf_array)+1,5,k)=w_array;
                                    MultiGaussian_data_all(k_position5(k),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,5,k)=E_array;
                                    EquivalentTriGaussian_data_all(k_position5(k),:,5,k)=EquivalentTriGaussian_parameter;
                                    TriGaussian_data_all(k_position5(k),:,5,k)=TriGaussian_parameter;
                                    wucha_best_all(k_position5(k),5,k)=wucha_best;
                                    k_position5(k)=k_position5(k)+1;
                                end

                            end
                            
                            %%% 拟合该细胞的另一个蛋白数据文件
                            if k_AP~=10  % 该细胞存在另一个蛋白数据文件
                                if k_AP==1
                                    filename_protein='lifeact';
                                elseif k_AP==2
                                    filename_protein='actin';
                                elseif k_AP==3
                                    filename_protein='MLC';
                                elseif k_AP==4
                                    filename_protein='RhoA';
                                elseif k_AP==5
                                    filename_protein='Anillin RBD';
                                elseif k_AP==6
                                    filename_protein='VD1';
                                elseif k_AP==7
                                    filename_protein='talinA';
                                elseif k_AP==8
                                    filename_protein='Cdc42WT';
                                elseif k_AP==9
                                    filename_protein='Cdc42D118A';
                                end
                                
                                %%% 对删减处理完成的细胞灰度值数据进行三高斯分布拟合
                                % 删去X轴和对应灰度值(y)上小于x_min和大于x_max的部分
                                x_delete_AP=zeros(round((X-1)*(x_max-x_min)+1),1);
                                y_delete_AP=zeros(round((X-1)*(x_max-x_min)+1),1);
                                for ii=1:round((X-1)*(x_max-x_min)+1)
                                    x_delete_AP(ii)=x_AP(ii+round((x_min+0.5)*(X-1)));
                                    y_delete_AP(ii)=y_AP(ii+round((x_min+0.5)*(X-1)));
                                end
                                x_AP=x_delete_AP;
                                y_AP=y_delete_AP;
                                % X轴重新平移至关于原点对称的区间，并拉伸至(-0.5,0.5)
                                bias=x_max+x_min; % 平移量
                                stretch=1/(x_max-x_min); % 拉伸量
                                x_AP=x_AP-bias/2; % x轴平移至关于坐标原点对称
                                x_AP=x_AP*stretch;

                                %%% 使用Levenberg-Marquardt method拟合
                                max_y=0;min_y=100;
                                for iy=floor(length(y_AP)*0.25):floor(length(y_AP)*0.75)
                                    if y_AP(iy)>max_y
                                        max_y=y_AP(iy);
                                    elseif y_AP(iy)<min_y
                                        min_y=y_AP(iy);
                                    end
                                end
                                if max_y/min_y<1.5  % 表示赤道轴收缩环已经几乎消失，lumen时期处于后期，为防止拟合出现误差，将收缩环位置固定在中央
                                    E_bottom=-0.05;E_top=0.05;
                                else  % 收缩环较强，其位置不固定
                                    E_bottom=-0.1;E_top=0.1;
                                end

                                % 多高斯拟合
                                % 调用MultiGaussian_fit子函数计算Multi-Gaussian distribution
                                cell_type=1;
                                [rho0_best,rho_inf_array,E_array,w_array,break_point,wucha_best,TriGaussian_parameter]=MultiGaussian_fit(x_AP,y_AP,data_changdubi(number_changdubi),anterior_lateral_length,basal_length,posterior_lateral_length,cell_type);
                                framepath_filename_AP
                                % 调用MultiGaussianTransfertoEquivalentTriGaussian子函数计算等效Tri-Gaussian distribution
                                [EquivalentTriGaussian_parameter,area_sumsum]=MultiGaussianTransfertoEquivalentTriGaussian(x_AP,y_AP,rho0_best,rho_inf_array,E_array,w_array,break_point);

                                % 判断细胞是否是左侧为anterior part，右侧为posterior part。若不是，将拟合参数、拟合曲线、原始数据关于原点做对称处理
                                if rho_inf_array(end-1)<rho_inf_array(end)
                                    % MultiGaussian参数做对称处理
                                    rho_inf_array_reverse=zeros(size(rho_inf_array));
                                    E_array_reverse=zeros(size(rho_inf_array));
                                    w_array_reverse=zeros(size(rho_inf_array));
                                    for i_Gaussian=1:length(rho_inf_array)-2
                                        rho_inf_array_reverse(i_Gaussian)=rho_inf_array(end-i_Gaussian-1);
                                        E_array_reverse(i_Gaussian)=E_array(end-i_Gaussian-1);
                                        w_array_reverse(i_Gaussian)=w_array(end-i_Gaussian-1);
                                    end
                                    rho_inf_array_reverse(end-1)=rho_inf_array(end);rho_inf_array_reverse(end)=rho_inf_array(end-1);
                                    E_array_reverse(end-1)=E_array(end);E_array_reverse(end)=E_array(end-1);
                                    w_array_reverse(end-1)=w_array(end);w_array_reverse(end)=w_array(end-1);

                                    rho_inf_array=rho_inf_array_reverse;
                                    E_array=-E_array_reverse;
                                    w_array=w_array_reverse;

                                    % EquivalentTriGaussian参数做对称处理
                                    anterior_posterior_ring_reverse=EquivalentTriGaussian_parameter(5:6);
                                    EquivalentTriGaussian_parameter(5:6)=EquivalentTriGaussian_parameter(8:9);
                                    EquivalentTriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                    EquivalentTriGaussian_parameter(4)=-EquivalentTriGaussian_parameter(4);

                                    % TriGaussian参数做对称处理
                                    anterior_posterior_ring_reverse=TriGaussian_parameter(5:6);
                                    TriGaussian_parameter(5:6)=TriGaussian_parameter(8:9);
                                    TriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                    TriGaussian_parameter(4)=-TriGaussian_parameter(4);

                                    y_reverse_AP=zeros(length(y_AP),1);  % 原始数据翻转
                                    for ii=1:length(y_AP)
                                        y_reverse_AP(length(y_AP)-ii+1)=y_AP(ii);
                                    end
                                    y_AP=y_reverse_AP;
                                end

                                % 荧光强度值(纵坐标)归一化
                                % 将拟合完成的多高斯分布在定义域上的数学期望变换为1，同时将原始测量数据和等效三高斯分布按等比例拉伸
                                if length(rho_inf_array)-2==1
                                    f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                        +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2);
                                elseif length(rho_inf_array)-2==2
                                    f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                        +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2);
                                elseif length(rho_inf_array)-2==3
                                    f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                        +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2)...
                                        +rho_inf_array(5)*exp(-1/2*((y-E_array(5))/w_array(5)).^2);
                                elseif length(rho_inf_array)-2==4
                                    f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                        +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2)+...
                                        rho_inf_array(5)*exp(-1/2*((y-E_array(5))/w_array(5)).^2)+rho_inf_array(6)*exp(-1/2*((y-E_array(6))/w_array(6)).^2);
                                elseif length(rho_inf_array)-2>=5
                                    f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                        +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2)...
                                        +rho_inf_array(5)*exp(-1/2*((y-E_array(5))/w_array(5)).^2)+rho_inf_array(6)*exp(-1/2*((y-E_array(6))/w_array(6)).^2)...
                                        +rho_inf_array(7)*exp(-1/2*((y-E_array(7))/w_array(7)).^2);
                                end

                                a0=integral(f,-1/2,1/2);
                                rho0_best=rho0_best/a0;
                                rho_inf_array=rho_inf_array./a0;
                                EquivalentTriGaussian_parameter(2:3:8)=EquivalentTriGaussian_parameter(2:3:8)./a0;EquivalentTriGaussian_parameter(1)=EquivalentTriGaussian_parameter(1)./a0;
                                TriGaussian_parameter(2:3:8)=TriGaussian_parameter(2:3:8)./a0;TriGaussian_parameter(1)=TriGaussian_parameter(1)./a0;
                                y_MultiGaussian_fit_AP=ones(size(x_AP))*rho0_best;  % 多高斯拟合曲线
                                for i_Gaussian=1:length(rho_inf_array)
                                    y_MultiGaussian_fit_AP=y_MultiGaussian_fit_AP+rho_inf_array(i_Gaussian)*exp(-1/2*((x_AP-E_array(i_Gaussian))/w_array(i_Gaussian)).^2);
                                end
                                y_EquivalentTriGaussian_fit_AP=ones(size(x_AP))*EquivalentTriGaussian_parameter(1);  % 等效三高斯拟合曲线
                                for i_Gaussian=1:3
                                    y_EquivalentTriGaussian_fit_AP=y_EquivalentTriGaussian_fit_AP+EquivalentTriGaussian_parameter(3*i_Gaussian-1)*exp(-1/2*((x_AP-EquivalentTriGaussian_parameter(3*i_Gaussian+1))/EquivalentTriGaussian_parameter(3*i_Gaussian)).^2);
                                end
                                y_AP=y_AP./a0;

                                %%% 绘图
                                filename_plot=[filename_photo filename_kongge filename1 filename_kongge filename_cell filename_kongge filename2 filename_kongge filename_protein];

                                figure
                                plot(x_AP,y_AP,'.','MarkerSize',15,'Color','k');  % 绘制数据点
                                hold on;
                                plot(x_AP,y_MultiGaussian_fit_AP,'b','linewidth',3);  % 绘制拟合的高斯曲线
                                plot(x_AP,y_EquivalentTriGaussian_fit_AP,'r','linewidth',3);  % 绘制拟合的高斯曲线
                                hold off;
                                axis ([-0.5 0.5 0 max(y_AP)+0.1]);
                                title(filename_plot)
                                xlabel('x')
                                ylabel('cortex thickness')
                                set(gca,'YTick',0:0.2:max(y_AP)+0.1);  % 修改y坐标刻度
                                set(gca,'XTick',-0.5:0.1:0.5);  % 修改y坐标刻度
                                legend('Experiments','Multi-Gaussian Curve','Tri-Gaussian Curve')
                                set(gca,'FontName','Arial','FontSize',12)
                                saveas(gcf,[framepath2,filename_plot,'.jpg']);
                                close

                                %%% 拟合参数汇总
                                for aa=1
                                    if data_changdubi(number_changdubi)<=0.2
                                        data_changdubi_all(k_position1(k_AP),1,k_AP)=data_changdubi(number_changdubi);
                                        filenamearray_all{k_position1(k_AP),1,k_AP}=filename_AP;
                                        MultiGaussian_data_all(k_position1(k_AP),1,1,k_AP)=rho0_best;
                                        MultiGaussian_data_all(k_position1(k_AP),2:length(rho_inf_array)+1,1,k_AP)=rho_inf_array;
                                        MultiGaussian_data_all(k_position1(k_AP),length(rho_inf_array)+2:2*length(rho_inf_array)+1,1,k_AP)=w_array;
                                        MultiGaussian_data_all(k_position1(k_AP),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,1,k_AP)=E_array;
                                        EquivalentTriGaussian_data_all(k_position1(k_AP),:,1,k_AP)=EquivalentTriGaussian_parameter;
                                        TriGaussian_data_all(k_position1(k_AP),:,1,k_AP)=TriGaussian_parameter;
                                        wucha_best_all(k_position1(k_AP),1,k_AP)=wucha_best;
                                        k_position1(k_AP)=k_position1(k_AP)+1;
                                    elseif data_changdubi(number_changdubi)>0.2 && data_changdubi(number_changdubi)<=0.4
                                        data_changdubi_all(k_position2(k_AP),2,k_AP)=data_changdubi(number_changdubi);
                                        filenamearray_all{k_position2(k_AP),2,k_AP}=filename_AP;
                                        MultiGaussian_data_all(k_position2(k_AP),1,2,k_AP)=rho0_best;
                                        MultiGaussian_data_all(k_position2(k_AP),2:length(rho_inf_array)+1,2,k_AP)=rho_inf_array;
                                        MultiGaussian_data_all(k_position2(k_AP),length(rho_inf_array)+2:2*length(rho_inf_array)+1,2,k_AP)=w_array;
                                        MultiGaussian_data_all(k_position2(k_AP),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,2,k_AP)=E_array;
                                        EquivalentTriGaussian_data_all(k_position2(k_AP),:,2,k_AP)=EquivalentTriGaussian_parameter;
                                        TriGaussian_data_all(k_position2(k_AP),:,2,k_AP)=TriGaussian_parameter;
                                        wucha_best_all(k_position2(k_AP),2,k_AP)=wucha_best;
                                        k_position2(k_AP)=k_position2(k_AP)+1;
                                    elseif data_changdubi(number_changdubi)>0.4 && data_changdubi(number_changdubi)<=0.6
                                        data_changdubi_all(k_position3(k_AP),3,k_AP)=data_changdubi(number_changdubi);
                                        filenamearray_all{k_position3(k_AP),3,k_AP}=filename_AP;
                                        MultiGaussian_data_all(k_position3(k_AP),1,3,k_AP)=rho0_best;
                                        MultiGaussian_data_all(k_position3(k_AP),2:length(rho_inf_array)+1,3,k_AP)=rho_inf_array;
                                        MultiGaussian_data_all(k_position3(k_AP),length(rho_inf_array)+2:2*length(rho_inf_array)+1,3,k_AP)=w_array;
                                        MultiGaussian_data_all(k_position3(k_AP),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,3,k_AP)=E_array;
                                        EquivalentTriGaussian_data_all(k_position3(k_AP),:,3,k_AP)=EquivalentTriGaussian_parameter;
                                        TriGaussian_data_all(k_position3(k_AP),:,3,k_AP)=TriGaussian_parameter;
                                        wucha_best_all(k_position3(k_AP),3,k_AP)=wucha_best;
                                        k_position3(k_AP)=k_position3(k_AP)+1;
                                    elseif data_changdubi(number_changdubi)>0.6 && data_changdubi(number_changdubi)<=0.8
                                        data_changdubi_all(k_position4(k_AP),4,k_AP)=data_changdubi(number_changdubi);
                                        filenamearray_all{k_position4(k_AP),4,k_AP}=filename_AP;
                                        MultiGaussian_data_all(k_position4(k_AP),1,4,k_AP)=rho0_best;
                                        MultiGaussian_data_all(k_position4(k_AP),2:length(rho_inf_array)+1,4,k_AP)=rho_inf_array;
                                        MultiGaussian_data_all(k_position4(k_AP),length(rho_inf_array)+2:2*length(rho_inf_array)+1,4,k_AP)=w_array;
                                        MultiGaussian_data_all(k_position4(k_AP),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,4,k_AP)=E_array;
                                        EquivalentTriGaussian_data_all(k_position4(k_AP),:,4,k_AP)=EquivalentTriGaussian_parameter;
                                        TriGaussian_data_all(k_position4(k_AP),:,4,k_AP)=TriGaussian_parameter;
                                        wucha_best_all(k_position4(k_AP),4,k_AP)=wucha_best;
                                        k_position4(k_AP)=k_position4(k_AP)+1;
                                    else
                                        data_changdubi_all(k_position5(k_AP),5,k_AP)=data_changdubi(number_changdubi);
                                        filenamearray_all{k_position5(k_AP),5,k_AP}=filename_AP;
                                        MultiGaussian_data_all(k_position5(k_AP),1,5,k_AP)=rho0_best;
                                        MultiGaussian_data_all(k_position5(k_AP),2:length(rho_inf_array)+1,5,k_AP)=rho_inf_array;
                                        MultiGaussian_data_all(k_position5(k_AP),length(rho_inf_array)+2:2*length(rho_inf_array)+1,5,k_AP)=w_array;
                                        MultiGaussian_data_all(k_position5(k_AP),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,5,k_AP)=E_array;
                                        EquivalentTriGaussian_data_all(k_position5(k_AP),:,5,k_AP)=EquivalentTriGaussian_parameter;
                                        TriGaussian_data_all(k_position5(k_AP),:,5,k_AP)=TriGaussian_parameter;
                                        wucha_best_all(k_position5(k_AP),5,k_AP)=wucha_best;
                                        k_position5(k_AP)=k_position5(k_AP)+1;
                                    end
                                end
                            elseif k_AP==10  % kk=10表示该细胞仅有一种蛋白的数据文件
                            end

                        elseif exist(framepath_filename_2,'file')
                            number_changdubi=number_changdubi+1;  % 跳过该细胞的长度比数据
                        else  % 若该张照片下，不存在该细胞的数据文件，说明该张照片中的细胞已经拟合完成，跳出j(细胞)循环
                            break
                        end
                        
                        % 覆盖该细胞拟合另一蛋白所用的kk循环参数值
                        filename3=num2str(k);
                        if k==1
                            filename_protein='lifeact';
                        elseif k==2
                            filename_protein='actin';
                        elseif k==3
                            filename_protein='MLC';
                        elseif k==4
                            filename_protein='RhoA';
                        elseif k==5
                            filename_protein='Anillin RBD';
                        elseif k==6
                            filename_protein='VD1';
                        elseif k==7
                            filename_protein='talinA';
                        elseif k==8
                            filename_protein='Cdc42WT';
                        elseif k==9
                            filename_protein='Cdc42D118A';
                        end
                        filename=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3];  % 如 0928_111
                        framepath_filename=[framepath2 filename];
                        framepath_filename_1=[framepath_filename filename_xiahuaxian filename_1 filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111_1.csv
                        framepath_filename_2=[framepath_filename filename_xiahuaxian filename_2 filename_csv];
                        framepath_filename_3=[framepath_filename filename_xiahuaxian filename_3 filename_csv];
                        framepath_filename_danyi=[framepath_filename filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111.csv
                    end
                    if j==1  % 若该张照片下不存在第一个细胞，则说明该日期下所有照片已经处理完成，跳出i(照片)循环
                        break
                    end
                end
            end
        end
    end
end


%%% 删除拟合结果不佳的细胞数据，和赤道轴收缩环特征参数过于偏离平均的细胞,不用于平均A-P轴强度比的计算
% 对wucha_best_all矩阵的每个时期、每种蛋白由小到大进行排序(使用冒泡排序法 bubble sort)，并将顺序同时作用于filenamearray_all data_all矩阵
% 将mix_all和filenamearray_all矩阵合并
mix_all=zeros(1000,25,5,9); % 储存某个时期、某种蛋白下的wucha_best filename data 1-9 changdubi

% mix_all矩阵赋值
for i=1:5 % i表示不同时期
    for j=1:9 % j表示不同蛋白
        for k=1:1000  % 将误差值、文件名、拟合数据输入至mix_all矩阵
            mix_all(k,1,i,j)=wucha_best_all(k,i,j);
            mix_all(k,2,i,j)=0;
            for ii=1:length(MultiGaussian_data_all(1,:,1,1))
                mix_all(k,ii+2,i,j)=MultiGaussian_data_all(k,ii,i,j);
            end
            mix_all(k,end,i,j)=data_changdubi_all(k,i,j);
        end
    end
end

% mix_sort矩阵赋值
mix_sort=mix_all;  % 用于排序的矩阵

% 使用冒泡排序法进行排序
% 首先根据误差大小排序，并将误差值大于8的细胞数据剔除
cell_number=zeros(5,9);
for i=1:5
    for j=1:9
        for k=1:1000  % 找出该矩阵中的实际数据组数
            if mix_sort(k,1,i,j)==0
                break
            end
        end
        cell_number(i,j)=k-1; % 将实际数组数(细胞个数)赋值给cell_number
        
        % 冒泡排序
        for ii=2:cell_number(i,j)
            for jj=cell_number(i,j):-1:ii
                if mix_sort(jj,1,i,j)<mix_sort(jj-1,1,i,j)  % 若前一个细胞的误差大于后一个细胞的误差，则调换顺序
                    temp=mix_sort(jj,:,i,j);
                    mix_sort(jj,:,i,j)=mix_sort(jj-1,:,i,j);
                    mix_sort(jj-1,:,i,j)=temp;
                end
            end
        end
        
        % 删除误差值大于8的细胞数据
        for k=1:cell_number(i,j)
            if mix_sort(k,1,i,j)>=8
                for ii=1:length(mix_sort(1,:,1,1))
                    mix_sort(k,ii,i,j)=0;  % 剔除误差大于8的细胞数据
                end
                cell_number(i,j)=cell_number(i,j)-1;
            end
        end
    end
end

% 计算basal domain多高斯的平均强度和聚集程度，储存到mix_sort_basalave矩阵中，用于后续排序查找需删除数据
mix_sort_basalave=zeros(1000,13,5,9);
for i=1:5
    for j=1:9
        for k=1:cell_number(i,j)
            for ii=4:length(mix_sort(1,:,1,1))
                if mix_sort(k,ii,i,j)==0  % 若第ii位为0，说明basal domain有(ii-10)/3个高斯
                    break
                end
            end
            if mix_sort(k,end-1,i,j)~=0  % 若最后一位不为0，说明basal domain有5(max)个高斯
                ii=25;
            end
            basal_Gaussian_number=(ii-10)/3;

            % 计算basal domain多高斯的平均强度和聚集程度
            basal_overactivity_average=sum(mix_all(k,4:basal_Gaussian_number+3,i,j))/basal_Gaussian_number;
            basal_signalwidth_average=sum(mix_all(k,basal_Gaussian_number+4:2*basal_Gaussian_number+3,i,j))/basal_Gaussian_number;
            basal_expectation_average=sum(mix_all(k,2*basal_Gaussian_number+4:3*basal_Gaussian_number+3,i,j))/basal_Gaussian_number;

            % 赋值给mix_sort_basalave
            for ii=1:length(mix_sort_basalave(1,:,1,1))
                if ii==4
                    mix_sort_basalave(k,ii,i,j)=basal_overactivity_average;
                elseif ii==5
                    mix_sort_basalave(k,ii,i,j)=basal_signalwidth_average;
                elseif ii==6
                    mix_sort_basalave(k,ii,i,j)=basal_expectation_average;
                elseif ii==13
                    mix_sort_basalave(k,ii,i,j)=mix_sort(k,end,i,j);
                elseif ii<4
                    mix_sort_basalave(k,ii,i,j)=mix_sort(k,ii,i,j);
                elseif ii>=7 && ii<=9
                    mix_sort_basalave(k,ii,i,j)=mix_sort(k,(ii-6)*(basal_Gaussian_number+2)+2,i,j);
                else
                    mix_sort_basalave(k,ii,i,j)=mix_sort(k,(ii-9)*(basal_Gaussian_number+2)+3,i,j);
                end
            end
        end
    end
end

% 分别按中心收缩环强度和聚集程度(basal domain多高斯的平均强度和聚集程度),和前侧收缩环的强度(一般两侧环拟合不佳的细胞其聚集程度和强度相关性强)排序，
% 剔除三者前10%和后10%的细胞数据的并集
memory_cell_wucha=zeros(200,5,9,3);  % 储存待删除细胞数据对应的误差值
k_cellname=ones(5,9,3);  % memory_cellname的数组长度
number_sort_position=[4 5 7];  % 特征参数对应在mix_all矩阵的位置
for k=1:3
    number_sort=number_sort_position(k);
    for i=1:5
        for j=1:9
            % 冒泡排序
            for ii=2:cell_number(i,j)
                for jj=cell_number(i,j):-1:ii
                    if mix_sort_basalave(jj,number_sort_position(k),i,j)<mix_sort_basalave(jj-1,number_sort_position(k),i,j)  % 若前一个细胞的特征参数大于后一个细胞的特征参数，则调换顺序
                        temp=mix_sort(jj,:,i,j);
                        mix_sort(jj,:,i,j)=mix_sort(jj-1,:,i,j);
                        mix_sort(jj-1,:,i,j)=temp;
                        temp2=mix_sort_basalave(jj,:,i,j);
                        mix_sort_basalave(jj,:,i,j)=mix_sort_basalave(jj-1,:,i,j);
                        mix_sort_basalave(jj-1,:,i,j)=temp2;
                    end
                end
            end
            
            if k~=3
                % 找出前10%的细胞，用误差值标记
                for jj=1:floor(cell_number(i,j)/10)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
                % 找出后10%的细胞，用误差值标记
                for jj=cell_number(i,j)-floor(cell_number(i,j)/10)+1:cell_number(i,j)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
            else
                % 找出前10%的细胞，用误差值标记
                for jj=1:floor(cell_number(i,j)/10)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
                % 找出后10%的细胞，用误差值标记
                for jj=cell_number(i,j)-floor(cell_number(i,j)/5)+1:cell_number(i,j)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
            end
        end
    end
end
% 取三次排序得到的细胞的并集，并在mix_all矩阵中删去
mix_sort2=mix_sort;  % 将mix_all数组内容赋值给mix_all2
mix_sort=zeros(1000,25,5,9);  % 清空数组，用于后续存放删除非正常细胞后的细胞拟合数据
for i=1:5
    for j=1:9
        kk2=1;  % 用于存储mix_all每个时期，每种蛋白的细胞数
        for k1=1:1000
            if mix_sort2(k1,1,i,j)~=0  % 存在第k1个细胞
                logic=0;  % 用于判断该细胞是否需要删除
                for k2=1:k_cellname(i,j,1)
                    if mix_sort2(k1,1,i,j)==memory_cell_wucha(k2,i,j,1)
                        logic=1;  % 若该细胞的wucha值为赤道轴收缩环中心强度,和前侧收缩环的强度的前10%和后10%的细胞，则说明该细胞需要删除
                        break
                    end
                end
                if logic==0  % 若第一次未找到有相同误差值的细胞，再进行第二次查找
                    for k2=1:k_cellname(i,j,2)
                        if mix_sort2(k1,1,i,j)==memory_cell_wucha(k2,i,j,2)
                            logic=1;  % 若该细胞的wucha值为赤道轴收缩环中心聚集程度,和前侧收缩环的强度的前10%和后10%的细胞，则说明该细胞需要删除
                            break
                        end
                    end
                end
                if logic==0  % 若第一、二次未找到有相同误差值的细胞，再进行第三次查找
                    for k2=1:k_cellname(i,j,3)
                        if mix_sort2(k1,1,i,j)==memory_cell_wucha(k2,i,j,3)
                            logic=1;  % 若该细胞的wucha值为赤道轴收缩环中心聚集程度,和前侧收缩环的强度的前10%和后10%的细胞，则说明该细胞需要删除
                            break
                        end
                    end
                end
                if logic==0  % 若logic为0，则表示该细胞不需删除，将其赋值给mix_sort filenamearray_all矩阵
                    mix_sort(kk2,:,i,j)=mix_sort2(k1,:,i,j);
                    kk2=kk2+1;
                else  % 若logic为1，则表示该细胞需删除，不赋值给mix_all矩阵，且cell_number减1
                    cell_number(i,j)=cell_number(i,j)-1;
                end
            else
                break
            end
        end
    end
end

% 将排完序并删减完成的mix_sort矩阵重新赋值给mix_all矩阵
mix_all=mix_sort;

% 重新计算basal domain平均后的拟合参数矩阵，并对wucha_best排序
mix_all_basalave=zeros(1000,13,5,9);
for i=1:5
    for j=1:9
        for k=1:cell_number(i,j)
            for ii=4:length(mix_all(1,:,1,1))
                if mix_all(k,ii,i,j)==0  % 若第ii位为0，说明basal domain有(ii-10)/3个高斯
                    break
                end
            end
            if mix_all(k,end-1,i,j)~=0  % 若最后一位不为0，说明basal domain有5(max)个高斯
                ii=25;
            end
            basal_Gaussian_number=(ii-10)/3;

            % 计算basal domain多高斯的平均强度和聚集程度
            basal_overactivity_average=sum(mix_all(k,4:basal_Gaussian_number+3,i,j))/basal_Gaussian_number;
            basal_signalwidth_average=sum(mix_all(k,basal_Gaussian_number+4:2*basal_Gaussian_number+3,i,j))/basal_Gaussian_number;
            basal_expectation_average=sum(mix_all(k,2*basal_Gaussian_number+4:3*basal_Gaussian_number+3,i,j))/basal_Gaussian_number;

            % 赋值给mix_sort_basalave
            for ii=1:length(mix_all_basalave(1,:,1,1))
                if ii==4
                    mix_all_basalave(k,ii,i,j)=basal_overactivity_average;
                elseif ii==5
                    mix_all_basalave(k,ii,i,j)=basal_signalwidth_average;
                elseif ii==6
                    mix_all_basalave(k,ii,i,j)=basal_expectation_average;
                elseif ii==13
                    mix_all_basalave(k,ii,i,j)=mix_all(k,end,i,j);
                elseif ii<4
                    mix_all_basalave(k,ii,i,j)=mix_all(k,ii,i,j);
                elseif ii>=7 && ii<=9
                    mix_all_basalave(k,ii,i,j)=mix_all(k,(ii-6)*(basal_Gaussian_number+2)+2,i,j);
                else
                    mix_all_basalave(k,ii,i,j)=mix_all(k,(ii-9)*(basal_Gaussian_number+2)+3,i,j);
                end
            end
        end
    end
end
% 冒泡排序
for i=1:5
    for j=1:9
        for ii=2:cell_number(i,j)
            for jj=cell_number(i,j):-1:ii
                if mix_all_basalave(jj,1,i,j)<mix_all_basalave(jj-1,1,i,j)  % 若前一个细胞的误差大于后一个细胞的误差，则调换顺序
                    temp=mix_all_basalave(jj,:,i,j);
                    mix_all_basalave(jj,:,i,j)=mix_all_basalave(jj-1,:,i,j);
                    mix_all_basalave(jj-1,:,i,j)=temp;
                    temp2=mix_all(jj,:,i,j);
                    mix_all(jj,:,i,j)=mix_all(jj-1,:,i,j);
                    mix_all(jj-1,:,i,j)=temp2;
                end
            end
        end
    end
end


%%% 计算A-P轴前后收缩环中心强度的比值(用于后续照片A-P polarity direction判断和叠加观测数据拆分)
AP_overactivity_ratio=zeros(5,9);  % 定义A-P轴收缩环强度比变量，每个时期、每种蛋白单独计算
max_AP_overactivity_ratio=zeros(5,9);  % 定义A-P轴收缩环强度比最大值变量
for i=1:5
    for j=1:9
        for k=1:1000
            if abs(mix_all_basalave(k,4,i,j))+abs(mix_all_basalave(k,7,i,j))+abs(mix_all_basalave(k,10,i,j))~=0  % 存在该细胞的拟合参数
                if mix_all_basalave(k,10,i,j)~=0  % 首先计算A-P比值存在的总数
                    AP_overactivity_ratio(i,j)=AP_overactivity_ratio(i,j)+mix_all_basalave(k,7,i,j)/mix_all_basalave(k,10,i,j);  % 找出最大的A-P轴强度比值
                    if max_AP_overactivity_ratio(i,j)<mix_all_basalave(k,7,i,j)/mix_all_basalave(k,10,i,j)
                        max_AP_overactivity_ratio(i,j)=mix_all_basalave(k,7,i,j)/mix_all_basalave(k,10,i,j);
                    end
                end
            else  % 不存在该细胞的拟合参数，跳出循环
                break
            end
            if wucha_best_all(k,i,j)>180  % 表示该细胞以及接下来的细胞拟合结果较差，不用于计算
                break
            end
        end
        k=k-1;
        if k==0
            AP_overactivity_ratio(i,j)=0;
        else
            for k_2=1:k  % 若A-P轴强度比值不存在，则将其替换为所测得的最大值计入
                if mix_all_basalave(k,10,i,j)==0
                    AP_overactivity_ratio(i,j)=AP_overactivity_ratio(i,j)+max_AP_overactivity_ratio(i,j);
                end
            end
            AP_overactivity_ratio(i,j)=AP_overactivity_ratio(i,j)/k;
        end
    end
end


%%% 计算A-P轴PCP强度比(C_AP) (C_AP定义：前侧收缩环的高斯强度除以后侧收缩环的高斯强度)
%%% protein_all_fit的第三维为rho_anterior和rho_posterior的关系，此关系即为C_AP
%%% 由于前侧收缩环和后侧收缩环的rho和w相关性较强，考虑将w表示为rho的函数，即C_AP表示为仅关于rho的函数(若前后两个细胞所处时期不同，则按C_AP_next/(C_AP_next+1)*(C_AP_now+1)计算)
%%% 并计算w关于rho回归函数的标准差

% 计算rho_anterior和rho_posterior比值随时间的变化
AP_rho_ratio_time=zeros(2000,9,2);  % 定义A-P polarity coefficient
for j=1:9
    if sum(cell_number(:,j))<=10
        continue
    end
    
    AP_rho_inf_array=zeros(sum(cell_number(:,j)),3);
    k_position=1;
    i_delete=0;
    for i=1:5
        for jj=1:cell_number(i,j)
            if mix_all_basalave(jj,10,i,j)<0.01
                i_delete=i_delete+1;
                continue
            end
            AP_rho_inf_array(k_position,1)=mix_all_basalave(jj,13,i,j);
            AP_rho_inf_array(k_position,2)=mix_all_basalave(jj,7,i,j);
            AP_rho_inf_array(k_position,3)=mix_all_basalave(jj,10,i,j);
            k_position=k_position+1;
        end
    end
    AP_rho_inf_array_fuben=AP_rho_inf_array(1:sum(cell_number(:,j))-i_delete,:);
    AP_rho_inf_array=AP_rho_inf_array_fuben;
    
    % 冒泡排序(根据TD/CD)
    for ii=2:sum(cell_number(:,j))-i_delete
        for jj=sum(cell_number(:,j))-i_delete:-1:ii
            if AP_rho_inf_array(jj,1)<AP_rho_inf_array(jj-1,1)
                temp=AP_rho_inf_array(jj,:);
                AP_rho_inf_array(jj,:)=AP_rho_inf_array(jj-1,:);
                AP_rho_inf_array(jj-1,:)=temp;
            end
        end
    end
    
    AP_rho_ratio_temp=zeros(sum(cell_number(:,j))-i_delete,2);
    AP_rho_ratio_MA=zeros(sum(cell_number(:,j))-i_delete,2);
    win=9;
    for i=1:sum(cell_number(:,j))-i_delete
        AP_rho_ratio_temp(i,1)=AP_rho_inf_array(i,1);
        AP_rho_ratio_temp(i,2)=AP_rho_inf_array(i,2)/AP_rho_inf_array(i,3);
    end
    for i=1:sum(cell_number(:,j))-i_delete
        if i<=(win-1)/2  % 若t为初始若干时刻
            AP_rho_ratio_MA(i,1)=sum(AP_rho_ratio_temp(1:2*i-1,1))/(2*i-1);
            AP_rho_ratio_MA(i,2)=sum(AP_rho_ratio_temp(1:2*i-1,2))/(2*i-1);
        elseif i<sum(cell_number(:,j))-i_delete-(win-1)/2+1  % 若t为中间时刻
            AP_rho_ratio_MA(i,1)=sum(AP_rho_ratio_temp(i-(win-1)/2:i+(win-1)/2,1))/win;
            AP_rho_ratio_MA(i,2)=sum(AP_rho_ratio_temp(i-(win-1)/2:i+(win-1)/2,2))/win;
        else  % 若t为末尾若干时刻
            temp=sum(cell_number(:,j))-i_delete-i+1;
            AP_rho_ratio_MA(i,1)=sum(AP_rho_ratio_temp(sum(cell_number(:,j))-i_delete-2*temp+2:sum(cell_number(:,j))-i_delete,1))/(2*temp-1);
            AP_rho_ratio_MA(i,2)=sum(AP_rho_ratio_temp(sum(cell_number(:,j))-i_delete-2*temp+2:sum(cell_number(:,j))-i_delete,2))/(2*temp-1);
        end
    end
    
    N=1;  % 拟合阶数
    [p,S]=polyfit(AP_rho_ratio_MA(:,1),AP_rho_ratio_MA(:,2),N);
    [yfit,delta]=polyconf(p,AP_rho_ratio_MA(:,1),S,'alpha',0.05,'predopt','curve');
    x_dummy=linspace(min(AP_rho_ratio_MA(:,1)),max(AP_rho_ratio_MA(:,1)),2000);  % 拟合函数的x坐标
    delta_interp=interp1(AP_rho_ratio_MA(:,1),delta,x_dummy,'linear');  % 拟合曲线置信区间插值
    AP_rho_ratio_time(:,j,1)=x_dummy;
    for i=1:N+1
        AP_rho_ratio_time(:,j,2)=AP_rho_ratio_time(:,j,2)+p(i).*x_dummy'.^(N-i+1);
    end
end

% 将一种蛋白的5个时期数据整合到一维数组中(protein_mix_all)
protein_mix_all=zeros(5000,13,9);  % 每种蛋白的所有时期特征数据集合矩阵
for j=1:9
    k_position=1;  % 赋值到protein_mix_all矩阵的位置
    for i=1:5
        protein_mix_all(k_position:k_position+cell_number(i,j)-1,:,j)=mix_all_basalave(1:cell_number(i,j),:,i,j);
        k_position=k_position+cell_number(i,j);
    end
end

% 拟合前侧收缩环的rho和w的函数关系，后侧收缩环的rho和w的函数关系，前侧收缩环的rho和后侧收缩环的rho的函数关系
% 首先对待拟合参数(对应x轴的)进行排序，然后取滑动平均并计算滑动标准差，将滑动平均值和滑动标准差数值插值到同一长度，最后拟合函数关系
parameter_x=[7 10 7];  % 三组待拟合函数关系的x轴参数位置
parameter_y=[8 11 10];  % 三组待拟合函数关系的y轴参数位置
win=1;  % 滑动平均窗口大小
protein_all_MA_interp=zeros(2000,9,2,3);  % 每种蛋白的滑动平均值插值矩阵，第一维为参数滑动平均值的插值结果，第二维为不同的蛋白，第三维为x轴和y轴，第四维为三组待拟合函数关系的y轴参数
std_protein_all_MA_interp=zeros(2000,9,2,3);  % 每种蛋白的滑动标准差插值矩阵，第一维为参数滑动标准差的插值结果，第二维为不同的蛋白，第三维为x轴和y轴，第四维为三组待拟合函数关系的y轴参数
protein_all_fit=zeros(2000,9,3);  % 每种蛋白的多项式拟合函数，第一维为拟合函数(y关于x的，例如w关于rho的函数)，第二维为不同的蛋白，第三维为三组待拟合函数关系
protein_all_delta=zeros(2000,9,3);  % 每种蛋白的多项式拟合函数的置信区间，第一维为拟合函数(y关于x的，例如w关于rho的函数)，第二维为不同的蛋白，第三维为三组待拟合函数关系

for k=1:3
    % 冒泡排序(依据待拟合函数关系的x轴参数)
    for j=1:9  % 9种不同蛋白
        for ii=2:sum(cell_number(:,j))
            for jj=sum(cell_number(:,j)):-1:ii
                if protein_mix_all(jj,parameter_x(k),j)<protein_mix_all(jj-1,parameter_x(k),j)  % 若前一个细胞的参数大于后一个细胞的参数，则调换顺序
                    temp=protein_mix_all(jj,:,j);
                    protein_mix_all(jj,:,j)=protein_mix_all(jj-1,:,j);
                    protein_mix_all(jj-1,:,j)=temp;
                end
            end
        end
    end
    
    % 将特征参数取滑动平均(过滤高频震荡)
    % 计算滑动平均值和滑动标准差
    protein_all_MA=zeros(sum(cell_number(:,i)),9,2);  % 每种蛋白的滑动平均值，第一列为x轴对应的拟合参数，第二列为y轴对应的拟合参数
    std_protein_all_MA=zeros(sum(cell_number(:,i)),9,2);  % 每种蛋白的滑动标准差，第一列储存x轴参数的滑动平均值，第二列储存y轴参数的滑动标准差
    for i=1:9  % 9种不同蛋白
        for j=1:sum(cell_number(:,i))
            if j<=(win-1)/2  % 若t为初始若干时刻
                protein_all_MA(j,i,1)=sum(protein_mix_all(1:2*j,parameter_x(k),i))/(2*j);
                protein_all_MA(j,i,2)=sum(protein_mix_all(1:2*j,parameter_y(k),i))/(2*j);
            elseif j<sum(cell_number(:,i))-(win-1)/2+1  % 若t为中间时刻
                protein_all_MA(j,i,1)=sum(protein_mix_all(j-(win-1)/2:j+(win-1)/2,parameter_x(k),i))/win;
                protein_all_MA(j,i,2)=sum(protein_mix_all(j-(win-1)/2:j+(win-1)/2,parameter_y(k),i))/win;
            else  % 若t为末尾若干时刻
                temp=sum(cell_number(:,i))-j+1;
                protein_all_MA(j,i,1)=sum(protein_mix_all(sum(cell_number(:,i))-2*temp-1:sum(cell_number(:,i)),parameter_x(k),i))/(2*temp+2);
                protein_all_MA(j,i,2)=sum(protein_mix_all(sum(cell_number(:,i))-2*temp-1:sum(cell_number(:,i)),parameter_y(k),i))/(2*temp+2);
            end
        end
    end
    for i=1:9  % 9种不同蛋白
        for j=1:sum(cell_number(:,i))
            if j<=(win-1)/2  % 若t为初始若干时刻
                std_protein_all_MA(j,i,1)=protein_all_MA(j,i,1);
                std_protein_all_MA(j,i,2)=std(protein_all_MA(1:2*j,i,2));
            elseif j<sum(cell_number(:,i))-(win-1)/2+1  % 若t为中间时刻
                std_protein_all_MA(j,i,1)=protein_all_MA(j,i,1);
                std_protein_all_MA(j,i,2)=std(protein_all_MA(j-(win-1)/2:j+(win-1)/2,i,2));
            else  % 若t为末尾若干时刻
                temp=sum(cell_number(:,i))-j+1;
                std_protein_all_MA(j,i,1)=protein_all_MA(j,i,1);
                std_protein_all_MA(j,i,2)=std(protein_all_MA(sum(cell_number(:,i))-2*temp-1:sum(cell_number(:,i)),i,2));
            end
        end
    end
    
    % 滑动标准差插值
    for i=1:9
        if sum(cell_number(:,i))-win+1<=0
            continue
        end
        
        x_interp=linspace(min(protein_all_MA(1:1:sum(cell_number(:,i)),i,1)),max(protein_all_MA(1:1:sum(cell_number(:,i)),i,1)),2000);
        protein_all_MA_interp(:,i,1,k)=interp1(protein_all_MA(1:1:sum(cell_number(:,i)),i,1),protein_all_MA(1:1:sum(cell_number(:,i)),i,1),x_interp,'linear');  % 拟合参数滑动平均标准差插值
        protein_all_MA_interp(:,i,2,k)=interp1(protein_all_MA(1:1:sum(cell_number(:,i)),i,1),protein_all_MA(1:1:sum(cell_number(:,i)),i,2),x_interp,'linear');  % 拟合参数滑动平均标准差插值
        std_protein_all_MA_interp(:,i,1,k)=interp1(std_protein_all_MA(1:1:sum(cell_number(:,i)),i,1),std_protein_all_MA(1:1:sum(cell_number(:,i)),i,1),x_interp,'linear');  % 拟合参数滑动平均标准差插值
        std_protein_all_MA_interp(:,i,2,k)=interp1(std_protein_all_MA(1:1:sum(cell_number(:,i)),i,1),std_protein_all_MA(1:1:sum(cell_number(:,i)),i,2),x_interp,'linear');  % 拟合参数滑动平均标准差插值
    end
    
    % 多项式拟合特征参数间的函数关系
    N=2;  % 拟合阶数
    for i=1:9
        if sum(cell_number(:,i))-win+1<=0
            continue
        end
        
        [p,S]=polyfit(protein_all_MA(1:1:sum(cell_number(:,i)),i,1),protein_all_MA(1:1:sum(cell_number(:,i)),i,2),N);
        [yfit,delta]=polyconf(p,protein_all_MA(1:1:sum(cell_number(:,i)),i,1),S,'alpha',0.05,'predopt','curve');  % 拟合函数的95%置信区间
        x_fit=linspace(min(protein_all_MA(1:1:sum(cell_number(:,i)),i,1)),max(protein_all_MA(1:1:sum(cell_number(:,i)),i,1)),2000);  % 拟合函数的x坐标
        delta_interp=interp1(protein_all_MA(1:1:sum(cell_number(:,i)),i,1),delta,x_fit,'linear');  % 拟合曲线置信区间插值
        protein_all_delta(:,i,k)=delta_interp;
        for ii=1:N+1
            protein_all_fit(:,i,k)=protein_all_fit(:,i,k)+p(ii).*x_fit'.^(N-ii+1);
        end
    end
end

f=fittype('a*x');
fit1=fit(protein_all_MA_interp(:,1,2,3),protein_all_MA_interp(:,1,1,3),f);
conf1=predint(fit1,protein_all_MA_interp(:,1,2,3),0.99,'functional','off');
C_AP=fit1.a;  % 正比例函数拟合的C_AP


% i=1;k=2;
% % 画图
% for ii=1:2000
%     plot(protein_all_MA_interp(ii,i,1,k),protein_all_MA_interp(ii,i,2,k),'.','MarkerSize',20,'color','r')  % 画均值点
%     hold on
% end
% 
% % 画滑动平均标准差
% for ii=1:2:2000
%     aaa=plot([std_protein_all_MA_interp(ii,i,1,k) std_protein_all_MA_interp(ii,i,1,k)],[protein_all_fit(ii,i,k)-std_protein_all_MA_interp(ii,i,2,k) protein_all_fit(ii,i,k)+std_protein_all_MA_interp(ii,i,2,k)],'-','Color','b','linewidth',1);  % 画TD误差bar竖线
%     aaa.Color(4) = 0.2;
%     hold on
% end
% plot(std_protein_all_MA_interp(:,i,1,k),protein_all_fit(:,i,k),'-b','linewidth',3)



%% Step 2
for k=1:9  % 1表示lifeact，2表示actin，3表示MLC，4表示RhoA，5表示Anillin RBD，6表示VD1,7表示talinA，8表示Cdc42WT，9表示Cdc42D118A
    filename3=num2str(k);
    
    if k==1
        filename_protein='lifeact';
    elseif k==2
        filename_protein='actin';
    elseif k==3
        filename_protein='MLC';
    elseif k==4
        filename_protein='RhoA';
    elseif k==5
        filename_protein='Anillin RBD';
    elseif k==6
        filename_protein='VD1';
    elseif k==7
        filename_protein='talinA';
    elseif k==8
        filename_protein='Cdc42WT';
    elseif k==9
        filename_protein='Cdc42D118A';
    end
    filename_photo='photo';
    filename_cell='cell';
    filename_kongge=' ';
    
    for month=1:12
        if month<10
            filename_month=[filename_0 num2str(month)];
        else
            filename_month=num2str(month);
        end
        for day=1:31
            if day<10
                filename_day=[filename_0 num2str(day)];
            else
                filename_day=num2str(day);
            end
            number_changdubi=0;  % 用于输出文件时指示长度比文件中的数据位置
            
            framepath2=[framepath num2str(k) filename_fanxiegang filename_month filename_day filename_fanxiegang];  % 子文件夹路径，如 'C:\Users\Desktop\1\0928'
            framepath_filename_changdubi=[framepath2 filename_changdubi filename_xlsx];
            framepath_filename_Radius=[framepath2 filename_Radius filename_xlsx];
            framepath_filename_ContactAngle=[framepath2 filename_ContactAngle filename_xlsx];

            if exist(framepath_filename_changdubi,'file')   % 长度比文件是否存在可用于快速筛选是否存在该天拍摄该蛋白的数据，缩短查找时间
                data_lumenradius=xlsread(framepath_filename_Radius,'sheet1','A1:A100');  % 读入lumen半径
                data_lumenradius=sqrt(data_lumenradius./pi);
                data_lumen_ContactAngle=xlsread(framepath_filename_ContactAngle,'sheet1','A1:A100');  % 读入lumen接触角
                data_lumenlength=xlsread(framepath_filename_changdubi,'sheet1','A1:A100');  % 读入该文件下的lumen长度
                data_celllength=xlsread(framepath_filename_changdubi,'sheet1','B1:B100');  % 读入该文件下的cell长度
                % 若Lumen length大于cell length，则cell length=Lumen length+0.5
                for i=1:length(data_lumenradius)
                    if 2*data_lumenradius(i)*sin(data_lumen_ContactAngle(i)*pi./180)>data_celllength(i)
                        data_celllength(i)=2*data_lumenradius(i)*sin(data_lumen_ContactAngle(i)*pi./180)+0.5;
                    end
                end
                data_changdubi=2.*data_lumenradius.*sin(data_lumen_ContactAngle.*pi./180)./data_celllength;
                
                for i=1:30  % 外侧循环用于查找目标文件是否存在(i表示第几张photo，j表示photo中的第几个细胞，k表示荧光通道)
                    filename1=num2str(i);
                    %%% 确定A-P轴极性方向
                    %%% 每张照片确定一个统一的A-P轴极性方向
                    %%% 步骤：1. 首先找出该张照片下的所有细胞名称，并把每个细胞的名称和原始数据分别保存在一个矩阵中
                    %%%       2. 查找该张照片中有无无相邻细胞的细胞数据
                    %%%         (2.1)若有，则根据无相邻细胞的细胞数据的A-P轴极性方向确定该张照片的极性方向(若有多个细胞
                    %%%            且极性方向不同，则根据大多数细胞的极性方向确定该张照片的极性方向)
                    %%%         (2.2)若无，分别尝试两种极性方向，选择该张照片下所有细胞按该种极性方向处理后，所得大部分
                    %%%            细胞的实际极性方向相同且平均A-P强度比更大的极性方向
                    filenamearray_PD=cell(10,1);  % PD：polarity determination
                    y_danyi=zeros(1500,10);  % 存储单一数据文件的细胞数据
                    y_part1=zeros(1500,10);  % 存储有多个数据文件的细胞数据的第一部分
                    y_part2=zeros(1500,10);  % 存储有多个数据文件的细胞数据的第二部分
                    y_part3=zeros(1500,10);  % 存储有多个数据文件的细胞数据的第三部分
                    number_monocell=0;  % 表示有单一数据文件的细胞个数
                    number_multicell=0;  % 表示有多个数据文件的细胞个数
                    polarity_direction=0;  % 用于储存A-P轴极性方向(1表示正向，2表示反向)
                    
                    %%% 1. 首先找出该张照片下的所有细胞名称，并把每个细胞的名称和原始数据分别保存在一个矩阵中
                    %%% 将该张照片下的所有细胞数据分别保存到y_danyi矩阵和y_part1/2/3矩阵中，文件名保存到filenamearray_PD矩阵中
                    for j=1:10
                        filename2=num2str(j);
                        filename=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3];  % 如 0928_111
                        framepath_filename=[framepath2 filename];
                        
                        % 拼接完整的文件路径
                        framepath_filename_1=[framepath_filename filename_xiahuaxian filename_1 filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111_1.csv
                        framepath_filename_2=[framepath_filename filename_xiahuaxian filename_2 filename_csv];
                        framepath_filename_3=[framepath_filename filename_xiahuaxian filename_3 filename_csv];
                        framepath_filename_danyi=[framepath_filename filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111.csv
                        
                        % 细胞数据和文件名矩阵赋值
                        if exist(framepath_filename_danyi,'file')   % 一个细胞的数据有一个文件
                            number_monocell=number_monocell+1;
                            filenamearray_PD{number_monocell+number_multicell}=filename;  % 文件名赋值
                            
                            data=xlsread(framepath_filename_danyi); % 导入数据
                            [X,Y]=size(data); % 丈量数据矩阵尺寸
                            for ii=1:X  % 数据矩阵保存到y_danyi矩阵中
                                y_danyi(ii,number_monocell)=data(ii,2);
                            end
                        elseif exist(framepath_filename_1,'file')   % 一个细胞的数据有多个文件
                            number_multicell=number_multicell+1;
                            filenamearray_PD{number_monocell+number_multicell}=filename;  % 文件名赋值
                            
                            data_1=xlsread(framepath_filename_1); % 导入数据
                            data_2=xlsread(framepath_filename_2);
                            data_3=xlsread(framepath_filename_3);
                            for ii=1:length(data_1)  % 数据矩阵保存到y_part1/2/3矩阵中
                                y_part1(ii,number_multicell)=data_1(ii,2);
                            end
                            for ii=1:length(data_2)  % 数据矩阵保存到y_part1/2/3矩阵中
                                y_part2(ii,number_multicell)=data_2(ii,2);
                            end
                            for ii=1:length(data_3)  % 数据矩阵保存到y_part1/2/3矩阵中
                                y_part3(ii,number_multicell)=data_3(ii,2);
                            end
                        end
                    end
                    
                    %%% 确定A-P轴极性方向
                    if number_multicell==0  % 该张照片下不存在有相邻细胞的细胞数据，直接结束该次计算，进入下一次计算
                        number_changdubi=number_changdubi+number_monocell;
                        continue
                    else
                        if number_monocell>=1
                         %%% 2. 查找该张照片中有无无相邻细胞的细胞数据
                         %%%    (2.1)若有，则根据无相邻细胞的细胞数据的A-P轴极性方向确定该张照片的极性方向(若有多个细胞
                         %%%         且极性方向不同，则根据大多数细胞的极性方向确定该张照片的极性方向)
                            AP_overactivity_ratio_PD=zeros(number_monocell,1);  % 存放所有无相邻细胞的细胞的A-P收缩环强度比
                            
                            % 确定无相邻细胞的细胞极性方向
                            % 首先通过lateral domain的收缩环高斯分布拟合获得A-P极性大小
                            for i_momocell=1:number_monocell
                                for X=length(y_danyi(:,i_momocell)):-1:1  % 确定y数据矩阵长度
                                    if y_danyi(X,i_momocell)~=0
                                        break
                                    end
                                end
                                X=X+1;
                                % X轴关于坐标原点对称
                                x=linspace(-0.5,0.5,X);
                                % Y轴数据归一化
                                y_danyi(:,i_momocell)=y_danyi(:,i_momocell)./max(y_danyi(:,i_momocell));
                                
                                Y_anterior=zeros(1,round((X-1)*0.2));  % 定义细胞两侧的收缩环灰度值数据，用于单高斯分布拟合确定极性方向
                                Y_posterior=zeros(1,round((X-1)*0.2));
                                
                                % Y_anterior Y_posterior矩阵赋值
                                for ii=1:round((X-1)*0.2)
                                    Y_anterior(ii)=y_danyi(ii,i_momocell);
                                end
                                for ii=X-round((X-1)*0.2)+1:X
                                    Y_posterior(ii-X+round((X-1)*0.2))=y_danyi(ii,i_momocell);
                                end
                                
                                % 最小二乘拟合两侧高斯分布，确定要删去的部分
                                % 拟合Y_anterior
                                wucha_best=100000000;
                                for rho0=0.01:(1-0.01)/25:1
                                    for rho_inf_anterior=0.01:(5-0.01)/25:5
                                        for w=0.01:(1-0.01)/25:1
                                            for E=-0.5:(x(round((X-1)*0.2))+0.5)/25:x(round((X-1)*0.2))
                                                wucha=0;
                                                for ii=1:round((X-1)*0.2)
                                                    wucha=wucha+(rho0+rho_inf_anterior*exp(-1/2*((x(ii)-E)/w)^2)-Y_anterior(ii))^2;
                                                end
                                                if wucha<wucha_best
                                                    wucha_best=wucha;
                                                    rho0_best=rho0;
                                                    rho_inf_anterior_best=rho_inf_anterior;
                                                    w_best=w;
                                                    E_best=E;
                                                end
                                            end
                                        end
                                    end
                                end
                                
                                % 拟合Y_posterior
                                wucha_best=100000000;
                                for rho0=0.01:(1-0.01)/25:1
                                    for rho_inf_posterior=0.01:(5-0.01)/25:5
                                        for w=0.01:(1-0.01)/25:1
                                            for E=x(X-round((X-1)*0.2)+1):(0.5-x(X-round((X-1)*0.2)+1))/25:0.5
                                                wucha=0;
                                                for ii=X-round((X-1)*0.2)+1:X
                                                    wucha=wucha+(rho0+rho_inf_posterior*exp(-1/2*((x(ii)-E)/w)^2)-Y_posterior(ii-X+round((X-1)*0.2)))^2;
                                                end
                                                if wucha<wucha_best
                                                    wucha_best=wucha;
                                                    rho0_best=rho0;
                                                    rho_inf_posterior_best=rho_inf_posterior;
                                                    w_best=w;
                                                    E_best=E;
                                                end
                                            end
                                        end
                                    end
                                end
                                if rho_inf_posterior_best~=0
                                    AP_overactivity_ratio_PD(i_momocell)=rho_inf_anterior_best/rho_inf_posterior_best;
                                else
                                    if data_changdubi(number_changdubi)<=0.2
                                        AP_overactivity_ratio_PD(i_momocell)=max_AP_overactivity_ratio(1,k);
                                    elseif data_changdubi(number_changdubi)>0.2 && data_changdubi(number_changdubi)<=0.4
                                        AP_overactivity_ratio_PD(i_momocell)=max_AP_overactivity_ratio(2,k);
                                    elseif data_changdubi(number_changdubi)>0.4 && data_changdubi(number_changdubi)<=0.6
                                        AP_overactivity_ratio_PD(i_momocell)=max_AP_overactivity_ratio(3,k);
                                    elseif data_changdubi(number_changdubi)>0.6 && data_changdubi(number_changdubi)<=0.8
                                        AP_overactivity_ratio_PD(i_momocell)=max_AP_overactivity_ratio(4,k);
                                    else
                                        AP_overactivity_ratio_PD(i_momocell)=max_AP_overactivity_ratio(5,k);
                                    end
                                end
                            end
                            
                            % 统计两个极性方向的细胞数
                            forward_cellnumber=0;  % anterior lateral domain收缩环强的细胞数
                            backward_cellnumber=0;  % posterior lateral domain收缩环强的细胞数
                            for ii=1:number_monocell
                                if AP_overactivity_ratio_PD(ii)>1
                                    forward_cellnumber=forward_cellnumber+1;
                                else
                                    backward_cellnumber=backward_cellnumber+1;
                                end
                            end
                            
                            % 判断该张照片的A-P轴极性方向
                            if forward_cellnumber==number_monocell  % 所有无相邻细胞的细胞均为正向
                                polarity_direction=1;
                            elseif backward_cellnumber==number_monocell  % 所有无相邻细胞的细胞均为反向
                                polarity_direction=2;
                            elseif forward_cellnumber~=0 && forward_cellnumber~=number_monocell  % 无相邻细胞的细胞有两种方向
                                if forward_cellnumber>=backward_cellnumber
                                    polarity_direction=1;
                                else
                                    polarity_direction=2;
                                end
                            end
                        else
                            %%% 2. 查找该张照片中有无无相邻细胞的细胞数据
                            %%%    (2.2)若无，分别尝试两种极性方向，选择该张照片下所有细胞按该种极性方向处理后，所得大部分
                            %%%         细胞的实际极性方向相同且平均A-P强度比更大的极性方向
                            AP_overactivity_ratio_PD=zeros(number_multicell,2);  % 存放所有有相邻细胞的细胞(即所有细胞)的A-P极性强度比
                            
                            %%% 获得按两种极性方向拆分原始数据计算得到的所有细胞A-P轴强度比
                            AP_overactivity_ratio_test=zeros(2);  % 第一个数据为正向极性，第二个数据为反向极性
                            k_sum=0;
                            for i_sum=1:5  % 正向极性为5个时期(有实际数值)的平均极性大小
                                AP_overactivity_ratio_test(1)=AP_overactivity_ratio_test(1)+AP_overactivity_ratio(i_sum,k);
                                if AP_overactivity_ratio(i_sum,k)~=0
                                    k_sum=k_sum+1;
                                end
                            end
                            AP_overactivity_ratio_test(1)=AP_overactivity_ratio_test(1)/k_sum;
                            AP_overactivity_ratio_test(2)=1/AP_overactivity_ratio_test(1);     % 反向极性为正向极性的倒数
                            
                            for i_test=1:2
                                % 确定所有细胞数据文件的长度，保存在X_part1/2/3矩阵中，并将数据归一化
                                X_part1=zeros(number_multicell);
                                X_part2=zeros(number_multicell);
                                X_part3=zeros(number_multicell);
                                for i_multicell=1:number_multicell
                                    for i_X_part1=length(y_part1(:,i_multicell)):-1:1  % 确定y数据矩阵长度
                                        if y_part1(i_X_part1,i_multicell)~=0
                                            break
                                        end
                                    end
                                    X_part1(i_multicell)=i_X_part1;
                                    for i_X_part2=length(y_part2(:,i_multicell)):-1:1  % 确定y数据矩阵长度
                                        if y_part1(i_X_part2,i_multicell)~=0
                                            break
                                        end
                                    end
                                    X_part2(i_multicell)=i_X_part2;
                                    for i_X_part3=length(y_part3(:,i_multicell)):-1:1  % 确定y数据矩阵长度
                                        if y_part1(i_X_part3,i_multicell)~=0
                                            break
                                        end
                                    end
                                    X_part3(i_multicell)=i_X_part3;
                                end
                                    
                                % Y轴数据归一化(同除以y_part1/2/3的总平均值)
                                average_y=0;
                                for i_sum=1:X_part1(i_multicell)
                                    average_y=average_y+y_part1(i_sum,i_multicell);
                                end
                                for i_sum=1:X_part2(i_multicell)
                                    average_y=average_y+y_part2(i_sum,i_multicell);
                                end
                                for i_sum=1:X_part3(i_multicell)
                                    average_y=average_y+y_part3(i_sum,i_multicell);
                                end
                                average_y=average_y/(X_part1(i_multicell)+X_part2(i_multicell)+X_part3(i_multicell));
                                
                                y_part1=y_part1./average_y;
                                y_part2=y_part2./average_y;
                                y_part3=y_part3./average_y;
                                
                                % 细胞灰度值数据的预处理
                                for i_multicell=1:number_multicell
                                    % X轴关于坐标原点对称
                                    x_part1=linspace(-0.5,0.5,X_part1(i_multicell));
                                    x_part2=linspace(-0.5,0.5,X_part2(i_multicell));
                                    x_part3=linspace(-0.5,0.5,X_part3(i_multicell));
                                    
                                    Y_anterior=zeros(1,X_part1(i_multicell));  % 定义细胞两侧和中央的收缩环灰度值数据，用于单高斯分布拟合
                                    Y_middle=zeros(1,X_part2(i_multicell));
                                    Y_posterior=zeros(1,X_part3(i_multicell));
                                    
                                    % 查找上一个细胞是否存在(若上一个细胞的第三个文件和该细胞的第一个文件相同，则说明两细胞相邻，对y_part1进行处理)
                                    if i_multicell>=2  % 表示存在上一个细胞
                                        if X_part3(i_multicell-1)==X_part1(i_multicell) && ~sum(abs(y_part3(:,i_multicell-1)-y_part1(:,i_multicell)))
                                            total_data_previous2=0;
                                            total_data_2=0;
                                            for i_sum=1:X_part2(i_multicell-1)  % 计算上一个细胞basal domain的总灰度值
                                                total_data_previous2=total_data_previous2+y_part2(i_sum,i_multicell-1);
                                            end
                                            total_data_previous2=total_data_previous2/X_part2(i_multicell-1);  % 上一个细胞的basal domain的总灰度值取平均
                                            for i_sum=1:X_part2(i_multicell)  % 计算该细胞basal domain的总灰度值
                                                total_data_2=total_data_2+y_part2(i_sum,i_multicell);
                                            end
                                            total_data_2=total_data_2/X_part2(i_multicell);  % 该细胞的basal domain的总灰度值取平均
                                            for i_ratio=1:X_part1(i_multicell)  % 该细胞anterior(左)侧lateral domain灰度值数据的拆分 (公式：y_part1=y_part1*[ave(y_part2)*A-P-ratio/(ave(y_part2_before)+ave(y_part2)*A-P-ratio)]
                                                Y_anterior(i_ratio)=y_part1(i_ratio,i_multicell)*total_data_2*AP_overactivity_ratio_test(i_test)/(total_data_previous2+total_data_2*AP_overactivity_ratio_test(i_test));
                                            end
                                            
                                            % 将该细胞第一个数据文件的数据进行颠倒
                                            Y_anterior_fuben=Y_anterior;
                                            for i_symmetry=1:length(Y_anterior)
                                                Y_anterior(i_symmetry)=Y_anterior_fuben(length(Y_anterior)-i_symmetry+1);
                                            end
                                        else  % 不存在上一个细胞，直接将原始y_part1的值赋值给Y_anterior
                                            for i_ratio=1:X_part1(i_multicell)
                                                Y_anterior(i_ratio)=y_part1(i_ratio,i_multicell);
                                            end
                                        end
                                    else  % 该张照片第一个细胞，直接将原始y_part1的值赋值给Y_anterior
                                        for i_ratio=1:X_part1(i_multicell)
                                            Y_anterior(i_ratio)=y_part1(i_ratio,i_multicell);
                                        end
                                    end
                                    
                                    % 查找下一个细胞是否存在(若下一个细胞的第一个文件和该细胞的第三个文件相同，则说明两细胞相邻，对y_part3进行处理)
                                    if i_multicell<=number_multicell-1  % 表示存在下一个细胞
                                        if X_part3(i_multicell)==X_part1(i_multicell+1) && ~sum(abs(y_part3(:,i_multicell)-y_part1(:,i_multicell+1)))
                                            total_data_next2=0;
                                            total_data_2=0;
                                            for i_sum=1:X_part2(i_multicell+1)  % 计算下一个细胞basal domain的总灰度值
                                                total_data_next2=total_data_next2+y_part2(i_sum,i_multicell+1);
                                            end
                                            total_data_next2=total_data_next2/X_part2(i_multicell+1);  % 下一个细胞basal domain的总灰度值取平均
                                            for i_sum=1:X_part2(i_multicell)  % 计算该细胞basal domain的总灰度值
                                                total_data_2=total_data_2+y_part2(i_sum,i_multicell);
                                            end
                                            total_data_2=total_data_2/X_part2(i_multicell);  % 该细胞basal domain的总灰度值取平均
                                            for i_ratio=1:X_part1(i_multicell)  % 该细胞posterior(右)侧lateral domain灰度值数据的拆分 (公式：y_part3=y_part3*[ave(y_part2)/(ave(y_part2)+ave(y_part2_after)*A-P-ratio)]
                                                Y_posterior(i_ratio)=y_part3(i_ratio,i_multicell)*total_data_2/(total_data_2+total_data_next2*AP_overactivity_ratio_test(i_test));
                                            end
                                        else  % 不存在下一个细胞，直接将原始y_part3的值赋值给Y_posterior
                                            for i_ratio=1:X_part1(i_multicell)
                                                Y_posterior(i_ratio)=y_part3(i_ratio,i_multicell);
                                            end
                                        end
                                    else  % 该张照片最后一个细胞，直接将原始y_part3的值赋值给Y_posterior
                                        for i_ratio=1:X_part1(i_multicell)
                                            Y_posterior(i_ratio)=y_part3(i_ratio,i_multicell);
                                        end
                                    end
                                    
                                    % 最小二乘拟合两侧高斯分布，获得两侧收缩环的强度，计算A-P强度比
                                    % 拟合Y_anterior
                                    wucha_best=100000000;
                                    for rho0=max(Y_anterior)*0.001:0.999*max(Y_anterior)/25:max(Y_anterior)
                                        for rho_inf_anterior=0.01:(5-0.01)/25:5
                                            for w=0.01:(1-0.01)/25:1
                                                for E=-0.5:(0.5+0.5)/25:0.5
                                                    wucha=0;
                                                    for ii=1:X_part1(i_multicell)
                                                        wucha=wucha+(rho0+rho_inf_anterior*exp(-1/2*((x_part1(ii)-E)/w)^2)-Y_anterior(ii))^2;
                                                    end
                                                    if wucha<wucha_best
                                                        wucha_best=wucha;
                                                        rho0_best=rho0;
                                                        rho_inf_anterior_best=rho_inf_anterior;
                                                        w_best=w;
                                                        E_best=E;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    
                                    % 拟合Y_posterior
                                    wucha_best=100000000;
                                    for rho0=max(Y_posterior)*0.001:0.999*max(Y_posterior)/25:max(Y_posterior)
                                        for rho_inf_posterior=0.01:(5-0.01)/25:5
                                            for w=0.01:(1-0.01)/25:1
                                                for E=-0.5:(0.5+0.5)/25:0.5
                                                    wucha=0;
                                                    for ii=1:X_part3(i_multicell)
                                                        wucha=wucha+(rho0+rho_inf_posterior*exp(-1/2*((x_part3(ii)-E)/w)^2)-Y_posterior(ii))^2;
                                                    end
                                                    if wucha<wucha_best
                                                        wucha_best=wucha;
                                                        rho0_best=rho0;
                                                        rho_inf_posterior_best=rho_inf_posterior;
                                                        w_best=w;
                                                        E_best=E;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    if i_test==1  % A-P极性强度比为rho_inf_anterior_best/rho_inf_posterior_best
                                        if rho_inf_posterior_best~=0
                                            AP_overactivity_ratio_PD(i_multicell,i_test)=rho_inf_anterior_best/rho_inf_posterior_best;
                                        else
                                            if data_changdubi(number_changdubi)<=0.2
                                                AP_overactivity_ratio_PD(i_multicell,i_test)=max_AP_overactivity_ratio(1,k);
                                            elseif data_changdubi(number_changdubi)>0.2 && data_changdubi(number_changdubi)<=0.4
                                                AP_overactivity_ratio_PD(i_multicell,i_test)=max_AP_overactivity_ratio(2,k);
                                            elseif data_changdubi(number_changdubi)>0.4 && data_changdubi(number_changdubi)<=0.6
                                                AP_overactivity_ratio_PD(i_multicell,i_test)=max_AP_overactivity_ratio(3,k);
                                            elseif data_changdubi(number_changdubi)>0.6 && data_changdubi(number_changdubi)<=0.8
                                                AP_overactivity_ratio_PD(i_multicell,i_test)=max_AP_overactivity_ratio(4,k);
                                            else
                                                AP_overactivity_ratio_PD(i_multicell,i_test)=max_AP_overactivity_ratio(5,k);
                                            end
                                        end
                                    elseif i_test==2  % A-P极性强度比为rho_inf_posterior_best/rho_inf_anterior_best
                                        if rho_inf_anterior_best~=0
                                            AP_overactivity_ratio_PD(i_multicell,i_test)=rho_inf_posterior_best/rho_inf_anterior_best;
                                        else
                                            if data_changdubi(number_changdubi)<=0.2
                                                AP_overactivity_ratio_PD(i_multicell,i_test)=max_AP_overactivity_ratio(1,k);
                                            elseif data_changdubi(number_changdubi)>0.2 && data_changdubi(number_changdubi)<=0.4
                                                AP_overactivity_ratio_PD(i_multicell,i_test)=max_AP_overactivity_ratio(2,k);
                                            elseif data_changdubi(number_changdubi)>0.4 && data_changdubi(number_changdubi)<=0.6
                                                AP_overactivity_ratio_PD(i_multicell,i_test)=max_AP_overactivity_ratio(3,k);
                                            elseif data_changdubi(number_changdubi)>0.6 && data_changdubi(number_changdubi)<=0.8
                                                AP_overactivity_ratio_PD(i_multicell,i_test)=max_AP_overactivity_ratio(4,k);
                                            else
                                                AP_overactivity_ratio_PD(i_multicell,i_test)=max_AP_overactivity_ratio(5,k);
                                            end
                                        end
                                    end
                                end
                            end
                            
                            % 统计两个极性方向的细胞数
                            forward_cellnumber=0;  % 正向极性拆分后，和假设极性方向相同(anterior lateral domain收缩环较强)的细胞数
                            backward_cellnumber=0;  % 反向极性拆分后，和假设极性方向相同(posterior lateral domain收缩环较强)的细胞数
                            for ii=1:number_multicell
                                if AP_overactivity_ratio_PD(ii,1)>1
                                    forward_cellnumber=forward_cellnumber+1;
                                end
                                if AP_overactivity_ratio_PD(ii,2)>1
                                    backward_cellnumber=backward_cellnumber+1;
                                end
                            end
                            
                            %%% 判断该张照片的A-P轴极性方向
                            if forward_cellnumber==number_multicell  % 正向极性拆分后，和假设极性方向相同(anterior lateral domain收缩环较强)的细胞数等于总细胞数
                                polarity_direction=1;
                            elseif backward_cellnumber==number_multicell  % 反向极性拆分后，和假设极性方向相同(posterior lateral domain收缩环较强)的细胞数等于总细胞数
                                polarity_direction=2;
                            elseif forward_cellnumber~=number_multicell  % 细胞极性方向无统一规律
                                % 根据 极性与该极性拆分方向相同的细胞数+所有细胞A-P ratio的积 作为评判标准
                                forward_standard=AP_overactivity_ratio_PD(1,1);  % 正向极性评判标准参数计算值
                                backward_standard=AP_overactivity_ratio_PD(1,2);  % 反向极性评判标准参数计算值
                                for i_multicell=2:number_multicell  % 计算所有细胞A-P ratio的积
                                    forward_standard=forward_standard*AP_overactivity_ratio_PD(i_multicell,1);
                                    backward_standard=backward_standard*AP_overactivity_ratio_PD(i_multicell,2);
                                end
                                forward_standard=forward_standard+forward_cellnumber;  % 加上极性与该极性拆分方向相同的细胞数
                                backward_standard=backward_standard+backward_cellnumber;
                                if forward_standard>=backward_standard
                                    polarity_direction=1;
                                else
                                    polarity_direction=2;
                                end
                            end
                        end
                    end
                    
                    
                    %%% 细胞特征参数的拟合(数据文件预处理、拼接、多余数据的删减、数据归一化、三高斯分布拟合)
                    for j=1:10
                        filename2=num2str(j);
                        filename=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3];  % 如 0928_111
                        framepath_filename=[framepath2 filename];
                        
                        % 拼接完整的文件路径，用于下一步exist函数的查找
                        framepath_filename_1=[framepath_filename filename_xiahuaxian filename_1 filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111_1.csv
                        framepath_filename_2=[framepath_filename filename_xiahuaxian filename_2 filename_csv];
                        framepath_filename_3=[framepath_filename filename_xiahuaxian filename_3 filename_csv];
                        framepath_filename_danyi=[framepath_filename filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111.csv
                        
                        if exist(framepath_filename_1,'file')   % 一个细胞的数据有多个文件
                            number_changdubi=number_changdubi+1;
                            for aa=1
                                %%% 确定该细胞的A-P轴强度比
                                A_P_ratio2=0;  % 定义该细胞的A-P轴极性强度比
                                if polarity_direction==1
                                    if data_changdubi(number_changdubi)<=0.2
                                        A_P_ratio2=AP_overactivity_ratio(1,k);
                                    elseif data_changdubi(number_changdubi)>0.2 && data_changdubi(number_changdubi)<=0.4
                                        A_P_ratio2=AP_overactivity_ratio(2,k);
                                    elseif data_changdubi(number_changdubi)>0.4 && data_changdubi(number_changdubi)<=0.6
                                        A_P_ratio2=AP_overactivity_ratio(3,k);
                                    elseif data_changdubi(number_changdubi)>0.6 && data_changdubi(number_changdubi)<=0.8
                                        A_P_ratio2=AP_overactivity_ratio(4,k);
                                    else
                                        A_P_ratio2=AP_overactivity_ratio(5,k);
                                    end
                                elseif polarity_direction==2
                                    if data_changdubi(number_changdubi)<=0.2
                                        A_P_ratio2=1/AP_overactivity_ratio(1,k);
                                    elseif data_changdubi(number_changdubi)>0.2 && data_changdubi(number_changdubi)<=0.4
                                        A_P_ratio2=1/AP_overactivity_ratio(2,k);
                                    elseif data_changdubi(number_changdubi)>0.4 && data_changdubi(number_changdubi)<=0.6
                                        A_P_ratio2=1/AP_overactivity_ratio(3,k);
                                    elseif data_changdubi(number_changdubi)>0.6 && data_changdubi(number_changdubi)<=0.8
                                        A_P_ratio2=1/AP_overactivity_ratio(4,k);
                                    else
                                        A_P_ratio2=1/AP_overactivity_ratio(5,k);
                                    end
                                end
                                
                                %%% 每个细胞的多个数据文件合并至一个矩阵中
                                data_1=xlsread(framepath_filename_1); % 导入数据
                                data_2=xlsread(framepath_filename_2);
                                data_3=xlsread(framepath_filename_3);
                                
                                % 计算细胞的anterior-lateral basal posterior-lateral domain的长度
                                anterior_lateral_length=data_1(end,1);
                                posterior_lateral_length=data_3(end,1);
                                basal_length=data_2(end,1);
                                
                                %%% 细胞原始数据拆分
                                protein_number=1;
                                [data_1,data_2,data_3] = CellOverlapDataSplitting(data_1,data_2,data_3,C_AP,protein_number,data_changdubi,number_changdubi,j);
                                data_split=[data_3_split1;data_3_split2];

                                % 将处理完成的一个细胞的三组数据合并
                                length_data_1=length(data_1);
                                length_data_2=length(data_2);
                                length_data_3=length(data_3);
                                data=zeros(length_data_1+length_data_2+length_data_3,2);
                                for ii=1:length_data_1
                                    data(ii,1)=ii;
                                    data(ii,2)=data_1(ii,2);
                                end
                                for ii=length_data_1+1:length_data_1+length_data_2
                                    data(ii,1)=ii;
                                    data(ii,2)=data_2(ii-length_data_1,2);
                                end
                                for ii=length_data_1+length_data_2+1:length_data_1+length_data_2+length_data_3
                                    data(ii,1)=ii;
                                    data(ii,2)=data_3(ii-length_data_1-length_data_2,2);
                                end

                                %%% 对原始细胞收缩环灰度值数据进行两侧多余数据的删减处理
                                % 根据测量数据拟合basal-lateral domain三高斯分布
                                % 数据点拉伸
                                [X,Y]=size(data); % 丈量数据矩阵尺寸
                                if mod(X,2)==0   % 若数据为偶数个，则进行线性插值，方便后续的平移工作中的对称性
                                    data2=zeros(X+1,Y);
                                    data2(:,1)=imresize(data(:,1),[X+1,1], 'bilinear'); % 将X轴进行线性插值
                                    data2(:,2)=imresize(data(:,2),[X+1,1], 'bilinear'); % 将Y轴进行线性插值
                                else
                                    data2=data;
                                end

                                % X轴平移至关于原点对称的区间，并归一化
                                bias=0.5; % 平移量
                                x=data2(:,1)/max(data2(:,1)); % 对x轴进行归一化
                                x=x-bias; % x轴平移至关于坐标原点对称
                                % Y轴数据归一化
                                y=data2(:,2)/max(data2(:,2));


                                Y_anterior=zeros(1,round((X-1)*0.2));  % 定义细胞两侧的收缩环灰度值数据，用于单高斯分布拟合
                                Y_posterior=zeros(1,round((X-1)*0.2));
                                Y_anterior_AP=zeros(1,round((X-1)*0.2));  % 定义细胞另一蛋白两侧的收缩环灰度值数据，用于单高斯分布拟合
                                Y_posterior_AP=zeros(1,round((X-1)*0.2));

                                if data_changdubi(number_changdubi)<=0.2  % lumen形成初期两侧收缩环极性未建立，不删除数据
                                    for k_AP=1:10
                                        if k_AP~=k
                                            filename3_AP=num2str(k_AP);  % AP: anothor protein
                                            filename_AP=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3_AP];  % 如 0928_111
                                            framepath_filename_AP=[framepath2 filename_AP];
                                            framepath_filename_1_AP=[framepath_filename_AP filename_xiahuaxian filename_1 filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111_1.csv
                                            framepath_filename_2_AP=[framepath_filename_AP filename_xiahuaxian filename_2 filename_csv];
                                            framepath_filename_3_AP=[framepath_filename_AP filename_xiahuaxian filename_3 filename_csv];
                                            if exist(framepath_filename_1_AP,'file')
                                                break
                                            end
                                        end
                                    end
                                    if k_AP~=10  % 表示该细胞有另一个蛋白数据
                                        %%% 每个细胞的多个数据文件合并至一个矩阵中
                                        data_1_AP=xlsread(framepath_filename_1_AP); % 导入数据
                                        data_2_AP=xlsread(framepath_filename_2_AP);
                                        data_3_AP=xlsread(framepath_filename_3_AP);

                                        %%% 细胞原始数据拆分
                                        protein_number=2;
                                        [data_1_AP,data_2_AP,data_3_AP] = CellOverlapDataSplitting(data_1_AP,data_2_AP,data_3_AP,C_AP,protein_number,data_changdubi,number_changdubi,j);
                                        data_AP_split=[data_3_AP_split1;data_3_AP_split2];

                                        % 将处理完成的一个细胞的三组数据合并
                                        length_data_1=length(data_1_AP);
                                        length_data_2=length(data_2_AP);
                                        length_data_3=length(data_3_AP);
                                        data_AP=zeros(length_data_1+length_data_2+length_data_3,2);
                                        for ii=1:length_data_1
                                            data_AP(ii,1)=ii;
                                            data_AP(ii,2)=data_1_AP(ii,2);
                                        end
                                        for ii=length_data_1+1:length_data_1+length_data_2
                                            data_AP(ii,1)=ii;
                                            data_AP(ii,2)=data_2_AP(ii-length_data_1,2);
                                        end
                                        for ii=length_data_1+length_data_2+1:length_data_1+length_data_2+length_data_3
                                            data_AP(ii,1)=ii;
                                            data_AP(ii,2)=data_3_AP(ii-length_data_1-length_data_2,2);
                                        end

                                        % 数据点拉伸
                                        [X_AP,Y_AP]=size(data_AP); % 丈量数据矩阵尺寸
                                        if mod(X_AP,2)==0   % 若数据为偶数个，则进行线性插值，方便后续的平移工作中的对称性
                                            data2_AP=zeros(X_AP+1,Y_AP);
                                            data2_AP(:,1)=imresize(data_AP(:,1),[X_AP+1,1], 'bilinear'); % 将X轴进行线性插值
                                            data2_AP(:,2)=imresize(data_AP(:,2),[X_AP+1,1], 'bilinear'); % 将Y轴进行线性插值
                                        else
                                            data2_AP=data_AP;
                                        end

                                        % X轴平移至关于原点对称的区间，并归一化
                                        bias=0.5; % 平移量
                                        x_AP=data2_AP(:,1)/max(data2_AP(:,1)); % 对x轴进行归一化
                                        x_AP=x_AP-bias; % x轴平移至关于坐标原点对称
                                        x=x_AP;
                                        % Y轴数据归一化
                                        y=data2(:,2)/max(data2(:,2));
                                        y_AP=data2_AP(:,2)/max(data2_AP(:,2));
                                    end
                                    x_min=-0.5;x_max=0.5;
                                else  % lumen形成中后期两侧收缩环极性已建立，需删除数据
                                    for k_AP=1:10
                                        if k_AP~=k
                                            filename3_AP=num2str(k_AP);  % AP: anothor protein
                                            filename_AP=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3_AP];  % 如 0928_111
                                            framepath_filename_AP=[framepath2 filename_AP];
                                            framepath_filename_1_AP=[framepath_filename_AP filename_xiahuaxian filename_1 filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111_1.csv
                                            framepath_filename_2_AP=[framepath_filename_AP filename_xiahuaxian filename_2 filename_csv];
                                            framepath_filename_3_AP=[framepath_filename_AP filename_xiahuaxian filename_3 filename_csv];
                                            if exist(framepath_filename_1_AP,'file')
                                                break
                                            end
                                        end
                                    end
                                    if k_AP~=10  % 表示该细胞有另一个蛋白数据
                                        %%% 每个细胞的多个数据文件合并至一个矩阵中
                                        data_1_AP=xlsread(framepath_filename_1_AP); % 导入数据
                                        data_2_AP=xlsread(framepath_filename_2_AP);
                                        data_3_AP=xlsread(framepath_filename_3_AP);

                                        %%% 细胞原始数据拆分
                                        protein_number=2;
                                        [data_1_AP,data_2_AP,data_3_AP] = CellOverlapDataSplitting(data_1_AP,data_2_AP,data_3_AP,C_AP,protein_number,data_changdubi,number_changdubi,j);
                                        data_AP_split=[data_3_AP_split1;data_3_AP_split2];

                                        % 将处理完成的一个细胞的三组数据合并
                                        length_data_1=length(data_1_AP);
                                        length_data_2=length(data_2_AP);
                                        length_data_3=length(data_3_AP);
                                        data_AP=zeros(length_data_1+length_data_2+length_data_3,2);
                                        for ii=1:length_data_1
                                            data_AP(ii,1)=ii;
                                            data_AP(ii,2)=data_1_AP(ii,2);
                                        end
                                        for ii=length_data_1+1:length_data_1+length_data_2
                                            data_AP(ii,1)=ii;
                                            data_AP(ii,2)=data_2_AP(ii-length_data_1,2);
                                        end
                                        for ii=length_data_1+length_data_2+1:length_data_1+length_data_2+length_data_3
                                            data_AP(ii,1)=ii;
                                            data_AP(ii,2)=data_3_AP(ii-length_data_1-length_data_2,2);
                                        end
                                        % 数据点拉伸
                                        [X_AP,Y_AP]=size(data_AP); % 丈量数据矩阵尺寸
                                        if mod(X_AP,2)==0   % 若数据为偶数个，则进行线性插值，方便后续的平移工作中的对称性
                                            data2_AP=zeros(X_AP+1,Y_AP);
                                            data2_AP(:,1)=imresize(data_AP(:,1),[X_AP+1,1], 'bilinear'); % 将X轴进行线性插值
                                            data2_AP(:,2)=imresize(data_AP(:,2),[X_AP+1,1], 'bilinear'); % 将Y轴进行线性插值
                                        else
                                            data2_AP=data_AP;
                                        end

                                        % X轴平移至关于原点对称的区间，并归一化
                                        bias=0.5; % 平移量
                                        x_AP=data2_AP(:,1)/max(data2_AP(:,1)); % 对x轴进行归一化
                                        x_AP=x_AP-bias; % x轴平移至关于坐标原点对称
                                        x=x_AP;
                                        % Y轴数据归一化
                                        y=data2(:,2)/max(data2(:,2));
                                        y_AP=data2_AP(:,2)/max(data2_AP(:,2));

                                        % Levenberg-Marquardt method拟合两侧高斯分布，确定要删去的部分
                                        % 该蛋白数据拟合
                                        for ab=1
                                            if data_changdubi(number_changdubi)<=0.2  % lumen形成早期，lateral domain收缩环未建立，共有9个特征参数拟合
                                                E_lateral_start=0.3;
                                                E_lateral_end=0.5;
                                            else  % lumen形成中后期，lateral domain收缩环建立，共有8个特征参数拟合
                                                E_lateral_start=0.5;
                                                E_lateral_end=0.5;
                                            end

                                            wucha_best=100000000;
                                            % 使用Levenberg-Marquardt method拟合
                                            for i_fit=1:5  % 循环10次，寻找最优拟合解
                                                ft=fittype('rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)','independent','x','dependent','y');
                                                opts=fitoptions('Method','NonlinearLeastSquares');
                                                opts.Algorithm='Levenberg-Marquardt';
                                                opts.Display='Off';
                                                opts.MaxFunEvals=600;
                                                opts.MaxIter=400;
                                                opts.Lower=[-0.2 -0.5 0.45 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
                                                opts.Upper=[0.2 -0.45 0.5 1 5 5 5 1 1 1];
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
                                                E_anterior=fitresult.E_anterior;
                                                E_posterior=fitresult.E_posterior;

                                                wucha=0;
                                                for ii=1:length(x)
                                                    wucha=wucha+(rho0+rho_inf*exp(-1/2*((x(ii)-E)/w)^2)+rho_inf_anterior*exp(-1/2*((x(ii)-E_anterior)/w_anterior)^2)+...
                                                        rho_inf_posterior*exp(-1/2*((x(ii)-E_posterior)/w_posterior)^2)-y(ii))^2;
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
                                                    E_anterior_best=E_anterior;
                                                    E_posterior_best=E_posterior;
                                                end
                                            end
                                        end
                                        % 另一个蛋白数据拟合
                                        for ab=1
                                            if data_changdubi(number_changdubi)<=0.2  % lumen形成早期，lateral domain收缩环未建立，共有9个特征参数拟合
                                                E_lateral_start=0.3;
                                                E_lateral_end=0.5;
                                            else  % lumen形成中后期，lateral domain收缩环建立，共有8个特征参数拟合
                                                E_lateral_start=0.5;
                                                E_lateral_end=0.5;
                                            end

                                            wucha_best=100000000;
                                            % 使用Levenberg-Marquardt method拟合
                                            for i_fit=1:5  % 循环10次，寻找最优拟合解
                                                ft=fittype('rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)','independent','x','dependent','y');
                                                opts=fitoptions('Method','NonlinearLeastSquares');
                                                opts.Algorithm='Levenberg-Marquardt';
                                                opts.Display='Off';
                                                opts.MaxFunEvals=600;
                                                opts.MaxIter=400;
                                                opts.Lower=[-0.2 -0.5 0.45 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
                                                opts.Upper=[0.2 -0.45 0.5 1 5 5 5 1 1 1];
                                                [xData,yData]=prepareCurveData(x_AP,y_AP);
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
                                                E_anterior=fitresult.E_anterior;
                                                E_posterior=fitresult.E_posterior;

                                                wucha=0;
                                                for ii=1:length(x_AP)
                                                    wucha=wucha+(rho0+rho_inf*exp(-1/2*((x_AP(ii)-E)/w)^2)+rho_inf_anterior*exp(-1/2*((x_AP(ii)-E_anterior)/w_anterior)^2)+...
                                                        rho_inf_posterior*exp(-1/2*((x_AP(ii)-E_posterior)/w_posterior)^2)-y_AP(ii))^2;
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
                                                    E_anterior_AP_best=E_anterior;
                                                    E_posterior_AP_best=E_posterior;
                                                end
                                            end
                                        end
                                        x_min=min(E_anterior_best,E_anterior_AP_best);  % 前侧边界取两个蛋白的拟合前边界中较小的
                                        x_max=max(E_posterior_best,E_posterior_AP_best);  % 后侧边界取两个蛋白的拟合后边界中较大的
                                    else  % kk=9表示该细胞仅有一种蛋白的数据文件
                                        % 读取该蛋白的数据
                                        % 数据点拉伸

                                        % Levenberg-Marquardt method拟合两侧高斯分布，确定要删去的部分
                                        for ab=1
                                            if data_changdubi(number_changdubi)<=0.2  % lumen形成早期，lateral domain收缩环未建立，共有9个特征参数拟合
                                                E_lateral_start=0.3;
                                                E_lateral_end=0.5;
                                            else  % lumen形成中后期，lateral domain收缩环建立，共有8个特征参数拟合
                                                E_lateral_start=0.5;
                                                E_lateral_end=0.5;
                                            end

                                            wucha_best=100000000;
                                            % 使用Levenberg-Marquardt method拟合
                                            for i_fit=1:5  % 循环10次，寻找最优拟合解
                                                ft=fittype('rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)','independent','x','dependent','y');
                                                opts=fitoptions('Method','NonlinearLeastSquares');
                                                opts.Algorithm='Levenberg-Marquardt';
                                                opts.Display='Off';
                                                opts.MaxFunEvals=600;
                                                opts.MaxIter=400;
                                                opts.Lower=[-0.2 -0.5 0.45 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
                                                opts.Upper=[0.2 -0.45 0.5 1 5 5 5 1 1 1];
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
                                                E_anterior=fitresult.E_anterior;
                                                E_posterior=fitresult.E_posterior;

                                                wucha=0;
                                                for ii=1:length(x)
                                                    wucha=wucha+(rho0+rho_inf*exp(-1/2*((x(ii)-E)/w)^2)+rho_inf_anterior*exp(-1/2*((x(ii)-E_anterior)/w_anterior)^2)+...
                                                        rho_inf_posterior*exp(-1/2*((x(ii)-E_posterior)/w_posterior)^2)-y(ii))^2;
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
                                                    E_anterior_best=E_anterior;
                                                    E_posterior_best=E_posterior;
                                                end
                                            end
                                        end
                                        x_min=E_anterior_best;
                                        x_max=E_posterior_best;
                                    end
                                end

                                % 若x_min或x_max偏离边界过多(超过总细胞长度的1/10)，说明数据拟合效果不好，仍取原始边界值
                                if (x_min+0.5)/0.5>0.1
                                    x_min=-0.5;
                                elseif (0.5-x_max)/0.5>0.1
                                    x_max=0.5;
                                end
                                
                                %%% 对删减处理完成的细胞灰度值数据进行三高斯分布拟合
                                % 找到原X轴上和x_min x_max最接近的点
                                wucha_best=100;
                                for ii=1:X
                                    wucha=abs(x(ii)-x_min);
                                    if wucha<wucha_best
                                        wucha_best=wucha;
                                        x_min_best=x(ii);  % 将最接近x_min的原X坐标值暂时保存在x_min_best变量中
                                    end
                                end
                                x_min=x_min_best;  % x_min改成X轴上最接近x_min点的值
                                wucha_best=100;
                                for ii=1:X
                                    wucha=abs(x(ii)-x_max);
                                    if wucha<wucha_best
                                        wucha_best=wucha;
                                        x_max_best=x(ii);  % 将最接近x_max的原X坐标值暂时保存在x_max_best变量中
                                    end
                                end
                                x_max=x_max_best;  % x_max改成X轴上最接近x_max点的值
                                
                                % 删去X轴和对应灰度值(y)上小于x_min和大于x_max的部分
                                x_delete=zeros(round((X-1)*(x_max-x_min)+1),1);
                                y_delete=zeros(round((X-1)*(x_max-x_min)+1),1);
                                for ii=1:round((X-1)*(x_max-x_min)+1)
                                    x_delete(ii)=x(ii+round((x_min+0.5)*(X-1)));
                                    y_delete(ii)=y(ii+round((x_min+0.5)*(X-1)));
                                end
                                x=x_delete;
                                y=y_delete;
                                % X轴重新平移至关于原点对称的区间，并拉伸至(-0.5,0.5)
                                bias=x_max+x_min; % 平移量
                                stretch=1/(x_max-x_min); % 拉伸量
                                x=x-bias/2; % x轴平移至关于坐标原点对称
                                x=x*stretch;
                                
                                %%% 使用Levenberg-Marquardt method拟合
                                max_y=0;min_y=100;
                                for iy=floor(length(y)*0.25):floor(length(y)*0.75)
                                    if y(iy)>max_y
                                        max_y=y(iy);
                                    elseif y(iy)<min_y
                                        min_y=y(iy);
                                    end
                                end
                                if max_y/min_y<1.5  % 表示赤道轴收缩环已经几乎消失，lumen时期处于后期，为防止拟合出现误差，将收缩环位置固定在中央
                                    E_bottom=-0.05;E_top=0.05;
                                else  % 收缩环较强，其位置不固定
                                    E_bottom=-0.1;E_top=0.1;
                                end

                                % 多高斯拟合
                                % 调用MultiGaussian_fit子函数计算Multi-Gaussian distribution
                                cell_type=2;
                                [rho0_best,rho_inf_array,E_array,w_array,break_point,wucha_best,TriGaussian_parameter]=MultiGaussian_fit(x,y,data_changdubi(number_changdubi),anterior_lateral_length,basal_length,posterior_lateral_length,cell_type);
                                framepath_filename
                                % 调用MultiGaussianTransfertoEquivalentTriGaussian子函数计算等效Tri-Gaussian distribution
                                [EquivalentTriGaussian_parameter,area_sumsum]=MultiGaussianTransfertoEquivalentTriGaussian(x,y,rho0_best,rho_inf_array,E_array,w_array,break_point);
                            end
                            
                            % 判断细胞是否是左侧为anterior part，右侧为posterior part。若不是，将拟合参数、拟合曲线、原始数据关于原点做对称处理
                            if rho_inf_array(end-1)<rho_inf_array(end)
                                % MultiGaussian参数做对称处理
                                rho_inf_array_reverse=zeros(size(rho_inf_array));
                                E_array_reverse=zeros(size(rho_inf_array));
                                w_array_reverse=zeros(size(rho_inf_array));
                                for i_Gaussian=1:length(rho_inf_array)-2
                                    rho_inf_array_reverse(i_Gaussian)=rho_inf_array(end-i_Gaussian-1);
                                    E_array_reverse(i_Gaussian)=E_array(end-i_Gaussian-1);
                                    w_array_reverse(i_Gaussian)=w_array(end-i_Gaussian-1);
                                end
                                rho_inf_array_reverse(end-1)=rho_inf_array(end);rho_inf_array_reverse(end)=rho_inf_array(end-1);
                                E_array_reverse(end-1)=E_array(end);E_array_reverse(end)=E_array(end-1);
                                w_array_reverse(end-1)=w_array(end);w_array_reverse(end)=w_array(end-1);
                                
                                rho_inf_array=rho_inf_array_reverse;
                                E_array=-E_array_reverse;
                                w_array=w_array_reverse;

                                % EquivalentTriGaussian参数做对称处理
                                anterior_posterior_ring_reverse=EquivalentTriGaussian_parameter(5:6);
                                EquivalentTriGaussian_parameter(5:6)=EquivalentTriGaussian_parameter(8:9);
                                EquivalentTriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                EquivalentTriGaussian_parameter(4)=-EquivalentTriGaussian_parameter(4);
                                
                                % TriGaussian参数做对称处理
                                anterior_posterior_ring_reverse=TriGaussian_parameter(5:6);
                                TriGaussian_parameter(5:6)=TriGaussian_parameter(8:9);
                                TriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                TriGaussian_parameter(4)=-TriGaussian_parameter(4);

                                y_reverse=zeros(length(y),1);  % 原始数据翻转
                                for ii=1:length(y)
                                    y_reverse(length(y)-ii+1)=y(ii);
                                end
                                y=y_reverse;
                            end

                            % 荧光强度值(纵坐标)归一化
                            % 将拟合完成的多高斯分布在定义域上的数学期望变换为1，同时将原始测量数据和等效三高斯分布按等比例拉伸
                            if length(rho_inf_array)-2==1
                                f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                    +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2);
                            elseif length(rho_inf_array)-2==2
                                f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                    +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2);
                            elseif length(rho_inf_array)-2==3
                                f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                    +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2)...
                                    +rho_inf_array(5)*exp(-1/2*((y-E_array(5))/w_array(5)).^2);
                            elseif length(rho_inf_array)-2==4
                                f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                    +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2)+...
                                    rho_inf_array(5)*exp(-1/2*((y-E_array(5))/w_array(5)).^2)+rho_inf_array(6)*exp(-1/2*((y-E_array(6))/w_array(6)).^2);
                            elseif length(rho_inf_array)-2>=5
                                f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                    +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2)...
                                    +rho_inf_array(5)*exp(-1/2*((y-E_array(5))/w_array(5)).^2)+rho_inf_array(6)*exp(-1/2*((y-E_array(6))/w_array(6)).^2)...
                                    +rho_inf_array(7)*exp(-1/2*((y-E_array(7))/w_array(7)).^2);
                            end

                            a0=integral(f,-1/2,1/2);
                            rho0_best=rho0_best/a0;
                            rho_inf_array=rho_inf_array./a0;
                            EquivalentTriGaussian_parameter(2:3:8)=EquivalentTriGaussian_parameter(2:3:8)./a0;EquivalentTriGaussian_parameter(1)=EquivalentTriGaussian_parameter(1)./a0;
                            TriGaussian_parameter(2:3:8)=TriGaussian_parameter(2:3:8)./a0;TriGaussian_parameter(1)=TriGaussian_parameter(1)./a0;
                            y_MultiGaussian_fit=ones(size(x))*rho0_best;  % 多高斯拟合曲线
                            for i_Gaussian=1:length(rho_inf_array)
                                y_MultiGaussian_fit=y_MultiGaussian_fit+rho_inf_array(i_Gaussian)*exp(-1/2*((x-E_array(i_Gaussian))/w_array(i_Gaussian)).^2);
                            end
                            y_EquivalentTriGaussian_fit=ones(size(x))*EquivalentTriGaussian_parameter(1);  % 等效三高斯拟合曲线
                            for i_Gaussian=1:3
                                y_EquivalentTriGaussian_fit=y_EquivalentTriGaussian_fit+EquivalentTriGaussian_parameter(3*i_Gaussian-1)*exp(-1/2*((x-EquivalentTriGaussian_parameter(3*i_Gaussian+1))/EquivalentTriGaussian_parameter(3*i_Gaussian)).^2);
                            end
                            y=y./a0;

                            %%% 绘图
                            filename_plot=[filename_photo filename_kongge filename1 filename_kongge filename_cell filename_kongge filename2 filename_kongge filename_protein];

                            figure
                            plot(x,y,'.','MarkerSize',15,'Color','k');  % 绘制数据点
                            hold on;
                            plot(x,y_MultiGaussian_fit,'b','linewidth',3);  % 绘制拟合的高斯曲线
                            plot(x,y_EquivalentTriGaussian_fit,'r','linewidth',3);  % 绘制拟合的高斯曲线
                            hold off;
                            axis ([-0.5 0.5 0 max(y)+0.1]);
                            title(filename_plot)
                            xlabel('x')
                            ylabel('cortex thickness')
                            set(gca,'YTick',0:0.2:max(y)+0.1);  % 修改y坐标刻度
                            set(gca,'XTick',-0.5:0.1:0.5);  % 修改y坐标刻度
                            legend('Experiments','Multi-Gaussian Curve','Tri-Gaussian Curve')
                            set(gca,'FontName','Arial','FontSize',12)
                            saveas(gcf,[framepath2,filename_plot,'.jpg']);
                            close

                            %%% 拟合参数汇总
                            for aa=1
                                if data_changdubi(number_changdubi)<=0.2
                                    data_changdubi_all(k_position1(k),1,k)=data_changdubi(number_changdubi);
                                    filenamearray_all{k_position1(k),1,k}=filename;
                                    MultiGaussian_data_all(k_position1(k),1,1,k)=rho0_best;
                                    MultiGaussian_data_all(k_position1(k),2:length(rho_inf_array)+1,1,k)=rho_inf_array;
                                    MultiGaussian_data_all(k_position1(k),length(rho_inf_array)+2:2*length(rho_inf_array)+1,1,k)=w_array;
                                    MultiGaussian_data_all(k_position1(k),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,1,k)=E_array;
                                    EquivalentTriGaussian_data_all(k_position1(k),:,1,k)=EquivalentTriGaussian_parameter;
                                    TriGaussian_data_all(k_position1(k),:,1,k)=TriGaussian_parameter;
                                    wucha_best_all(k_position1(k),1,k)=wucha_best;
                                    k_position1(k)=k_position1(k)+1;
                                elseif data_changdubi(number_changdubi)>0.2 && data_changdubi(number_changdubi)<=0.4
                                    data_changdubi_all(k_position2(k),2,k)=data_changdubi(number_changdubi);
                                    filenamearray_all{k_position2(k),2,k}=filename;
                                    MultiGaussian_data_all(k_position2(k),1,2,k)=rho0_best;
                                    MultiGaussian_data_all(k_position2(k),2:length(rho_inf_array)+1,2,k)=rho_inf_array;
                                    MultiGaussian_data_all(k_position2(k),length(rho_inf_array)+2:2*length(rho_inf_array)+1,2,k)=w_array;
                                    MultiGaussian_data_all(k_position2(k),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,2,k)=E_array;
                                    EquivalentTriGaussian_data_all(k_position2(k),:,2,k)=EquivalentTriGaussian_parameter;
                                    TriGaussian_data_all(k_position2(k),:,2,k)=TriGaussian_parameter;
                                    wucha_best_all(k_position2(k),2,k)=wucha_best;
                                    k_position2(k)=k_position2(k)+1;
                                elseif data_changdubi(number_changdubi)>0.4 && data_changdubi(number_changdubi)<=0.6
                                    data_changdubi_all(k_position3(k),3,k)=data_changdubi(number_changdubi);
                                    filenamearray_all{k_position3(k),3,k}=filename;
                                    MultiGaussian_data_all(k_position3(k),1,3,k)=rho0_best;
                                    MultiGaussian_data_all(k_position3(k),2:length(rho_inf_array)+1,3,k)=rho_inf_array;
                                    MultiGaussian_data_all(k_position3(k),length(rho_inf_array)+2:2*length(rho_inf_array)+1,3,k)=w_array;
                                    MultiGaussian_data_all(k_position3(k),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,3,k)=E_array;
                                    EquivalentTriGaussian_data_all(k_position3(k),:,3,k)=EquivalentTriGaussian_parameter;
                                    TriGaussian_data_all(k_position3(k),:,3,k)=TriGaussian_parameter;
                                    wucha_best_all(k_position3(k),3,k)=wucha_best;
                                    k_position3(k)=k_position3(k)+1;
                                elseif data_changdubi(number_changdubi)>0.6 && data_changdubi(number_changdubi)<=0.8
                                    data_changdubi_all(k_position4(k),4,k)=data_changdubi(number_changdubi);
                                    filenamearray_all{k_position4(k),4,k}=filename;
                                    MultiGaussian_data_all(k_position4(k),1,4,k)=rho0_best;
                                    MultiGaussian_data_all(k_position4(k),2:length(rho_inf_array)+1,4,k)=rho_inf_array;
                                    MultiGaussian_data_all(k_position4(k),length(rho_inf_array)+2:2*length(rho_inf_array)+1,4,k)=w_array;
                                    MultiGaussian_data_all(k_position4(k),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,4,k)=E_array;
                                    EquivalentTriGaussian_data_all(k_position4(k),:,4,k)=EquivalentTriGaussian_parameter;
                                    TriGaussian_data_all(k_position4(k),:,4,k)=TriGaussian_parameter;
                                    wucha_best_all(k_position4(k),4,k)=wucha_best;
                                    k_position4(k)=k_position4(k)+1;
                                else
                                    data_changdubi_all(k_position5(k),5,k)=data_changdubi(number_changdubi);
                                    filenamearray_all{k_position5(k),5,k}=filename;
                                    MultiGaussian_data_all(k_position5(k),1,5,k)=rho0_best;
                                    MultiGaussian_data_all(k_position5(k),2:length(rho_inf_array)+1,5,k)=rho_inf_array;
                                    MultiGaussian_data_all(k_position5(k),length(rho_inf_array)+2:2*length(rho_inf_array)+1,5,k)=w_array;
                                    MultiGaussian_data_all(k_position5(k),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,5,k)=E_array;
                                    EquivalentTriGaussian_data_all(k_position5(k),:,5,k)=EquivalentTriGaussian_parameter;
                                    TriGaussian_data_all(k_position5(k),:,5,k)=TriGaussian_parameter;
                                    wucha_best_all(k_position5(k),5,k)=wucha_best;
                                    k_position5(k)=k_position5(k)+1;
                                end
                            end
                            
                            %%% 拟合该细胞的另一个蛋白数据文件
                            if k_AP~=10  % 该细胞存在另一个蛋白数据文件
                                if k_AP==1
                                    filename_protein='lifeact';
                                elseif k_AP==2
                                    filename_protein='actin';
                                elseif k_AP==3
                                    filename_protein='MLC';
                                elseif k_AP==4
                                    filename_protein='RhoA';
                                elseif k_AP==5
                                    filename_protein='Anillin RBD';
                                elseif k_AP==6
                                    filename_protein='VD1';
                                elseif k_AP==7
                                    filename_protein='talinA';
                                elseif k_AP==8
                                    filename_protein='Cdc42WT';
                                elseif k_AP==9
                                    filename_protein='Cdc42D118A';
                                end

                                %%% 对删减处理完成的细胞灰度值数据进行三高斯分布拟合
                                % 删去X轴和对应灰度值(y)上小于x_min和大于x_max的部分
                                x_delete_AP=zeros(round((X-1)*(x_max-x_min)+1),1);
                                y_delete_AP=zeros(round((X-1)*(x_max-x_min)+1),1);
                                for ii=1:round((X-1)*(x_max-x_min)+1)
                                    x_delete_AP(ii)=x_AP(ii+round((x_min+0.5)*(X-1)));
                                    y_delete_AP(ii)=y_AP(ii+round((x_min+0.5)*(X-1)));
                                end
                                x_AP=x_delete_AP;
                                y_AP=y_delete_AP;
                                % X轴重新平移至关于原点对称的区间，并拉伸至(-0.5,0.5)
                                bias=x_max+x_min; % 平移量
                                stretch=1/(x_max-x_min); % 拉伸量
                                x_AP=x_AP-bias/2; % x轴平移至关于坐标原点对称
                                x_AP=x_AP*stretch;

                                %%% 使用Levenberg-Marquardt method拟合
                                max_y=0;min_y=100;
                                for iy=floor(length(y_AP)*0.25):floor(length(y_AP)*0.75)
                                    if y_AP(iy)>max_y
                                        max_y=y_AP(iy);
                                    elseif y_AP(iy)<min_y
                                        min_y=y_AP(iy);
                                    end
                                end
                                if max_y/min_y<1.5  % 表示赤道轴收缩环已经几乎消失，lumen时期处于后期，为防止拟合出现误差，将收缩环位置固定在中央
                                    E_bottom=-0.05;E_top=0.05;
                                else  % 收缩环较强，其位置不固定
                                    E_bottom=-0.1;E_top=0.1;
                                end

                                % 多高斯拟合
                                % 调用MultiGaussian_fit子函数计算Multi-Gaussian distribution
                                cell_type=2;
                                [rho0_best,rho_inf_array,E_array,w_array,break_point,wucha_best,TriGaussian_parameter]=MultiGaussian_fit(x_AP,y_AP,data_changdubi(number_changdubi),anterior_lateral_length,basal_length,posterior_lateral_length,cell_type);
                                framepath_filename_AP
                                % 调用MultiGaussianTransfertoEquivalentTriGaussian子函数计算等效Tri-Gaussian distribution
                                [EquivalentTriGaussian_parameter,area_sumsum]=MultiGaussianTransfertoEquivalentTriGaussian(x_AP,y_AP,rho0_best,rho_inf_array,E_array,w_array,break_point);

                                % 判断细胞是否是左侧为anterior part，右侧为posterior part。若不是，将拟合参数、拟合曲线、原始数据关于原点做对称处理
                                if rho_inf_array(end-1)<rho_inf_array(end)
                                    % MultiGaussian参数做对称处理
                                    rho_inf_array_reverse=zeros(size(rho_inf_array));
                                    E_array_reverse=zeros(size(rho_inf_array));
                                    w_array_reverse=zeros(size(rho_inf_array));
                                    for i_Gaussian=1:length(rho_inf_array)-2
                                        rho_inf_array_reverse(i_Gaussian)=rho_inf_array(end-i_Gaussian-1);
                                        E_array_reverse(i_Gaussian)=E_array(end-i_Gaussian-1);
                                        w_array_reverse(i_Gaussian)=w_array(end-i_Gaussian-1);
                                    end
                                    rho_inf_array_reverse(end-1)=rho_inf_array(end);rho_inf_array_reverse(end)=rho_inf_array(end-1);
                                    E_array_reverse(end-1)=E_array(end);E_array_reverse(end)=E_array(end-1);
                                    w_array_reverse(end-1)=w_array(end);w_array_reverse(end)=w_array(end-1);

                                    rho_inf_array=rho_inf_array_reverse;
                                    E_array=-E_array_reverse;
                                    w_array=w_array_reverse;

                                    % EquivalentTriGaussian参数做对称处理
                                    anterior_posterior_ring_reverse=EquivalentTriGaussian_parameter(5:6);
                                    EquivalentTriGaussian_parameter(5:6)=EquivalentTriGaussian_parameter(8:9);
                                    EquivalentTriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                    EquivalentTriGaussian_parameter(4)=-EquivalentTriGaussian_parameter(4);

                                    % TriGaussian参数做对称处理
                                    anterior_posterior_ring_reverse=TriGaussian_parameter(5:6);
                                    TriGaussian_parameter(5:6)=TriGaussian_parameter(8:9);
                                    TriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                    TriGaussian_parameter(4)=-TriGaussian_parameter(4);

                                    y_reverse_AP=zeros(length(y_AP),1);  % 原始数据翻转
                                    for ii=1:length(y_AP)
                                        y_reverse_AP(length(y_AP)-ii+1)=y_AP(ii);
                                    end
                                    y_AP=y_reverse_AP;
                                end

                                % 荧光强度值(纵坐标)归一化
                                % 将拟合完成的多高斯分布在定义域上的数学期望变换为1，同时将原始测量数据和等效三高斯分布按等比例拉伸
                                if length(rho_inf_array)-2==1
                                    f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                        +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2);
                                elseif length(rho_inf_array)-2==2
                                    f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                        +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2);
                                elseif length(rho_inf_array)-2==3
                                    f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                        +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2)...
                                        +rho_inf_array(5)*exp(-1/2*((y-E_array(5))/w_array(5)).^2);
                                elseif length(rho_inf_array)-2==4
                                    f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                        +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2)+...
                                        rho_inf_array(5)*exp(-1/2*((y-E_array(5))/w_array(5)).^2)+rho_inf_array(6)*exp(-1/2*((y-E_array(6))/w_array(6)).^2);
                                elseif length(rho_inf_array)-2>=5
                                    f=@(y) rho0_best+rho_inf_array(1)*exp(-1/2*((y-E_array(1))/w_array(1)).^2)+rho_inf_array(2)*exp(-1/2*((y-E_array(2))/w_array(2)).^2)...
                                        +rho_inf_array(3)*exp(-1/2*((y-E_array(3))/w_array(3)).^2)+rho_inf_array(4)*exp(-1/2*((y-E_array(4))/w_array(4)).^2)...
                                        +rho_inf_array(5)*exp(-1/2*((y-E_array(5))/w_array(5)).^2)+rho_inf_array(6)*exp(-1/2*((y-E_array(6))/w_array(6)).^2)...
                                        +rho_inf_array(7)*exp(-1/2*((y-E_array(7))/w_array(7)).^2);
                                end

                                a0=integral(f,-1/2,1/2);
                                rho0_best=rho0_best/a0;
                                rho_inf_array=rho_inf_array./a0;
                                EquivalentTriGaussian_parameter(2:3:8)=EquivalentTriGaussian_parameter(2:3:8)./a0;EquivalentTriGaussian_parameter(1)=EquivalentTriGaussian_parameter(1)./a0;
                                TriGaussian_parameter(2:3:8)=TriGaussian_parameter(2:3:8)./a0;TriGaussian_parameter(1)=TriGaussian_parameter(1)./a0;
                                y_MultiGaussian_fit_AP=ones(size(x_AP))*rho0_best;  % 多高斯拟合曲线
                                for i_Gaussian=1:length(rho_inf_array)
                                    y_MultiGaussian_fit_AP=y_MultiGaussian_fit_AP+rho_inf_array(i_Gaussian)*exp(-1/2*((x_AP-E_array(i_Gaussian))/w_array(i_Gaussian)).^2);
                                end
                                y_EquivalentTriGaussian_fit_AP=ones(size(x_AP))*EquivalentTriGaussian_parameter(1);  % 等效三高斯拟合曲线
                                for i_Gaussian=1:3
                                    y_EquivalentTriGaussian_fit_AP=y_EquivalentTriGaussian_fit_AP+EquivalentTriGaussian_parameter(3*i_Gaussian-1)*exp(-1/2*((x_AP-EquivalentTriGaussian_parameter(3*i_Gaussian+1))/EquivalentTriGaussian_parameter(3*i_Gaussian)).^2);
                                end
                                y_AP=y_AP./a0;

                                %%% 绘图
                                filename_plot=[filename_photo filename_kongge filename1 filename_kongge filename_cell filename_kongge filename2 filename_kongge filename_protein];

                                figure
                                plot(x_AP,y_AP,'.','MarkerSize',15,'Color','k');  % 绘制数据点
                                hold on;
                                plot(x_AP,y_MultiGaussian_fit_AP,'b','linewidth',3);  % 绘制拟合的高斯曲线
                                plot(x_AP,y_EquivalentTriGaussian_fit_AP,'r','linewidth',3);  % 绘制拟合的高斯曲线
                                hold off;
                                axis ([-0.5 0.5 0 max(y_AP)+0.1]);
                                title(filename_plot)
                                xlabel('x')
                                ylabel('cortex thickness')
                                set(gca,'YTick',0:0.2:max(y_AP)+0.1);  % 修改y坐标刻度
                                set(gca,'XTick',-0.5:0.1:0.5);  % 修改y坐标刻度
                                legend('Experiments','Multi-Gaussian Curve','Tri-Gaussian Curve')
                                set(gca,'FontName','Arial','FontSize',14)
                                saveas(gcf,[framepath2,filename_plot,'.jpg']);
                                close

                                %%% 拟合参数汇总
                                for aa=1
                                    if data_changdubi(number_changdubi)<=0.2
                                        data_changdubi_all(k_position1(k_AP),1,k_AP)=data_changdubi(number_changdubi);
                                        filenamearray_all{k_position1(k_AP),1,k_AP}=filename_AP;
                                        MultiGaussian_data_all(k_position1(k_AP),1,1,k_AP)=rho0_best;
                                        MultiGaussian_data_all(k_position1(k_AP),2:length(rho_inf_array)+1,1,k_AP)=rho_inf_array;
                                        MultiGaussian_data_all(k_position1(k_AP),length(rho_inf_array)+2:2*length(rho_inf_array)+1,1,k_AP)=w_array;
                                        MultiGaussian_data_all(k_position1(k_AP),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,1,k_AP)=E_array;
                                        EquivalentTriGaussian_data_all(k_position1(k_AP),:,1,k_AP)=EquivalentTriGaussian_parameter;
                                        TriGaussian_data_all(k_position1(k_AP),:,1,k_AP)=TriGaussian_parameter;
                                        wucha_best_all(k_position1(k_AP),1,k_AP)=wucha_best;
                                        k_position1(k_AP)=k_position1(k_AP)+1;
                                    elseif data_changdubi(number_changdubi)>0.2 && data_changdubi(number_changdubi)<=0.4
                                        data_changdubi_all(k_position2(k_AP),2,k_AP)=data_changdubi(number_changdubi);
                                        filenamearray_all{k_position2(k_AP),2,k_AP}=filename_AP;
                                        MultiGaussian_data_all(k_position2(k_AP),1,2,k_AP)=rho0_best;
                                        MultiGaussian_data_all(k_position2(k_AP),2:length(rho_inf_array)+1,2,k_AP)=rho_inf_array;
                                        MultiGaussian_data_all(k_position2(k_AP),length(rho_inf_array)+2:2*length(rho_inf_array)+1,2,k_AP)=w_array;
                                        MultiGaussian_data_all(k_position2(k_AP),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,2,k_AP)=E_array;
                                        EquivalentTriGaussian_data_all(k_position2(k_AP),:,2,k_AP)=EquivalentTriGaussian_parameter;
                                        TriGaussian_data_all(k_position2(k_AP),:,2,k_AP)=TriGaussian_parameter;
                                        wucha_best_all(k_position2(k_AP),2,k_AP)=wucha_best;
                                        k_position2(k_AP)=k_position2(k_AP)+1;
                                    elseif data_changdubi(number_changdubi)>0.4 && data_changdubi(number_changdubi)<=0.6
                                        data_changdubi_all(k_position3(k_AP),3,k_AP)=data_changdubi(number_changdubi);
                                        filenamearray_all{k_position3(k_AP),3,k_AP}=filename_AP;
                                        MultiGaussian_data_all(k_position3(k_AP),1,3,k_AP)=rho0_best;
                                        MultiGaussian_data_all(k_position3(k_AP),2:length(rho_inf_array)+1,3,k_AP)=rho_inf_array;
                                        MultiGaussian_data_all(k_position3(k_AP),length(rho_inf_array)+2:2*length(rho_inf_array)+1,3,k_AP)=w_array;
                                        MultiGaussian_data_all(k_position3(k_AP),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,3,k_AP)=E_array;
                                        EquivalentTriGaussian_data_all(k_position3(k_AP),:,3,k_AP)=EquivalentTriGaussian_parameter;
                                        TriGaussian_data_all(k_position3(k_AP),:,3,k_AP)=TriGaussian_parameter;
                                        wucha_best_all(k_position3(k_AP),3,k_AP)=wucha_best;
                                        k_position3(k_AP)=k_position3(k_AP)+1;
                                    elseif data_changdubi(number_changdubi)>0.6 && data_changdubi(number_changdubi)<=0.8
                                        data_changdubi_all(k_position4(k_AP),4,k_AP)=data_changdubi(number_changdubi);
                                        filenamearray_all{k_position4(k_AP),4,k_AP}=filename_AP;
                                        MultiGaussian_data_all(k_position4(k_AP),1,4,k_AP)=rho0_best;
                                        MultiGaussian_data_all(k_position4(k_AP),2:length(rho_inf_array)+1,4,k_AP)=rho_inf_array;
                                        MultiGaussian_data_all(k_position4(k_AP),length(rho_inf_array)+2:2*length(rho_inf_array)+1,4,k_AP)=w_array;
                                        MultiGaussian_data_all(k_position4(k_AP),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,4,k_AP)=E_array;
                                        EquivalentTriGaussian_data_all(k_position4(k_AP),:,4,k_AP)=EquivalentTriGaussian_parameter;
                                        TriGaussian_data_all(k_position4(k_AP),:,4,k_AP)=TriGaussian_parameter;
                                        wucha_best_all(k_position4(k_AP),4,k_AP)=wucha_best;
                                        k_position4(k_AP)=k_position4(k_AP)+1;
                                    else
                                        data_changdubi_all(k_position5(k_AP),5,k_AP)=data_changdubi(number_changdubi);
                                        filenamearray_all{k_position5(k_AP),5,k_AP}=filename_AP;
                                        MultiGaussian_data_all(k_position5(k_AP),1,5,k_AP)=rho0_best;
                                        MultiGaussian_data_all(k_position5(k_AP),2:length(rho_inf_array)+1,5,k_AP)=rho_inf_array;
                                        MultiGaussian_data_all(k_position5(k_AP),length(rho_inf_array)+2:2*length(rho_inf_array)+1,5,k_AP)=w_array;
                                        MultiGaussian_data_all(k_position5(k_AP),2*length(rho_inf_array)+2:3*length(rho_inf_array)+1,5,k_AP)=E_array;
                                        EquivalentTriGaussian_data_all(k_position5(k_AP),:,5,k_AP)=EquivalentTriGaussian_parameter;
                                        TriGaussian_data_all(k_position5(k_AP),:,5,k_AP)=TriGaussian_parameter;
                                        wucha_best_all(k_position5(k_AP),5,k_AP)=wucha_best;
                                        k_position5(k_AP)=k_position5(k_AP)+1;
                                    end
                                end
                            elseif k_AP==10  % kk=9表示该细胞仅有一种蛋白的数据文件
                            end
                            
                        elseif exist(framepath_filename_danyi,'file')
                            number_changdubi=number_changdubi+1;  % 跳过该细胞的长度比数据
                        else  % 若该张照片下，不存在该细胞的数据文件，说明该张照片中的细胞已经拟合完成，跳出j(细胞)循环
                            break
                        end
                        
                        % 覆盖该细胞拟合另一蛋白所用的kk循环参数值
                        filename3=num2str(k);
                        if k==1
                            filename_protein='lifeact';
                        elseif k==2
                            filename_protein='actin';
                        elseif k==3
                            filename_protein='MLC';
                        elseif k==4
                            filename_protein='RhoA';
                        elseif k==5
                            filename_protein='Anillin RBD';
                        elseif k==6
                            filename_protein='VD1';
                        elseif k==7
                            filename_protein='talinA';
                        elseif k==8
                            filename_protein='Cdc42WT';
                        elseif k==9
                            filename_protein='Cdc42D118A';
                        end
                        filename=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3];  % 如 0928_111
                        framepath_filename=[framepath2 filename];
                        framepath_filename_1=[framepath_filename filename_xiahuaxian filename_1 filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111_1.csv
                        framepath_filename_2=[framepath_filename filename_xiahuaxian filename_2 filename_csv];
                        framepath_filename_3=[framepath_filename filename_xiahuaxian filename_3 filename_csv];
                        framepath_filename_danyi=[framepath_filename filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111.csv
                    end
                    if j==1  % 若该张照片下不存在第一个细胞，则说明该日期下所有照片已经处理完成，跳出i(照片)循环
                        break
                    end
                end
            end
        end
    end
end



