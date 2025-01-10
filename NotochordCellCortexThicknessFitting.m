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

% �����ļ���char�ͱ���
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
filename_changdubi='���ȱ�';
filename_Radius='�뾶 (Radius, �������Ϊlumen���)';
filename_ContactAngle='�Ӵ��� (Contact Angle, �������ΪԲ�ܽ�)';

% ���岻ͬ����5��ʱ�ڵ��ļ�������(����д�뵽excel�ļ�ʱʹ��)
filenamearray_all=cell(1000,5,9); % �����ļ����ķ�����ܣ����е�һά��ʾĳһ������ĳһʱ�ڵ��������ڶ�ά��ʾʱ�ڣ�����ά��ʾ��������

% ���岻ͬ����5��ʱ�ڵ������������(����д�뵽excel�ļ�ʱʹ��)
EquivalentTriGaussian_data_all=zeros(1000,10,5,9); % ��Ч����˹��ϲ����ķ�����ܣ����е�һά��ʾĳһ������ĳһʱ�ڵ��������ڶ�ά��ʾ��ϲ���������ά��ʾʱ�ڣ�����ά��ʾ��������
TriGaussian_data_all=zeros(1000,10,5,9);  % ֱ������˹��ϲ����ķ������
MultiGaussian_data_all=zeros(1000,22,5,9);  % ���˹��ϲ����ķ������
data_changdubi_all=zeros(1000,5,9);  % ����ϸ��lumen/cell���ȱ����ݻ���
wucha_best_all=zeros(1000,5,9);
k_position1=ones(1,9);k_position2=ones(1,9);k_position3=ones(1,9);k_position4=ones(1,9);k_position5=ones(1,9);  % ����ͳ��ÿ��ʱ�ڲ��������data_all������ÿ��ʱ�ڶ�Ӧ��λ��

data_FittingParameterName=cell(1,10);
data_FittingParameterName{1,1}='�������';
data_FittingParameterName{1,2}='�����������ǿ��';
data_FittingParameterName{1,3}='������������ۼ��̶�';
data_FittingParameterName{1,4}='�������������ѧ����';
data_FittingParameterName{1,5}='ǰ��������ǿ��';
data_FittingParameterName{1,6}='ǰ���������ۼ��̶�';
data_FittingParameterName{1,7}='ǰ����������ѧ����';
data_FittingParameterName{1,8}='���������ǿ��';
data_FittingParameterName{1,9}='����������ۼ��̶�';
data_FittingParameterName{1,10}='�����������ѧ����';


%%% �����������ļ���·��������Ԥд��
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
    % lifeact xlsx�ļ�Ԥд��
    xlswrite(framepath_filename_lifeact,data_FittingParameterName,'sheet1','B1');
    % writematrix(M,'M.xls','Sheet',2,'Range','A3:E8')
    xlswrite(framepath_filename_lifeact,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_lifeact,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_lifeact,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_lifeact,data_FittingParameterName,'sheet5','B1');
    % actin xlsx�ļ�Ԥд��
    xlswrite(framepath_filename_actin,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_actin,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_actin,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_actin,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_actin,data_FittingParameterName,'sheet5','B1');
    % MLC xlsx�ļ�Ԥд��
    xlswrite(framepath_filename_MLC,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_MLC,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_MLC,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_MLC,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_MLC,data_FittingParameterName,'sheet5','B1');
    % RhoA xlsx�ļ�Ԥд��
    xlswrite(framepath_filename_RhoA,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_RhoA,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_RhoA,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_RhoA,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_RhoA,data_FittingParameterName,'sheet5','B1');
    % RBD xlsx�ļ�Ԥд��
    xlswrite(framepath_filename_RBD,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_RBD,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_RBD,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_RBD,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_RBD,data_FittingParameterName,'sheet5','B1');
    % VD1 xlsx�ļ�Ԥд��
    xlswrite(framepath_filename_VD1,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_VD1,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_VD1,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_VD1,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_VD1,data_FittingParameterName,'sheet5','B1');
    % talinA xlsx�ļ�Ԥд��
    xlswrite(framepath_filename_talinA,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_talinA,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_talinA,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_talinA,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_talinA,data_FittingParameterName,'sheet5','B1');
    % Cdc42WT xlsx�ļ�Ԥд��
    xlswrite(framepath_filename_Cdc42WT,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_Cdc42WT,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_Cdc42WT,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_Cdc42WT,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_Cdc42WT,data_FittingParameterName,'sheet5','B1');
    % Cdc42D118A xlsx�ļ�Ԥд��
    xlswrite(framepath_filename_Cdc42D118A,data_FittingParameterName,'sheet1','B1');
    xlswrite(framepath_filename_Cdc42D118A,data_FittingParameterName,'sheet2','B1');
    xlswrite(framepath_filename_Cdc42D118A,data_FittingParameterName,'sheet3','B1');
    xlswrite(framepath_filename_Cdc42D118A,data_FittingParameterName,'sheet4','B1');
    xlswrite(framepath_filename_Cdc42D118A,data_FittingParameterName,'sheet5','B1');
end

%%% �������ǰ�����ж��ļ�����ϸ�����Ƿ�һ��
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
            number_changdubi=0;  % ��������ļ�ʱָʾ���ȱ��ļ��е�����λ��
            
            filename_changdubi='���ȱ�';
            framepath2=[framepath num2str(k) filename_fanxiegang filename_month filename_day filename_fanxiegang];  % ���ļ���·������ 'C:\Users\Desktop\1\0928'
            framepath_filename_changdubi=[framepath2 filename_changdubi filename_xlsx];
            
            if exist(framepath_filename_changdubi,'file')   % ���ȱ��ļ��Ƿ���ڿ����ڿ���ɸѡ�Ƿ���ڸ�������õ��׵����ݣ����̲���ʱ��
                file_number=0;
                data_lumenlength=xlsread(framepath_filename_changdubi,'sheet1','A1:A100');  % ������ļ��µ�lumen����
                data_celllength=xlsread(framepath_filename_changdubi,'sheet1','B1:B100');  % ������ļ��µ�cell����
                data_changdubi=data_lumenlength./data_celllength;
                cell_number=length(data_changdubi);
                
                for i=1:30  % ���ѭ�����ڲ���Ŀ���ļ��Ƿ����(i��ʾ�ڼ���photo��j��ʾphoto�еĵڼ���ϸ����k��ʾӫ��ͨ��)
                    filename1=num2str(i);
                    for j=1:10
                        filename2=num2str(j);
                        filename=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3];  % �� 0928_111
                        framepath_filename=[framepath2 filename];
                        
                        % ƴ���������ļ�·����������һ��exist�����Ĳ���
                        framepath_filename_1=[framepath_filename filename_xiahuaxian filename_1 filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111_1.csv
                        framepath_filename_2=[framepath_filename filename_xiahuaxian filename_2 filename_csv];
                        framepath_filename_3=[framepath_filename filename_xiahuaxian filename_3 filename_csv];
                        framepath_filename_danyi=[framepath_filename filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111.csv
                        
                        if exist(framepath_filename_danyi,'file') || exist(framepath_filename_1,'file')
                            file_number=file_number+1;
                        end
                    end
                end
                if cell_number~=file_number
                    fprintf('ȱ�������ļ�');
                    framepath2
                    abc
                end
            end
        end
    end
end


%% Step 1
for k=1:9  % 1��ʾlifeact��2��ʾactin��3��ʾMLC��4��ʾRhoA��5��ʾAnillin RBD��6��ʾVD1,7��ʾtalinA��8��ʾCdc42WT��9��ʾCdc42D118A
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
            number_changdubi=0;  % ��������ļ�ʱָʾ���ȱ��ļ��е�����λ��
            
            framepath2=[framepath num2str(k) filename_fanxiegang filename_month filename_day filename_fanxiegang];  % ���ļ���·������ 'C:\Users\Desktop\1\0928'
            framepath_filename_changdubi=[framepath2 filename_changdubi filename_xlsx];
            framepath_filename_Radius=[framepath2 filename_Radius filename_xlsx];
            framepath_filename_ContactAngle=[framepath2 filename_ContactAngle filename_xlsx];
            
            if exist(framepath_filename_changdubi,'file')   % ���ȱ��ļ��Ƿ���ڿ����ڿ���ɸѡ�Ƿ���ڸ�������õ��׵����ݣ����̲���ʱ��
                data_lumenradius=xlsread(framepath_filename_Radius,'sheet1','A1:A100');  % ����lumen�뾶
                data_lumenradius=sqrt(data_lumenradius./pi);
                data_lumen_ContactAngle=xlsread(framepath_filename_ContactAngle,'sheet1','A1:A100');  % ����lumen�Ӵ���
                data_lumenlength=xlsread(framepath_filename_changdubi,'sheet1','A1:A100');  % ������ļ��µ�lumen����
                data_celllength=xlsread(framepath_filename_changdubi,'sheet1','B1:B100');  % ������ļ��µ�cell����
                % ��Lumen length����cell length����cell length=Lumen length+0.5
                for i=1:length(data_lumenradius)
                    if 2*data_lumenradius(i)*sin(data_lumen_ContactAngle(i)*pi./180)>data_celllength(i)
                        data_celllength(i)=2*data_lumenradius(i)*sin(data_lumen_ContactAngle(i)*pi./180)+0.5;
                    end
                end
                data_changdubi=2.*data_lumenradius.*sin(data_lumen_ContactAngle.*pi./180)./data_celllength;

                for i=1:30  % ���ѭ�����ڲ���Ŀ���ļ��Ƿ����(i��ʾ�ڼ���photo��j��ʾphoto�еĵڼ���ϸ����k��ʾӫ��ͨ��)
                    filename1=num2str(i);
                    for j=1:10
                        filename2=num2str(j);
                        filename=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3];  % �� 0928_111
                        framepath_filename=[framepath2 filename];
                        
                        % ƴ���������ļ�·����������һ��exist�����Ĳ���
                        framepath_filename_1=[framepath_filename filename_xiahuaxian filename_1 filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111_1.csv
                        framepath_filename_2=[framepath_filename filename_xiahuaxian filename_2 filename_csv];
                        framepath_filename_3=[framepath_filename filename_xiahuaxian filename_3 filename_csv];
                        framepath_filename_danyi=[framepath_filename filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111.csv
                        
                        if exist(framepath_filename_danyi,'file')   % һ��ϸ���������е����ļ�
                            number_changdubi=number_changdubi+1;
                            for aa=1
                                %%% ��ԭʼϸ���������Ҷ�ֵ���ݽ�������������ݵ�ɾ������
                                % ���ݲ����������basal-lateral domain˫��˹�ֲ�
                                % ���ݵ�����
                                data=xlsread(framepath_filename_danyi); % ��������
                                [X,Y]=size(data); % �������ݾ���ߴ�
                                if mod(X,2)==0   % ������Ϊż��������������Բ�ֵ�����������ƽ�ƹ����еĶԳ���
                                    data2=zeros(X+1,Y);
                                    data2(:,1)=imresize(data(:,1),[X+1,1], 'bilinear'); % ��X��������Բ�ֵ
                                    data2(:,2)=imresize(data(:,2),[X+1,1], 'bilinear'); % ��Y��������Բ�ֵ
                                else
                                    data2=data;
                                end
                                
                                % ����ϸ����anterior-lateral basal posterior-lateral domain�ĳ���
                                basal_lateral_length=data(end,1);  % ϸ����basal-lateral domain�ܳ���
                                anterior_lateral_length=data_celllength(number_changdubi)-data_lumenlength(number_changdubi);
                                posterior_lateral_length=anterior_lateral_length;
                                basal_length=basal_lateral_length-anterior_lateral_length-posterior_lateral_length;
                                
                                % X��ƽ��������ԭ��ԳƵ����䣬����һ��
                                bias=0.5; % ƽ����
                                x=data2(:,1)/max(data2(:,1)); % ��x����й�һ��
                                x=x-bias; % x��ƽ������������ԭ��Գ�
                                % Y�����ݹ�һ��
                                y=data2(:,2)/max(data2(:,2));
                                
                                
                                Y_anterior=zeros(1,round((X-1)*0.2));  % ����ϸ��������������Ҷ�ֵ���ݣ����ڵ���˹�ֲ����
                                Y_posterior=zeros(1,round((X-1)*0.2));
                                Y_anterior_AP=zeros(1,round((X-1)*0.2));  % ����ϸ����һ����������������Ҷ�ֵ���ݣ����ڵ���˹�ֲ����
                                Y_posterior_AP=zeros(1,round((X-1)*0.2));

                                if data_changdubi(number_changdubi)<=0.2  % lumen�γɳ�����������������δ��������ɾ������
                                    for k_AP=1:10
                                        if k_AP~=k
                                            filename3_AP=num2str(k_AP);  % AP: anothor protein
                                            filename_AP=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3_AP];  % �� 0928_111
                                            framepath_filename_AP=[framepath2 filename_AP];
                                            framepath_filename_danyi_AP=[framepath_filename_AP filename_csv];
                                            if exist(framepath_filename_danyi_AP,'file')
                                                break
                                            end
                                        end
                                    end
                                    if k_AP~=10  % ��ʾ��ϸ������һ����������
                                        % ��ȡ�õ��׵�����
                                        % ���ݵ�����
                                        data_AP=xlsread(framepath_filename_danyi_AP); % ��������
                                        [X_AP,Y_AP]=size(data_AP); % �������ݾ���ߴ�
                                        if mod(X_AP,2)==0   % ������Ϊż��������������Բ�ֵ�����������ƽ�ƹ����еĶԳ���
                                            data2_AP=zeros(X_AP+1,Y_AP);
                                            data2_AP(:,1)=imresize(data_AP(:,1),[X_AP+1,1], 'bilinear'); % ��X��������Բ�ֵ
                                            data2_AP(:,2)=imresize(data_AP(:,2),[X_AP+1,1], 'bilinear'); % ��Y��������Բ�ֵ
                                        else
                                            data2_AP=data_AP;
                                        end

                                        % X��ƽ��������ԭ��ԳƵ����䣬����һ��
                                        bias=0.5; % ƽ����
                                        x_AP=data2_AP(:,1)/max(data2_AP(:,1)); % ��x����й�һ��
                                        x_AP=x_AP-bias; % x��ƽ������������ԭ��Գ�
                                        % Y�����ݹ�һ��
                                        y_AP=data2_AP(:,2)/max(data2_AP(:,2));
                                    end
                                    x_min=-0.5;x_max=0.5;
                                else  % lumen�γ��к������������������ѽ�������ɾ������
                                    for k_AP=1:10
                                        if k_AP~=k
                                            filename3_AP=num2str(k_AP);  % AP: anothor protein
                                            filename_AP=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3_AP];  % �� 0928_111
                                            framepath_filename_AP=[framepath2 filename_AP];
                                            framepath_filename_danyi_AP=[framepath_filename_AP filename_csv];
                                            if exist(framepath_filename_danyi_AP,'file')
                                                break
                                            end
                                        end
                                    end
                                    if k_AP~=10  % ��ʾ��ϸ������һ����������
                                        % ��ȡ�õ��׵�����
                                        % ���ݵ�����
                                        data_AP=xlsread(framepath_filename_danyi_AP); % ��������
                                        [X_AP,Y_AP]=size(data_AP); % �������ݾ���ߴ�
                                        if mod(X_AP,2)==0   % ������Ϊż��������������Բ�ֵ�����������ƽ�ƹ����еĶԳ���
                                            data2_AP=zeros(X_AP+1,Y_AP);
                                            data2_AP(:,1)=imresize(data_AP(:,1),[X_AP+1,1], 'bilinear'); % ��X��������Բ�ֵ
                                            data2_AP(:,2)=imresize(data_AP(:,2),[X_AP+1,1], 'bilinear'); % ��Y��������Բ�ֵ
                                        else
                                            data2_AP=data_AP;
                                        end

                                        % X��ƽ��������ԭ��ԳƵ����䣬����һ��
                                        bias=0.5; % ƽ����
                                        x_AP=data2_AP(:,1)/max(data2_AP(:,1)); % ��x����й�һ��
                                        x_AP=x_AP-bias; % x��ƽ������������ԭ��Գ�
                                        % Y�����ݹ�һ��
                                        y_AP=data2_AP(:,2)/max(data2_AP(:,2));

                                        % Levenberg-Marquardt method��������˹�ֲ���ȷ��Ҫɾȥ�Ĳ���
                                        for ab=1
                                            if data_changdubi(number_changdubi)<=0.2  % lumen�γ����ڣ�lateral domain������δ����������9�������������
                                                E_lateral_start=0.3;
                                                E_lateral_end=0.5;
                                            else  % lumen�γ��к��ڣ�lateral domain����������������8�������������
                                                E_lateral_start=0.5;
                                                E_lateral_end=0.5;
                                            end

                                            wucha_best=100000000;
                                            % ʹ��Levenberg-Marquardt method���
                                            for i_fit=1:5  % ѭ��10�Σ�Ѱ��������Ͻ�
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
                                        % ��һ�������������
                                        for ab=1
                                            if data_changdubi(number_changdubi)<=0.2  % lumen�γ����ڣ�lateral domain������δ����������9�������������
                                                E_lateral_start=0.3;
                                                E_lateral_end=0.5;
                                            else  % lumen�γ��к��ڣ�lateral domain����������������8�������������
                                                E_lateral_start=0.5;
                                                E_lateral_end=0.5;
                                            end

                                            wucha_best=100000000;
                                            % ʹ��Levenberg-Marquardt method���
                                            for i_fit=1:5  % ѭ��10�Σ�Ѱ��������Ͻ�
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
                                        x_min=min(E_anterior_best,E_anterior_AP_best);  % ǰ��߽�ȡ�������׵����ǰ�߽��н�С��
                                        x_max=max(E_posterior_best,E_posterior_AP_best);  % ���߽�ȡ�������׵���Ϻ�߽��нϴ��
                                    else  % kk=9��ʾ��ϸ������һ�ֵ��׵������ļ�
                                        % ��ȡ�õ��׵�����
                                        % ���ݵ�����

                                        % Levenberg-Marquardt method��������˹�ֲ���ȷ��Ҫɾȥ�Ĳ���
                                        for ab=1
                                            if data_changdubi(number_changdubi)<=0.2  % lumen�γ����ڣ�lateral domain������δ����������9�������������
                                                E_lateral_start=0.3;
                                                E_lateral_end=0.5;
                                            else  % lumen�γ��к��ڣ�lateral domain����������������8�������������
                                                E_lateral_start=0.5;
                                                E_lateral_end=0.5;
                                            end

                                            wucha_best=100000000;
                                            % ʹ��Levenberg-Marquardt method���
                                            for i_fit=1:5  % ѭ��10�Σ�Ѱ��������Ͻ�
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
                                % ��x_min��x_maxƫ��߽����(������ϸ�����ȵ�1/10)��˵���������Ч�����ã���ȡԭʼ�߽�ֵ
                                if (x_min+0.5)/0.5>0.1
                                    x_min=-0.5;
                                elseif (0.5-x_max)/0.5>0.1
                                    x_max=0.5;
                                end
                                
                                %%% ��ɾ��������ɵ�ϸ���Ҷ�ֵ���ݽ�������˹�ֲ����
                                % �ҵ�ԭX���Ϻ�x_min x_max��ӽ��ĵ�
                                wucha_best=100;
                                for ii=1:X
                                    wucha=abs(x(ii)-x_min);
                                    if wucha<wucha_best
                                        wucha_best=wucha;
                                        x_min_best=x(ii);  % ����ӽ�x_min��ԭX����ֵ��ʱ������x_min_best������
                                    end
                                end
                                x_min=x_min_best;  % x_min�ĳ�X������ӽ�x_min���ֵ
                                wucha_best=100;
                                for ii=1:X
                                    wucha=abs(x(ii)-x_max);
                                    if wucha<wucha_best
                                        wucha_best=wucha;
                                        x_max_best=x(ii);  % ����ӽ�x_max��ԭX����ֵ��ʱ������x_max_best������
                                    end
                                end
                                x_max=x_max_best;  % x_max�ĳ�X������ӽ�x_max���ֵ
                                
                                % ɾȥX��Ͷ�Ӧ�Ҷ�ֵ(y)��С��x_min�ʹ���x_max�Ĳ���
                                x_delete=zeros(round((X-1)*(x_max-x_min)+1),1);
                                y_delete=zeros(round((X-1)*(x_max-x_min)+1),1);
                                for ii=1:round((X-1)*(x_max-x_min)+1)
                                    x_delete(ii)=x(ii+round((x_min+0.5)*(X-1)));
                                    y_delete(ii)=y(ii+round((x_min+0.5)*(X-1)));
                                end
                                x=x_delete;
                                y=y_delete;
                                % X������ƽ��������ԭ��ԳƵ����䣬��������(-0.5,0.5)
                                bias=x_max+x_min; % ƽ����
                                stretch=1/(x_max-x_min); % ������
                                x=x-bias/2; % x��ƽ������������ԭ��Գ�
                                x=x*stretch;
                                
                                %%% ʹ��Levenberg-Marquardt method���
                                max_y=0;min_y=100;
                                for iy=floor(length(y)*0.25):floor(length(y)*0.75)
                                    if y(iy)>max_y
                                        max_y=y(iy);
                                    elseif y(iy)<min_y
                                        min_y=y(iy);
                                    end
                                end
                                if max_y/min_y<1.5  % ��ʾ������������Ѿ�������ʧ��lumenʱ�ڴ��ں��ڣ�Ϊ��ֹ��ϳ�������������λ�ù̶�������
                                    E_bottom=-0.05;E_top=0.05;
                                else  % ��������ǿ����λ�ò��̶�
                                    E_bottom=-0.1;E_top=0.1;
                                end

                                % ���˹���
                                % ����MultiGaussian_fit�Ӻ�������Multi-Gaussian distribution
                                cell_type=1;
                                [rho0_best,rho_inf_array,E_array,w_array,break_point,wucha_best,TriGaussian_parameter]=MultiGaussian_fit(x,y,data_changdubi(number_changdubi),anterior_lateral_length,basal_length,posterior_lateral_length,cell_type);
                                framepath_filename
                                % ����MultiGaussianTransfertoEquivalentTriGaussian�Ӻ��������ЧTri-Gaussian distribution
                                [EquivalentTriGaussian_parameter,area_sumsum]=MultiGaussianTransfertoEquivalentTriGaussian(x,y,rho0_best,rho_inf_array,E_array,w_array,break_point);
                            end
                            
                            % �ж�ϸ���Ƿ������Ϊanterior part���Ҳ�Ϊposterior part�������ǣ�����ϲ�����������ߡ�ԭʼ���ݹ���ԭ�����Գƴ���
                            if rho_inf_array(end-1)<rho_inf_array(end)
                                % MultiGaussian�������Գƴ���
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

                                % EquivalentTriGaussian�������Գƴ���
                                anterior_posterior_ring_reverse=EquivalentTriGaussian_parameter(5:6);
                                EquivalentTriGaussian_parameter(5:6)=EquivalentTriGaussian_parameter(8:9);
                                EquivalentTriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                EquivalentTriGaussian_parameter(4)=-EquivalentTriGaussian_parameter(4);

                                % TriGaussian�������Գƴ���
                                anterior_posterior_ring_reverse=TriGaussian_parameter(5:6);
                                TriGaussian_parameter(5:6)=TriGaussian_parameter(8:9);
                                TriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                TriGaussian_parameter(4)=-TriGaussian_parameter(4);
                                
                                y_reverse=zeros(length(y),1);  % ԭʼ���ݷ�ת
                                for ii=1:length(y)
                                    y_reverse(length(y)-ii+1)=y(ii);
                                end
                                y=y_reverse;
                            end

                            % ӫ��ǿ��ֵ(������)��һ��
                            % �������ɵĶ��˹�ֲ��ڶ������ϵ���ѧ�����任Ϊ1��ͬʱ��ԭʼ�������ݺ͵�Ч����˹�ֲ����ȱ�������
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
                            y_MultiGaussian_fit=ones(size(x))*rho0_best;  % ���˹�������
                            for i_Gaussian=1:length(rho_inf_array)
                                y_MultiGaussian_fit=y_MultiGaussian_fit+rho_inf_array(i_Gaussian)*exp(-1/2*((x-E_array(i_Gaussian))/w_array(i_Gaussian)).^2);
                            end
                            y_EquivalentTriGaussian_fit=ones(size(x))*EquivalentTriGaussian_parameter(1);  % ��Ч����˹�������
                            for i_Gaussian=1:3
                                y_EquivalentTriGaussian_fit=y_EquivalentTriGaussian_fit+EquivalentTriGaussian_parameter(3*i_Gaussian-1)*exp(-1/2*((x-EquivalentTriGaussian_parameter(3*i_Gaussian+1))/EquivalentTriGaussian_parameter(3*i_Gaussian)).^2);
                            end
                            y=y./a0;
                            

                            %%% ��ͼ
                            filename_plot=[filename_photo filename_kongge filename1 filename_kongge filename_cell filename_kongge filename2 filename_kongge filename_protein];
                            
                            figure
                            plot(x,y,'.','MarkerSize',15,'Color','k');  % �������ݵ�
                            hold on;
                            plot(x,y_MultiGaussian_fit,'b','linewidth',3);  % ������ϵĸ�˹����
                            plot(x,y_EquivalentTriGaussian_fit,'r','linewidth',3);  % ������ϵĸ�˹����
                            hold off;
                            axis ([-0.5 0.5 0 max(y)+0.1]);
                            title(filename_plot)
                            xlabel('x')
                            ylabel('cortex thickness')
                            set(gca,'YTick',0:0.2:max(y)+0.1);  % �޸�y����̶�
                            set(gca,'XTick',-0.5:0.1:0.5);  % �޸�y����̶�
                            legend('Experiments','Multi-Gaussian Curve','Tri-Gaussian Curve')
                            set(gca,'FontName','Arial','FontSize',12)
                            saveas(gcf,[framepath2,filename_plot,'.jpg']);
                            close

                            %%% ��ϲ�������
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
                            
                            %%% ��ϸ�ϸ������һ�����������ļ�
                            if k_AP~=10  % ��ϸ��������һ�����������ļ�
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
                                
                                %%% ��ɾ��������ɵ�ϸ���Ҷ�ֵ���ݽ�������˹�ֲ����
                                % ɾȥX��Ͷ�Ӧ�Ҷ�ֵ(y)��С��x_min�ʹ���x_max�Ĳ���
                                x_delete_AP=zeros(round((X-1)*(x_max-x_min)+1),1);
                                y_delete_AP=zeros(round((X-1)*(x_max-x_min)+1),1);
                                for ii=1:round((X-1)*(x_max-x_min)+1)
                                    x_delete_AP(ii)=x_AP(ii+round((x_min+0.5)*(X-1)));
                                    y_delete_AP(ii)=y_AP(ii+round((x_min+0.5)*(X-1)));
                                end
                                x_AP=x_delete_AP;
                                y_AP=y_delete_AP;
                                % X������ƽ��������ԭ��ԳƵ����䣬��������(-0.5,0.5)
                                bias=x_max+x_min; % ƽ����
                                stretch=1/(x_max-x_min); % ������
                                x_AP=x_AP-bias/2; % x��ƽ������������ԭ��Գ�
                                x_AP=x_AP*stretch;

                                %%% ʹ��Levenberg-Marquardt method���
                                max_y=0;min_y=100;
                                for iy=floor(length(y_AP)*0.25):floor(length(y_AP)*0.75)
                                    if y_AP(iy)>max_y
                                        max_y=y_AP(iy);
                                    elseif y_AP(iy)<min_y
                                        min_y=y_AP(iy);
                                    end
                                end
                                if max_y/min_y<1.5  % ��ʾ������������Ѿ�������ʧ��lumenʱ�ڴ��ں��ڣ�Ϊ��ֹ��ϳ�������������λ�ù̶�������
                                    E_bottom=-0.05;E_top=0.05;
                                else  % ��������ǿ����λ�ò��̶�
                                    E_bottom=-0.1;E_top=0.1;
                                end

                                % ���˹���
                                % ����MultiGaussian_fit�Ӻ�������Multi-Gaussian distribution
                                cell_type=1;
                                [rho0_best,rho_inf_array,E_array,w_array,break_point,wucha_best,TriGaussian_parameter]=MultiGaussian_fit(x_AP,y_AP,data_changdubi(number_changdubi),anterior_lateral_length,basal_length,posterior_lateral_length,cell_type);
                                framepath_filename_AP
                                % ����MultiGaussianTransfertoEquivalentTriGaussian�Ӻ��������ЧTri-Gaussian distribution
                                [EquivalentTriGaussian_parameter,area_sumsum]=MultiGaussianTransfertoEquivalentTriGaussian(x_AP,y_AP,rho0_best,rho_inf_array,E_array,w_array,break_point);

                                % �ж�ϸ���Ƿ������Ϊanterior part���Ҳ�Ϊposterior part�������ǣ�����ϲ�����������ߡ�ԭʼ���ݹ���ԭ�����Գƴ���
                                if rho_inf_array(end-1)<rho_inf_array(end)
                                    % MultiGaussian�������Գƴ���
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

                                    % EquivalentTriGaussian�������Գƴ���
                                    anterior_posterior_ring_reverse=EquivalentTriGaussian_parameter(5:6);
                                    EquivalentTriGaussian_parameter(5:6)=EquivalentTriGaussian_parameter(8:9);
                                    EquivalentTriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                    EquivalentTriGaussian_parameter(4)=-EquivalentTriGaussian_parameter(4);

                                    % TriGaussian�������Գƴ���
                                    anterior_posterior_ring_reverse=TriGaussian_parameter(5:6);
                                    TriGaussian_parameter(5:6)=TriGaussian_parameter(8:9);
                                    TriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                    TriGaussian_parameter(4)=-TriGaussian_parameter(4);

                                    y_reverse_AP=zeros(length(y_AP),1);  % ԭʼ���ݷ�ת
                                    for ii=1:length(y_AP)
                                        y_reverse_AP(length(y_AP)-ii+1)=y_AP(ii);
                                    end
                                    y_AP=y_reverse_AP;
                                end

                                % ӫ��ǿ��ֵ(������)��һ��
                                % �������ɵĶ��˹�ֲ��ڶ������ϵ���ѧ�����任Ϊ1��ͬʱ��ԭʼ�������ݺ͵�Ч����˹�ֲ����ȱ�������
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
                                y_MultiGaussian_fit_AP=ones(size(x_AP))*rho0_best;  % ���˹�������
                                for i_Gaussian=1:length(rho_inf_array)
                                    y_MultiGaussian_fit_AP=y_MultiGaussian_fit_AP+rho_inf_array(i_Gaussian)*exp(-1/2*((x_AP-E_array(i_Gaussian))/w_array(i_Gaussian)).^2);
                                end
                                y_EquivalentTriGaussian_fit_AP=ones(size(x_AP))*EquivalentTriGaussian_parameter(1);  % ��Ч����˹�������
                                for i_Gaussian=1:3
                                    y_EquivalentTriGaussian_fit_AP=y_EquivalentTriGaussian_fit_AP+EquivalentTriGaussian_parameter(3*i_Gaussian-1)*exp(-1/2*((x_AP-EquivalentTriGaussian_parameter(3*i_Gaussian+1))/EquivalentTriGaussian_parameter(3*i_Gaussian)).^2);
                                end
                                y_AP=y_AP./a0;

                                %%% ��ͼ
                                filename_plot=[filename_photo filename_kongge filename1 filename_kongge filename_cell filename_kongge filename2 filename_kongge filename_protein];

                                figure
                                plot(x_AP,y_AP,'.','MarkerSize',15,'Color','k');  % �������ݵ�
                                hold on;
                                plot(x_AP,y_MultiGaussian_fit_AP,'b','linewidth',3);  % ������ϵĸ�˹����
                                plot(x_AP,y_EquivalentTriGaussian_fit_AP,'r','linewidth',3);  % ������ϵĸ�˹����
                                hold off;
                                axis ([-0.5 0.5 0 max(y_AP)+0.1]);
                                title(filename_plot)
                                xlabel('x')
                                ylabel('cortex thickness')
                                set(gca,'YTick',0:0.2:max(y_AP)+0.1);  % �޸�y����̶�
                                set(gca,'XTick',-0.5:0.1:0.5);  % �޸�y����̶�
                                legend('Experiments','Multi-Gaussian Curve','Tri-Gaussian Curve')
                                set(gca,'FontName','Arial','FontSize',12)
                                saveas(gcf,[framepath2,filename_plot,'.jpg']);
                                close

                                %%% ��ϲ�������
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
                            elseif k_AP==10  % kk=10��ʾ��ϸ������һ�ֵ��׵������ļ�
                            end

                        elseif exist(framepath_filename_2,'file')
                            number_changdubi=number_changdubi+1;  % ������ϸ���ĳ��ȱ�����
                        else  % ��������Ƭ�£������ڸ�ϸ���������ļ���˵��������Ƭ�е�ϸ���Ѿ������ɣ�����j(ϸ��)ѭ��
                            break
                        end
                        
                        % ���Ǹ�ϸ�������һ�������õ�kkѭ������ֵ
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
                        filename=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3];  % �� 0928_111
                        framepath_filename=[framepath2 filename];
                        framepath_filename_1=[framepath_filename filename_xiahuaxian filename_1 filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111_1.csv
                        framepath_filename_2=[framepath_filename filename_xiahuaxian filename_2 filename_csv];
                        framepath_filename_3=[framepath_filename filename_xiahuaxian filename_3 filename_csv];
                        framepath_filename_danyi=[framepath_filename filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111.csv
                    end
                    if j==1  % ��������Ƭ�²����ڵ�һ��ϸ������˵����������������Ƭ�Ѿ�������ɣ�����i(��Ƭ)ѭ��
                        break
                    end
                end
            end
        end
    end
end


%%% ɾ����Ͻ�����ѵ�ϸ�����ݣ��ͳ����������������������ƫ��ƽ����ϸ��,������ƽ��A-P��ǿ�ȱȵļ���
% ��wucha_best_all�����ÿ��ʱ�ڡ�ÿ�ֵ�����С�����������(ʹ��ð������ bubble sort)������˳��ͬʱ������filenamearray_all data_all����
% ��mix_all��filenamearray_all����ϲ�
mix_all=zeros(1000,25,5,9); % ����ĳ��ʱ�ڡ�ĳ�ֵ����µ�wucha_best filename data 1-9 changdubi

% mix_all����ֵ
for i=1:5 % i��ʾ��ͬʱ��
    for j=1:9 % j��ʾ��ͬ����
        for k=1:1000  % �����ֵ���ļ������������������mix_all����
            mix_all(k,1,i,j)=wucha_best_all(k,i,j);
            mix_all(k,2,i,j)=0;
            for ii=1:length(MultiGaussian_data_all(1,:,1,1))
                mix_all(k,ii+2,i,j)=MultiGaussian_data_all(k,ii,i,j);
            end
            mix_all(k,end,i,j)=data_changdubi_all(k,i,j);
        end
    end
end

% mix_sort����ֵ
mix_sort=mix_all;  % ��������ľ���

% ʹ��ð�����򷨽�������
% ���ȸ�������С���򣬲������ֵ����8��ϸ�������޳�
cell_number=zeros(5,9);
for i=1:5
    for j=1:9
        for k=1:1000  % �ҳ��þ����е�ʵ����������
            if mix_sort(k,1,i,j)==0
                break
            end
        end
        cell_number(i,j)=k-1; % ��ʵ��������(ϸ������)��ֵ��cell_number
        
        % ð������
        for ii=2:cell_number(i,j)
            for jj=cell_number(i,j):-1:ii
                if mix_sort(jj,1,i,j)<mix_sort(jj-1,1,i,j)  % ��ǰһ��ϸ���������ں�һ��ϸ�����������˳��
                    temp=mix_sort(jj,:,i,j);
                    mix_sort(jj,:,i,j)=mix_sort(jj-1,:,i,j);
                    mix_sort(jj-1,:,i,j)=temp;
                end
            end
        end
        
        % ɾ�����ֵ����8��ϸ������
        for k=1:cell_number(i,j)
            if mix_sort(k,1,i,j)>=8
                for ii=1:length(mix_sort(1,:,1,1))
                    mix_sort(k,ii,i,j)=0;  % �޳�������8��ϸ������
                end
                cell_number(i,j)=cell_number(i,j)-1;
            end
        end
    end
end

% ����basal domain���˹��ƽ��ǿ�Ⱥ;ۼ��̶ȣ����浽mix_sort_basalave�����У����ں������������ɾ������
mix_sort_basalave=zeros(1000,13,5,9);
for i=1:5
    for j=1:9
        for k=1:cell_number(i,j)
            for ii=4:length(mix_sort(1,:,1,1))
                if mix_sort(k,ii,i,j)==0  % ����iiλΪ0��˵��basal domain��(ii-10)/3����˹
                    break
                end
            end
            if mix_sort(k,end-1,i,j)~=0  % �����һλ��Ϊ0��˵��basal domain��5(max)����˹
                ii=25;
            end
            basal_Gaussian_number=(ii-10)/3;

            % ����basal domain���˹��ƽ��ǿ�Ⱥ;ۼ��̶�
            basal_overactivity_average=sum(mix_all(k,4:basal_Gaussian_number+3,i,j))/basal_Gaussian_number;
            basal_signalwidth_average=sum(mix_all(k,basal_Gaussian_number+4:2*basal_Gaussian_number+3,i,j))/basal_Gaussian_number;
            basal_expectation_average=sum(mix_all(k,2*basal_Gaussian_number+4:3*basal_Gaussian_number+3,i,j))/basal_Gaussian_number;

            % ��ֵ��mix_sort_basalave
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

% �ֱ�����������ǿ�Ⱥ;ۼ��̶�(basal domain���˹��ƽ��ǿ�Ⱥ;ۼ��̶�),��ǰ����������ǿ��(һ�����໷��ϲ��ѵ�ϸ����ۼ��̶Ⱥ�ǿ�������ǿ)����
% �޳�����ǰ10%�ͺ�10%��ϸ�����ݵĲ���
memory_cell_wucha=zeros(200,5,9,3);  % �����ɾ��ϸ�����ݶ�Ӧ�����ֵ
k_cellname=ones(5,9,3);  % memory_cellname�����鳤��
number_sort_position=[4 5 7];  % ����������Ӧ��mix_all�����λ��
for k=1:3
    number_sort=number_sort_position(k);
    for i=1:5
        for j=1:9
            % ð������
            for ii=2:cell_number(i,j)
                for jj=cell_number(i,j):-1:ii
                    if mix_sort_basalave(jj,number_sort_position(k),i,j)<mix_sort_basalave(jj-1,number_sort_position(k),i,j)  % ��ǰһ��ϸ���������������ں�һ��ϸ�������������������˳��
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
                % �ҳ�ǰ10%��ϸ���������ֵ���
                for jj=1:floor(cell_number(i,j)/10)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
                % �ҳ���10%��ϸ���������ֵ���
                for jj=cell_number(i,j)-floor(cell_number(i,j)/10)+1:cell_number(i,j)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
            else
                % �ҳ�ǰ10%��ϸ���������ֵ���
                for jj=1:floor(cell_number(i,j)/10)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
                % �ҳ���10%��ϸ���������ֵ���
                for jj=cell_number(i,j)-floor(cell_number(i,j)/5)+1:cell_number(i,j)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
            end
        end
    end
end
% ȡ��������õ���ϸ���Ĳ���������mix_all������ɾȥ
mix_sort2=mix_sort;  % ��mix_all�������ݸ�ֵ��mix_all2
mix_sort=zeros(1000,25,5,9);  % ������飬���ں������ɾ��������ϸ�����ϸ���������
for i=1:5
    for j=1:9
        kk2=1;  % ���ڴ洢mix_allÿ��ʱ�ڣ�ÿ�ֵ��׵�ϸ����
        for k1=1:1000
            if mix_sort2(k1,1,i,j)~=0  % ���ڵ�k1��ϸ��
                logic=0;  % �����жϸ�ϸ���Ƿ���Ҫɾ��
                for k2=1:k_cellname(i,j,1)
                    if mix_sort2(k1,1,i,j)==memory_cell_wucha(k2,i,j,1)
                        logic=1;  % ����ϸ����wuchaֵΪ���������������ǿ��,��ǰ����������ǿ�ȵ�ǰ10%�ͺ�10%��ϸ������˵����ϸ����Ҫɾ��
                        break
                    end
                end
                if logic==0  % ����һ��δ�ҵ�����ͬ���ֵ��ϸ�����ٽ��еڶ��β���
                    for k2=1:k_cellname(i,j,2)
                        if mix_sort2(k1,1,i,j)==memory_cell_wucha(k2,i,j,2)
                            logic=1;  % ����ϸ����wuchaֵΪ��������������ľۼ��̶�,��ǰ����������ǿ�ȵ�ǰ10%�ͺ�10%��ϸ������˵����ϸ����Ҫɾ��
                            break
                        end
                    end
                end
                if logic==0  % ����һ������δ�ҵ�����ͬ���ֵ��ϸ�����ٽ��е����β���
                    for k2=1:k_cellname(i,j,3)
                        if mix_sort2(k1,1,i,j)==memory_cell_wucha(k2,i,j,3)
                            logic=1;  % ����ϸ����wuchaֵΪ��������������ľۼ��̶�,��ǰ����������ǿ�ȵ�ǰ10%�ͺ�10%��ϸ������˵����ϸ����Ҫɾ��
                            break
                        end
                    end
                end
                if logic==0  % ��logicΪ0�����ʾ��ϸ������ɾ�������丳ֵ��mix_sort filenamearray_all����
                    mix_sort(kk2,:,i,j)=mix_sort2(k1,:,i,j);
                    kk2=kk2+1;
                else  % ��logicΪ1�����ʾ��ϸ����ɾ��������ֵ��mix_all������cell_number��1
                    cell_number(i,j)=cell_number(i,j)-1;
                end
            else
                break
            end
        end
    end
end

% ��������ɾ����ɵ�mix_sort�������¸�ֵ��mix_all����
mix_all=mix_sort;

% ���¼���basal domainƽ�������ϲ������󣬲���wucha_best����
mix_all_basalave=zeros(1000,13,5,9);
for i=1:5
    for j=1:9
        for k=1:cell_number(i,j)
            for ii=4:length(mix_all(1,:,1,1))
                if mix_all(k,ii,i,j)==0  % ����iiλΪ0��˵��basal domain��(ii-10)/3����˹
                    break
                end
            end
            if mix_all(k,end-1,i,j)~=0  % �����һλ��Ϊ0��˵��basal domain��5(max)����˹
                ii=25;
            end
            basal_Gaussian_number=(ii-10)/3;

            % ����basal domain���˹��ƽ��ǿ�Ⱥ;ۼ��̶�
            basal_overactivity_average=sum(mix_all(k,4:basal_Gaussian_number+3,i,j))/basal_Gaussian_number;
            basal_signalwidth_average=sum(mix_all(k,basal_Gaussian_number+4:2*basal_Gaussian_number+3,i,j))/basal_Gaussian_number;
            basal_expectation_average=sum(mix_all(k,2*basal_Gaussian_number+4:3*basal_Gaussian_number+3,i,j))/basal_Gaussian_number;

            % ��ֵ��mix_sort_basalave
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
% ð������
for i=1:5
    for j=1:9
        for ii=2:cell_number(i,j)
            for jj=cell_number(i,j):-1:ii
                if mix_all_basalave(jj,1,i,j)<mix_all_basalave(jj-1,1,i,j)  % ��ǰһ��ϸ���������ں�һ��ϸ�����������˳��
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


%%% ����A-P��ǰ������������ǿ�ȵı�ֵ(���ں�����ƬA-P polarity direction�жϺ͵��ӹ۲����ݲ��)
AP_overactivity_ratio=zeros(5,9);  % ����A-P��������ǿ�ȱȱ�����ÿ��ʱ�ڡ�ÿ�ֵ��׵�������
max_AP_overactivity_ratio=zeros(5,9);  % ����A-P��������ǿ�ȱ����ֵ����
for i=1:5
    for j=1:9
        for k=1:1000
            if abs(mix_all_basalave(k,4,i,j))+abs(mix_all_basalave(k,7,i,j))+abs(mix_all_basalave(k,10,i,j))~=0  % ���ڸ�ϸ������ϲ���
                if mix_all_basalave(k,10,i,j)~=0  % ���ȼ���A-P��ֵ���ڵ�����
                    AP_overactivity_ratio(i,j)=AP_overactivity_ratio(i,j)+mix_all_basalave(k,7,i,j)/mix_all_basalave(k,10,i,j);  % �ҳ�����A-P��ǿ�ȱ�ֵ
                    if max_AP_overactivity_ratio(i,j)<mix_all_basalave(k,7,i,j)/mix_all_basalave(k,10,i,j)
                        max_AP_overactivity_ratio(i,j)=mix_all_basalave(k,7,i,j)/mix_all_basalave(k,10,i,j);
                    end
                end
            else  % �����ڸ�ϸ������ϲ���������ѭ��
                break
            end
            if wucha_best_all(k,i,j)>180  % ��ʾ��ϸ���Լ���������ϸ����Ͻ���ϲ�����ڼ���
                break
            end
        end
        k=k-1;
        if k==0
            AP_overactivity_ratio(i,j)=0;
        else
            for k_2=1:k  % ��A-P��ǿ�ȱ�ֵ�����ڣ������滻Ϊ����õ����ֵ����
                if mix_all_basalave(k,10,i,j)==0
                    AP_overactivity_ratio(i,j)=AP_overactivity_ratio(i,j)+max_AP_overactivity_ratio(i,j);
                end
            end
            AP_overactivity_ratio(i,j)=AP_overactivity_ratio(i,j)/k;
        end
    end
end


%%% ����A-P��PCPǿ�ȱ�(C_AP) (C_AP���壺ǰ���������ĸ�˹ǿ�ȳ��Ժ���������ĸ�˹ǿ��)
%%% protein_all_fit�ĵ���άΪrho_anterior��rho_posterior�Ĺ�ϵ���˹�ϵ��ΪC_AP
%%% ����ǰ���������ͺ����������rho��w����Խ�ǿ�����ǽ�w��ʾΪrho�ĺ�������C_AP��ʾΪ������rho�ĺ���(��ǰ������ϸ������ʱ�ڲ�ͬ����C_AP_next/(C_AP_next+1)*(C_AP_now+1)����)
%%% ������w����rho�ع麯���ı�׼��

% ����rho_anterior��rho_posterior��ֵ��ʱ��ı仯
AP_rho_ratio_time=zeros(2000,9,2);  % ����A-P polarity coefficient
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
    
    % ð������(����TD/CD)
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
        if i<=(win-1)/2  % ��tΪ��ʼ����ʱ��
            AP_rho_ratio_MA(i,1)=sum(AP_rho_ratio_temp(1:2*i-1,1))/(2*i-1);
            AP_rho_ratio_MA(i,2)=sum(AP_rho_ratio_temp(1:2*i-1,2))/(2*i-1);
        elseif i<sum(cell_number(:,j))-i_delete-(win-1)/2+1  % ��tΪ�м�ʱ��
            AP_rho_ratio_MA(i,1)=sum(AP_rho_ratio_temp(i-(win-1)/2:i+(win-1)/2,1))/win;
            AP_rho_ratio_MA(i,2)=sum(AP_rho_ratio_temp(i-(win-1)/2:i+(win-1)/2,2))/win;
        else  % ��tΪĩβ����ʱ��
            temp=sum(cell_number(:,j))-i_delete-i+1;
            AP_rho_ratio_MA(i,1)=sum(AP_rho_ratio_temp(sum(cell_number(:,j))-i_delete-2*temp+2:sum(cell_number(:,j))-i_delete,1))/(2*temp-1);
            AP_rho_ratio_MA(i,2)=sum(AP_rho_ratio_temp(sum(cell_number(:,j))-i_delete-2*temp+2:sum(cell_number(:,j))-i_delete,2))/(2*temp-1);
        end
    end
    
    N=1;  % ��Ͻ���
    [p,S]=polyfit(AP_rho_ratio_MA(:,1),AP_rho_ratio_MA(:,2),N);
    [yfit,delta]=polyconf(p,AP_rho_ratio_MA(:,1),S,'alpha',0.05,'predopt','curve');
    x_dummy=linspace(min(AP_rho_ratio_MA(:,1)),max(AP_rho_ratio_MA(:,1)),2000);  % ��Ϻ�����x����
    delta_interp=interp1(AP_rho_ratio_MA(:,1),delta,x_dummy,'linear');  % ����������������ֵ
    AP_rho_ratio_time(:,j,1)=x_dummy;
    for i=1:N+1
        AP_rho_ratio_time(:,j,2)=AP_rho_ratio_time(:,j,2)+p(i).*x_dummy'.^(N-i+1);
    end
end

% ��һ�ֵ��׵�5��ʱ���������ϵ�һά������(protein_mix_all)
protein_mix_all=zeros(5000,13,9);  % ÿ�ֵ��׵�����ʱ���������ݼ��Ͼ���
for j=1:9
    k_position=1;  % ��ֵ��protein_mix_all�����λ��
    for i=1:5
        protein_mix_all(k_position:k_position+cell_number(i,j)-1,:,j)=mix_all_basalave(1:cell_number(i,j),:,i,j);
        k_position=k_position+cell_number(i,j);
    end
end

% ���ǰ����������rho��w�ĺ�����ϵ�������������rho��w�ĺ�����ϵ��ǰ����������rho�ͺ����������rho�ĺ�����ϵ
% ���ȶԴ���ϲ���(��Ӧx���)��������Ȼ��ȡ����ƽ�������㻬����׼�������ƽ��ֵ�ͻ�����׼����ֵ��ֵ��ͬһ���ȣ������Ϻ�����ϵ
parameter_x=[7 10 7];  % �������Ϻ�����ϵ��x�����λ��
parameter_y=[8 11 10];  % �������Ϻ�����ϵ��y�����λ��
win=1;  % ����ƽ�����ڴ�С
protein_all_MA_interp=zeros(2000,9,2,3);  % ÿ�ֵ��׵Ļ���ƽ��ֵ��ֵ���󣬵�һάΪ��������ƽ��ֵ�Ĳ�ֵ������ڶ�άΪ��ͬ�ĵ��ף�����άΪx���y�ᣬ����άΪ�������Ϻ�����ϵ��y�����
std_protein_all_MA_interp=zeros(2000,9,2,3);  % ÿ�ֵ��׵Ļ�����׼���ֵ���󣬵�һάΪ����������׼��Ĳ�ֵ������ڶ�άΪ��ͬ�ĵ��ף�����άΪx���y�ᣬ����άΪ�������Ϻ�����ϵ��y�����
protein_all_fit=zeros(2000,9,3);  % ÿ�ֵ��׵Ķ���ʽ��Ϻ�������һάΪ��Ϻ���(y����x�ģ�����w����rho�ĺ���)���ڶ�άΪ��ͬ�ĵ��ף�����άΪ�������Ϻ�����ϵ
protein_all_delta=zeros(2000,9,3);  % ÿ�ֵ��׵Ķ���ʽ��Ϻ������������䣬��һάΪ��Ϻ���(y����x�ģ�����w����rho�ĺ���)���ڶ�άΪ��ͬ�ĵ��ף�����άΪ�������Ϻ�����ϵ

for k=1:3
    % ð������(���ݴ���Ϻ�����ϵ��x�����)
    for j=1:9  % 9�ֲ�ͬ����
        for ii=2:sum(cell_number(:,j))
            for jj=sum(cell_number(:,j)):-1:ii
                if protein_mix_all(jj,parameter_x(k),j)<protein_mix_all(jj-1,parameter_x(k),j)  % ��ǰһ��ϸ���Ĳ������ں�һ��ϸ���Ĳ����������˳��
                    temp=protein_mix_all(jj,:,j);
                    protein_mix_all(jj,:,j)=protein_mix_all(jj-1,:,j);
                    protein_mix_all(jj-1,:,j)=temp;
                end
            end
        end
    end
    
    % ����������ȡ����ƽ��(���˸�Ƶ��)
    % ���㻬��ƽ��ֵ�ͻ�����׼��
    protein_all_MA=zeros(sum(cell_number(:,i)),9,2);  % ÿ�ֵ��׵Ļ���ƽ��ֵ����һ��Ϊx���Ӧ����ϲ������ڶ���Ϊy���Ӧ����ϲ���
    std_protein_all_MA=zeros(sum(cell_number(:,i)),9,2);  % ÿ�ֵ��׵Ļ�����׼���һ�д���x������Ļ���ƽ��ֵ���ڶ��д���y������Ļ�����׼��
    for i=1:9  % 9�ֲ�ͬ����
        for j=1:sum(cell_number(:,i))
            if j<=(win-1)/2  % ��tΪ��ʼ����ʱ��
                protein_all_MA(j,i,1)=sum(protein_mix_all(1:2*j,parameter_x(k),i))/(2*j);
                protein_all_MA(j,i,2)=sum(protein_mix_all(1:2*j,parameter_y(k),i))/(2*j);
            elseif j<sum(cell_number(:,i))-(win-1)/2+1  % ��tΪ�м�ʱ��
                protein_all_MA(j,i,1)=sum(protein_mix_all(j-(win-1)/2:j+(win-1)/2,parameter_x(k),i))/win;
                protein_all_MA(j,i,2)=sum(protein_mix_all(j-(win-1)/2:j+(win-1)/2,parameter_y(k),i))/win;
            else  % ��tΪĩβ����ʱ��
                temp=sum(cell_number(:,i))-j+1;
                protein_all_MA(j,i,1)=sum(protein_mix_all(sum(cell_number(:,i))-2*temp-1:sum(cell_number(:,i)),parameter_x(k),i))/(2*temp+2);
                protein_all_MA(j,i,2)=sum(protein_mix_all(sum(cell_number(:,i))-2*temp-1:sum(cell_number(:,i)),parameter_y(k),i))/(2*temp+2);
            end
        end
    end
    for i=1:9  % 9�ֲ�ͬ����
        for j=1:sum(cell_number(:,i))
            if j<=(win-1)/2  % ��tΪ��ʼ����ʱ��
                std_protein_all_MA(j,i,1)=protein_all_MA(j,i,1);
                std_protein_all_MA(j,i,2)=std(protein_all_MA(1:2*j,i,2));
            elseif j<sum(cell_number(:,i))-(win-1)/2+1  % ��tΪ�м�ʱ��
                std_protein_all_MA(j,i,1)=protein_all_MA(j,i,1);
                std_protein_all_MA(j,i,2)=std(protein_all_MA(j-(win-1)/2:j+(win-1)/2,i,2));
            else  % ��tΪĩβ����ʱ��
                temp=sum(cell_number(:,i))-j+1;
                std_protein_all_MA(j,i,1)=protein_all_MA(j,i,1);
                std_protein_all_MA(j,i,2)=std(protein_all_MA(sum(cell_number(:,i))-2*temp-1:sum(cell_number(:,i)),i,2));
            end
        end
    end
    
    % ������׼���ֵ
    for i=1:9
        if sum(cell_number(:,i))-win+1<=0
            continue
        end
        
        x_interp=linspace(min(protein_all_MA(1:1:sum(cell_number(:,i)),i,1)),max(protein_all_MA(1:1:sum(cell_number(:,i)),i,1)),2000);
        protein_all_MA_interp(:,i,1,k)=interp1(protein_all_MA(1:1:sum(cell_number(:,i)),i,1),protein_all_MA(1:1:sum(cell_number(:,i)),i,1),x_interp,'linear');  % ��ϲ�������ƽ����׼���ֵ
        protein_all_MA_interp(:,i,2,k)=interp1(protein_all_MA(1:1:sum(cell_number(:,i)),i,1),protein_all_MA(1:1:sum(cell_number(:,i)),i,2),x_interp,'linear');  % ��ϲ�������ƽ����׼���ֵ
        std_protein_all_MA_interp(:,i,1,k)=interp1(std_protein_all_MA(1:1:sum(cell_number(:,i)),i,1),std_protein_all_MA(1:1:sum(cell_number(:,i)),i,1),x_interp,'linear');  % ��ϲ�������ƽ����׼���ֵ
        std_protein_all_MA_interp(:,i,2,k)=interp1(std_protein_all_MA(1:1:sum(cell_number(:,i)),i,1),std_protein_all_MA(1:1:sum(cell_number(:,i)),i,2),x_interp,'linear');  % ��ϲ�������ƽ����׼���ֵ
    end
    
    % ����ʽ�������������ĺ�����ϵ
    N=2;  % ��Ͻ���
    for i=1:9
        if sum(cell_number(:,i))-win+1<=0
            continue
        end
        
        [p,S]=polyfit(protein_all_MA(1:1:sum(cell_number(:,i)),i,1),protein_all_MA(1:1:sum(cell_number(:,i)),i,2),N);
        [yfit,delta]=polyconf(p,protein_all_MA(1:1:sum(cell_number(:,i)),i,1),S,'alpha',0.05,'predopt','curve');  % ��Ϻ�����95%��������
        x_fit=linspace(min(protein_all_MA(1:1:sum(cell_number(:,i)),i,1)),max(protein_all_MA(1:1:sum(cell_number(:,i)),i,1)),2000);  % ��Ϻ�����x����
        delta_interp=interp1(protein_all_MA(1:1:sum(cell_number(:,i)),i,1),delta,x_fit,'linear');  % ����������������ֵ
        protein_all_delta(:,i,k)=delta_interp;
        for ii=1:N+1
            protein_all_fit(:,i,k)=protein_all_fit(:,i,k)+p(ii).*x_fit'.^(N-ii+1);
        end
    end
end

f=fittype('a*x');
fit1=fit(protein_all_MA_interp(:,1,2,3),protein_all_MA_interp(:,1,1,3),f);
conf1=predint(fit1,protein_all_MA_interp(:,1,2,3),0.99,'functional','off');
C_AP=fit1.a;  % ������������ϵ�C_AP


% i=1;k=2;
% % ��ͼ
% for ii=1:2000
%     plot(protein_all_MA_interp(ii,i,1,k),protein_all_MA_interp(ii,i,2,k),'.','MarkerSize',20,'color','r')  % ����ֵ��
%     hold on
% end
% 
% % ������ƽ����׼��
% for ii=1:2:2000
%     aaa=plot([std_protein_all_MA_interp(ii,i,1,k) std_protein_all_MA_interp(ii,i,1,k)],[protein_all_fit(ii,i,k)-std_protein_all_MA_interp(ii,i,2,k) protein_all_fit(ii,i,k)+std_protein_all_MA_interp(ii,i,2,k)],'-','Color','b','linewidth',1);  % ��TD���bar����
%     aaa.Color(4) = 0.2;
%     hold on
% end
% plot(std_protein_all_MA_interp(:,i,1,k),protein_all_fit(:,i,k),'-b','linewidth',3)



%% Step 2
for k=1:9  % 1��ʾlifeact��2��ʾactin��3��ʾMLC��4��ʾRhoA��5��ʾAnillin RBD��6��ʾVD1,7��ʾtalinA��8��ʾCdc42WT��9��ʾCdc42D118A
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
            number_changdubi=0;  % ��������ļ�ʱָʾ���ȱ��ļ��е�����λ��
            
            framepath2=[framepath num2str(k) filename_fanxiegang filename_month filename_day filename_fanxiegang];  % ���ļ���·������ 'C:\Users\Desktop\1\0928'
            framepath_filename_changdubi=[framepath2 filename_changdubi filename_xlsx];
            framepath_filename_Radius=[framepath2 filename_Radius filename_xlsx];
            framepath_filename_ContactAngle=[framepath2 filename_ContactAngle filename_xlsx];

            if exist(framepath_filename_changdubi,'file')   % ���ȱ��ļ��Ƿ���ڿ����ڿ���ɸѡ�Ƿ���ڸ�������õ��׵����ݣ����̲���ʱ��
                data_lumenradius=xlsread(framepath_filename_Radius,'sheet1','A1:A100');  % ����lumen�뾶
                data_lumenradius=sqrt(data_lumenradius./pi);
                data_lumen_ContactAngle=xlsread(framepath_filename_ContactAngle,'sheet1','A1:A100');  % ����lumen�Ӵ���
                data_lumenlength=xlsread(framepath_filename_changdubi,'sheet1','A1:A100');  % ������ļ��µ�lumen����
                data_celllength=xlsread(framepath_filename_changdubi,'sheet1','B1:B100');  % ������ļ��µ�cell����
                % ��Lumen length����cell length����cell length=Lumen length+0.5
                for i=1:length(data_lumenradius)
                    if 2*data_lumenradius(i)*sin(data_lumen_ContactAngle(i)*pi./180)>data_celllength(i)
                        data_celllength(i)=2*data_lumenradius(i)*sin(data_lumen_ContactAngle(i)*pi./180)+0.5;
                    end
                end
                data_changdubi=2.*data_lumenradius.*sin(data_lumen_ContactAngle.*pi./180)./data_celllength;
                
                for i=1:30  % ���ѭ�����ڲ���Ŀ���ļ��Ƿ����(i��ʾ�ڼ���photo��j��ʾphoto�еĵڼ���ϸ����k��ʾӫ��ͨ��)
                    filename1=num2str(i);
                    %%% ȷ��A-P�Ἣ�Է���
                    %%% ÿ����Ƭȷ��һ��ͳһ��A-P�Ἣ�Է���
                    %%% ���裺1. �����ҳ�������Ƭ�µ�����ϸ�����ƣ�����ÿ��ϸ�������ƺ�ԭʼ���ݷֱ𱣴���һ��������
                    %%%       2. ���Ҹ�����Ƭ������������ϸ����ϸ������
                    %%%         (2.1)���У������������ϸ����ϸ�����ݵ�A-P�Ἣ�Է���ȷ��������Ƭ�ļ��Է���(���ж��ϸ��
                    %%%            �Ҽ��Է���ͬ������ݴ����ϸ���ļ��Է���ȷ��������Ƭ�ļ��Է���)
                    %%%         (2.2)���ޣ��ֱ������ּ��Է���ѡ�������Ƭ������ϸ�������ּ��Է���������ô󲿷�
                    %%%            ϸ����ʵ�ʼ��Է�����ͬ��ƽ��A-Pǿ�ȱȸ���ļ��Է���
                    filenamearray_PD=cell(10,1);  % PD��polarity determination
                    y_danyi=zeros(1500,10);  % �洢��һ�����ļ���ϸ������
                    y_part1=zeros(1500,10);  % �洢�ж�������ļ���ϸ�����ݵĵ�һ����
                    y_part2=zeros(1500,10);  % �洢�ж�������ļ���ϸ�����ݵĵڶ�����
                    y_part3=zeros(1500,10);  % �洢�ж�������ļ���ϸ�����ݵĵ�������
                    number_monocell=0;  % ��ʾ�е�һ�����ļ���ϸ������
                    number_multicell=0;  % ��ʾ�ж�������ļ���ϸ������
                    polarity_direction=0;  % ���ڴ���A-P�Ἣ�Է���(1��ʾ����2��ʾ����)
                    
                    %%% 1. �����ҳ�������Ƭ�µ�����ϸ�����ƣ�����ÿ��ϸ�������ƺ�ԭʼ���ݷֱ𱣴���һ��������
                    %%% ��������Ƭ�µ�����ϸ�����ݷֱ𱣴浽y_danyi�����y_part1/2/3�����У��ļ������浽filenamearray_PD������
                    for j=1:10
                        filename2=num2str(j);
                        filename=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3];  % �� 0928_111
                        framepath_filename=[framepath2 filename];
                        
                        % ƴ���������ļ�·��
                        framepath_filename_1=[framepath_filename filename_xiahuaxian filename_1 filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111_1.csv
                        framepath_filename_2=[framepath_filename filename_xiahuaxian filename_2 filename_csv];
                        framepath_filename_3=[framepath_filename filename_xiahuaxian filename_3 filename_csv];
                        framepath_filename_danyi=[framepath_filename filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111.csv
                        
                        % ϸ�����ݺ��ļ�������ֵ
                        if exist(framepath_filename_danyi,'file')   % һ��ϸ����������һ���ļ�
                            number_monocell=number_monocell+1;
                            filenamearray_PD{number_monocell+number_multicell}=filename;  % �ļ�����ֵ
                            
                            data=xlsread(framepath_filename_danyi); % ��������
                            [X,Y]=size(data); % �������ݾ���ߴ�
                            for ii=1:X  % ���ݾ��󱣴浽y_danyi������
                                y_danyi(ii,number_monocell)=data(ii,2);
                            end
                        elseif exist(framepath_filename_1,'file')   % һ��ϸ���������ж���ļ�
                            number_multicell=number_multicell+1;
                            filenamearray_PD{number_monocell+number_multicell}=filename;  % �ļ�����ֵ
                            
                            data_1=xlsread(framepath_filename_1); % ��������
                            data_2=xlsread(framepath_filename_2);
                            data_3=xlsread(framepath_filename_3);
                            for ii=1:length(data_1)  % ���ݾ��󱣴浽y_part1/2/3������
                                y_part1(ii,number_multicell)=data_1(ii,2);
                            end
                            for ii=1:length(data_2)  % ���ݾ��󱣴浽y_part1/2/3������
                                y_part2(ii,number_multicell)=data_2(ii,2);
                            end
                            for ii=1:length(data_3)  % ���ݾ��󱣴浽y_part1/2/3������
                                y_part3(ii,number_multicell)=data_3(ii,2);
                            end
                        end
                    end
                    
                    %%% ȷ��A-P�Ἣ�Է���
                    if number_multicell==0  % ������Ƭ�²�����������ϸ����ϸ�����ݣ�ֱ�ӽ����ôμ��㣬������һ�μ���
                        number_changdubi=number_changdubi+number_monocell;
                        continue
                    else
                        if number_monocell>=1
                         %%% 2. ���Ҹ�����Ƭ������������ϸ����ϸ������
                         %%%    (2.1)���У������������ϸ����ϸ�����ݵ�A-P�Ἣ�Է���ȷ��������Ƭ�ļ��Է���(���ж��ϸ��
                         %%%         �Ҽ��Է���ͬ������ݴ����ϸ���ļ��Է���ȷ��������Ƭ�ļ��Է���)
                            AP_overactivity_ratio_PD=zeros(number_monocell,1);  % �������������ϸ����ϸ����A-P������ǿ�ȱ�
                            
                            % ȷ��������ϸ����ϸ�����Է���
                            % ����ͨ��lateral domain����������˹�ֲ���ϻ��A-P���Դ�С
                            for i_momocell=1:number_monocell
                                for X=length(y_danyi(:,i_momocell)):-1:1  % ȷ��y���ݾ��󳤶�
                                    if y_danyi(X,i_momocell)~=0
                                        break
                                    end
                                end
                                X=X+1;
                                % X���������ԭ��Գ�
                                x=linspace(-0.5,0.5,X);
                                % Y�����ݹ�һ��
                                y_danyi(:,i_momocell)=y_danyi(:,i_momocell)./max(y_danyi(:,i_momocell));
                                
                                Y_anterior=zeros(1,round((X-1)*0.2));  % ����ϸ��������������Ҷ�ֵ���ݣ����ڵ���˹�ֲ����ȷ�����Է���
                                Y_posterior=zeros(1,round((X-1)*0.2));
                                
                                % Y_anterior Y_posterior����ֵ
                                for ii=1:round((X-1)*0.2)
                                    Y_anterior(ii)=y_danyi(ii,i_momocell);
                                end
                                for ii=X-round((X-1)*0.2)+1:X
                                    Y_posterior(ii-X+round((X-1)*0.2))=y_danyi(ii,i_momocell);
                                end
                                
                                % ��С������������˹�ֲ���ȷ��Ҫɾȥ�Ĳ���
                                % ���Y_anterior
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
                                
                                % ���Y_posterior
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
                            
                            % ͳ���������Է����ϸ����
                            forward_cellnumber=0;  % anterior lateral domain������ǿ��ϸ����
                            backward_cellnumber=0;  % posterior lateral domain������ǿ��ϸ����
                            for ii=1:number_monocell
                                if AP_overactivity_ratio_PD(ii)>1
                                    forward_cellnumber=forward_cellnumber+1;
                                else
                                    backward_cellnumber=backward_cellnumber+1;
                                end
                            end
                            
                            % �жϸ�����Ƭ��A-P�Ἣ�Է���
                            if forward_cellnumber==number_monocell  % ����������ϸ����ϸ����Ϊ����
                                polarity_direction=1;
                            elseif backward_cellnumber==number_monocell  % ����������ϸ����ϸ����Ϊ����
                                polarity_direction=2;
                            elseif forward_cellnumber~=0 && forward_cellnumber~=number_monocell  % ������ϸ����ϸ�������ַ���
                                if forward_cellnumber>=backward_cellnumber
                                    polarity_direction=1;
                                else
                                    polarity_direction=2;
                                end
                            end
                        else
                            %%% 2. ���Ҹ�����Ƭ������������ϸ����ϸ������
                            %%%    (2.2)���ޣ��ֱ������ּ��Է���ѡ�������Ƭ������ϸ�������ּ��Է���������ô󲿷�
                            %%%         ϸ����ʵ�ʼ��Է�����ͬ��ƽ��A-Pǿ�ȱȸ���ļ��Է���
                            AP_overactivity_ratio_PD=zeros(number_multicell,2);  % �������������ϸ����ϸ��(������ϸ��)��A-P����ǿ�ȱ�
                            
                            %%% ��ð����ּ��Է�����ԭʼ���ݼ���õ�������ϸ��A-P��ǿ�ȱ�
                            AP_overactivity_ratio_test=zeros(2);  % ��һ������Ϊ�����ԣ��ڶ�������Ϊ������
                            k_sum=0;
                            for i_sum=1:5  % ������Ϊ5��ʱ��(��ʵ����ֵ)��ƽ�����Դ�С
                                AP_overactivity_ratio_test(1)=AP_overactivity_ratio_test(1)+AP_overactivity_ratio(i_sum,k);
                                if AP_overactivity_ratio(i_sum,k)~=0
                                    k_sum=k_sum+1;
                                end
                            end
                            AP_overactivity_ratio_test(1)=AP_overactivity_ratio_test(1)/k_sum;
                            AP_overactivity_ratio_test(2)=1/AP_overactivity_ratio_test(1);     % ������Ϊ�����Եĵ���
                            
                            for i_test=1:2
                                % ȷ������ϸ�������ļ��ĳ��ȣ�������X_part1/2/3�����У��������ݹ�һ��
                                X_part1=zeros(number_multicell);
                                X_part2=zeros(number_multicell);
                                X_part3=zeros(number_multicell);
                                for i_multicell=1:number_multicell
                                    for i_X_part1=length(y_part1(:,i_multicell)):-1:1  % ȷ��y���ݾ��󳤶�
                                        if y_part1(i_X_part1,i_multicell)~=0
                                            break
                                        end
                                    end
                                    X_part1(i_multicell)=i_X_part1;
                                    for i_X_part2=length(y_part2(:,i_multicell)):-1:1  % ȷ��y���ݾ��󳤶�
                                        if y_part1(i_X_part2,i_multicell)~=0
                                            break
                                        end
                                    end
                                    X_part2(i_multicell)=i_X_part2;
                                    for i_X_part3=length(y_part3(:,i_multicell)):-1:1  % ȷ��y���ݾ��󳤶�
                                        if y_part1(i_X_part3,i_multicell)~=0
                                            break
                                        end
                                    end
                                    X_part3(i_multicell)=i_X_part3;
                                end
                                    
                                % Y�����ݹ�һ��(ͬ����y_part1/2/3����ƽ��ֵ)
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
                                
                                % ϸ���Ҷ�ֵ���ݵ�Ԥ����
                                for i_multicell=1:number_multicell
                                    % X���������ԭ��Գ�
                                    x_part1=linspace(-0.5,0.5,X_part1(i_multicell));
                                    x_part2=linspace(-0.5,0.5,X_part2(i_multicell));
                                    x_part3=linspace(-0.5,0.5,X_part3(i_multicell));
                                    
                                    Y_anterior=zeros(1,X_part1(i_multicell));  % ����ϸ�������������������Ҷ�ֵ���ݣ����ڵ���˹�ֲ����
                                    Y_middle=zeros(1,X_part2(i_multicell));
                                    Y_posterior=zeros(1,X_part3(i_multicell));
                                    
                                    % ������һ��ϸ���Ƿ����(����һ��ϸ���ĵ������ļ��͸�ϸ���ĵ�һ���ļ���ͬ����˵����ϸ�����ڣ���y_part1���д���)
                                    if i_multicell>=2  % ��ʾ������һ��ϸ��
                                        if X_part3(i_multicell-1)==X_part1(i_multicell) && ~sum(abs(y_part3(:,i_multicell-1)-y_part1(:,i_multicell)))
                                            total_data_previous2=0;
                                            total_data_2=0;
                                            for i_sum=1:X_part2(i_multicell-1)  % ������һ��ϸ��basal domain���ܻҶ�ֵ
                                                total_data_previous2=total_data_previous2+y_part2(i_sum,i_multicell-1);
                                            end
                                            total_data_previous2=total_data_previous2/X_part2(i_multicell-1);  % ��һ��ϸ����basal domain���ܻҶ�ֵȡƽ��
                                            for i_sum=1:X_part2(i_multicell)  % �����ϸ��basal domain���ܻҶ�ֵ
                                                total_data_2=total_data_2+y_part2(i_sum,i_multicell);
                                            end
                                            total_data_2=total_data_2/X_part2(i_multicell);  % ��ϸ����basal domain���ܻҶ�ֵȡƽ��
                                            for i_ratio=1:X_part1(i_multicell)  % ��ϸ��anterior(��)��lateral domain�Ҷ�ֵ���ݵĲ�� (��ʽ��y_part1=y_part1*[ave(y_part2)*A-P-ratio/(ave(y_part2_before)+ave(y_part2)*A-P-ratio)]
                                                Y_anterior(i_ratio)=y_part1(i_ratio,i_multicell)*total_data_2*AP_overactivity_ratio_test(i_test)/(total_data_previous2+total_data_2*AP_overactivity_ratio_test(i_test));
                                            end
                                            
                                            % ����ϸ����һ�������ļ������ݽ��еߵ�
                                            Y_anterior_fuben=Y_anterior;
                                            for i_symmetry=1:length(Y_anterior)
                                                Y_anterior(i_symmetry)=Y_anterior_fuben(length(Y_anterior)-i_symmetry+1);
                                            end
                                        else  % ��������һ��ϸ����ֱ�ӽ�ԭʼy_part1��ֵ��ֵ��Y_anterior
                                            for i_ratio=1:X_part1(i_multicell)
                                                Y_anterior(i_ratio)=y_part1(i_ratio,i_multicell);
                                            end
                                        end
                                    else  % ������Ƭ��һ��ϸ����ֱ�ӽ�ԭʼy_part1��ֵ��ֵ��Y_anterior
                                        for i_ratio=1:X_part1(i_multicell)
                                            Y_anterior(i_ratio)=y_part1(i_ratio,i_multicell);
                                        end
                                    end
                                    
                                    % ������һ��ϸ���Ƿ����(����һ��ϸ���ĵ�һ���ļ��͸�ϸ���ĵ������ļ���ͬ����˵����ϸ�����ڣ���y_part3���д���)
                                    if i_multicell<=number_multicell-1  % ��ʾ������һ��ϸ��
                                        if X_part3(i_multicell)==X_part1(i_multicell+1) && ~sum(abs(y_part3(:,i_multicell)-y_part1(:,i_multicell+1)))
                                            total_data_next2=0;
                                            total_data_2=0;
                                            for i_sum=1:X_part2(i_multicell+1)  % ������һ��ϸ��basal domain���ܻҶ�ֵ
                                                total_data_next2=total_data_next2+y_part2(i_sum,i_multicell+1);
                                            end
                                            total_data_next2=total_data_next2/X_part2(i_multicell+1);  % ��һ��ϸ��basal domain���ܻҶ�ֵȡƽ��
                                            for i_sum=1:X_part2(i_multicell)  % �����ϸ��basal domain���ܻҶ�ֵ
                                                total_data_2=total_data_2+y_part2(i_sum,i_multicell);
                                            end
                                            total_data_2=total_data_2/X_part2(i_multicell);  % ��ϸ��basal domain���ܻҶ�ֵȡƽ��
                                            for i_ratio=1:X_part1(i_multicell)  % ��ϸ��posterior(��)��lateral domain�Ҷ�ֵ���ݵĲ�� (��ʽ��y_part3=y_part3*[ave(y_part2)/(ave(y_part2)+ave(y_part2_after)*A-P-ratio)]
                                                Y_posterior(i_ratio)=y_part3(i_ratio,i_multicell)*total_data_2/(total_data_2+total_data_next2*AP_overactivity_ratio_test(i_test));
                                            end
                                        else  % ��������һ��ϸ����ֱ�ӽ�ԭʼy_part3��ֵ��ֵ��Y_posterior
                                            for i_ratio=1:X_part1(i_multicell)
                                                Y_posterior(i_ratio)=y_part3(i_ratio,i_multicell);
                                            end
                                        end
                                    else  % ������Ƭ���һ��ϸ����ֱ�ӽ�ԭʼy_part3��ֵ��ֵ��Y_posterior
                                        for i_ratio=1:X_part1(i_multicell)
                                            Y_posterior(i_ratio)=y_part3(i_ratio,i_multicell);
                                        end
                                    end
                                    
                                    % ��С������������˹�ֲ������������������ǿ�ȣ�����A-Pǿ�ȱ�
                                    % ���Y_anterior
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
                                    
                                    % ���Y_posterior
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
                                    if i_test==1  % A-P����ǿ�ȱ�Ϊrho_inf_anterior_best/rho_inf_posterior_best
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
                                    elseif i_test==2  % A-P����ǿ�ȱ�Ϊrho_inf_posterior_best/rho_inf_anterior_best
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
                            
                            % ͳ���������Է����ϸ����
                            forward_cellnumber=0;  % �����Բ�ֺ󣬺ͼ��輫�Է�����ͬ(anterior lateral domain��������ǿ)��ϸ����
                            backward_cellnumber=0;  % �����Բ�ֺ󣬺ͼ��輫�Է�����ͬ(posterior lateral domain��������ǿ)��ϸ����
                            for ii=1:number_multicell
                                if AP_overactivity_ratio_PD(ii,1)>1
                                    forward_cellnumber=forward_cellnumber+1;
                                end
                                if AP_overactivity_ratio_PD(ii,2)>1
                                    backward_cellnumber=backward_cellnumber+1;
                                end
                            end
                            
                            %%% �жϸ�����Ƭ��A-P�Ἣ�Է���
                            if forward_cellnumber==number_multicell  % �����Բ�ֺ󣬺ͼ��輫�Է�����ͬ(anterior lateral domain��������ǿ)��ϸ����������ϸ����
                                polarity_direction=1;
                            elseif backward_cellnumber==number_multicell  % �����Բ�ֺ󣬺ͼ��輫�Է�����ͬ(posterior lateral domain��������ǿ)��ϸ����������ϸ����
                                polarity_direction=2;
                            elseif forward_cellnumber~=number_multicell  % ϸ�����Է�����ͳһ����
                                % ���� ������ü��Բ�ַ�����ͬ��ϸ����+����ϸ��A-P ratio�Ļ� ��Ϊ���б�׼
                                forward_standard=AP_overactivity_ratio_PD(1,1);  % ���������б�׼��������ֵ
                                backward_standard=AP_overactivity_ratio_PD(1,2);  % ���������б�׼��������ֵ
                                for i_multicell=2:number_multicell  % ��������ϸ��A-P ratio�Ļ�
                                    forward_standard=forward_standard*AP_overactivity_ratio_PD(i_multicell,1);
                                    backward_standard=backward_standard*AP_overactivity_ratio_PD(i_multicell,2);
                                end
                                forward_standard=forward_standard+forward_cellnumber;  % ���ϼ�����ü��Բ�ַ�����ͬ��ϸ����
                                backward_standard=backward_standard+backward_cellnumber;
                                if forward_standard>=backward_standard
                                    polarity_direction=1;
                                else
                                    polarity_direction=2;
                                end
                            end
                        end
                    end
                    
                    
                    %%% ϸ���������������(�����ļ�Ԥ����ƴ�ӡ��������ݵ�ɾ�������ݹ�һ��������˹�ֲ����)
                    for j=1:10
                        filename2=num2str(j);
                        filename=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3];  % �� 0928_111
                        framepath_filename=[framepath2 filename];
                        
                        % ƴ���������ļ�·����������һ��exist�����Ĳ���
                        framepath_filename_1=[framepath_filename filename_xiahuaxian filename_1 filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111_1.csv
                        framepath_filename_2=[framepath_filename filename_xiahuaxian filename_2 filename_csv];
                        framepath_filename_3=[framepath_filename filename_xiahuaxian filename_3 filename_csv];
                        framepath_filename_danyi=[framepath_filename filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111.csv
                        
                        if exist(framepath_filename_1,'file')   % һ��ϸ���������ж���ļ�
                            number_changdubi=number_changdubi+1;
                            for aa=1
                                %%% ȷ����ϸ����A-P��ǿ�ȱ�
                                A_P_ratio2=0;  % �����ϸ����A-P�Ἣ��ǿ�ȱ�
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
                                
                                %%% ÿ��ϸ���Ķ�������ļ��ϲ���һ��������
                                data_1=xlsread(framepath_filename_1); % ��������
                                data_2=xlsread(framepath_filename_2);
                                data_3=xlsread(framepath_filename_3);
                                
                                % ����ϸ����anterior-lateral basal posterior-lateral domain�ĳ���
                                anterior_lateral_length=data_1(end,1);
                                posterior_lateral_length=data_3(end,1);
                                basal_length=data_2(end,1);
                                
                                %%% ϸ��ԭʼ���ݲ��
                                protein_number=1;
                                [data_1,data_2,data_3] = CellOverlapDataSplitting(data_1,data_2,data_3,C_AP,protein_number,data_changdubi,number_changdubi,j);
                                data_split=[data_3_split1;data_3_split2];

                                % ��������ɵ�һ��ϸ�����������ݺϲ�
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

                                %%% ��ԭʼϸ���������Ҷ�ֵ���ݽ�������������ݵ�ɾ������
                                % ���ݲ����������basal-lateral domain����˹�ֲ�
                                % ���ݵ�����
                                [X,Y]=size(data); % �������ݾ���ߴ�
                                if mod(X,2)==0   % ������Ϊż��������������Բ�ֵ�����������ƽ�ƹ����еĶԳ���
                                    data2=zeros(X+1,Y);
                                    data2(:,1)=imresize(data(:,1),[X+1,1], 'bilinear'); % ��X��������Բ�ֵ
                                    data2(:,2)=imresize(data(:,2),[X+1,1], 'bilinear'); % ��Y��������Բ�ֵ
                                else
                                    data2=data;
                                end

                                % X��ƽ��������ԭ��ԳƵ����䣬����һ��
                                bias=0.5; % ƽ����
                                x=data2(:,1)/max(data2(:,1)); % ��x����й�һ��
                                x=x-bias; % x��ƽ������������ԭ��Գ�
                                % Y�����ݹ�һ��
                                y=data2(:,2)/max(data2(:,2));


                                Y_anterior=zeros(1,round((X-1)*0.2));  % ����ϸ��������������Ҷ�ֵ���ݣ����ڵ���˹�ֲ����
                                Y_posterior=zeros(1,round((X-1)*0.2));
                                Y_anterior_AP=zeros(1,round((X-1)*0.2));  % ����ϸ����һ����������������Ҷ�ֵ���ݣ����ڵ���˹�ֲ����
                                Y_posterior_AP=zeros(1,round((X-1)*0.2));

                                if data_changdubi(number_changdubi)<=0.2  % lumen�γɳ�����������������δ��������ɾ������
                                    for k_AP=1:10
                                        if k_AP~=k
                                            filename3_AP=num2str(k_AP);  % AP: anothor protein
                                            filename_AP=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3_AP];  % �� 0928_111
                                            framepath_filename_AP=[framepath2 filename_AP];
                                            framepath_filename_1_AP=[framepath_filename_AP filename_xiahuaxian filename_1 filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111_1.csv
                                            framepath_filename_2_AP=[framepath_filename_AP filename_xiahuaxian filename_2 filename_csv];
                                            framepath_filename_3_AP=[framepath_filename_AP filename_xiahuaxian filename_3 filename_csv];
                                            if exist(framepath_filename_1_AP,'file')
                                                break
                                            end
                                        end
                                    end
                                    if k_AP~=10  % ��ʾ��ϸ������һ����������
                                        %%% ÿ��ϸ���Ķ�������ļ��ϲ���һ��������
                                        data_1_AP=xlsread(framepath_filename_1_AP); % ��������
                                        data_2_AP=xlsread(framepath_filename_2_AP);
                                        data_3_AP=xlsread(framepath_filename_3_AP);

                                        %%% ϸ��ԭʼ���ݲ��
                                        protein_number=2;
                                        [data_1_AP,data_2_AP,data_3_AP] = CellOverlapDataSplitting(data_1_AP,data_2_AP,data_3_AP,C_AP,protein_number,data_changdubi,number_changdubi,j);
                                        data_AP_split=[data_3_AP_split1;data_3_AP_split2];

                                        % ��������ɵ�һ��ϸ�����������ݺϲ�
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

                                        % ���ݵ�����
                                        [X_AP,Y_AP]=size(data_AP); % �������ݾ���ߴ�
                                        if mod(X_AP,2)==0   % ������Ϊż��������������Բ�ֵ�����������ƽ�ƹ����еĶԳ���
                                            data2_AP=zeros(X_AP+1,Y_AP);
                                            data2_AP(:,1)=imresize(data_AP(:,1),[X_AP+1,1], 'bilinear'); % ��X��������Բ�ֵ
                                            data2_AP(:,2)=imresize(data_AP(:,2),[X_AP+1,1], 'bilinear'); % ��Y��������Բ�ֵ
                                        else
                                            data2_AP=data_AP;
                                        end

                                        % X��ƽ��������ԭ��ԳƵ����䣬����һ��
                                        bias=0.5; % ƽ����
                                        x_AP=data2_AP(:,1)/max(data2_AP(:,1)); % ��x����й�һ��
                                        x_AP=x_AP-bias; % x��ƽ������������ԭ��Գ�
                                        x=x_AP;
                                        % Y�����ݹ�һ��
                                        y=data2(:,2)/max(data2(:,2));
                                        y_AP=data2_AP(:,2)/max(data2_AP(:,2));
                                    end
                                    x_min=-0.5;x_max=0.5;
                                else  % lumen�γ��к������������������ѽ�������ɾ������
                                    for k_AP=1:10
                                        if k_AP~=k
                                            filename3_AP=num2str(k_AP);  % AP: anothor protein
                                            filename_AP=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3_AP];  % �� 0928_111
                                            framepath_filename_AP=[framepath2 filename_AP];
                                            framepath_filename_1_AP=[framepath_filename_AP filename_xiahuaxian filename_1 filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111_1.csv
                                            framepath_filename_2_AP=[framepath_filename_AP filename_xiahuaxian filename_2 filename_csv];
                                            framepath_filename_3_AP=[framepath_filename_AP filename_xiahuaxian filename_3 filename_csv];
                                            if exist(framepath_filename_1_AP,'file')
                                                break
                                            end
                                        end
                                    end
                                    if k_AP~=10  % ��ʾ��ϸ������һ����������
                                        %%% ÿ��ϸ���Ķ�������ļ��ϲ���һ��������
                                        data_1_AP=xlsread(framepath_filename_1_AP); % ��������
                                        data_2_AP=xlsread(framepath_filename_2_AP);
                                        data_3_AP=xlsread(framepath_filename_3_AP);

                                        %%% ϸ��ԭʼ���ݲ��
                                        protein_number=2;
                                        [data_1_AP,data_2_AP,data_3_AP] = CellOverlapDataSplitting(data_1_AP,data_2_AP,data_3_AP,C_AP,protein_number,data_changdubi,number_changdubi,j);
                                        data_AP_split=[data_3_AP_split1;data_3_AP_split2];

                                        % ��������ɵ�һ��ϸ�����������ݺϲ�
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
                                        % ���ݵ�����
                                        [X_AP,Y_AP]=size(data_AP); % �������ݾ���ߴ�
                                        if mod(X_AP,2)==0   % ������Ϊż��������������Բ�ֵ�����������ƽ�ƹ����еĶԳ���
                                            data2_AP=zeros(X_AP+1,Y_AP);
                                            data2_AP(:,1)=imresize(data_AP(:,1),[X_AP+1,1], 'bilinear'); % ��X��������Բ�ֵ
                                            data2_AP(:,2)=imresize(data_AP(:,2),[X_AP+1,1], 'bilinear'); % ��Y��������Բ�ֵ
                                        else
                                            data2_AP=data_AP;
                                        end

                                        % X��ƽ��������ԭ��ԳƵ����䣬����һ��
                                        bias=0.5; % ƽ����
                                        x_AP=data2_AP(:,1)/max(data2_AP(:,1)); % ��x����й�һ��
                                        x_AP=x_AP-bias; % x��ƽ������������ԭ��Գ�
                                        x=x_AP;
                                        % Y�����ݹ�һ��
                                        y=data2(:,2)/max(data2(:,2));
                                        y_AP=data2_AP(:,2)/max(data2_AP(:,2));

                                        % Levenberg-Marquardt method��������˹�ֲ���ȷ��Ҫɾȥ�Ĳ���
                                        % �õ����������
                                        for ab=1
                                            if data_changdubi(number_changdubi)<=0.2  % lumen�γ����ڣ�lateral domain������δ����������9�������������
                                                E_lateral_start=0.3;
                                                E_lateral_end=0.5;
                                            else  % lumen�γ��к��ڣ�lateral domain����������������8�������������
                                                E_lateral_start=0.5;
                                                E_lateral_end=0.5;
                                            end

                                            wucha_best=100000000;
                                            % ʹ��Levenberg-Marquardt method���
                                            for i_fit=1:5  % ѭ��10�Σ�Ѱ��������Ͻ�
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
                                        % ��һ�������������
                                        for ab=1
                                            if data_changdubi(number_changdubi)<=0.2  % lumen�γ����ڣ�lateral domain������δ����������9�������������
                                                E_lateral_start=0.3;
                                                E_lateral_end=0.5;
                                            else  % lumen�γ��к��ڣ�lateral domain����������������8�������������
                                                E_lateral_start=0.5;
                                                E_lateral_end=0.5;
                                            end

                                            wucha_best=100000000;
                                            % ʹ��Levenberg-Marquardt method���
                                            for i_fit=1:5  % ѭ��10�Σ�Ѱ��������Ͻ�
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
                                        x_min=min(E_anterior_best,E_anterior_AP_best);  % ǰ��߽�ȡ�������׵����ǰ�߽��н�С��
                                        x_max=max(E_posterior_best,E_posterior_AP_best);  % ���߽�ȡ�������׵���Ϻ�߽��нϴ��
                                    else  % kk=9��ʾ��ϸ������һ�ֵ��׵������ļ�
                                        % ��ȡ�õ��׵�����
                                        % ���ݵ�����

                                        % Levenberg-Marquardt method��������˹�ֲ���ȷ��Ҫɾȥ�Ĳ���
                                        for ab=1
                                            if data_changdubi(number_changdubi)<=0.2  % lumen�γ����ڣ�lateral domain������δ����������9�������������
                                                E_lateral_start=0.3;
                                                E_lateral_end=0.5;
                                            else  % lumen�γ��к��ڣ�lateral domain����������������8�������������
                                                E_lateral_start=0.5;
                                                E_lateral_end=0.5;
                                            end

                                            wucha_best=100000000;
                                            % ʹ��Levenberg-Marquardt method���
                                            for i_fit=1:5  % ѭ��10�Σ�Ѱ��������Ͻ�
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

                                % ��x_min��x_maxƫ��߽����(������ϸ�����ȵ�1/10)��˵���������Ч�����ã���ȡԭʼ�߽�ֵ
                                if (x_min+0.5)/0.5>0.1
                                    x_min=-0.5;
                                elseif (0.5-x_max)/0.5>0.1
                                    x_max=0.5;
                                end
                                
                                %%% ��ɾ��������ɵ�ϸ���Ҷ�ֵ���ݽ�������˹�ֲ����
                                % �ҵ�ԭX���Ϻ�x_min x_max��ӽ��ĵ�
                                wucha_best=100;
                                for ii=1:X
                                    wucha=abs(x(ii)-x_min);
                                    if wucha<wucha_best
                                        wucha_best=wucha;
                                        x_min_best=x(ii);  % ����ӽ�x_min��ԭX����ֵ��ʱ������x_min_best������
                                    end
                                end
                                x_min=x_min_best;  % x_min�ĳ�X������ӽ�x_min���ֵ
                                wucha_best=100;
                                for ii=1:X
                                    wucha=abs(x(ii)-x_max);
                                    if wucha<wucha_best
                                        wucha_best=wucha;
                                        x_max_best=x(ii);  % ����ӽ�x_max��ԭX����ֵ��ʱ������x_max_best������
                                    end
                                end
                                x_max=x_max_best;  % x_max�ĳ�X������ӽ�x_max���ֵ
                                
                                % ɾȥX��Ͷ�Ӧ�Ҷ�ֵ(y)��С��x_min�ʹ���x_max�Ĳ���
                                x_delete=zeros(round((X-1)*(x_max-x_min)+1),1);
                                y_delete=zeros(round((X-1)*(x_max-x_min)+1),1);
                                for ii=1:round((X-1)*(x_max-x_min)+1)
                                    x_delete(ii)=x(ii+round((x_min+0.5)*(X-1)));
                                    y_delete(ii)=y(ii+round((x_min+0.5)*(X-1)));
                                end
                                x=x_delete;
                                y=y_delete;
                                % X������ƽ��������ԭ��ԳƵ����䣬��������(-0.5,0.5)
                                bias=x_max+x_min; % ƽ����
                                stretch=1/(x_max-x_min); % ������
                                x=x-bias/2; % x��ƽ������������ԭ��Գ�
                                x=x*stretch;
                                
                                %%% ʹ��Levenberg-Marquardt method���
                                max_y=0;min_y=100;
                                for iy=floor(length(y)*0.25):floor(length(y)*0.75)
                                    if y(iy)>max_y
                                        max_y=y(iy);
                                    elseif y(iy)<min_y
                                        min_y=y(iy);
                                    end
                                end
                                if max_y/min_y<1.5  % ��ʾ������������Ѿ�������ʧ��lumenʱ�ڴ��ں��ڣ�Ϊ��ֹ��ϳ�������������λ�ù̶�������
                                    E_bottom=-0.05;E_top=0.05;
                                else  % ��������ǿ����λ�ò��̶�
                                    E_bottom=-0.1;E_top=0.1;
                                end

                                % ���˹���
                                % ����MultiGaussian_fit�Ӻ�������Multi-Gaussian distribution
                                cell_type=2;
                                [rho0_best,rho_inf_array,E_array,w_array,break_point,wucha_best,TriGaussian_parameter]=MultiGaussian_fit(x,y,data_changdubi(number_changdubi),anterior_lateral_length,basal_length,posterior_lateral_length,cell_type);
                                framepath_filename
                                % ����MultiGaussianTransfertoEquivalentTriGaussian�Ӻ��������ЧTri-Gaussian distribution
                                [EquivalentTriGaussian_parameter,area_sumsum]=MultiGaussianTransfertoEquivalentTriGaussian(x,y,rho0_best,rho_inf_array,E_array,w_array,break_point);
                            end
                            
                            % �ж�ϸ���Ƿ������Ϊanterior part���Ҳ�Ϊposterior part�������ǣ�����ϲ�����������ߡ�ԭʼ���ݹ���ԭ�����Գƴ���
                            if rho_inf_array(end-1)<rho_inf_array(end)
                                % MultiGaussian�������Գƴ���
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

                                % EquivalentTriGaussian�������Գƴ���
                                anterior_posterior_ring_reverse=EquivalentTriGaussian_parameter(5:6);
                                EquivalentTriGaussian_parameter(5:6)=EquivalentTriGaussian_parameter(8:9);
                                EquivalentTriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                EquivalentTriGaussian_parameter(4)=-EquivalentTriGaussian_parameter(4);
                                
                                % TriGaussian�������Գƴ���
                                anterior_posterior_ring_reverse=TriGaussian_parameter(5:6);
                                TriGaussian_parameter(5:6)=TriGaussian_parameter(8:9);
                                TriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                TriGaussian_parameter(4)=-TriGaussian_parameter(4);

                                y_reverse=zeros(length(y),1);  % ԭʼ���ݷ�ת
                                for ii=1:length(y)
                                    y_reverse(length(y)-ii+1)=y(ii);
                                end
                                y=y_reverse;
                            end

                            % ӫ��ǿ��ֵ(������)��һ��
                            % �������ɵĶ��˹�ֲ��ڶ������ϵ���ѧ�����任Ϊ1��ͬʱ��ԭʼ�������ݺ͵�Ч����˹�ֲ����ȱ�������
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
                            y_MultiGaussian_fit=ones(size(x))*rho0_best;  % ���˹�������
                            for i_Gaussian=1:length(rho_inf_array)
                                y_MultiGaussian_fit=y_MultiGaussian_fit+rho_inf_array(i_Gaussian)*exp(-1/2*((x-E_array(i_Gaussian))/w_array(i_Gaussian)).^2);
                            end
                            y_EquivalentTriGaussian_fit=ones(size(x))*EquivalentTriGaussian_parameter(1);  % ��Ч����˹�������
                            for i_Gaussian=1:3
                                y_EquivalentTriGaussian_fit=y_EquivalentTriGaussian_fit+EquivalentTriGaussian_parameter(3*i_Gaussian-1)*exp(-1/2*((x-EquivalentTriGaussian_parameter(3*i_Gaussian+1))/EquivalentTriGaussian_parameter(3*i_Gaussian)).^2);
                            end
                            y=y./a0;

                            %%% ��ͼ
                            filename_plot=[filename_photo filename_kongge filename1 filename_kongge filename_cell filename_kongge filename2 filename_kongge filename_protein];

                            figure
                            plot(x,y,'.','MarkerSize',15,'Color','k');  % �������ݵ�
                            hold on;
                            plot(x,y_MultiGaussian_fit,'b','linewidth',3);  % ������ϵĸ�˹����
                            plot(x,y_EquivalentTriGaussian_fit,'r','linewidth',3);  % ������ϵĸ�˹����
                            hold off;
                            axis ([-0.5 0.5 0 max(y)+0.1]);
                            title(filename_plot)
                            xlabel('x')
                            ylabel('cortex thickness')
                            set(gca,'YTick',0:0.2:max(y)+0.1);  % �޸�y����̶�
                            set(gca,'XTick',-0.5:0.1:0.5);  % �޸�y����̶�
                            legend('Experiments','Multi-Gaussian Curve','Tri-Gaussian Curve')
                            set(gca,'FontName','Arial','FontSize',12)
                            saveas(gcf,[framepath2,filename_plot,'.jpg']);
                            close

                            %%% ��ϲ�������
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
                            
                            %%% ��ϸ�ϸ������һ�����������ļ�
                            if k_AP~=10  % ��ϸ��������һ�����������ļ�
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

                                %%% ��ɾ��������ɵ�ϸ���Ҷ�ֵ���ݽ�������˹�ֲ����
                                % ɾȥX��Ͷ�Ӧ�Ҷ�ֵ(y)��С��x_min�ʹ���x_max�Ĳ���
                                x_delete_AP=zeros(round((X-1)*(x_max-x_min)+1),1);
                                y_delete_AP=zeros(round((X-1)*(x_max-x_min)+1),1);
                                for ii=1:round((X-1)*(x_max-x_min)+1)
                                    x_delete_AP(ii)=x_AP(ii+round((x_min+0.5)*(X-1)));
                                    y_delete_AP(ii)=y_AP(ii+round((x_min+0.5)*(X-1)));
                                end
                                x_AP=x_delete_AP;
                                y_AP=y_delete_AP;
                                % X������ƽ��������ԭ��ԳƵ����䣬��������(-0.5,0.5)
                                bias=x_max+x_min; % ƽ����
                                stretch=1/(x_max-x_min); % ������
                                x_AP=x_AP-bias/2; % x��ƽ������������ԭ��Գ�
                                x_AP=x_AP*stretch;

                                %%% ʹ��Levenberg-Marquardt method���
                                max_y=0;min_y=100;
                                for iy=floor(length(y_AP)*0.25):floor(length(y_AP)*0.75)
                                    if y_AP(iy)>max_y
                                        max_y=y_AP(iy);
                                    elseif y_AP(iy)<min_y
                                        min_y=y_AP(iy);
                                    end
                                end
                                if max_y/min_y<1.5  % ��ʾ������������Ѿ�������ʧ��lumenʱ�ڴ��ں��ڣ�Ϊ��ֹ��ϳ�������������λ�ù̶�������
                                    E_bottom=-0.05;E_top=0.05;
                                else  % ��������ǿ����λ�ò��̶�
                                    E_bottom=-0.1;E_top=0.1;
                                end

                                % ���˹���
                                % ����MultiGaussian_fit�Ӻ�������Multi-Gaussian distribution
                                cell_type=2;
                                [rho0_best,rho_inf_array,E_array,w_array,break_point,wucha_best,TriGaussian_parameter]=MultiGaussian_fit(x_AP,y_AP,data_changdubi(number_changdubi),anterior_lateral_length,basal_length,posterior_lateral_length,cell_type);
                                framepath_filename_AP
                                % ����MultiGaussianTransfertoEquivalentTriGaussian�Ӻ��������ЧTri-Gaussian distribution
                                [EquivalentTriGaussian_parameter,area_sumsum]=MultiGaussianTransfertoEquivalentTriGaussian(x_AP,y_AP,rho0_best,rho_inf_array,E_array,w_array,break_point);

                                % �ж�ϸ���Ƿ������Ϊanterior part���Ҳ�Ϊposterior part�������ǣ�����ϲ�����������ߡ�ԭʼ���ݹ���ԭ�����Գƴ���
                                if rho_inf_array(end-1)<rho_inf_array(end)
                                    % MultiGaussian�������Գƴ���
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

                                    % EquivalentTriGaussian�������Գƴ���
                                    anterior_posterior_ring_reverse=EquivalentTriGaussian_parameter(5:6);
                                    EquivalentTriGaussian_parameter(5:6)=EquivalentTriGaussian_parameter(8:9);
                                    EquivalentTriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                    EquivalentTriGaussian_parameter(4)=-EquivalentTriGaussian_parameter(4);

                                    % TriGaussian�������Գƴ���
                                    anterior_posterior_ring_reverse=TriGaussian_parameter(5:6);
                                    TriGaussian_parameter(5:6)=TriGaussian_parameter(8:9);
                                    TriGaussian_parameter(8:9)=anterior_posterior_ring_reverse;
                                    TriGaussian_parameter(4)=-TriGaussian_parameter(4);

                                    y_reverse_AP=zeros(length(y_AP),1);  % ԭʼ���ݷ�ת
                                    for ii=1:length(y_AP)
                                        y_reverse_AP(length(y_AP)-ii+1)=y_AP(ii);
                                    end
                                    y_AP=y_reverse_AP;
                                end

                                % ӫ��ǿ��ֵ(������)��һ��
                                % �������ɵĶ��˹�ֲ��ڶ������ϵ���ѧ�����任Ϊ1��ͬʱ��ԭʼ�������ݺ͵�Ч����˹�ֲ����ȱ�������
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
                                y_MultiGaussian_fit_AP=ones(size(x_AP))*rho0_best;  % ���˹�������
                                for i_Gaussian=1:length(rho_inf_array)
                                    y_MultiGaussian_fit_AP=y_MultiGaussian_fit_AP+rho_inf_array(i_Gaussian)*exp(-1/2*((x_AP-E_array(i_Gaussian))/w_array(i_Gaussian)).^2);
                                end
                                y_EquivalentTriGaussian_fit_AP=ones(size(x_AP))*EquivalentTriGaussian_parameter(1);  % ��Ч����˹�������
                                for i_Gaussian=1:3
                                    y_EquivalentTriGaussian_fit_AP=y_EquivalentTriGaussian_fit_AP+EquivalentTriGaussian_parameter(3*i_Gaussian-1)*exp(-1/2*((x_AP-EquivalentTriGaussian_parameter(3*i_Gaussian+1))/EquivalentTriGaussian_parameter(3*i_Gaussian)).^2);
                                end
                                y_AP=y_AP./a0;

                                %%% ��ͼ
                                filename_plot=[filename_photo filename_kongge filename1 filename_kongge filename_cell filename_kongge filename2 filename_kongge filename_protein];

                                figure
                                plot(x_AP,y_AP,'.','MarkerSize',15,'Color','k');  % �������ݵ�
                                hold on;
                                plot(x_AP,y_MultiGaussian_fit_AP,'b','linewidth',3);  % ������ϵĸ�˹����
                                plot(x_AP,y_EquivalentTriGaussian_fit_AP,'r','linewidth',3);  % ������ϵĸ�˹����
                                hold off;
                                axis ([-0.5 0.5 0 max(y_AP)+0.1]);
                                title(filename_plot)
                                xlabel('x')
                                ylabel('cortex thickness')
                                set(gca,'YTick',0:0.2:max(y_AP)+0.1);  % �޸�y����̶�
                                set(gca,'XTick',-0.5:0.1:0.5);  % �޸�y����̶�
                                legend('Experiments','Multi-Gaussian Curve','Tri-Gaussian Curve')
                                set(gca,'FontName','Arial','FontSize',14)
                                saveas(gcf,[framepath2,filename_plot,'.jpg']);
                                close

                                %%% ��ϲ�������
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
                            elseif k_AP==10  % kk=9��ʾ��ϸ������һ�ֵ��׵������ļ�
                            end
                            
                        elseif exist(framepath_filename_danyi,'file')
                            number_changdubi=number_changdubi+1;  % ������ϸ���ĳ��ȱ�����
                        else  % ��������Ƭ�£������ڸ�ϸ���������ļ���˵��������Ƭ�е�ϸ���Ѿ������ɣ�����j(ϸ��)ѭ��
                            break
                        end
                        
                        % ���Ǹ�ϸ�������һ�������õ�kkѭ������ֵ
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
                        filename=[filename_month filename_day filename_xiahuaxian filename1 filename2 filename3];  % �� 0928_111
                        framepath_filename=[framepath2 filename];
                        framepath_filename_1=[framepath_filename filename_xiahuaxian filename_1 filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111_1.csv
                        framepath_filename_2=[framepath_filename filename_xiahuaxian filename_2 filename_csv];
                        framepath_filename_3=[framepath_filename filename_xiahuaxian filename_3 filename_csv];
                        framepath_filename_danyi=[framepath_filename filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111.csv
                    end
                    if j==1  % ��������Ƭ�²����ڵ�һ��ϸ������˵����������������Ƭ�Ѿ�������ɣ�����i(��Ƭ)ѭ��
                        break
                    end
                end
            end
        end
    end
end



