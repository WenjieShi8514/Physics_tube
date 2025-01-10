%%% Relative lateral edge signal calculation

%%% Steps:
%%% 1. Read the mix_all data file (xlsx) that stores the fitted Tri-Gaussian distribution characteristic parameters 
%%%    and the protein file (e.g., lifeact) (xlsx) that stores file names corresponding to the characteristic parameters.
%%% 2. Exclude cells with excessive fitting errors or cells where equatorial axis contractile ring (or lateral 
%%%    contractile ring) characteristic parameters significantly deviate from the average.
%%% 3. Calculate the average thickness and average length of the lateral domain cortex before lumen formation 
%%%    (lumen length/cell length < 0.08).
%%%    (The average thickness is calculated using both cells with three data files and cells with a single data 
%%%    file, while the average length is calculated only from cells with a single data file.)
%%% 4. Calculate the cortex thickness and length of the anterior- and posterior-lateral domains of each cell at 
%%%    different stages after lumen formation.
%%%    (1) Calculate the cortex thickness and length of the anterior- and posterior-lateral domains of cells with 
%%%        three data files and fit the functional relationship between anterior/posterior-lateral domain length/lateral-basal-lateral 
%%%        whole length and lumen length/cell length.
%%%    (2) Calculate the cortex thickness of the anterior- and posterior-lateral domains of cells with a single data
%%%        file. The division of the three cell parts is based on the lateral-basal-lateral whole length and lumen length/cell length.
%%% 5. Calculate the contractile force exerted on each lumen by the contractile rings of notochord cells on both sides.
%%%    (1) First, calculate the contractile force exerted on lumens with grayscale data from cells on both sides. 
%%%        The contractile force is equal to the sum of the product of the posterior-lateral domain cortex thickness, length, 
%%%        and unit thickness-generated pressure from the previous cell and the anterior-lateral domain cortex thickness, length, 
%%%        and unit thickness-generated pressure from the next cell, minus the contractile force corresponding to the basal thickness.
%%%    (2) Then, calculate the contractile force exerted on lumens with grayscale data from only one side. The 
%%%        contractile force is equal to the sum of the product of the anterior-lateral domain cortex thickness, length, 
%%%        and unit thickness-generated pressure and the posterior-lateral domain cortex thickness, length, and unit 
%%%        thickness-generated pressure from the same cell, minus the contractile force corresponding to the basal thickness.
%%%    (3) Fit the contractile force exerted on lumens by the contractile rings of notochord cells on both sides,
%%%        sorted by lumen length/cell length.
%%%    (4) Fit the relationship between normalized TD (lumen length/cell length) and time (t) using observational
%%%        data. Fit the functional relationship between the contractile force exerted on lumens by the contractile rings 
%%%        of notochord cells on both sides and time (t).



%% step 1
framepath='C:\Users\';

filename_xiahuaxian='_';
filename_fanxiegang='\';
filename_1='1';
filename_2='2';
filename_3='3';
filename_xlsx='.xlsx';
filename_csv='.csv';

filename_mix_MultiGaussian_lifeact='mix_lifeact_MultiGaussian';
filename_mix_MultiGaussian_actin='mix_actin_MultiGaussian';
filename_mix_MultiGaussian_MLC='mix_MLC_MultiGaussian';
filename_mix_MultiGaussian_RhoA='mix_RhoA_MultiGaussian';
filename_mix_MultiGaussian_RBD='mix_RBD_MultiGaussian';
filename_mix_MultiGaussian_VD1='mix_VD1_MultiGaussian';
filename_mix_MultiGaussian_talinA='mix_talinA_MultiGaussian';
filename_mix_MultiGaussian_Cdc42WT='mix_Cdc42WT_MultiGaussian';
filename_mix_MultiGaussian_Cdc42D118A='mix_Cdc42D118A_MultiGaussian';

filename_mix_EquivalentTriGaussian_lifeact='mix_lifeact_EquivalentTriGaussian';
filename_mix_EquivalentTriGaussian_actin='mix_actin_EquivalentTriGaussian';
filename_mix_EquivalentTriGaussian_MLC='mix_MLC_EquivalentTriGaussian';
filename_mix_EquivalentTriGaussian_RhoA='mix_RhoA_EquivalentTriGaussian';
filename_mix_EquivalentTriGaussian_RBD='mix_RBD_EquivalentTriGaussian';
filename_mix_EquivalentTriGaussian_VD1='mix_VD1_EquivalentTriGaussian';
filename_mix_EquivalentTriGaussian_talinA='mix_talinA_EquivalentTriGaussian';
filename_mix_EquivalentTriGaussian_Cdc42WT='mix_Cdc42WT_EquivalentTriGaussian';
filename_mix_EquivalentTriGaussian_Cdc42D118A='mix_Cdc42D118A_EquivalentTriGaussian';

filename_mix_TriGaussian_lifeact='mix_lifeact_TriGaussian';
filename_mix_TriGaussian_actin='mix_actin_TriGaussian';
filename_mix_TriGaussian_MLC='mix_MLC_TriGaussian';
filename_mix_TriGaussian_RhoA='mix_RhoA_TriGaussian';
filename_mix_TriGaussian_RBD='mix_RBD_TriGaussian';
filename_mix_TriGaussian_VD1='mix_VD1_TriGaussian';
filename_mix_TriGaussian_talinA='mix_talinA_TriGaussian';
filename_mix_TriGaussian_Cdc42WT='mix_Cdc42WT_TriGaussian';
filename_mix_TriGaussian_Cdc42D118A='mix_Cdc42D118A_TriGaussian';

filename_lifeact='lifeact';
filename_actin='actin';
filename_MLC='MLC';
filename_RhoA='RhoA';
filename_RBD='RBD';
filename_VD1='VD1';
filename_talinA='talinA';
filename_Cdc42WT='Cdc42WT';
filename_Cdc42D118A='Cdc42D118A';

filename_changdubi='长度比';
filename_Radius='半径 (Radius, 测量结果为lumen面积)';
filename_ContactAngle='接触角 (Contact Angle, 测量结果为圆周角)';

framepath_filename_mix_MultiGaussian_lifeact=[framepath filename_mix_MultiGaussian_lifeact filename_xlsx];
framepath_filename_mix_MultiGaussian_actin=[framepath filename_mix_MultiGaussian_actin filename_xlsx];
framepath_filename_mix_MultiGaussian_MLC=[framepath filename_mix_MultiGaussian_MLC filename_xlsx];
framepath_filename_mix_MultiGaussian_RhoA=[framepath filename_mix_MultiGaussian_RhoA filename_xlsx];
framepath_filename_mix_MultiGaussian_RBD=[framepath filename_mix_MultiGaussian_RBD filename_xlsx];
framepath_filename_mix_MultiGaussian_VD1=[framepath filename_mix_MultiGaussian_VD1 filename_xlsx];
framepath_filename_mix_MultiGaussian_talinA=[framepath filename_mix_MultiGaussian_talinA filename_xlsx];
framepath_filename_mix_MultiGaussian_Cdc42WT=[framepath filename_mix_MultiGaussian_Cdc42WT filename_xlsx];
framepath_filename_mix_MultiGaussian_Cdc42D118A=[framepath filename_mix_MultiGaussian_Cdc42D118A filename_xlsx];
framepath_filename_mix_EquivalentTriGaussian_lifeact=[framepath filename_mix_EquivalentTriGaussian_lifeact filename_xlsx];
framepath_filename_mix_EquivalentTriGaussian_actin=[framepath filename_mix_EquivalentTriGaussian_actin filename_xlsx];
framepath_filename_mix_EquivalentTriGaussian_MLC=[framepath filename_mix_EquivalentTriGaussian_MLC filename_xlsx];
framepath_filename_mix_EquivalentTriGaussian_RhoA=[framepath filename_mix_EquivalentTriGaussian_RhoA filename_xlsx];
framepath_filename_mix_EquivalentTriGaussian_RBD=[framepath filename_mix_EquivalentTriGaussian_RBD filename_xlsx];
framepath_filename_mix_EquivalentTriGaussian_VD1=[framepath filename_mix_EquivalentTriGaussian_VD1 filename_xlsx];
framepath_filename_mix_EquivalentTriGaussian_talinA=[framepath filename_mix_EquivalentTriGaussian_talinA filename_xlsx];
framepath_filename_mix_EquivalentTriGaussian_Cdc42WT=[framepath filename_mix_EquivalentTriGaussian_Cdc42WT filename_xlsx];
framepath_filename_mix_EquivalentTriGaussian_Cdc42D118A=[framepath filename_mix_EquivalentTriGaussian_Cdc42D118A filename_xlsx];
framepath_filename_mix_TriGaussian_lifeact=[framepath filename_mix_TriGaussian_lifeact filename_xlsx];
framepath_filename_mix_TriGaussian_actin=[framepath filename_mix_TriGaussian_actin filename_xlsx];
framepath_filename_mix_TriGaussian_MLC=[framepath filename_mix_TriGaussian_MLC filename_xlsx];
framepath_filename_mix_TriGaussian_RhoA=[framepath filename_mix_TriGaussian_RhoA filename_xlsx];
framepath_filename_mix_TriGaussian_RBD=[framepath filename_mix_TriGaussian_RBD filename_xlsx];
framepath_filename_mix_TriGaussian_VD1=[framepath filename_mix_TriGaussian_VD1 filename_xlsx];
framepath_filename_mix_TriGaussian_talinA=[framepath filename_mix_TriGaussian_talinA filename_xlsx];
framepath_filename_mix_TriGaussian_Cdc42WT=[framepath filename_mix_TriGaussian_Cdc42WT filename_xlsx];
framepath_filename_mix_TriGaussian_Cdc42D118A=[framepath filename_mix_TriGaussian_Cdc42D118A filename_xlsx];



%%% 读取存储拟合三(多)高斯分布特征参数的mix_all数据文件
for a=1
    % 读取mix_lifeact变量矩阵
    mix_MultiGaussian_lifeact_stage1=xlsread(framepath_filename_mix_MultiGaussian_lifeact,'sheet1','A1:Y1000');
    mix_MultiGaussian_lifeact_stage2=xlsread(framepath_filename_mix_MultiGaussian_lifeact,'sheet2','A1:Y1000');
    mix_MultiGaussian_lifeact_stage3=xlsread(framepath_filename_mix_MultiGaussian_lifeact,'sheet3','A1:Y1000');
    mix_MultiGaussian_lifeact_stage4=xlsread(framepath_filename_mix_MultiGaussian_lifeact,'sheet4','A1:Y1000');
    mix_MultiGaussian_lifeact_stage5=xlsread(framepath_filename_mix_MultiGaussian_lifeact,'sheet5','A1:Y1000');
    mix_EquivalentTriGaussian_lifeact_stage1=xlsread(framepath_filename_mix_EquivalentTriGaussian_lifeact,'sheet1','A1:M1000');
    mix_EquivalentTriGaussian_lifeact_stage2=xlsread(framepath_filename_mix_EquivalentTriGaussian_lifeact,'sheet2','A1:M1000');
    mix_EquivalentTriGaussian_lifeact_stage3=xlsread(framepath_filename_mix_EquivalentTriGaussian_lifeact,'sheet3','A1:M1000');
    mix_EquivalentTriGaussian_lifeact_stage4=xlsread(framepath_filename_mix_EquivalentTriGaussian_lifeact,'sheet4','A1:M1000');
    mix_EquivalentTriGaussian_lifeact_stage5=xlsread(framepath_filename_mix_EquivalentTriGaussian_lifeact,'sheet5','A1:M1000');
    mix_TriGaussian_lifeact_stage1=xlsread(framepath_filename_mix_TriGaussian_lifeact,'sheet1','A1:M1000');
    mix_TriGaussian_lifeact_stage2=xlsread(framepath_filename_mix_TriGaussian_lifeact,'sheet2','A1:M1000');
    mix_TriGaussian_lifeact_stage3=xlsread(framepath_filename_mix_TriGaussian_lifeact,'sheet3','A1:M1000');
    mix_TriGaussian_lifeact_stage4=xlsread(framepath_filename_mix_TriGaussian_lifeact,'sheet4','A1:M1000');
    mix_TriGaussian_lifeact_stage5=xlsread(framepath_filename_mix_TriGaussian_lifeact,'sheet5','A1:M1000');
    
    for aa=1
        for k1=1:1000
            if mix_MultiGaussian_lifeact_stage1(k1,1)==0
                k1_end=k1-1;
                break
            end
        end
        for k2=1:1000
            if mix_MultiGaussian_lifeact_stage2(k2,1)==0
                k2_end=k2-1;
                break
            end
        end
        if mix_MultiGaussian_lifeact_stage2(1000,1)~=0
            k2_end=1000;
        end
        for k3=1:1000
            if mix_MultiGaussian_lifeact_stage3(k3,1)==0
                k3_end=k3-1;
                break
            end
        end
        if mix_MultiGaussian_lifeact_stage3(1000,1)~=0
            k3_end=1000;
        end
        for k4=1:1000
            if mix_MultiGaussian_lifeact_stage4(k4,1)==0
                k4_end=k4-1;
                break
            end
        end
        for k5=1:1000
            if mix_MultiGaussian_lifeact_stage5(k5,1)==0
                k5_end=k5-1;
                break
            end
        end
    end
    mix_MultiGaussian_lifeact=[mix_MultiGaussian_lifeact_stage1(1:k1_end,:);mix_MultiGaussian_lifeact_stage2(1:k2_end,:);mix_MultiGaussian_lifeact_stage3(1:k3_end,:);mix_MultiGaussian_lifeact_stage4(1:k4_end,:);mix_MultiGaussian_lifeact_stage5(1:k5_end,:)];
    mix_EquivalentTriGaussian_lifeact=[mix_EquivalentTriGaussian_lifeact_stage1(1:k1_end,:);mix_EquivalentTriGaussian_lifeact_stage2(1:k2_end,:);mix_EquivalentTriGaussian_lifeact_stage3(1:k3_end,:);mix_EquivalentTriGaussian_lifeact_stage4(1:k4_end,:);mix_EquivalentTriGaussian_lifeact_stage5(1:k5_end,:)];
    mix_TriGaussian_lifeact=[mix_TriGaussian_lifeact_stage1(1:k1_end,:);mix_TriGaussian_lifeact_stage2(1:k2_end,:);mix_TriGaussian_lifeact_stage3(1:k3_end,:);mix_TriGaussian_lifeact_stage4(1:k4_end,:);mix_TriGaussian_lifeact_stage5(1:k5_end,:)];
    
    % 读取mix_actin变量矩阵
    mix_MultiGaussian_actin_stage1=xlsread(framepath_filename_mix_MultiGaussian_actin,'sheet1','A1:Y1000');
    mix_MultiGaussian_actin_stage2=xlsread(framepath_filename_mix_MultiGaussian_actin,'sheet2','A1:Y1000');
    mix_MultiGaussian_actin_stage3=xlsread(framepath_filename_mix_MultiGaussian_actin,'sheet3','A1:Y1000');
    mix_MultiGaussian_actin_stage4=xlsread(framepath_filename_mix_MultiGaussian_actin,'sheet4','A1:Y1000');
    mix_MultiGaussian_actin_stage5=xlsread(framepath_filename_mix_MultiGaussian_actin,'sheet5','A1:Y1000');
    mix_EquivalentTriGaussian_actin_stage1=xlsread(framepath_filename_mix_EquivalentTriGaussian_actin,'sheet1','A1:M1000');
    mix_EquivalentTriGaussian_actin_stage2=xlsread(framepath_filename_mix_EquivalentTriGaussian_actin,'sheet2','A1:M1000');
    mix_EquivalentTriGaussian_actin_stage3=xlsread(framepath_filename_mix_EquivalentTriGaussian_actin,'sheet3','A1:M1000');
    mix_EquivalentTriGaussian_actin_stage4=xlsread(framepath_filename_mix_EquivalentTriGaussian_actin,'sheet4','A1:M1000');
    mix_EquivalentTriGaussian_actin_stage5=xlsread(framepath_filename_mix_EquivalentTriGaussian_actin,'sheet5','A1:M1000');
    mix_TriGaussian_actin_stage1=xlsread(framepath_filename_mix_TriGaussian_actin,'sheet1','A1:M1000');
    mix_TriGaussian_actin_stage2=xlsread(framepath_filename_mix_TriGaussian_actin,'sheet2','A1:M1000');
    mix_TriGaussian_actin_stage3=xlsread(framepath_filename_mix_TriGaussian_actin,'sheet3','A1:M1000');
    mix_TriGaussian_actin_stage4=xlsread(framepath_filename_mix_TriGaussian_actin,'sheet4','A1:M1000');
    mix_TriGaussian_actin_stage5=xlsread(framepath_filename_mix_TriGaussian_actin,'sheet5','A1:M1000');
    
    for aa=1
        for k1=1:1000
            if mix_MultiGaussian_actin_stage1(k1,1)==0
                k1_end=k1-1;
                break
            end
        end
        for k2=1:1000
            if mix_MultiGaussian_actin_stage2(k2,1)==0
                k2_end=k2-1;
                break
            end
        end
        if mix_MultiGaussian_actin_stage2(1000,1)~=0
            k2_end=1000;
        end
        for k3=1:1000
            if mix_MultiGaussian_actin_stage3(k3,1)==0
                k3_end=k3-1;
                break
            end
        end
        if mix_MultiGaussian_actin_stage3(1000,1)~=0
            k3_end=1000;
        end
        for k4=1:1000
            if mix_MultiGaussian_actin_stage4(k4,1)==0
                k4_end=k4-1;
                break
            end
        end
        for k5=1:1000
            if mix_MultiGaussian_actin_stage5(k5,1)==0
                k5_end=k5-1;
                break
            end
        end
    end
    mix_MultiGaussian_actin=[mix_MultiGaussian_actin_stage1(1:k1_end,:);mix_MultiGaussian_actin_stage2(1:k2_end,:);mix_MultiGaussian_actin_stage3(1:k3_end,:);mix_MultiGaussian_actin_stage4(1:k4_end,:);mix_MultiGaussian_actin_stage5(1:k5_end,:)];
    mix_EquivalentTriGaussian_actin=[mix_EquivalentTriGaussian_actin_stage1(1:k1_end,:);mix_EquivalentTriGaussian_actin_stage2(1:k2_end,:);mix_EquivalentTriGaussian_actin_stage3(1:k3_end,:);mix_EquivalentTriGaussian_actin_stage4(1:k4_end,:);mix_EquivalentTriGaussian_actin_stage5(1:k5_end,:)];
    mix_TriGaussian_actin=[mix_TriGaussian_actin_stage1(1:k1_end,:);mix_TriGaussian_actin_stage2(1:k2_end,:);mix_TriGaussian_actin_stage3(1:k3_end,:);mix_TriGaussian_actin_stage4(1:k4_end,:);mix_TriGaussian_actin_stage5(1:k5_end,:)];
    
    % 读取mix_MLC变量矩阵
    mix_MultiGaussian_MLC_stage1=xlsread(framepath_filename_mix_MultiGaussian_MLC,'sheet1','A1:Y1000');
    mix_MultiGaussian_MLC_stage2=xlsread(framepath_filename_mix_MultiGaussian_MLC,'sheet2','A1:Y1000');
    mix_MultiGaussian_MLC_stage3=xlsread(framepath_filename_mix_MultiGaussian_MLC,'sheet3','A1:Y1000');
    mix_MultiGaussian_MLC_stage4=xlsread(framepath_filename_mix_MultiGaussian_MLC,'sheet4','A1:Y1000');
    mix_MultiGaussian_MLC_stage5=xlsread(framepath_filename_mix_MultiGaussian_MLC,'sheet5','A1:Y1000');
    mix_EquivalentTriGaussian_MLC_stage1=xlsread(framepath_filename_mix_EquivalentTriGaussian_MLC,'sheet1','A1:M1000');
    mix_EquivalentTriGaussian_MLC_stage2=xlsread(framepath_filename_mix_EquivalentTriGaussian_MLC,'sheet2','A1:M1000');
    mix_EquivalentTriGaussian_MLC_stage3=xlsread(framepath_filename_mix_EquivalentTriGaussian_MLC,'sheet3','A1:M1000');
    mix_EquivalentTriGaussian_MLC_stage4=xlsread(framepath_filename_mix_EquivalentTriGaussian_MLC,'sheet4','A1:M1000');
    mix_EquivalentTriGaussian_MLC_stage5=xlsread(framepath_filename_mix_EquivalentTriGaussian_MLC,'sheet5','A1:M1000');
    mix_TriGaussian_MLC_stage1=xlsread(framepath_filename_mix_TriGaussian_MLC,'sheet1','A1:M1000');
    mix_TriGaussian_MLC_stage2=xlsread(framepath_filename_mix_TriGaussian_MLC,'sheet2','A1:M1000');
    mix_TriGaussian_MLC_stage3=xlsread(framepath_filename_mix_TriGaussian_MLC,'sheet3','A1:M1000');
    mix_TriGaussian_MLC_stage4=xlsread(framepath_filename_mix_TriGaussian_MLC,'sheet4','A1:M1000');
    mix_TriGaussian_MLC_stage5=xlsread(framepath_filename_mix_TriGaussian_MLC,'sheet5','A1:M1000');
    
    for aa=1
        for k1=1:1000
            if mix_MultiGaussian_MLC_stage1(k1,1)==0
                k1_end=k1-1;
                break
            end
        end
        for k2=1:1000
            if mix_MultiGaussian_MLC_stage2(k2,1)==0
                k2_end=k2-1;
                break
            end
        end
        if mix_MultiGaussian_MLC_stage2(1000,1)~=0
            k2_end=1000;
        end
        for k3=1:1000
            if mix_MultiGaussian_MLC_stage3(k3,1)==0
                k3_end=k3-1;
                break
            end
        end
        if mix_MultiGaussian_MLC_stage3(1000,1)~=0
            k3_end=1000;
        end
        for k4=1:1000
            if mix_MultiGaussian_MLC_stage4(k4,1)==0
                k4_end=k4-1;
                break
            end
        end
        for k5=1:1000
            if mix_MultiGaussian_MLC_stage5(k5,1)==0
                k5_end=k5-1;
                break
            end
        end
    end
    mix_MultiGaussian_MLC=[mix_MultiGaussian_MLC_stage1(1:k1_end,:);mix_MultiGaussian_MLC_stage2(1:k2_end,:);mix_MultiGaussian_MLC_stage3(1:k3_end,:);mix_MultiGaussian_MLC_stage4(1:k4_end,:);mix_MultiGaussian_MLC_stage5(1:k5_end,:)];
    mix_EquivalentTriGaussian_MLC=[mix_EquivalentTriGaussian_MLC_stage1(1:k1_end,:);mix_EquivalentTriGaussian_MLC_stage2(1:k2_end,:);mix_EquivalentTriGaussian_MLC_stage3(1:k3_end,:);mix_EquivalentTriGaussian_MLC_stage4(1:k4_end,:);mix_EquivalentTriGaussian_MLC_stage5(1:k5_end,:)];
    mix_TriGaussian_MLC=[mix_TriGaussian_MLC_stage1(1:k1_end,:);mix_TriGaussian_MLC_stage2(1:k2_end,:);mix_TriGaussian_MLC_stage3(1:k3_end,:);mix_TriGaussian_MLC_stage4(1:k4_end,:);mix_TriGaussian_MLC_stage5(1:k5_end,:)];
    
    % 读取mix_RhoA变量矩阵
    mix_MultiGaussian_RhoA_stage1=xlsread(framepath_filename_mix_MultiGaussian_RhoA,'sheet1','A1:Y1000');
    mix_MultiGaussian_RhoA_stage2=xlsread(framepath_filename_mix_MultiGaussian_RhoA,'sheet2','A1:Y1000');
    mix_MultiGaussian_RhoA_stage3=xlsread(framepath_filename_mix_MultiGaussian_RhoA,'sheet3','A1:Y1000');
    mix_MultiGaussian_RhoA_stage4=xlsread(framepath_filename_mix_MultiGaussian_RhoA,'sheet4','A1:Y1000');
    mix_MultiGaussian_RhoA_stage5=xlsread(framepath_filename_mix_MultiGaussian_RhoA,'sheet5','A1:Y1000');
    mix_EquivalentTriGaussian_RhoA_stage1=xlsread(framepath_filename_mix_EquivalentTriGaussian_RhoA,'sheet1','A1:M1000');
    mix_EquivalentTriGaussian_RhoA_stage2=xlsread(framepath_filename_mix_EquivalentTriGaussian_RhoA,'sheet2','A1:M1000');
    mix_EquivalentTriGaussian_RhoA_stage3=xlsread(framepath_filename_mix_EquivalentTriGaussian_RhoA,'sheet3','A1:M1000');
    mix_EquivalentTriGaussian_RhoA_stage4=xlsread(framepath_filename_mix_EquivalentTriGaussian_RhoA,'sheet4','A1:M1000');
    mix_EquivalentTriGaussian_RhoA_stage5=xlsread(framepath_filename_mix_EquivalentTriGaussian_RhoA,'sheet5','A1:M1000');
    mix_TriGaussian_RhoA_stage1=xlsread(framepath_filename_mix_TriGaussian_RhoA,'sheet1','A1:M1000');
    mix_TriGaussian_RhoA_stage2=xlsread(framepath_filename_mix_TriGaussian_RhoA,'sheet2','A1:M1000');
    mix_TriGaussian_RhoA_stage3=xlsread(framepath_filename_mix_TriGaussian_RhoA,'sheet3','A1:M1000');
    mix_TriGaussian_RhoA_stage4=xlsread(framepath_filename_mix_TriGaussian_RhoA,'sheet4','A1:M1000');
    mix_TriGaussian_RhoA_stage5=xlsread(framepath_filename_mix_TriGaussian_RhoA,'sheet5','A1:M1000');
    
    for aa=1
        for k1=1:1000
            if mix_MultiGaussian_RhoA_stage1(k1,1)==0
                k1_end=k1-1;
                break
            end
        end
        for k2=1:1000
            if mix_MultiGaussian_RhoA_stage2(k2,1)==0
                k2_end=k2-1;
                break
            end
        end
        if mix_MultiGaussian_RhoA_stage2(1000,1)~=0
            k2_end=1000;
        end
        for k3=1:1000
            if mix_MultiGaussian_RhoA_stage3(k3,1)==0
                k3_end=k3-1;
                break
            end
        end
        if mix_MultiGaussian_RhoA_stage3(1000,1)~=0
            k3_end=1000;
        end
        for k4=1:1000
            if mix_MultiGaussian_RhoA_stage4(k4,1)==0
                k4_end=k4-1;
                break
            end
        end
        for k5=1:1000
            if mix_MultiGaussian_RhoA_stage5(k5,1)==0
                k5_end=k5-1;
                break
            end
        end
    end
    mix_MultiGaussian_RhoA=[mix_MultiGaussian_RhoA_stage1(1:k1_end,:);mix_MultiGaussian_RhoA_stage2(1:k2_end,:);mix_MultiGaussian_RhoA_stage3(1:k3_end,:);mix_MultiGaussian_RhoA_stage4(1:k4_end,:);mix_MultiGaussian_RhoA_stage5(1:k5_end,:)];
    mix_EquivalentTriGaussian_RhoA=[mix_EquivalentTriGaussian_RhoA_stage1(1:k1_end,:);mix_EquivalentTriGaussian_RhoA_stage2(1:k2_end,:);mix_EquivalentTriGaussian_RhoA_stage3(1:k3_end,:);mix_EquivalentTriGaussian_RhoA_stage4(1:k4_end,:);mix_EquivalentTriGaussian_RhoA_stage5(1:k5_end,:)];
    mix_TriGaussian_RhoA=[mix_TriGaussian_RhoA_stage1(1:k1_end,:);mix_TriGaussian_RhoA_stage2(1:k2_end,:);mix_TriGaussian_RhoA_stage3(1:k3_end,:);mix_TriGaussian_RhoA_stage4(1:k4_end,:);mix_TriGaussian_RhoA_stage5(1:k5_end,:)];
    
    % 读取mix_RBD变量矩阵
    mix_MultiGaussian_RBD_stage1=xlsread(framepath_filename_mix_MultiGaussian_RBD,'sheet1','A1:Y1000');
    mix_MultiGaussian_RBD_stage2=xlsread(framepath_filename_mix_MultiGaussian_RBD,'sheet2','A1:Y1000');
    mix_MultiGaussian_RBD_stage3=xlsread(framepath_filename_mix_MultiGaussian_RBD,'sheet3','A1:Y1000');
    mix_MultiGaussian_RBD_stage4=xlsread(framepath_filename_mix_MultiGaussian_RBD,'sheet4','A1:Y1000');
    mix_MultiGaussian_RBD_stage5=xlsread(framepath_filename_mix_MultiGaussian_RBD,'sheet5','A1:Y1000');
    mix_EquivalentTriGaussian_RBD_stage1=xlsread(framepath_filename_mix_EquivalentTriGaussian_RBD,'sheet1','A1:M1000');
    mix_EquivalentTriGaussian_RBD_stage2=xlsread(framepath_filename_mix_EquivalentTriGaussian_RBD,'sheet2','A1:M1000');
    mix_EquivalentTriGaussian_RBD_stage3=xlsread(framepath_filename_mix_EquivalentTriGaussian_RBD,'sheet3','A1:M1000');
    mix_EquivalentTriGaussian_RBD_stage4=xlsread(framepath_filename_mix_EquivalentTriGaussian_RBD,'sheet4','A1:M1000');
    mix_EquivalentTriGaussian_RBD_stage5=xlsread(framepath_filename_mix_EquivalentTriGaussian_RBD,'sheet5','A1:M1000');
    mix_TriGaussian_RBD_stage1=xlsread(framepath_filename_mix_TriGaussian_RBD,'sheet1','A1:M1000');
    mix_TriGaussian_RBD_stage2=xlsread(framepath_filename_mix_TriGaussian_RBD,'sheet2','A1:M1000');
    mix_TriGaussian_RBD_stage3=xlsread(framepath_filename_mix_TriGaussian_RBD,'sheet3','A1:M1000');
    mix_TriGaussian_RBD_stage4=xlsread(framepath_filename_mix_TriGaussian_RBD,'sheet4','A1:M1000');
    mix_TriGaussian_RBD_stage5=xlsread(framepath_filename_mix_TriGaussian_RBD,'sheet5','A1:M1000');
    
    for aa=1
        for k1=1:1000
            if mix_MultiGaussian_RBD_stage1(k1,1)==0
                k1_end=k1-1;
                break
            end
        end
        for k2=1:1000
            if mix_MultiGaussian_RBD_stage2(k2,1)==0
                k2_end=k2-1;
                break
            end
        end
        if mix_MultiGaussian_RBD_stage2(1000,1)~=0
            k2_end=1000;
        end
        for k3=1:1000
            if mix_MultiGaussian_RBD_stage3(k3,1)==0
                k3_end=k3-1;
                break
            end
        end
        if mix_MultiGaussian_RBD_stage3(1000,1)~=0
            k3_end=1000;
        end
        for k4=1:1000
            if mix_MultiGaussian_RBD_stage4(k4,1)==0
                k4_end=k4-1;
                break
            end
        end
        for k5=1:1000
            if mix_MultiGaussian_RBD_stage5(k5,1)==0
                k5_end=k5-1;
                break
            end
        end
    end
    mix_MultiGaussian_RBD=[mix_MultiGaussian_RBD_stage1(1:k1_end,:);mix_MultiGaussian_RBD_stage2(1:k2_end,:);mix_MultiGaussian_RBD_stage3(1:k3_end,:);mix_MultiGaussian_RBD_stage4(1:k4_end,:);mix_MultiGaussian_RBD_stage5(1:k5_end,:)];
    mix_EquivalentTriGaussian_RBD=[mix_EquivalentTriGaussian_RBD_stage1(1:k1_end,:);mix_EquivalentTriGaussian_RBD_stage2(1:k2_end,:);mix_EquivalentTriGaussian_RBD_stage3(1:k3_end,:);mix_EquivalentTriGaussian_RBD_stage4(1:k4_end,:);mix_EquivalentTriGaussian_RBD_stage5(1:k5_end,:)];
    mix_TriGaussian_RBD=[mix_TriGaussian_RBD_stage1(1:k1_end,:);mix_TriGaussian_RBD_stage2(1:k2_end,:);mix_TriGaussian_RBD_stage3(1:k3_end,:);mix_TriGaussian_RBD_stage4(1:k4_end,:);mix_TriGaussian_RBD_stage5(1:k5_end,:)];
    
    % 读取mix_VD1变量矩阵
    mix_MultiGaussian_VD1_stage1=xlsread(framepath_filename_mix_MultiGaussian_VD1,'sheet1','A1:Y1000');
    mix_MultiGaussian_VD1_stage2=xlsread(framepath_filename_mix_MultiGaussian_VD1,'sheet2','A1:Y1000');
    mix_MultiGaussian_VD1_stage3=xlsread(framepath_filename_mix_MultiGaussian_VD1,'sheet3','A1:Y1000');
    mix_MultiGaussian_VD1_stage4=xlsread(framepath_filename_mix_MultiGaussian_VD1,'sheet4','A1:Y1000');
    mix_MultiGaussian_VD1_stage5=xlsread(framepath_filename_mix_MultiGaussian_VD1,'sheet5','A1:Y1000');
    mix_EquivalentTriGaussian_VD1_stage1=xlsread(framepath_filename_mix_EquivalentTriGaussian_VD1,'sheet1','A1:M1000');
    mix_EquivalentTriGaussian_VD1_stage2=xlsread(framepath_filename_mix_EquivalentTriGaussian_VD1,'sheet2','A1:M1000');
    mix_EquivalentTriGaussian_VD1_stage3=xlsread(framepath_filename_mix_EquivalentTriGaussian_VD1,'sheet3','A1:M1000');
    mix_EquivalentTriGaussian_VD1_stage4=xlsread(framepath_filename_mix_EquivalentTriGaussian_VD1,'sheet4','A1:M1000');
    mix_EquivalentTriGaussian_VD1_stage5=xlsread(framepath_filename_mix_EquivalentTriGaussian_VD1,'sheet5','A1:M1000');
    mix_TriGaussian_VD1_stage1=xlsread(framepath_filename_mix_TriGaussian_VD1,'sheet1','A1:M1000');
    mix_TriGaussian_VD1_stage2=xlsread(framepath_filename_mix_TriGaussian_VD1,'sheet2','A1:M1000');
    mix_TriGaussian_VD1_stage3=xlsread(framepath_filename_mix_TriGaussian_VD1,'sheet3','A1:M1000');
    mix_TriGaussian_VD1_stage4=xlsread(framepath_filename_mix_TriGaussian_VD1,'sheet4','A1:M1000');
    mix_TriGaussian_VD1_stage5=xlsread(framepath_filename_mix_TriGaussian_VD1,'sheet5','A1:M1000');
    
    for aa=1
        for k1=1:1000
            if mix_MultiGaussian_VD1_stage1(k1,1)==0
                k1_end=k1-1;
                break
            end
        end
        for k2=1:1000
            if mix_MultiGaussian_VD1_stage2(k2,1)==0
                k2_end=k2-1;
                break
            end
        end
        if mix_MultiGaussian_VD1_stage2(1000,1)~=0
            k2_end=1000;
        end
        for k3=1:1000
            if mix_MultiGaussian_VD1_stage3(k3,1)==0
                k3_end=k3-1;
                break
            end
        end
        if mix_MultiGaussian_VD1_stage3(1000,1)~=0
            k3_end=1000;
        end
        for k4=1:1000
            if mix_MultiGaussian_VD1_stage4(k4,1)==0
                k4_end=k4-1;
                break
            end
        end
        for k5=1:1000
            if mix_MultiGaussian_VD1_stage5(k5,1)==0
                k5_end=k5-1;
                break
            end
        end
    end
    mix_MultiGaussian_VD1=[mix_MultiGaussian_VD1_stage1(1:k1_end,:);mix_MultiGaussian_VD1_stage2(1:k2_end,:);mix_MultiGaussian_VD1_stage3(1:k3_end,:);mix_MultiGaussian_VD1_stage4(1:k4_end,:);mix_MultiGaussian_VD1_stage5(1:k5_end,:)];
    mix_EquivalentTriGaussian_VD1=[mix_EquivalentTriGaussian_VD1_stage1(1:k1_end,:);mix_EquivalentTriGaussian_VD1_stage2(1:k2_end,:);mix_EquivalentTriGaussian_VD1_stage3(1:k3_end,:);mix_EquivalentTriGaussian_VD1_stage4(1:k4_end,:);mix_EquivalentTriGaussian_VD1_stage5(1:k5_end,:)];
    mix_TriGaussian_VD1=[mix_TriGaussian_VD1_stage1(1:k1_end,:);mix_TriGaussian_VD1_stage2(1:k2_end,:);mix_TriGaussian_VD1_stage3(1:k3_end,:);mix_TriGaussian_VD1_stage4(1:k4_end,:);mix_TriGaussian_VD1_stage5(1:k5_end,:)];
    
    % 读取mix_talinA变量矩阵
    mix_MultiGaussian_talinA_stage1=xlsread(framepath_filename_mix_MultiGaussian_talinA,'sheet1','A1:Y1000');
    mix_MultiGaussian_talinA_stage2=xlsread(framepath_filename_mix_MultiGaussian_talinA,'sheet2','A1:Y1000');
    mix_MultiGaussian_talinA_stage3=xlsread(framepath_filename_mix_MultiGaussian_talinA,'sheet3','A1:Y1000');
    mix_MultiGaussian_talinA_stage4=xlsread(framepath_filename_mix_MultiGaussian_talinA,'sheet4','A1:Y1000');
    mix_MultiGaussian_talinA_stage5=xlsread(framepath_filename_mix_MultiGaussian_talinA,'sheet5','A1:Y1000');
    mix_EquivalentTriGaussian_talinA_stage1=xlsread(framepath_filename_mix_EquivalentTriGaussian_talinA,'sheet1','A1:M1000');
    mix_EquivalentTriGaussian_talinA_stage2=xlsread(framepath_filename_mix_EquivalentTriGaussian_talinA,'sheet2','A1:M1000');
    mix_EquivalentTriGaussian_talinA_stage3=xlsread(framepath_filename_mix_EquivalentTriGaussian_talinA,'sheet3','A1:M1000');
    mix_EquivalentTriGaussian_talinA_stage4=xlsread(framepath_filename_mix_EquivalentTriGaussian_talinA,'sheet4','A1:M1000');
    mix_EquivalentTriGaussian_talinA_stage5=xlsread(framepath_filename_mix_EquivalentTriGaussian_talinA,'sheet5','A1:M1000');
    mix_TriGaussian_talinA_stage1=xlsread(framepath_filename_mix_TriGaussian_talinA,'sheet1','A1:M1000');
    mix_TriGaussian_talinA_stage2=xlsread(framepath_filename_mix_TriGaussian_talinA,'sheet2','A1:M1000');
    mix_TriGaussian_talinA_stage3=xlsread(framepath_filename_mix_TriGaussian_talinA,'sheet3','A1:M1000');
    mix_TriGaussian_talinA_stage4=xlsread(framepath_filename_mix_TriGaussian_talinA,'sheet4','A1:M1000');
    mix_TriGaussian_talinA_stage5=xlsread(framepath_filename_mix_TriGaussian_talinA,'sheet5','A1:M1000');
    
    for aa=1
        for k1=1:1000
            if mix_MultiGaussian_talinA_stage1(k1,1)==0
                k1_end=k1-1;
                break
            end
        end
        for k2=1:1000
            if mix_MultiGaussian_talinA_stage2(k2,1)==0
                k2_end=k2-1;
                break
            end
        end
        if mix_MultiGaussian_talinA_stage2(1000,1)~=0
            k2_end=1000;
        end
        for k3=1:1000
            if mix_MultiGaussian_talinA_stage3(k3,1)==0
                k3_end=k3-1;
                break
            end
        end
        if mix_MultiGaussian_talinA_stage3(1000,1)~=0
            k3_end=1000;
        end
        for k4=1:1000
            if mix_MultiGaussian_talinA_stage4(k4,1)==0
                k4_end=k4-1;
                break
            end
        end
        for k5=1:1000
            if mix_MultiGaussian_talinA_stage5(k5,1)==0
                k5_end=k5-1;
                break
            end
        end
    end
    mix_MultiGaussian_talinA=[mix_MultiGaussian_talinA_stage1(1:k1_end,:);mix_MultiGaussian_talinA_stage2(1:k2_end,:);mix_MultiGaussian_talinA_stage3(1:k3_end,:);mix_MultiGaussian_talinA_stage4(1:k4_end,:);mix_MultiGaussian_talinA_stage5(1:k5_end,:)];
    mix_EquivalentTriGaussian_talinA=[mix_EquivalentTriGaussian_talinA_stage1(1:k1_end,:);mix_EquivalentTriGaussian_talinA_stage2(1:k2_end,:);mix_EquivalentTriGaussian_talinA_stage3(1:k3_end,:);mix_EquivalentTriGaussian_talinA_stage4(1:k4_end,:);mix_EquivalentTriGaussian_talinA_stage5(1:k5_end,:)];
    mix_TriGaussian_talinA=[mix_TriGaussian_talinA_stage1(1:k1_end,:);mix_TriGaussian_talinA_stage2(1:k2_end,:);mix_TriGaussian_talinA_stage3(1:k3_end,:);mix_TriGaussian_talinA_stage4(1:k4_end,:);mix_TriGaussian_talinA_stage5(1:k5_end,:)];
    
    % 读取mix_Cdc42WT变量矩阵
    mix_MultiGaussian_Cdc42WT_stage1=xlsread(framepath_filename_mix_MultiGaussian_Cdc42WT,'sheet1','A1:Y1000');
    mix_MultiGaussian_Cdc42WT_stage2=xlsread(framepath_filename_mix_MultiGaussian_Cdc42WT,'sheet2','A1:Y1000');
    mix_MultiGaussian_Cdc42WT_stage3=xlsread(framepath_filename_mix_MultiGaussian_Cdc42WT,'sheet3','A1:Y1000');
    mix_MultiGaussian_Cdc42WT_stage4=xlsread(framepath_filename_mix_MultiGaussian_Cdc42WT,'sheet4','A1:Y1000');
    mix_MultiGaussian_Cdc42WT_stage5=xlsread(framepath_filename_mix_MultiGaussian_Cdc42WT,'sheet5','A1:Y1000');
    mix_EquivalentTriGaussian_Cdc42WT_stage1=xlsread(framepath_filename_mix_EquivalentTriGaussian_Cdc42WT,'sheet1','A1:M1000');
    mix_EquivalentTriGaussian_Cdc42WT_stage2=xlsread(framepath_filename_mix_EquivalentTriGaussian_Cdc42WT,'sheet2','A1:M1000');
    mix_EquivalentTriGaussian_Cdc42WT_stage3=xlsread(framepath_filename_mix_EquivalentTriGaussian_Cdc42WT,'sheet3','A1:M1000');
    mix_EquivalentTriGaussian_Cdc42WT_stage4=xlsread(framepath_filename_mix_EquivalentTriGaussian_Cdc42WT,'sheet4','A1:M1000');
    mix_EquivalentTriGaussian_Cdc42WT_stage5=xlsread(framepath_filename_mix_EquivalentTriGaussian_Cdc42WT,'sheet5','A1:M1000');
    mix_TriGaussian_Cdc42WT_stage1=xlsread(framepath_filename_mix_TriGaussian_Cdc42WT,'sheet1','A1:M1000');
    mix_TriGaussian_Cdc42WT_stage2=xlsread(framepath_filename_mix_TriGaussian_Cdc42WT,'sheet2','A1:M1000');
    mix_TriGaussian_Cdc42WT_stage3=xlsread(framepath_filename_mix_TriGaussian_Cdc42WT,'sheet3','A1:M1000');
    mix_TriGaussian_Cdc42WT_stage4=xlsread(framepath_filename_mix_TriGaussian_Cdc42WT,'sheet4','A1:M1000');
    mix_TriGaussian_Cdc42WT_stage5=xlsread(framepath_filename_mix_TriGaussian_Cdc42WT,'sheet5','A1:M1000');
    
    for aa=1
        for k1=1:1000
            if mix_MultiGaussian_Cdc42WT_stage1(k1,1)==0
                k1_end=k1-1;
                break
            end
        end
        for k2=1:1000
            if mix_MultiGaussian_Cdc42WT_stage2(k2,1)==0
                k2_end=k2-1;
                break
            end
        end
        if mix_MultiGaussian_Cdc42WT_stage2(1000,1)~=0
            k2_end=1000;
        end
        for k3=1:1000
            if mix_MultiGaussian_Cdc42WT_stage3(k3,1)==0
                k3_end=k3-1;
                break
            end
        end
        if mix_MultiGaussian_Cdc42WT_stage3(1000,1)~=0
            k3_end=1000;
        end
        for k4=1:1000
            if mix_MultiGaussian_Cdc42WT_stage4(k4,1)==0
                k4_end=k4-1;
                break
            end
        end
        for k5=1:1000
            if mix_MultiGaussian_Cdc42WT_stage5(k5,1)==0
                k5_end=k5-1;
                break
            end
        end
    end
    mix_MultiGaussian_Cdc42WT=[mix_MultiGaussian_Cdc42WT_stage1(1:k1_end,:);mix_MultiGaussian_Cdc42WT_stage2(1:k2_end,:);mix_MultiGaussian_Cdc42WT_stage3(1:k3_end,:);mix_MultiGaussian_Cdc42WT_stage4(1:k4_end,:);mix_MultiGaussian_Cdc42WT_stage5(1:k5_end,:)];
    mix_EquivalentTriGaussian_Cdc42WT=[mix_EquivalentTriGaussian_Cdc42WT_stage1(1:k1_end,:);mix_EquivalentTriGaussian_Cdc42WT_stage2(1:k2_end,:);mix_EquivalentTriGaussian_Cdc42WT_stage3(1:k3_end,:);mix_EquivalentTriGaussian_Cdc42WT_stage4(1:k4_end,:);mix_EquivalentTriGaussian_Cdc42WT_stage5(1:k5_end,:)];
    mix_TriGaussian_Cdc42WT=[mix_TriGaussian_Cdc42WT_stage1(1:k1_end,:);mix_TriGaussian_Cdc42WT_stage2(1:k2_end,:);mix_TriGaussian_Cdc42WT_stage3(1:k3_end,:);mix_TriGaussian_Cdc42WT_stage4(1:k4_end,:);mix_TriGaussian_Cdc42WT_stage5(1:k5_end,:)];
    
    % 读取mix_Cdc42D118A变量矩阵
    mix_MultiGaussian_Cdc42D118A_stage1=xlsread(framepath_filename_mix_MultiGaussian_Cdc42D118A,'sheet1','A1:Y1000');
    mix_MultiGaussian_Cdc42D118A_stage2=xlsread(framepath_filename_mix_MultiGaussian_Cdc42D118A,'sheet2','A1:Y1000');
    mix_MultiGaussian_Cdc42D118A_stage3=xlsread(framepath_filename_mix_MultiGaussian_Cdc42D118A,'sheet3','A1:Y1000');
    mix_MultiGaussian_Cdc42D118A_stage4=xlsread(framepath_filename_mix_MultiGaussian_Cdc42D118A,'sheet4','A1:Y1000');
    mix_MultiGaussian_Cdc42D118A_stage5=xlsread(framepath_filename_mix_MultiGaussian_Cdc42D118A,'sheet5','A1:Y1000');
    mix_EquivalentTriGaussian_Cdc42D118A_stage1=xlsread(framepath_filename_mix_EquivalentTriGaussian_Cdc42D118A,'sheet1','A1:M1000');
    mix_EquivalentTriGaussian_Cdc42D118A_stage2=xlsread(framepath_filename_mix_EquivalentTriGaussian_Cdc42D118A,'sheet2','A1:M1000');
    mix_EquivalentTriGaussian_Cdc42D118A_stage3=xlsread(framepath_filename_mix_EquivalentTriGaussian_Cdc42D118A,'sheet3','A1:M1000');
    mix_EquivalentTriGaussian_Cdc42D118A_stage4=xlsread(framepath_filename_mix_EquivalentTriGaussian_Cdc42D118A,'sheet4','A1:M1000');
    mix_EquivalentTriGaussian_Cdc42D118A_stage5=xlsread(framepath_filename_mix_EquivalentTriGaussian_Cdc42D118A,'sheet5','A1:M1000');
    mix_TriGaussian_Cdc42D118A_stage1=xlsread(framepath_filename_mix_TriGaussian_Cdc42D118A,'sheet1','A1:M1000');
    mix_TriGaussian_Cdc42D118A_stage2=xlsread(framepath_filename_mix_TriGaussian_Cdc42D118A,'sheet2','A1:M1000');
    mix_TriGaussian_Cdc42D118A_stage3=xlsread(framepath_filename_mix_TriGaussian_Cdc42D118A,'sheet3','A1:M1000');
    mix_TriGaussian_Cdc42D118A_stage4=xlsread(framepath_filename_mix_TriGaussian_Cdc42D118A,'sheet4','A1:M1000');
    mix_TriGaussian_Cdc42D118A_stage5=xlsread(framepath_filename_mix_TriGaussian_Cdc42D118A,'sheet5','A1:M1000');
    
    for aa=1
        for k1=1:1000
            if mix_MultiGaussian_Cdc42D118A_stage1(k1,1)==0
                k1_end=k1-1;
                break
            end
        end
        for k2=1:1000
            if mix_MultiGaussian_Cdc42D118A_stage2(k2,1)==0
                k2_end=k2-1;
                break
            end
        end
        if mix_MultiGaussian_Cdc42D118A_stage2(1000,1)~=0
            k2_end=1000;
        end
        for k3=1:1000
            if mix_MultiGaussian_Cdc42D118A_stage3(k3,1)==0
                k3_end=k3-1;
                break
            end
        end
        if mix_MultiGaussian_Cdc42D118A_stage3(1000,1)~=0
            k3_end=1000;
        end
        for k4=1:1000
            if mix_MultiGaussian_Cdc42D118A_stage4(k4,1)==0
                k4_end=k4-1;
                break
            end
        end
        for k5=1:1000
            if mix_MultiGaussian_Cdc42D118A_stage5(k5,1)==0
                k5_end=k5-1;
                break
            end
        end
    end
    mix_MultiGaussian_Cdc42D118A=[mix_MultiGaussian_Cdc42D118A_stage1(1:k1_end,:);mix_MultiGaussian_Cdc42D118A_stage2(1:k2_end,:);mix_MultiGaussian_Cdc42D118A_stage3(1:k3_end,:);mix_MultiGaussian_Cdc42D118A_stage4(1:k4_end,:);mix_MultiGaussian_Cdc42D118A_stage5(1:k5_end,:)];
    mix_EquivalentTriGaussian_Cdc42D118A=[mix_EquivalentTriGaussian_Cdc42D118A_stage1(1:k1_end,:);mix_EquivalentTriGaussian_Cdc42D118A_stage2(1:k2_end,:);mix_EquivalentTriGaussian_Cdc42D118A_stage3(1:k3_end,:);mix_EquivalentTriGaussian_Cdc42D118A_stage4(1:k4_end,:);mix_EquivalentTriGaussian_Cdc42D118A_stage5(1:k5_end,:)];
    mix_TriGaussian_Cdc42D118A=[mix_TriGaussian_Cdc42D118A_stage1(1:k1_end,:);mix_TriGaussian_Cdc42D118A_stage2(1:k2_end,:);mix_TriGaussian_Cdc42D118A_stage3(1:k3_end,:);mix_TriGaussian_Cdc42D118A_stage4(1:k4_end,:);mix_TriGaussian_Cdc42D118A_stage5(1:k5_end,:)];
    
    
    % 将各个时期、各种蛋白的数据赋值给mix_all矩阵
    mix_MultiGaussian_all=zeros(1000,25,5,9); % 储存某个时期、某种蛋白下的wucha_best filename data 1-25 changdubi
    mix_EquivalentTriGaussian_all=zeros(1000,13,5,9); % 储存某个时期、某种蛋白下的wucha_best filename data 1-25 changdubi
    mix_TriGaussian_all=zeros(1000,13,5,9); % 储存某个时期、某种蛋白下的wucha_best filename data 1-25 changdubi
    mix_MultiGaussian_all(:,:,1,1)=mix_MultiGaussian_lifeact_stage1;mix_MultiGaussian_all(:,:,2,1)=mix_MultiGaussian_lifeact_stage2;mix_MultiGaussian_all(:,:,3,1)=mix_MultiGaussian_lifeact_stage3;mix_MultiGaussian_all(:,:,4,1)=mix_MultiGaussian_lifeact_stage4;mix_MultiGaussian_all(:,:,5,1)=mix_MultiGaussian_lifeact_stage5;
    mix_MultiGaussian_all(:,:,1,2)=mix_MultiGaussian_actin_stage1;mix_MultiGaussian_all(:,:,2,2)=mix_MultiGaussian_actin_stage2;mix_MultiGaussian_all(:,:,3,2)=mix_MultiGaussian_actin_stage3;mix_MultiGaussian_all(:,:,4,2)=mix_MultiGaussian_actin_stage4;mix_MultiGaussian_all(:,:,5,2)=mix_MultiGaussian_actin_stage5;
    mix_MultiGaussian_all(:,:,1,3)=mix_MultiGaussian_MLC_stage1;mix_MultiGaussian_all(:,:,2,3)=mix_MultiGaussian_MLC_stage2;mix_MultiGaussian_all(:,:,3,3)=mix_MultiGaussian_MLC_stage3;mix_MultiGaussian_all(:,:,4,3)=mix_MultiGaussian_MLC_stage4;mix_MultiGaussian_all(:,:,5,3)=mix_MultiGaussian_MLC_stage5;
    mix_MultiGaussian_all(:,:,1,4)=mix_MultiGaussian_RhoA_stage1;mix_MultiGaussian_all(:,:,2,4)=mix_MultiGaussian_RhoA_stage2;mix_MultiGaussian_all(:,:,3,4)=mix_MultiGaussian_RhoA_stage3;mix_MultiGaussian_all(:,:,4,4)=mix_MultiGaussian_RhoA_stage4;mix_MultiGaussian_all(:,:,5,4)=mix_MultiGaussian_RhoA_stage5;
    mix_MultiGaussian_all(:,:,1,5)=mix_MultiGaussian_RBD_stage1;mix_MultiGaussian_all(:,:,2,5)=mix_MultiGaussian_RBD_stage2;mix_MultiGaussian_all(:,:,3,5)=mix_MultiGaussian_RBD_stage3;mix_MultiGaussian_all(:,:,4,5)=mix_MultiGaussian_RBD_stage4;mix_MultiGaussian_all(:,:,5,5)=mix_MultiGaussian_RBD_stage5;
    mix_MultiGaussian_all(:,:,1,6)=mix_MultiGaussian_VD1_stage1;mix_MultiGaussian_all(:,:,2,6)=mix_MultiGaussian_VD1_stage2;mix_MultiGaussian_all(:,:,3,6)=mix_MultiGaussian_VD1_stage3;mix_MultiGaussian_all(:,:,4,6)=mix_MultiGaussian_VD1_stage4;mix_MultiGaussian_all(:,:,5,6)=mix_MultiGaussian_VD1_stage5;
    mix_MultiGaussian_all(:,:,1,7)=mix_MultiGaussian_talinA_stage1;mix_MultiGaussian_all(:,:,2,7)=mix_MultiGaussian_talinA_stage2;mix_MultiGaussian_all(:,:,3,7)=mix_MultiGaussian_talinA_stage3;mix_MultiGaussian_all(:,:,4,7)=mix_MultiGaussian_talinA_stage4;mix_MultiGaussian_all(:,:,5,7)=mix_MultiGaussian_talinA_stage5;
    mix_MultiGaussian_all(:,:,1,8)=mix_MultiGaussian_Cdc42WT_stage1;mix_MultiGaussian_all(:,:,2,8)=mix_MultiGaussian_Cdc42WT_stage2;mix_MultiGaussian_all(:,:,3,8)=mix_MultiGaussian_Cdc42WT_stage3;mix_MultiGaussian_all(:,:,4,8)=mix_MultiGaussian_Cdc42WT_stage4;mix_MultiGaussian_all(:,:,5,8)=mix_MultiGaussian_Cdc42WT_stage5;
    mix_MultiGaussian_all(:,:,1,9)=mix_MultiGaussian_Cdc42D118A_stage1;mix_MultiGaussian_all(:,:,2,9)=mix_MultiGaussian_Cdc42D118A_stage2;mix_MultiGaussian_all(:,:,3,9)=mix_MultiGaussian_Cdc42D118A_stage3;mix_MultiGaussian_all(:,:,4,9)=mix_MultiGaussian_Cdc42D118A_stage4;mix_MultiGaussian_all(:,:,5,9)=mix_MultiGaussian_Cdc42D118A_stage5;
    
    mix_EquivalentTriGaussian_all(:,:,1,1)=mix_EquivalentTriGaussian_lifeact_stage1;mix_EquivalentTriGaussian_all(:,:,2,1)=mix_EquivalentTriGaussian_lifeact_stage2;mix_EquivalentTriGaussian_all(:,:,3,1)=mix_EquivalentTriGaussian_lifeact_stage3;mix_EquivalentTriGaussian_all(:,:,4,1)=mix_EquivalentTriGaussian_lifeact_stage4;mix_EquivalentTriGaussian_all(:,:,5,1)=mix_EquivalentTriGaussian_lifeact_stage5;
    mix_EquivalentTriGaussian_all(:,:,1,2)=mix_EquivalentTriGaussian_actin_stage1;mix_EquivalentTriGaussian_all(:,:,2,2)=mix_EquivalentTriGaussian_actin_stage2;mix_EquivalentTriGaussian_all(:,:,3,2)=mix_EquivalentTriGaussian_actin_stage3;mix_EquivalentTriGaussian_all(:,:,4,2)=mix_EquivalentTriGaussian_actin_stage4;mix_EquivalentTriGaussian_all(:,:,5,2)=mix_EquivalentTriGaussian_actin_stage5;
    mix_EquivalentTriGaussian_all(:,:,1,3)=mix_EquivalentTriGaussian_MLC_stage1;mix_EquivalentTriGaussian_all(:,:,2,3)=mix_EquivalentTriGaussian_MLC_stage2;mix_EquivalentTriGaussian_all(:,:,3,3)=mix_EquivalentTriGaussian_MLC_stage3;mix_EquivalentTriGaussian_all(:,:,4,3)=mix_EquivalentTriGaussian_MLC_stage4;mix_EquivalentTriGaussian_all(:,:,5,3)=mix_EquivalentTriGaussian_MLC_stage5;
    mix_EquivalentTriGaussian_all(:,:,1,4)=mix_EquivalentTriGaussian_RhoA_stage1;mix_EquivalentTriGaussian_all(:,:,2,4)=mix_EquivalentTriGaussian_RhoA_stage2;mix_EquivalentTriGaussian_all(:,:,3,4)=mix_EquivalentTriGaussian_RhoA_stage3;mix_EquivalentTriGaussian_all(:,:,4,4)=mix_EquivalentTriGaussian_RhoA_stage4;mix_EquivalentTriGaussian_all(:,:,5,4)=mix_EquivalentTriGaussian_RhoA_stage5;
    mix_EquivalentTriGaussian_all(:,:,1,5)=mix_EquivalentTriGaussian_RBD_stage1;mix_EquivalentTriGaussian_all(:,:,2,5)=mix_EquivalentTriGaussian_RBD_stage2;mix_EquivalentTriGaussian_all(:,:,3,5)=mix_EquivalentTriGaussian_RBD_stage3;mix_EquivalentTriGaussian_all(:,:,4,5)=mix_EquivalentTriGaussian_RBD_stage4;mix_EquivalentTriGaussian_all(:,:,5,5)=mix_EquivalentTriGaussian_RBD_stage5;
    mix_EquivalentTriGaussian_all(:,:,1,6)=mix_EquivalentTriGaussian_VD1_stage1;mix_EquivalentTriGaussian_all(:,:,2,6)=mix_EquivalentTriGaussian_VD1_stage2;mix_EquivalentTriGaussian_all(:,:,3,6)=mix_EquivalentTriGaussian_VD1_stage3;mix_EquivalentTriGaussian_all(:,:,4,6)=mix_EquivalentTriGaussian_VD1_stage4;mix_EquivalentTriGaussian_all(:,:,5,6)=mix_EquivalentTriGaussian_VD1_stage5;
    mix_EquivalentTriGaussian_all(:,:,1,7)=mix_EquivalentTriGaussian_talinA_stage1;mix_EquivalentTriGaussian_all(:,:,2,7)=mix_EquivalentTriGaussian_talinA_stage2;mix_EquivalentTriGaussian_all(:,:,3,7)=mix_EquivalentTriGaussian_talinA_stage3;mix_EquivalentTriGaussian_all(:,:,4,7)=mix_EquivalentTriGaussian_talinA_stage4;mix_EquivalentTriGaussian_all(:,:,5,7)=mix_EquivalentTriGaussian_talinA_stage5;
    mix_EquivalentTriGaussian_all(:,:,1,8)=mix_EquivalentTriGaussian_Cdc42WT_stage1;mix_EquivalentTriGaussian_all(:,:,2,8)=mix_EquivalentTriGaussian_Cdc42WT_stage2;mix_EquivalentTriGaussian_all(:,:,3,8)=mix_EquivalentTriGaussian_Cdc42WT_stage3;mix_EquivalentTriGaussian_all(:,:,4,8)=mix_EquivalentTriGaussian_Cdc42WT_stage4;mix_EquivalentTriGaussian_all(:,:,5,8)=mix_EquivalentTriGaussian_Cdc42WT_stage5;
    mix_EquivalentTriGaussian_all(:,:,1,9)=mix_EquivalentTriGaussian_Cdc42D118A_stage1;mix_EquivalentTriGaussian_all(:,:,2,9)=mix_EquivalentTriGaussian_Cdc42D118A_stage2;mix_EquivalentTriGaussian_all(:,:,3,9)=mix_EquivalentTriGaussian_Cdc42D118A_stage3;mix_EquivalentTriGaussian_all(:,:,4,9)=mix_EquivalentTriGaussian_Cdc42D118A_stage4;mix_EquivalentTriGaussian_all(:,:,5,9)=mix_EquivalentTriGaussian_Cdc42D118A_stage5;
    
    mix_TriGaussian_all(:,:,1,1)=mix_TriGaussian_lifeact_stage1;mix_TriGaussian_all(:,:,2,1)=mix_TriGaussian_lifeact_stage2;mix_TriGaussian_all(:,:,3,1)=mix_TriGaussian_lifeact_stage3;mix_TriGaussian_all(:,:,4,1)=mix_TriGaussian_lifeact_stage4;mix_TriGaussian_all(:,:,5,1)=mix_TriGaussian_lifeact_stage5;
    mix_TriGaussian_all(:,:,1,2)=mix_TriGaussian_actin_stage1;mix_TriGaussian_all(:,:,2,2)=mix_TriGaussian_actin_stage2;mix_TriGaussian_all(:,:,3,2)=mix_TriGaussian_actin_stage3;mix_TriGaussian_all(:,:,4,2)=mix_TriGaussian_actin_stage4;mix_TriGaussian_all(:,:,5,2)=mix_TriGaussian_actin_stage5;
    mix_TriGaussian_all(:,:,1,3)=mix_TriGaussian_MLC_stage1;mix_TriGaussian_all(:,:,2,3)=mix_TriGaussian_MLC_stage2;mix_TriGaussian_all(:,:,3,3)=mix_TriGaussian_MLC_stage3;mix_TriGaussian_all(:,:,4,3)=mix_TriGaussian_MLC_stage4;mix_TriGaussian_all(:,:,5,3)=mix_TriGaussian_MLC_stage5;
    mix_TriGaussian_all(:,:,1,4)=mix_TriGaussian_RhoA_stage1;mix_TriGaussian_all(:,:,2,4)=mix_TriGaussian_RhoA_stage2;mix_TriGaussian_all(:,:,3,4)=mix_TriGaussian_RhoA_stage3;mix_TriGaussian_all(:,:,4,4)=mix_TriGaussian_RhoA_stage4;mix_TriGaussian_all(:,:,5,4)=mix_TriGaussian_RhoA_stage5;
    mix_TriGaussian_all(:,:,1,5)=mix_TriGaussian_RBD_stage1;mix_TriGaussian_all(:,:,2,5)=mix_TriGaussian_RBD_stage2;mix_TriGaussian_all(:,:,3,5)=mix_TriGaussian_RBD_stage3;mix_TriGaussian_all(:,:,4,5)=mix_TriGaussian_RBD_stage4;mix_TriGaussian_all(:,:,5,5)=mix_TriGaussian_RBD_stage5;
    mix_TriGaussian_all(:,:,1,6)=mix_TriGaussian_VD1_stage1;mix_TriGaussian_all(:,:,2,6)=mix_TriGaussian_VD1_stage2;mix_TriGaussian_all(:,:,3,6)=mix_TriGaussian_VD1_stage3;mix_TriGaussian_all(:,:,4,6)=mix_TriGaussian_VD1_stage4;mix_TriGaussian_all(:,:,5,6)=mix_TriGaussian_VD1_stage5;
    mix_TriGaussian_all(:,:,1,7)=mix_TriGaussian_talinA_stage1;mix_TriGaussian_all(:,:,2,7)=mix_TriGaussian_talinA_stage2;mix_TriGaussian_all(:,:,3,7)=mix_TriGaussian_talinA_stage3;mix_TriGaussian_all(:,:,4,7)=mix_TriGaussian_talinA_stage4;mix_TriGaussian_all(:,:,5,7)=mix_TriGaussian_talinA_stage5;
    mix_TriGaussian_all(:,:,1,8)=mix_TriGaussian_Cdc42WT_stage1;mix_TriGaussian_all(:,:,2,8)=mix_TriGaussian_Cdc42WT_stage2;mix_TriGaussian_all(:,:,3,8)=mix_TriGaussian_Cdc42WT_stage3;mix_TriGaussian_all(:,:,4,8)=mix_TriGaussian_Cdc42WT_stage4;mix_TriGaussian_all(:,:,5,8)=mix_TriGaussian_Cdc42WT_stage5;
    mix_TriGaussian_all(:,:,1,9)=mix_TriGaussian_Cdc42D118A_stage1;mix_TriGaussian_all(:,:,2,9)=mix_TriGaussian_Cdc42D118A_stage2;mix_TriGaussian_all(:,:,3,9)=mix_TriGaussian_Cdc42D118A_stage3;mix_TriGaussian_all(:,:,4,9)=mix_TriGaussian_Cdc42D118A_stage4;mix_TriGaussian_all(:,:,5,9)=mix_TriGaussian_Cdc42D118A_stage5;
end


%%% 读取各蛋白数据文件中保存的文件名，保存到filenamearray_all矩阵中，用于后续查找该细胞的原始灰度值文件
filenamearray_all=cell(1000,5,9);  % lumen形成前的文件名矩阵

% 文件路径拼接
framepath_filename_lifeact=[framepath filename_lifeact filename_xlsx];
framepath_filename_actin=[framepath filename_actin filename_xlsx];
framepath_filename_MLC=[framepath filename_MLC filename_xlsx];
framepath_filename_RhoA=[framepath filename_RhoA filename_xlsx];
framepath_filename_RBD=[framepath filename_RBD filename_xlsx];
framepath_filename_VD1=[framepath filename_VD1 filename_xlsx];
framepath_filename_talinA=[framepath filename_talinA filename_xlsx];
framepath_filename_Cdc42WT=[framepath filename_Cdc42WT filename_xlsx];
framepath_filename_Cdc42D118A=[framepath filename_Cdc42D118A filename_xlsx];

% 文件名读取，保存到filenamearray_all矩阵
for aa=1
    % 读取lifeact对应的文件名
    [~, ~, filenamearray_all(:,1,1)] = xlsread(framepath_filename_lifeact,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,1)] = xlsread(framepath_filename_lifeact,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,1)] = xlsread(framepath_filename_lifeact,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,1)] = xlsread(framepath_filename_lifeact,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,1)] = xlsread(framepath_filename_lifeact,'sheet5','A2:A1001');
    
    % 读取actin对应的文件名
    [~, ~, filenamearray_all(:,1,2)] = xlsread(framepath_filename_actin,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,2)] = xlsread(framepath_filename_actin,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,2)] = xlsread(framepath_filename_actin,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,2)] = xlsread(framepath_filename_actin,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,2)] = xlsread(framepath_filename_actin,'sheet5','A2:A1001');
    
    % 读取MLC对应的文件名
    [~, ~, filenamearray_all(:,1,3)] = xlsread(framepath_filename_MLC,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,3)] = xlsread(framepath_filename_MLC,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,3)] = xlsread(framepath_filename_MLC,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,3)] = xlsread(framepath_filename_MLC,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,3)] = xlsread(framepath_filename_MLC,'sheet5','A2:A1001');
    
    % 读取RhoA对应的文件名
    [~, ~, filenamearray_all(:,1,4)] = xlsread(framepath_filename_RhoA,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,4)] = xlsread(framepath_filename_RhoA,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,4)] = xlsread(framepath_filename_RhoA,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,4)] = xlsread(framepath_filename_RhoA,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,4)] = xlsread(framepath_filename_RhoA,'sheet5','A2:A1001');
    
    % 读取RBD对应的文件名
    [~, ~, filenamearray_all(:,1,5)] = xlsread(framepath_filename_RBD,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,5)] = xlsread(framepath_filename_RBD,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,5)] = xlsread(framepath_filename_RBD,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,5)] = xlsread(framepath_filename_RBD,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,5)] = xlsread(framepath_filename_RBD,'sheet5','A2:A1001');
    
    % 读取VD1对应的文件名
    [~, ~, filenamearray_all(:,1,6)] = xlsread(framepath_filename_VD1,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,6)] = xlsread(framepath_filename_VD1,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,6)] = xlsread(framepath_filename_VD1,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,6)] = xlsread(framepath_filename_VD1,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,6)] = xlsread(framepath_filename_VD1,'sheet5','A2:A1001');
    
    % 读取talinA对应的文件名
    [~, ~, filenamearray_all(:,1,7)] = xlsread(framepath_filename_talinA,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,7)] = xlsread(framepath_filename_talinA,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,7)] = xlsread(framepath_filename_talinA,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,7)] = xlsread(framepath_filename_talinA,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,7)] = xlsread(framepath_filename_talinA,'sheet5','A2:A1001');
    
    % 读取Cdc42WT对应的文件名
    [~, ~, filenamearray_all(:,1,8)] = xlsread(framepath_filename_Cdc42WT,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,8)] = xlsread(framepath_filename_Cdc42WT,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,8)] = xlsread(framepath_filename_Cdc42WT,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,8)] = xlsread(framepath_filename_Cdc42WT,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,8)] = xlsread(framepath_filename_Cdc42WT,'sheet5','A2:A1001');
    
    % 读取Cdc42D118A对应的文件名
    [~, ~, filenamearray_all(:,1,9)] = xlsread(framepath_filename_Cdc42D118A,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,9)] = xlsread(framepath_filename_Cdc42D118A,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,9)] = xlsread(framepath_filename_Cdc42D118A,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,9)] = xlsread(framepath_filename_Cdc42D118A,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,9)] = xlsread(framepath_filename_Cdc42D118A,'sheet5','A2:A1001');
end


%% step 2

% mix_sort矩阵赋值
mix_MultiGaussian_sort=mix_MultiGaussian_all;  % 用于排序
mix_EquivalentTriGaussian_sort=mix_EquivalentTriGaussian_all;
mix_TriGaussian_sort=mix_TriGaussian_all;

% 使用冒泡排序法进行排序
% 首先根据误差大小排序，并将误差值大于180的细胞数据剔除
cell_number=zeros(5,9);
for i=1:5
    for j=1:9
        for k=1:1000  % 找出该矩阵中的实际数据组数
            if mix_TriGaussian_sort(k,1,i,j)==0
                break
            end
        end
        cell_number(i,j)=k-1; % 将实际数组数(细胞个数)赋值给cell_number
        
        % 冒泡排序
        for ii=2:cell_number(i,j)
            for jj=cell_number(i,j):-1:ii
                if mix_TriGaussian_sort(jj,1,i,j)<mix_TriGaussian_sort(jj-1,1,i,j)  % 若前一个细胞的误差大于后一个细胞的误差，则调换顺序
                    temp=mix_TriGaussian_sort(jj,:,i,j);
                    mix_TriGaussian_sort(jj,:,i,j)=mix_TriGaussian_sort(jj-1,:,i,j);
                    mix_TriGaussian_sort(jj-1,:,i,j)=temp;
                    temp2=mix_MultiGaussian_sort(jj,:,i,j);
                    mix_MultiGaussian_sort(jj,:,i,j)=mix_MultiGaussian_sort(jj-1,:,i,j);
                    mix_MultiGaussian_sort(jj-1,:,i,j)=temp2;
                    temp3=mix_EquivalentTriGaussian_sort(jj,:,i,j);
                    mix_EquivalentTriGaussian_sort(jj,:,i,j)=mix_EquivalentTriGaussian_sort(jj-1,:,i,j);
                    mix_EquivalentTriGaussian_sort(jj-1,:,i,j)=temp3;
                    temp4=filenamearray_all(jj,i,j);
                    filenamearray_all(jj,i,j)=filenamearray_all(jj-1,i,j);
                    filenamearray_all(jj-1,i,j)=temp4;
                end
            end
        end
        
        % 删除误差值大于180的细胞数据
        for k=1:cell_number(i,j)
            if mix_TriGaussian_sort(k,1,i,j)>=180
                mix_TriGaussian_sort(k,:,i,j)=0;  % 剔除误差大于8的细胞数据
                mix_MultiGaussian_sort(k,:,i,j)=0;
                mix_EquivalentTriGaussian_sort(k,:,i,j)=0;
                filenamearray_all{k,i,j}='0';
                cell_number(i,j)=cell_number(i,j)-1;
            end
        end
    end
end

% 删除赤道轴收缩环特征参数和两侧收缩环特征参数过于偏离平均的细胞，剔除三者前10%和后10%的细胞数据的并集
memory_cell_wucha=zeros(200,5,9,6);  % 储存待删除细胞数据对应的误差值
k_cellname=ones(5,9,6);  % memory_cellname的数组长度

% 计算basal domain多高斯的平均强度和聚集程度，储存到mix_sort_basalave矩阵中，用于后续排序查找需删除数据
mix_sort_basalave=zeros(1000,13,5,9);
for i=1:5
    for j=1:9
        for k=1:cell_number(i,j)
            for ii=4:length(mix_MultiGaussian_sort(1,:,1,1))
                if mix_MultiGaussian_sort(k,ii,i,j)==0  % 若第ii位为0，说明basal domain有(ii-10)/3个高斯
                    break
                end
            end
            if mix_MultiGaussian_sort(k,end-1,i,j)~=0  % 若最后一位不为0，说明basal domain有5(max)个高斯
                ii=25;
            end
            basal_Gaussian_number=(ii-10)/3;

            % 计算basal domain多高斯的平均强度和聚集程度
            basal_overactivity_average=sum(mix_MultiGaussian_all(k,4:basal_Gaussian_number+3,i,j))/basal_Gaussian_number;
            basal_signalwidth_average=sum(mix_MultiGaussian_all(k,basal_Gaussian_number+4:2*basal_Gaussian_number+3,i,j))/basal_Gaussian_number;
            basal_expectation_average=sum(mix_MultiGaussian_all(k,2*basal_Gaussian_number+4:3*basal_Gaussian_number+3,i,j))/basal_Gaussian_number;

            % 赋值给mix_sort_basalave
            for ii=1:length(mix_sort_basalave(1,:,1,1))
                if ii==4
                    mix_sort_basalave(k,ii,i,j)=basal_overactivity_average;
                elseif ii==5
                    mix_sort_basalave(k,ii,i,j)=basal_signalwidth_average;
                elseif ii==6
                    mix_sort_basalave(k,ii,i,j)=basal_expectation_average;
                elseif ii==13
                    mix_sort_basalave(k,ii,i,j)=mix_TriGaussian_sort(k,end,i,j);
                elseif ii<4
                    mix_sort_basalave(k,ii,i,j)=mix_TriGaussian_sort(k,ii,i,j);
                elseif ii>=7 && ii<=9
                    mix_sort_basalave(k,ii,i,j)=mix_MultiGaussian_sort(k,(ii-6)*(basal_Gaussian_number+2)+2,i,j);
                else
                    mix_sort_basalave(k,ii,i,j)=mix_MultiGaussian_sort(k,(ii-9)*(basal_Gaussian_number+2)+3,i,j);
                end
            end
        end
    end
end

% 两侧收缩环特征参数使用多高斯拟合的结果为依据
number_sort_position=[7 8 10 11];  % 特征参数对应在mix_all矩阵的位置
for k=1:4
    number_sort=number_sort_position(k);
    for i=1:5
        for j=1:9
            % 冒泡排序
            for ii=2:cell_number(i,j)
                for jj=cell_number(i,j):-1:ii
                    if mix_sort_basalave(jj,number_sort,i,j)<mix_sort_basalave(jj-1,number_sort,i,j)  % 若前一个细胞的特征参数大于后一个细胞的特征参数，则调换顺序
                        temp5=mix_sort_basalave(jj,:,i,j);
                        mix_sort_basalave(jj,:,i,j)=mix_sort_basalave(jj-1,:,i,j);
                        mix_sort_basalave(jj-1,:,i,j)=temp5;
                        temp=mix_TriGaussian_sort(jj,:,i,j);
                        mix_TriGaussian_sort(jj,:,i,j)=mix_TriGaussian_sort(jj-1,:,i,j);
                        mix_TriGaussian_sort(jj-1,:,i,j)=temp;
                        temp2=mix_MultiGaussian_sort(jj,:,i,j);
                        mix_MultiGaussian_sort(jj,:,i,j)=mix_MultiGaussian_sort(jj-1,:,i,j);
                        mix_MultiGaussian_sort(jj-1,:,i,j)=temp2;
                        temp3=mix_EquivalentTriGaussian_sort(jj,:,i,j);
                        mix_EquivalentTriGaussian_sort(jj,:,i,j)=mix_EquivalentTriGaussian_sort(jj-1,:,i,j);
                        mix_EquivalentTriGaussian_sort(jj-1,:,i,j)=temp3;
                        temp4=filenamearray_all(jj,i,j);
                        filenamearray_all(jj,i,j)=filenamearray_all(jj-1,i,j);
                        filenamearray_all(jj-1,i,j)=temp4;
                    end
                end
            end

            % 找出前5%的细胞，用误差值标记
            if number_sort==7 || number_sort==10
                for jj=1:floor(cell_number(i,j)/20)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_TriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
                % 找出后5%的细胞，用误差值标记
                for jj=cell_number(i,j)-floor(cell_number(i,j)/20)+1:cell_number(i,j)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_TriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
            else
                for jj=1:floor(cell_number(i,j)/10)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_TriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
                % 找出后5%的细胞，用误差值标记
                for jj=cell_number(i,j)-floor(cell_number(i,j)/10)+1:cell_number(i,j)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_TriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
            end
        end
    end
end

% 赤道轴收缩环特征参数使用直接三高斯和等效三高斯拟合的比较数据(若两种方法的拟合结果差异超过20%，则使用直接三高斯拟合参数，反之使用等效三高斯拟合参数)
number_sort_position=[4 5];  % 特征参数对应在mix_all矩阵的位置
for k=5:6
    number_sort=number_sort_position(k-4);
    for i=1:5
        for j=1:9
            % 冒泡排序
            for ii=2:cell_number(i,j)
                for jj=cell_number(i,j):-1:ii
                    if mix_TriGaussian_sort(jj,number_sort,i,j)<mix_TriGaussian_sort(jj-1,number_sort,i,j)  % 若前一个细胞的特征参数大于后一个细胞的特征参数，则调换顺序
                        temp=mix_TriGaussian_sort(jj,:,i,j);
                        mix_TriGaussian_sort(jj,:,i,j)=mix_TriGaussian_sort(jj-1,:,i,j);
                        mix_TriGaussian_sort(jj-1,:,i,j)=temp;
                        temp2=mix_MultiGaussian_sort(jj,:,i,j);
                        mix_MultiGaussian_sort(jj,:,i,j)=mix_MultiGaussian_sort(jj-1,:,i,j);
                        mix_MultiGaussian_sort(jj-1,:,i,j)=temp2;
                        temp3=mix_EquivalentTriGaussian_sort(jj,:,i,j);
                        mix_EquivalentTriGaussian_sort(jj,:,i,j)=mix_EquivalentTriGaussian_sort(jj-1,:,i,j);
                        mix_EquivalentTriGaussian_sort(jj-1,:,i,j)=temp3;
                        temp4=filenamearray_all(jj,i,j);
                        filenamearray_all(jj,i,j)=filenamearray_all(jj-1,i,j);
                        filenamearray_all(jj-1,i,j)=temp4;
                    end
                end
            end

            % 找出前5%的细胞，用误差值标记
            for jj=1:floor(cell_number(i,j)/20)
                if abs(mix_TriGaussian_sort(jj,number_sort,i,j)-mix_EquivalentTriGaussian_sort(jj,number_sort,i,j))/mix_TriGaussian_sort(jj,number_sort,i,j)>0.2
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_TriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                else
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_EquivalentTriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
            end
            % 找出后5%的细胞，用误差值标记
            for jj=cell_number(i,j)-floor(cell_number(i,j)/20)+1:cell_number(i,j)
                if abs(mix_TriGaussian_sort(jj,number_sort,i,j)-mix_EquivalentTriGaussian_sort(jj,number_sort,i,j))/mix_TriGaussian_sort(jj,number_sort,i,j)>0.2
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_TriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                else
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_EquivalentTriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
            end
        end
    end
end


% 取三次排序得到的细胞的并集，并在mix_all矩阵中删去
mix_TriGaussian_sort2=mix_TriGaussian_sort;  % 将mix_TriGaussian_sort数组内容赋值给mix_TriGaussian_sort2
mix_EquivalentTriGaussian_sort2=mix_EquivalentTriGaussian_sort;  % 将mix_EquivalentTriGaussian_sort数组内容赋值给mix_EquivalentTriGaussian_sort2
mix_MultiGaussian_sort2=mix_MultiGaussian_sort;  % 将mix_MultiGaussian_sort数组内容赋值给mix_MultiGaussian_sort2

mix_TriGaussian_sort=zeros(1000,13,5,9);  % 清空数组，用于后续存放删除非正常细胞后的细胞拟合数据
mix_EquivalentTriGaussian_sort=zeros(1000,13,5,9);  % 清空数组，用于后续存放删除非正常细胞后的细胞拟合数据
mix_MultiGaussian_sort=zeros(1000,25,5,9);  % 清空数组，用于后续存放删除非正常细胞后的细胞拟合数据
filenamearray_all_sort=cell(1000,5,9);
for i=1:5
    for j=1:9
        kk2=1;  % 用于存储mix_all每个时期，每种蛋白的细胞数
        for k1=1:1000
            if mix_TriGaussian_sort2(k1,1,i,j)~=0  % 存在第k1个细胞
                logic=0;  % 用于判断该细胞是否需要删除
                for k=1:6
                    for k2=1:k_cellname(i,j,k)
                        if mix_TriGaussian_sort2(k1,1,i,j)==memory_cell_wucha(k2,i,j,k)
                            logic=1;  % 若该细胞的wucha值,赤道轴收缩环参数,和前侧收缩环参数为前5%和后5%的细胞，则说明该细胞需要删除
                            break
                        end
                    end
                end
                if logic==0  % 若logic为0，则表示该细胞不需删除，将其赋值给mix_sort filenamearray_all矩阵
                    mix_TriGaussian_sort(kk2,:,i,j)=mix_TriGaussian_sort2(k1,:,i,j);
                    mix_EquivalentTriGaussian_sort(kk2,:,i,j)=mix_EquivalentTriGaussian_sort2(k1,:,i,j);
                    mix_MultiGaussian_sort(kk2,:,i,j)=mix_MultiGaussian_sort2(k1,:,i,j);
                    filenamearray_all_sort(kk2,i,j)=filenamearray_all(k1,i,j);
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
mix_TriGaussian_all=mix_TriGaussian_sort;
mix_EquivalentTriGaussian_all=mix_EquivalentTriGaussian_sort;
mix_MultiGaussian_all=mix_MultiGaussian_sort;
filenamearray_all=filenamearray_all_sort;



%%% 计算各时期各拟合参数的均值
% 三高斯特征参数汇总
TriGaussian_summarize=mix_EquivalentTriGaussian_sort;



%% step 3
% 首先将每种蛋白的细胞特征值根据lumen length/cell length由小到大排序
for i=1:5
    for j=1:9
        % 冒泡排序
        for ii=2:cell_number(i,j)
            for jj=cell_number(i,j):-1:ii
                if TriGaussian_summarize(jj,end,i,j)<TriGaussian_summarize(jj-1,end,i,j)  % 若前一个细胞的特征参数大于后一个细胞的特征参数，则调换顺序
                    temp=TriGaussian_summarize(jj,:,i,j);
                    TriGaussian_summarize(jj,:,i,j)=TriGaussian_summarize(jj-1,:,i,j);
                    TriGaussian_summarize(jj-1,:,i,j)=temp;
                    temp2=filenamearray_all(jj,i,j);
                    filenamearray_all(jj,i,j)=filenamearray_all(jj-1,i,j);
                    filenamearray_all(jj-1,i,j)=temp2;
                end
            end
        end
    end
end


pre_cell_number=zeros(9,1);  % 各蛋白的lumen形成前的细胞数
pre_cell_number_monocell=zeros(9,1);  % 各蛋白的lumen形成前的有单个数据文件的细胞数
pre_cell_number_multicell=zeros(9,1);  % 各蛋白的lumen形成前的有三个数据文件的细胞数
pre_AL_thickness=zeros(200,9,2);  % lumen形成前细胞的anterior-laterior thickness，第二个维度表示不同种蛋白，第三个维度表示有三个数据文件的细胞和有单个数据文件的细胞
pre_AL_length=zeros(200,9,2);  % lumen形成前细胞的\anterior-laterior length
pre_PL_thickness=zeros(200,9,2);  % lumen形成前细胞的\posterior-laterior thickness
pre_PL_length=zeros(200,9,2);  % lumen形成前细胞的\posterior-laterior length
pre_average_anteriorlength=zeros(1,9);pre_average_posteriorlength=zeros(1,9);  % 定义anterior-和posterior-lateral domain的cortex平均长度变量
pre_average_anteriorthickness=zeros(1,9);pre_average_posteriorthickness=zeros(1,9);  % 定义anterior-和posterior-lateral domain的cortex平均厚度变量

% 找出其中lumen length/cell length<0.08的细胞，按文件数分类计算其anterior-和posterior-lateral domain的cortex平均厚度和平均长度
for i=1:9
    % 统计lumen形成前的细胞数
    for j=1:cell_number(1,i)
        if TriGaussian_summarize(j,end,1,i)<0.08
            pre_cell_number(i)=pre_cell_number(i)+1;
        else
            break
        end
    end
    
    if pre_cell_number(i)==0  % 判断是否存在lumen形成前的细胞
        continue
    end
    
    for j=1:pre_cell_number(i)
        filename=mat2str(cell2mat(filenamearray_all(j,1,i)));  % 读取文件名
        filename=filename(2:length(filename)-1);  % 删除两侧的单引号
        S=regexp(filename,'_','split');  % 拆分文件名字符串
        month_day=mat2str(cell2mat(S(1)));  % S矩阵的第一个值为该细胞的日期(month+day)
        month_day=month_day(2:length(month_day)-1);  % 删除两侧的单引号
        framepath_filename_1=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_xiahuaxian filename_1 filename_csv];  % 拼接完整的细胞灰度值数据文件原始路径(多个文件形式)
        framepath_filename_2=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_xiahuaxian filename_2 filename_csv];
        framepath_filename_3=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_xiahuaxian filename_3 filename_csv];
        framepath_filename_danyi=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_csv];  % 拼接完整的细胞灰度值数据文件原始路径(单个文件形式)

        % 根据该细胞的收缩环特征参数描绘tri-Gaussian distribution
        rho0=TriGaussian_summarize(j,3,1,i);rho_inf=TriGaussian_summarize(j,4,1,i);w=TriGaussian_summarize(j,5,1,i);E=TriGaussian_summarize(j,6,1,i);  % 赤道轴收缩环特征参数赋值
        rho_inf_anterior=TriGaussian_summarize(j,7,1,i);w_anterior=TriGaussian_summarize(j,8,1,i);E_anterior=TriGaussian_summarize(j,9,1,i);
        rho_inf_posterior=TriGaussian_summarize(j,10,1,i);w_posterior=TriGaussian_summarize(j,11,1,i);E_posterior=TriGaussian_summarize(j,12,1,i);  % 两侧收缩环特征参数赋值
        f=@(x) rho0+rho_inf.*exp(-1/2.*((x-E)./w).^2)+rho_inf_anterior.*exp(-1/2.*((x-E_anterior)./w_anterior).^2)+rho_inf_posterior.*exp(-1/2.*((x-E_posterior)./w_posterior).^2);
        
        % 计算anterior-和posterior-lateral domain的cortex平均厚度和平均长度
        if exist(framepath_filename_1,'file')  % 统计lumen形成前的细胞中有三个数据文件的细胞数，并计算其anterior-和posterior-lateral domain的cortex平均厚度和平均长度
            pre_cell_number_multicell(i)=pre_cell_number_multicell(i)+1;  % 三个数据文件的细胞数统计
            
            data_1=xlsread(framepath_filename_1); % 导入数据
            data_2=xlsread(framepath_filename_2);
            data_3=xlsread(framepath_filename_3);
            pre_AL_length(pre_cell_number_multicell(i),i,1)=length(data_1)/(length(data_1)+length(data_2)+length(data_3));  % anterior-laterior length统计
            pre_PL_length(pre_cell_number_multicell(i),i,1)=length(data_3)/(length(data_1)+length(data_2)+length(data_3));  % posterior-laterior length统计
            
            a_anterior=integral(f,-1/2,-1/2+pre_AL_length(pre_cell_number_multicell(i),i,1));  % anterior-laterior cortex thickness积分
            a_anterior=a_anterior/pre_AL_length(pre_cell_number_multicell(i),i,1);  % 积分结果除以x轴长度，得到平均cortex thickness
            pre_AL_thickness(pre_cell_number_multicell(i),i,1)=a_anterior;  % 平均cortex thickness值赋值给pre_average_AT矩阵
            a_posterior=integral(f,1/2-pre_PL_length(pre_cell_number_multicell(i),i,1),1/2);  % posterior-laterior cortex thickness积分
            a_posterior=a_posterior/pre_PL_length(pre_cell_number_multicell(i),i,1);  % 积分结果除以x轴长度，得到平均cortex thickness
            pre_PL_thickness(pre_cell_number_multicell(i),i,1)=a_posterior;  % 平均cortex thickness值赋值给pre_average_PT矩阵
        elseif exist(framepath_filename_danyi,'file')  % 统计lumen形成前的细胞中有单个数据文件的细胞数，并计算其anterior-和posterior-lateral domain的cortex平均厚度
            pre_cell_number_monocell(i)=pre_cell_number_monocell(i)+1;  % 三个数据文件的细胞数统计
            
            % 根据文件名确定该细胞在该日期下的位置
            number_changdubi=0;
            for ii=1:30  % 外侧循环用于查找目标文件是否存在(i表示第几张photo，j表示photo中的第几个细胞，k表示荧光通道)
                filename1=num2str(ii);
                for jj=1:10
                    filename2=num2str(jj);
                    filename3='1';
                    filename_temp=[month_day filename_xiahuaxian filename1 filename2 filename3];  % 如 0928_111
                    
                    % 拼接完整的文件路径，用于下一步exist函数的查找
                    framepath_filename_1_temp=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_temp filename_xiahuaxian filename_1 filename_csv];  % 拼接完整的细胞灰度值数据文件原始路径(多个文件形式)
                    framepath_filename_danyi_temp=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_temp filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111.csv
                    
                    if exist(framepath_filename_danyi_temp,'file')   % 一个细胞的数据有单个文件
                        number_changdubi=number_changdubi+1;
                    elseif exist(framepath_filename_1_temp,'file')
                        number_changdubi=number_changdubi+1;
                    elseif strcmp(filename_temp,filename)  % 若该张照片下，不存在该细胞的数据文件，说明该张照片中的细胞已经拟合完成，跳出j(细胞)循环
                        number_changdubi=number_changdubi+1;
                        break
                    end
                end
            end
            
            % 查找该细胞对应的lumen length和cell length (从data_changdubi中寻找)
            framepath2=[framepath num2str(k) filename_fanxiegang month_day filename_fanxiegang];  % 子文件夹路径，如 'C:\Users\Desktop\1\0928'
            framepath_filename_changdubi=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_changdubi filename_xlsx];
            framepath_filename_Radius=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_Radius filename_xlsx];
            framepath_filename_ContactAngle=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_ContactAngle filename_xlsx];

            Lumen_Radius=xlsread(framepath_filename_Radius,'sheet1','A1:A100');  % 读入lumen半径
            Lumen_Radius=sqrt(Lumen_Radius./pi);
            Lumen_ContactAngle=xlsread(framepath_filename_ContactAngle,'sheet1','A1:A100');  % 读入lumen接触角
            data_lumenlength=xlsread(framepath_filename_changdubi,'sheet1','A1:A100');  % 读入该文件下的lumen长度
            data_celllength=xlsread(framepath_filename_changdubi,'sheet1','B1:B100');  % 读入该文件下的cell长度

            data=xlsread(framepath_filename_danyi); % 导入数据
            anterior_lateral_length=(data_celllength(number_changdubi)-Lumen_Radius(number_changdubi)*sin(Lumen_ContactAngle(number_changdubi)*pi./180))/2;  % anterior-lateral domain长度
            posterior_lateral_length=(data_celllength(number_changdubi)-Lumen_Radius(number_changdubi)*sin(Lumen_ContactAngle(number_changdubi)*pi./180))/2;  % posterior-lateral domain长度
            basal_length=data(end,1)-anterior_lateral_length-posterior_lateral_length;  % basal domain长度


            a_anterior=integral(f,-1/2,-1/2+anterior_lateral_length/data(end,1));  % anterior-lateral cortex thickness积分
            a_anterior=a_anterior/(anterior_lateral_length/data(end,1));  % 积分结果除以x轴长度，得到平均cortex thickness
            pre_AL_length(pre_cell_number_monocell(i),i,2)=anterior_lateral_length/data(end,1);  % anterior length赋值给pre_AL_length矩阵
            pre_AL_thickness(pre_cell_number_monocell(i),i)=a_anterior;  % 平均cortex thickness值赋值给pre_AL_thickness矩阵
            a_posterior=integral(f,1/2-posterior_lateral_length/data(end,1),1/2);  % posterior-lateral cortex thickness积分
            a_posterior=a_posterior/(posterior_lateral_length/data(end,1));  % 积分结果除以x轴长度，得到平均cortex thickness
            pre_PL_length(pre_cell_number_monocell(i),i,2)=posterior_lateral_length/data(end,1);  % posterior length赋值给pre_LL_length矩阵
            pre_PL_thickness(pre_cell_number_monocell(i),i)=a_posterior;  % 平均cortex thickness值赋值给pre_PL_thickness矩阵
        end
    end

    % 计算lumen形成前的anterior-和posterior-lateral domain的cortex平均厚度
    for j=1:pre_cell_number_multicell(i)
        pre_average_anteriorlength(i)=pre_average_anteriorlength(i)+pre_AL_length(j,i,1);
        pre_average_posteriorlength(i)=pre_average_posteriorlength(i)+pre_PL_length(j,i,1);
        pre_average_anteriorthickness(i)=pre_average_anteriorthickness(i)+pre_AL_thickness(j,i,1);
        pre_average_posteriorthickness(i)=pre_average_posteriorthickness(i)+pre_PL_thickness(j,i,1);
    end
    for j=1:pre_cell_number_monocell(i)
        pre_average_anteriorlength(i)=pre_average_anteriorlength(i)+pre_AL_length(j,i,2);
        pre_average_posteriorlength(i)=pre_average_posteriorlength(i)+pre_PL_length(j,i,2);
        pre_average_anteriorthickness(i)=pre_average_anteriorthickness(i)+pre_AL_thickness(j,i,2);
        pre_average_posteriorthickness(i)=pre_average_posteriorthickness(i)+pre_PL_thickness(j,i,2);
    end
    pre_average_anteriorlength(i)=pre_average_anteriorlength(i)/(pre_cell_number_multicell(i)+pre_cell_number_monocell(i));
    pre_average_posteriorlength(i)=pre_average_posteriorlength(i)/(pre_cell_number_multicell(i)+pre_cell_number_monocell(i));
    pre_average_anteriorthickness(i)=pre_average_anteriorthickness(i)/(pre_cell_number_multicell(i)+pre_cell_number_monocell(i));
    pre_average_posteriorthickness(i)=pre_average_posteriorthickness(i)/(pre_cell_number_multicell(i)+pre_cell_number_monocell(i));
end


%% step 4
%  首先将mix_all矩阵中每种蛋白的5个时期的数据合并，保存在TriGaussian_summarize2矩阵中
TriGaussian_summarize2=zeros(2000,13,9);  % 保存9种蛋白的数据，5个时期合并
filenamearray_all2=cell(2000,9);
cell_number2=zeros(9,1);  % 保存9种蛋白的细胞数，5个时期合并
for i=1:9
    k2=0;  % k2指示mix_all2的位置，循环结束后为5个时期细胞数的总和
    for j=1:5
        for k=1:cell_number(j,i)  % k为mix_all的位置
            k2=k2+1;
            TriGaussian_summarize2(k2,:,i)=TriGaussian_summarize(k,:,j,i);
            filenamearray_all2(k2,i)=filenamearray_all(k,j,i);
        end
    end
    cell_number2(i)=k2;
end

cell_number_monocell=zeros(9,1);  % 各蛋白的lumen形成后的有单个数据文件的细胞数
filenamearray_monocell=cell(1000,9);  % 储存有单个数据文件的细胞计算的anterior-/posterior-lateral thickness对应的文件名和路径
anterior_thickness_monocell=zeros(1000,9,3);posterior_thickness_monocell=zeros(1000,9,3);  % 定义anterior-和posterior-lateral domain的cortex厚度变量，第三个维度的两层分别存放厚度、对应的长度和lumen length/cell length
cell_number_multicell=zeros(9,1);  % 各蛋白的lumen形成后的有三个数据文件的细胞数
filenamearray_multicell_anterior=cell(1000,9);  % 储存有多个数据文件的细胞计算的anterior-lateral thickness对应的文件名和路径
filenamearray_multicell_posterior=cell(1000,9);  % 储存有多个数据文件的细胞计算的posterior-lateral thickness对应的文件名和路径
anterior_length_summarize=zeros(2000,2);posterior_length_summarize=zeros(2000,2);  % 定义anterior-和posterior-lateral domain的cortex长度变量，第二个维度的两层分别存放长度和对应的lumen length/cell length值
anterior_thickness_multicell=zeros(1000,9,3);posterior_thickness_multicell=zeros(1000,9,3);  % 定义anterior-和posterior-lateral domain的cortex厚度变量，第三个维度的三层分别存放厚度、对应的长度和lumen length/cell length
k_length=0;  % 指示anterior_length_summarize赋值位置
for i=1:9
    for j=pre_cell_number(i)+1:cell_number2(i)
        filename=mat2str(cell2mat(filenamearray_all2(j,i)));  % 读取文件名
        filename=filename(2:length(filename)-1);  % 删除两侧的单引号
        S=regexp(filename,'_','split');  % 拆分文件名字符串
        month_day=mat2str(cell2mat(S(1)));  % S矩阵的第一个值为该细胞的日期(month+day)
        month_day=month_day(2:length(month_day)-1);  % 删除两侧的单引号
        framepath_filename_1=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_xiahuaxian filename_1 filename_csv];  % 拼接完整的细胞灰度值数据文件原始路径(多个文件形式)
        framepath_filename_2=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_xiahuaxian filename_2 filename_csv];
        framepath_filename_3=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_xiahuaxian filename_3 filename_csv];
        framepath_filename_danyi=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_csv];  % 拼接完整的细胞灰度值数据文件原始路径(单个文件形式)

        % 根据该细胞的收缩环特征参数描绘tri-Gaussian distribution
        rho0=TriGaussian_summarize2(j,3,i);rho_inf=TriGaussian_summarize2(j,4,i);w=TriGaussian_summarize2(j,5,i);E=TriGaussian_summarize2(j,6,i);  % 赤道轴收缩环特征参数赋值
        rho_inf_anterior=TriGaussian_summarize2(j,7,i);w_anterior=TriGaussian_summarize2(j,8,i);E_anterior=TriGaussian_summarize2(j,9,i);
        rho_inf_posterior=TriGaussian_summarize2(j,10,i);w_posterior=TriGaussian_summarize2(j,11,i);E_posterior=TriGaussian_summarize2(j,12,i);  % 两侧收缩环特征参数赋值
        f=@(x) rho0+rho_inf.*exp(-1/2.*((x-E)./w).^2)+rho_inf_anterior.*exp(-1/2.*((x-E_anterior)./w_anterior).^2)+rho_inf_posterior.*exp(-1/2.*((x-E_posterior)./w_posterior).^2);
        
        % 计算anterior-和posterior-lateral domain的cortex平均厚度和平均长度
        if exist(framepath_filename_1,'file')
            %% Step 4.1：计算有三个数据文件的细胞的anterior-和posterior-lateral domain的cortex厚度和长度，并拟合anterior/posterior-
            %%           lateral domain length/lateral-basal-lateral whole length与lumen length/cell length的函数关系
            cell_number_multicell(i)=cell_number_multicell(i)+1;  % 三个数据文件的细胞数统计
            
            data_1=xlsread(framepath_filename_1); % 导入数据
            data_2=xlsread(framepath_filename_2);
            data_3=xlsread(framepath_filename_3);
            
            % anterior-/posterior- cortex thickness和对应的lateral domain length赋值(因不同蛋白有不同的分布形式，顾分开计算)
%             a_anterior=integral(f,-1/2,-1/2+length(data_1)/(length(data_1)+length(data_2)+length(data_3)));  % anterior-lateral cortex thickness积分
%             a_anterior=a_anterior/length(data_1)*(length(data_1)+length(data_2)+length(data_3));  % 积分结果除以x轴长度，得到平均cortex thickness
            a_anterior=integral(f,-1/2,-1/2+0.05);  % anterior-lateral cortex thickness积分
            a_anterior=a_anterior/0.05;  % 积分结果除以x轴长度，得到平均cortex thickness
            anterior_thickness_multicell(cell_number_multicell(i),i,1)=a_anterior;  % cortex thickness值赋值给anterior_thickness_multifile矩阵
            anterior_thickness_multicell(cell_number_multicell(i),i,2)=length(data_1)/(length(data_1)+length(data_2)+length(data_3));  % 该细胞cortex thickness对应的anterior-lateral length
            anterior_thickness_multicell(cell_number_multicell(i),i,3)=TriGaussian_summarize2(j,end,i);  % 该细胞的lumen length/cell length
            filenamearray_multicell_anterior{cell_number_multicell(i),i}=framepath_filename_1;
            
%             a_posterior=integral(f,1/2-length(data_3)/(length(data_1)+length(data_2)+length(data_3)),1/2);  % posterior-lateral cortex thickness积分
%             a_posterior=a_posterior/length(data_3)*(length(data_1)+length(data_2)+length(data_3));  % 积分结果除以x轴长度，得到平均cortex thickness
            a_posterior=integral(f,1/2-0.05,1/2);  % posterior-lateral cortex thickness积分
            a_posterior=a_posterior/0.05;  % 积分结果除以x轴长度，得到平均cortex thickness
            posterior_thickness_multicell(cell_number_multicell(i),i,1)=a_posterior;  % cortex thickness值赋值给posterior_thickness_multifile矩阵
            posterior_thickness_multicell(cell_number_multicell(i),i,2)=length(data_3)/(length(data_1)+length(data_2)+length(data_3));  % 该细胞cortex thickness对应的posterior-lateral length
            posterior_thickness_multicell(cell_number_multicell(i),i,3)=TriGaussian_summarize2(j,end,i);  % 该细胞的lumen length/cell length
            filenamearray_multicell_posterior{cell_number_multicell(i),i}=framepath_filename_3;
            
            % anterior-/posterior-lateral length和对应的lumen length/cell length赋值(因不同蛋白不影响anterior-/posterior-lateral length，顾合并计算以增加样本量)
            logic=0;
            for k=1:k_length  % 判断该细胞是否在anterior_length矩阵中已经存在
                if anterior_length_summarize(k,2)==TriGaussian_summarize2(j,end,i)
                    logic=1;
                    break
                end
            end
            if logic==0  % 该细胞在anterior_length矩阵中不存在，则赋值计算该细胞的anterior-/posterior-lateral length
                k_length=k_length+1;
                anterior_length_summarize(k_length,1)=length(data_1)/(length(data_1)+length(data_2)+length(data_3));  % anterior-lateral length统计
                posterior_length_summarize(k_length,1)=length(data_3)/(length(data_1)+length(data_2)+length(data_3));  % posterior-lateral length统计
                anterior_length_summarize(k_length,2)=TriGaussian_summarize2(j,end,i);  % anterior-lateral length对应的lumen length/cell length
                posterior_length_summarize(k_length,2)=TriGaussian_summarize2(j,end,i);
            end
        elseif exist(framepath_filename_danyi,'file')
            %% Step 4.2：计算有一个数据文件的细胞的anterior-和posterior-lateral domain的cortex厚度，细胞三部分的拆分依据
            %%           lateral-basal-lateral whole length与lumen length/cell length
            cell_number_monocell(i)=cell_number_monocell(i)+1;  % 三个数据文件的细胞数统计
            
            % 根据文件名确定该细胞在该日期下的位置
            number_changdubi=0;
            for ii=1:30  % 外侧循环用于查找目标文件是否存在(i表示第几张photo，j表示photo中的第几个细胞，k表示荧光通道)
                filename1=num2str(ii);
                for jj=1:10
                    filename2=num2str(jj);
                    filename3='1';
                    filename_temp=[month_day filename_xiahuaxian filename1 filename2 filename3];  % 如 0928_111
                    
                    % 拼接完整的文件路径，用于下一步exist函数的查找
                    framepath_filename_1_temp=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_temp filename_xiahuaxian filename_1 filename_csv];  % 拼接完整的细胞灰度值数据文件原始路径(多个文件形式)
                    framepath_filename_danyi_temp=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_temp filename_csv];  % 如 'C:\Users\Desktop\1\0928\0928_111.csv
                    
                    if exist(framepath_filename_danyi_temp,'file')   % 一个细胞的数据有单个文件
                        number_changdubi=number_changdubi+1;
                    elseif exist(framepath_filename_1_temp,'file')
                        number_changdubi=number_changdubi+1;
                    elseif strcmp(filename_temp,filename)  % 若该张照片下，不存在该细胞的数据文件，说明该张照片中的细胞已经拟合完成，跳出j(细胞)循环
                        number_changdubi=number_changdubi+1;
                        break
                    end
                end
            end
            
            % 查找该细胞对应的lumen length和cell length (从data_changdubi中寻找)
            framepath_filename_changdubi=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_changdubi filename_xlsx];
            framepath_filename_Radius=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_Radius filename_xlsx];
            framepath_filename_ContactAngle=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_ContactAngle filename_xlsx];

            Lumen_Radius=xlsread(framepath_filename_Radius,'sheet1','A1:A100');  % 读入lumen半径
            Lumen_Radius=sqrt(Lumen_Radius./pi);
            Lumen_ContactAngle=xlsread(framepath_filename_ContactAngle,'sheet1','A1:A100');  % 读入lumen接触角
            data_lumenlength=xlsread(framepath_filename_changdubi,'sheet1','A1:A100');  % 读入该文件下的lumen长度
            data_celllength=xlsread(framepath_filename_changdubi,'sheet1','B1:B100');  % 读入该文件下的cell长度

            data=xlsread(framepath_filename_danyi); % 导入数据
            anterior_lateral_length=(data_celllength(number_changdubi)-Lumen_Radius(number_changdubi)*sin(Lumen_ContactAngle(number_changdubi)*pi./180))/2;  % anterior-lateral domain长度
            posterior_lateral_length=(data_celllength(number_changdubi)-Lumen_Radius(number_changdubi)*sin(Lumen_ContactAngle(number_changdubi)*pi./180))/2;  % posterior-lateral domain长度
            basal_length=data(end,1)-anterior_lateral_length-posterior_lateral_length;  % basal domain长度

            % anterior-/posterior- cortex thickness和对应的lateral domain length赋值(因不同蛋白有不同的分布形式，顾分开计算)
%             a_anterior=integral(f,-1/2,-1/2+anterior_lateral_length/data(end,1));  % anterior-lateral cortex thickness积分
%             a_anterior=a_anterior/(anterior_lateral_length/data(end,1));  % 积分结果除以x轴长度，得到平均cortex thickness
            a_anterior=integral(f,-1/2,-1/2+0.05);  % anterior-lateral cortex thickness积分
            a_anterior=a_anterior/0.05;  % 积分结果除以x轴长度，得到平均cortex thickness
            anterior_thickness_monocell(cell_number_monocell(i),i,1)=a_anterior;  % cortex thickness值赋值给anterior_thickness_multifile矩阵
            anterior_thickness_monocell(cell_number_monocell(i),i,2)=anterior_lateral_length/data(end,1);  % 该细胞cortex thickness对应的anterior-lateral length
            anterior_thickness_monocell(cell_number_monocell(i),i,3)=TriGaussian_summarize2(j,end,i);  % 该细胞的lumen length/cell length
            
%             a_posterior=integral(f,1/2-posterior_lateral_length/data(end,1),1/2);  % posterior-lateral cortex thickness积分
%             a_posterior=a_posterior/(posterior_lateral_length/data(end,1));  % 积分结果除以x轴长度，得到平均cortex thickness
            a_posterior=integral(f,1/2-0.05,1/2);  % posterior-lateral cortex thickness积分
            a_posterior=a_posterior/0.05;  % 积分结果除以x轴长度，得到平均cortex thickness
            posterior_thickness_monocell(cell_number_monocell(i),i,1)=a_posterior;  % cortex thickness值赋值给posterior_thickness_multifile矩阵
            posterior_thickness_monocell(cell_number_monocell(i),i,2)=posterior_lateral_length/data(end,1);  % 该细胞cortex thickness对应的posterior-lateral length
            posterior_thickness_monocell(cell_number_monocell(i),i,3)=TriGaussian_summarize2(j,end,i);  % 该细胞的lumen length/cell length
            filenamearray_monocell{cell_number_monocell(i),i}=framepath_filename_danyi;

            % anterior-/posterior-lateral length和对应的lumen length/cell length赋值(因不同蛋白不影响anterior-/posterior-lateral length，顾合并计算以增加样本量)
            logic=0;
            for k=1:k_length  % 判断该细胞是否在anterior_length矩阵中已经存在
                if anterior_length_summarize(k,2)==TriGaussian_summarize2(j,end,i)
                    logic=1;
                    break
                end
            end
            if logic==0  % 该细胞在anterior_length矩阵中不存在，则赋值计算该细胞的anterior-/posterior-lateral length
                k_length=k_length+1;
                anterior_length_summarize(k_length,1)=anterior_lateral_length/data(end,1);  % anterior-lateral length统计
                posterior_length_summarize(k_length,1)=posterior_lateral_length/data(end,1);  % posterior-lateral length统计
                anterior_length_summarize(k_length,2)=TriGaussian_summarize2(j,end,i);  % anterior-lateral length对应的lumen length/cell length
                posterior_length_summarize(k_length,2)=TriGaussian_summarize2(j,end,i);
            end
        end
    end
end

% 对anterior/posterior_length矩阵按照lumen length/cell length进行排序(冒泡排序)
for i=k_length:-1:2
    % anterior_length矩阵排序
    for j=2:i
        if anterior_length_summarize(j,2)<anterior_length_summarize(j-1,2)
            temp2=anterior_length_summarize(j,1);
            anterior_length_summarize(j,1)=anterior_length_summarize(j-1,1);
            anterior_length_summarize(j-1,1)=temp2;
            temp=anterior_length_summarize(j,2);
            anterior_length_summarize(j,2)=anterior_length_summarize(j-1,2);
            anterior_length_summarize(j-1,2)=temp;
        end
    end
    
    % posterior_length矩阵排序
    for j=2:i
        if posterior_length_summarize(j,2)<posterior_length_summarize(j-1,2)
            temp2=posterior_length_summarize(j,1);
            posterior_length_summarize(j,1)=posterior_length_summarize(j-1,1);
            posterior_length_summarize(j-1,1)=temp2;
            temp=posterior_length_summarize(j,2);
            posterior_length_summarize(j,2)=posterior_length_summarize(j-1,2);
            posterior_length_summarize(j-1,2)=temp;
        end
    end
end

% 拟合anterior/posterior-lateral domain length/lateral-basal-lateral whole length与t的函数关系
% 找出anterior/posterior_length_summarize矩阵中lumen length/cell length对应的时间t
% 根据实际观测拟合transverse diameter(TD)和longitudinal radius(LR)随时间的变化函数
TD=[3.248130325	3.697268921	4.042205695	4.875820261	5.623190217	5.105764223	5.795637773	5.996833531	6.456763119	6.571755932	6.485511322	6.370518509	6.399266712	6.485511322	6.370518509	6.600504135	6.284315564	6.226819157	6.226819157	6.830448097	7.37658063	7.577818053	8.411390955	8.267691603	8.900068746	9.474949482	10.02108202	10.4235152	11.14213695	11.77451409	12.57938046	12.86682082	13.38420515	13.64293898	14.21781972	14.36151907	14.5915047	15.10888903	15.14842822	15.05139262	16];
LR=[0.8114597 0.88031998 1.012697123 0.632377143 0.966522926 1.167760348 1.523446451 1.943878508 2.08757786	2.288815282	2.202570673	2.949940629	3.381080349	3.668520717	4.214694915	4.818323855	4.703331042	5.174051622	5.192008833	5.651938421	5.795637773	5.968126992	6.054329938	6.140574547	6.169322751	5.910630586	5.853134179	6.111826344	6.284315564	5.996833531	6.083078141	5.996833531	5.824385976	5.863883507	5.853134179	5.853134179	5.795637773	5.939378789	5.738141366	6.111826344	5.996833531];
TD2=[TD 16.1 16.1 16.1 16.1 16.1 16.1 16.1 16.1 16.1 16.1];
t=linspace(0,120,length(TD))+6;  % 观测数据的时间t
t2=[t t(end) t(end) t(end) t(end) t(end) t(end) t(end) t(end) t(end) t(end)];
tlength=1000;  % 拟合时间轴的精度
t_fit=linspace(0,126,tlength);  % 拟合时间轴
N1=5;N2=6;
TD_fit=zeros(1,length(t_fit));
LR_fit=zeros(1,length(t_fit));
p=polyfit(t2,TD2,N1);
p2=polyfit(t,LR,N2);
for k=1:N1+1
    TD_fit=TD_fit+p(k).*t_fit.^(N1-k+1);
end
for k=1:N2+1
    LR_fit=LR_fit+p2(k).*t_fit.^(N2-k+1);
end
TD_fit(1)=0.01;TD_fit(2)=0.02;
% figure
% plot(t,TD,'.','Color','b','MarkerSize',20)
% hold on
% plot(t,LR,'.','Color','r','MarkerSize',20)
% plot(t_fit,TD_fit,'Color','b','LineWidth',3)
% plot(t_fit,LR_fit,'Color','r','LineWidth',3)

% 根据TD和LR与R和theta的关系，通过联立方程组求出R和theta随时间的变化
R=zeros(1,length(t_fit));
theta=zeros(1,length(t_fit));
for i=1:length(t_fit)
    wucha_best=1000;  % 最小误差值
    for theta_search=0:pi/1000:pi
        R_search=TD_fit(i)/2/sin(theta_search);
        if abs(R_search+R_search*cos(theta_search)-LR_fit(i))<wucha_best
            wucha_best=abs(R_search+R_search*cos(theta_search)-LR_fit(i));
            R(i)=R_search;
            theta(i)=theta_search*180/pi;
        end
    end
end


% 将lumen length/cell length(lateral_tension_all(:,:,2))转化为对应的时间t，储存在t_rectify矩阵中
t_rectify=zeros(1,tlength);  % 矫正文献观测数据后的时间t
t_fit_rectify=linspace(0,126,tlength);  % 矫正文献观测数据后的拟合曲线横坐标t
TD_rectify=zeros(1,tlength);  % 矫正后的TD
TD_fit_rectify=zeros(1,tlength);  % 矫正后的TD拟合函数
TD_jump=6.2;  % lumen length/cell length的跳跃值

i_rectify_position=0;  % 赋值到rectify矩阵中的位置
for i=1:tlength  % 去除TD_jump周围若干数据
    if TD_fit(i)>=TD_jump+1 || TD_fit(i)<=TD_jump-1
        i_rectify_position=i_rectify_position+1;
        t_rectify(i_rectify_position)=t_fit(i);
        TD_rectify(i_rectify_position)=TD_fit(i);
    end
end
% 对TD_rectify进行拟合
% N_rectify=6;
% p=polyfit(t_rectify(1:i_rectify_position),TD_rectify(1:i_rectify_position),N_rectify);
% for k=1:N_rectify+1
%     TD_fit_rectify=TD_fit_rectify+p(k).*t_fit_rectify.^(N_rectify-k+1);
% end
% 对TD_rectify进行插值
TD_fit_rectify=interp1(t_rectify(1:i_rectify_position),TD_rectify(1:i_rectify_position),t_fit_rectify,'pchip');  % 三次样条插值

tlength=1000;
t_observation_AL=zeros(1000,1);  % 实验观测的lumen length/cell length对应的时间t
t_observation_PL=zeros(1000,1);  % 实验观测的lumen length/cell length对应的时间t
for i=1:k_length
    wucha_best=1000;
    for j=1:tlength
        wucha=abs(anterior_length_summarize(i,2)*16-TD_fit_rectify(j));
        if wucha<wucha_best
            wucha_best=wucha;
            t_observation_AL(i)=t_fit_rectify(j);
        end
    end

    wucha_best=1000;
    for j=1:tlength
        wucha=abs(posterior_length_summarize(i,2)*16-TD_fit_rectify(j));
        if wucha<wucha_best
            wucha_best=wucha;
            t_observation_PL(i)=t_fit_rectify(j);
        end
    end
end

N=3;  % 拟合阶数
lateral_lumen_length_fit=zeros(N+1,2);  % 储存拟合多项式系数变量，第一层为anterior的函数关系，第二层为posterior的函数关系
p_fit=polyfit(t_observation_AL(1:k_length),anterior_length_summarize(1:k_length,1),N);  % anterior-lateral domain length的多项式拟合
for k=1:N+1
    lateral_lumen_length_fit(k,1)=p_fit(k);
end
p_fit=polyfit(t_observation_PL(1:k_length),posterior_length_summarize(1:k_length,1),N);  % posterior-lateral domain length的多项式拟合
for k=1:N+1
    lateral_lumen_length_fit(k,2)=p_fit(k);
end

% x_fit=linspace(min(t_observation_AL(1:k_length)),max(t_observation_AL(1:k_length)),1000);
% y_fit=0;y_fit2=0;
% for k=1:N+1
%     y_fit=y_fit+lateral_lumen_length_fit(k,1).*x_fit.^(N-k+1);
%     y_fit2=y_fit2+lateral_lumen_length_fit(k,2).*x_fit.^(N-k+1);
% end
% plot(t_observation_PL(1:k_length),anterior_length_summarize(1:k_length,1))
% hold on
% plot(x_fit,y_fit)
% plot(t_observation_PL(1:k_length),posterior_length_summarize(1:k_length,1))
% hold on
% plot(x_fit,y_fit2)





%% step 5
%% Step 5.1：首先计算两侧都有细胞灰度值数据的lumen所受到的来自脊索细胞两侧收缩环的收缩力大小，收缩力大小等于前一个细
%%           胞的posterior-lateral domain的cortex厚度、长度和单位厚度产生的压强的乘积与后一个细胞的anterior-lateral 
%%           domain的cortex厚度、长度和单位厚度产生的压强的乘积的和减去基础厚度对应的收缩力
lateral_tension_situation1=zeros(1000,9,2);  % 储存两侧都有细胞灰度值数据的lumen所受到的来自脊索细胞两侧收缩环的收缩力值，第二个维度储存每个lumen对应的lumen length/cell length
filenamearray_lateraltension_situation1=cell(1000,9,2);  % 储存计算每一个lateral_tension对应的两个细胞的路径文件名
i_tension_position=zeros(9,1);  % 指示lateral_tension_situation1矩阵的写入位置
for i=1:9
    for j=1:cell_number_multicell(i)
        framepath_filename_1=mat2str(cell2mat(filenamearray_multicell_anterior(j,i)));  % 读取anterior-lateral domain的路径和文件名
        framepath_filename_1=framepath_filename_1(2:length(framepath_filename_1)-1);  % 删除两侧的单引号
        

        % 提取路径和文件名中包含的该文件的日期(month+day)、照片数和细胞数
        S=regexp(framepath_filename_1,'\','split');  % 拆分文件名字符串
        anterior_date=mat2str(cell2mat(S(7)));  % 提取路径文件名中的日期
        anterior_date=anterior_date(2:length(anterior_date)-1);  % 删除两侧的单引号
        filename=mat2str(cell2mat(S(end)));  % 提取路径文件名中的完整文件名
        if length(filename)==17  % filename长度为17表示照片数为两位数
            anterior_photo_number=filename(7:8);
            anterior_cell_number=filename(9);
        elseif length(filename)==16  % filename长度为16表示照片数为一位数           
            anterior_photo_number=filename(7);
            anterior_cell_number=filename(8);
        end
        
        % 查找filenamearray_multicell_posterior矩阵中有无与该anterior_lateral domain所属的细胞对应的posterior_lateral domain
        logic_filenamematch=0;  % 是否存在对应细胞的逻辑值
        for k=1:cell_number_multicell(i)
            framepath_filename_3=mat2str(cell2mat(filenamearray_multicell_posterior(k,i)));  % 读取anterior-lateral domain的路径和文件名
            framepath_filename_3=framepath_filename_3(2:length(framepath_filename_3)-1);  % 删除两侧的单引号
            
            % 提取路径和文件名中包含的该文件的日期(month+day)、照片数和细胞数
            S=regexp(framepath_filename_3,'\','split');  % 拆分文件名字符串
            posterior_date=mat2str(cell2mat(S(7)));  % 提取路径文件名中的日期
            posterior_date=posterior_date(2:length(posterior_date)-1);  % 删除两侧的单引号
            filename=mat2str(cell2mat(S(end)));  % 提取路径文件名中的完整文件名
            if length(filename)==17  % filename长度为17表示照片数为两位数
                posterior_photo_number=filename(7:8);
                posterior_cell_number=filename(9);
            elseif length(filename)==16  % filename长度为16表示照片数为一位数
                posterior_photo_number=filename(7);
                posterior_cell_number=filename(8);
            end
            
            % 判断日期(month+day)、照片数和细胞数是否对应
            if ~sum(abs(anterior_date-posterior_date)) && ~sum(abs(anterior_photo_number-posterior_photo_number))  % 日期和照片数是否相同
                if anterior_cell_number==posterior_cell_number+1  % anterior part是否是posterior part的前一个细胞
                    logic_filenamematch=1;
                    break
                end
            end
        end
        
        % 计算该lumen所受到的来自两侧lateral domain的收缩力合力
        if logic_filenamematch==1  % 表示有匹配的一组文件名
            i_tension_position(i)=i_tension_position(i)+1;
            lateral_tension_situation1(i_tension_position(i),i,1)=(anterior_thickness_multicell(j,i,1)-pre_average_anteriorthickness(i))*anterior_thickness_multicell(j,i,2)+(posterior_thickness_multicell(k,i,1)-pre_average_posteriorthickness(i))*posterior_thickness_multicell(k,i,2);  % 计算收缩力并赋值
            lateral_tension_situation1(i_tension_position(i),i,1)=lateral_tension_situation1(i_tension_position(i),i,1)/(pre_average_anteriorthickness(i)*pre_average_anteriorlength(i)+pre_average_posteriorthickness(i)*pre_average_posteriorlength(i))*2;  % 得到相对厚度的收缩力
%             lateral_tension_situation1(i_tension_position(i),i,1)=(anterior_thickness_multicell(j,i,1)-pre_average_anteriorthickness(i))+(posterior_thickness_multicell(k,i,1)-pre_average_posteriorthickness(i));  % 计算收缩力并赋值
%             lateral_tension_situation1(i_tension_position(i),i,1)=lateral_tension_situation1(i_tension_position(i),i,1)/(pre_average_anteriorthickness(i)+pre_average_posteriorthickness(i))*2;  % 得到相对厚度的收缩力
            if lateral_tension_situation1(i_tension_position(i),i,1)<0
                lateral_tension_situation1(i_tension_position(i),i,1)=0;
            end
            lateral_tension_situation1(i_tension_position(i),i,2)=anterior_thickness_multicell(j,i,3);  % 对应的lumen length/cell length赋值
            filenamearray_lateraltension_situation1{i_tension_position(i),i,1}=framepath_filename_1;
            filenamearray_lateraltension_situation1{i_tension_position(i),i,2}=framepath_filename_3;
        end
    end
end




%% Step 5.2：再计算仅有单侧有细胞灰度值数据的lumen所受到的来自脊索细胞两侧收缩环的收缩力大小，收缩力大小等于该细胞自
%%           身的anterior-lateral domain的cortex厚度、长度和单位厚度产生的压强的乘积和posterior-lateral domain的cor-
%%           tex厚度、长度和单位厚度产生的压强的乘积的和减去基础厚度对应的收缩力
lateral_tension_situation2=zeros(1000,9,2);  % 储存一侧有细胞灰度值数据的lumen所受到的来自脊索细胞两侧收缩环的收缩力值，第二个维度储存每个lumen对应的lumen length/cell length
filenamearray_lateraltension_situation2=cell(1000,9);  % 储存计算每一个lateral_force对应的两个细胞的路径文件名
i_tension2_position=zeros(9,1);  % 指示lateral_tension_situation2矩阵的写入位置

% 计算有三个数据文件的细胞中，未计算过lateral_tension的细胞所对应的lumen受到的收缩力值，并减去基础厚度对应的收缩力
for i=1:9
    % 计算anterior-/posterior-lateral domain所对应的还未计算过的细胞
    for j=1:cell_number_multicell(i)
        framepath_filename_1=mat2str(cell2mat(filenamearray_multicell_anterior(j,i)));  % 读取anterior-lateral domain的路径和文件名
        framepath_filename_1=framepath_filename_1(2:length(framepath_filename_1)-1);  % 删除两侧的单引号
        
        % 判断该文件名是否在filenamearray_lateraltension_situation1矩阵中不存在
        logic_filenamematch=1;  % 文件名是否存在的逻辑值
        for k=1:i_tension_position(i)
            if length(framepath_filename_1)==length(filenamearray_lateraltension_situation1{k,i,1}) && ~sum(abs(framepath_filename_1-filenamearray_lateraltension_situation1{k,i,1}))
                logic_filenamematch=0;
                break
            end
        end
        
        % 计算该lumen所受到的来自两侧lateral domain的收缩力合力
        if logic_filenamematch==1
            i_tension2_position(i)=i_tension2_position(i)+1;
            lateral_tension_situation2(i_tension2_position(i),i,1)=(anterior_thickness_multicell(j,i,1)-pre_average_anteriorthickness(i))*anterior_thickness_multicell(j,i,2)+(posterior_thickness_multicell(j,i,1)-pre_average_posteriorthickness(i))*posterior_thickness_multicell(j,i,2);  % 计算收缩力并赋值
            lateral_tension_situation2(i_tension2_position(i),i,1)=lateral_tension_situation2(i_tension2_position(i),i,1)/(pre_average_anteriorthickness(i)*pre_average_anteriorlength(i)+pre_average_posteriorthickness(i)*pre_average_posteriorlength(i))*2;  % 得到相对厚度的收缩力
%             lateral_tension_situation2(i_tension2_position(i),i,1)=(anterior_thickness_multicell(j,i,1)-pre_average_anteriorthickness(i))+(posterior_thickness_multicell(j,i,1)-pre_average_posteriorthickness(i));  % 计算收缩力并赋值
%             lateral_tension_situation2(i_tension2_position(i),i,1)=lateral_tension_situation2(i_tension2_position(i),i,1)/(pre_average_anteriorthickness(i)+pre_average_posteriorthickness(i))*2;  % 得到相对厚度的收缩力
            if lateral_tension_situation2(i_tension2_position(i),i,1)<0
                lateral_tension_situation2(i_tension2_position(i),i,1)=0;
            end
            lateral_tension_situation2(i_tension2_position(i),i,2)=anterior_thickness_multicell(j,i,3);  % 对应的lumen length/cell length赋值
            filenamearray_lateraltension_situation2{i_tension2_position(i),i}=framepath_filename_1;
        end
    end
end

% 计算有单个数据文件的细胞所对应的lumen受到的收缩力值，并减去基础厚度对应的收缩力
for i=1:9
    for j=1:cell_number_monocell(i)
        framepath_filename_danyi=mat2str(cell2mat(filenamearray_monocell(j,i)));  % 读取anterior-lateral domain的路径和文件名
        framepath_filename_danyi=framepath_filename_danyi(2:length(framepath_filename_danyi)-1);  % 删除两侧的单引号
        
        i_tension2_position(i)=i_tension2_position(i)+1;
        lateral_tension_situation2(i_tension2_position(i),i,1)=(anterior_thickness_monocell(j,i,1)-pre_average_anteriorthickness(i))*anterior_thickness_monocell(j,i,2)+(posterior_thickness_monocell(j,i,1)-pre_average_posteriorthickness(i))*posterior_thickness_monocell(j,i,2);  % 计算收缩力并赋值
        lateral_tension_situation2(i_tension2_position(i),i,1)=lateral_tension_situation2(i_tension2_position(i),i,1)/(pre_average_anteriorthickness(i)*pre_average_anteriorlength(i)+pre_average_posteriorthickness(i)*pre_average_posteriorlength(i))*2;  % 得到相对厚度的收缩力
        if lateral_tension_situation2(i_tension2_position(i),i,1)<0
            lateral_tension_situation2(i_tension2_position(i),i,1)=0;
        end
        lateral_tension_situation2(i_tension2_position(i),i,2)=anterior_thickness_monocell(j,i,3);  % 对应的lumen length/cell length赋值
        filenamearray_lateraltension_situation2{i_tension2_position(i),i}=framepath_filename_danyi;
    end
end



%% Step 5.3：拟合lumen所受到的来自脊索细胞两侧收缩环的收缩力按lumen length/cell length大小进行排序
% 将两侧都有细胞数据的lumen和单侧有细胞数据的lumen收缩力结果合并
lateral_tension_all=zeros(1000,9,2);  % 收缩力计算结果合并矩阵
for i=1:9
    for j=1:i_tension_position(i)
        lateral_tension_all(j,i,:)=lateral_tension_situation1(j,i,:);
    end
    for j=1:i_tension2_position(i)
        lateral_tension_all(j+i_tension_position(i),i,:)=lateral_tension_situation2(j,i,:);
    end
end
lateral_tension_all(:,:,1)=lateral_tension_all(:,:,1)./2;

% 对lateral_tension_all矩阵按lumen length/cell length大小进行排序(冒泡排序法)
for i=1:9
    for j=i_tension_position(i)+i_tension2_position(i):-1:2
        for k=2:j
            if lateral_tension_all(k,i,2)<lateral_tension_all(k-1,i,2)
                temp2=lateral_tension_all(k,i,1);
                lateral_tension_all(k,i,1)=lateral_tension_all(k-1,i,1);
                lateral_tension_all(k-1,i,1)=temp2;
                temp=lateral_tension_all(k,i,2);
                lateral_tension_all(k,i,2)=lateral_tension_all(k-1,i,2);
                lateral_tension_all(k-1,i,2)=temp;
            end
        end
    end
end



%% Step 5.4：根据归一化的TD(lumen length/cell length)与时间t的观测数据拟合变化关系，拟合lumen所受到的来自脊索细胞两侧
%%           收缩环的收缩力与时间t的函数关系
% 根据观测获得R和theta随时间的实际变化关系
% 根据实际观测拟合transverse diameter(TD)和longitudinal radius(LR)随时间的变化函数
TD=[3.248130325	3.697268921	4.042205695	4.875820261	5.623190217	5.105764223	5.795637773	5.996833531	6.456763119	6.571755932	6.485511322	6.370518509	6.399266712	6.485511322	6.370518509	6.600504135	6.284315564	6.226819157	6.226819157	6.830448097	7.37658063	7.577818053	8.411390955	8.267691603	8.900068746	9.474949482	10.02108202	10.4235152	11.14213695	11.77451409	12.57938046	12.86682082	13.38420515	13.64293898	14.21781972	14.36151907	14.5915047	15.10888903	15.14842822	15.05139262	16];
LR=[0.8114597 0.88031998 1.012697123 0.632377143 0.966522926 1.167760348 1.523446451 1.943878508 2.08757786	2.288815282	2.202570673	2.949940629	3.381080349	3.668520717	4.214694915	4.818323855	4.703331042	5.174051622	5.192008833	5.651938421	5.795637773	5.968126992	6.054329938	6.140574547	6.169322751	5.910630586	5.853134179	6.111826344	6.284315564	5.996833531	6.083078141	5.996833531	5.824385976	5.863883507	5.853134179	5.853134179	5.795637773	5.939378789	5.738141366	6.111826344	5.996833531];
TD2=[TD 16.1 16.1 16.1 16.1 16.1 16.1 16.1 16.1 16.1 16.1];
t=linspace(0,120,length(TD))+6;  % 观测数据的时间t
t2=[t t(end) t(end) t(end) t(end) t(end) t(end) t(end) t(end) t(end) t(end)];
tlength=1000;  % 拟合时间轴的精度
t_fit=linspace(0,126,tlength);  % 拟合时间轴
N1=5;N2=6;
TD_fit=zeros(1,length(t_fit));
LR_fit=zeros(1,length(t_fit));
p=polyfit(t2,TD2,N1);
p2=polyfit(t,LR,N2);
for k=1:N1+1
    TD_fit=TD_fit+p(k).*t_fit.^(N1-k+1);
end
for k=1:N2+1
    LR_fit=LR_fit+p2(k).*t_fit.^(N2-k+1);
end
TD_fit(1)=0.01;TD_fit(2)=0.02;
% figure
% plot(t,TD,'.','Color','b','MarkerSize',20)
% hold on
% plot(t,LR,'.','Color','r','MarkerSize',20)
% plot(t_fit,TD_fit,'Color','b','LineWidth',3)
% plot(t_fit,LR_fit,'Color','r','LineWidth',3)

% 根据TD和LR与R和theta的关系，通过联立方程组求出R和theta随时间的变化
R=zeros(1,length(t_fit));
theta=zeros(1,length(t_fit));
for i=1:length(t_fit)
    wucha_best=1000;  % 最小误差值
    for theta_search=0:pi/1000:pi
        R_search=TD_fit(i)/2/sin(theta_search);
        if abs(R_search+R_search*cos(theta_search)-LR_fit(i))<wucha_best
            wucha_best=abs(R_search+R_search*cos(theta_search)-LR_fit(i));
            R(i)=R_search;
            theta(i)=theta_search*180/pi;
        end
    end
end


% 将lumen length/cell length(lateral_tension_all(:,:,2))转化为对应的时间t，储存在t_rectify矩阵中
t_rectify=zeros(1,tlength);  % 矫正文献观测数据后的时间t
t_fit_rectify=linspace(0,126,tlength);  % 矫正文献观测数据后的拟合曲线横坐标t
TD_rectify=zeros(1,tlength);  % 矫正后的TD
TD_fit_rectify=zeros(1,tlength);  % 矫正后的TD拟合函数
TD_jump=6.2;  % lumen length/cell length的跳跃值

i_rectify_position=0;  % 赋值到rectify矩阵中的位置
for i=1:tlength  % 去除TD_jump周围若干数据
    if TD_fit(i)>=TD_jump+1 || TD_fit(i)<=TD_jump-1
        i_rectify_position=i_rectify_position+1;
        t_rectify(i_rectify_position)=t_fit(i);
        TD_rectify(i_rectify_position)=TD_fit(i);
    end
end
% 对TD_rectify进行拟合
% N_rectify=6;
% p=polyfit(t_rectify(1:i_rectify_position),TD_rectify(1:i_rectify_position),N_rectify);
% for k=1:N_rectify+1
%     TD_fit_rectify=TD_fit_rectify+p(k).*t_fit_rectify.^(N_rectify-k+1);
% end
% 对TD_rectify进行插值
TD_fit_rectify=interp1(t_rectify(1:i_rectify_position),TD_rectify(1:i_rectify_position),t_fit_rectify,'pchip');  % 三次样条插值

% plot(t_rectify(1:i_rectify_position),TD_rectify(1:i_rectify_position),'.','MarkerSize',20)
% hold on
% plot(t_fit_rectify,TD_fit_rectify,'LineWidth',3)


% 找出lateral_tension_all矩阵中lumen length/cell length对应的时间t
t_observation=zeros(1000,1);  % 实验观测的lumen length/cell length对应的时间t
for i=1:i_tension_position(1)+i_tension2_position(1)
    wucha_best=1000;
    for j=1:tlength
        wucha=abs(lateral_tension_all(i,1,2)*16-TD_fit_rectify(j));
        if wucha<wucha_best
            wucha_best=wucha;
            t_observation(i)=t_fit_rectify(j);
        end
    end
end


% 拟合lumen所受到的来自脊索细胞两侧收缩环的收缩力与时间t的函数关系
N=3;  % 拟合阶数
lateral_tension_fit=zeros(N+1,9);  % 储存拟合多项式系数变量
delta_tensionfit=zeros(1000,9);  % 收缩压强拟合的标准差
i=1;
[p_fit,S]=polyfit(t_observation(1:i_tension_position(i)+i_tension2_position(i)),lateral_tension_all(1:i_tension_position(i)+i_tension2_position(i),i,1),N);  % anterior-lateral domain length的多项式拟合
[y_polyval,delta_tensionfit(1:i_tension_position(i)+i_tension2_position(i),i)]=polyval(p_fit,t_observation(1:i_tension_position(i)+i_tension2_position(i)),S);
for k=1:N+1
    lateral_tension_fit(k,i)=p_fit(k);
end



%% 计算Lumen Contact Angle和lateral tension的相关性
win=21;  % 滑动平均值的跨度
ContactAngle_LateralTension_MA=zeros(cell_number_match,3);  % 储存lateral tension的滑动平均值
ContactAngle_LateralTension_Mstd=zeros(cell_number_match,3);  % 储存lateral tension的滑动标准差
ContactAngle_LateralTension_MA2=zeros(cell_number_match,3);  % 储存lateral tension的二次滑动平均值
ContactAngle_LateralTension_Mstd2=zeros(cell_number_match,3);  % 储存lateral tension的二次滑动标准差

for i=1:cell_number_match
    for j=1:2
        if i<=(win-1)/2  % 若t为初始若干时刻
            ContactAngle_LateralTension_MA(i,j)=sum(ContactAngle_LateralTension(1:2*i-1,j))/(2*i-1);
            ContactAngle_LateralTension_Mstd(i,j)=std(ContactAngle_LateralTension(1:2*i-1,j));
        elseif i<i_tension_position(1)+i_tension2_position(1)-(win-1)/2+1  % 若t为中间时刻
            ContactAngle_LateralTension_MA(i,j)=sum(ContactAngle_LateralTension(i-(win-1)/2:i+(win-1)/2,j))/win;
            ContactAngle_LateralTension_Mstd(i,j)=std(ContactAngle_LateralTension(i-(win-1)/2:i+(win-1)/2,j));
        else  % 若t为末尾若干时刻
            temp=i_tension_position(1)+i_tension2_position(1)-i+1;
            ContactAngle_LateralTension_MA(i,j)=sum(ContactAngle_LateralTension(cell_number_match-2*temp+2:cell_number_match,j))/(2*temp-1);
            ContactAngle_LateralTension_Mstd(i,j)=std(ContactAngle_LateralTension(cell_number_match-2*temp+2:cell_number_match,j));
        end
    end
end
tension_sort=ContactAngle_LateralTension_MA(1:cell_number_match,2);
[tension_sort,I]=sort(tension_sort);
ContactAngle_LateralTension_MA(1:cell_number_match,:)=ContactAngle_LateralTension_MA(I,:);


for i=1:cell_number_match
    for j=1:2
        if i<=(win-1)/2  % 若t为初始若干时刻
            ContactAngle_LateralTension_MA2(i,j)=sum(ContactAngle_LateralTension_MA(1:2*i-1,j))/(2*i-1);
            ContactAngle_LateralTension_Mstd2(i,j)=std(ContactAngle_LateralTension_MA(1:2*i-1,j));
        elseif i<i_tension_position(1)+i_tension2_position(1)-(win-1)/2+1  % 若t为中间时刻
            ContactAngle_LateralTension_MA2(i,j)=sum(ContactAngle_LateralTension_MA(i-(win-1)/2:i+(win-1)/2,j))/win;
            ContactAngle_LateralTension_Mstd2(i,j)=std(ContactAngle_LateralTension_MA(i-(win-1)/2:i+(win-1)/2,j));
        else  % 若t为末尾若干时刻
            temp=i_tension_position(1)+i_tension2_position(1)-i+1;
            ContactAngle_LateralTension_MA2(i,j)=sum(ContactAngle_LateralTension_MA(cell_number_match-2*temp+2:cell_number_match,j))/(2*temp-1);
            ContactAngle_LateralTension_Mstd2(i,j)=std(ContactAngle_LateralTension_MA(cell_number_match-2*temp+2:cell_number_match,j));
        end
    end
end
ContactAngle_LateralTension_MA2=ContactAngle_LateralTension_MA2((win-1)/2+1:cell_number_match-(win-1)/2,:);
ContactAngle_LateralTension_Mstd2=ContactAngle_LateralTension_Mstd2((win-1)/2+1:cell_number_match-(win-1)/2,:);
cell_number_match=cell_number_match-(win-1);

% 拟合Lumen Contact Angle和lateral tension的线性关系
N=1;
[p,S]=polyfit(ContactAngle_LateralTension_MA2(1:cell_number_match,2),ContactAngle_LateralTension_MA2(1:cell_number_match,1),N); 
x_fit=linspace(ContactAngle_LateralTension_MA2(1,2)-0.2,ContactAngle_LateralTension_MA2(cell_number_match,2)+0.2,cell_number_match);
[y_fit,delta]=polyval(p,x_fit,S);

% 计算Lumen Contact Angle和lateral tension的拟合优度
x_fit_interp=linspace(x_fit(1),x_fit(end),10000);
y_fit_interp=interp1(x_fit,y_fit,x_fit_interp,'linear');  % 三次样条插值
sum_raw=sum((ContactAngle_LateralTension_MA2(1:cell_number_match,1)-sum(ContactAngle_LateralTension_MA2(1:cell_number_match,1))/cell_number_match).^2);
sum_fit=0;
for i=1:cell_number_match
    wucha_best=10000;
    for j=1:length(x_fit_interp)
        wucha=abs(x_fit_interp(j)-ContactAngle_LateralTension_MA2(i,2));

        if wucha<wucha_best
            j_best=j;
            wucha_best=wucha;
        end
    end
    sum_fit=(y_fit_interp(j_best)-sum(ContactAngle_LateralTension_MA2(1:cell_number_match,1))/cell_number_match)^2+sum_fit;
end
r_ContactAngle_LateralTension=sqrt(sum_fit/sum_raw);

% 计算Lumen Contact Angle和lateral tension的相关系数
A=corrcoef(ContactAngle_LateralTension_MA2(1:cell_number_match,1),ContactAngle_LateralTension_MA2(1:cell_number_match,2));
ContactAngle_LateralTension_time_corrcoef=A(1,2);

