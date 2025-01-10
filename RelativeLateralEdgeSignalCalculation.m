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

filename_changdubi='���ȱ�';
filename_Radius='�뾶 (Radius, �������Ϊlumen���)';
filename_ContactAngle='�Ӵ��� (Contact Angle, �������ΪԲ�ܽ�)';

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



%%% ��ȡ�洢�����(��)��˹�ֲ�����������mix_all�����ļ�
for a=1
    % ��ȡmix_lifeact��������
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
    
    % ��ȡmix_actin��������
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
    
    % ��ȡmix_MLC��������
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
    
    % ��ȡmix_RhoA��������
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
    
    % ��ȡmix_RBD��������
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
    
    % ��ȡmix_VD1��������
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
    
    % ��ȡmix_talinA��������
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
    
    % ��ȡmix_Cdc42WT��������
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
    
    % ��ȡmix_Cdc42D118A��������
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
    
    
    % ������ʱ�ڡ����ֵ��׵����ݸ�ֵ��mix_all����
    mix_MultiGaussian_all=zeros(1000,25,5,9); % ����ĳ��ʱ�ڡ�ĳ�ֵ����µ�wucha_best filename data 1-25 changdubi
    mix_EquivalentTriGaussian_all=zeros(1000,13,5,9); % ����ĳ��ʱ�ڡ�ĳ�ֵ����µ�wucha_best filename data 1-25 changdubi
    mix_TriGaussian_all=zeros(1000,13,5,9); % ����ĳ��ʱ�ڡ�ĳ�ֵ����µ�wucha_best filename data 1-25 changdubi
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


%%% ��ȡ�����������ļ��б�����ļ��������浽filenamearray_all�����У����ں������Ҹ�ϸ����ԭʼ�Ҷ�ֵ�ļ�
filenamearray_all=cell(1000,5,9);  % lumen�γ�ǰ���ļ�������

% �ļ�·��ƴ��
framepath_filename_lifeact=[framepath filename_lifeact filename_xlsx];
framepath_filename_actin=[framepath filename_actin filename_xlsx];
framepath_filename_MLC=[framepath filename_MLC filename_xlsx];
framepath_filename_RhoA=[framepath filename_RhoA filename_xlsx];
framepath_filename_RBD=[framepath filename_RBD filename_xlsx];
framepath_filename_VD1=[framepath filename_VD1 filename_xlsx];
framepath_filename_talinA=[framepath filename_talinA filename_xlsx];
framepath_filename_Cdc42WT=[framepath filename_Cdc42WT filename_xlsx];
framepath_filename_Cdc42D118A=[framepath filename_Cdc42D118A filename_xlsx];

% �ļ�����ȡ�����浽filenamearray_all����
for aa=1
    % ��ȡlifeact��Ӧ���ļ���
    [~, ~, filenamearray_all(:,1,1)] = xlsread(framepath_filename_lifeact,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,1)] = xlsread(framepath_filename_lifeact,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,1)] = xlsread(framepath_filename_lifeact,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,1)] = xlsread(framepath_filename_lifeact,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,1)] = xlsread(framepath_filename_lifeact,'sheet5','A2:A1001');
    
    % ��ȡactin��Ӧ���ļ���
    [~, ~, filenamearray_all(:,1,2)] = xlsread(framepath_filename_actin,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,2)] = xlsread(framepath_filename_actin,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,2)] = xlsread(framepath_filename_actin,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,2)] = xlsread(framepath_filename_actin,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,2)] = xlsread(framepath_filename_actin,'sheet5','A2:A1001');
    
    % ��ȡMLC��Ӧ���ļ���
    [~, ~, filenamearray_all(:,1,3)] = xlsread(framepath_filename_MLC,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,3)] = xlsread(framepath_filename_MLC,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,3)] = xlsread(framepath_filename_MLC,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,3)] = xlsread(framepath_filename_MLC,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,3)] = xlsread(framepath_filename_MLC,'sheet5','A2:A1001');
    
    % ��ȡRhoA��Ӧ���ļ���
    [~, ~, filenamearray_all(:,1,4)] = xlsread(framepath_filename_RhoA,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,4)] = xlsread(framepath_filename_RhoA,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,4)] = xlsread(framepath_filename_RhoA,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,4)] = xlsread(framepath_filename_RhoA,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,4)] = xlsread(framepath_filename_RhoA,'sheet5','A2:A1001');
    
    % ��ȡRBD��Ӧ���ļ���
    [~, ~, filenamearray_all(:,1,5)] = xlsread(framepath_filename_RBD,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,5)] = xlsread(framepath_filename_RBD,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,5)] = xlsread(framepath_filename_RBD,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,5)] = xlsread(framepath_filename_RBD,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,5)] = xlsread(framepath_filename_RBD,'sheet5','A2:A1001');
    
    % ��ȡVD1��Ӧ���ļ���
    [~, ~, filenamearray_all(:,1,6)] = xlsread(framepath_filename_VD1,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,6)] = xlsread(framepath_filename_VD1,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,6)] = xlsread(framepath_filename_VD1,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,6)] = xlsread(framepath_filename_VD1,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,6)] = xlsread(framepath_filename_VD1,'sheet5','A2:A1001');
    
    % ��ȡtalinA��Ӧ���ļ���
    [~, ~, filenamearray_all(:,1,7)] = xlsread(framepath_filename_talinA,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,7)] = xlsread(framepath_filename_talinA,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,7)] = xlsread(framepath_filename_talinA,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,7)] = xlsread(framepath_filename_talinA,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,7)] = xlsread(framepath_filename_talinA,'sheet5','A2:A1001');
    
    % ��ȡCdc42WT��Ӧ���ļ���
    [~, ~, filenamearray_all(:,1,8)] = xlsread(framepath_filename_Cdc42WT,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,8)] = xlsread(framepath_filename_Cdc42WT,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,8)] = xlsread(framepath_filename_Cdc42WT,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,8)] = xlsread(framepath_filename_Cdc42WT,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,8)] = xlsread(framepath_filename_Cdc42WT,'sheet5','A2:A1001');
    
    % ��ȡCdc42D118A��Ӧ���ļ���
    [~, ~, filenamearray_all(:,1,9)] = xlsread(framepath_filename_Cdc42D118A,'sheet1','A2:A1001');
    [~, ~, filenamearray_all(:,2,9)] = xlsread(framepath_filename_Cdc42D118A,'sheet2','A2:A1001');
    [~, ~, filenamearray_all(:,3,9)] = xlsread(framepath_filename_Cdc42D118A,'sheet3','A2:A1001');
    [~, ~, filenamearray_all(:,4,9)] = xlsread(framepath_filename_Cdc42D118A,'sheet4','A2:A1001');
    [~, ~, filenamearray_all(:,5,9)] = xlsread(framepath_filename_Cdc42D118A,'sheet5','A2:A1001');
end


%% step 2

% mix_sort����ֵ
mix_MultiGaussian_sort=mix_MultiGaussian_all;  % ��������
mix_EquivalentTriGaussian_sort=mix_EquivalentTriGaussian_all;
mix_TriGaussian_sort=mix_TriGaussian_all;

% ʹ��ð�����򷨽�������
% ���ȸ�������С���򣬲������ֵ����180��ϸ�������޳�
cell_number=zeros(5,9);
for i=1:5
    for j=1:9
        for k=1:1000  % �ҳ��þ����е�ʵ����������
            if mix_TriGaussian_sort(k,1,i,j)==0
                break
            end
        end
        cell_number(i,j)=k-1; % ��ʵ��������(ϸ������)��ֵ��cell_number
        
        % ð������
        for ii=2:cell_number(i,j)
            for jj=cell_number(i,j):-1:ii
                if mix_TriGaussian_sort(jj,1,i,j)<mix_TriGaussian_sort(jj-1,1,i,j)  % ��ǰһ��ϸ���������ں�һ��ϸ�����������˳��
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
        
        % ɾ�����ֵ����180��ϸ������
        for k=1:cell_number(i,j)
            if mix_TriGaussian_sort(k,1,i,j)>=180
                mix_TriGaussian_sort(k,:,i,j)=0;  % �޳�������8��ϸ������
                mix_MultiGaussian_sort(k,:,i,j)=0;
                mix_EquivalentTriGaussian_sort(k,:,i,j)=0;
                filenamearray_all{k,i,j}='0';
                cell_number(i,j)=cell_number(i,j)-1;
            end
        end
    end
end

% ɾ���������������������������������������������ƫ��ƽ����ϸ�����޳�����ǰ10%�ͺ�10%��ϸ�����ݵĲ���
memory_cell_wucha=zeros(200,5,9,6);  % �����ɾ��ϸ�����ݶ�Ӧ�����ֵ
k_cellname=ones(5,9,6);  % memory_cellname�����鳤��

% ����basal domain���˹��ƽ��ǿ�Ⱥ;ۼ��̶ȣ����浽mix_sort_basalave�����У����ں������������ɾ������
mix_sort_basalave=zeros(1000,13,5,9);
for i=1:5
    for j=1:9
        for k=1:cell_number(i,j)
            for ii=4:length(mix_MultiGaussian_sort(1,:,1,1))
                if mix_MultiGaussian_sort(k,ii,i,j)==0  % ����iiλΪ0��˵��basal domain��(ii-10)/3����˹
                    break
                end
            end
            if mix_MultiGaussian_sort(k,end-1,i,j)~=0  % �����һλ��Ϊ0��˵��basal domain��5(max)����˹
                ii=25;
            end
            basal_Gaussian_number=(ii-10)/3;

            % ����basal domain���˹��ƽ��ǿ�Ⱥ;ۼ��̶�
            basal_overactivity_average=sum(mix_MultiGaussian_all(k,4:basal_Gaussian_number+3,i,j))/basal_Gaussian_number;
            basal_signalwidth_average=sum(mix_MultiGaussian_all(k,basal_Gaussian_number+4:2*basal_Gaussian_number+3,i,j))/basal_Gaussian_number;
            basal_expectation_average=sum(mix_MultiGaussian_all(k,2*basal_Gaussian_number+4:3*basal_Gaussian_number+3,i,j))/basal_Gaussian_number;

            % ��ֵ��mix_sort_basalave
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

% ������������������ʹ�ö��˹��ϵĽ��Ϊ����
number_sort_position=[7 8 10 11];  % ����������Ӧ��mix_all�����λ��
for k=1:4
    number_sort=number_sort_position(k);
    for i=1:5
        for j=1:9
            % ð������
            for ii=2:cell_number(i,j)
                for jj=cell_number(i,j):-1:ii
                    if mix_sort_basalave(jj,number_sort,i,j)<mix_sort_basalave(jj-1,number_sort,i,j)  % ��ǰһ��ϸ���������������ں�һ��ϸ�������������������˳��
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

            % �ҳ�ǰ5%��ϸ���������ֵ���
            if number_sort==7 || number_sort==10
                for jj=1:floor(cell_number(i,j)/20)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_TriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
                % �ҳ���5%��ϸ���������ֵ���
                for jj=cell_number(i,j)-floor(cell_number(i,j)/20)+1:cell_number(i,j)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_TriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
            else
                for jj=1:floor(cell_number(i,j)/10)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_TriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
                % �ҳ���5%��ϸ���������ֵ���
                for jj=cell_number(i,j)-floor(cell_number(i,j)/10)+1:cell_number(i,j)
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_TriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
            end
        end
    end
end

% �������������������ʹ��ֱ������˹�͵�Ч����˹��ϵıȽ�����(�����ַ�������Ͻ�����쳬��20%����ʹ��ֱ������˹��ϲ�������֮ʹ�õ�Ч����˹��ϲ���)
number_sort_position=[4 5];  % ����������Ӧ��mix_all�����λ��
for k=5:6
    number_sort=number_sort_position(k-4);
    for i=1:5
        for j=1:9
            % ð������
            for ii=2:cell_number(i,j)
                for jj=cell_number(i,j):-1:ii
                    if mix_TriGaussian_sort(jj,number_sort,i,j)<mix_TriGaussian_sort(jj-1,number_sort,i,j)  % ��ǰһ��ϸ���������������ں�һ��ϸ�������������������˳��
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

            % �ҳ�ǰ5%��ϸ���������ֵ���
            for jj=1:floor(cell_number(i,j)/20)
                if abs(mix_TriGaussian_sort(jj,number_sort,i,j)-mix_EquivalentTriGaussian_sort(jj,number_sort,i,j))/mix_TriGaussian_sort(jj,number_sort,i,j)>0.2
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_TriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                else
                    memory_cell_wucha(k_cellname(i,j,k),i,j,k)=mix_EquivalentTriGaussian_sort(jj,1,i,j);
                    k_cellname(i,j,k)=k_cellname(i,j,k)+1;
                end
            end
            % �ҳ���5%��ϸ���������ֵ���
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


% ȡ��������õ���ϸ���Ĳ���������mix_all������ɾȥ
mix_TriGaussian_sort2=mix_TriGaussian_sort;  % ��mix_TriGaussian_sort�������ݸ�ֵ��mix_TriGaussian_sort2
mix_EquivalentTriGaussian_sort2=mix_EquivalentTriGaussian_sort;  % ��mix_EquivalentTriGaussian_sort�������ݸ�ֵ��mix_EquivalentTriGaussian_sort2
mix_MultiGaussian_sort2=mix_MultiGaussian_sort;  % ��mix_MultiGaussian_sort�������ݸ�ֵ��mix_MultiGaussian_sort2

mix_TriGaussian_sort=zeros(1000,13,5,9);  % ������飬���ں������ɾ��������ϸ�����ϸ���������
mix_EquivalentTriGaussian_sort=zeros(1000,13,5,9);  % ������飬���ں������ɾ��������ϸ�����ϸ���������
mix_MultiGaussian_sort=zeros(1000,25,5,9);  % ������飬���ں������ɾ��������ϸ�����ϸ���������
filenamearray_all_sort=cell(1000,5,9);
for i=1:5
    for j=1:9
        kk2=1;  % ���ڴ洢mix_allÿ��ʱ�ڣ�ÿ�ֵ��׵�ϸ����
        for k1=1:1000
            if mix_TriGaussian_sort2(k1,1,i,j)~=0  % ���ڵ�k1��ϸ��
                logic=0;  % �����жϸ�ϸ���Ƿ���Ҫɾ��
                for k=1:6
                    for k2=1:k_cellname(i,j,k)
                        if mix_TriGaussian_sort2(k1,1,i,j)==memory_cell_wucha(k2,i,j,k)
                            logic=1;  % ����ϸ����wuchaֵ,���������������,��ǰ������������Ϊǰ5%�ͺ�5%��ϸ������˵����ϸ����Ҫɾ��
                            break
                        end
                    end
                end
                if logic==0  % ��logicΪ0�����ʾ��ϸ������ɾ�������丳ֵ��mix_sort filenamearray_all����
                    mix_TriGaussian_sort(kk2,:,i,j)=mix_TriGaussian_sort2(k1,:,i,j);
                    mix_EquivalentTriGaussian_sort(kk2,:,i,j)=mix_EquivalentTriGaussian_sort2(k1,:,i,j);
                    mix_MultiGaussian_sort(kk2,:,i,j)=mix_MultiGaussian_sort2(k1,:,i,j);
                    filenamearray_all_sort(kk2,i,j)=filenamearray_all(k1,i,j);
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
mix_TriGaussian_all=mix_TriGaussian_sort;
mix_EquivalentTriGaussian_all=mix_EquivalentTriGaussian_sort;
mix_MultiGaussian_all=mix_MultiGaussian_sort;
filenamearray_all=filenamearray_all_sort;



%%% �����ʱ�ڸ���ϲ����ľ�ֵ
% ����˹������������
TriGaussian_summarize=mix_EquivalentTriGaussian_sort;



%% step 3
% ���Ƚ�ÿ�ֵ��׵�ϸ������ֵ����lumen length/cell length��С��������
for i=1:5
    for j=1:9
        % ð������
        for ii=2:cell_number(i,j)
            for jj=cell_number(i,j):-1:ii
                if TriGaussian_summarize(jj,end,i,j)<TriGaussian_summarize(jj-1,end,i,j)  % ��ǰһ��ϸ���������������ں�һ��ϸ�������������������˳��
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


pre_cell_number=zeros(9,1);  % �����׵�lumen�γ�ǰ��ϸ����
pre_cell_number_monocell=zeros(9,1);  % �����׵�lumen�γ�ǰ���е��������ļ���ϸ����
pre_cell_number_multicell=zeros(9,1);  % �����׵�lumen�γ�ǰ�������������ļ���ϸ����
pre_AL_thickness=zeros(200,9,2);  % lumen�γ�ǰϸ����anterior-laterior thickness���ڶ���ά�ȱ�ʾ��ͬ�ֵ��ף�������ά�ȱ�ʾ�����������ļ���ϸ�����е��������ļ���ϸ��
pre_AL_length=zeros(200,9,2);  % lumen�γ�ǰϸ����\anterior-laterior length
pre_PL_thickness=zeros(200,9,2);  % lumen�γ�ǰϸ����\posterior-laterior thickness
pre_PL_length=zeros(200,9,2);  % lumen�γ�ǰϸ����\posterior-laterior length
pre_average_anteriorlength=zeros(1,9);pre_average_posteriorlength=zeros(1,9);  % ����anterior-��posterior-lateral domain��cortexƽ�����ȱ���
pre_average_anteriorthickness=zeros(1,9);pre_average_posteriorthickness=zeros(1,9);  % ����anterior-��posterior-lateral domain��cortexƽ����ȱ���

% �ҳ�����lumen length/cell length<0.08��ϸ�������ļ������������anterior-��posterior-lateral domain��cortexƽ����Ⱥ�ƽ������
for i=1:9
    % ͳ��lumen�γ�ǰ��ϸ����
    for j=1:cell_number(1,i)
        if TriGaussian_summarize(j,end,1,i)<0.08
            pre_cell_number(i)=pre_cell_number(i)+1;
        else
            break
        end
    end
    
    if pre_cell_number(i)==0  % �ж��Ƿ����lumen�γ�ǰ��ϸ��
        continue
    end
    
    for j=1:pre_cell_number(i)
        filename=mat2str(cell2mat(filenamearray_all(j,1,i)));  % ��ȡ�ļ���
        filename=filename(2:length(filename)-1);  % ɾ������ĵ�����
        S=regexp(filename,'_','split');  % ����ļ����ַ���
        month_day=mat2str(cell2mat(S(1)));  % S����ĵ�һ��ֵΪ��ϸ��������(month+day)
        month_day=month_day(2:length(month_day)-1);  % ɾ������ĵ�����
        framepath_filename_1=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_xiahuaxian filename_1 filename_csv];  % ƴ��������ϸ���Ҷ�ֵ�����ļ�ԭʼ·��(����ļ���ʽ)
        framepath_filename_2=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_xiahuaxian filename_2 filename_csv];
        framepath_filename_3=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_xiahuaxian filename_3 filename_csv];
        framepath_filename_danyi=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_csv];  % ƴ��������ϸ���Ҷ�ֵ�����ļ�ԭʼ·��(�����ļ���ʽ)

        % ���ݸ�ϸ���������������������tri-Gaussian distribution
        rho0=TriGaussian_summarize(j,3,1,i);rho_inf=TriGaussian_summarize(j,4,1,i);w=TriGaussian_summarize(j,5,1,i);E=TriGaussian_summarize(j,6,1,i);  % ���������������������ֵ
        rho_inf_anterior=TriGaussian_summarize(j,7,1,i);w_anterior=TriGaussian_summarize(j,8,1,i);E_anterior=TriGaussian_summarize(j,9,1,i);
        rho_inf_posterior=TriGaussian_summarize(j,10,1,i);w_posterior=TriGaussian_summarize(j,11,1,i);E_posterior=TriGaussian_summarize(j,12,1,i);  % ��������������������ֵ
        f=@(x) rho0+rho_inf.*exp(-1/2.*((x-E)./w).^2)+rho_inf_anterior.*exp(-1/2.*((x-E_anterior)./w_anterior).^2)+rho_inf_posterior.*exp(-1/2.*((x-E_posterior)./w_posterior).^2);
        
        % ����anterior-��posterior-lateral domain��cortexƽ����Ⱥ�ƽ������
        if exist(framepath_filename_1,'file')  % ͳ��lumen�γ�ǰ��ϸ���������������ļ���ϸ��������������anterior-��posterior-lateral domain��cortexƽ����Ⱥ�ƽ������
            pre_cell_number_multicell(i)=pre_cell_number_multicell(i)+1;  % ���������ļ���ϸ����ͳ��
            
            data_1=xlsread(framepath_filename_1); % ��������
            data_2=xlsread(framepath_filename_2);
            data_3=xlsread(framepath_filename_3);
            pre_AL_length(pre_cell_number_multicell(i),i,1)=length(data_1)/(length(data_1)+length(data_2)+length(data_3));  % anterior-laterior lengthͳ��
            pre_PL_length(pre_cell_number_multicell(i),i,1)=length(data_3)/(length(data_1)+length(data_2)+length(data_3));  % posterior-laterior lengthͳ��
            
            a_anterior=integral(f,-1/2,-1/2+pre_AL_length(pre_cell_number_multicell(i),i,1));  % anterior-laterior cortex thickness����
            a_anterior=a_anterior/pre_AL_length(pre_cell_number_multicell(i),i,1);  % ���ֽ������x�᳤�ȣ��õ�ƽ��cortex thickness
            pre_AL_thickness(pre_cell_number_multicell(i),i,1)=a_anterior;  % ƽ��cortex thicknessֵ��ֵ��pre_average_AT����
            a_posterior=integral(f,1/2-pre_PL_length(pre_cell_number_multicell(i),i,1),1/2);  % posterior-laterior cortex thickness����
            a_posterior=a_posterior/pre_PL_length(pre_cell_number_multicell(i),i,1);  % ���ֽ������x�᳤�ȣ��õ�ƽ��cortex thickness
            pre_PL_thickness(pre_cell_number_multicell(i),i,1)=a_posterior;  % ƽ��cortex thicknessֵ��ֵ��pre_average_PT����
        elseif exist(framepath_filename_danyi,'file')  % ͳ��lumen�γ�ǰ��ϸ�����е��������ļ���ϸ��������������anterior-��posterior-lateral domain��cortexƽ�����
            pre_cell_number_monocell(i)=pre_cell_number_monocell(i)+1;  % ���������ļ���ϸ����ͳ��
            
            % �����ļ���ȷ����ϸ���ڸ������µ�λ��
            number_changdubi=0;
            for ii=1:30  % ���ѭ�����ڲ���Ŀ���ļ��Ƿ����(i��ʾ�ڼ���photo��j��ʾphoto�еĵڼ���ϸ����k��ʾӫ��ͨ��)
                filename1=num2str(ii);
                for jj=1:10
                    filename2=num2str(jj);
                    filename3='1';
                    filename_temp=[month_day filename_xiahuaxian filename1 filename2 filename3];  % �� 0928_111
                    
                    % ƴ���������ļ�·����������һ��exist�����Ĳ���
                    framepath_filename_1_temp=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_temp filename_xiahuaxian filename_1 filename_csv];  % ƴ��������ϸ���Ҷ�ֵ�����ļ�ԭʼ·��(����ļ���ʽ)
                    framepath_filename_danyi_temp=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_temp filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111.csv
                    
                    if exist(framepath_filename_danyi_temp,'file')   % һ��ϸ���������е����ļ�
                        number_changdubi=number_changdubi+1;
                    elseif exist(framepath_filename_1_temp,'file')
                        number_changdubi=number_changdubi+1;
                    elseif strcmp(filename_temp,filename)  % ��������Ƭ�£������ڸ�ϸ���������ļ���˵��������Ƭ�е�ϸ���Ѿ������ɣ�����j(ϸ��)ѭ��
                        number_changdubi=number_changdubi+1;
                        break
                    end
                end
            end
            
            % ���Ҹ�ϸ����Ӧ��lumen length��cell length (��data_changdubi��Ѱ��)
            framepath2=[framepath num2str(k) filename_fanxiegang month_day filename_fanxiegang];  % ���ļ���·������ 'C:\Users\Desktop\1\0928'
            framepath_filename_changdubi=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_changdubi filename_xlsx];
            framepath_filename_Radius=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_Radius filename_xlsx];
            framepath_filename_ContactAngle=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_ContactAngle filename_xlsx];

            Lumen_Radius=xlsread(framepath_filename_Radius,'sheet1','A1:A100');  % ����lumen�뾶
            Lumen_Radius=sqrt(Lumen_Radius./pi);
            Lumen_ContactAngle=xlsread(framepath_filename_ContactAngle,'sheet1','A1:A100');  % ����lumen�Ӵ���
            data_lumenlength=xlsread(framepath_filename_changdubi,'sheet1','A1:A100');  % ������ļ��µ�lumen����
            data_celllength=xlsread(framepath_filename_changdubi,'sheet1','B1:B100');  % ������ļ��µ�cell����

            data=xlsread(framepath_filename_danyi); % ��������
            anterior_lateral_length=(data_celllength(number_changdubi)-Lumen_Radius(number_changdubi)*sin(Lumen_ContactAngle(number_changdubi)*pi./180))/2;  % anterior-lateral domain����
            posterior_lateral_length=(data_celllength(number_changdubi)-Lumen_Radius(number_changdubi)*sin(Lumen_ContactAngle(number_changdubi)*pi./180))/2;  % posterior-lateral domain����
            basal_length=data(end,1)-anterior_lateral_length-posterior_lateral_length;  % basal domain����


            a_anterior=integral(f,-1/2,-1/2+anterior_lateral_length/data(end,1));  % anterior-lateral cortex thickness����
            a_anterior=a_anterior/(anterior_lateral_length/data(end,1));  % ���ֽ������x�᳤�ȣ��õ�ƽ��cortex thickness
            pre_AL_length(pre_cell_number_monocell(i),i,2)=anterior_lateral_length/data(end,1);  % anterior length��ֵ��pre_AL_length����
            pre_AL_thickness(pre_cell_number_monocell(i),i)=a_anterior;  % ƽ��cortex thicknessֵ��ֵ��pre_AL_thickness����
            a_posterior=integral(f,1/2-posterior_lateral_length/data(end,1),1/2);  % posterior-lateral cortex thickness����
            a_posterior=a_posterior/(posterior_lateral_length/data(end,1));  % ���ֽ������x�᳤�ȣ��õ�ƽ��cortex thickness
            pre_PL_length(pre_cell_number_monocell(i),i,2)=posterior_lateral_length/data(end,1);  % posterior length��ֵ��pre_LL_length����
            pre_PL_thickness(pre_cell_number_monocell(i),i)=a_posterior;  % ƽ��cortex thicknessֵ��ֵ��pre_PL_thickness����
        end
    end

    % ����lumen�γ�ǰ��anterior-��posterior-lateral domain��cortexƽ�����
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
%  ���Ƚ�mix_all������ÿ�ֵ��׵�5��ʱ�ڵ����ݺϲ���������TriGaussian_summarize2������
TriGaussian_summarize2=zeros(2000,13,9);  % ����9�ֵ��׵����ݣ�5��ʱ�ںϲ�
filenamearray_all2=cell(2000,9);
cell_number2=zeros(9,1);  % ����9�ֵ��׵�ϸ������5��ʱ�ںϲ�
for i=1:9
    k2=0;  % k2ָʾmix_all2��λ�ã�ѭ��������Ϊ5��ʱ��ϸ�������ܺ�
    for j=1:5
        for k=1:cell_number(j,i)  % kΪmix_all��λ��
            k2=k2+1;
            TriGaussian_summarize2(k2,:,i)=TriGaussian_summarize(k,:,j,i);
            filenamearray_all2(k2,i)=filenamearray_all(k,j,i);
        end
    end
    cell_number2(i)=k2;
end

cell_number_monocell=zeros(9,1);  % �����׵�lumen�γɺ���е��������ļ���ϸ����
filenamearray_monocell=cell(1000,9);  % �����е��������ļ���ϸ�������anterior-/posterior-lateral thickness��Ӧ���ļ�����·��
anterior_thickness_monocell=zeros(1000,9,3);posterior_thickness_monocell=zeros(1000,9,3);  % ����anterior-��posterior-lateral domain��cortex��ȱ�����������ά�ȵ�����ֱ��ź�ȡ���Ӧ�ĳ��Ⱥ�lumen length/cell length
cell_number_multicell=zeros(9,1);  % �����׵�lumen�γɺ�������������ļ���ϸ����
filenamearray_multicell_anterior=cell(1000,9);  % �����ж�������ļ���ϸ�������anterior-lateral thickness��Ӧ���ļ�����·��
filenamearray_multicell_posterior=cell(1000,9);  % �����ж�������ļ���ϸ�������posterior-lateral thickness��Ӧ���ļ�����·��
anterior_length_summarize=zeros(2000,2);posterior_length_summarize=zeros(2000,2);  % ����anterior-��posterior-lateral domain��cortex���ȱ������ڶ���ά�ȵ�����ֱ��ų��ȺͶ�Ӧ��lumen length/cell lengthֵ
anterior_thickness_multicell=zeros(1000,9,3);posterior_thickness_multicell=zeros(1000,9,3);  % ����anterior-��posterior-lateral domain��cortex��ȱ�����������ά�ȵ�����ֱ��ź�ȡ���Ӧ�ĳ��Ⱥ�lumen length/cell length
k_length=0;  % ָʾanterior_length_summarize��ֵλ��
for i=1:9
    for j=pre_cell_number(i)+1:cell_number2(i)
        filename=mat2str(cell2mat(filenamearray_all2(j,i)));  % ��ȡ�ļ���
        filename=filename(2:length(filename)-1);  % ɾ������ĵ�����
        S=regexp(filename,'_','split');  % ����ļ����ַ���
        month_day=mat2str(cell2mat(S(1)));  % S����ĵ�һ��ֵΪ��ϸ��������(month+day)
        month_day=month_day(2:length(month_day)-1);  % ɾ������ĵ�����
        framepath_filename_1=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_xiahuaxian filename_1 filename_csv];  % ƴ��������ϸ���Ҷ�ֵ�����ļ�ԭʼ·��(����ļ���ʽ)
        framepath_filename_2=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_xiahuaxian filename_2 filename_csv];
        framepath_filename_3=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_xiahuaxian filename_3 filename_csv];
        framepath_filename_danyi=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename filename_csv];  % ƴ��������ϸ���Ҷ�ֵ�����ļ�ԭʼ·��(�����ļ���ʽ)

        % ���ݸ�ϸ���������������������tri-Gaussian distribution
        rho0=TriGaussian_summarize2(j,3,i);rho_inf=TriGaussian_summarize2(j,4,i);w=TriGaussian_summarize2(j,5,i);E=TriGaussian_summarize2(j,6,i);  % ���������������������ֵ
        rho_inf_anterior=TriGaussian_summarize2(j,7,i);w_anterior=TriGaussian_summarize2(j,8,i);E_anterior=TriGaussian_summarize2(j,9,i);
        rho_inf_posterior=TriGaussian_summarize2(j,10,i);w_posterior=TriGaussian_summarize2(j,11,i);E_posterior=TriGaussian_summarize2(j,12,i);  % ��������������������ֵ
        f=@(x) rho0+rho_inf.*exp(-1/2.*((x-E)./w).^2)+rho_inf_anterior.*exp(-1/2.*((x-E_anterior)./w_anterior).^2)+rho_inf_posterior.*exp(-1/2.*((x-E_posterior)./w_posterior).^2);
        
        % ����anterior-��posterior-lateral domain��cortexƽ����Ⱥ�ƽ������
        if exist(framepath_filename_1,'file')
            %% Step 4.1�����������������ļ���ϸ����anterior-��posterior-lateral domain��cortex��Ⱥͳ��ȣ������anterior/posterior-
            %%           lateral domain length/lateral-basal-lateral whole length��lumen length/cell length�ĺ�����ϵ
            cell_number_multicell(i)=cell_number_multicell(i)+1;  % ���������ļ���ϸ����ͳ��
            
            data_1=xlsread(framepath_filename_1); % ��������
            data_2=xlsread(framepath_filename_2);
            data_3=xlsread(framepath_filename_3);
            
            % anterior-/posterior- cortex thickness�Ͷ�Ӧ��lateral domain length��ֵ(��ͬ�����в�ͬ�ķֲ���ʽ���˷ֿ�����)
%             a_anterior=integral(f,-1/2,-1/2+length(data_1)/(length(data_1)+length(data_2)+length(data_3)));  % anterior-lateral cortex thickness����
%             a_anterior=a_anterior/length(data_1)*(length(data_1)+length(data_2)+length(data_3));  % ���ֽ������x�᳤�ȣ��õ�ƽ��cortex thickness
            a_anterior=integral(f,-1/2,-1/2+0.05);  % anterior-lateral cortex thickness����
            a_anterior=a_anterior/0.05;  % ���ֽ������x�᳤�ȣ��õ�ƽ��cortex thickness
            anterior_thickness_multicell(cell_number_multicell(i),i,1)=a_anterior;  % cortex thicknessֵ��ֵ��anterior_thickness_multifile����
            anterior_thickness_multicell(cell_number_multicell(i),i,2)=length(data_1)/(length(data_1)+length(data_2)+length(data_3));  % ��ϸ��cortex thickness��Ӧ��anterior-lateral length
            anterior_thickness_multicell(cell_number_multicell(i),i,3)=TriGaussian_summarize2(j,end,i);  % ��ϸ����lumen length/cell length
            filenamearray_multicell_anterior{cell_number_multicell(i),i}=framepath_filename_1;
            
%             a_posterior=integral(f,1/2-length(data_3)/(length(data_1)+length(data_2)+length(data_3)),1/2);  % posterior-lateral cortex thickness����
%             a_posterior=a_posterior/length(data_3)*(length(data_1)+length(data_2)+length(data_3));  % ���ֽ������x�᳤�ȣ��õ�ƽ��cortex thickness
            a_posterior=integral(f,1/2-0.05,1/2);  % posterior-lateral cortex thickness����
            a_posterior=a_posterior/0.05;  % ���ֽ������x�᳤�ȣ��õ�ƽ��cortex thickness
            posterior_thickness_multicell(cell_number_multicell(i),i,1)=a_posterior;  % cortex thicknessֵ��ֵ��posterior_thickness_multifile����
            posterior_thickness_multicell(cell_number_multicell(i),i,2)=length(data_3)/(length(data_1)+length(data_2)+length(data_3));  % ��ϸ��cortex thickness��Ӧ��posterior-lateral length
            posterior_thickness_multicell(cell_number_multicell(i),i,3)=TriGaussian_summarize2(j,end,i);  % ��ϸ����lumen length/cell length
            filenamearray_multicell_posterior{cell_number_multicell(i),i}=framepath_filename_3;
            
            % anterior-/posterior-lateral length�Ͷ�Ӧ��lumen length/cell length��ֵ(��ͬ���ײ�Ӱ��anterior-/posterior-lateral length���˺ϲ�����������������)
            logic=0;
            for k=1:k_length  % �жϸ�ϸ���Ƿ���anterior_length�������Ѿ�����
                if anterior_length_summarize(k,2)==TriGaussian_summarize2(j,end,i)
                    logic=1;
                    break
                end
            end
            if logic==0  % ��ϸ����anterior_length�����в����ڣ���ֵ�����ϸ����anterior-/posterior-lateral length
                k_length=k_length+1;
                anterior_length_summarize(k_length,1)=length(data_1)/(length(data_1)+length(data_2)+length(data_3));  % anterior-lateral lengthͳ��
                posterior_length_summarize(k_length,1)=length(data_3)/(length(data_1)+length(data_2)+length(data_3));  % posterior-lateral lengthͳ��
                anterior_length_summarize(k_length,2)=TriGaussian_summarize2(j,end,i);  % anterior-lateral length��Ӧ��lumen length/cell length
                posterior_length_summarize(k_length,2)=TriGaussian_summarize2(j,end,i);
            end
        elseif exist(framepath_filename_danyi,'file')
            %% Step 4.2��������һ�������ļ���ϸ����anterior-��posterior-lateral domain��cortex��ȣ�ϸ�������ֵĲ������
            %%           lateral-basal-lateral whole length��lumen length/cell length
            cell_number_monocell(i)=cell_number_monocell(i)+1;  % ���������ļ���ϸ����ͳ��
            
            % �����ļ���ȷ����ϸ���ڸ������µ�λ��
            number_changdubi=0;
            for ii=1:30  % ���ѭ�����ڲ���Ŀ���ļ��Ƿ����(i��ʾ�ڼ���photo��j��ʾphoto�еĵڼ���ϸ����k��ʾӫ��ͨ��)
                filename1=num2str(ii);
                for jj=1:10
                    filename2=num2str(jj);
                    filename3='1';
                    filename_temp=[month_day filename_xiahuaxian filename1 filename2 filename3];  % �� 0928_111
                    
                    % ƴ���������ļ�·����������һ��exist�����Ĳ���
                    framepath_filename_1_temp=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_temp filename_xiahuaxian filename_1 filename_csv];  % ƴ��������ϸ���Ҷ�ֵ�����ļ�ԭʼ·��(����ļ���ʽ)
                    framepath_filename_danyi_temp=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_temp filename_csv];  % �� 'C:\Users\Desktop\1\0928\0928_111.csv
                    
                    if exist(framepath_filename_danyi_temp,'file')   % һ��ϸ���������е����ļ�
                        number_changdubi=number_changdubi+1;
                    elseif exist(framepath_filename_1_temp,'file')
                        number_changdubi=number_changdubi+1;
                    elseif strcmp(filename_temp,filename)  % ��������Ƭ�£������ڸ�ϸ���������ļ���˵��������Ƭ�е�ϸ���Ѿ������ɣ�����j(ϸ��)ѭ��
                        number_changdubi=number_changdubi+1;
                        break
                    end
                end
            end
            
            % ���Ҹ�ϸ����Ӧ��lumen length��cell length (��data_changdubi��Ѱ��)
            framepath_filename_changdubi=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_changdubi filename_xlsx];
            framepath_filename_Radius=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_Radius filename_xlsx];
            framepath_filename_ContactAngle=[framepath filename_1 filename_fanxiegang month_day filename_fanxiegang filename_ContactAngle filename_xlsx];

            Lumen_Radius=xlsread(framepath_filename_Radius,'sheet1','A1:A100');  % ����lumen�뾶
            Lumen_Radius=sqrt(Lumen_Radius./pi);
            Lumen_ContactAngle=xlsread(framepath_filename_ContactAngle,'sheet1','A1:A100');  % ����lumen�Ӵ���
            data_lumenlength=xlsread(framepath_filename_changdubi,'sheet1','A1:A100');  % ������ļ��µ�lumen����
            data_celllength=xlsread(framepath_filename_changdubi,'sheet1','B1:B100');  % ������ļ��µ�cell����

            data=xlsread(framepath_filename_danyi); % ��������
            anterior_lateral_length=(data_celllength(number_changdubi)-Lumen_Radius(number_changdubi)*sin(Lumen_ContactAngle(number_changdubi)*pi./180))/2;  % anterior-lateral domain����
            posterior_lateral_length=(data_celllength(number_changdubi)-Lumen_Radius(number_changdubi)*sin(Lumen_ContactAngle(number_changdubi)*pi./180))/2;  % posterior-lateral domain����
            basal_length=data(end,1)-anterior_lateral_length-posterior_lateral_length;  % basal domain����

            % anterior-/posterior- cortex thickness�Ͷ�Ӧ��lateral domain length��ֵ(��ͬ�����в�ͬ�ķֲ���ʽ���˷ֿ�����)
%             a_anterior=integral(f,-1/2,-1/2+anterior_lateral_length/data(end,1));  % anterior-lateral cortex thickness����
%             a_anterior=a_anterior/(anterior_lateral_length/data(end,1));  % ���ֽ������x�᳤�ȣ��õ�ƽ��cortex thickness
            a_anterior=integral(f,-1/2,-1/2+0.05);  % anterior-lateral cortex thickness����
            a_anterior=a_anterior/0.05;  % ���ֽ������x�᳤�ȣ��õ�ƽ��cortex thickness
            anterior_thickness_monocell(cell_number_monocell(i),i,1)=a_anterior;  % cortex thicknessֵ��ֵ��anterior_thickness_multifile����
            anterior_thickness_monocell(cell_number_monocell(i),i,2)=anterior_lateral_length/data(end,1);  % ��ϸ��cortex thickness��Ӧ��anterior-lateral length
            anterior_thickness_monocell(cell_number_monocell(i),i,3)=TriGaussian_summarize2(j,end,i);  % ��ϸ����lumen length/cell length
            
%             a_posterior=integral(f,1/2-posterior_lateral_length/data(end,1),1/2);  % posterior-lateral cortex thickness����
%             a_posterior=a_posterior/(posterior_lateral_length/data(end,1));  % ���ֽ������x�᳤�ȣ��õ�ƽ��cortex thickness
            a_posterior=integral(f,1/2-0.05,1/2);  % posterior-lateral cortex thickness����
            a_posterior=a_posterior/0.05;  % ���ֽ������x�᳤�ȣ��õ�ƽ��cortex thickness
            posterior_thickness_monocell(cell_number_monocell(i),i,1)=a_posterior;  % cortex thicknessֵ��ֵ��posterior_thickness_multifile����
            posterior_thickness_monocell(cell_number_monocell(i),i,2)=posterior_lateral_length/data(end,1);  % ��ϸ��cortex thickness��Ӧ��posterior-lateral length
            posterior_thickness_monocell(cell_number_monocell(i),i,3)=TriGaussian_summarize2(j,end,i);  % ��ϸ����lumen length/cell length
            filenamearray_monocell{cell_number_monocell(i),i}=framepath_filename_danyi;

            % anterior-/posterior-lateral length�Ͷ�Ӧ��lumen length/cell length��ֵ(��ͬ���ײ�Ӱ��anterior-/posterior-lateral length���˺ϲ�����������������)
            logic=0;
            for k=1:k_length  % �жϸ�ϸ���Ƿ���anterior_length�������Ѿ�����
                if anterior_length_summarize(k,2)==TriGaussian_summarize2(j,end,i)
                    logic=1;
                    break
                end
            end
            if logic==0  % ��ϸ����anterior_length�����в����ڣ���ֵ�����ϸ����anterior-/posterior-lateral length
                k_length=k_length+1;
                anterior_length_summarize(k_length,1)=anterior_lateral_length/data(end,1);  % anterior-lateral lengthͳ��
                posterior_length_summarize(k_length,1)=posterior_lateral_length/data(end,1);  % posterior-lateral lengthͳ��
                anterior_length_summarize(k_length,2)=TriGaussian_summarize2(j,end,i);  % anterior-lateral length��Ӧ��lumen length/cell length
                posterior_length_summarize(k_length,2)=TriGaussian_summarize2(j,end,i);
            end
        end
    end
end

% ��anterior/posterior_length������lumen length/cell length��������(ð������)
for i=k_length:-1:2
    % anterior_length��������
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
    
    % posterior_length��������
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

% ���anterior/posterior-lateral domain length/lateral-basal-lateral whole length��t�ĺ�����ϵ
% �ҳ�anterior/posterior_length_summarize������lumen length/cell length��Ӧ��ʱ��t
% ����ʵ�ʹ۲����transverse diameter(TD)��longitudinal radius(LR)��ʱ��ı仯����
TD=[3.248130325	3.697268921	4.042205695	4.875820261	5.623190217	5.105764223	5.795637773	5.996833531	6.456763119	6.571755932	6.485511322	6.370518509	6.399266712	6.485511322	6.370518509	6.600504135	6.284315564	6.226819157	6.226819157	6.830448097	7.37658063	7.577818053	8.411390955	8.267691603	8.900068746	9.474949482	10.02108202	10.4235152	11.14213695	11.77451409	12.57938046	12.86682082	13.38420515	13.64293898	14.21781972	14.36151907	14.5915047	15.10888903	15.14842822	15.05139262	16];
LR=[0.8114597 0.88031998 1.012697123 0.632377143 0.966522926 1.167760348 1.523446451 1.943878508 2.08757786	2.288815282	2.202570673	2.949940629	3.381080349	3.668520717	4.214694915	4.818323855	4.703331042	5.174051622	5.192008833	5.651938421	5.795637773	5.968126992	6.054329938	6.140574547	6.169322751	5.910630586	5.853134179	6.111826344	6.284315564	5.996833531	6.083078141	5.996833531	5.824385976	5.863883507	5.853134179	5.853134179	5.795637773	5.939378789	5.738141366	6.111826344	5.996833531];
TD2=[TD 16.1 16.1 16.1 16.1 16.1 16.1 16.1 16.1 16.1 16.1];
t=linspace(0,120,length(TD))+6;  % �۲����ݵ�ʱ��t
t2=[t t(end) t(end) t(end) t(end) t(end) t(end) t(end) t(end) t(end) t(end)];
tlength=1000;  % ���ʱ����ľ���
t_fit=linspace(0,126,tlength);  % ���ʱ����
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

% ����TD��LR��R��theta�Ĺ�ϵ��ͨ���������������R��theta��ʱ��ı仯
R=zeros(1,length(t_fit));
theta=zeros(1,length(t_fit));
for i=1:length(t_fit)
    wucha_best=1000;  % ��С���ֵ
    for theta_search=0:pi/1000:pi
        R_search=TD_fit(i)/2/sin(theta_search);
        if abs(R_search+R_search*cos(theta_search)-LR_fit(i))<wucha_best
            wucha_best=abs(R_search+R_search*cos(theta_search)-LR_fit(i));
            R(i)=R_search;
            theta(i)=theta_search*180/pi;
        end
    end
end


% ��lumen length/cell length(lateral_tension_all(:,:,2))ת��Ϊ��Ӧ��ʱ��t��������t_rectify������
t_rectify=zeros(1,tlength);  % �������׹۲����ݺ��ʱ��t
t_fit_rectify=linspace(0,126,tlength);  % �������׹۲����ݺ��������ߺ�����t
TD_rectify=zeros(1,tlength);  % �������TD
TD_fit_rectify=zeros(1,tlength);  % �������TD��Ϻ���
TD_jump=6.2;  % lumen length/cell length����Ծֵ

i_rectify_position=0;  % ��ֵ��rectify�����е�λ��
for i=1:tlength  % ȥ��TD_jump��Χ��������
    if TD_fit(i)>=TD_jump+1 || TD_fit(i)<=TD_jump-1
        i_rectify_position=i_rectify_position+1;
        t_rectify(i_rectify_position)=t_fit(i);
        TD_rectify(i_rectify_position)=TD_fit(i);
    end
end
% ��TD_rectify�������
% N_rectify=6;
% p=polyfit(t_rectify(1:i_rectify_position),TD_rectify(1:i_rectify_position),N_rectify);
% for k=1:N_rectify+1
%     TD_fit_rectify=TD_fit_rectify+p(k).*t_fit_rectify.^(N_rectify-k+1);
% end
% ��TD_rectify���в�ֵ
TD_fit_rectify=interp1(t_rectify(1:i_rectify_position),TD_rectify(1:i_rectify_position),t_fit_rectify,'pchip');  % ����������ֵ

tlength=1000;
t_observation_AL=zeros(1000,1);  % ʵ��۲��lumen length/cell length��Ӧ��ʱ��t
t_observation_PL=zeros(1000,1);  % ʵ��۲��lumen length/cell length��Ӧ��ʱ��t
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

N=3;  % ��Ͻ���
lateral_lumen_length_fit=zeros(N+1,2);  % ������϶���ʽϵ����������һ��Ϊanterior�ĺ�����ϵ���ڶ���Ϊposterior�ĺ�����ϵ
p_fit=polyfit(t_observation_AL(1:k_length),anterior_length_summarize(1:k_length,1),N);  % anterior-lateral domain length�Ķ���ʽ���
for k=1:N+1
    lateral_lumen_length_fit(k,1)=p_fit(k);
end
p_fit=polyfit(t_observation_PL(1:k_length),posterior_length_summarize(1:k_length,1),N);  % posterior-lateral domain length�Ķ���ʽ���
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
%% Step 5.1�����ȼ������඼��ϸ���Ҷ�ֵ���ݵ�lumen���ܵ������Լ���ϸ����������������������С����������С����ǰһ��ϸ
%%           ����posterior-lateral domain��cortex��ȡ����Ⱥ͵�λ��Ȳ�����ѹǿ�ĳ˻����һ��ϸ����anterior-lateral 
%%           domain��cortex��ȡ����Ⱥ͵�λ��Ȳ�����ѹǿ�ĳ˻��ĺͼ�ȥ������ȶ�Ӧ��������
lateral_tension_situation1=zeros(1000,9,2);  % �������඼��ϸ���Ҷ�ֵ���ݵ�lumen���ܵ������Լ���ϸ��������������������ֵ���ڶ���ά�ȴ���ÿ��lumen��Ӧ��lumen length/cell length
filenamearray_lateraltension_situation1=cell(1000,9,2);  % �������ÿһ��lateral_tension��Ӧ������ϸ����·���ļ���
i_tension_position=zeros(9,1);  % ָʾlateral_tension_situation1�����д��λ��
for i=1:9
    for j=1:cell_number_multicell(i)
        framepath_filename_1=mat2str(cell2mat(filenamearray_multicell_anterior(j,i)));  % ��ȡanterior-lateral domain��·�����ļ���
        framepath_filename_1=framepath_filename_1(2:length(framepath_filename_1)-1);  % ɾ������ĵ�����
        

        % ��ȡ·�����ļ����а����ĸ��ļ�������(month+day)����Ƭ����ϸ����
        S=regexp(framepath_filename_1,'\','split');  % ����ļ����ַ���
        anterior_date=mat2str(cell2mat(S(7)));  % ��ȡ·���ļ����е�����
        anterior_date=anterior_date(2:length(anterior_date)-1);  % ɾ������ĵ�����
        filename=mat2str(cell2mat(S(end)));  % ��ȡ·���ļ����е������ļ���
        if length(filename)==17  % filename����Ϊ17��ʾ��Ƭ��Ϊ��λ��
            anterior_photo_number=filename(7:8);
            anterior_cell_number=filename(9);
        elseif length(filename)==16  % filename����Ϊ16��ʾ��Ƭ��Ϊһλ��           
            anterior_photo_number=filename(7);
            anterior_cell_number=filename(8);
        end
        
        % ����filenamearray_multicell_posterior�������������anterior_lateral domain������ϸ����Ӧ��posterior_lateral domain
        logic_filenamematch=0;  % �Ƿ���ڶ�Ӧϸ�����߼�ֵ
        for k=1:cell_number_multicell(i)
            framepath_filename_3=mat2str(cell2mat(filenamearray_multicell_posterior(k,i)));  % ��ȡanterior-lateral domain��·�����ļ���
            framepath_filename_3=framepath_filename_3(2:length(framepath_filename_3)-1);  % ɾ������ĵ�����
            
            % ��ȡ·�����ļ����а����ĸ��ļ�������(month+day)����Ƭ����ϸ����
            S=regexp(framepath_filename_3,'\','split');  % ����ļ����ַ���
            posterior_date=mat2str(cell2mat(S(7)));  % ��ȡ·���ļ����е�����
            posterior_date=posterior_date(2:length(posterior_date)-1);  % ɾ������ĵ�����
            filename=mat2str(cell2mat(S(end)));  % ��ȡ·���ļ����е������ļ���
            if length(filename)==17  % filename����Ϊ17��ʾ��Ƭ��Ϊ��λ��
                posterior_photo_number=filename(7:8);
                posterior_cell_number=filename(9);
            elseif length(filename)==16  % filename����Ϊ16��ʾ��Ƭ��Ϊһλ��
                posterior_photo_number=filename(7);
                posterior_cell_number=filename(8);
            end
            
            % �ж�����(month+day)����Ƭ����ϸ�����Ƿ��Ӧ
            if ~sum(abs(anterior_date-posterior_date)) && ~sum(abs(anterior_photo_number-posterior_photo_number))  % ���ں���Ƭ���Ƿ���ͬ
                if anterior_cell_number==posterior_cell_number+1  % anterior part�Ƿ���posterior part��ǰһ��ϸ��
                    logic_filenamematch=1;
                    break
                end
            end
        end
        
        % �����lumen���ܵ�����������lateral domain������������
        if logic_filenamematch==1  % ��ʾ��ƥ���һ���ļ���
            i_tension_position(i)=i_tension_position(i)+1;
            lateral_tension_situation1(i_tension_position(i),i,1)=(anterior_thickness_multicell(j,i,1)-pre_average_anteriorthickness(i))*anterior_thickness_multicell(j,i,2)+(posterior_thickness_multicell(k,i,1)-pre_average_posteriorthickness(i))*posterior_thickness_multicell(k,i,2);  % ��������������ֵ
            lateral_tension_situation1(i_tension_position(i),i,1)=lateral_tension_situation1(i_tension_position(i),i,1)/(pre_average_anteriorthickness(i)*pre_average_anteriorlength(i)+pre_average_posteriorthickness(i)*pre_average_posteriorlength(i))*2;  % �õ���Ժ�ȵ�������
%             lateral_tension_situation1(i_tension_position(i),i,1)=(anterior_thickness_multicell(j,i,1)-pre_average_anteriorthickness(i))+(posterior_thickness_multicell(k,i,1)-pre_average_posteriorthickness(i));  % ��������������ֵ
%             lateral_tension_situation1(i_tension_position(i),i,1)=lateral_tension_situation1(i_tension_position(i),i,1)/(pre_average_anteriorthickness(i)+pre_average_posteriorthickness(i))*2;  % �õ���Ժ�ȵ�������
            if lateral_tension_situation1(i_tension_position(i),i,1)<0
                lateral_tension_situation1(i_tension_position(i),i,1)=0;
            end
            lateral_tension_situation1(i_tension_position(i),i,2)=anterior_thickness_multicell(j,i,3);  % ��Ӧ��lumen length/cell length��ֵ
            filenamearray_lateraltension_situation1{i_tension_position(i),i,1}=framepath_filename_1;
            filenamearray_lateraltension_situation1{i_tension_position(i),i,2}=framepath_filename_3;
        end
    end
end




%% Step 5.2���ټ�����е�����ϸ���Ҷ�ֵ���ݵ�lumen���ܵ������Լ���ϸ����������������������С����������С���ڸ�ϸ����
%%           ���anterior-lateral domain��cortex��ȡ����Ⱥ͵�λ��Ȳ�����ѹǿ�ĳ˻���posterior-lateral domain��cor-
%%           tex��ȡ����Ⱥ͵�λ��Ȳ�����ѹǿ�ĳ˻��ĺͼ�ȥ������ȶ�Ӧ��������
lateral_tension_situation2=zeros(1000,9,2);  % ����һ����ϸ���Ҷ�ֵ���ݵ�lumen���ܵ������Լ���ϸ��������������������ֵ���ڶ���ά�ȴ���ÿ��lumen��Ӧ��lumen length/cell length
filenamearray_lateraltension_situation2=cell(1000,9);  % �������ÿһ��lateral_force��Ӧ������ϸ����·���ļ���
i_tension2_position=zeros(9,1);  % ָʾlateral_tension_situation2�����д��λ��

% ���������������ļ���ϸ���У�δ�����lateral_tension��ϸ������Ӧ��lumen�ܵ���������ֵ������ȥ������ȶ�Ӧ��������
for i=1:9
    % ����anterior-/posterior-lateral domain����Ӧ�Ļ�δ�������ϸ��
    for j=1:cell_number_multicell(i)
        framepath_filename_1=mat2str(cell2mat(filenamearray_multicell_anterior(j,i)));  % ��ȡanterior-lateral domain��·�����ļ���
        framepath_filename_1=framepath_filename_1(2:length(framepath_filename_1)-1);  % ɾ������ĵ�����
        
        % �жϸ��ļ����Ƿ���filenamearray_lateraltension_situation1�����в�����
        logic_filenamematch=1;  % �ļ����Ƿ���ڵ��߼�ֵ
        for k=1:i_tension_position(i)
            if length(framepath_filename_1)==length(filenamearray_lateraltension_situation1{k,i,1}) && ~sum(abs(framepath_filename_1-filenamearray_lateraltension_situation1{k,i,1}))
                logic_filenamematch=0;
                break
            end
        end
        
        % �����lumen���ܵ�����������lateral domain������������
        if logic_filenamematch==1
            i_tension2_position(i)=i_tension2_position(i)+1;
            lateral_tension_situation2(i_tension2_position(i),i,1)=(anterior_thickness_multicell(j,i,1)-pre_average_anteriorthickness(i))*anterior_thickness_multicell(j,i,2)+(posterior_thickness_multicell(j,i,1)-pre_average_posteriorthickness(i))*posterior_thickness_multicell(j,i,2);  % ��������������ֵ
            lateral_tension_situation2(i_tension2_position(i),i,1)=lateral_tension_situation2(i_tension2_position(i),i,1)/(pre_average_anteriorthickness(i)*pre_average_anteriorlength(i)+pre_average_posteriorthickness(i)*pre_average_posteriorlength(i))*2;  % �õ���Ժ�ȵ�������
%             lateral_tension_situation2(i_tension2_position(i),i,1)=(anterior_thickness_multicell(j,i,1)-pre_average_anteriorthickness(i))+(posterior_thickness_multicell(j,i,1)-pre_average_posteriorthickness(i));  % ��������������ֵ
%             lateral_tension_situation2(i_tension2_position(i),i,1)=lateral_tension_situation2(i_tension2_position(i),i,1)/(pre_average_anteriorthickness(i)+pre_average_posteriorthickness(i))*2;  % �õ���Ժ�ȵ�������
            if lateral_tension_situation2(i_tension2_position(i),i,1)<0
                lateral_tension_situation2(i_tension2_position(i),i,1)=0;
            end
            lateral_tension_situation2(i_tension2_position(i),i,2)=anterior_thickness_multicell(j,i,3);  % ��Ӧ��lumen length/cell length��ֵ
            filenamearray_lateraltension_situation2{i_tension2_position(i),i}=framepath_filename_1;
        end
    end
end

% �����е��������ļ���ϸ������Ӧ��lumen�ܵ���������ֵ������ȥ������ȶ�Ӧ��������
for i=1:9
    for j=1:cell_number_monocell(i)
        framepath_filename_danyi=mat2str(cell2mat(filenamearray_monocell(j,i)));  % ��ȡanterior-lateral domain��·�����ļ���
        framepath_filename_danyi=framepath_filename_danyi(2:length(framepath_filename_danyi)-1);  % ɾ������ĵ�����
        
        i_tension2_position(i)=i_tension2_position(i)+1;
        lateral_tension_situation2(i_tension2_position(i),i,1)=(anterior_thickness_monocell(j,i,1)-pre_average_anteriorthickness(i))*anterior_thickness_monocell(j,i,2)+(posterior_thickness_monocell(j,i,1)-pre_average_posteriorthickness(i))*posterior_thickness_monocell(j,i,2);  % ��������������ֵ
        lateral_tension_situation2(i_tension2_position(i),i,1)=lateral_tension_situation2(i_tension2_position(i),i,1)/(pre_average_anteriorthickness(i)*pre_average_anteriorlength(i)+pre_average_posteriorthickness(i)*pre_average_posteriorlength(i))*2;  % �õ���Ժ�ȵ�������
        if lateral_tension_situation2(i_tension2_position(i),i,1)<0
            lateral_tension_situation2(i_tension2_position(i),i,1)=0;
        end
        lateral_tension_situation2(i_tension2_position(i),i,2)=anterior_thickness_monocell(j,i,3);  % ��Ӧ��lumen length/cell length��ֵ
        filenamearray_lateraltension_situation2{i_tension2_position(i),i}=framepath_filename_danyi;
    end
end



%% Step 5.3�����lumen���ܵ������Լ���ϸ����������������������lumen length/cell length��С��������
% �����඼��ϸ�����ݵ�lumen�͵�����ϸ�����ݵ�lumen����������ϲ�
lateral_tension_all=zeros(1000,9,2);  % �������������ϲ�����
for i=1:9
    for j=1:i_tension_position(i)
        lateral_tension_all(j,i,:)=lateral_tension_situation1(j,i,:);
    end
    for j=1:i_tension2_position(i)
        lateral_tension_all(j+i_tension_position(i),i,:)=lateral_tension_situation2(j,i,:);
    end
end
lateral_tension_all(:,:,1)=lateral_tension_all(:,:,1)./2;

% ��lateral_tension_all����lumen length/cell length��С��������(ð������)
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



%% Step 5.4�����ݹ�һ����TD(lumen length/cell length)��ʱ��t�Ĺ۲�������ϱ仯��ϵ�����lumen���ܵ������Լ���ϸ������
%%           ����������������ʱ��t�ĺ�����ϵ
% ���ݹ۲���R��theta��ʱ���ʵ�ʱ仯��ϵ
% ����ʵ�ʹ۲����transverse diameter(TD)��longitudinal radius(LR)��ʱ��ı仯����
TD=[3.248130325	3.697268921	4.042205695	4.875820261	5.623190217	5.105764223	5.795637773	5.996833531	6.456763119	6.571755932	6.485511322	6.370518509	6.399266712	6.485511322	6.370518509	6.600504135	6.284315564	6.226819157	6.226819157	6.830448097	7.37658063	7.577818053	8.411390955	8.267691603	8.900068746	9.474949482	10.02108202	10.4235152	11.14213695	11.77451409	12.57938046	12.86682082	13.38420515	13.64293898	14.21781972	14.36151907	14.5915047	15.10888903	15.14842822	15.05139262	16];
LR=[0.8114597 0.88031998 1.012697123 0.632377143 0.966522926 1.167760348 1.523446451 1.943878508 2.08757786	2.288815282	2.202570673	2.949940629	3.381080349	3.668520717	4.214694915	4.818323855	4.703331042	5.174051622	5.192008833	5.651938421	5.795637773	5.968126992	6.054329938	6.140574547	6.169322751	5.910630586	5.853134179	6.111826344	6.284315564	5.996833531	6.083078141	5.996833531	5.824385976	5.863883507	5.853134179	5.853134179	5.795637773	5.939378789	5.738141366	6.111826344	5.996833531];
TD2=[TD 16.1 16.1 16.1 16.1 16.1 16.1 16.1 16.1 16.1 16.1];
t=linspace(0,120,length(TD))+6;  % �۲����ݵ�ʱ��t
t2=[t t(end) t(end) t(end) t(end) t(end) t(end) t(end) t(end) t(end) t(end)];
tlength=1000;  % ���ʱ����ľ���
t_fit=linspace(0,126,tlength);  % ���ʱ����
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

% ����TD��LR��R��theta�Ĺ�ϵ��ͨ���������������R��theta��ʱ��ı仯
R=zeros(1,length(t_fit));
theta=zeros(1,length(t_fit));
for i=1:length(t_fit)
    wucha_best=1000;  % ��С���ֵ
    for theta_search=0:pi/1000:pi
        R_search=TD_fit(i)/2/sin(theta_search);
        if abs(R_search+R_search*cos(theta_search)-LR_fit(i))<wucha_best
            wucha_best=abs(R_search+R_search*cos(theta_search)-LR_fit(i));
            R(i)=R_search;
            theta(i)=theta_search*180/pi;
        end
    end
end


% ��lumen length/cell length(lateral_tension_all(:,:,2))ת��Ϊ��Ӧ��ʱ��t��������t_rectify������
t_rectify=zeros(1,tlength);  % �������׹۲����ݺ��ʱ��t
t_fit_rectify=linspace(0,126,tlength);  % �������׹۲����ݺ��������ߺ�����t
TD_rectify=zeros(1,tlength);  % �������TD
TD_fit_rectify=zeros(1,tlength);  % �������TD��Ϻ���
TD_jump=6.2;  % lumen length/cell length����Ծֵ

i_rectify_position=0;  % ��ֵ��rectify�����е�λ��
for i=1:tlength  % ȥ��TD_jump��Χ��������
    if TD_fit(i)>=TD_jump+1 || TD_fit(i)<=TD_jump-1
        i_rectify_position=i_rectify_position+1;
        t_rectify(i_rectify_position)=t_fit(i);
        TD_rectify(i_rectify_position)=TD_fit(i);
    end
end
% ��TD_rectify�������
% N_rectify=6;
% p=polyfit(t_rectify(1:i_rectify_position),TD_rectify(1:i_rectify_position),N_rectify);
% for k=1:N_rectify+1
%     TD_fit_rectify=TD_fit_rectify+p(k).*t_fit_rectify.^(N_rectify-k+1);
% end
% ��TD_rectify���в�ֵ
TD_fit_rectify=interp1(t_rectify(1:i_rectify_position),TD_rectify(1:i_rectify_position),t_fit_rectify,'pchip');  % ����������ֵ

% plot(t_rectify(1:i_rectify_position),TD_rectify(1:i_rectify_position),'.','MarkerSize',20)
% hold on
% plot(t_fit_rectify,TD_fit_rectify,'LineWidth',3)


% �ҳ�lateral_tension_all������lumen length/cell length��Ӧ��ʱ��t
t_observation=zeros(1000,1);  % ʵ��۲��lumen length/cell length��Ӧ��ʱ��t
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


% ���lumen���ܵ������Լ���ϸ����������������������ʱ��t�ĺ�����ϵ
N=3;  % ��Ͻ���
lateral_tension_fit=zeros(N+1,9);  % ������϶���ʽϵ������
delta_tensionfit=zeros(1000,9);  % ����ѹǿ��ϵı�׼��
i=1;
[p_fit,S]=polyfit(t_observation(1:i_tension_position(i)+i_tension2_position(i)),lateral_tension_all(1:i_tension_position(i)+i_tension2_position(i),i,1),N);  % anterior-lateral domain length�Ķ���ʽ���
[y_polyval,delta_tensionfit(1:i_tension_position(i)+i_tension2_position(i),i)]=polyval(p_fit,t_observation(1:i_tension_position(i)+i_tension2_position(i)),S);
for k=1:N+1
    lateral_tension_fit(k,i)=p_fit(k);
end



%% ����Lumen Contact Angle��lateral tension�������
win=21;  % ����ƽ��ֵ�Ŀ��
ContactAngle_LateralTension_MA=zeros(cell_number_match,3);  % ����lateral tension�Ļ���ƽ��ֵ
ContactAngle_LateralTension_Mstd=zeros(cell_number_match,3);  % ����lateral tension�Ļ�����׼��
ContactAngle_LateralTension_MA2=zeros(cell_number_match,3);  % ����lateral tension�Ķ��λ���ƽ��ֵ
ContactAngle_LateralTension_Mstd2=zeros(cell_number_match,3);  % ����lateral tension�Ķ��λ�����׼��

for i=1:cell_number_match
    for j=1:2
        if i<=(win-1)/2  % ��tΪ��ʼ����ʱ��
            ContactAngle_LateralTension_MA(i,j)=sum(ContactAngle_LateralTension(1:2*i-1,j))/(2*i-1);
            ContactAngle_LateralTension_Mstd(i,j)=std(ContactAngle_LateralTension(1:2*i-1,j));
        elseif i<i_tension_position(1)+i_tension2_position(1)-(win-1)/2+1  % ��tΪ�м�ʱ��
            ContactAngle_LateralTension_MA(i,j)=sum(ContactAngle_LateralTension(i-(win-1)/2:i+(win-1)/2,j))/win;
            ContactAngle_LateralTension_Mstd(i,j)=std(ContactAngle_LateralTension(i-(win-1)/2:i+(win-1)/2,j));
        else  % ��tΪĩβ����ʱ��
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
        if i<=(win-1)/2  % ��tΪ��ʼ����ʱ��
            ContactAngle_LateralTension_MA2(i,j)=sum(ContactAngle_LateralTension_MA(1:2*i-1,j))/(2*i-1);
            ContactAngle_LateralTension_Mstd2(i,j)=std(ContactAngle_LateralTension_MA(1:2*i-1,j));
        elseif i<i_tension_position(1)+i_tension2_position(1)-(win-1)/2+1  % ��tΪ�м�ʱ��
            ContactAngle_LateralTension_MA2(i,j)=sum(ContactAngle_LateralTension_MA(i-(win-1)/2:i+(win-1)/2,j))/win;
            ContactAngle_LateralTension_Mstd2(i,j)=std(ContactAngle_LateralTension_MA(i-(win-1)/2:i+(win-1)/2,j));
        else  % ��tΪĩβ����ʱ��
            temp=i_tension_position(1)+i_tension2_position(1)-i+1;
            ContactAngle_LateralTension_MA2(i,j)=sum(ContactAngle_LateralTension_MA(cell_number_match-2*temp+2:cell_number_match,j))/(2*temp-1);
            ContactAngle_LateralTension_Mstd2(i,j)=std(ContactAngle_LateralTension_MA(cell_number_match-2*temp+2:cell_number_match,j));
        end
    end
end
ContactAngle_LateralTension_MA2=ContactAngle_LateralTension_MA2((win-1)/2+1:cell_number_match-(win-1)/2,:);
ContactAngle_LateralTension_Mstd2=ContactAngle_LateralTension_Mstd2((win-1)/2+1:cell_number_match-(win-1)/2,:);
cell_number_match=cell_number_match-(win-1);

% ���Lumen Contact Angle��lateral tension�����Թ�ϵ
N=1;
[p,S]=polyfit(ContactAngle_LateralTension_MA2(1:cell_number_match,2),ContactAngle_LateralTension_MA2(1:cell_number_match,1),N); 
x_fit=linspace(ContactAngle_LateralTension_MA2(1,2)-0.2,ContactAngle_LateralTension_MA2(cell_number_match,2)+0.2,cell_number_match);
[y_fit,delta]=polyval(p,x_fit,S);

% ����Lumen Contact Angle��lateral tension������Ŷ�
x_fit_interp=linspace(x_fit(1),x_fit(end),10000);
y_fit_interp=interp1(x_fit,y_fit,x_fit_interp,'linear');  % ����������ֵ
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

% ����Lumen Contact Angle��lateral tension�����ϵ��
A=corrcoef(ContactAngle_LateralTension_MA2(1:cell_number_match,1),ContactAngle_LateralTension_MA2(1:cell_number_match,2));
ContactAngle_LateralTension_time_corrcoef=A(1,2);

