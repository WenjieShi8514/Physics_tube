%% Assign the values of Transverse Diameter (radial) and Longitudinal Radius (axial)
Cr_ZO1_WT=zeros(10,2,20,2);  % Lumen parameters of wild-type ZO1 overexpression. The first dimension represents data of different lumens in one image; the second dimension represents different parameters of a single lumen (vertical length, transverse length),
                             % the third dimension represents different images; the fourth dimension indicates experimental or control group (1 for lumens with fluorescence, 2 for lumens without fluorescence).
Cr_ZO1_GuK=zeros(10,2,20,2);  % ZO1 U5-GuK缺失突变体过表达lumen参数
Cr_ZO1_ABR=zeros(10,2,20,2);  % ZO1 ABR缺失突变体过表达lumen参数



%% Solve normalized volume ratio and lumen opening length ratio
lumen_volume=zeros(200,3);  % Normalized lumen volume. The second dimension represents WT, GuK, ABR
lumen_volume_control=zeros(200,3);  % Lumen volume of cells without fluorescence (control group)
lumen_TD=zeros(200,3);  % Normalized lumen opening length. The second dimension represents WT, GuK, ABR
lumen_TD_control=zeros(200,3);  % Lumen opening length of cells without fluorescence (control group)
lumen_LR=zeros(200,3);  % Normalized transverse (A-P) radius of the lumen. The second dimension represents ZO1 WT overexpression group and cells without fluorescence (control group).
lumen_LR_control=zeros(200,3);
i_experiment_array=zeros(20,3);  % Number of experimental samples per image
i_control_array=zeros(20,3);  % Number of control samples per image

i_number1=1;  % Indicates the position when inputting to lumen_volume and lumen_openinglength matrices
for k=1:length(Cr_ZO1_WT)
    % First, check if the image exists
    if Cr_ZO1_WT(1,1,k,1)==0
        break
    end
    
    % Determine the number of lumens with fluorescence and without fluorescence in the experimental and control groups for this image.
    for i=1:10
        if Cr_ZO1_WT(i,1,k,1)==0
            i_experiment_array(k,1)=i-1;  % Number of experimental samples per image
            break
        end
    end
    for i=1:10
        if Cr_ZO1_WT(i,1,k,2)==0
            i_control_array(k,1)=i-1;  % Number of control samples per image
            break
        end
    end
    
    % 计算Average lumen volume of the control group.比和lumen opening length比
    opening_length_control=sum(Cr_ZO1_WT(1:i_control_array(k,1),1,k,2))/i_control_array(k,1);  % Average lumen opening length of the control group
    volume_control=0;
    for i=1:i_control_array(k,1)
        wucha_best=1000;
        for theta=0:pi/1000:pi*999/1000
            wucha=abs(Cr_ZO1_WT(i,2,k,2)*sin(theta)-Cr_ZO1_WT(i,1,k,2)*(1+cos(theta)));
            if wucha<wucha_best
                wucha_best=wucha;
                theta_best=theta;
            end
        end
        r_best=Cr_ZO1_WT(i,1,k,2)/sin(theta_best);
        volume_control=volume_control+4/3*pi*r_best^3*(1+3/2*cos(theta_best)-1/2*cos(theta_best)^3);
        lumen_volume_control(i+sum(i_control_array(:,1))-i_control_array(k,1),1)=4/3*pi*r_best^3*(1+3/2*cos(theta_best)-1/2*cos(theta_best)^3);  % Lumen volume of cells without fluorescence (control group)比赋值
        lumen_TD_control(i+sum(i_control_array(:,1))-i_control_array(k,1),1)=Cr_ZO1_WT(i,1,k,2);  % Lumen opening length of cells without fluorescence (control group)比
        lumen_LR_control(i+sum(i_control_array(:,1))-i_control_array(k,1),1)=Cr_ZO1_WT(i,2,k,2);  % Lumen opening length of cells without fluorescence (control group)比
    end
    volume_control=volume_control/i_control_array(k,1);  % Average lumen volume of the control group.
    
    % 计算各实验样本的归一化体积比和lumen opening length比
    for i=1:i_experiment_array(k,1)
        lumen_TD(i_number1,1)=Cr_ZO1_WT(i,1,k,1)/opening_length_control;  % Lumen opening length ratio of experimental samples.
        lumen_LR(i_number1,1)=Cr_ZO1_WT(i,2,k,1)/opening_length_control;  % Lumen opening length ratio of experimental samples.
        wucha_best=1000;
        for theta=0:pi/1000:pi*999/1000
            wucha=abs(Cr_ZO1_WT(i,2,k,1)*sin(theta)-Cr_ZO1_WT(i,1,k,1)*(1+cos(theta)));
            if wucha<wucha_best
                wucha_best=wucha;
                theta_best=theta;
            end
        end
        r_best=Cr_ZO1_WT(i,1,k,1)/sin(theta_best);
        lumen_volume(i_number1,1)=4/3*pi*r_best^3*(1+3/2*cos(theta_best)-1/2*cos(theta_best)^3)/volume_control;  % Volume ratio of experimental samples.
        i_number1=i_number1+1;
    end
end

% Solve the ratio for ZO1 U5-GuK
i_number2=1;  % Indicates the position when inputting to lumen_volume and lumen_openinglength matrices
for k=1:length(Cr_ZO1_GuK)
    % First, check if the image exists
    if Cr_ZO1_GuK(1,1,k,1)==0
        break
    end
    
    % Determine the number of lumens with fluorescence and without fluorescence in the experimental and control groups for this image.
    for i=1:10
        if Cr_ZO1_GuK(i,1,k,1)==0
            i_experiment_array(k,2)=i-1;  % 实验组样本数
            break
        end
    end
    for i=1:10
        if Cr_ZO1_GuK(i,1,k,2)==0
            i_control_array(k,2)=i-1;  % 对照组样本数
            break
        end
    end
    
    % 计算Average lumen volume of the control group.比和lumen opening length比
    opening_length_control=sum(Cr_ZO1_GuK(1:i_control_array(k,2),1,k,2))/i_control_array(k,2);  % Average lumen opening length of the control group
    lumen_longitundialradius_control=sum(Cr_ZO1_GuK(1:i_control_array(k,2),2,k,2))/i_control_array(k,2);  % 对照组平均lumen longitundial radius
    volume_control=0;
    for i=1:i_control_array(k,2)
        wucha_best=1000;
        for theta=pi/1000:pi/1000:pi*999/1000
            wucha=abs(Cr_ZO1_GuK(i,1,k,2)/sin(theta)-Cr_ZO1_GuK(i,2,k,2)/(1+cos(theta)));
            if wucha<wucha_best
                wucha_best=wucha;
                theta_best=theta;
            end
        end
        r_best=Cr_ZO1_GuK(i,1,k,2)/sin(theta_best)/2;
        volume_control=volume_control+4/3*pi*r_best^3*(1+3/2*cos(theta_best)-1/2*cos(theta_best)^3);
        lumen_volume_control(i+sum(i_control_array(:,2))-i_control_array(k,2),2)=4/3*pi*r_best^3*(1+3/2*cos(theta_best)-1/2*cos(theta_best)^3);  % Lumen volume of cells without fluorescence (control group)比赋值
        lumen_TD_control(i+sum(i_control_array(:,2))-i_control_array(k,2),2)=Cr_ZO1_GuK(i,1,k,2);  % Lumen opening length of cells without fluorescence (control group)比
        lumen_LR_control(i+sum(i_control_array(:,2))-i_control_array(k,2),2)=Cr_ZO1_GuK(i,2,k,2);  % Lumen opening length of cells without fluorescence (control group)比
    end
    volume_control=volume_control/i_control_array(k,2);  % Average lumen volume of the control group.
    
    % 计算各实验样本的归一化体积比和lumen opening length比
    for i=1:i_experiment_array(k,2)
        lumen_TD(i_number2,2)=Cr_ZO1_GuK(i,1,k,1)/opening_length_control;  % Lumen opening length ratio of experimental samples.
        lumen_LR(i_number2,2)=Cr_ZO1_GuK(i,2,k,1)/lumen_longitundialradius_control;  % 实验样本的lumen longitundial radius比
        wucha_best=1000;
        for theta=pi/1000:pi/1000:pi*999/1000
            wucha=abs(Cr_ZO1_GuK(i,1,k,1)/sin(theta)-Cr_ZO1_GuK(i,2,k,1)/(1+cos(theta)));
            if wucha<wucha_best
                wucha_best=wucha;
                theta_best=theta;
            end
        end
        r_best=Cr_ZO1_GuK(i,1,k,1)/sin(theta_best)/2;
        lumen_volume(i_number2,2)=4/3*pi*r_best^3*(1+3/2*cos(theta_best)-1/2*cos(theta_best)^3)/volume_control;  % Volume ratio of experimental samples.
        i_number2=i_number2+1;
    end
end

% 求解ZO1 ABR比值
i_number3=1;  % Indicates the position when inputting to lumen_volume and lumen_openinglength matrices
for k=1:length(Cr_ZO1_ABR)
    % First, check if the image exists
    if Cr_ZO1_ABR(1,1,k,1)==0
        break
    end
    
    % Determine the number of lumens with fluorescence and without fluorescence in the experimental and control groups for this image.
    for i=1:10
        if Cr_ZO1_ABR(i,1,k,1)==0
            i_experiment_array(k,3)=i-1;  % 实验组样本数
            break
        end
    end
    for i=1:10
        if Cr_ZO1_ABR(i,1,k,2)==0
            i_control_array(k,3)=i-1;  % 对照组样本数
            break
        end
        if Cr_ZO1_ABR(10,1,k,2)~=0
            i_control_array(k,3)=10;  % 对照组样本数
        end
    end
    
    % 计算Average lumen volume of the control group.比和lumen opening length比
    opening_length_control=sum(Cr_ZO1_ABR(1:i_control_array(k,3),1,k,2))/i_control_array(k,3);  % Average lumen opening length of the control group
    lumen_longitundialradius_control=sum(Cr_ZO1_ABR(1:i_control_array(k,3),2,k,2))/i_control_array(k,3);  % 对照组平均lumen longitundial radius
    volume_control=0;
    for i=1:i_control_array(k,3)
        wucha_best=1000;
        for theta=pi/1000:pi/1000:pi*999/1000
            wucha=abs(Cr_ZO1_ABR(i,1,k,2)/sin(theta)-Cr_ZO1_ABR(i,2,k,2)/(1+cos(theta)));
            if wucha<wucha_best
                wucha_best=wucha;
                theta_best=theta;
            end
        end
        r_best=Cr_ZO1_ABR(i,1,k,2)/sin(theta_best)/2;
        volume_control=volume_control+4/3*pi*r_best^3*(1+3/2*cos(theta_best)-1/2*cos(theta_best)^3);
        lumen_volume_control(i+sum(i_control_array(:,3))-i_control_array(k,3),3)=4/3*pi*r_best^3*(1+3/2*cos(theta_best)-1/2*cos(theta_best)^3);  % Lumen volume of cells without fluorescence (control group)比赋值
        lumen_TD_control(i+sum(i_control_array(:,3))-i_control_array(k,3),3)=Cr_ZO1_ABR(i,1,k,2);  % Lumen opening length of cells without fluorescence (control group)比
        lumen_LR_control(i+sum(i_control_array(:,3))-i_control_array(k,3),3)=Cr_ZO1_ABR(i,2,k,2);  % Lumen opening length of cells without fluorescence (control group)比
    end
    volume_control=volume_control/i_control_array(k,3);  % Average lumen volume of the control group.
    
    % 计算各实验样本的归一化体积比和lumen opening length比
    for i=1:i_experiment_array(k,3)
        lumen_TD(i_number3,3)=Cr_ZO1_ABR(i,1,k,1)/opening_length_control;  % Lumen opening length ratio of experimental samples.
        lumen_LR(i_number3,3)=Cr_ZO1_ABR(i,2,k,1)/lumen_longitundialradius_control;  % 实验样本的lumen longitundial radius比
        wucha_best=1000;
        for theta=pi/1000:pi/1000:pi*999/1000
            wucha=abs(Cr_ZO1_ABR(i,1,k,1)/sin(theta)-Cr_ZO1_ABR(i,2,k,1)/(1+cos(theta)));
            if wucha<wucha_best
                wucha_best=wucha;
                theta_best=theta;
            end
        end
        r_best=Cr_ZO1_ABR(i,1,k,1)/sin(theta_best)/2;
        lumen_volume(i_number3,3)=4/3*pi*r_best^3*(1+3/2*cos(theta_best)-1/2*cos(theta_best)^3)/volume_control;  % Volume ratio of experimental samples.
        i_number3=i_number3+1;
    end
end



%% Calculate the mean and standard deviation of the lumen volume for each individual in the ZO1 WT overexpression group and the control group
individual_TD=zeros(20,3,2);  % ZO1 WT组的ZO1 WT overexpression和control group的lumen opening length均值
individual_LR=zeros(20,3,2);  % ZO1 WT组的ZO1 WT overexpression和control group的LR均值
individual_volume=zeros(20,3,2);  % ZO1 WT组的ZO1 WT overexpression和control group的lumen volume均值
std_individual_TD=zeros(20,3,2);  % ZO1 WT组的ZO1 WT overexpression和control group的lumen opening length标准差
std_individual_LR=zeros(20,3,2);  % ZO1 WT组的ZO1 WT overexpression和control group的LR标准差
std_individual_volume=zeros(20,3,2);  % ZO1 WT组的ZO1 WT overexpression和control group的lumen volume标准差
k_final_array=zeros(1,3);
for i_ZO1_type=1:3
for k=1:20  % 分别表示ZO1 WT和control group
    if i_experiment_array(k,i_ZO1_type)==0  % First, check if the image exists
        k_final=k-1;  % Record the number of images.
        break
    end
    
    % 计算ZO1 WT/DN overexpression每张照片的归一化lumen opening length和lumen volume
    individual_TD(k,i_ZO1_type,1)=sum(lumen_TD_control(1+sum(i_control_array(1:k,i_ZO1_type))-i_control_array(k,i_ZO1_type):sum(i_control_array(1:k,i_ZO1_type)),i_ZO1_type))/i_control_array(k,i_ZO1_type);
    individual_LR(k,i_ZO1_type,1)=sum(lumen_LR_control(1+sum(i_control_array(1:k,i_ZO1_type))-i_control_array(k,i_ZO1_type):sum(i_control_array(1:k,i_ZO1_type)),i_ZO1_type))/i_control_array(k,i_ZO1_type);
    individual_volume(k,i_ZO1_type,1)=sum(lumen_volume_control(1+sum(i_control_array(1:k,i_ZO1_type))-i_control_array(k,i_ZO1_type):sum(i_control_array(1:k,i_ZO1_type)),i_ZO1_type))/i_control_array(k,i_ZO1_type);
    std_individual_TD_temp=std(lumen_TD_control(1+sum(i_control_array(1:k,i_ZO1_type))-i_control_array(k,i_ZO1_type):sum(i_control_array(1:k,i_ZO1_type)),i_ZO1_type));
    std_individual_LR_temp=std(lumen_LR_control(1+sum(i_control_array(1:k,i_ZO1_type))-i_control_array(k,i_ZO1_type):sum(i_control_array(1:k,i_ZO1_type)),i_ZO1_type));
    std_individual_volume_temp=std(lumen_volume_control(1+sum(i_control_array(1:k,i_ZO1_type))-i_control_array(k,i_ZO1_type):sum(i_control_array(1:k,i_ZO1_type)),i_ZO1_type));
    
    individual_TD(k,i_ZO1_type,2)=sum(lumen_TD(1+sum(i_experiment_array(1:k,i_ZO1_type))-i_experiment_array(k,i_ZO1_type):sum(i_experiment_array(1:k,i_ZO1_type)),i_ZO1_type))/i_experiment_array(k,i_ZO1_type);
    individual_LR(k,i_ZO1_type,2)=sum(lumen_LR(1+sum(i_experiment_array(1:k,i_ZO1_type))-i_experiment_array(k,i_ZO1_type):sum(i_experiment_array(1:k,i_ZO1_type)),i_ZO1_type))/i_experiment_array(k,i_ZO1_type);
    std_individual_TD(k,i_ZO1_type,2)=std(lumen_TD(1+sum(i_experiment_array(1:k,i_ZO1_type))-i_experiment_array(k,i_ZO1_type):sum(i_experiment_array(1:k,i_ZO1_type)),i_ZO1_type))*individual_TD(k,i_ZO1_type,1);
    std_individual_LR(k,i_ZO1_type,2)=std(lumen_LR(1+sum(i_experiment_array(1:k,i_ZO1_type))-i_experiment_array(k,i_ZO1_type):sum(i_experiment_array(1:k,i_ZO1_type)),i_ZO1_type))*individual_LR(k,i_ZO1_type,1);
    individual_volume(k,i_ZO1_type,2)=sum(lumen_volume(1+sum(i_experiment_array(1:k,i_ZO1_type))-i_experiment_array(k,i_ZO1_type):sum(i_experiment_array(1:k,i_ZO1_type)),i_ZO1_type))/i_experiment_array(k,i_ZO1_type);
    std_individual_volume(k,i_ZO1_type,2)=std(lumen_volume(1+sum(i_experiment_array(1:k,i_ZO1_type))-i_experiment_array(k,i_ZO1_type):sum(i_experiment_array(1:k,i_ZO1_type)),i_ZO1_type))*individual_volume(k,i_ZO1_type,1);
end
k_final_array(i_ZO1_type)=k_final;
end


%% Calculate the mean and standard deviation of volume and length for the control group and the experimental group
average_TD=zeros(1,3);  % 长度均值
average_LR=zeros(1,3);  % LR均值
average_volume=zeros(1,3);  % 体积均值
std_TD=zeros(1,3);  % 长度标准差
std_LR=zeros(1,3);  % LR标准差
std_volume=zeros(1,3);  % 体积标准差
i_number=[i_number1-1 i_number2-1 i_number3-1];

for i=1:3
    average_TD(i)=sum(individual_TD(:,i,2))/k_final_array(i);
    average_LR(i)=sum(individual_LR(:,i,2))/k_final_array(i);
    average_volume(i)=sum(individual_volume(:,i,2))/k_final_array(i);
    std_TD(i)=std(individual_TD(1:k_final_array(i),i,2).*individual_TD(1:k_final_array(i),i,1))/(sum(individual_TD(1:k_final_array(i),i,1))/k_final_array(i));
    std_LR(i)=std(individual_LR(1:k_final_array(i),i,2).*individual_LR(1:k_final_array(i),i,1))/(sum(individual_LR(1:k_final_array(i),i,1))/k_final_array(i));
    std_volume(i)=std(individual_volume(1:k_final_array(i),i,2).*individual_volume(1:k_final_array(i),i,1))/(sum(individual_volume(1:k_final_array(i),i,1))/k_final_array(i));
end


%% Perform statistical tests
% 首先通过D’Agostino–Pearson normality test判断实验样本是否服从正态分布。若服从，则使用t-test进行检验，若不服从，则使用Mann–Whitney U test进行检验
TD_test=zeros(1,3);  % 第一维表示WT比GuK，第二维表示WT比ABR，第三维表示GuK比ABR
LR_test=zeros(1,3);
volume_test=zeros(1,3);

% lumen TD的统计检验
TD_normality_test=zeros(1,3);
for i=1:3
    % 首先通过D’Agostino–Pearson normality test判断实验样本是否服从正态分布
    X=lumen_TD(1:i_number(i),1);
    X=sort(X);
    X=X(:);
    n=length(X);

    [c,v]=hist(X,X);  % record of data in a frequency table form
    nc=find(c~=0);
    c=[v(nc) c(nc)'];

    x=c(:,1);
    f=c(:,2);
    s1=f'*x;
    s2=f'*x.^2;
    s3=f'*x.^3;
    s4=f'*x.^4;
    SS=s2-(s1^2/n);
    v=SS/(n-1);
    k3=((n*s3)-(3*s1*s2)+((2*(s1^3))/n))/((n-1)*(n-2));
    g1=k3/sqrt(v^3);
    k4=((n+1)*((n*s4)-(4*s1*s3)+(6*(s1^2)*(s2/n))-((3*(s1^4))/(n^2)))/((n-1)*(n-2)*(n-3)))-((3*(SS^2))/((n-2)*(n-3)));
    g2=k4/v^2;
    eg1=((n-2)*g1)/sqrt(n*(n-1));  % measure of skewness
    eg2=((n-2)*(n-3)*g2)/((n+1)*(n-1))+((3*(n-1))/(n+1));  % measure of kurtosis

    A=eg1*sqrt(((n+1)*(n+3))/(6*(n-2)));
    B=(3*((n^2)+(27*n)-70)*((n+1)*(n+3)))/((n-2)*(n+5)*(n+7)*(n+9));
    C=sqrt(2*(B-1))-1;
    D=sqrt(C);
    E=1/sqrt(log(D));
    F=A/sqrt(2/(C-1));
    Zg1=E*log(F+sqrt(F^2+1));

    G=(24*n*(n-2)*(n-3))/((n+1)^2*(n+3)*(n+5));
    H=((n-2)*(n-3)*g2)/((n+1)*(n-1)*sqrt(G));
    % H=((n-2)*(n-3)*abs(g2))/((n+1)*(n-1)*sqrt(G));
    J=((6*(n^2-(5*n)+2))/((n+7)*(n+9)))*sqrt((6*(n+3)*(n+5))/((n*(n-2)*(n-3))));
    K=6+((8/J)*((2/J)+sqrt(1+(4/J^2))));
    L=(1-(2/K))/(1+H*sqrt(2/(K-4)));
    Zg2=(1-(2/(9*K))-L^(1/3))/sqrt(2/(9*K));

    K2=Zg1^2+Zg2^2;  % D'Agostino-Pearson statistic
    X2=K2;  % approximation to chi-distribution
    df=2;  % degrees of freedom

    p=1-chi2cdf(X2,df);  % probability associated to the chi-squared statistic
    TD_normality_test(i)=p;
end

% 若服从正态分布，则使用t-test进行检验，若不服从正态分布，则使用Mann–Whitney U test进行检验
for i=1:3
    if i<3
        if TD_normality_test(1)>0.05 && TD_normality_test(i+1)>0.05
            [h1,p]=ttest2(lumen_TD(1:i_number(1),1),lumen_TD(1:i_number(i+1),i+1),'Alpha',0.05,'Tail','both');
            TD_test(i)=p;
        else
            p=ranksum(lumen_TD(1:i_number(1),1),lumen_TD(1:i_number(i+1),i+1),'Alpha',0.05,'Tail','both');
            TD_test(i)=p;
        end
    else
        if TD_normality_test(2)>0.05 && TD_normality_test(3)>0.05
            [h1,p]=ttest2(lumen_TD(1:i_number(2),2),lumen_TD(1:i_number(3),3),'Alpha',0.05,'Tail','both');
            TD_test(3)=p;
        else
            p=ranksum(lumen_TD(1:i_number(2),2),lumen_TD(1:i_number(3),3),'Alpha',0.05,'Tail','both');
            TD_test(3)=p;
        end
    end
end


% lumen LR的统计检验
LR_normality_test=zeros(1,3);
for i=1:3
    % 首先通过D’Agostino–Pearson normality test判断实验样本是否服从正态分布
    X=lumen_LR(1:i_number(i),1);
    X=sort(X);
    X=X(:);
    n=length(X);

    [c,v]=hist(X,X);  % record of data in a frequency table form
    nc=find(c~=0);
    c=[v(nc) c(nc)'];

    x=c(:,1);
    f=c(:,2);
    s1=f'*x;
    s2=f'*x.^2;
    s3=f'*x.^3;
    s4=f'*x.^4;
    SS=s2-(s1^2/n);
    v=SS/(n-1);
    k3=((n*s3)-(3*s1*s2)+((2*(s1^3))/n))/((n-1)*(n-2));
    g1=k3/sqrt(v^3);
    k4=((n+1)*((n*s4)-(4*s1*s3)+(6*(s1^2)*(s2/n))-((3*(s1^4))/(n^2)))/((n-1)*(n-2)*(n-3)))-((3*(SS^2))/((n-2)*(n-3)));
    g2=k4/v^2;
    eg1=((n-2)*g1)/sqrt(n*(n-1));  % measure of skewness
    eg2=((n-2)*(n-3)*g2)/((n+1)*(n-1))+((3*(n-1))/(n+1));  % measure of kurtosis

    A=eg1*sqrt(((n+1)*(n+3))/(6*(n-2)));
    B=(3*((n^2)+(27*n)-70)*((n+1)*(n+3)))/((n-2)*(n+5)*(n+7)*(n+9));
    C=sqrt(2*(B-1))-1;
    D=sqrt(C);
    E=1/sqrt(log(D));
    F=A/sqrt(2/(C-1));
    Zg1=E*log(F+sqrt(F^2+1));

    G=(24*n*(n-2)*(n-3))/((n+1)^2*(n+3)*(n+5));
    H=((n-2)*(n-3)*g2)/((n+1)*(n-1)*sqrt(G));
    % H=((n-2)*(n-3)*abs(g2))/((n+1)*(n-1)*sqrt(G));
    J=((6*(n^2-(5*n)+2))/((n+7)*(n+9)))*sqrt((6*(n+3)*(n+5))/((n*(n-2)*(n-3))));
    K=6+((8/J)*((2/J)+sqrt(1+(4/J^2))));
    L=(1-(2/K))/(1+H*sqrt(2/(K-4)));
    Zg2=(1-(2/(9*K))-L^(1/3))/sqrt(2/(9*K));

    K2=Zg1^2+Zg2^2;  % D'Agostino-Pearson statistic
    X2=K2;  % approximation to chi-distribution
    df=2;  % degrees of freedom

    p=1-chi2cdf(X2,df);  % probability associated to the chi-squared statistic
    LR_normality_test(i)=p;
end

% 若服从正态分布，则使用t-test进行检验，若不服从正态分布，则使用Mann–Whitney U test进行检验
for i=1:3
    if i<3
        if LR_normality_test(1)>0.05 && LR_normality_test(i+1)>0.05
            [h1,p]=ttest2(lumen_LR(1:i_number(1),1),lumen_LR(1:i_number(i+1),i+1),'Alpha',0.05,'Tail','both');
            LR_test(i)=p;
        else
            p=ranksum(lumen_LR(1:i_number(1),1),lumen_LR(1:i_number(i+1),i+1),'Alpha',0.05,'Tail','both');
            LR_test(i)=p;
        end
    else
        if LR_normality_test(2)>0.05 && LR_normality_test(3)>0.05
            [h1,p]=ttest2(lumen_LR(1:i_number(2),2),lumen_LR(1:i_number(3),3),'Alpha',0.05,'Tail','both');
            LR_test(3)=p;
        else
            p=ranksum(lumen_LR(1:i_number(2),2),lumen_LR(1:i_number(3),3),'Alpha',0.05,'Tail','both');
            LR_test(3)=p;
        end
    end
end


% lumen volume的统计检验
volume_normality_test=zeros(1,3);
for i=1:3
    % 首先通过D’Agostino–Pearson normality test判断实验样本是否服从正态分布
    X=lumen_volume(1:i_number(i),1);
    X=sort(X);
    X=X(:);
    n=length(X);

    [c,v]=hist(X,X);  % record of data in a frequency table form
    nc=find(c~=0);
    c=[v(nc) c(nc)'];

    x=c(:,1);
    f=c(:,2);
    s1=f'*x;
    s2=f'*x.^2;
    s3=f'*x.^3;
    s4=f'*x.^4;
    SS=s2-(s1^2/n);
    v=SS/(n-1);
    k3=((n*s3)-(3*s1*s2)+((2*(s1^3))/n))/((n-1)*(n-2));
    g1=k3/sqrt(v^3);
    k4=((n+1)*((n*s4)-(4*s1*s3)+(6*(s1^2)*(s2/n))-((3*(s1^4))/(n^2)))/((n-1)*(n-2)*(n-3)))-((3*(SS^2))/((n-2)*(n-3)));
    g2=k4/v^2;
    eg1=((n-2)*g1)/sqrt(n*(n-1));  % measure of skewness
    eg2=((n-2)*(n-3)*g2)/((n+1)*(n-1))+((3*(n-1))/(n+1));  % measure of kurtosis

    A=eg1*sqrt(((n+1)*(n+3))/(6*(n-2)));
    B=(3*((n^2)+(27*n)-70)*((n+1)*(n+3)))/((n-2)*(n+5)*(n+7)*(n+9));
    C=sqrt(2*(B-1))-1;
    D=sqrt(C);
    E=1/sqrt(log(D));
    F=A/sqrt(2/(C-1));
    Zg1=E*log(F+sqrt(F^2+1));

    G=(24*n*(n-2)*(n-3))/((n+1)^2*(n+3)*(n+5));
    H=((n-2)*(n-3)*g2)/((n+1)*(n-1)*sqrt(G));
    % H=((n-2)*(n-3)*abs(g2))/((n+1)*(n-1)*sqrt(G));
    J=((6*(n^2-(5*n)+2))/((n+7)*(n+9)))*sqrt((6*(n+3)*(n+5))/((n*(n-2)*(n-3))));
    K=6+((8/J)*((2/J)+sqrt(1+(4/J^2))));
    L=(1-(2/K))/(1+H*sqrt(2/(K-4)));
    Zg2=(1-(2/(9*K))-L^(1/3))/sqrt(2/(9*K));

    K2=Zg1^2+Zg2^2;  % D'Agostino-Pearson statistic
    X2=K2;  % approximation to chi-distribution
    df=2;  % degrees of freedom

    p=1-chi2cdf(X2,df);  % probability associated to the chi-squared statistic
    volume_normality_test(i)=p;
end

% 若服从正态分布，则使用t-test进行检验，若不服从正态分布，则使用Mann–Whitney U test进行检验
for i=1:3
    if i<3
        if volume_normality_test(1)>0.05 && volume_normality_test(i+1)>0.05
            [h1,p]=ttest2(lumen_volume(1:i_number(1),1),lumen_volume(1:i_number(i+1),i+1),'Alpha',0.05,'Tail','both');
            volume_test(i)=p;
        else
            p=ranksum(lumen_volume(1:i_number(1),1),lumen_volume(1:i_number(i+1),i+1),'Alpha',0.05,'Tail','both');
            volume_test(i)=p;
        end
    else
        if volume_normality_test(2)>0.05 && volume_normality_test(3)>0.05
            [h1,p]=ttest2(lumen_volume(1:i_number(2),2),lumen_volume(1:i_number(3),3),'Alpha',0.05,'Tail','both');
            volume_test(3)=p;
        else
            p=ranksum(lumen_volume(1:i_number(2),2),lumen_volume(1:i_number(3),3),'Alpha',0.05,'Tail','both');
            volume_test(3)=p;
        end
    end
end



