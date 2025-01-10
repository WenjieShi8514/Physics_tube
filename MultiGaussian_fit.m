function [rho0_best,rho_inf_array,E_array,w_array,break_point,wucha_best,TriGaussian_parameter] = MultiGaussian_fit(x,y,data_changdubi,anterior_lateral_length,basal_length,posterior_lateral_length,cell_type)
warning('off','curvefit:fit:usingTrustRegion')
warning('off','curvefit:fit:noStartPoint')

%%% Cell Grayscale Data Fitting with Multi-Gaussian Functions  
%%% Steps:  
%%% 1. Identify the boundary point between the basal domain and lateral domain as the segmentation basis, 
%%%    and extract the middle section.  
%%%   (1.1) For non-adjacent cells, first use a Tri-Gaussian fit to locate the approximate region of 
%%%         the breakpoint. Around the breakpoint, take a smaller range and use cubic spline fitting to find 
%%%         the minimum value point within the breakpoint range, which will serve as the boundary point.  
%%%   (1.2) For adjacent cells, directly determine the position of the boundary point based on the lengths
%%%         of the three data files.  
%%% 2. Fit the lateral domain section with a mono-Gaussian function, and use spline regression for the 
%%%    basal domain. Determine the number of Gaussian functions based on the extrema.  
%%% 3. Perform multi-Gaussian fittings for the basal domain using the number of Gaussians determined 
%%%    by the spline regression.  
%%% 4. Fix the mean values of each Gaussian in the basal domain's multi-Gaussian fit results and perform
%%%    a global multi-Gaussian fit for the basal-lateral domain.  

%% step 1
if cell_type==1
    % step (1.1)
    if data_changdubi<=0.2  % lumen形成早期，lateral domain收缩环未建立，共有9个特征参数拟合
        E_anterior_start=-0.5;
        E_anterior_end=-0.3;
        E_posterior_start=0.3;
        E_posterior_end=0.5;
    else  % lumen形成中后期，lateral domain收缩环建立，共有8个特征参数拟合
        E_anterior_start=-0.5;
        E_anterior_end=-0.5;
        E_posterior_start=0.5;
        E_posterior_end=0.5;
    end
    
    % 首先通过三高斯拟合找到断点的大致区域
    ft=fittype(['rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)+' ...
        'rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)'],'independent','x','dependent','y');
    opts=fitoptions('Method','NonlinearLeastSquares');
    opts.Algorithm='Levenberg-Marquardt';
    opts.Display='Off';
    opts.MaxFunEvals=600;
    opts.MaxIter=400;
    opts.Lower=[-0.1 E_anterior_start E_posterior_start 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];
    opts.Upper=[0.1 E_anterior_end E_posterior_end max(y)/2 max(y) max(y) max(y) 1 1 1];
    [xData,yData]=prepareCurveData(x,y);
    % Fit model to data.
    wucha_best=1000000000;
    for i_fit=1:10  % 循环10次，寻找最优拟合解
        [fitresult,~]=fit(xData,yData,ft,opts);
        % Fit parameter result to variables
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
            E_basal_best=E;
            rho_inf_anterior_best=rho_inf_anterior;
            w_anterior_best=w_anterior;
            rho_inf_posterior_best=rho_inf_posterior;
            w_posterior_best=w_posterior;
            E_anterior_best=E_anterior;
            E_posterior_best=E_posterior;
        end
    end
    TriGaussian_parameter=[rho0_best rho_inf_best w_best E_basal_best rho_inf_anterior_best w_anterior_best E_anterior_best rho_inf_posterior_best w_posterior_best E_posterior_best];
    
    % 寻找极小值，作为basal domain和lateral domain断点
    y_fit=rho0_best+rho_inf_best*exp(-1/2*((x-E_basal_best)/w_best).^2)+rho_inf_anterior_best*exp(-1/2*((x-E_anterior_best)/w_anterior_best).^2)+...
        rho_inf_posterior_best*exp(-1/2*((x-E_posterior_best)/w_posterior_best).^2);
    break_point=islocalmin(y_fit);
    break_point=x(break_point);
    if length(break_point)<=1
        [~,index]=min(abs(x-anterior_lateral_length/(anterior_lateral_length+basal_length+posterior_lateral_length)+0.5));
        break_point(1)=max(x(index),x(5));
        [~,index]=min(abs(x-(anterior_lateral_length+basal_length)/(anterior_lateral_length+basal_length+posterior_lateral_length)+0.5));
        break_point(2)=min(x(index),x(end-5));
    end
    if length(x(1:find(x==break_point(1))))<=5
        break_point(1)=x(5);
    end
    if length(x(find(x==break_point(1))+1:find(x==break_point(2))))<=5
        break_point(2)=x(end-5);
    end

    % 在断点处左右取较小范围
    break_point_range=zeros(2,2);  % 两个断点范围的左右边界
    break_point_range(1,1)=max(break_point(1)-0.1,-0.5);break_point_range(1,2)=min(break_point(1)+0.1,0.5);
    break_point_range(2,1)=max(break_point(2)-0.1,-0.5);break_point_range(2,2)=min(break_point(2)+0.1,0.5);

    % 通过三次样条拟合找出断点范围内的最小值点，作为交界点
    k=50;
    b=linspace(min(x),max(x),k+1);  % cubic spline断点
    spline1=spline(b,y'/spline(b,eye(length(b)),x'));  % cubic spline拟合
    y_cubicspline=fnval(spline1,x);

    local_min_x=islocalmin(y_cubicspline);
    local_min_y=y(local_min_x);  % 所有极小值点对应的y坐标
    local_min_x=x(local_min_x);  % 所有极小值点对应的x坐标
    y_minimum=[1000 1000];
    for ii=1:length(local_min_x)  % 找出断点范围内的最小值点
        if local_min_x(ii)>=break_point_range(1,1) && local_min_x(ii)<=break_point_range(1,2)
            if local_min_y(ii)<y_minimum(1)
                y_minimum(1)=local_min_y(ii);
                break_point(1)=local_min_x(ii);
            end
        end

        if local_min_x(ii)>=break_point_range(2,1) && local_min_x(ii)<=break_point_range(2,2)
            if local_min_y(ii)<y_minimum(2)
                y_minimum(2)=local_min_y(ii);
                break_point(2)=local_min_x(ii);
            end
        end
    end
    
    if break_point(1)<x(5)
        break_point(1)=x(5);
    end
    if break_point(2)>x(end-5)
        break_point(2)=x(end-5);
    end
    if length(x(find(x==break_point(1)):find(x==break_point(2))))<=15
        if length(x(find(x==break_point(2)):end))>20
            break_point(2)=x(find(x==break_point(2))+15);
        elseif length(x(1:find(x==break_point(1))))>20
            break_point(1)=x(find(x==break_point(1))-15);
        end
    end
elseif cell_type==2
    % step (1.2)
    [~,index]=min(abs(x-anterior_lateral_length/(anterior_lateral_length+basal_length+posterior_lateral_length)+0.5));
    break_point(1)=x(index);
    [~,index]=min(abs(x-(anterior_lateral_length+basal_length)/(anterior_lateral_length+basal_length+posterior_lateral_length)+0.5));
    break_point(2)=x(index);
    
    if break_point(1)<x(5)
        break_point(1)=x(5);
    end
    if break_point(2)>x(end-5)
        break_point(2)=x(end-5);
    end
    if length(x(find(x==break_point(1)):find(x==break_point(2))))<=15
        if length(x(find(x==break_point(2)):end))>20
            break_point(2)=x(find(x==break_point(2))+15);
        elseif length(x(1:find(x==break_point(1))))>20
            break_point(1)=x(find(x==break_point(1))-15);
        end
    end

    % 拟合直接三高斯分布特征参数
    if data_changdubi<=0.2  % lumen形成早期，lateral domain收缩环未建立，共有9个特征参数拟合
        E_anterior_start=-0.5;
        E_anterior_end=break_point(1);
        E_posterior_start=break_point(2);
        E_posterior_end=0.5;
    else  % lumen形成中后期，lateral domain收缩环建立，共有8个特征参数拟合
        E_anterior_start=-0.5;
        E_anterior_end=-0.5;
        E_posterior_start=0.5;
        E_posterior_end=0.5;
    end

    ft=fittype(['rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)+' ...
        'rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)'],'independent','x','dependent','y');
    opts=fitoptions('Method','NonlinearLeastSquares');
    opts.Algorithm='Levenberg-Marquardt';
    opts.Display='Off';
    opts.MaxFunEvals=600;
    opts.MaxIter=400;
    opts.Lower=[-0.1 E_anterior_start E_posterior_start 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];
    opts.Upper=[0.1 E_anterior_end E_posterior_end max(y)/2 max(y) max(y) max(y) 1 1 1];
    [xData,yData]=prepareCurveData(x,y);
    % Fit model to data.
    wucha_best=1000000000;
    for i_fit=1:10  % 循环10次，寻找最优拟合解
        [fitresult,~]=fit(xData,yData,ft,opts);
        % Fit parameter result to variables
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
            E_basal_best=E;
            rho_inf_anterior_best=rho_inf_anterior;
            w_anterior_best=w_anterior;
            rho_inf_posterior_best=rho_inf_posterior;
            w_posterior_best=w_posterior;
            E_anterior_best=E_anterior;
            E_posterior_best=E_posterior;
        end
    end
    TriGaussian_parameter=[rho0_best rho_inf_best w_best E_basal_best rho_inf_anterior_best w_anterior_best E_anterior_best rho_inf_posterior_best w_posterior_best E_posterior_best];
end

%% step 2
x_anterior_lateral=x(1:find(x==break_point(1)));  % anterior-lateral domain区域
y_anterior_lateral=y(1:find(x==break_point(1)));
x_basal=x(find(x==break_point(1))+1:find(x==break_point(2)));  % basal domain区域
y_basal=y(find(x==break_point(1))+1:find(x==break_point(2)));
x_posterior_lateral=x(find(x==break_point(2))+1:end);  % posterior-lateral domain区域
y_posterior_lateral=y(find(x==break_point(2))+1:end);

% lateral domain单高斯拟合
for ab=1
    ft_anterior_lateral=fittype('rho0+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)','independent','x','dependent','y');
    opts=fitoptions('Method','NonlinearLeastSquares');
    opts.Algorithm='Levenberg-Marquardt';
    opts.Display='Off';
    opts.MaxFunEvals=600;
    opts.MaxIter=400;
    opts.Lower=[E_anterior_start 0.0001 0.0001 0.0001];
    opts.Upper=[E_anterior_end max(y_anterior_lateral)/2 max(y_anterior_lateral) 1];
    [xData_anterior_lateral,yData_anterior_lateral]=prepareCurveData(x_anterior_lateral,y_anterior_lateral);
    % Fit model to data.
    wucha_best1=100000000;
    for i_fit=1:10  % 循环10次，寻找最优拟合解
        [fitresult_anterior_lateral,~]=fit(xData_anterior_lateral,yData_anterior_lateral,ft_anterior_lateral,opts);
        % Fit parameter result to variables
        rho0_anterior=fitresult_anterior_lateral.rho0;
        rho_inf_anterior=fitresult_anterior_lateral.rho_inf_anterior;
        w_anterior=fitresult_anterior_lateral.w_anterior;
        E_anterior=fitresult_anterior_lateral.E_anterior;

        wucha1=0;
        for ii=1:length(x_anterior_lateral)
            wucha1=wucha1+(rho0_anterior+rho_inf_anterior*exp(-1/2*((x_anterior_lateral(ii)-E_anterior)/w_anterior)^2)-y_anterior_lateral(ii))^2;
        end
        if wucha1<wucha_best1
            wucha_best1=wucha1;
            rho0_anterior_best=rho0_anterior;
            rho_inf_anterior_best=rho_inf_anterior;
            w_anterior_best=w_anterior;
            E_anterior_best=E_anterior;
        end
    end
    ft_posterior_lateral=fittype('rho0+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)','independent','x','dependent','y');
    opts=fitoptions('Method','NonlinearLeastSquares');
    opts.Algorithm='Levenberg-Marquardt';
    opts.Display='Off';
    opts.MaxFunEvals=600;
    opts.MaxIter=400;
    opts.Lower=[E_posterior_start 0.0001 0.0001 0.0001];
    opts.Upper=[E_posterior_end max(y_posterior_lateral)/2 max(y_posterior_lateral) 1];
    [xData_posterior_lateral,yData_posterior_lateral]=prepareCurveData(x_posterior_lateral,y_posterior_lateral);
    % Fit model to data.
    wucha_best2=100000000;
    for i_fit=1:10  % 循环10次，寻找最优拟合解
        [fitresult_posterior_lateral,~]=fit(xData_posterior_lateral,yData_posterior_lateral,ft_posterior_lateral,opts);
        % Fit parameter result to variables
        rho0_posterior=fitresult_posterior_lateral.rho0;
        rho_inf_posterior=fitresult_posterior_lateral.rho_inf_posterior;
        w_posterior=fitresult_posterior_lateral.w_posterior;
        E_posterior=fitresult_posterior_lateral.E_posterior;

        wucha2=0;
        for ii=1:length(x_posterior_lateral)
            wucha2=wucha2+(rho0_posterior+rho_inf_posterior*exp(-1/2*((x_posterior_lateral(ii)-E_posterior)/w_posterior)^2)-y_posterior_lateral(ii))^2;
        end
        if wucha2<wucha_best2
            wucha_best2=wucha2;
            rho0_posterior_best=rho0_posterior;
            rho_inf_posterior_best=rho_inf_posterior;
            w_posterior_best=w_posterior;
            E_posterior_best=E_posterior;
        end
    end
end

% basal domain的部分使用样条回归
k=round((length(x_basal))^(1/3)*1.6)+1;
spline1=spap2(k,4,x_basal,y_basal);  % 4 order B-spline拟合
% b=linspace(min(x_basal),1,k+1);  % cubic spline断点
% spline1=spline(b,y_basal'/spline(b,eye(length(b)),x_basal'));  % cubic spline拟合
y_basal_bspline=fnval(spline1,x_basal);

break_point_basal=islocalmin(y_basal_bspline);
break_point_basal=x_basal(break_point_basal);

Gaussian_number=length(break_point_basal)+1;  % basal domain的高斯个数


%% step 3
% 根据高斯个数确定拟合函数形式
opts=fitoptions('Method','NonlinearLeastSquares');
opts.Algorithm='Levenberg-Marquardt';
opts.Display='Off';
opts.MaxFunEvals=600;
opts.MaxIter=400;
if Gaussian_number==1
    ft=fittype('rho0+rho_inf1*exp(-1/2*((x-E1)/w1).^2)','independent','x','dependent','y');
    opts.Lower=[min(x_basal) 0.0001 0.0001 0.0001];
    opts.Upper=[max(x_basal) max(y_basal)/2 max(y_basal) 1];
elseif Gaussian_number==2
    ft=fittype('rho0+rho_inf1*exp(-1/2*((x-E1)/w1).^2)+rho_inf2*exp(-1/2*((x-E2)/w2).^2)','independent','x','dependent','y');
    opts.Lower=[min(x_basal) min(x_basal) 0.0001 0.0001 0.0001 0.0001 0.0001];
    opts.Upper=[max(x_basal) max(x_basal) max(y_basal)/2 max(y_basal) max(y_basal) 1 1];
elseif Gaussian_number==3
    ft=fittype(['rho0+rho_inf1*exp(-1/2*((x-E1)/w1).^2)+rho_inf2*exp(-1/2*((x-E2)/w2).^2)' ...
        '+rho_inf3*exp(-1/2*((x-E3)/w3).^2)'],'independent','x','dependent','y');
    opts.Lower=[min(x_basal) min(x_basal) min(x_basal) 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];
    opts.Upper=[max(x_basal) max(x_basal) max(x_basal) max(y_basal)/2 max(y_basal) max(y_basal) max(y_basal) 1 1 1];
elseif Gaussian_number==4
    ft=fittype(['rho0+rho_inf1*exp(-1/2*((x-E1)/w1).^2)+rho_inf2*exp(-1/2*((x-E2)/w2).^2)' ...
        '+rho_inf3*exp(-1/2*((x-E3)/w3).^2)+rho_inf4*exp(-1/2*((x-E4)/w4).^2)'],'independent','x','dependent','y');
    opts.Lower=[min(x_basal) min(x_basal) min(x_basal) min(x_basal) 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];
    opts.Upper=[max(x_basal) max(x_basal) max(x_basal) max(x_basal) max(y_basal)/2 max(y_basal) max(y_basal) max(y_basal) max(y_basal) 1 1 1 1];
elseif Gaussian_number>=5
    ft=fittype(['rho0+rho_inf1*exp(-1/2*((x-E1)/w1).^2)+rho_inf2*exp(-1/2*((x-E2)/w2).^2)' ...
        '+rho_inf3*exp(-1/2*((x-E3)/w3).^2)+rho_inf4*exp(-1/2*((x-E4)/w4).^2)+rho_inf5*exp(-1/2*((x-E5)/w5).^2)'],'independent','x','dependent','y');
    opts.Lower=[min(x_basal) min(x_basal) min(x_basal) min(x_basal) min(x_basal) 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];
    opts.Upper=[max(x_basal) max(x_basal) max(x_basal) max(x_basal) max(x_basal) max(y_basal)/2 max(y_basal) max(y_basal) max(y_basal) max(y_basal) max(y_basal) 1 1 1 1 1];
end
[xData,yData]=prepareCurveData(x_basal,y_basal);
% Fit model to data.
wucha_best=1000000000;
for i_fit=1:20  % 循环50次，寻找最优拟合解
    wucha=0;
    [fitresult,~]=fit(xData,yData,ft,opts);
    % Fit parameter result to variables
    rho0=fitresult.rho0;
    if Gaussian_number==1
        rho_inf1=fitresult.rho_inf1;
        w1=fitresult.w1;
        E1=fitresult.E1;
        for ii=1:length(x_basal)
            wucha=wucha+(rho0+rho_inf1*exp(-1/2*((x_basal(ii)-E1)/w1).^2)-y_basal(ii))^2;
        end
    elseif Gaussian_number==2
        rho_inf1=fitresult.rho_inf1;rho_inf2=fitresult.rho_inf2;
        w1=fitresult.w1;w2=fitresult.w2;
        E1=fitresult.E1;E2=fitresult.E2;
        for ii=1:length(x_basal)
            wucha=wucha+(rho0+rho_inf1*exp(-1/2*((x_basal(ii)-E1)/w1).^2)+rho_inf2*exp(-1/2*((x_basal(ii)-E2)/w2).^2)-y_basal(ii))^2;
        end
    elseif Gaussian_number==3
        rho_inf1=fitresult.rho_inf1;rho_inf2=fitresult.rho_inf2;rho_inf3=fitresult.rho_inf3;
        w1=fitresult.w1;w2=fitresult.w2;w3=fitresult.w3;
        E1=fitresult.E1;E2=fitresult.E2;E3=fitresult.E3;
        for ii=1:length(x_basal)
            wucha=wucha+(rho0+rho_inf1*exp(-1/2*((x_basal(ii)-E1)/w1).^2)+rho_inf2*exp(-1/2*((x_basal(ii)-E2)/w2).^2)...
                +rho_inf3*exp(-1/2*((x_basal(ii)-E3)/w3).^2)-y_basal(ii))^2;
        end
    elseif Gaussian_number==4
        rho_inf1=fitresult.rho_inf1;rho_inf2=fitresult.rho_inf2;rho_inf3=fitresult.rho_inf3;rho_inf4=fitresult.rho_inf4;
        w1=fitresult.w1;w2=fitresult.w2;w3=fitresult.w3;w4=fitresult.w4;
        E1=fitresult.E1;E2=fitresult.E2;E3=fitresult.E3;E4=fitresult.E4;
        for ii=1:length(x_basal)
            wucha=wucha+(rho0+rho_inf1*exp(-1/2*((x_basal(ii)-E1)/w1).^2)+rho_inf2*exp(-1/2*((x_basal(ii)-E2)/w2).^2)...
                +rho_inf3*exp(-1/2*((x_basal(ii)-E3)/w3).^2)+rho_inf4*exp(-1/2*((x_basal(ii)-E4)/w4).^2)-y_basal(ii))^2;
        end
    elseif Gaussian_number>=5
        rho_inf1=fitresult.rho_inf1;rho_inf2=fitresult.rho_inf2;rho_inf3=fitresult.rho_inf3;rho_inf4=fitresult.rho_inf4;rho_inf5=fitresult.rho_inf5;
        w1=fitresult.w1;w2=fitresult.w2;w3=fitresult.w3;w4=fitresult.w4;w5=fitresult.w5;
        E1=fitresult.E1;E2=fitresult.E2;E3=fitresult.E3;E4=fitresult.E4;E5=fitresult.E5;
        for ii=1:length(x_basal)
            wucha=wucha+(rho0+rho_inf1*exp(-1/2*((x_basal(ii)-E1)/w1).^2)+rho_inf2*exp(-1/2*((x_basal(ii)-E2)/w2).^2)...
                +rho_inf3*exp(-1/2*((x_basal(ii)-E3)/w3).^2)+rho_inf4*exp(-1/2*((x_basal(ii)-E4)/w4).^2) ...
                +rho_inf5*exp(-1/2*((x_basal(ii)-E5)/w5).^2)-y_basal(ii))^2;
        end
    end

    if wucha<wucha_best
        wucha_best=wucha;
        rho0_best=rho0;
        if Gaussian_number==1
            rho_inf_best(1)=fitresult.rho_inf1;
            w_best(1)=fitresult.w1;
            E_basal_best(1)=fitresult.E1;
        elseif Gaussian_number==2
            rho_inf_best(1)=fitresult.rho_inf1;rho_inf_best(2)=fitresult.rho_inf2;
            w_best(1)=fitresult.w1;w_best(2)=fitresult.w2;
            E_basal_best(1)=fitresult.E1;E_basal_best(2)=fitresult.E2;
        elseif Gaussian_number==3
            rho_inf_best(1)=fitresult.rho_inf1;rho_inf_best(2)=fitresult.rho_inf2;rho_inf_best(3)=fitresult.rho_inf3;
            w_best(1)=fitresult.w1;w_best(2)=fitresult.w2;w_best(3)=fitresult.w3;
            E_basal_best(1)=fitresult.E1;E_basal_best(2)=fitresult.E2;E_basal_best(3)=fitresult.E3;
        elseif Gaussian_number==4
            rho_inf_best(1)=fitresult.rho_inf1;rho_inf_best(2)=fitresult.rho_inf2;rho_inf_best(3)=fitresult.rho_inf3;rho_inf_best(4)=fitresult.rho_inf4;
            w_best(1)=fitresult.w1;w_best(2)=fitresult.w2;w_best(3)=fitresult.w3;w_best(4)=fitresult.w4;
            E_basal_best(1)=fitresult.E1;E_basal_best(2)=fitresult.E2;E_basal_best(3)=fitresult.E3;E_basal_best(4)=fitresult.E4;
        elseif Gaussian_number>=5
            rho_inf_best(1)=fitresult.rho_inf1;rho_inf_best(2)=fitresult.rho_inf2;rho_inf_best(3)=fitresult.rho_inf3;rho_inf_best(4)=fitresult.rho_inf4;rho_inf_best(5)=fitresult.rho_inf5;
            w_best(1)=fitresult.w1;w_best(2)=fitresult.w2;w_best(3)=fitresult.w3;w_best(4)=fitresult.w4;w_best(5)=fitresult.w5;
            E_basal_best(1)=fitresult.E1;E_basal_best(2)=fitresult.E2;E_basal_best(3)=fitresult.E3;E_basal_best(4)=fitresult.E4;E_basal_best(5)=fitresult.E5;
        end
    end
end

% 最佳basal domain Multi-Gaussian拟合结果赋值给y_basal_fit
% y_basal_fit=ones(size(y_basal))*rho0_best;
% for i_Gaussian=1:Gaussian_number
%     y_basal_fit=y_basal_fit+rho_inf_best(i_Gaussian)*exp(-1/2*((x_basal-E_basal_best(i_Gaussian))/w_best(i_Gaussian)).^2);
% end


%% step 4
[E_basal_best,I]=sort(E_basal_best);
E_best=[E_basal_best E_anterior_best E_posterior_best];
rho_inf_best=[rho_inf_best(I) rho_inf_anterior_best rho_inf_posterior_best];

min_y=max(0,min(y));

opts=fitoptions('Method','NonlinearLeastSquares');
opts.Algorithm='Levenberg-Marquardt';
opts.Display='Off';
opts.MaxFunEvals=600;
opts.MaxIter=400;
if Gaussian_number==1
    ft=fittype(['rho0+rho_inf1*exp(-1/2*((x-E1)/w1).^2)+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)' ...
        '+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)'],'independent','x','dependent','y');
    opts.Lower=[E_best(1) E_best(2) E_best(3) 0 rho_inf_best(1)*0.7 rho_inf_best(2)*0.7 rho_inf_best(3)*0.7 0.0001 0.0001 0.0001];
    opts.Upper=[E_best(1) E_best(2) E_best(3) min_y rho_inf_best(1)*1.5 rho_inf_best(2)*1.5 rho_inf_best(3)*1.5 1 1 1];
elseif Gaussian_number==2
    ft=fittype(['rho0+rho_inf1*exp(-1/2*((x-E1)/w1).^2)+rho_inf2*exp(-1/2*((x-E2)/w2).^2)' ...
        '+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)'],'independent','x','dependent','y');
    opts.Lower=[E_best(1) E_best(2) E_best(3) E_best(4)...
        0 rho_inf_best(1)*0.7 rho_inf_best(2)*0.7 rho_inf_best(3)*0.7 rho_inf_best(4)*0.7 0.0001 0.0001 0.0001 0.0001];
    opts.Upper=[E_best(1) E_best(2) E_best(3) E_best(4)...
        min_y rho_inf_best(1)*1.5 rho_inf_best(2)*1.5 rho_inf_best(3)*1.5 rho_inf_best(4)*1.5 1 1 1 1];
elseif Gaussian_number==3
    ft=fittype(['rho0+rho_inf1*exp(-1/2*((x-E1)/w1).^2)+rho_inf2*exp(-1/2*((x-E2)/w2).^2)+rho_inf3*exp(-1/2*((x-E3)/w3).^2)' ...
        '+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)' ...
        '+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)'],'independent','x','dependent','y');
    opts.Lower=[E_best(1) E_best(2) E_best(3) E_best(4) E_best(5)...
        0 rho_inf_best(1)*0.7 rho_inf_best(2)*0.7 rho_inf_best(3)*0.7 rho_inf_best(4)*0.7 rho_inf_best(5)*0.7 0.0001 0.0001 0.0001 0.0001 0.0001];
    opts.Upper=[E_best(1) E_best(2) E_best(3) E_best(4) E_best(5)...
        min_y rho_inf_best(1)*1.5 rho_inf_best(2)*1.5 rho_inf_best(3)*1.5 rho_inf_best(4)*1.5 rho_inf_best(5)*1.5 1 1 1 1 1];
elseif Gaussian_number==4
    ft=fittype(['rho0+rho_inf1*exp(-1/2*((x-E1)/w1).^2)+rho_inf2*exp(-1/2*((x-E2)/w2).^2)+rho_inf3*exp(-1/2*((x-E3)/w3).^2)' ...
        '+rho_inf4*exp(-1/2*((x-E4)/w4).^2)+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)' ...
        '+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)'],'independent','x','dependent','y');
    opts.Lower=[E_best(1) E_best(2) E_best(3) E_best(4) E_best(5) E_best(6)...
        0 rho_inf_best(1)*0.7 rho_inf_best(2)*0.7 rho_inf_best(3)*0.7 rho_inf_best(4)*0.7 rho_inf_best(5)*0.7 rho_inf_best(6)*0.7 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];
    opts.Upper=[E_best(1) E_best(2) E_best(3) E_best(4) E_best(5) E_best(6)...
        min_y rho_inf_best(1)*1.5 rho_inf_best(2)*1.5 rho_inf_best(3)*1.5 rho_inf_best(4)*1.5 rho_inf_best(5)*1.5 rho_inf_best(6)*1.5 1 1 1 1 1 1];
elseif Gaussian_number>=5
    ft=fittype(['rho0+rho_inf1*exp(-1/2*((x-E1)/w1).^2)+rho_inf2*exp(-1/2*((x-E2)/w2).^2)+rho_inf3*exp(-1/2*((x-E3)/w3).^2)' ...
        '+rho_inf4*exp(-1/2*((x-E4)/w4).^2)+rho_inf5*exp(-1/2*((x-E5)/w5).^2)+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)' ...
        '+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)'],'independent','x','dependent','y');
    opts.Lower=[E_best(1) E_best(2) E_best(3) E_best(4) E_best(5) E_best(6) E_best(7)...
        0 rho_inf_best(1)*0.7 rho_inf_best(2)*0.7 rho_inf_best(3)*0.7 rho_inf_best(4)*0.7 rho_inf_best(5)*0.7 rho_inf_best(6)*0.7 rho_inf_best(7)*0.7 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];
    opts.Upper=[E_best(1) E_best(2) E_best(3) E_best(4) E_best(5) E_best(6) E_best(7)...
        min_y rho_inf_best(1)*1.5 rho_inf_best(2)*1.5 rho_inf_best(3)*1.5 rho_inf_best(4)*1.5 rho_inf_best(5)*1.5 rho_inf_best(6)*1.5 rho_inf_best(7)*1.5 1 1 1 1 1 1 1];
end
[xData,yData]=prepareCurveData(x,y);
% Fit model to data.
wucha_best=1000000000;
fit_logic=1;  % 判断是否需要重复拟合的逻辑值
fit_cycle_number=1;  % 重复拟合的循环次数
while fit_logic==1  % 逻辑值为1时，重复拟合
    for i_fit=1:20  % 循环100次，寻找最优拟合解
        wucha=0;
        [fitresult,~]=fit(xData,yData,ft,opts);
        % Fit parameter result to variables
        rho0=fitresult.rho0;
        rho_inf_anterior=fitresult.rho_inf_anterior;
        rho_inf_posterior=fitresult.rho_inf_posterior;
        w_anterior=fitresult.w_anterior;
        w_posterior=fitresult.w_posterior;
        for ii=1:length(x)
            wucha=wucha+(rho_inf_anterior*exp(-1/2*((x(ii)-E_anterior)/w_anterior)^2)+rho_inf_posterior*exp(-1/2*((x(ii)-E_posterior)/w_posterior)^2)-y(ii))^2;
        end
        if Gaussian_number==1
            rho_inf1=fitresult.rho_inf1;
            w1=fitresult.w1;
            E1=fitresult.E1;
            for ii=1:length(x)
                wucha=wucha+(rho0+rho_inf1*exp(-1/2*((x(ii)-E1)/w1).^2)-y(ii))^2;
            end
        elseif Gaussian_number==2
            rho_inf1=fitresult.rho_inf1;rho_inf2=fitresult.rho_inf2;
            w1=fitresult.w1;w2=fitresult.w2;
            E1=fitresult.E1;E2=fitresult.E2;
            for ii=1:length(x)
                wucha=wucha+(rho0+rho_inf1*exp(-1/2*((x(ii)-E1)/w1).^2)+rho_inf2*exp(-1/2*((x(ii)-E2)/w2).^2)-y(ii))^2;
            end
        elseif Gaussian_number==3
            rho_inf1=fitresult.rho_inf1;rho_inf2=fitresult.rho_inf2;rho_inf3=fitresult.rho_inf3;
            w1=fitresult.w1;w2=fitresult.w2;w3=fitresult.w3;
            E1=fitresult.E1;E2=fitresult.E2;E3=fitresult.E3;
            for ii=1:length(x)
                wucha=wucha+(rho0+rho_inf1*exp(-1/2*((x(ii)-E1)/w1).^2)+rho_inf2*exp(-1/2*((x(ii)-E2)/w2).^2)...
                    +rho_inf3*exp(-1/2*((x(ii)-E3)/w3).^2)-y(ii))^2;
            end
        elseif Gaussian_number==4
            rho_inf1=fitresult.rho_inf1;rho_inf2=fitresult.rho_inf2;rho_inf3=fitresult.rho_inf3;rho_inf4=fitresult.rho_inf4;
            w1=fitresult.w1;w2=fitresult.w2;w3=fitresult.w3;w4=fitresult.w4;
            E1=fitresult.E1;E2=fitresult.E2;E3=fitresult.E3;E4=fitresult.E4;
            for ii=1:length(x)
                wucha=wucha+(rho0+rho_inf1*exp(-1/2*((x(ii)-E1)/w1).^2)+rho_inf2*exp(-1/2*((x(ii)-E2)/w2).^2)...
                    +rho_inf3*exp(-1/2*((x(ii)-E3)/w3).^2)+rho_inf4*exp(-1/2*((x(ii)-E4)/w4).^2)-y(ii))^2;
            end
        elseif Gaussian_number>=5
            rho_inf1=fitresult.rho_inf1;rho_inf2=fitresult.rho_inf2;rho_inf3=fitresult.rho_inf3;rho_inf4=fitresult.rho_inf4;rho_inf5=fitresult.rho_inf5;
            w1=fitresult.w1;w2=fitresult.w2;w3=fitresult.w3;w4=fitresult.w4;w5=fitresult.w5;
            E1=fitresult.E1;E2=fitresult.E2;E3=fitresult.E3;E4=fitresult.E4;E5=fitresult.E5;
            for ii=1:length(x)
                wucha=wucha+(rho0+rho_inf1*exp(-1/2*((x(ii)-E1)/w1).^2)+rho_inf2*exp(-1/2*((x(ii)-E2)/w2).^2)...
                    +rho_inf3*exp(-1/2*((x(ii)-E3)/w3).^2)+rho_inf4*exp(-1/2*((x(ii)-E4)/w4).^2) ...
                    +rho_inf5*exp(-1/2*((x(ii)-E5)/w5).^2)-y(ii))^2;
            end
        end

        if wucha<wucha_best
            wucha_best=wucha;
            rho0_best=rho0;
            rho_inf_anterior_best=rho_inf_anterior;
            w_anterior_best=w_anterior;
            rho_inf_posterior_best=rho_inf_posterior;
            w_posterior_best=w_posterior;
            if Gaussian_number==1
                rho_inf_basal_best(1)=fitresult.rho_inf1;
                w_best(1)=fitresult.w1;
                E_basal_best(1)=fitresult.E1;
            elseif Gaussian_number==2
                rho_inf_basal_best(1)=fitresult.rho_inf1;rho_inf_basal_best(2)=fitresult.rho_inf2;
                w_best(1)=fitresult.w1;w_best(2)=fitresult.w2;
                E_basal_best(1)=fitresult.E1;E_basal_best(2)=fitresult.E2;
            elseif Gaussian_number==3
                rho_inf_basal_best(1)=fitresult.rho_inf1;rho_inf_basal_best(2)=fitresult.rho_inf2;rho_inf_basal_best(3)=fitresult.rho_inf3;
                w_best(1)=fitresult.w1;w_best(2)=fitresult.w2;w_best(3)=fitresult.w3;
                E_basal_best(1)=fitresult.E1;E_basal_best(2)=fitresult.E2;E_basal_best(3)=fitresult.E3;
            elseif Gaussian_number==4
                rho_inf_basal_best(1)=fitresult.rho_inf1;rho_inf_basal_best(2)=fitresult.rho_inf2;rho_inf_basal_best(3)=fitresult.rho_inf3;rho_inf_basal_best(4)=fitresult.rho_inf4;
                w_best(1)=fitresult.w1;w_best(2)=fitresult.w2;w_best(3)=fitresult.w3;w_best(4)=fitresult.w4;
                E_basal_best(1)=fitresult.E1;E_basal_best(2)=fitresult.E2;E_basal_best(3)=fitresult.E3;E_basal_best(4)=fitresult.E4;
            elseif Gaussian_number>=5
                rho_inf_basal_best(1)=fitresult.rho_inf1;rho_inf_basal_best(2)=fitresult.rho_inf2;rho_inf_basal_best(3)=fitresult.rho_inf3;rho_inf_basal_best(4)=fitresult.rho_inf4;rho_inf_basal_best(5)=fitresult.rho_inf5;
                w_best(1)=fitresult.w1;w_best(2)=fitresult.w2;w_best(3)=fitresult.w3;w_best(4)=fitresult.w4;w_best(5)=fitresult.w5;
                E_basal_best(1)=fitresult.E1;E_basal_best(2)=fitresult.E2;E_basal_best(3)=fitresult.E3;E_basal_best(4)=fitresult.E4;E_basal_best(5)=fitresult.E5;
            end
        end
    end

    % 若拟合后的E存在小于0.02的情况，说明掉入局部最优解，重复循环
    for ii=1:length(E_basal_best)
        if w_best(ii)<=0.02
            fit_logic=0;
            break
        end
    end
    if w_anterior_best<=0.02 || w_posterior_best<=0.02
        fit_logic=0;
    end
    if fit_logic==0  % 前面判断E存在小于0.02的情况后，调转fit_logic数值，重复while循环
        fit_logic=1;
    elseif fit_logic==1  % 前面判断E不存在小于0.02的情况，调转fit_logic数值，跳出while循环
        fit_logic=0;
    end
    fit_cycle_number=fit_cycle_number+1;
    if fit_cycle_number>=5  % 若while循环重复超过5次，则跳出循环
        fit_logic=0;
    end
end

% 拟合参数汇总
rho_inf_array=[rho_inf_basal_best rho_inf_anterior_best rho_inf_posterior_best];
E_array=[E_basal_best E_anterior_best E_posterior_best];
w_array=[w_best w_anterior_best w_posterior_best];

% % 最佳basal_lateral domain Multi-Gaussian拟合结果赋值给y_fit
% y_fit=rho0_best+rho_inf_anterior_best*exp(-1/2*((x-E_anterior_best)/w_anterior_best).^2)+rho_inf_posterior_best*exp(-1/2*((x-E_posterior_best)/w_posterior_best).^2);
% for i_Gaussian=1:Gaussian_number
%     y_fit=y_fit+rho_inf_basal_best(i_Gaussian)*exp(-1/2*((x-E_basal_best(i_Gaussian))/w_best(i_Gaussian)).^2);
% end
