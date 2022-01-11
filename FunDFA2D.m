%function X=FunDFA2D(I,div,SHOW)
function X=FunDFA2D(I,SHOW)
[U,V,W]=size(I);
% temp=min([U,V])/4;
% s=6:floor((temp-6)/div):temp;
% stop
% s=[6 128];
%s=[6 8 10 13 17 22 29 38 49 64];
s=[6 8 10 16 32 64];
% s=[8 16 32 64];
% s=[8 64];
for rgb=1:W
    F2=zeros(1,length(s));
    for iter=1:length(s)
        Us=floor(U/s(iter));
        Vs=floor(V/s(iter));
        % fig=figure;
        % DefaultScreenSize;
        Iuv=zeros(s(iter),s(iter),Us*Vs);
        temp=1;
        for u=1:Us
            for v=1:Vs
                Iuv(:,:,temp)=I(s(iter)*(u-1)+1:s(iter)*u,s(iter)*(v-1)+1:s(iter)*v,rgb);
                temp=temp+1;
            end
        end
        Puv=zeros(s(iter),s(iter),Us*Vs);
        for k=1:Us*Vs
            It=Iuv(:,:,k);
            It=It';
            It=It(:)';
            MIuv=mean(It);
            temp=cumsum(It-MIuv);
            Puv(:,:,k)=vec2mat(temp,s(iter));
        end
        
        [x,y]=meshgrid(1:s(iter),1:s(iter));
        x=x(:);
        y=y(:);
        gPuv=zeros(s(iter),s(iter),Us*Vs);
        for k=1:Us*Vs
            z=Puv(:,:,k);
            z=z';
            z=z(:);
            a = [x,y,ones(s(iter)*s(iter),1)] \ z;
            for i=1:s(iter)
                for j=1:s(iter)
                    gPuv(i,j,k)=a(1)*i+a(2)*j+a(3);
                end
            end
        end
        
        temp=(Puv-gPuv).^2;
        F1=sum(temp,3)/(s(iter)*s(iter));
        
        F2(iter)=sqrt((sum(sum(F1))/(Us*Vs)));
    end
    
    x=log(s);
    y=log(F2);
    
    a=(x*y'-((sum(x)*sum(y))/length(x)))/(x*x'-((sum(x)^2)/length(x)));
%     temp=(length(x)*(x*y')-sum(x)*sum(y))/(length(x)*(x*x')-sum(x)*sum(x))
    b=mean(y)-a*mean(x);
%     temp=(sum(y)*(x*x')-sum(x)*(x*y'))/(length(x)*(x*x')-sum(x)*sum(x))
    
    if SHOW==1
        %     FP=[20 60 1880 920];
        %     FS=20;
        fig=figure;
%         DefaultScreenSize;
        FS2=14;
        FN='Times New Roman';
        subplot(1,2,1);
        imshow(I,[]);
        sp=subplot(1,2,2);
        plot1=plot(x,a*x+b);
        hold on;
        plot2=plot(x,y,'r');
        plot1.Color=[0 0 1];
        plot1.LineWidth=2;
        plot2.Marker='o';
        plot2.MarkerFaceColor=[1 0 0];
        plot2.LineStyle='none';
        
        sp.FontName=FN;
        sp.FontSize=FS2;
        sp.XLabel.String='ln(s)';
        sp.YLabel.String='ln(F_2(s))';
        %     sp.XLim=[1.5 4.5];
        %     sp.YLim=[4 10];
        lgd=legend(['\alpha=',num2str(a,5)],'Location','southeast');
        
        %     sp.XLim=[1.5 4.3];
        %     sp.YLim=[4 9.1];
    end
    
    X.s=s;
    X.F2(rgb,:)=F2;
    X.a(rgb,:)=a;
    X.b(rgb,:)=b;
end


if SHOW==2 && W==3
    %     FP=[20 60 1880 920];
    %     FS=20;
    fig=figure;
%     DefaultScreenSize;
    FS2=14;
    FN='Times New Roman';
    
    subplot(1,2,1);
    imshow(uint8(I),[]);
    
    rgb=1;
    sp=subplot(3,2,2);
    plot1=plot(log(X.s),X.a(rgb,:)*log(X.s)+X.b(rgb,:));
    hold on;
    plot2=plot(log(X.s),log(X.F2(rgb,:)),'r');
    plot1.Color=[0 0 1];
    plot1.LineWidth=2;
    plot2.Marker='o';
    plot2.MarkerFaceColor=[1 0 0];
    plot2.LineStyle='none';
    sp.FontName=FN;
    sp.FontSize=FS2;
    sp.XLabel.String='ln(s)';
    legend(['\alpha=',num2str(X.a(rgb,:),5)],'Location','southeast');
    sp.YLabel.String='ln(F_2(s))';
    
    rgb=2;
    sp=subplot(3,2,4);
    plot1=plot(log(X.s),X.a(rgb,:)*log(X.s)+X.b(rgb,:));
    hold on;
    plot2=plot(log(X.s),log(X.F2(rgb,:)),'r');
    plot1.Color=[0 0 1];
    plot1.LineWidth=2;
    plot2.Marker='o';
    plot2.MarkerFaceColor=[1 0 0];
    plot2.LineStyle='none';
    sp.FontName=FN;
    sp.FontSize=FS2;
    sp.XLabel.String='ln(s)';
    legend(['\alpha=',num2str(X.a(rgb,:),5)],'Location','southeast');
    sp.YLabel.String='ln(F_2(s))';
    
    rgb=3;
    sp=subplot(3,2,6);
    plot1=plot(log(X.s),X.a(rgb,:)*log(X.s)+X.b(rgb,:));
    hold on;
    plot2=plot(log(X.s),log(X.F2(rgb,:)),'r');
    plot1.Color=[0 0 1];
    plot1.LineWidth=2;
    plot2.Marker='o';
    plot2.MarkerFaceColor=[1 0 0];
    plot2.LineStyle='none';
    sp.FontName=FN;
    sp.FontSize=FS2;
    sp.XLabel.String='ln(s)';
    legend(['\alpha=',num2str(X.a(rgb,:),5)],'Location','southeast');
    sp.YLabel.String='ln(F_2(s))';
    
end