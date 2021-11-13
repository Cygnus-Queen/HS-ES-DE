clc 
clear all;
format shorte 
rootpath="C:\Users\Lenovo\Desktop\最终数据\result_2014_100d\";%root path
%for fun=[4 7 12 15 22:30]%需要画的函数
for fun=[4 7 12 15 16 22:30]

fg=figure;
fp1= strcat(rootpath, "1L-SHADE-EpSin\F", num2str(fun),"_Eps_100_2014.txt");  %各算法路径      
fp2 = strcat(rootpath, "2UMOEAsII\F", num2str(fun),"_original.txt");
fp3 =  strcat(rootpath, "3MPEDE\f", num2str(fun),"_original.txt");      
fp4= strcat(rootpath, "4ETI-JADE\f", num2str(fun),"_original.txt");       
fp5 =  strcat(rootpath, "5HSES\100_f",num2str(fun),".txt");   %JADE
fp6 = strcat(rootpath, "6AMECoDEs\f", num2str(fun),"_original.txt");   
fp7 =  strcat(rootpath ,"7EDEV\f", num2str(fun),"_30D.txt");        
fp8 =   strcat(rootpath ,"8MLCC-SI\f",num2str(fun),"_original.txt");    %MLCC-S
fp9 = strcat(rootpath ,"9NDE\f",num2str(fun),"_original_100D.txt");   %NDE 
fp10 =  strcat(rootpath,"9ZHS-ES-DE\100_F", num2str(fun),"_9HS-ES-DE_2014.txt");   
%set(gca,'position',[0.1 0.1 0.9 0.9])
%set (gcf,'unit','normalized','position',[0.2,0.2, 0.64,0.32])

set(0,'DefaultFigureVisible', 'off')
%x = 0:1000; 
%y = log(x);             %1
%color=[255 0 0;255 0 0;0 191 191;0 191 191;0 128 0;0 128 0;191 191 0;191 191 0;0 0 255;0 0 255];
color=[255 0 0;255 0 0;0 128 0;0 128 0;128 128 128;128 128 128;255 215 0;255 215 0;0 0 255;0 0 255];
color=color/255;
LINE=['--' ;' -'; '--' ;'- ' ;'--'; '- '; '--'; '- ' ;'--' ; '- ' ];
for alg=1:10
 h(alg)=semilogy(0:3000:300000,gettxt(eval(strcat("fp",num2str(alg)))),LINE(alg,:),'color',color(alg,:),'LineWidth',1.5);%gettxt为自定函数，提取txt文件内容
    
 %h(alg)=semilogy(0:10000:1000000,gettxt(eval(strcat("fp",num2str(alg)))),LINE(alg,:),'color',color(alg,:),'LineWidth',1.5);%gettxt为自定函数，提取txt文件内容
 hold on
end
%-\uppercase\expandafter{\romannumeral2}
 %legend ('L-SHADE-EpSin','UMOEAs-Ⅱ','MPEDE','ETI-JADE',...
 %    'HS-ES','AMECoDEs','EDEV','MLCC-SI','NDE','HS-ES-DE');
 										
 legend ('jSO','ETI-JADE','L-SHADE-RSP','EDEV',...
      'NDE','L-ETI-EDEV');
 
 
% legend([h(1),h(2)],'L-SHADE','L-SHADE-EpSin','orientation','horizontal');
% ah1=axes('position',get(gca,'position'),'visible','off');
% legend(ah1,[h(3),h(4)],'ETI-JADE','L-SHADE-RSP','orientation','horizontal');
% ah2=axes('position',get(gca,'position'),'visible','off');
% legend(ah2,[h(5),h(6)],'MPEDE','AMECoDEs','orientation','horizontal');
% ah3=axes('position',get(gca,'position'),'visible','off');
% legend(ah3,[h(7),h(8)],'EDEV','MLCC-SI','orientation','horizontal');
% ah4=axes('position',get(gca,'position'),'visible','off');
% legend(ah4,[h(9),h(10)],'NDE','L-SHADE-E','orientation','horizontal');

y_val=get(gca,'XTick');   
y_str=num2str(y_val');     
%set(gca,'XTickLabel',['0      '; '       ' ; '       ' ; '       ' ; '       ' ; '       ' ; '       ' ; '       ' ; '       ' ; '       ' ; '1000000']);
set(gca,'XTickLabel',['0      '; '       ' ; '       ' ; '       ' ; '       ' ; '       ' ; '       ' ; '       ' ; '       ' ; '       ' ; '300000']);


ylabel("The best fitness")
xlabel("Number of function evaluations")
set(gca,'fontsize',10,'linewidth',1.1)
%axis([0 1000000 35 50]);
set(gca, 'LooseInset', [0,0,0,0]);
figure(fg)
saveas(fg,char(strcat(rootpath,'figure\f',num2str(fun),'.eps')),'psc2');%图片存入根路径下的figure文件
fun
end