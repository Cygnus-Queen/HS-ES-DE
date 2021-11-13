function [   meannum ] = gettxt(fn)
%提取txt中的内容
times=30;
hangshu=101;
fp1 =fopen (fn, 'rt');   
  count=1;
  while count~=hangshu*30+1%设置读取数据总行数+1
      str=strsplit( fgetl(fp1))  ;%按列分隔

            num(count)=str2num(str{2})  ;%取第五列数据
  
       count=count+1;
  end
  meannum=[];

  for line=1:101
       num2=[];
     for i=1:times                                   
          num2=[num2 num(line+hangshu*(i-1)) ];
     end
     meannum(line)=mean(num2);
  end
   fclose(fp1);
end

