function dimreduc(fl)
inp=strcat(fl,'_fea.txt');
tobe=load(inp);


x=tobe;
y=evalin('base','pm');
 %y=load('pcaMatrix.txt');

[m,n]=size(x);
x_mean=mean(x,2);
x_var=(x-repmat(x_mean,1,n));
x_data=x_var*y;


fid=fopen(sprintf('%s_feaDR.txt',fl),'wb');
n=size(x_data,1);
for i=1:n
     fprintf(fid,'%f ',x_data(i,:));
    fprintf(fid,'\n');
end
fclose(fid);

end