function forcpp(fl)

feas=load(sprintf('%s_pts.txt',fl));

fid=fopen(sprintf('%s_ptscpp.txt',fl),'wb');
n=size(feas,1);
fprintf(fid,'%d %d\n',size(feas,1),size(feas,2));
for i=1:n
    fprintf(fid,'%f ',feas(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
end