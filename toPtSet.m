function toPtSet(fl)

inp1=load(sprintf('%s_feaDR.txt',fl));
inp2=load(sprintf('%s_pos.txt',fl));

n=size(inp1,1);

oup=fopen(sprintf('%s_pts.txt',fl),'w');
for i=1:n
    fprintf(oup,'%f ',inp1(i,:));
    fprintf(oup,'%d ',inp2(i,:));
    fprintf(oup,'\n');
end
fclose(oup);


end