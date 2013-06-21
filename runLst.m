function runLst(fl,cmd)
inp=fopen(fl)
num=fscanf(inp,'%d')

for i=1:num
    tline = fgetl(inp);
    rcmd=strcat(cmd,'(',char(39),tline,char(39),')')
    eval(rcmd);
end
fclose(inp);
end