fid = fopen('gmsd1660_cut_4sat_2.13o', 'r');
gid = fopen('gmsd1660_cut_4sat_3.13o', 'w');
tline = fgetl(fid);

while ischar(tline)
    if tline(1) == '>'
        fprintf(gid, tline(1:33));
        fprintf(gid, '4\n');
    else
        fprintf(gid, tline);
        fprintf(gid, '\n');
    end
    tline = fgetl(fid);
end
fclose(fid);
fclose(gid);
     