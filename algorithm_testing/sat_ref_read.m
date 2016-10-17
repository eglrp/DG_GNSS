
fid=fopen('data.sp3','r');

line = fgetl(fid);

n = 0;

while ~feof(fid)
    if line(1) == '*'
        n = n + 1;  
        [d] = sscanf(line, '%c %d %d %d %d %d %lf');
        date(n,1:5) = d(2:6);
        line = fgetl(fid);
    elseif line(1) == 'P'
        [p] = sscanf(line, '%c %c %d %lf %lf %lf %f');
        sat_ind = find(p(3) == k);
        if ~isempty(sat_ind)
            XS(sat_ind,:,n) = 1000*p(4:6);
            dtS(sat_ind,n) = (10^-6)*p(7);
        end
        line = fgetl(fid);
    end
    %line = fgetl(fid);
end

    
