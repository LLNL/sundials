function pstatus = get_pstatus()

fid = fopen('/proc/self/status');

while ~feof(fid)
  line = fgetl(fid);
  sp = isspace(line);
  indx = find(sp);
  fname = line(1:indx(1)-2);
  fval = line(indx(1)+1:length(line));
  pstatus.(fname)=fval;
end
