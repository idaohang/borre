s = serial('COM1');
fopen(s);
fprintf(s,'*IDN?');
idn = fscanf(s);

%proper code

fclose(s);