function print_switch_table (p8, filename)

file=fopen(filename, 'w');
fprintf(file, "#include <ogg/os_types.h>\n");
fprintf(file, "const ogg_uint16_t od_switch_size8_cdf[][16] = {\n");
for k=1:length(p8)
   p = max(1,round(p8(:,k)*32768));
   s = sum(p);
   [a,b]=max(p);
   p(b) = p(b)+32768-s;
   if (sum(p)!=32768)
      fprintf(stderr, "failed to renormalize\n");
   end
   p = filter(1, [1,-1], p);
   fprintf(file, "  {");
   fprintf(file, "%d, ", p(1:8));
   fprintf(file, "\n   ");
   fprintf(file, "%d, ", p(9:end-1));
   fprintf(file, "%d", p(end));
   fprintf(file, "}");
   if (k!=length(p8))
      fprintf(file,",");
   end
   fprintf(file, "\n");
end
fprintf(file, "};\n");
fclose(file);
end
