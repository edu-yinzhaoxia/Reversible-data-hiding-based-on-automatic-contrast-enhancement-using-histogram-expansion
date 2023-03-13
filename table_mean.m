function [table,ref_mean,index_max,index_min]= table_mean(cover)
[m,n] = size(cover);
cover_hor = reshape(cover,m*n,1);
ref_mean = mean(cover_hor);

table = [zeros(256,1) transpose(0:255)];
for i=1:length(cover_hor)
    table(cover_hor(i)+1)=table(cover_hor(i)+1)+1;
end

[ref_max,index_max] = max(table(:,1));
 
[ref_min,index_min] = min(table(:,1));
 
end
