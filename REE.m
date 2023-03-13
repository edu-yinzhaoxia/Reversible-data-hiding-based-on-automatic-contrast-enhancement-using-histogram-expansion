function ree = REE(I1,I2)

   L=256;
   
   ree=(info_entropies(I1)-info_entropies(I2))/(2*log2(L))+0.5;

end   