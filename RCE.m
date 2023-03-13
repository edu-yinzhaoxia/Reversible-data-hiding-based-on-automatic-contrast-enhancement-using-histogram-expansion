function rce =RCE(I1,I2)

   L=256;
   rce =(std2(I1)-std2(I2))/(L-1)+0.5; 
   
   end