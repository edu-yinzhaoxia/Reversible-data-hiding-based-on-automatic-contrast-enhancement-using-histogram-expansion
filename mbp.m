function [rdh_image, iteration_max, EC_list, LM_size_list, embedding_capacy_left]=mbp(image,payload_length,iteration_max,m,n)

if isempty(iteration_max)
        iteration_max = 1000;
end
max_contrast_bypass_mode = 0;
% if isempty(max_contrast_bypass_mode)
    % max_contrast_bypass_mode = 0;
% end

%处理
image = reshape(image,m,n);
image_size=size(image);


% rng(1) 
% actual_payload=randi([0,1],payload_length,1);


actual_payload = payload_length;


%有效载荷
payload_length_max=2*ceil(log2(image_size(1)*image_size(2)+1));
actual_payload=[de2bi(length(actual_payload),payload_length_max)'; actual_payload];

P_s_list=zeros(1,iteration_max);
P_c_list=zeros(1,iteration_max);
EC_list=zeros(1,iteration_max);
LM_size_list=zeros(1,iteration_max);
ref_image_hor = zeros(image_size(1)*image_size(2),iteration_max);
ref_image_hor(:,1) = reshape(image,image_size(1)*image_size(2),1);

image_hor=ref_image_hor(:,1);
%original_brightness = mean(ref_image_hor(:,1));

P_s=0;
P_c=0;
iteration=1;

payload_total=[];

while true

P_s_previous=P_s;
P_c_previous=P_c;
[P_s, P_c]=adaptive_peak_selection(image_hor);
 
if P_s < P_c %RHS
        d = 1;
else %LHS
        d = -1;
end
	
LM=double(image_hor(image_hor==P_c | image_hor==(P_c-d))==(P_c-d));
H_P_s=sum(image_hor==P_s);
% save H_P_s.mat

	if H_P_s-sum(image_hor(1:16)==P_s) < length(LM)+32 || iteration == iteration_max || (max_contrast_bypass_mode==0)&&length(payload_total) > length(actual_payload) 
		P_s=P_s_previous;
        P_c=P_c_previous;
        P_s_previous=P_s_list(iteration);
        P_c_previous=P_c_list(iteration);
        first_16_pixels=ref_image_hor(1:16,iteration);
        original_16_lsb=mod(first_16_pixels,2);
        
        if P_s < P_c %RHS
            d = 1;
        else %LHS
            d = -1;
        end
        

        image_hor = ref_image_hor(17:end,iteration);
        H_P_s=sum(image_hor==P_s);
        LM=(image_hor(image_hor==P_c | image_hor==P_c-d)==P_c-d);
        
        payload_total(end-payload_length_last+1:end)=[];
        if length(payload_total) < length(actual_payload)
            payload_left_over=length(actual_payload)-length(payload_total);
            if payload_left_over < H_P_s-length(LM)-32
                synthetic_payload =randi([0,1],H_P_s-length(LM)-32-payload_left_over,1);
                payload = [actual_payload(length(payload_total)+1:end); synthetic_payload];
            else
                payload = actual_payload(length(payload_total)+1:length(payload_total)+H_P_s-length(LM)-32);
            end
        else
            payload =randi([0,1],H_P_s-length(LM)-32,1);
        end

        message=[LM ; de2bi(P_s_previous,8)'; de2bi(P_c_previous,8)';original_16_lsb;payload];
        EC_list(iteration)=H_P_s-length(LM)-32;
        LM_size_list(iteration)=length(LM);
        message_whole=zeros(length(image_hor),1);
        message_whole(image_hor==P_s)=message;



        image_hor(image_hor==P_c-d)=image_hor(image_hor==P_c-d)+d;
      
        if d == 1
            image_hor(image_hor > P_s & image_hor < P_c)=image_hor(image_hor > P_s & image_hor < P_c)+d; %RHS
        else
            image_hor(image_hor < P_s & image_hor > P_c)=image_hor(image_hor < P_s & image_hor > P_c)+d; %LHS
        end


        image_hor(image_hor==P_s & message_whole)=image_hor(image_hor==P_s & message_whole)+d;
        
   
        image_hor=[bitxor(bitxor(first_16_pixels,mod(first_16_pixels,2)),[de2bi(P_s,8)'; de2bi(P_c,8)']) ;image_hor];
        iteration_max = iteration;
        rdh_image=reshape(image_hor,image_size(1)*image_size(2),1);
        
        EC_list(iteration_max+1:end)=[];
        LM_size_list(iteration_max+1:end)=[];
        ref_image_hor(:,iteration_max+1:end) = [];
        payload_total=[payload_total; payload];
        embedding_capacy_left=length(payload_total)-length(actual_payload);
        break

    else
			if length(payload_total) < length(actual_payload)
				payload_left_over=length(actual_payload)-length(payload_total);
				if payload_left_over < H_P_s-length(LM)-16
					synthetic_payload =randi([0,1],H_P_s-length(LM)-16-payload_left_over,1);
					payload = [actual_payload(length(payload_total)+1:end); synthetic_payload];
				else
					payload = actual_payload(length(payload_total)+1:length(payload_total)+H_P_s-length(LM)-16);
				end
			else
				payload =randi([0,1],H_P_s-length(LM)-16,1);
			end
			
			payload_length_last=length(payload);
			payload_total=[payload_total; payload];
				
			message=[LM ; de2bi(P_s_previous,8)'; de2bi(P_c_previous,8)';payload];
			iteration=iteration+1;
			EC_list(iteration)=H_P_s-length(LM)-16;
			LM_size_list(iteration)=length(LM);
			
			ref_image_hor(:,iteration)=image_hor;
			P_s_list(iteration)=P_s_previous;
			P_c_list(iteration)=P_c_previous;

	end
    message_whole=zeros(length(image_hor),1);
    message_whole(image_hor==P_s)=message;

    
  
    image_hor(image_hor==P_c-d)=image_hor(image_hor==P_c-d)+d;
    

    if d == 1
        image_hor(image_hor > P_s & image_hor < P_c)=image_hor(image_hor > P_s & image_hor < P_c)+d; %RHS
    else
        image_hor(image_hor < P_s & image_hor > P_c)=image_hor(image_hor < P_s & image_hor > P_c)+d; %LHS
    end
  
    image_hor(image_hor==P_s & message_whole)=image_hor(image_hor==P_s & message_whole)+d;
    
end

 end



function [P_s, P_c]=adaptive_peak_selection(image_hor)
 
[table,current_brightness,index_max,~]= table_mean(image_hor);
%save table.mat
%划定范围
[P_s,P_c] = find_range(table,current_brightness,index_max);
end



%

function [P_s,P_c] = find_range(table,current_brightness,index_max)
percent_range = 0.1;
current_brightness = fix(current_brightness);
ref_left = fix((current_brightness-index_max)*percent_range);
ref_right = fix((index_max-current_brightness)*percent_range);

right_dec = current_brightness+ref_left;
left_dec = current_brightness-ref_right;

 % per_left = 0.2;
 % per_right =per_left;
 % ref_dec_left = fix(per_left *dis_left);
 % ref_dec_right = fix(per_right*dis_right);
 % left_dec  = current_brightness-ref_dec_left;
 % right_dec = current_brightness+ref_dec_right;
 

if (left_dec>253) left_dec = 253;end
if (right_dec>253) right_dec = 253;end
if (left_dec<2) left_dec = 2;end
if(right_dec<2) right_dec = 2;end

if(current_brightness>index_max)
    
    range_min = index_max;
    range_max = right_dec;
    if(range_min<2)range_min = 2;end
    table_Ps(:,1) = table(range_min:range_max,1);
    table_Ps(:,2) = table(range_min:range_max,2);
    [~,ps_index] = max(table_Ps(:,1));
    P_s=table_Ps(ps_index,2)


    table_Pc(:,1) = table(2:range_min+ref_left,1);
    table_Pc(:,2) = table(2:range_min+ref_left,2);
    %处理min 选择离均值最远的
    combine_table = table_Pc;
    [sort_combine_table, sort_combine_table_index]=sort(combine_table(:,1));
    sort_combine_table = [sort_combine_table combine_table(sort_combine_table_index,2)];
    list_P_c=sort_combine_table(sort_combine_table(:,1)==sort_combine_table(1,1),:);
    [~,index]=min(abs(list_P_c(:,2)-current_brightness));
    P_c=list_P_c(index,2);
    
else(current_brightness<=index_max)
  
    range_min = left_dec;
    range_max = index_max;
	  if(range_max>253) range_max = 253;end
    table_Ps(:,1) = table(range_min:range_max,1);
    table_Ps(:,2) = table(range_min:range_max,2);
    table_Pc(:,1) = table(range_max-ref_right:253,1);
    table_Pc(:,2) = table(range_max-ref_right:253,2);
	
    [~,ps_index] = max(table_Ps(:,1));
    P_s=table_Ps(ps_index,2)
    combine_table = table_Pc;
    [sort_combine_table, sort_combine_table_index]=sort(combine_table(:,1));
    sort_combine_table = [sort_combine_table combine_table(sort_combine_table_index,2)];
    list_P_c=sort_combine_table(sort_combine_table(:,1)==sort_combine_table(1,1),:);
    [~,index]=min(abs(list_P_c(:,2)-current_brightness));
    P_c=list_P_c(index,2);
end 

end