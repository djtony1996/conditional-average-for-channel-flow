function [u] = get_new_3d_box(u,kx,ky,kx_middle,ky_middle)

if ky<ky_middle
    temp_u = u(:,1:end-ky_middle+ky,:);
    u = cat(2,u(:,(end-ky_middle+ky+1):end,:),temp_u);
else
    temp_u = u(:,(1+ky-ky_middle):end,:);
    u = cat(2,temp_u,u(:,1:(ky-ky_middle),:));
end

if kx<kx_middle
    temp_u = u(:,:,1:end-kx_middle+kx);
    u = cat(3,u(:,:,(end-kx_middle+kx+1):end),temp_u);
else 
    temp_u = u(:,:,(1+kx-kx_middle):end);
    u = cat(3,temp_u,u(:,:,1:(kx-kx_middle)));
end

end