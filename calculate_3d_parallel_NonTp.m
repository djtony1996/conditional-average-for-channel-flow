function calculate_3d_parallel_NonTp(k_z, k_rms2, kz_max, kx_detection_array, ky_detection_array, Retau, read_array, jobid, workers)
myCluster = parcluster('local');
myCluster.NumWorkers = workers;
parpool(workers)

kx_middle = 56;
ky_middle = 56;

if Retau == 180
    kz_detection = 111;
    loadname1 = '180/112x112x150';
elseif Retau == 395
    kz_detection = 232;
    loadname1 = '395/256x256x300';
elseif Retau == 590
    kz_detection = 237;
    loadname1 = '590/384x384x500';
end

temp = load(strcat('full',num2str(Retau),'_mean.mat'));
nx = temp.nx;
ny = temp.ny;
nz = temp.nz;
nzDNS = temp.nzDNS;
zp = temp.zp;
xu = temp.xu;
xp = temp.xp;
yv = temp.yv;
yp = temp.yp;
zw = temp.zw;
dkx = temp.dkx;
dky = temp.dky;
kx = [0:dkx:(50*dkx),-(50*dkx):dkx:-dkx];
ky = [0:dky:(50*dky)];
kx_array = repelem(kx,length(ky));
ky_array = repmat(ky,1,length(kx));
kx_array = kx_array.';
ky_array = ky_array.';

[Diff,zc] = cheb(nz);
Diff = Diff(2:end-1,2:end-1);

[~,~,~,NonTp_slice] = get_uvwNonTp_z(k_z,Retau, read_array);
disp("all good now")

u_cd_posi_all     = zeros(nz-1,ny,nx,length(ky_detection_array));
v_cd_posi_all     = zeros(nz-1,ny,nx,length(ky_detection_array));
w_cd_posi_all     = zeros(nz-1,ny,nx,length(ky_detection_array));
uw_cd_posi_all    = zeros(nz-1,ny,nx,length(ky_detection_array));
swirling_posi_all = zeros(nz-1,ny,nx,length(ky_detection_array));
u_cd_nega_all     = zeros(nz-1,ny,nx,length(ky_detection_array));
v_cd_nega_all     = zeros(nz-1,ny,nx,length(ky_detection_array));
w_cd_nega_all     = zeros(nz-1,ny,nx,length(ky_detection_array));
uw_cd_nega_all    = zeros(nz-1,ny,nx,length(ky_detection_array));
swirling_nega_all = zeros(nz-1,ny,nx,length(ky_detection_array));
number_posi       = zeros(length(ky_detection_array),length(kx_detection_array)); 
number_nega       = zeros(length(ky_detection_array),length(kx_detection_array));

parfor kx_detection_index = 1: length(kx_detection_array)
    kx_detection = kx_detection_array(kx_detection_index);
    disp(kx_detection)
    temp_u_cd_posi_all     = zeros(nz-1,ny,nx); 
    temp_v_cd_posi_all     = zeros(nz-1,ny,nx);
    temp_w_cd_posi_all     = zeros(nz-1,ny,nx); 
    temp_uw_cd_posi_all    = zeros(nz-1,ny,nx);
    temp_swirling_posi_all = zeros(nz-1,ny,nx);
    temp_u_cd_nega_all     = zeros(nz-1,ny,nx);
    temp_v_cd_nega_all     = zeros(nz-1,ny,nx);
    temp_w_cd_nega_all     = zeros(nz-1,ny,nx);
    temp_uw_cd_nega_all    = zeros(nz-1,ny,nx);
    temp_swirling_nega_all = zeros(nz-1,ny,nx);
    temp_number_posi       = zeros(length(ky_detection_array),1);
    temp_number_nega       = zeros(length(ky_detection_array),1);
   
    for ky_detection_index = 1: length(ky_detection_array)
        ky_detection = ky_detection_array(ky_detection_index);
        NonTp_slice2 = squeeze(NonTp_slice(ky_detection,kx_detection,:));

        pick_NonTp = get_detection_events(NonTp_slice2,k_rms2,1,1);
        
        NonTp_posi = NonTp_slice2>0;
        NonTp_nega = NonTp_slice2<0;
        
        pick_NonTp_posi = NonTp_posi .* pick_NonTp;
        pick_NonTp_nega = NonTp_nega .* pick_NonTp .* (-1);
        
        pick_NonTp_posi_index = find(pick_NonTp_posi);
        pick_NonTp_nega_index = find(pick_NonTp_nega);
        
        temp_number_posi(ky_detection_index) = length(pick_NonTp_posi_index);
        temp_number_nega(ky_detection_index) = length(pick_NonTp_nega_index);

        u_cd_posi     = zeros(nz-1,ny,nx);
        v_cd_posi     = zeros(nz-1,ny,nx);
        w_cd_posi     = zeros(nz-1,ny,nx);
        uw_cd_posi    = zeros(nz-1,ny,nx);
        swirling_posi = zeros(nz-1,ny,nx);
        u_cd_nega     = zeros(nz-1,ny,nx);
        v_cd_nega     = zeros(nz-1,ny,nx);
        w_cd_nega     = zeros(nz-1,ny,nx);
        uw_cd_nega    = zeros(nz-1,ny,nx);
        swirling_nega = zeros(nz-1,ny,nx);
        for k_array_index = 1: length(pick_NonTp_posi_index)
            k_array = pick_NonTp_posi_index(k_array_index);
            u = zeros(nzDNS+2,ny,nx);
            v = zeros(nzDNS+2,ny,nx);
            w = zeros(nzDNS+1,ny,nx);
            loadname = strcat('../ChanFast/grid_',loadname1,'/outputdir/u_it',num2str(read_array(k_array),'%.0f'),'.dat');
            u(2:end-1,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
            u = interp1(zp,u,zc);
            u = u - mean(mean(u,3),2);
            u = u(2:end-1,:,:);
            u = permute(u,[3 2 1]);
            u = interp1(xu,u,xp,'linear','extrap');
            u = permute(u,[3 2 1]);
            loadname = strcat('../ChanFast/grid_',loadname1,'/outputdir/v_it',num2str(read_array(k_array),'%.0f'),'.dat');
            v(2:end-1,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
            v = interp1(zp,v,zc);
            v = v - mean(mean(v,3),2);
            v = v(2:end-1,:,:);
            v = permute(v,[2 1 3]);
            v = interp1(yv,v,yp,'linear','extrap');
            v = permute(v, [2 1 3]);
            loadname = strcat('../ChanFast/grid_',loadname1,'/outputdir/w_it',num2str(read_array(k_array),'%.0f'),'.dat');
            w(2:end,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
            w = interp1(zw,w,zc);
            w = w - mean(mean(w,3),2);
            w = w(2:end-1,:,:);
        
            [u] = get_new_3d_box(u,kx_detection,ky_detection,kx_middle,ky_middle);
            [v] = get_new_3d_box(v,kx_detection,ky_detection,kx_middle,ky_middle);
            [w] = get_new_3d_box(w,kx_detection,ky_detection,kx_middle,ky_middle);
        
            [velocity_tensor] = real(get_velocity_tensor(u,v,w,kx_array,ky_array,nx,ny,nz,dkx,dky));
            velocity_tensor = permute(velocity_tensor,[4,1,2,3]); 
            [swirling_strength] = get_swirling_strength(velocity_tensor,nx,ny,nz);
        
            u_cd_posi = (k_array_index-1).*u_cd_posi./k_array_index + u./k_array_index;
            v_cd_posi = (k_array_index-1).*v_cd_posi./k_array_index + v./k_array_index;
            w_cd_posi = (k_array_index-1).*w_cd_posi./k_array_index + w./k_array_index;
            uw_cd_posi= (k_array_index-1).*uw_cd_posi./k_array_index + (u.*w)./k_array_index;
            swirling_posi = (k_array_index-1).*swirling_posi./k_array_index + (swirling_strength)./k_array_index;
        end

        temp_u_cd_posi_all     = length(pick_NonTp_posi_index) .* u_cd_posi + temp_u_cd_posi_all;
        temp_v_cd_posi_all     = length(pick_NonTp_posi_index) .* v_cd_posi + temp_v_cd_posi_all;
        temp_w_cd_posi_all     = length(pick_NonTp_posi_index) .* w_cd_posi + temp_w_cd_posi_all;
        temp_uw_cd_posi_all    = length(pick_NonTp_posi_index) .* uw_cd_posi + temp_uw_cd_posi_all;
        temp_swirling_posi_all = length(pick_NonTp_posi_index) .* swirling_posi + temp_swirling_posi_all;
        
        for k_array_index = 1: length(pick_NonTp_nega_index)
            k_array = pick_NonTp_nega_index(k_array_index);
            u = zeros(nzDNS+2,ny,nx);
            v = zeros(nzDNS+2,ny,nx);
            w = zeros(nzDNS+1,ny,nx);
            loadname = strcat('../ChanFast/grid_',loadname1,'/outputdir/u_it',num2str(read_array(k_array),'%.0f'),'.dat');
            u(2:end-1,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
            u = interp1(zp,u,zc);
            u = u - mean(mean(u,3),2);
            u = u(2:end-1,:,:);
            u = permute(u,[3 2 1]);
            u = interp1(xu,u,xp,'linear','extrap');
            u = permute(u,[3 2 1]);
            loadname = strcat('../ChanFast/grid_',loadname1,'/outputdir/v_it',num2str(read_array(k_array),'%.0f'),'.dat');
            v(2:end-1,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
            v = interp1(zp,v,zc);
            v = v - mean(mean(v,3),2);
            v = v(2:end-1,:,:);
            v = permute(v,[2 1 3]);
            v = interp1(yv,v,yp,'linear','extrap');
            v = permute(v, [2 1 3]);
            loadname = strcat('../ChanFast/grid_',loadname1,'/outputdir/w_it',num2str(read_array(k_array),'%.0f'),'.dat');
            w(2:end,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
            w = interp1(zw,w,zc);
            w = w - mean(mean(w,3),2);
            w = w(2:end-1,:,:);
        
            [u] = get_new_3d_box(u,kx_detection,ky_detection,kx_middle,ky_middle);
            [v] = get_new_3d_box(v,kx_detection,ky_detection,kx_middle,ky_middle);
            [w] = get_new_3d_box(w,kx_detection,ky_detection,kx_middle,ky_middle);
        
            [velocity_tensor] = real(get_velocity_tensor(u,v,w,kx_array,ky_array,nx,ny,nz,dkx,dky));
            velocity_tensor = permute(velocity_tensor,[4,1,2,3]); 
            [swirling_strength] = get_swirling_strength(velocity_tensor,nx,ny,nz);
        
            u_cd_nega = (k_array_index-1).*u_cd_nega./k_array_index + u./k_array_index;
            v_cd_nega = (k_array_index-1).*v_cd_nega./k_array_index + v./k_array_index;
            w_cd_nega = (k_array_index-1).*w_cd_nega./k_array_index + w./k_array_index;
            uw_cd_nega= (k_array_index-1).*uw_cd_nega./k_array_index + (u.*w)./k_array_index;
            swirling_nega = (k_array_index-1).*swirling_nega./k_array_index + (swirling_strength)./k_array_index;
        end

        temp_u_cd_nega_all     = length(pick_NonTp_nega_index) .* u_cd_nega + temp_u_cd_nega_all;
        temp_v_cd_nega_all     = length(pick_NonTp_nega_index) .* v_cd_nega + temp_v_cd_nega_all;
        temp_w_cd_nega_all     = length(pick_NonTp_nega_index) .* w_cd_nega + temp_w_cd_nega_all;
        temp_uw_cd_nega_all    = length(pick_NonTp_nega_index) .* uw_cd_nega + temp_uw_cd_nega_all;
        temp_swirling_nega_all = length(pick_NonTp_nega_index) .* swirling_nega + temp_swirling_nega_all;
    end


    u_cd_posi_all(:,:,:,kx_detection_index)        = single(temp_u_cd_posi_all);
    v_cd_posi_all(:,:,:,kx_detection_index)        = single(temp_v_cd_posi_all);
    w_cd_posi_all(:,:,:,kx_detection_index)        = single(temp_w_cd_posi_all);
    uw_cd_posi_all(:,:,:,kx_detection_index)       = single(temp_uw_cd_posi_all);
    swirling_posi_all(:,:,:,kx_detection_index)    = single(temp_swirling_posi_all);
    u_cd_nega_all(:,:,:,kx_detection_index)        = single(temp_u_cd_nega_all);
    v_cd_nega_all(:,:,:,kx_detection_index)        = single(temp_v_cd_nega_all);
    w_cd_nega_all(:,:,:,kx_detection_index)        = single(temp_w_cd_nega_all);
    uw_cd_nega_all(:,:,:,kx_detection_index)       = single(temp_uw_cd_nega_all);
    swirling_nega_all(:,:,:,kx_detection_index)    = single(temp_swirling_nega_all);
    number_posi(:,kx_detection_index)              = temp_number_posi;
    number_nega(:,kx_detection_index)              = temp_number_nega;
end

delete(gcp('nocreate'))

u_cd_posi_all     = sum(u_cd_posi_all,4) ./ sum(number_posi,'all');
v_cd_posi_all     = sum(v_cd_posi_all,4) ./ sum(number_posi,'all');
w_cd_posi_all     = sum(w_cd_posi_all,4) ./ sum(number_posi,'all');
uw_cd_posi_all    = sum(uw_cd_posi_all,4) ./ sum(number_posi,'all');
swirling_posi_all = sum(swirling_posi_all,4) ./ sum(number_posi,'all');
u_cd_nega_all     = sum(u_cd_nega_all,4) ./ sum(number_nega,'all');
v_cd_nega_all     = sum(v_cd_nega_all,4) ./ sum(number_nega,'all');
w_cd_nega_all     = sum(w_cd_nega_all,4) ./ sum(number_nega,'all');
uw_cd_nega_all    = sum(uw_cd_nega_all,4) ./ sum(number_nega,'all');
swirling_nega_all = sum(swirling_nega_all,4) ./ sum(number_nega,'all');

[velocity_tensor] = real(get_velocity_tensor(u_cd_nega_all,v_cd_nega_all,w_cd_nega_all,kx_array,ky_array,nx,ny,nz,dkx,dky));
velocity_tensor = permute(velocity_tensor,[4,1,2,3]); 
[swirling_nega_aftercd] = get_swirling_strength(velocity_tensor,nx,ny,nz);

[velocity_tensor] = real(get_velocity_tensor(u_cd_posi_all,v_cd_posi_all,w_cd_posi_all,kx_array,ky_array,nx,ny,nz,dkx,dky));
velocity_tensor = permute(velocity_tensor,[4,1,2,3]); 
[swirling_posi_aftercd] = get_swirling_strength(velocity_tensor,nx,ny,nz);

savename = strcat('store_pick/store_pick_3d_',num2str(k_z),'_all',num2str(k_rms2),'_',num2str(Retau),'_',num2str(jobid),'.mat');
save(savename,'u_cd_posi_all','v_cd_posi_all','w_cd_posi_all','uw_cd_posi_all','swirling_posi_all','u_cd_nega_all','v_cd_nega_all','w_cd_nega_all','uw_cd_nega_all','swirling_nega_all','number_posi','number_nega','swirling_nega_aftercd','swirling_posi_aftercd','-v7.3');

end
