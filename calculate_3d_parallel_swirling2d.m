function calculate_3d_parallel_swirling2d(which_component, k_z, k_rms2, kz_max, kx_detection_array, ky_detection_array, Retau, read_array, jobid, workers)
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

[~,~,~,swirling_spanwise_slice] = get_uvwswirl2d_z(which_component,k_z,Retau,read_array);
disp("all good now")

u_cd_all     = zeros(nz-1,ny,nx,length(ky_detection_array));
v_cd_all     = zeros(nz-1,ny,nx,length(ky_detection_array));
w_cd_all     = zeros(nz-1,ny,nx,length(ky_detection_array));
uw_cd_all    = zeros(nz-1,ny,nx,length(ky_detection_array));
swirling_all = zeros(nz-1,ny,nx,length(ky_detection_array));
number       = zeros(length(ky_detection_array),length(kx_detection_array)); 

parfor kx_detection_index = 1: length(kx_detection_array)
    kx_detection = kx_detection_array(kx_detection_index);
    disp(kx_detection)
    temp_u_cd_all     = zeros(nz-1,ny,nx); 
    temp_v_cd_all     = zeros(nz-1,ny,nx);
    temp_w_cd_all     = zeros(nz-1,ny,nx); 
    temp_uw_cd_all    = zeros(nz-1,ny,nx);
    temp_swirling_all = zeros(nz-1,ny,nx);
    temp_number       = zeros(length(ky_detection_array),1);
   
    for ky_detection_index = 1: length(ky_detection_array)
        ky_detection = ky_detection_array(ky_detection_index);
        swirling_spanwise_slice2 = squeeze(swirling_spanwise_slice(ky_detection,kx_detection,:));

        pick_swirling_spanwise = get_detection_events(swirling_spanwise_slice2,k_rms2,4,1);
        
        pick_swirling_spanwise_index = find(pick_swirling_spanwise);
        
        temp_number(ky_detection_index) = length(pick_swirling_spanwise_index);

        u_cd     = zeros(nz-1,ny,nx);
        v_cd     = zeros(nz-1,ny,nx);
        w_cd     = zeros(nz-1,ny,nx);
        uw_cd    = zeros(nz-1,ny,nx);
        swirling = zeros(nz-1,ny,nx);
        for k_array_index = 1: length(pick_swirling_spanwise_index)
            k_array = pick_swirling_spanwise_index(k_array_index);
            u = zeros(nzDNS+2,ny,nx);
            v = zeros(nzDNS+2,ny,nx);
            w = zeros(nzDNS+1,ny,nx);
            loadname = strcat('../ChanFast/grid_',loadname1,'/outputdir/u_it',num2str(read_array(k_array),'%.0f'),'.dat');
            u(2:end-1,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
            u = interp1(zp,u,zc);
            u = u(2:end-1,:,:);
            u = permute(u,[3 2 1]);
            u = interp1(xu,u,xp,'linear','extrap');
            u = permute(u,[3 2 1]);
            loadname = strcat('../ChanFast/grid_',loadname1,'/outputdir/v_it',num2str(read_array(k_array),'%.0f'),'.dat');
            v(2:end-1,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
            v = interp1(zp,v,zc);
            v = v(2:end-1,:,:);
            v = permute(v,[2 1 3]);
            v = interp1(yv,v,yp,'linear','extrap');
            v = permute(v, [2 1 3]);
            loadname = strcat('../ChanFast/grid_',loadname1,'/outputdir/w_it',num2str(read_array(k_array),'%.0f'),'.dat');
            w(2:end,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
            w = interp1(zw,w,zc);
            w = w(2:end-1,:,:);
        
            [u] = get_new_3d_box(u,kx_detection,ky_detection,kx_middle,ky_middle);
            [v] = get_new_3d_box(v,kx_detection,ky_detection,kx_middle,ky_middle);
            [w] = get_new_3d_box(w,kx_detection,ky_detection,kx_middle,ky_middle);
        
            [velocity_tensor] = real(get_velocity_tensor(u,v,w,kx_array,ky_array,nx,ny,nz,dkx,dky));
            velocity_tensor = permute(velocity_tensor,[4,1,2,3]); 
            [swirling_strength] = get_swirling_strength(velocity_tensor,nx,ny,nz);
        
            u_cd = (k_array_index-1).*u_cd./k_array_index + u./k_array_index;
            v_cd = (k_array_index-1).*v_cd./k_array_index + v./k_array_index;
            w_cd = (k_array_index-1).*w_cd./k_array_index + w./k_array_index;
            uw_cd= (k_array_index-1).*uw_cd./k_array_index + (u.*w)./k_array_index;
            swirling = (k_array_index-1).*swirling./k_array_index + (swirling_strength)./k_array_index;
        end

        temp_u_cd_all     = length(pick_swirling_spanwise_index) .* u_cd + temp_u_cd_all;
        temp_v_cd_all     = length(pick_swirling_spanwise_index) .* v_cd + temp_v_cd_all;
        temp_w_cd_all     = length(pick_swirling_spanwise_index) .* w_cd + temp_w_cd_all;
        temp_uw_cd_all    = length(pick_swirling_spanwise_index) .* uw_cd + temp_uw_cd_all;
        temp_swirling_all = length(pick_swirling_spanwise_index) .* swirling + temp_swirling_all;
    end


    u_cd_all(:,:,:,kx_detection_index)        = single(temp_u_cd_all);
    v_cd_all(:,:,:,kx_detection_index)        = single(temp_v_cd_all);
    w_cd_all(:,:,:,kx_detection_index)        = single(temp_w_cd_all);
    uw_cd_all(:,:,:,kx_detection_index)       = single(temp_uw_cd_all);
    swirling_all(:,:,:,kx_detection_index)    = single(temp_swirling_all);
    number(:,kx_detection_index)              = temp_number;
end

delete(gcp('nocreate'))

u_cd_all     = sum(u_cd_all,4) ./ sum(number,'all');
v_cd_all     = sum(v_cd_all,4) ./ sum(number,'all');
w_cd_all     = sum(w_cd_all,4) ./ sum(number,'all');
uw_cd_all    = sum(uw_cd_all,4) ./ sum(number,'all');
swirling_all = sum(swirling_all,4) ./ sum(number,'all');

[velocity_tensor] = real(get_velocity_tensor(u_cd_all,v_cd_all,w_cd_all,kx_array,ky_array,nx,ny,nz,dkx,dky));
velocity_tensor = permute(velocity_tensor,[4,1,2,3]); 
[swirling_aftercd] = get_swirling_strength(velocity_tensor,nx,ny,nz);

savename = strcat('store_pick/store_pick_3d_swirling2d',num2str(k_z),'_all',num2str(k_rms2),'_',num2str(Retau),'_',num2str(jobid),'.mat');
save(savename,'u_cd_all','v_cd_all','w_cd_all','uw_cd_all','swirling_all','swirling_aftercd','-v7.3');

end
