function [u_slice,v_slice,w_slice,swirling_spanwise_slice] = get_uvwswirl2d_z(which_component,k_z,Retau,read_array) 
fft_x = @(x) fft(x, [], 3) / size(x, 3);
fft_y = @(x) fft(x, [], 2) / size(x, 2);
fft_xy = @(x) fft_x(fft_y(x));
ifft_x = @(x) ifft(x, [], 3) * size(x, 3);
ifft_y = @(x) ifft(x, [], 2) * size(x, 2);
ifft_xy = @(x) ifft_x(ifft_y(x));

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
dUdz = temp.channelRe.Up_diff1(2:end-1);

kx = [0:dkx:(50*dkx),-(50*dkx):dkx:-dkx];
ky = [0:dky:(50*dky)];
kx_array = repelem(kx,length(ky));
ky_array = repmat(ky,1,length(kx));
kx_array = kx_array.';
ky_array = ky_array.';

[Diff,zc] = cheb(nz);
Diff = Diff(2:end-1,2:end-1);

swirling_spanwise_slice = zeros(ny,nx,length(read_array));
u_slice     = zeros(ny,nx,length(read_array));
v_slice     = zeros(ny,nx,length(read_array));
w_slice     = zeros(ny,nx,length(read_array));

if Retau == 180
    loadname1 = '180/112x112x150';
elseif Retau == 395
    loadname1 = '395/256x256x300';
elseif Retau == 590
    loadname1 = '590/384x384x500';
end

for k_array = 1:length(read_array)
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

    u_slice(:,:,k_array) = squeeze(u(k_z,:,:));
    v_slice(:,:,k_array) = squeeze(v(k_z,:,:));
    w_slice(:,:,k_array) = squeeze(w(k_z,:,:));

    [velocity_tensor]   = real(get_velocity_tensor(u,v,w,kx_array,ky_array,nx,ny,nz,dkx,dky));
    velocity_tensor     = permute(velocity_tensor,[4,1,2,3]); 
    [swirling_spanwise] = get_swirling_strength_2d(velocity_tensor,nx,ny,nz,which_component);

    swirling_spanwise_slice(:,:,k_array) = squeeze(swirling_spanwise(k_z,:,:));
end
u_slice = single(u_slice);
v_slice = single(v_slice);
w_slice = single(w_slice);
swirling_spanwise_slice = single(real(swirl2d_slice));
end
