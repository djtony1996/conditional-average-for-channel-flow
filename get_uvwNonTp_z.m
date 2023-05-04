function [u_slice,v_slice,w_slice,NonTp_slice,Prop_slice,Dissp_slice] = get_uvwNonTp_z(k_z,Retau,read_array) 
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
kx_array = temp.kx_array;
ky_array = temp.ky_array;
dkx = temp.dkx;
dky = temp.dky;
dUdz = temp.channelRe.Up_diff1(2:end-1);

[Diff,zc] = cheb(nz);
Diff = Diff(2:end-1,2:end-1);

NonTp_slice = zeros(ny,nx,length(read_array));
Prop_slice  = zeros(ny,nx,length(read_array));
Dissp_slice = zeros(ny,nx,length(read_array));
u_slice     = zeros(ny,nx,length(read_array));
v_slice     = zeros(ny,nx,length(read_array));
w_slice     = zeros(ny,nx,length(read_array));
% p_slice     = zeros(ny,nx,length(read_array));

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
    % loadname = strcat('../ChanFast/grid_',loadname1,'/outputdir/p_it',num2str(read_array(k_array),'%.0f'),'.dat');
    % p = read_bin(loadname, [nzDNS, ny, nx]);
    % p = interp1(zp(2:end-1),p,zc(2:end-1));

    u_slice(:,:,k_array) = squeeze(u(k_z,:,:));
    v_slice(:,:,k_array) = squeeze(v(k_z,:,:));
    w_slice(:,:,k_array) = squeeze(w(k_z,:,:));
    % p_slice(:,:,k_array) = squeeze(p(k_z,:,:));

    u_F  = fft_xy(u);
    v_F  = fft_xy(v);
    w_F  = fft_xy(w);

    [kx_m,ky_m] = getkm(kx_array,ky_array,nx,ny,dkx,dky);
    [duF_dx,duF_dy,duF_dz]    = get_3d(u_F,Diff,kx_m,ky_m);
    [dvF_dx,dvF_dy,dvF_dz]    = get_3d(v_F,Diff,kx_m,ky_m);
    [dwF_dx,dwF_dy,dwF_dz]    = get_3d(w_F,Diff,kx_m,ky_m);

    NonTxp = -u.*u.*ifft_xy(duF_dx) - u.*v.*ifft_xy(duF_dy) - u.*w.*ifft_xy(duF_dz);
    NonTyp = -v.*u.*ifft_xy(dvF_dx) - v.*v.*ifft_xy(dvF_dy) - v.*w.*ifft_xy(dvF_dz);
    NonTzp = -w.*u.*ifft_xy(dwF_dx) - w.*v.*ifft_xy(dwF_dy) - w.*w.*ifft_xy(dwF_dz);
    NonTp  = NonTxp + NonTyp + NonTzp;
    NonTp_slice(:,:,k_array) = squeeze(NonTp(k_z,:,:));
    
    Dissxp = -ifft_xy(duF_dx).*ifft_xy(duF_dx) - ifft_xy(duF_dy).*ifft_xy(duF_dy) - ifft_xy(duF_dz).*ifft_xy(duF_dz);
    Dissyp = -ifft_xy(dvF_dx).*ifft_xy(dvF_dx) - ifft_xy(dvF_dy).*ifft_xy(dvF_dy) - ifft_xy(dvF_dz).*ifft_xy(dvF_dz);
    Disszp = -ifft_xy(dwF_dx).*ifft_xy(dwF_dx) - ifft_xy(dwF_dy).*ifft_xy(dwF_dy) - ifft_xy(dwF_dz).*ifft_xy(dwF_dz);
    Dissp  = Dissxp + Dissyp + Disszp;
    Dissp_slice(:,:,k_array) = squeeze(Dissp(k_z,:,:));

    Proxp = -u.*w.*dUdz;
    Prop_slice(:,:,k_array) = squeeze(Proxp(k_z,:,:));
end
u_slice = single(u_slice);
v_slice = single(v_slice);
w_slice = single(w_slice);
% p_slice = single(p_slice);
NonTp_slice = single(real(NonTp_slice));
Prop_slice  = single(Prop_slice);
Dissp_slice = single(Dissp_slice);
end
