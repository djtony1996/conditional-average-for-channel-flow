function [velocity_tensor] = get_velocity_tensor(u,v,w,kx_array,ky_array,nx,ny,nz,dkx,dky)

fft_x = @(x) fft(x, [], 3) / size(x, 3);
fft_y = @(x) fft(x, [], 2) / size(x, 2);
fft_xy = @(x) fft_x(fft_y(x));

ifft_x = @(x) ifft(x, [], 3) * size(x, 3);
ifft_y = @(x) ifft(x, [], 2) * size(x, 2);
ifft_xy = @(x) ifft_x(ifft_y(x));

[Diff,zc] = cheb(nz);
Diff = Diff(2:end-1,2:end-1);

u_F  = fft_xy(u);
v_F  = fft_xy(v);
w_F  = fft_xy(w);

[kx_m,ky_m] = getkm(kx_array,ky_array,nx,ny,dkx,dky);

[duF_dx,duF_dy,duF_dz]    = get_3d(u_F,Diff,kx_m,ky_m);
[dvF_dx,dvF_dy,dvF_dz]    = get_3d(v_F,Diff,kx_m,ky_m);
[dwF_dx,dwF_dy,dwF_dz]    = get_3d(w_F,Diff,kx_m,ky_m);

du_dx = ifft_xy(duF_dx);    du_dy = ifft_xy(duF_dy);    du_dz = ifft_xy(duF_dz);
dv_dx = ifft_xy(dvF_dx);    dv_dy = ifft_xy(dvF_dy);    dv_dz = ifft_xy(dvF_dz);
dw_dx = ifft_xy(dwF_dx);    dw_dy = ifft_xy(dwF_dy);    dw_dz = ifft_xy(dwF_dz);

velocity_tensor = cat(4,du_dx,dv_dx,dw_dx,du_dy,dv_dy,dw_dy,du_dz,dv_dz,dw_dz);
end