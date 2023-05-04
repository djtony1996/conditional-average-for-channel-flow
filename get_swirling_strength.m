function [swirling_strength] = get_swirling_strength(velocity_tensor,nx,ny,nz)

swirling_strength = zeros(nz-1,ny,nx);
for kz = 1: nz-1
    for ky = 1: ny
        for kx = 1: nx
            velocity_tensor_onepoint = reshape(velocity_tensor(:,kz,ky,kx),[3,3]);
            swirling_strength(kz,ky,kx) = max(imag(eig(velocity_tensor_onepoint))); %[z,y,x]
        end
    end
end

end