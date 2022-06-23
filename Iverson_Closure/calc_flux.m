function flux = calc_flux(u_vec)
    dz=1/(length(u_vec)-0.5);
    flux = u_vec(1)+u_vec(end);
    for i=2:(length(u_vec)-1)
        flux = flux + 2*u_vec(i);
    end
    flux = flux*dz/2;
end