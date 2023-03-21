function flux = get_flux(xi,h,u)
    flux = 0;
    for i=1:(max(size(xi))-1)
        flux = flux + (u(i)*h(i)+u(i+1)*h(i+1))*(xi(i+1)-xi(i))/2;
    end
    flux = flux/xi(end);
end