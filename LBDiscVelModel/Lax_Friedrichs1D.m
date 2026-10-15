function flux = Lax_Friedrichs1D(x_Abl,fun_var,alpha,dx)

flux = WENO3rec_x(x_Abl,fun_var,dx,alpha);

end


