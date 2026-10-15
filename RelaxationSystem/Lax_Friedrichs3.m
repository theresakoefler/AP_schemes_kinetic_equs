function flux = Lax_Friedrichs3(x_Abl,y_Abl,fun_var,alpha,dx,dy)
   

    flux = WENO3rec_x(x_Abl,fun_var,dx,alpha)+WENO3rec_y(y_Abl,fun_var,dy,alpha);
end