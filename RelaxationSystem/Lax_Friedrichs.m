function flux = Lax_Friedrichs(x_Abl,y_Abl,fun_var,alpha,dx,dy)
   

    flux = WENO5rec_x(x_Abl,fun_var,dx,alpha)+WENO5rec_y(y_Abl,fun_var,dy,alpha);
end

