clear; 
k_r = 0.4;
final_measured = load('fermi_surface.mat','I_thetax_thetay','kx','ky');
standard_measured = load('fermi_surface_0mA.mat','I_thetax_thetay','kx','ky');
standard_field = load('standard_field.mat','X','Y','Z','Bx','By','Bz');
initial_guess = [0;0;0;0;0;0;0.025];
tic
parameters = optimized_parameters(standard_measured, final_measured, k_r, standard_field, initial_guess);
elapsed_time = toc;
save('optimized_p.mat','parameters','elapsed_time');