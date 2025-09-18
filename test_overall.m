clear; 
k_r = 0.4;
final_measured = load('fermi_surface_30mA.mat','I_thetax_thetay','kx','ky');
standard_measured = load('fermi_surface_0mA.mat','I_thetax_thetay','kx','ky');
initial_guess = [0;0;0;0;0;0;0.025];
calc_standard_field()
tic
parameters = optimized_parameters(standard_measured, final_measured, k_r, standard_field, initial_guess);
elapsed_time = toc;
save('optimized_p.mat','parameters','elapsed_time');