% simple isequal(a,b) function to compare results from an individual
% packages with results from "merged" functions
% 
% control data ... results from an individual packages (sofi2d, sofi3d,
%                  and multicolor_sofi)
% tested data .... results from the "utils" after merging functions
%
clc, clear
addpath('experiments');
addpath(genpath('utils'));

%% sofi2d
clear
input_path = 'data\sofi2d\test_main\20210615';
test_main_script(input_path, 'experiment_sofi2d')

%% sofi3d
clear
input_path = 'data\sofi3d\test_main\20210615\';
test_main_script(input_path, 'experiment_sofi3d')

%% multicolor sofi
clear
input_path = 'data\sofimc\test_main\20210615\';
test_main_script(input_path, 'experiment_sofimc')

%% biplane sofi
clear
input_path = 'data\sofibiplane\test_main\2021618\';
test_main_script(input_path, 'experiment_sofibp')
