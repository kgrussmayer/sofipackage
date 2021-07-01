%%% test helpers


function test_main_script(input_path, main_script)
% Helper function to test and evaluate script given by main_script
% input_path  ...  path to input test data              [str]
% main-script ...  name of the main stript to be tested [str]
    % control data
    folder_expected = 'expected';
    % test output
    folder_result = 'result';
    % remove results
    output_path = [input_path, filesep , folder_result, filesep];
    if isdir(output_path)
        rmdir([input_path, filesep , folder_result, filesep], 's')
    end
    
    % run test experiment
    eval(main_script)

    files_expected = getnamesdir(...
        [input_path, filesep , folder_expected], '.tif', true);
    files_result = getnamesdir(...
        [input_path, filesep , folder_result], '.tif', true);

    F = length(files_expected);
    for f = 1:F
        file_expected = imread(files_expected{f});
        file_result = imread(files_result{f});
        
        [ssimval,ssimmap] = ssim(file_result,file_expected);
        
        if ssimval==1
            disp(['TEST PASSED - SSIM = 1.0000 ' files_result{f}])
        elseif ssimval>=9999e-4
            disp(['TEST PASSED - SSIM > 0.9999 ' files_result{f}])
        else
            disp(['TEST FAILED - SSIM < 0.9999 ' files_result{f}])
        end
    end

end
