import os
import shutil

def test_launch_job():
    # Test case 1: Test with 'qe' program and input_label
    folder = '/path/to/folder'
    program = 'qe'
    jobscript_file = '/path/to/jobscript'
    input_label = 'input'
    
    launch_job(folder, program, jobscript_file, input_label)
    
    # Assert that the jobscript file is copied to the specified folder
    assert os.path.exists(f'{folder}/jobscript')
    
    # Assert that the jobscript is executed with the correct command
    expected_command = f'sbatch jobscript {input_label}.pwi >> {input_label}.pwo'
    assert os.system.call_count == 1
    assert os.system.call_args[0][0] == expected_command
    
    # Test case 2: Test with 'vasp' program and no input_label
    folder = '/path/to/folder'
    program = 'vasp'
    jobscript_file = '/path/to/jobscript'
    
    launch_job(folder, program, jobscript_file)
    
    # Assert that the jobscript file is copied to the specified folder
    assert os.path.exists(f'{folder}/jobscript')
    
    # Assert that the jobscript is executed with the correct command
    expected_command = 'sbatch jobscript'
    assert os.system.call_count == 2
    assert os.system.call_args[0][0] == expected_command

test_launch_job()