% This is the superscript that has the paths for all of the patient data.

% Paths.
setenv ( 'PATH22' , '/home/ignis/ISMRM_patients/Patient0002/000/I_O/');
path22 = getenv ( 'PATH22' );

cd (path22)
load index.txt

[metric] = simple_model_obj_fxn ( path22 );

index = index + 1;
csvwrite ('index.txt' , index);