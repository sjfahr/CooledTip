% This is the superscript that has the paths for all of the patient data.

% Paths.
setenv ( 'PATH22' , '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/000');
path22 = getenv ( 'PATH22' );

[metric] = simple_model_obj_fxn ( path22 );