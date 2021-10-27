function D = initializeDataStructure
% INITIALIZEDATASTRUCTURE generate a structure array to store the temp data
% used during the analysis such as the data processing images.
% Initialization speed up the computation as matlab use strcture array as
% pointers and does not generate a new variable at each step of the for
% loop.

% Initialise data structure  
D = struct;

end