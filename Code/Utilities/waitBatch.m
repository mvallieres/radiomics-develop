function waitBatch(pathCheck,time,nBatch)
% -------------------------------------------------------------------------
% function waitBatch(pathCheck,time,nBatch)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function implements a waiting loop ensuring that all the
% computations from all parallel batch are done.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathCheck: Full path to the directory where the 'batch1_end',
%   'batch2_end', etc. (parallel checkpoints) are saved.
% 2. time: Number of seconds to wait before checking if parallel 
%          computations are done.
%          --> Ex: 60
% 3. nBatch: Number of parallel batch.
%            --> Ex: 8
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: April 2017
% -------------------------------------------------------------------------
% DISCLAIMER:
% "I'm not a programmer, I'm just a scientist doing stuff!"
% -------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics-develop/>, 
% a private repository dedicated to the development of programming code for
% new radiomics applications.
% --> Copyright (C) 2017  Martin Vallieres
%     All rights reserved.
%
% This file is written on the basis of a scientific collaboration for the 
% "radiomics-develop" team.
%
% By using this file, all members of the team acknowledge that it is to be 
% kept private until public release. Other scientists willing to join the 
% "radiomics-develop" team is however highly encouraged. Please contact 
% Martin Vallieres for this matter.
% -------------------------------------------------------------------------

startpath = pwd;

cd(pathCheck)
if nBatch
    while 1
        pause(time);
        check = zeros(nBatch,1);
        for i = 1:nBatch
            check(i) = exist(['batch',num2str(i),'_end'],'file');
        end
        if sum(check) == nBatch*2
            break
        end
    end
end

cd(startpath)
end