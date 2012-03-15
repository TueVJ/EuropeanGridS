% csolve  Solves a custom quadratic program very rapidly.
%
% [vars, status] = csolve(params, settings)
%
% solves the convex optimization problem
%
%   minimize(sum(max(-(Delta - K*F), 0)))
%   subject to
%     h_mns <= F
%     F <= h_pls
%
% with variables
%        F  44 x 1
%
% and parameters
%    Delta  27 x 1
%        K  27 x 44
%    h_mns  44 x 1
%    h_pls  44 x 1
%
% Note:
%   - Check status.converged, which will be 1 if optimization succeeded.
%   - You don't have to specify settings if you don't want to.
%   - To hide output, use settings.verbose = 0.
%   - To change iterations, use settings.max_iters = 20.
%   - You may wish to compare with cvxsolve to check the solver is correct.
%
% Specify params.Delta, ..., params.h_pls, then run
%   [vars, status] = csolve(params, settings)


% Produced by CVXGEN, 2012-03-15 05:13:46 -0700.
% CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2011 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: csolve.m.
% Description: Help file for the Matlab solver interface.
