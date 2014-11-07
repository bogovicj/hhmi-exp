function [acc, errors, fpc, fnc, tpc, tnc ] = binaryClassifierEval( truth, pred )
% [acc, errors, fpc, fnc, tpc, tnc ] = binaryClassifierEval( truth, pred )

errors = ( pred - truth );

fpc    = nnz( errors ==  1 );
fnc    = nnz( errors == -1 );
tpc    = nnz( errors == 0 & truth == 1 );
tnc    = nnz( errors == 0 & truth == 0 );

acc = nnz( errors == 0 )./length(truth);
