function feats = dawmr_feat_comp( dm_obj_or_fn, xs, ys, zs, data_size, ...
                                 flag, cube_size )
% Usage:
%   feats = dawmr_feat_comp( dm_obj_or_fn, xs, ys, zs, data_size, flag, cube_size )
%
% 

if( isa( dm_obj_or_fn, 'dawmr' ))
    this =  dm_obj_or_fn
elseif( exist( dm_obj_or_fn, 'file' ))
    this_struct = load(obj_fn,'this');
    this = this_struct.this
end

recur_min = this.ds.get_recur_min_offset()
recur_max = this.ds.get_recur_max_offset()

for xx=xs(1):cube_size:xs(2)
    for yy=ys(1):cube_size:ys(2)
        for zz=zs(1):cube_size:zs(2)
            
            xx_end = min(xx+cube_size-1,xs(2));
            yy_end = min(yy+cube_size-1,ys(2));
            zz_end = min(zz+cube_size-1,zs(2));
            xi = xx:xx_end;
            yi = yy:yy_end;
            zi = zz:zz_end;
            
            [xn, yn, zn] = ndgrid(xi,yi,zi);
            xn = xn(:); yn = yn(:); zn = zn(:);
            inds = sub2ind(data_size, xn, yn, zn);
            
            
            for i=1:length(this.dds)
                dd = this.dds(i);
                res = dd.scaling;
                if(dd.downsampling==1)
                    res = [1 1 1];
                end
                inds
                inds = dd.filter_indices(inds, ...
                    ceil(data_size./res), ...
                    recur_min, recur_max);
                inds
            end
            if(isempty(inds))
                continue
            end
            
            [xn, yn, zn] = ind2sub(data_size, inds);
            
            feats = this.get_features(inds, xn, yn, zn, flag, 1);
            feats
            switch(this.svm_normalization)
                case 0
                case 1
                    assert(~isempty(this.sn1_max{1}), ...
                        'DAWMRLIB:AssertionFailed', ...
                        'normalization parameter not set');
                    feats = bsxfun(@rdivide, feats, this.sn1_max{1});
                case 2
                    assert(~isempty(this.sn2_min{1}) && ...
                        ~isempty(this.sn2_max{1}), ...
                        'DAWMRLIB:AssertionFailed', ...
                        'normalization parameter not set');
                    feats = 2*bsxfun(@rdivide, ...
                        bsxfun(@minus, feats, ...
                        this.sn2_min{1}), ...
                        this.sn2_max{1}) - 1;
                case 3
                    assert(~isempty(this.sn3_mn{1}) && ...
                        ~isempty(this.sn3_std{1}), ...
                        'DAWMRLIB:AssertionFailed', ...
                        'normalization parameter not set');
                    feats = bsxfun(@rdivide, ...
                        bsxfun(@minus, feats, ...
                        this.sn3_mn{1}), ...
                        this.sn3_std{1});
                otherwise
                    assert(0, 'DAWMRLIB:AssertionFailed', ...
                        'unknown svm normalization');
            end
            
        end
    end
end

end
