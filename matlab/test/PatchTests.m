classdef PatchTests < matlab.unittest.TestCase
% tc = PatchTests; run( tc, 'testGrabPatches' );
% tc = PatchTests; run( tc, 'testReconIm' );

    properties
        im_grX;
        im_grY;
    end
    
    methods(TestMethodSetup)
        function createInstance(tc)
            [grX, grY] = ndgrid( 1:9, 1:9, 1:9 );
            tc.im_grX = grX;
            tc.im_grY = grY;
        end

    end
    
    methods(TestMethodSetup)
        function destroyInstance(tc)
            tc.im_grX  = [];
            tc.im_grY = [];
        end
    end
    
    methods (Test)
        
        function testGrabPatches( tc )
            N = 5;
            sz = [3 3 3];
            pRand = grabPatchesSimple( tc.im_grX, sz, N );
            
            tc.verifySize( pRand, [N prod(sz)], ...
                        'patch matrix size')
                    
        end
        
        function testReconIm( tc )
            sz = [3 3 3];
            [cx,cy,cz] = ndgrid( 2:3:8, 2:3:8, 2:3:8 );
            coords = {cx(:),cy(:),cz(:)};
            pIm = grabPatchesSimple( tc.im_grX, sz, [], coords );
            imre = imageFromPatches( pIm, size(tc.im_grX), sz, coords );
            
            tc.verifyEqual( imre, tc.im_grX, ...
                            'check that reconstruction works');
        end
    end
    
end
