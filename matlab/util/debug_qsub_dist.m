% function out = debug_qsub_dist( fn );

fn = '/nobackup/saalfeld/john/tmp/job_20140701T141608_599456705/task_9.mat';
load( fn );

out = run_obj_method_dist( bundle.args{1}, bundle.args{2}, bundle.args{3}, bundle.args{4}, ...
                             bundle.args{5:end});
