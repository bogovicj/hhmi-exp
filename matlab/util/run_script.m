function run_script(script_name, summary)
%RUN_SCRIPT   runs script with automated numbering and logging
%   RUN_SCRIPT(SCRIPT_NAME, SUMMARY) creates a new directory
%   under the EXPPATH directory with the next experiment number.
%   The script SCRIPT_NAME is copied into this folder, the
%   current mercurial changeset information is printed to file, 
%   and the matlab command output is logged to file.
%   The string SUMMARY is echoed to the prompt and captured by
%   the log file.
  
  full_name = which(script_name);
  assert(~isempty(full_name), 'JANLIB:AssertionFailed', ...
         sprintf('could not find script %s', script_name));
  
  [f,n] = get_next_exp_number(script_name);
  global SAVEPATH
  SAVEPATH = f;
  global SAVEPREFIX
  global DATDIR
  SAVEPREFIX = sprintf('exp%04d', n);
  
  [~,sn_name,sn_ext] = fileparts(full_name);
  
  % copy over file
  system(sprintf('cp %s %s/%s_%s%s', full_name, f, ...
                 SAVEPREFIX, sn_name, sn_ext));

  % log output
  diary(sprintf('%s/output.txt', f));

  fprintf('running %s: %s\n', SAVEPREFIX, script_name);
  fprintf('summary: %s\n\n', summary);
  % write out repo information
  %system(sprintf(['hg log -l 1 -f -R %s > %s/code_ver.txt; ' ...
  %                'hg log -l 1 -f -R %s'], ...
  %               DANVILLIBPATH, f, DANVILLIBPATH));
  %fprintf('\n');
  
  eval(script_name);
  
  diary off
                 
end

function [f,n] = get_next_exp_number(script_name)
  global EXPPATH
  d = sprintf('%s/saved_exp', EXPPATH);
  if(~exist(d,'dir'))
    system(sprintf('mkdir %s', d));
  end
  de = dir(sprintf('%s/exp*', d));
  if(isempty(de))
    n = 1;
  else
    n = max(arrayfun(@(x)(str2double(x.name(4:7))),de)) + 1;
  end
  f = sprintf('%s/exp%04d_%s', d, n, script_name);
  system(sprintf('mkdir %s', f));
end
