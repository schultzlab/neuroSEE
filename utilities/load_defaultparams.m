function params_def = load_defaultparams( params )

    params_def = load( 'default_params.mat' );
    
    % Remove irrelevant parameters 
    if strcmpi(params.methods.mcorr_method,'normcorre')
        params_def = rmfield(params_def.mcorr,'fftRigid');
        fields = {'df_prctile','df_medfilt1'};
        params_def.ROIsegment = rmfield(params_def.ROIsegment,fields);
    elseif strcmpi(params.methods.mcorr_method,'normcorre-r')
        fields = {'normcorre_nr','fftRigid'};
        params_def.mcorr = rmfield(params_def.mcorr,fields);
    elseif strcmpi(params.methods.mcorr_method,'normcorre-nr')
        fields = {'normcorre_r','fftRigid'};
        params_def.mcorr = rmfield(params_def.mcorr,fields);
    elseif fftRigid
        fields = {'normcorre_r','normcorre_nr'};
        params_def.mcorr = rmfield(params_def.mcorr,fields);
    end
    
    if ~params.methods.dofissa
        params_def = rmfield(params_def,'fissa');
    end
    
    params_def.methods = params.methods;
