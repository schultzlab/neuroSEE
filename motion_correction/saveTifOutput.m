function saveTifOutput(out_g, out_r, shifts, col_shift, template, imG, imR, template_g, template_r, params_mcorr,...
                file, fname_mat_mcorr, fname_tif_gr_mcorr, fname_tif_red_mcorr)
    if nargin<14, fname_tif_red_mcorr = []; end

    % Save output
    mcorr_output.green = out_g;
    if ~isempty(out_r)
        mcorr_output.red = out_r;
    end
    mcorr_output.shifts = shifts;
    mcorr_output.col_shift = col_shift;
    mcorr_output.template = template;
    mcorr_output.params = params_mcorr;
    if ~isempty(template_g), mcorr_output.template_g = template_g; end
    if ~isempty(template_r), mcorr_output.template_r = template_r; end
    save(fname_mat_mcorr,'-struct','mcorr_output');

    % Save motion corrected tif images
    if isempty(reffile)
        prevstr = sprintf( '%s: Saving motion corrected tif images...\n', file );
    else
        prevstr = sprintf( '%s: Saving registered tif images...\n', file );
    end
    cprintf('Text',prevstr);
        writeTifStack( imG,fname_tif_gr_mcorr );
        if ~isempty(fname_tif_red_mcorr), writeTifStack( imR, fname_tif_red_mcorr ); end
    if isempty(reffile)
        str = sprintf( '%s: Motion corrected tif images saved\n', file );
    else
        str = sprintf( '%s: Registered tif images saved\n', file );
    end
    refreshdisp( str, prevstr );
end