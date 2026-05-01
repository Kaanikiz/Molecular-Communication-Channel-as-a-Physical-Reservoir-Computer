function generateSimConfig_MCConfirmation(N_i, input_time_points, txMode, r_tx, smoldyn_visualization_on)

    filename = 'MCpointTXsphRX_generated.txt';
    fileID = fopen(filename, 'w');

    fprintf(fileID, '# Monte Carlo confirmation for point approximation\n\n');

    if smoldyn_visualization_on == 1
        fprintf(fileID, 'graphics opengl\n');
        fprintf(fileID, 'graphic_iter 1\n\n');
    end

    fprintf(fileID, 'dim 3\n');
    fprintf(fileID, 'boundaries x -BOUNDARYLENGTH BOUNDARYLENGTH\n');
    fprintf(fileID, 'boundaries y -BOUNDARYLENGTH BOUNDARYLENGTH\n');
    fprintf(fileID, 'boundaries z -BOUNDARYLENGTH BOUNDARYLENGTH\n\n');

    fprintf(fileID, 'species Receptor Messenger ReceptorActive\n\n');

    fprintf(fileID, 'difc Receptor DIFFRECEPTOR\n');
    fprintf(fileID, 'difc ReceptorActive DIFFRECEPTORACT\n');
    fprintf(fileID, 'difc Messenger DIFFMESSENGER\n\n');

    if smoldyn_visualization_on == 1
        fprintf(fileID, 'color Receptor(all) maroon\n');
        fprintf(fileID, 'color Messenger(all) orange\n');
        fprintf(fileID, 'color ReceptorActive(up) fuchsia\n\n');

        fprintf(fileID, 'display_size all(all) 3\n');
        fprintf(fileID, 'display_size *Active(all) 7\n\n');
    end

    fprintf(fileID, 'time_start START_TIME\n');
    fprintf(fileID, 'time_stop STOP_TIME\n');
    fprintf(fileID, 'time_step TIME_STEP\n\n');

    fprintf(fileID, 'frame_thickness 0\n\n');

    fprintf(fileID, 'start_surface outsides\n');
    fprintf(fileID, '  unbounded_emitter front Messenger 1 TXPOSITION 0 0\n');
    fprintf(fileID, '  panel rect +x -BOUNDARYLENGTH -BOUNDARYLENGTH -BOUNDARYLENGTH BOUNDARYLENGTH*2 BOUNDARYLENGTH*2\n');
    fprintf(fileID, '  panel rect -x BOUNDARYLENGTH -BOUNDARYLENGTH -BOUNDARYLENGTH BOUNDARYLENGTH*2 BOUNDARYLENGTH*2\n');
    fprintf(fileID, '  panel rect +y -BOUNDARYLENGTH -BOUNDARYLENGTH -BOUNDARYLENGTH BOUNDARYLENGTH*2 BOUNDARYLENGTH*2\n');
    fprintf(fileID, '  panel rect -y -BOUNDARYLENGTH BOUNDARYLENGTH -BOUNDARYLENGTH BOUNDARYLENGTH*2 BOUNDARYLENGTH*2\n');
    fprintf(fileID, '  panel rect +z -BOUNDARYLENGTH -BOUNDARYLENGTH -BOUNDARYLENGTH BOUNDARYLENGTH*2 BOUNDARYLENGTH*2\n');
    fprintf(fileID, '  panel rect -z -BOUNDARYLENGTH -BOUNDARYLENGTH BOUNDARYLENGTH BOUNDARYLENGTH*2 BOUNDARYLENGTH*2\n');

    if smoldyn_visualization_on == 1
        fprintf(fileID, '  polygon both edge\n');
        fprintf(fileID, '  color both green\n');
    end

    fprintf(fileID, 'end_surface\n\n');

    fprintf(fileID, 'start_surface receiver\n');
    fprintf(fileID, 'action both all reflect\n');

    if smoldyn_visualization_on == 1
        fprintf(fileID, 'color both grey\n');
        fprintf(fileID, 'polygon both edge\n');
    end

    fprintf(fileID, 'panel sphere TXPOSITION+TXRXDISTANCE 0 0 RXRADIUS 30 30\n');
    fprintf(fileID, 'end_surface\n\n');

    fprintf(fileID, 'start_surface receiverspace\n');
    fprintf(fileID, 'action both all transmit\n');

    if smoldyn_visualization_on == 1
        fprintf(fileID, 'polygon both edge\n');
    end

    fprintf(fileID, 'panel sphere TXPOSITION+TXRXDISTANCE 0 0 RXRADIUS+RXRECEPTIONSPACETHICKNESS 30 30\n');
    fprintf(fileID, 'end_surface\n\n');

    fprintf(fileID, 'reaction binding Messenger(fsoln) + Receptor(up) -> ReceptorActive(up) KON\n');
    fprintf(fileID, 'reaction unbinding ReceptorActive(up) -> Messenger(fsoln) + Receptor(up) KOFF\n\n');

    fprintf(fileID, 'start_compartment receiver\n');
    fprintf(fileID, 'surface receiver\n');
    fprintf(fileID, 'point TXPOSITION+TXRXDISTANCE 0 0\n');
    fprintf(fileID, 'end_compartment\n\n');

    fprintf(fileID, 'start_compartment outreceiver\n');
    fprintf(fileID, 'surface receiverspace\n');
    fprintf(fileID, 'point TXPOSITION+TXRXDISTANCE 0 0\n');
    fprintf(fileID, 'end_compartment\n\n');

    fprintf(fileID, 'start_compartment receiverspace\n');
    fprintf(fileID, 'compartment equal outreceiver\n');
    fprintf(fileID, 'compartment andnot receiver\n');
    fprintf(fileID, 'end_compartment\n\n');

    fprintf(fileID, 'surface_mol numRECEPTOR Receptor(up) receiver all all\n\n');

    fprintf(fileID, 'reaction_serialnum binding r2\n');
    fprintf(fileID, 'reaction_serialnum unbinding new + r1\n');
    fprintf(fileID, 'product_placement unbinding pgemmax 0.2\n\n');

    fprintf(fileID, 'output_files receivedsignal_varNo_variableNo_varValue_variableValue_iter_index.txt\n');
    fprintf(fileID, 'output_files allmolecules_varNo_variableNo_varValue_variableValue_iter_index.txt\n');
    fprintf(fileID, 'cmd b overwrite receivedsignal_varNo_variableNo_varValue_variableValue_iter_index.txt\n');
    fprintf(fileID, 'cmd b overwrite allmolecules_varNo_variableNo_varValue_variableValue_iter_index.txt\n\n');

    %% Release molecules

    for i = 1:length(N_i)

        moleculeCount = round(N_i(i));

        if moleculeCount <= 0
            continue
        end

        tPulse = input_time_points(i);

        if txMode == "point"

            fprintf(fileID, 'cmd @ %.2f set mol %d Messenger TXPOSITION 0 0\n', ...
                tPulse, moleculeCount);

        elseif txMode == "finite"

            xyz = randomPointsInSphere(moleculeCount, r_tx);

            for m = 1:moleculeCount
                fprintf(fileID, 'cmd @ %.2f set mol 1 Messenger TXPOSITION%+.8f %.8f %.8f\n', ...
                    tPulse, xyz(m,1), xyz(m,2), xyz(m,3));
            end

        else
            error('Unknown txMode.')
        end
    end

    fprintf(fileID, '\ncmd b molcountheader allmolecules_varNo_variableNo_varValue_variableValue_iter_index.txt\n');
    fprintf(fileID, 'cmd N SAMPLING_PERIOD molcountincmpt2 receiverspace soln receivedsignal_varNo_variableNo_varValue_variableValue_iter_index.txt\n');
    fprintf(fileID, 'cmd N SAMPLING_PERIOD molcount allmolecules_varNo_variableNo_varValue_variableValue_iter_index.txt\n\n');

    fprintf(fileID, 'text_display time Receptor(up) Messenger ReceptorActive(up)\n\n');

    fprintf(fileID, 'end_file\n');

    fclose(fileID);
    fclose("all");

    fprintf('Smoldyn config generated: %s, txMode = %s, r_tx = %.2f, visualization = %d\n', ...
        filename, txMode, r_tx, smoldyn_visualization_on)
end
