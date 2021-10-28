function [SIR_table,SIR_table_latex,time_table,time_table_latex]=results_print(Y,SIR,time_elapsed,Structure_data)
% function for pirintig the results
%[SIR_table,time_table]=results_print(Y,SIR,time_elapsed,Structure_data)

Methods=Structure_data.Methods;
labels_methods_NMF=Structure_data.labels;

%Print images

if ~isempty(Y)
    y=length(Y)+1;
    if size(Structure_data.X,3)<=3
        
        figure('Name','Original and reconstructed images')
        
        %Show original image
        
        subplot((ceil((y)/3)),3,1)
        imshow(uint8(Structure_data.X))
        title('Original image')
        
        
        for i=2:y
            subplot((ceil(y/3)),3,i)
            imshow(uint8(Y{i-1}))
            title(labels_methods_NMF{Methods(i-1)})
            
        end
        
    else
        
        figure('Name','Original and reconstructed images')
        
        %Show original image
        
        subplot((ceil((y)/3)),3,1)
        imshow(uint8(Structure_data.X(:,:,round(size(Structure_data.X,3)/2))))
        title('Original image')
        
        
        for i=2:y
            subplot((ceil(y/3)),3,i)
            z=Y{i-1};
            imshow(uint8(z(:,:,round(size(Structure_data.X,3)/2))))
            title(labels_methods_NMF{Methods(i-1)})
            
        end
        
    end
end

%Print SIRs results plots and table


if Structure_data.MC==1
    
    disp('To show boxplot more iterations needed!')
    
else
    figure('Name','SIR results for selected method')
    boxplot(SIR')
    xlabel('Algorithms','FontName','Times New Roman','FontSize',24)
    ylabel('SIR [dB]','FontName','Times New Roman','FontSize',24)
    set(gca,'FontName','Times New Roman','FontSize',20)
    set(gcf,'Color',[1 1 1])
    set(gca,'XTickLabel',labels_methods_NMF(Methods))
end

SIR_mean=mean(SIR,2);
SIR_std=std(SIR,0,2);

SIR_table = table(SIR_mean ,SIR_std ,'RowNames',labels_methods_NMF(Methods));

SIR_table_latex=latex((sym(vpa(round(table2array(SIR_table),2)))));

time_mean=mean(time_elapsed,2);
time_std=std(time_elapsed,0,2);

time_table = table(time_mean ,time_std ,'RowNames',labels_methods_NMF(Methods));

time_table_latex=latex((sym(vpa(round(table2array(time_table),2)))));



end
