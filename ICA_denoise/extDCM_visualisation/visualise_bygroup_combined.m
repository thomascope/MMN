function visualise_bygroup(dirname_DCM,diagnosis_list,source_names,conductances)
% A simple function for combining the extDCM circuit diagrams

dirname_circuit = [pwd '/circuit_diagrams'];

addpath(['/group/language/data/thomascope/MMN/ICA_denoise/extDCM_visualisation/ojwoodford-export_fig-216b30e'])

thisdir = pwd;

cleanupObj = onCleanup(@()cd(thisdir));
cd([dirname_DCM 'PEB_secondlevel'])

%Load any input DCM and PEB of PEBs as template
this_DCM = load(['PEB_H_' source_names{1} '_' conductances{1} '_' diagnosis_list{1} '.mat']);
this_DCM = this_DCM.DCM{1};

template_PEB = load(['PEB_H_' source_names{1} '_' conductances{1} '_Overall_combined.mat']);
template_PEB = template_PEB.PEB_Overall;

assert(all(template_PEB.M.X(:,1)==1),'The first column of the PEB of PEBs contrast should be all ones, check please.')

xw = 400; %X_width
yw = 600; %Y_width

for this_contrast = 1:size(template_PEB.M.X,2)
    this_source = 1;
    combined_diagram = figure;
    if this_contrast == 1
        this_title = ['All positives ' source_names{this_source}];
    else
        this_title = [diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' minus ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' ' source_names{this_source}];
    end
    [image_to_add, ~, image_alphas] = imread([thisdir '/circuit_diagrams/' this_title '.png']);
    imagearray = uint8(zeros(size(image_to_add,1), size(image_to_add,2)*length(source_names), 3));
    imagearray_alphas = uint8(zeros(size(image_to_add,1), size(image_to_add,2)*length(source_names), 1));
    imagearray(1:size(image_to_add,1),1:size(image_to_add,2),:) = image_to_add;
    imagearray_alphas(1:size(image_alphas,1),1:size(image_alphas,2),1) = image_alphas;
    for this_source = 2:length(source_names)
        if this_contrast == 1
            this_title = ['All positives ' source_names{this_source}];
        else
            this_title = [diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' minus ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' ' source_names{this_source}];
        end
        [image_to_add, ~, image_alphas] = imread([thisdir '/circuit_diagrams/' this_title '.png']);
        imagearray(1:size(image_to_add,1),1+(this_source-1)*size(image_to_add,2):this_source*size(image_to_add,2),:) = image_to_add;
        %size(image_to_add,2) %Bugchecking
        imagearray_alphas(1:size(image_alphas,1),1+(this_source-1)*size(image_alphas,2):this_source*size(image_alphas,2),1) = image_alphas;
        
    end
    f = imshow(imagearray);
    set(f, 'AlphaData', imagearray_alphas);
    if this_contrast == 1
        this_title = ['All positives combined'];
    else
        this_title = [diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' minus ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' combined'];
    end
    drawnow
    
    %eval(['export_fig ' thisdir '/circuit_diagrams/' this_title '.pdf -transparent'])
    %     saveas(combined_diagram,[thisdir '/circuit_diagrams/' this_title '.png']);
    saveas(combined_diagram,[thisdir '/circuit_diagrams/' this_title '.pdf']);
    eval(['export_fig ''' thisdir '/circuit_diagrams/' this_title '''.png -transparent'])
    close all
end
end