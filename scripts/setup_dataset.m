vox_thresh_scale = 3;
vox_thresh_num = 1;

mets = {'pyr', 'ala', 'lac', 'urea'};

switch dataset
    % Normal
     case 1
        dirpath = '/data/vig1/C13/animal/Normal_Mice/FVB_Normal_100421/';
        imname = 't2fse_ax';
        exppath = {'experiment2-3d/', 'experiment1-steprep/'}; %steprep_TM5ms
        specname = {'3d', '3dsteprep'};
        pol = [25.8 33.2];
        Imets = [228 148 55];
       
    
     case 2
        dirpath = '/data/vig1/C13/animal/Mice_liver/MSL57_Apr26_10/';
        imname = 't2fse_ax';
        exppath = {'experiment2-3d/', 'experiment1-steprep/'};
        specname = {'3d', '3dsteprep'};
        pol = [31.9 31.9];
        Imets = [228 148 55 77];
       
     case 3
        dirpath = '/data/vig1/C13/animal/Normal_Mice/M013_100922/';
        imname = 't2fse_ax';
        exppath = {'experiment1-3dcsi/', 'experiment2-steprep/'};
        specname = {'3dcsi', '3dsteprep'};
        pol = [25 24.8];
        Imets = [228 148 55] -5;
    
      case 4
        dirpath = '/data/vig1/C13/animal/Normal_Mice/M013_100923/';
        imname = 't2fse_ax';
        exppath = {'exp3-3dcsi/', 'exp2-steprep/'};
        specname = {'3dcsi', 'steprep'};
        pol = [21.2 23.6];
        Imets = [228 148 55];
   
    
    %TRAMPs
    
    case 101
        dirpath = '/data/vig1/C13/animal/TRAMP/MS204_100325/';
        imname = 't2fse_ax';
        exppath = {'experiment2-3dCS/', 'experiment1-steprep_3dCS/'};
        specname = {'3d', '3d_steprep'};
        pol = [14.5 15.4];
        Imets = [235 158 61 84];
        
    case 102
        dirpath = '/data/vig1/C13/animal/TRAMP/MS206_100331/';
        imname = 't2fse_ax';
        exppath = {'experiment2-3d/', 'experiment1-3d_steprep/'};
        specname = {'3d', '3d_steprep'};
        pol = [21 23.7];
        Imets = [229 147 56];
        
     case 103
        dirpath = '/data/vig1/C13/animal/TRAMP/MS217_100625/';
        imname = 't2fse_ax';
        exppath = {'experiment2-3d/', 'experiment1-3dsteprep/'};
        specname = {'3dnoprep', '3dsteprep'};
        pol = [32.4 31.9];
        Imets = [223 144 50];
        Idisp = 9;  axzoom = [75 230 80 175];
        
      case 104
        dirpath = '/data/vig1/C13/animal/TRAMP/MS215_100813/';
        imname = 't2fse_ax';
        exppath = {'experiment2-csi/', 'experiment1-steprep/'};
        specname = {'csi', 'steprep'};
        pol = [23.1 28.2];
        Imets = [229 147 56];
        
      
        % Liver tumors
        
     case 201
        dirpath = '/data/vig1/C13/animal/Mice_liver/MSL58_Apr12_10/';
        imname = 't2fse_ax';
        exppath = {'experiment2-3d/', 'experiment1-ste_prep/'};
        specname = {'3d', '3dsteprep'};
        pol = [23.7 24.7];
        Imets = [228 148 55 77];
        Idisp = 9;  axzoom = [75 230 70 170];
       
     case 202
        dirpath = '/data/vig1/C13/animal/Mice_liver/MSL58_Apr13_10/';
        imname = 't2fse_ax';
        exppath = {'experiment2-3d/', 'experiment1-steprep/'};
        specname = {'3d_noprep', '3d_steprep'};
        pol = [21.4 24.8];
        Imets = [228 148 55 77];
       
       
        
end

Nmets = length(Imets);