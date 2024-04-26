% WFSORTER() - Alfred Manual Waveform Sorter - Matlab app/graphical 
%   user interface for manual sorting of multivariate data. Designed for 
%   spike detection and sorting, but useful for manual exclusion and 
%   clustering of any data type.
% 
%   Usage:
%       wfapp = wfsorter(X,thr,L,pr,align);
%   
%   Inputs:
%       X = data(observations,variables). If X is a vector (signal),
%           performs threshold-based waveform detection. If X is array
%           (waveforms), goes directly to sorting and uses "thr" and "L" 
%           inputs differently.
%       thr = threshold for spike detection. [default: Quiroga et al., 2004].
%           If thr = 0, creates array of duplicate X vector for sorting.
%           Tip: change 'X' to 'Time' and use for sorting actual signal.
%           If X is array, "thr" is used as samples indices or timestamps.
%       L = length (samples) of waveforms. Some following clustering
%           analysis needs it to be a power of two. [default: 32].
%           if X is array, "L" is used as cluster ID vector
%       pr = refractory period. Minimum samples for next spike detection.
%           [default: L].
%       align = method to align waveforms. {'thr','mod','max','min'}
%           [default: thr]
% 
%   Outputs:
%       wfapp = handle of app.
%       sortresults = structure array with sorting results
%           idx_valid = indices of waveforms not excluded.
%           idx_invalid = indices of waveforms excluded.
%           idx_units = numerical cluster ID vector of valid waveforms.
%           sortstats = statistical estimates of sorting quality.
% 
%   Buttons and callbacks:
%       Restart = Restart sorting with original data, with no units. You
%           can define a new amplitude threshold for waveform detection 
%           and waveforms alignment method {'thr','max','min','mod'}. 
% 
%       Select Unit[ ] = Select unit to select or remove waveforms.
% 
%       Select WF 
%           (dot) = Draw any shape on map to select inside waveforms for 
%               current unit
%           (line) = Draw line on plot to selecct intersecting waveforms 
%               for current unit
% 
%       Remove WF 
%           (dot) = Draw any shape on map to remove inside waveforms from 
%               current unit
%       	(line) = Draw line on plot to remove intersecting waveforms 
%               from current unit
%           (unit) = Remove all waveforms of current unit
%  
%       Exclude WF 
%           (dot) = Select dots to exclude on map.
%           (line) = Draw line to exclude on graph.
%           (unit) = Exclude all waveforms from current unit.
% 
%       Merge units = Merge two units
%           Units[ ][ ] = Units to merge. Merge to the lowest cluster ID value.
% 
%       Check
%           Visualize WF = Click next to waveform of interest on map, then
%               press enter. Will show waveform in thick yellow in the waveform
%               graph.
%           Average WF = Plot average waveform and standard deviation for
%               each unit.
%           Visualize dist. = Plot 2D histogram of current
%              map data including all units. Tip: if dim.1 is 'Time' it is
%              equivalent to firing counts.
%           Visualize ISI = Plot Inter-spike interval for each unit. Does
%               not perform if X is array.
%           Check in signal = Plot original signal and timestamps of spikes
%               for each unit. Does not perform if X is array.
%       
%       Show(Item) = Selection of what waveforms to show
%           All units = Show all units and unsorted
%           Single unit = Show current unit and unsorted
%           None = Show unsorted
% 
%       Recalculate map = Performs dimensionality reduction or select
%           specific variables and plots two-dimensional map.
%           ( )PCA = Uses PC scores. [default]
%           ( )Var# = Uses variable number#. 
%           X[ ] = Uses PC# or Var# for dim. 1.
%           Y[ ] = Uses PC# or Var# for dim. 2.
%           ( )t-SNE [ = Calculates t-SNE with [] number of PC scores.
%           [ ]X (Item) = Select specific variable for dim. 1
%           [ ]Y (Item) = Select specific variable for dim. 2
%               Index = Index of observations
%               Time = Index of observations of detected spikes. If X was
%                   already concatenated it is the same as Index
%               Peak = Max. ampl. of waveform
%               Valley = Min. ampl. of waveform
%               Energy = Sum of squares of waveform
%               Template = Euclidean distance from the average of the
%                   current selected unit
% 
%       Select outliers = Identifies outliers and add to current unit
%           ( )line = Identifies in any variable of waveforms
%           ( )dot = Identifies in any dimensionn of current feature space
%           Method (Item) = Method for detecting outliers. Uses 'isoutlier'
%               median = 3 scaled MAD
%               mean =  3 standard deviations from the mean
%               quartiles =  tukey fences with 1.5*(IQR)
%               density = density-based outlier detection based on map.
%                   Calculates 100-bin 2D histogram, then locates the knee
%                   point using 'findchangepts' of proportions across bins.
%                   The "threshold" value is the ratio to the knee point.
%                   [default = 1/100].
%       
%       Clustering
%           Import clustering = Import clustering index variable from
%               worspace. Type the name inside the text box and push
%               button. Do not use ''. Allows cluster ID vector or logical
%               array (observations,clusters). Must have the same number of
%               rows (obervations) than waveforms.
% 
%           Automatic cluster = Performs automatic clustering.
%               Method (Items): 
%                   GMM = Performs automatic Gaussian Mixture Model
%                       clustering. Does 100 replicates and selects waveforms
%                       with higher than 0.95 posterior probability.
%                   k-means = Performs automatic k-means clustering. Does
%                       100 replicates with 1000 maximum iterations.
%                   linkage = Performs automatic agglomerative hierarchical
%                       clustering. 
%                   dbscan = Performs automatic density-based clustering.
%                       Allows manual selection of epsilon on k graph.
%                   wave_clus = Performs automatic clustering using the
%                       wave_clus toolbox. Must have it installed.
% 
%               k[ ] = Number of clusters
%               PC[ ] = Number of PC scores. If PC = 0, uses the current
%                   map data
%     
%       Statistics = displays measures of sorting quality statistics. Uses
%           'sortstatsfcn.m'. Measures are returned within 'sortresults' as 
%           field 'sortstats' presenting:
%               Fmanova = F-value of one-way multivariate analysis of
%                   variance (MANOVA)
%               pvalue = p-value of significant dimension of MANOVA
%               J3 = J3 index
%               PseudoF = Pseudo-F value
%               DBindex = Davies-Bouldin index
%               Dunn = Dunn indexz
% 
%           all( ) = Includes all waveforms for calculations
%           sorted( ) = Ignores "unsorted" for calculations
% 
%       End sorting = Closes sorter app and returns the results structure
%           variables on workspace:
% 
%           wfresults
%               thr = threshdold used for spike detection
%               idx_wf = sample index of detected spikes
%               wf = waveforms
%           sortresults
%               idx_valid = index of valid (non-excluded) spikes based on
%                   the initially concatenated waveforms
%               idx_units = cluster ID vector for each sorted units.
%               sortstats = sorting quality statistics
% 
%   Author: Danilo Benette Marques, 2023
% 