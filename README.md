# extract the AP
To analysis the action potential, run the "code_FIcurve_extract_abf.m" to extract the AP.
The code will extract the ".abf" documents, rename the documents to "N*-*.abf" (for example, C1-1.abf)
"fidx=dir('*.abf')" is the path of your folder, change the "current folder" to your destination folder in MATLAB
modify the "StimAmp_DP" to [initial stimulate current: Incremental current: final stimulate current]
modify the "TimeStim" to [start time : end time]
run the code, you will get the ".mat" documents
# analysis the AP
run the "Code_FIcurve_statistic_abf.m"
modify the "Fidx{1}=dir('DATA_C*.mat')" to "Fidx{1}=dir('DATA_"the group 1 name"*.mat')"
modify the "Fidx{2}=dir('DATA_C*.mat')" to "Fidx{1}=dir('DATA_"the group 2 name"*.mat')"
modify the "GroupName" to match the "Fidx{1}" and "Fidx{2}"
"Num" and "NPlusMinus" identify which curve you want to analyze, the trial evoking "Num-NPlusMinus to Num+NPlusMinus" APs will be analyzed. For example, if the Num=10, the NPlusMinus=4, the trial evoking 6-14 APs will be analyzed.
"thAP" identify which AP in the trial you want to analyze. For example, if thAP=1, the first AP in the trial will be analyzed.
run the code, results will be established to "Results_Statistic.xls".
