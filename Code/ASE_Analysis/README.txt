You've arrived at the final, and interesting part, of the analysis!
Assuming you've done the ASEReadCounter
    ASEReadCounter_OutputProcessing will take the ASEReadCounter output and prepare it
    for ASE score calculations, analyses and visualizations with the following scripts:

There are a few options here:
    You want to see what we specifically did with our dataset? --> ASE_Analysis.R
    You don't have genotype data? --> ASE_NoGeno.R
    You don't have prior genes of interest? --> ASE_NoTop.R
    You don't have groups to compare? --> ASE_NoGroup.R
    You do have all of that? --> ASE_Analysis.R
    You just want ASE analysis for your non-genotyped samples? --> ASE_Crude.R

All of these take the output from ASEReadCounter_OutputProcessing and analyse and visualize data!

PS: In case you're interested, some of the functions created in the ASE_Analysis script(s) 
    are in the Functions folder for when you want to do your own thing completely or just orientate!
