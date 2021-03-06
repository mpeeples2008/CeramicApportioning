# CeramicApportioning

This is an R implementation of the procedure outlined by John M. Roberts, Jr. et al. (2012) "A method for chronologically apportioning of ceramic assemblages" published in the Journal of Archaeological Science

Roberts, John M., Jr., Barbara J. Mills, Jeffery J. Clark, W. Randall Haas, Jr., Deborah L. Huntley, and Meaghan A. Trowbridge
2012  A method for chronologically apportioning of ceramic assemblages. _Journal of Archaeological Science_ 39(5):1513-1520.

Abstract:

Artifact assemblages from long-inhabited sites may include ceramic types and wares from multiple time periods, making temporal comparisons between sites difficult. This is especially problematic in macro-regional data sets compiled from multiple sources with varying degrees of chronological control. We present a method for chronological apportioning of ceramic assemblages that considers site occupation dates, ceramic production dates, and popularity distribution curves. The chronological apportioning can also be adjusted to take into account different population sizes during the site occupation span. Our method is illustrated with ceramic data from late prehispanic sites in the San Pedro Valley and Tonto Basin, Arizona, U.S.A., compiled as part of the Southwest Social Networks Project. The accuracy of the apportioning method is evaluated by comparing apportioned assemblages with those from nearby contemporaneous single component sites.

The main script to reproduce the apportioning analysis is in the "Apportion.R" file and the csv file "preapportion.csv" provides sample data from the San Pedro valley of Arizona that can be used to apply this method. The Rmd and associated html output provides additional details guiding you through the use of the technique. I have also now added instructions for the iterative proportional fitting procedure and added another scirpt that performs that analysis on the output generated from the Apportion.R function.
