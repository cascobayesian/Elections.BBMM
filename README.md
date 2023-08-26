# Elections.BBMM
The Elections.BBMM package provides a set of tools for applying the Beta-binomial mixture model to longitudinal referendum data. The purpose of this model is to identify voting blocs -- long-term patterns of vote similarity -- across a set of political districts (e.g. precincts, municipalities, counties). The model results can be used to infer:
* the distribution on the number of voting blocs
* the proportion of each voting bloc in each district
* the mean and standard deviaation for each question for each voting bloc
* identifications of unusual questions from the perspective of the model
* an information-theoretic visualization of which questions drive distinctions
* the topological structure of the voting blocs

The functions in this package let you execute all of the steps necessary to complete these analyses, including the BD-MCMC algorithm that generate samples from the posterior distribution. It also include the scripts necessary to generate the various posterior visualizations that might be of interest. 

This form of data analysis might appeal to political scientists, sociologists, or historians. Why might you want to perform this analysis:
* Identify political voting blocs and their presence in a set of communities;
* Observe and understand long-term political/cultural structures and dynamics;
* Identify questions that provide strong distinction between voting blocs; and,
* Identify questions that run counter to the inferred voting blocs (questions that behave strangely relative to the model's expectations).
