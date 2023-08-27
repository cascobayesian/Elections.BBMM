# Identifying voting blocs from referendum data

The Elections.BBMM package provides a set of tools for applying the Beta-binomial mixture model to longitudinal referendum data. The purpose of this model is to identify voting blocs -- long-term patterns of vote similarity -- across a set of political districts (e.g. precincts, municipalities, counties). The model results can be used to infer:
* the distribution on the number of voting blocs
* the proportion of each voting bloc in each district
* the mean and standard deviaation for each question for each voting bloc
* identifications of unusual questions from the perspective of the model
* an information-theoretic visualization of which questions drive distinctions
* the topological structure of the voting blocs

The functions in this package let you execute all of the steps necessary to complete these analyses, including the BD-MCMC algorithm that generate samples from the posterior distribution. It also includes the scripts necessary to generate the various posterior visualizations that might be of interest. 

This form of data analysis might appeal to political scientists, sociologists, or historians. Why might you want to perform this analysis:
* Identify political voting blocs and their presence in a set of communities;
* Use voting blocs to identify long-term political/cultural structures and dynamics;
* Identify questions that provide strong distinction between voting blocs; and,
* Identify questions that run counter to the inferred voting blocs (questions that behave strangely relative to the model's expectations).

## A NOTE ON INSTALLATION
- The R package here should install without difficulty on most machines. However, as the packaage has a dependency with Rarmadillo, you may need to tweak the install flags to ensure that the gcc compile has access to the right paths. To aid with this process, I've provided a tar ball as well that should allow for an easier installation. 
