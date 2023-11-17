# PLM-perturbations

As it says in the code itself, it comes from a paper (also attached) with some explanation in supplementary document (also attached). The code needs a dataset to run and that is also attached (DATA_PLM.csv).

In short, this is a “fitting” problem to describe lactation curves (of goats in this case). We have done it and even have the “curve fitting” tool online that is somehow like it. The main difference is that in this paper they introduce “perturbations” as part of the fitting problem. At this point what I would like to be able to do is to reproduce what the paper reports. Then later we could use to apply to the development of a tool.

## Notes

to run the PLM code we need to specify the number of maximum perturbations we want the model to fit for. I dug into the code little,  played with it, and understood it. 

'N' is directly controlled by the loop index and does not depend on other conditions in the code.
There's no condition in the loop that would prematurely terminate it or reduce the value of 'N' below 'NMAX'.
'N' will not terminate being less than 'NMAX'. It will exactly follow the sequence from 1 to 'NMAX'. In my understanding of the code nowhere the number od perturbations is being calculated by any function. 

In the word document too, the "Two-step algorithm for PLM fitting data of one lactation" also specifies that the perturbation number 'n' will range starting from 1 till 'N'. Also, in the algorithm diagram, the authors specify that the number of perturbations 'N' or the 'NMAX' (both are same) should be either be equal to 15 or can be less than 15. I anticipate that if we set N more than 15 the model would severely overfit. 

In the page 6 of the Microsoft word document that the authors provide, there is these two diagrams a) and b)- 

![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

In the diagrams a) and b) how is the fitting (the red line) changing when we change the RMSE (0.11 to 0.41) like what does it mean? How is dependent on the RMSE (kg)
