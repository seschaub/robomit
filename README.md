# robomit

[![CRAN](https://www.r-pkg.org/badges/version/robomit)](https://cran.r-project.org/package=robomit)

## Summary 
In nonexperimental work, omitted variables bias can cause flawed research conclusions. **robomit** implements the recent framework by [Oster](#References) (2019) in R, which assesses the potential severity of the omitted variable bias for the research conclusion. For this assessment, we can estimate i) the bias-adjusted treatment effect or correlation <img src="https://render.githubusercontent.com/render/math?math=\beta^{*}"> and ii) the degree of selection on unobservables relative to observables that would be necessary to eliminate the result <img src="https://render.githubusercontent.com/render/math?math=\delta^{*}"> based on the Oster framework, using standard regression output. **robomit** implements an easy-to-use estimation of both of these variables. Additionally, **robomit** includes sensitivity analyses of these variables and their visualization.
Find introduction to the package [here](https://sites.google.com/view/sergeischaub/robomit). 

## Install package 
```
# from CRAN
install.packages("robomit")
```

## References
Oster. Unobservable selection and coefficient stability: Theory and evidence. *Journal of Business & Economic Statistics*, 37(2):187–204, 2019. URL https://doi.org/10.1080/07350015.2016.1227711.
