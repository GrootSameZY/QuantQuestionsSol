# QuantQuestionsSolutions
## Q1
### (a)
Dan's salary is greater than his manager Sally.
Sally's salary is greater than her manager Jane.
Joe's salary is greater than her manager Phil.
### (b)
If there is no subsidy drawn to their managers, the average salary should be 
John:300
Mike:200
Joe:600
Dan:600
Equal-weighted Avg = 425;
Weighted Avg = 500.

## Q2 
I think it depends on the type of the variable. For instance, the ```const``` variable in ```C++``` must be defined and initialized. Otherwise, there is the compilation error. Most undefined variables can be detected by the compiler. For those cannot be detected, e.g. the incomplete initialization of ```pointer``` and ```int```, can be detected by checking their address and values (i.e. whether the pointer points to a void address, whether the int value is random and extremely positive/negative, or whether there is a memory leak at the runtime).

## Q3
See ```mainQ3.cpp```.

![image](https://user-images.githubusercontent.com/89716697/131442265-848eed66-2a67-4b16-8aad-53a85be8a85a.png)

## Q4
The daily adjusted prices of the 7 stocks from 01/01/2016 to 12/31/2016 are accessed from Yahoo Finance. Calculate the log daily returns of each stock and the weighted return of the portfolio.
### (a)
Used historical simulation, （+: loss, -: return）
|  	| 95% VaR 	| 95% CVaR 	|
|---	|---	|---	|
| Historical Simulation 	| 0.01641343 	| 0.02300513 	|

### (b)
Use ```density``` to find the mean ```mean = -0.00057 (return)```and standard deviation ```std = 0.00950``` of the portfolio return.
The p.d.f below shows left-skewed and left heavy-tailed, which means the greater possibility of the occurence of huge loss.
![image](https://user-images.githubusercontent.com/89716697/131369628-41bf533a-82a9-41e2-b73d-030e44c6a470.png)

The correlation matrix is shown below.
| Corr 	| AAPL 	| IBM 	| GOOG 	| BP 	| XOM 	| COST 	| GS 	|
|---	|---	|---	|---	|---	|---	|---	|---	|
| AAPL 	| 1.0000000 	|  	|  	|  	|  	|  	|  	|
| IBM 	| 0.3063460  	| 1.0000000 	|  	|  	|  	|  	|  	|
| GOOG 	| 0.4727447  	| 0.3401606  	| 1.0000000 	|  	|  	|  	|  	|
| BP 	| 0.2907049  	| 0.4129107  	| 0.2070709  	| 1.0000000 	|  	|  	|  	|
| XOM 	| 0.2716068  	| 0.4342248  	| 0.1924154  	| 0.6814065  	| 1.0000000 	|  	|  	|
| COST 	| 0.3164456 	| 0.1754480 	| 0.3552544 	| 0.1309359 	| 0.1622474 	| 1.0000000 	|  	|
| GS 	| 0.3620238 	| 0.4276626 	| 0.3095373 	| 0.5398358 	| 0.4167982 	| 0.2491354 	| 1.0000000 	|
p.s: I don't have too much idea to solve it by its correlation matrix, so I just leave the correlation matrix here.

Use t distribution (df = 251-2-1=248,<a href="https://www.codecogs.com/eqnedit.php?latex=\lambda&space;=&space;\frac{std}{\sqrt{\frac{\lambda}{\lambda-2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda&space;=&space;\frac{std}{\sqrt{\frac{\lambda}{\lambda-2}}}" title="\lambda = \frac{std}{\sqrt{\frac{\lambda}{\lambda-2}}}" /></a>), （+: loss, -: return）
|  	| 95% VaR 	| 95% CVaR 	|
|---	|---	|---	|
|t distribution	| 0.01504312 	| 0.02019171 	|

### (c)
We aim to maximize the weighted daily return,

<a href="https://www.codecogs.com/eqnedit.php?latex=\text{Maximize}&space;\sum_{i=1}^{7}w_i\cdot&space;r_i,\\&space;\text{subject&space;to}&space;\sum_{i=1}^{7}w_i=1." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\text{Maximize}&space;\sum_{i=1}^{7}w_i\cdot&space;r_i,\\&space;\text{subject&space;to}&space;\sum_{i=1}^{7}w_i=1." title="\text{Maximize} \sum_{i=1}^{7}w_i\cdot r_i,\\ \text{subject to} \sum_{i=1}^{7}w_i=1." /></a>

See ```Q4.r```. Use ```lpSolve```. The result seems not correct since no negative weights appear, and the returns of each day seem not optimal.

## Q5

## Q6
I am not clear about how to extract the different format of date. I just have some ideas about splitting the text into categories which was done previously. Part of the code is shown in ```mainQ6.cpp```.
