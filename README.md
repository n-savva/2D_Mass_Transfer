# 2D Droplet Dynamics on Heterogeneous surfaces with mass transfer

This repository contains an implementation of the asymptotic model in Python (v.3) derived in

* [Droplet motion on chemically heterogeneous substrates with mass transfer. I. Two-dimensional dynamics](https://arxiv.org/abs/2007.07004)

In this work, we develop a reduced model based on matched asymptotic expansions that allows us to simulate the dynamics of 2D droplets on heterogeneous surfaces as they undergo changes in their mass. Extended comparisons are offered with simulations to the long-wave governing equation and highlight the generally excellent agreement with the governing long-wave equation.

## Files

The provided python scripts reproduce the figures in the above mentioned paper, namedly accordingly and within their corresponding folders. Generating the figures relies on the data and codes in the directory `main`, and, whenever appropriate on PDE data, which was obtained using MATLAB scripts (not currently made openly available). The `main` folder includes

* `ODEdrop2D.py`, which implements the solver of the asymptotic model
* `pdeloader.py` is a python script that loads matlab files which contain PDE solutions.

## Acknowledgement

This work has received funding through the Engineering and Physical Sciences Research Council of the United Kingdom through Grant No. EP/P505453/1 and from the European 
Union’s Horizon 2020 research and innovation program under Grant No. 810660.
