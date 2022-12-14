Flexible MSTPP
================

[![Status](https://www.repostatus.org/badges/latest/active.svg)](https://github.com/ForeStats/flexible-msttp-football)
[![GitHub](https://img.shields.io/github/license/ForeStats/flexible-msttp-football)](https://opensource.org/licenses/GPL-3.0)
[![arXiv](https://img.shields.io/badge/arXiv-2103.04647-blue.svg)](https://arxiv.org/abs/2103.04647)

## Introduction

This GitHub repository provides all the computer code used to carry out the analyses presented 
in our paper "Flexible marked spatio-temporal point processes with applications to event sequences 
from association football". The arXiv link to the paper is [here](https://arxiv.org/abs/2103.04647).

The raw data that motivated our work has been provided by Stratagem Technologies Ltd and consists 
of all touch-ball events for the 2013/14 season of the English Premier League. The raw data cannot 
be disclosed as the authors do not have the license to do so. However, the code can be used to apply 
our methods to similar data sets, like the publicly-available 2020/21 FA Women's Super League Data 
provided by StatsBomb Inc. available [here](https://github.com/statsbomb/open-data).

## Data format

Our methods provide a general framework to analyse marked spatio-temporal event sequences where 
events are charaterised by their occurrence times (t\_i), locations (z\_i), and marks (m\_i). In 
applications like football, we also have the team information and the game period in addition to 
the core dimensions. The expected data structure is given below

| id | period | team\_id | time (t\_i) | zone (z\_i) | mark (m\_i) |
| --- | --- | --- | --- | --- | --- |
|  101 | 1 | 1 | 0 | 2 | Away\_Pass\_S |
|  101 | 1 | 1 | 1 | 2 | Away\_Pass\_U |
|  101 | 1 | 2 | 3 | 1 | Home\_Clear |
|  101 | 1 | 1 | 6 | 3 | Away\_Win | 
|  101 | 1 | 1 | 8 | 3 | Away\_Pass\_S |
|  101 | 1 | 1 | 15 | 2 | Away\_Pass\_S |
|  101 | 1 | 1 | 16 | 1 | Away\_Pass\_U |
|  101 | 1 | 2 | 19 | 1 | Home\_Out\_Throw |
