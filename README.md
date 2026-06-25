# Introduction
I made this package about separation processes graphical calculation in Python to distribute freely for educational purpose (targeted for my friends and juniors in Department of Chemical Engineering, Chulalongkorn University). This package is free, light-duty, and white-boxed; therefore, you can run this in your potato computer, or even just a phone while able to know how it works visually. (not like some 20GB program which is heavy, expensive, and black-boxed)

Let's enjoy with this subject.

*-Supphawit Sripusitto (CU INTANIA107)*
> [!NOTE]
> I just have too much time on my hands, so I begin this project a week before final exam to the time during internship to boost my cognitive load.

> [!CAUTION]
> This package works with only linear equilibrium line (constant partitioning coefficient / equilibrium constant) 

# (Extremely) Brief guide
|Method|What do it do|
|--|:--|
|AbsorberOperation|Calculate the operating conditions from spec. & design|
|AbsorberDesign|Calculate number of stages from spec.|

# AbsorberOperation

## Syntax
~~~
import SepUnit

SepUnit.AbsorberOperation(X0, YN1, Y1, V, L, K, N, (report), (graph))
~~~
## Parameters

|Parameter|Data Type|Meaning|
|--|--|:--|
|X0 |float| Inlet liquid mole ratio|
|YN1 |float| Inlet gas mole ratio|
|Y1 |float| Outlet gas spec. (put negative value for recovery fraction)|
|V |float| Molar gas (w/o solute) flow rate|
|L |float| Molar liquid (w/o solute) flow rate (put neg|for times of minimum liquid flow rate)|
|K |float| Equilibrium constant / Partitioning coefficient|
|N |integer| Number of stages|
|report *(optional)*| boolean| show report or not|
|graph *(optional)*| boolean| show graph or not|

## Description

AdsorptionOperation works by guessing the outlet liquid mole ratio. Then, it calculate the actual outlet liquid mole ratio to compare with the guessing one. If the guessing one differ from the actual by less than the tolerance (1e-7 relative fraction), the program stop.

## Example

This code
~~~
import SepUnit

SepUnit.AbsorberOperation(0.0000067,0.0067,-0.67,1,-1.5,67,3)
~~~
will give this report
~~~
===== Calculation Report =====
Gas outlet mole ratio:        0.0018263607316187501
Liquid outlet mole ratio:     7.375502817242288e-05
Solute recovery in absorbent: 0.7274088460270522
===== End of the report =====
~~~
and this plot

# AbsorberDesign
## Syntax
~~~
import SepUnit

SepUnit.AbsorberDesign(X0, YN1, Y1, V, L, K)
SepUnit.AbsorberDesign(X0, YN1, Y1, V, L, K, (Nm), (report), (graph))
~~~
## Parameters

|Parameter|Data Type|Meaning|
|--|--|:--|
|X0 |float| Inlet liquid mole ratio|
|YN1 |float| Inlet gas mole ratio|
|Y1 |float| Outlet gas spec. (put negative value for recovery fraction)|
|V |float| Molar gas (w/o solute) flow rate|
|L |float| Molar liquid (w/o solute) flow rate (put neg|for times of minimum liquid flow rate)|
|K |float| Equilibrium constant / Partitioning coefficient|
|Nm *(optional)*|integer| Maximum number of stages|
|report *(optional)*| boolean| show report or not|
|graph *(optional)*| boolean| show graph or not|
## Description
AbsorberDesign is the original algorithm as McCabe-Thiele to design absorber with iconic staircase shape graph between operating and equilibrium line. The program run until the outlet liquid mole ratio is not less than the spec.
## Example

This code
~~~
import SepUnit

SepUnit.AbsorberDesign(0.0000067,0.0067,-0.67,1,-1.5,67)
~~~
will give this report
~~~
===== Calculation Report =====
Number of stages:                  3
Outlet liquid mole raio:           9.24793347065299e-05
Maximum feed mole ratio capable:   0.008445544156276913
Liquid to Feed ratio:              72.68119037769026
===== End of the report =====
~~~
and this plot