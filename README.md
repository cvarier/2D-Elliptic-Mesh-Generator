# 2D Elliptic Mesh Generator
This is a powerful <b>2D orthogonal elliptic mesh generator</b> which uses the <b>Winslow (modified Laplace) partial differential equations</b>.
It also uses <b>univariate stretching functions</b> and a <b>tilted parabola tangent line fitter</b> (original discovery). The grid generator is packaged as a Java program which can be compiled and executed via the command line. 

The program allows one to choose
from six different boundary types: rectangular, Gaussian, absolute value, greatest-integer, forwards step and semi-ellipse. Then one must
specify the coordinates of the grid domain (<b><i>warning: the domain must be perfectly square</i></b>). Finally, one can choose to add refinements
to the grid, such as <b>orthogonality</b> adjustment and <b>stretching functions</b>. The program will then generate an initial course grid and iteratively refine it to produce a <b>smooth grid</b> with the given parameters and refinement options. 

A distinct feature of the elliptic grid solver is that it <b>corrects overlapping and misplaced gridlines</b> very well. A detailed analysis of the quality of the resulting grid will also be provided. 

## Screenshots
Here are some examples of grids generated with the program (<b><i>initial</i></b> in <b>blue</b> and <b><i>final</i></b> in <b>green</b>):

<p align="center"><img src ="https://cloud.githubusercontent.com/assets/16710726/22316883/234eed46-e33e-11e6-8913-56a7a99bd28b.png" /></p>
<p align="center"><img src ="https://cloud.githubusercontent.com/assets/16710726/22316880/234d9f72-e33e-11e6-8e5b-9222c5ea2d6d.png" /></p>
<p align="center"><img src ="https://cloud.githubusercontent.com/assets/16710726/22316879/234c80c4-e33e-11e6-9d17-83fc904a4cd4.png" /></p>
<p align="center"><img src ="https://cloud.githubusercontent.com/assets/16710726/22316882/234ec23a-e33e-11e6-84b4-b82c847caa85.png" /></p>
<p align="center"><img src ="https://cloud.githubusercontent.com/assets/16710726/22316974/d329dbf4-e33e-11e6-88b3-5b16058221a7.png" /></p>
<p align="center"><img src ="https://cloud.githubusercontent.com/assets/16710726/22316973/d329c68c-e33e-11e6-8de2-8b475ad00936.png" /></p>

A more complete collection can be found within the `Screenshots` folder.

## Elliptic Grid Generation Algorithm
Firstly to construct an initial grid, the <b>Transfinite Interpolation algorithm</b> is applied to the given domain constrained by the
specified boundary conditions. This algorithm is implemented by mapping each point within the domain (regardless of the boundaries) to a new domain existing within the boundaries. This algorithm works by iteratively solving the parametric vector equation

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31158969-99999b48-a893-11e7-8c8f-f625ff310829.gif" /></p>

where <img src ="https://user-images.githubusercontent.com/16710726/31159037-0afef594-a894-11e7-9bda-406151b5590b.gif" /> and <img src ="https://user-images.githubusercontent.com/16710726/31159055-3326d938-a894-11e7-95c2-97329c6cd80f.gif" /> represent parameters in the original domain and <img src="https://user-images.githubusercontent.com/16710726/31159084-6cb8a92e-a894-11e7-8ec6-b92841c48b54.gif" />, <img src="https://user-images.githubusercontent.com/16710726/31159087-810fd820-a894-11e7-8db0-22c5a7300895.gif" />, <img src="https://user-images.githubusercontent.com/16710726/31159103-9722f26e-a894-11e7-83de-ff8bb44f5a03.gif" /> and <img src="https://user-images.githubusercontent.com/16710726/31159111-a60e2816-a894-11e7-803b-ede51eeac4b0.gif" /> represent the curves defining the left, top, bottom and right boundaries. *P<sub>ij</sub>* represents the point of intersection between curve <img src="https://user-images.githubusercontent.com/16710726/31159401-b654631e-a896-11e7-811f-64737a10f14c.gif" /> and <img src="https://user-images.githubusercontent.com/16710726/31159403-bbffb994-a896-11e7-9832-66573661a963.gif" />.

At the heart of the solver is the grid smoothing algorithm, which at a high level, works by solving the pair of Laplace equations

<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31159558-f3981bd4-a897-11e7-9d1c-40bc4e531f6b.gif" />    and    <img src="https://user-images.githubusercontent.com/16710726/31159563-fbd4eaca-a897-11e7-9c79-5c77eb501134.gif" />.</p>

Next, the Winslow equations are applied to the grid using the method of <b>mixed-order finite differences</b>, thereby generating a system of equations for each one-dimensional line of nodes in the grid. 

This system of equations is then modeled in matrix representation, resulting in a tri-diagonal matrix. This matrix is then solved using the <b>Thomas Tri-Diagonal Matrix
Algorithm</b>. 

The solution to the current iteration is then further processed by the orthogonality adjustment algorithm and stretching
function methods as necessary. The solver then calculates the solution for all other node lines and repeats the process until the difference between adjacent nodes meets a threshold convergence criteria.

## Orthogonality Adjustment Algorithm
In several <b>computational fluid dynamics</b> applications, an orthogonal mesh is necessary in certain regions to ensure a high enough accuracy when performing calculations. However, it is not always possible to achieve a fully orthogonal solution, and thus the problem becomes finding a nearly-orthogonal solution to an arbitrarily defined domain. 

The implemented solution uses an iterative approach to find the angles of intersection and adjust the position of the nodes until their respective angles of intersection converge to a reasonable threshold value from 90 degrees. The exact method makes use of the <b>linear approximation</b> of the grid lines intersecting at each node within the grid. 

A remarkable result from the research was the development of an accurate method for obtaining these linear approximations. This method consists of fitting a tilted parabola to three adjacent nodes using <b>coordinate transformations</b>. The resulting trigonometrically transformed equation of the parabola is solved using the bisection method for the angular position of the parabola (can be improved with Newton's method). 

The same process is applied to the three oppositely adjacent nodes. From this, a suitable linear approximation is obtained, and the adjustment is obtained by plugging the slopes of the two linear functions into the linear equation relating the two obtained by basic analytical geometry.

## Stretching Functions
In order to further improve the quality of the grid, one can introduce <b>univariate stretching functions</b> to either compress or expand grid lines in order to correct non-uniformity where grid lines are more or less dense. These functions are arbitrarily chosen and only reflect the distribution of grid lines. Upon implementation, the Winslow system becomes a Poisson system, thereby slightly modifying the solution process by changing the values of the matrix coefficients.

## Grid Quality Analysis Report
In order to determine the quality of the resulting mesh, it was necessary to construct an objective means of quality measurement. Therefore, several <b>statistical procedures</b> were implemented in the program to produce a <b>meaningful grid quality analysis report</b>. The metrics which are presented are divided into the following categories:

* Orthogonality Metrics
  * Standard deviation of angles
  * Mean angle
  * Maximum deviation from 90 degrees
  * Percentage of angles within x degrees from 90 degrees (x can be set as a constant in the code)
  
* Cell Quality Metrics
  * Average aspect ratio of all cells
  * Standard deviation of all aspect ratios
  
In the above list, "angle" refers to the angle of intersection of grid lines at a node and "aspect ratio" refers to the skewness of a grid cell measured as a ratio of the height to the length of the cell.

## Libraries Used
* JMathPlot by Yann Richet

## License
© Chaitanya Varier 2017. This software is protected under the MIT License.

## Context/Acknowledgement
This software was independently written for Prof. Jonathon Sumner at Dawson College as part of a research internship. A research paper detailing the mathematical results of the project is currently in the works.

## References
I gleaned all the necessary mathematical knowledge for completing this project from the following sources:

* Farrashkhalvat, M., and J. P. Miles. Basic structured grid generation with an introduction to unstructured grid generation. Oxford: Butterworth Heinemann, 2003. Print.

* Versteeg, H.K, and W. Malalasekera. An introduction to computational fluid dynamics: the finite volume method. Harlow: Pearson Education, 2011. Print.

* Knupp, Patrick M., and Stanly Steinberg. Fundamentals of grid generation. Boca Raton: CRC Press, 1994. Print.

