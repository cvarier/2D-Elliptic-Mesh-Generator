# 2D Elliptic Mesh Generator
This is a powerful <b>2D orthogonal elliptic mesh (grid) generator</b> which works by solving the <b>Winslow partial differential equations</b>.
It is capable of modifying the meshes with <b>stretching functions</b> and an <b>orthogonality adjustment</b> algorithm. This algorithm works by calculating curve slopes using a tilted parabola tangent line fitter (original discovery). The mesh generator is packaged as a Java program which can be compiled and executed via the command line. 

The program allows one to choose
from six different boundary types: rectangular, Gaussian, absolute value, greatest-integer, forwards step and semi-ellipse. One can then
specify the coordinates of the mesh domain (<b><i>warning: the domain must be perfectly square</i></b>). Finally, one can choose to add refinements
to the mesh, such as <b>orthogonality</b> adjustment and <b>stretching functions</b>. The program will then generate an initial coarse mesh and iteratively refine it to produce a <b>smooth mesh</b> with the given parameters and refinement options. 

A distinct feature of the elliptic mesh solver is that it <b>corrects overlapping and misplaced grid lines</b> very well. A detailed analysis of the quality of the resulting mesh will also be provided. 

## Screenshots
Here are some examples of meshes generated with the program (<b><i>initial</i></b> in <b>blue</b> and <b><i>final</i></b> in <b>green</b>):

<p align="center"><img src ="https://cloud.githubusercontent.com/assets/16710726/22316883/234eed46-e33e-11e6-8913-56a7a99bd28b.png" width = "400" height = "400"/><img src ="https://cloud.githubusercontent.com/assets/16710726/22316880/234d9f72-e33e-11e6-8e5b-9222c5ea2d6d.png" width = "400" height = "400"/></p>

<p align="center"><img src ="https://cloud.githubusercontent.com/assets/16710726/22316879/234c80c4-e33e-11e6-9d17-83fc904a4cd4.png" width= "410" height = "410"/><img src ="https://cloud.githubusercontent.com/assets/16710726/22316882/234ec23a-e33e-11e6-84b4-b82c847caa85.png" width= "400" height = "400"/></p>

<p align="center">&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp<img src ="https://cloud.githubusercontent.com/assets/16710726/22316974/d329dbf4-e33e-11e6-88b3-5b16058221a7.png" width= "320" height = "320"/>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp<img src ="https://cloud.githubusercontent.com/assets/16710726/22316973/d329c68c-e33e-11e6-8de2-8b475ad00936.png" width = "329" height = "320"/><br><br></p>


A more complete collection can be found within the `Screenshots` folder.

## Mathematical Framework
The algorithms implemented in this project are mainly founded upon the principles of differential geometry and tensor calculus. The most fundamental concept behind the mathematics governing this project is the transition between coordinate systems in order to obtain the solution to partial differential equations in the most efficient manner possible. More specifically, in the context of this project, this implies transforming a set of equations written in Cartesian coordinates to curvilinear coordinates. This concept can be extended to any number of spatial dimensions, which will later be shown. However, a two-dimensional solution was developed in order to demonstrate the feasibility of generating smooth elliptic grids with the prescribed specifications.

Consider a system of <i>n</i> dimensions which can be represented by the set of Cartesian coordinates <img src ="https://user-images.githubusercontent.com/16710726/31331568-d764704e-acb0-11e7-9bfe-f731ab10bc6a.gif" /> and the set of curvilinear coordinates <img src = "https://user-images.githubusercontent.com/16710726/31331564-d3a35088-acb0-11e7-9b25-68c3bcfdec27.gif"/>. We define the covariant metric tensor <i>g<sub>ij</sub></i> to be:

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31331597-00423ea6-acb1-11e7-896d-48ba46dd996d.gif"/></p>

where each of the partial derivatives comes from the definition of the covariant base vectors. Each of these base vectors describe how one coordinate system changes with respect to another, when any particular coordinate is held fixed. For our two-dimensional problem, we can expand these sums to yield the following expressions for the metric tensors:

<p align="center">&nbsp&nbsp&nbsp<img src="https://user-images.githubusercontent.com/16710726/31332560-d5ca47aa-acb4-11e7-9524-27911f357c6b.gif" /></p>
<p align="center">&nbsp&nbsp&nbsp<img src="https://user-images.githubusercontent.com/16710726/31332553-d0a63a04-acb4-11e7-80c7-b2340ca48d0f.gif" /></p>
<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31332548-cc300bb2-acb4-11e7-8c28-8129ebcd557d.gif" /></p>

where (*x*, *y*) represents our Cartesian coordinate system and <img src="https://user-images.githubusercontent.com/16710726/31160211-cb4fb862-a89c-11e7-8550-4873dc518b45.gif" /> represents our curvilinear coordinate system. We could also define the contravariant base vectors and metric tensors if we wished to solve the inverse problem, which would be transitioning from curvilinear to Cartesian coordinates.

## Elliptic Mesh Generation Algorithm
Firstly to construct an initial mesh, the <b>Transfinite Interpolation algorithm</b> is applied to the given domain constrained by the
specified boundary conditions. This algorithm is implemented by mapping each point within the domain (regardless of the boundaries) to a new domain existing within the boundaries. This algorithm works by iteratively solving the parametric vector equation

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31158969-99999b48-a893-11e7-8c8f-f625ff310829.gif" /></p>

where <img src ="https://user-images.githubusercontent.com/16710726/31159037-0afef594-a894-11e7-9bda-406151b5590b.gif" /> and <img src ="https://user-images.githubusercontent.com/16710726/31159055-3326d938-a894-11e7-95c2-97329c6cd80f.gif" /> represent parameters in the original domain and <img src="https://user-images.githubusercontent.com/16710726/31159084-6cb8a92e-a894-11e7-8ec6-b92841c48b54.gif" />, <img src="https://user-images.githubusercontent.com/16710726/31159087-810fd820-a894-11e7-8db0-22c5a7300895.gif" />, <img src="https://user-images.githubusercontent.com/16710726/31159103-9722f26e-a894-11e7-83de-ff8bb44f5a03.gif" /> and <img src="https://user-images.githubusercontent.com/16710726/31159111-a60e2816-a894-11e7-803b-ede51eeac4b0.gif" /> represent the curves defining the left, top, bottom and right boundaries. *P<sub>ij</sub>* represents the point of intersection between curve <img src="https://user-images.githubusercontent.com/16710726/31159401-b654631e-a896-11e7-811f-64737a10f14c.gif" /> and <img src="https://user-images.githubusercontent.com/16710726/31159403-bbffb994-a896-11e7-9832-66573661a963.gif" />.

At the heart of the solver is the mesh smoothing algorithm, which at a high level, works by solving the pair of Laplace equations

<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31159558-f3981bd4-a897-11e7-9d1c-40bc4e531f6b.gif" />    and    <img src="https://user-images.githubusercontent.com/16710726/31159563-fbd4eaca-a897-11e7-9c79-5c77eb501134.gif" /></p>

where <img src="https://user-images.githubusercontent.com/16710726/31160309-7a2a7b2e-a89d-11e7-8b7a-f7fd86db0e0d.gif" /> and <img src="https://user-images.githubusercontent.com/16710726/31159710-f8d219a0-a898-11e7-9195-7c403297e18f.gif" /> represent the *x* and *y* coordinates of every point in the target domain, mapped to a transformed, computational space using the change of variables method. This renders the calculations simpler and faster to compute. However, we wish to solve the inverse problem, where we transition from the computational space to the curvilinear solution space. Using tensor mathematics, it can be shown that this problem entails solving the equations

<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31159981-f7c7ad2a-a89a-11e7-9981-37efc81753c5.gif" /></p>
<p align="center">and</p>
<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31160011-2c9d829a-a89b-11e7-8efb-b28433e13e8f.gif" /></p>

where *g<sub>ij</sub>* is the covariant metric tensor at entry (*i*,*j*) within the matrix of covariant tensor components defining the mapping of the computational space coordinates <img src="https://user-images.githubusercontent.com/16710726/31160211-cb4fb862-a89c-11e7-8550-4873dc518b45.gif" /> onto the physical solution space coordinates (*x*,*y*). In this model, *x* and *y* are computed as functions of <img src="https://user-images.githubusercontent.com/16710726/31160309-7a2a7b2e-a89d-11e7-8b7a-f7fd86db0e0d.gif" /> and <img src="https://user-images.githubusercontent.com/16710726/31159710-f8d219a0-a898-11e7-9195-7c403297e18f.gif" />.

This set of equations are the elliptic PDEs known as the Winslow equations. These are applied to the mesh using the method of <b>mixed-order finite differences</b> on the partial derivatives (and tensor coefficients, as they are a function of these derivatives), thereby resulting in the equations (for a single node):

<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31333974-28d59be8-acba-11e7-8038-faa0ee406884.gif" /></p>
<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31333975-28d72648-acba-11e7-96b1-fb4ad8830a02.gif" /></p>
<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31333976-28d6ab5a-acba-11e7-8702-d83acfa84786.gif" /></p>
<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31333546-8439ea40-acb8-11e7-8a73-d6475b8cadbb.gif" /></p>
<p align="center">and</p>
<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31333563-8dee7fb0-acb8-11e7-8dd2-680b72ad5227.gif" />,</p>

where *i* and *j* are the coordinates of a node in the mesh in computational space. Here <img src="https://user-images.githubusercontent.com/16710726/31160826-788660d6-a8a1-11e7-9088-d36b7147147f.gif" /> and <img src="https://user-images.githubusercontent.com/16710726/31160836-8b33a8ba-a8a1-11e7-91ac-7c4a3c63fbf6.gif" /> are equal increments in <img src="https://user-images.githubusercontent.com/16710726/31160309-7a2a7b2e-a89d-11e7-8b7a-f7fd86db0e0d.gif" /> and <img src="https://user-images.githubusercontent.com/16710726/31159710-f8d219a0-a898-11e7-9195-7c403297e18f.gif" /> respectively.

The coefficients for these equations can be generated for each point to form a system of linear equations, which is then modeled in matrix representation, resulting in a tri-diagonal matrix. This matrix is then solved iteratively using the <b>Thomas Tri-Diagonal Matrix
Algorithm</b> line-by-line by traversing from the bottom up.

The solution to the matrix generated from a single iteration can then be further processed by the orthogonality adjustment algorithm and stretching
function methods as necessary. The solver then calculates the solution for all other node lines and repeats the process until the difference between adjacent nodes meets a threshold convergence criteria.

## Orthogonality Adjustment Algorithm
In several <b>computational fluid dynamics</b> applications, an orthogonal mesh is necessary in certain regions to ensure a high enough accuracy when performing calculations. However, it is not always possible to achieve a fully orthogonal solution, and thus the problem becomes finding a nearly-orthogonal solution to an arbitrarily defined domain. 

The implemented solution uses an iterative approach to find the angles of intersection and adjust the position of the nodes until their respective angles of intersection converge to a reasonable threshold value from 90 degrees. The exact method makes use of the <b>linear approximation</b> of the grid lines intersecting at each node within the mesh. 

A remarkable result from the research was the development of an accurate method for obtaining these linear approximations. This method consists of fitting a tilted parabola to three adjacent nodes, which are defined as (*x<sub>1</sub>*, *y<sub>1</sub>*), (*x<sub>2</sub>*, *y<sub>2</sub>*) and the node in between these two. By applying <b>coordinate transformations</b>, we can obtain the trigonometric function

<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31161603-28f37b1c-a8a6-11e7-9b0c-08bddb27ed6a.gif" /></p>

whose roots can be solved for using the bisection method, which represent the angular position of the parabola (can be improved with Newton's method). 

The same process is applied to the three oppositely adjacent nodes. From this, a suitable linear approximation is obtained, and the adjustment is determined by plugging the slopes of the two linear functions into the linear equation relating the two derived using basic analytical geometry. This describes the system of equations

<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31161932-6a369012-a8a8-11e7-994c-ba75d709aedd.gif" /></p>
<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31162043-214b9144-a8a9-11e7-9905-252cfce64f09.gif" /></p>

for the vertical grid line *V* and the horizontal grid line *H* at a given iteration *k*. Thus, this system is solved iteratively for each mesh node in the same line-by-line fashion as the Winslow system solver.

## Stretching Functions
In order to further improve the quality of the mesh, one can introduce <b>univariate stretching functions</b> to either compress or expand grid lines in order to correct non-uniformity where grid lines are more or less dense. These functions are arbitrarily chosen and only reflect the distribution of grid lines. 

We can derive a new set of equations by combining our previously established differential model for grid generation and a set of univariate stretching functions of our choice. In order to do so in a straightforward manner, we can transform our Cartesian coordinates (*x*, *y*) to a new set of coordinates (*χ*, *σ*) which exist in a different space, *χσ*, called the parameter space. We define our stretching functions, *f<sub>1</sub>* and *f<sub>2</sub>* as onto and one-to-one univariate functions of *ξ* and *η* respectively. They are described as follows: 

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31361649-47fd6178-ad22-11e7-8721-ff9e2a12096b.gif" /></p>
<p align="center">and</p>
<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31361650-47fe6f82-ad22-11e7-8287-0efd31004c6a.gif" /></p>

If we define the mapping from physical space to parameter space as the same elliptic system as before, we get the equations

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31362584-45ee28b8-ad27-11e7-935d-4384b68feb87.gif" /></p>
<p align="center">and</p>
<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31362576-3bf1b852-ad27-11e7-9ad8-4c4f02d03026.gif" /></p>

To obtain an expression for this system, we can get the first order partial derivative with respect to *x* as so:

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31362595-519fe3ea-ad27-11e7-8ae7-43dc57faf68e.gif" /></p>

Using the chain rule, we can get the second order partial derivative:

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31362604-5b86f588-ad27-11e7-9119-93f83a4c48f5.gif" /></p>

Similarly for *y*, we get:

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31362625-6eed8c68-ad27-11e7-83d2-81542f8a89c0.gif" /></p>

Adding the previous two equations, and rearranging, we get:

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31362634-8456d7f8-ad27-11e7-9a99-9a097dcbc75c.gif" /></p>

Similarly for *η*, we get:

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31362664-b087ab86-ad27-11e7-9b18-1d69d0bf99f8.gif" /></p>

However, we wish to solve the inverse problem as before, and thus in finding the inverse of these equations, we obtain the Poisson system 

<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31162397-29a18d7e-a8ab-11e7-910c-34050762fdf6.gif" /></p>
<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31162408-40ed191c-a8ab-11e7-86bb-f82e3f92fe87.gif" /></p>

However, if we wish to maintain orthogonality everywhere in the grid, then we must eliminate the *g<sub>12</sub>* terms to get the equations

<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31363152-768f1600-ad2a-11e7-9353-deadd6ddb6db.gif" /></p>
<p align="center"><img src="https://user-images.githubusercontent.com/16710726/31363168-91f09108-ad2a-11e7-86aa-1407ba271ec5.gif" /></p>

This slightly modifies the solution process by changing the values of the matrix coefficients in the TDMA setup.

The discretized version of the new system is

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31335263-c913f136-acbf-11e7-97a6-53dd0fa2a5d4.gif" /></p>
<p align="center">and</p>
<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31335260-c6cf7332-acbf-11e7-9243-5a16ac75e1e7.gif" /></p>

where the metric tensor coefficients are the same as before. 

The default stretching functions used in the program are

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31335710-bfe3db2e-acc1-11e7-9e86-f5010fcbf04f.gif" /></p>
<p align="center">and</p>
<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31335711-bfe475f2-acc1-11e7-8573-2032fa4f471c.gif" /></p>

where |*∝*| << 1 and |*β*| << 1. The parameters *∝* and *β*, when positive, indicate how much stretching occurs in the *x* and *y* directions respectively. When these values are negative, compression of grid lines will occur instead.

## Extension to Higher Dimensions
If we wished to extend the elliptic solver to 3D, we would need to develop equations for transitioning from a curvilinear coordinate system (*ξ*, *η*, *ς*) to a 3D Cartesian coordinate system (*x*, *y*, *z*). Using the same elliptic model as before, we get the 3D version of the Winslow equations

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31380152-7e1f4d78-ad7e-11e7-8302-54131024b983.gif" /></p>
<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31380154-7e26f85c-ad7e-11e7-973a-a818ada37008.gif" /></p>
<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/31380153-7e1fcb7c-ad7e-11e7-9052-8287e4de8d3b.gif" /></p>

where each of the metric tensor coefficients is determined by taking the cofactors of the contravariant tensor matrix. The contravariant tensor matrix is used to obtain the coefficients for the Winslow equations, which are the inverse of the Laplace equations as stated before.

In general, if we wish to extend our elliptic mesh solver to *n* dimensions, then we will have *n* sets of equations each with *n*!/(2(*n*-2)!) + *n* terms. This renders the problem gradually more and more difficult to solve for higher dimensions with the existing elliptic scheme, implying that a different type of PDE might be needed in these cases.

Another complication that arises in higher dimensions is adjusting grid lines to enforce orthogonality. Using the aforementioned algorithm for adjusting grid lines to achieve either complete or partial orthogonality on the boundary, we would need to iteratively solve three sets of two linear equations for each node in the mesh, as well as solve three trigonometric equations per iteration to compute the tangents.

For two dimensions, if our grid domain was of size *n* x *n*, and we assumed an average of *k* iterations per mesh node, then our runtime complexity would be approximately <img src ="https://user-images.githubusercontent.com/16710726/31381680-3d5c9110-ad83-11e7-8bcd-f7d9f98d69a4.gif" />. However for three dimensions, our runtime complexity would be <img src ="https://user-images.githubusercontent.com/16710726/31381681-3d5ee9ba-ad83-11e7-95b4-fac0c46d5f9a.gif" />. If our domain had a side length of 100 and for each node, the solver took 10 iterations in the 2D version and 25 iterations in the 3D version, then the 3D version would take over 1000 times longer to return a solution. This means we would definitely need to seek a more efficient and powerful algorithm. However, another issue that may arise is that for three dimensions, orthogonal solutions to smooth meshes are very rare, and even getting a solution with the most efficient algorithm might be impossible for many cases.

## Mesh Quality Analysis Report
In order to determine the quality of the resulting mesh, it was necessary to construct an objective means of quality measurement. Therefore, several <b>statistical procedures</b> were implemented in the program to produce a <b>meaningful mesh quality analysis report</b>. The metrics which are presented are divided into the following categories:

* Orthogonality Metrics
  * Standard deviation of angles
  * Mean angle
  * Maximum deviation from 90 degrees
  * Percentage of angles within x degrees from 90 degrees (x can be set as a constant in the code)
  
* Cell Quality Metrics
  * Average aspect ratio of all cells
  * Standard deviation of all aspect ratios
  
In the above list, "angle" refers to the angle of intersection of grid lines at a node and "aspect ratio" refers to the skewness of a grid cell measured as a ratio of the cell's longest side to its shortest side.

## Libraries Used
* JMathPlot by Yann Richet

## License
© Chaitanya Varier 2017. This software is protected under the MIT License.

## Context and Acknowledgement
I wrote this software independently for Prof. Jonathon Sumner at Dawson College, who provided me with invaluable mentorship and advice throughout the course of my research.

## References
I gleaned all the necessary mathematical knowledge for completing this project from the following sources:

* Farrashkhalvat, M., and J. P. Miles. Basic structured grid generation with an introduction to unstructured grid generation. Oxford: Butterworth Heinemann, 2003. Print.

* Versteeg, H.K, and W. Malalasekera. An introduction to computational fluid dynamics: the finite volume method. Harlow: Pearson Education, 2011. Print.

* Knupp, Patrick M., and Stanly Steinberg. Fundamentals of grid generation. Boca Raton: CRC Press, 1994. Print.

