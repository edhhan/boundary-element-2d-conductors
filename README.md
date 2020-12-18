# Boundary-element-2d-conductors

In this project, we wanted to familiarize ourselves with the boundary element method (BEM) for 2d coaxial conductors (the depth is supposed to be infinite) in a context of electromagnetism. To do so, we have solved a well-known problem in physics, the 2D Laplace potential with various geometries. 

We used the simplest form of BEM : each element has a constant value over his domain and the element have a linear geometry. 

The project was done in a context of a summer internship at Polytechnique Montréal University, with Pr. Frédiric Sirois as the supervisor. 

## Difference with classical finite element method
The BEM only requires a discretization of the boundary. Also, only one type of condition is imposed on the boundary : either (exclusive) Dirichlet or Neumann. 

### Advantages
1. The interior or outer domain is not discretized with elements. The solution can be easily reconstructed in all space when the BEM has solved entirely the boundary.
2. Less elements are required.

### Disadvantages
1. Each boundary element are treated as a physical source. The interaction between the elements is treated with a Green-integral-source formalism. Thus, the problem becomes more complex with the physical potential. Plus, a considerable amount of mathematical work is required before implementing the BEM. 

## Exemples and results

- Coaxial conductor 
2D Laplace with Dirichlet condition of 0 for the outer circle and 1 for the inner circle.
![Elements](https://github.com/edhhan/Boundary-element-2d-conductors/main/results/coaxial_elements.png?raw=true)
![Domain](https://github.com/edhhan/Boundary-element-2d-conductors/main/results/coaxial_domain.png?raw=true)

- Hole conductor (outer square with inner circle)
2D Laplace with Dirichlet condition of 0 for the bottom and 1 for the top (opposite sides). Other sides are Neumann with 0.


- Slope 
2D Laplace with Dirichlet condition of


## Author
Edward H-Hannan

## Acknowledgments
Pr. Frédéric Sirois 
