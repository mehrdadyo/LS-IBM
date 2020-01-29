# LS-IBM
## This code is under developement, i.e. not final version.
A level set immersed boundary method for reactive transport. The level set method is used to find the exact location of the evolving solid-fluid interface. Then with the knowledge of the interface we can get the normals and use the immersed boundary method (IBM) to enforce the no-slip BC for velocity as well as Robin or Neumann BC for reactive transport.


The IBM is explanied in the following paper:
Mehrdad Yousefzadeh, Ilenia Battiato, High order ghost-cell immersed boundary method for generalized boundary conditions,
International Journal of Heat and Mass Transfer, Volume 137, 2019, Pages 585-598, ISSN 0017-9310

https://doi.org/10.1016/j.ijheatmasstransfer.2019.03.061.

(http://www.sciencedirect.com/science/article/pii/S0017931018357636)

# Abstract: 
Flow and reactive transport problems in engineering, medical and environmental applications often involve complex geometries. Grid based methods (e.g. finite volume, finite element, etc.) are a vital tool for studying such problems. Cartesian grids are one of the most attractive options as they possess simple discretization stencils and are usually straightforward to generate at roughly no computational cost. The Immersed Boundary Method, a Cartesian based methodology, maintains most of the useful features of structured grids, while it exhibits a great resilience in dealing with complex geometries. These features make it increasingly more attractive to model transport in evolving porous media as the cost of grid generation reduces greatly. Yet, stability issues due to the geometry of the interpolation stencil combined with limited studies on the implementation of Neumann (constant flux) and linear Robin (e.g. reaction) boundary conditions have significantly limited its applicability to transport in complex topologies. We develop a high-order compact Cartesian model based on ghost cell immersed boundary method for incompressible flow and scalar transport subject to different boundary conditions. The accuracy test shows at least second order of accuracy in L1,L2 and L∞ norms of error. The proposed method is capable of accurately capturing the transport physics near the boundaries for Dirichlet, Neumann and Robin boundary conditions. We tested the method for several transport and flow scenarios, including heat transfer close to an immersed object and mass transport over reactive surfaces.

# Keywords: 
Immersed boundary method; Fluid-solid interaction; High order discretization; Robin boundary condition
