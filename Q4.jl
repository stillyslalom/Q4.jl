# Calculate basis functions for arbitrary hx, hy ===============================
using SymPy
function ϕQ4()
  @syms x y hx hy   # initialize symbolic variables
  pts = [(0, 0), (hx, 0), (hx, hy), (0, hy)]  # list elemental nodes
  # evaluate shape functions at each elemental node
  dxy = Sym.([f(p) for f in [p->1, p->p[1], p->p[2], p->p[1]*p[2]], p in pts])

  λ = (dxy\I) * [1, x, y, x*y]  # compute shape functions
  ∇λ = diff.(λ, [x y])          # compute gradient of shape functions
  ∇²λ = ∇λ*∇λ.'                 # outer product of gradient matrices
  ∫λ   = [lambdify(integrate(λ[i]*λ[j], (x,0,hx), (y,0,hy))) for i=1:4, j=1:4]
  ∫∇²λ = [lambdify(integrate(∇²λ[i,j],  (x,0,hx), (y,0,hy))) for i=1:4, j=1:4]
  return ∫λ, ∫∇²λ   # return two matrices of basis integrals w/ arguments hx, hy
end
const ∫λ, ∫∇²λ = ϕQ4()

function Q4(f, xx::StepRangeLen, yy = xx)
  nx, ny = length(xx), length(yy)
  nn = nx*ny              # Number of nodes
  ne = (nx - 1)*(ny - 1)  # Number of elements
  Zx, Zy = repmat(xx, ny), vec(repmat(yy', nx)) # Nodal coordinates
  T = zeros(Int, 4, ne)   # Element member nodes
  k = 1
  for j = 1:ny - 1, i = 1:nx - 1  # for each element, find corresponding nodes
    T[1, k] = k + (j-1)             # bottom left
    T[2, k] = k + (j-1) + 1         # bottom right
    T[3, k] = k + (j-1) + nx + 1    # top right
    T[4, k] = k + (j-1) + nx        # top left
    k += 1
  end

  # Basis function matrix instantiation
  ϕϕ   =  Float64[∫λ[i,j](step(xx), step(yy))   for i = 1:4, j = 1:4]
  ϕ′ϕ′ =  Float64[∫∇²λ[i,j](step(xx), step(yy)) for i = 1:4, j = 1:4]
  a = (ϕϕ + ϕ′ϕ′)

  # Stiffness matrix assembly
  A = spzeros(nn, nn) # Global stiffness matrix
  b = zeros(nn)       # Global load vector
  w = zeros(4)        # Element load vector
  for k = 1:ne, ic = 1:4
    i = T[ic, k]
    for jc = 1:4
      j = T[jc, k]
      A[i,j] += a[ic,jc]
      w[jc] = f(Zx[j], Zy[j])
    end
    b[i] += ϕϕ[ic, :]⋅w
  end
  return A, b
end
