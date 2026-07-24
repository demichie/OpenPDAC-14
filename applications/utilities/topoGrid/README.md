# Mesh Deformation Procedure in **topoGridNew**

The **topoGridNew** utility deforms the computational mesh to conform to a
given topography. The process uses a combination of bilinear interpolation,
inverse distance weighting (IDW), and blending techniques to compute vertical
displacements and areas. The deformation accounts for contributions from the
base ($z = 0$) and the top boundary to compute the deformation at the internal
points of the mesh.

---

## 1. Bilinear Interpolation for Face Centers with $z = 0$

- **Objective**: Compute vertical displacements ($\\Delta z$) and areas for
  the centers of mesh faces located at $z = 0$, based on the input topography.

- **Procedure**:
  1. The $x, y$ coordinates of each face center are mapped to the
     corresponding cell in the input topography grid.

  1. **Bilinear interpolation** is performed using the four surrounding grid
     nodes ($v\_{00}, v\_{01}, v\_{10}, v\_{11}$) to estimate the interpolated
     elevation ($z\_{\\text{interpolated}}$):

  $$
  z\_{\\text{interpolated}} = v\_{00}(1 - x\_{\\text{lerp}})(1 - y\_{\\text{lerp}})
  \+ v\_{01}x\_{\\text{lerp}}(1 - y\_{\\text{lerp}})
  \+ v\_{10}(1 - x\_{\\text{lerp}})y\_{\\text{lerp}}
  \+ v\_{11}x\_{\\text{lerp}}y\_{\\text{lerp}}
  $$
  1. The vertical displacement ($\\Delta z$) is calculated as:

  $$
  \\Delta z\_{\\text{face}} = z\_{\\text{interpolated}} - z\_{\\text{original}}
  $$
  1. The area of each face ($\\text{Area}\_{\\text{face}}$) is computed based
     on its geometry and stored for later weighting.

- **Output**:
  - Each face center at $z = 0$ is associated with:
    - Vertical displacement ($\\Delta z\_{\\text{face}}$).
    - Area ($\\text{Area}\_{\\text{face}}$).

---

## 2. First Inverse Distance Weighting (IDW) for Mesh Points at $z = 0$

- **Objective**: Interpolate the vertical displacements
  ($\\Delta z\_{\\text{face}}$) and areas ($\\text{Area}\_{\\text{face}}$)
  from the face centers to the mesh points at $z = 0$.

- **Procedure**:
  1. For each mesh point at $z = 0$, calculate its distances ($d_i$) to all
     $z = 0$ face centers.

  1. Compute the **minimum distance** ($d\_{\\text{min}}$) to the face
     centers and define a threshold distance:

  $$
  d\_{\\text{threshold}} = \\text{interpRelRadius} \\cdot d\_{\\text{min}}
  $$

  $\\text{interpRelRadius}$: A user-defined multiplier (default = 4).
  1. Include only face centers within $d\_{\\text{threshold}}$ in the summation.

  1. Compute weights ($w_i$) for the included face centers:

  $$
  w_i = \\frac{1}{d_i}, \\quad d_i \\leq d\_{\\text{threshold}}
  $$
  1. Interpolate the vertical displacement ($\\Delta z$) as:

  $$
  \\Delta z\_{\\text{p}} = \\frac{\\sum_i w_i \\Delta z\_{\\text{face},i}}{\\sum_i w_i}
  $$
  1. Similarly, interpolate the area ($\\text{Area}\_{\\text{p}}$):

  $$Area_p = \\frac{\\sum_i w_i \\text{Area}\_{\\text{face},i}}{\\sum_i w_i}$$

- **Output**:
  - Each mesh point at $z = 0$ is associated with:
    - Vertical displacement ($\\Delta z\_{p}$).
    - Interpolated area ($\\text{Area}\_{p}$).

---

## 3. Merging Displacements and Areas

- **Objective**: Combine the interpolated values from mesh points at $z = 0$
  with the fixed displacements and areas associated with the centers of the
  top-face boundary.

- **Procedure**:
  1. Prescribed vertical displacements ($\\Delta z\_{\\text{top}}$) and
     areas are assigned to the top-face centers, typically based on the
     maximum topography height.
  1. Combine these top-face values with the $z = 0$ mesh-point values to
     form a single global list:
       - $\\Delta z\_{\\text{global}}$: Merged list of vertical
         displacements.
       - $\\text{Area}\_{\\text{global}}$: Merged list of areas.

- **Output**:
  - A single global list of:
    - Vertical displacements ($\\Delta z\_{\\text{global}}$).
    - Areas ($\\text{Area}\_{\\text{global}}$).

---

## 4. Second Inverse Distance Weighting (IDW) for Internal Mesh Points

- **Objective**: Compute the vertical deformation for all internal mesh points
  using the global list of displacements and areas.

### **3D IDW Interpolation**:

1. For each internal mesh point, compute the weights ($w_i$) for all global
   points using:

   $$
   w_i = \\text{Area}\_{i} \\cdot \\left[\\left(\\frac{L}{d_i}\\right)^3 + \\left(\\alpha \\cdot \\frac{L}{d_i}\\right)^5\\right]
   $$
   - $L$: Estimated length of the deformation region.
   - $\\alpha$: Fraction of $L$, representing the near-body influence region.
   - $d_i$: Euclidean distance to the global point.

1. Interpolate the vertical deformation ($\\Delta z\_{\\text{3D}}$):

   $$
   \\Delta z\_{\\text{3D}} = \\frac{\\sum_i w_i \\Delta z_i}{\\sum_i w_i}
   $$

### **Bottom-Up Interpolation**:

- This follows the same procedure as **Step 2**, using only $z = 0$ face
  centers.

### **Blending**:

1. Compute the relative height of the mesh point:

   $$
   z\_{\\text{rel}} = \\frac{z\_{\\text{mesh}} - z\_{\\text{min}}}{z\_{\\text{max}} - z\_{\\text{min}}}
   $$

1. Blend the deformations:

   $$
   \\Delta z\_{\\text{mesh}} = z\_{\\text{rel}} \\cdot \\Delta z\_{\\text{3D}} + (1 - z\_{\\text{rel}}) \\cdot \\Delta z\_{\\text{2D}}
   $$

- **Output**:
  - Vertical deformation ($\\Delta z\_{\\text{mesh}}$) for all internal mesh points.

---

This procedure ensures smooth deformation, accurately capturing local
influences at $z = 0$ and global effects across the mesh.
