# Bone Volume Solver: Linear Volume Distribution for Skeletal Meshes

A solver for computing the per-bone volume of a non-manifold mesh for procedural animation use cases. 

 * [Terms](#terms)
 * [Algorithm Overview](#algorithm-overview)
 * [Volume Function Computation](#volume-function-computation)
 * [Hierarchical Error Correction](#hierarchical-error-correction)
 * [Final Volume Extraction](#final-volume-extraction)
 * [Extensions](#extensions)


## Terms

- **Volume Function**: A linear function `ax + by + cz + d` that gives the volume when evaluated at any 3D point, for a manifold mesh a,b, and c will all be 0. 
- **Manifold Mesh**: Best understood through non-manifold topology: such as the eye holes in a character model, or the way that the teeth are missing faces inside gums, or that sleeves are just double sided faces, all of these are NOT manifold. Manifold means it has no holes, all the surfaces are closed, there's no one sided surfaces, and a bunch of other constraints. In practice: almost no meshes are manifold.
- **Subtree Volume**: The total volume of a bone plus all its children
- **Tetrahedron Volume**: Volume between a triangle and an origin point, inherently linear in the origin position
- **Tetrahedra Summation**: We create a tetrahedra for each tri on the mesh and some origin point, for tris facing away from the origin the contribution is positive, and for tris facing towards it the contribution is negative. So if you imagine 2 tris, one facing up, and below that one facing down, the first one's tetrahedra is all the space between it and you, and the second is negative and all the space between it and you, so the second *removes* the part that is outside of the area bound by the 2 tris. So we just get the area between the two.
For a watertight mesh with no holes, this means that the summation will be the same *regardless of what point you pick for the origin.*

---

## Algorithm Overview

This solver tackles a fundamental problem: how do you fairly distribute mesh volume across bones in a skeleton? Most approaches use voxels or geometric decomposition, but this algorithm exploits the mathematical fact that tetrahedron volumes are linear functions of position.

**The key insight**: Since tetrahedron volume = `(1/6) * triangle_area * dot(triangle_normal, origin_offset)`, and we're summing weighted tetrahedra, the result is always a linear function `ax + by + cz + d`. This means we can solve exact linear systems instead of doing expensive iterative approximations.

**What makes this good**:
- **Fast**: O(triangles + bones), no iterations needed
- **Exact**: R² = 1.0 (within floating point), not an approximation  
- **Reusable**: Same error correction works for volume, center of mass, moments of inertia
- **Robust**: Handles non-manifold meshes gracefully via error redistribution,

  For example: lets imagine a player character from the sega dreamcast, where they have a hand but no arm and a sleeve that's just a cylinder with no end cap. This is tragically non-manifold!

   So taking what I said before, if we compute the volume based on the origin, typically for a player character thats the ground between the feet, it's going to have this massive un-canceled contribution from between the interior of the sleeve to the ground! but if we instead measure it relative to the elbow, now it's a fairly sensible contribution that is visually about as much as the arm should have, problem solved!

**Limitations**:
- Requires skinned meshes with bone weights
- Assumes linear volume distribution (which is mathematically correct for tetrahedra)
- As written the algorithm assumes nodes are in strict DFS order.
- The per-vertex weights of each bone must sum to 1 or the sum of the volume functions is not equivalent to the volume function of the entire mesh. 
 
---

## Volume Function Computation

First, we compute a linear volume function for each bone. The math here is: tetrahedron volumes are inherently linear, so when you sum them with weights, you get another linear function.

```
GetVectorsPerBone(mesh) -> vec4[]
{
    // Get the 8 corners of the mesh bounding box
    corners[8] = {
        {min.x, min.y, min.z}, {min.x, min.y, max.z}, 
        {min.x, max.y, min.z}, {min.x, max.y, max.z},
        {max.x, min.y, min.z}, {max.x, min.y, max.z}, 
        {max.x, max.y, min.z}, {max.x, max.y, max.z}
    };
    
    // Storage for volume functions - one per bone
    volumeFunctions : vec4[numBones];
    volumeAccumulators : double[numBones][8];  // 8 values per bone, one per corner
    
    // Process each triangle in the mesh
    for(triangle in mesh.triangles)
    {
        // Calculate tetrahedron volume from triangle to origin point
        calculateVolume(origin) -> double
        {
            // Form vectors from origin to triangle vertices
            vectorAB = triangle.vertex[0] - origin;
            vectorAC = triangle.vertex[1] - origin;  
            vectorAD = triangle.vertex[2] - origin;
            
            // Signed volume = scalar triple product / 6
            return dot(vectorAB, cross(vectorAC, vectorAD)) / 6.0;
        }
        
        // Find which bones influence this triangle and by how much
        for(bone in triangle.influencingBones)
        {
            // Get skinning weights for this triangle's vertices on this bone
            weights = getSkinningWeights(triangle, bone);
            avgWeight = (weights.x + weights.y + weights.z) / 3.0;
            
            // Accumulate weighted volume contributions for each bounding box corner
            for(cornerIndex = 0; cornerIndex < 8; cornerIndex++)
            {
                volume = calculateVolume(corners[cornerIndex]);
                volumeAccumulators[bone][cornerIndex] += avgWeight * volume;
            }
        }
    }
    
    // For each bone, solve the linear system to get volume function coefficients
    for(bone = 0; bone < numBones; bone++)
    {
        // We now have 8 equations: ax + by + cz + d = volume for each corner
        // Solve for coefficients [a, b, c, d] using least squares
        // x = (A^T * A)^(-1) * A^T * b
        // I just used eigen for this
        volumeFunctions[bone] = solveLeastSquares(corners, volumeAccumulators[bone]);
    }
    
    return volumeFunctions;
}


```

The `GetVolumeVector` function solves an exact 8x4 linear system using QR decomposition. This isn't least-squares - it's an exact solution because the volume relationship is genuinely linear.

**Why 8 corners?** We need enough constraints to uniquely determine the 4 coefficients, and the bounding box corners give us a well-distributed set of points that spans the mesh space nicely.

---

## Hierarchical Error Correction

Here's where it gets clever. For a manifold mesh, the spatial coefficients `{a, b, c}` at the root should be `{0, 0, 0}` - meaning the total volume function has no spatial dependence. If they're not zero, we have an error that needs to be redistributed down the hierarchy.

```
// STEP 1: Accumulate volume functions up the tree to get subtree totals
subtreeFunctions = originalVolumeFunctions.copy();
for(boneIndex = numBones-1; boneIndex >= 0; boneIndex--)  // leaves to root
{
    parentIndex = parents[boneIndex];
    if(parentIndex >= 0)
        subtreeFunctions[parentIndex] += subtreeFunctions[boneIndex];
}

// STEP 2: Initialize error correction - root spatial coefficients should be {0,0,0}
errorCorrections : vec4[numBones] = {0};
rootError = subtreeFunctions[0].xyz;  // Extract just the spatial coefficients [a,b,c]
errorCorrections[0] = vec4(rootError, 0);  // Initialize root with spatial error, d=0

// STEP 3: Propagate error down the tree (root to leaves)
for(boneIndex = 0; boneIndex < numBones; boneIndex++)
{
    errorToHandle = errorCorrections[boneIndex];
    
    // Each bone tries to "absorb" as much error as it can handle
    for(axis = 0; axis < 3; axis++)  // x, y, z axes
    {
        boneCoeff = originalVolumeFunctions[boneIndex][axis];
        errorAmount = errorToHandle[axis];
        
        // Can only fix error if bone coefficient has same sign as error
        if(boneCoeff * errorAmount <= 0)
        {
            errorCorrections[boneIndex][axis] = 0;  // Can't help with this axis
            continue;
        }
        
        // Absorb up to the bone's own magnitude, no more
        maxCanAbsorb = min(abs(boneCoeff), abs(errorAmount));
        errorCorrections[boneIndex][axis] = maxCanAbsorb * sign(errorAmount);
        
        // Recalculate constant term: d = -dot(origin, [a,b,c])
        // The way to think about this, is imagine we extrude all the non-manifold parts of the mesh,
        // and join to a single vertex, that we place where the root of the bone affecting these vertices is. 
        // as in the prior example of the sleeve: the elbow becomes the root.
        // so the non-manifold issues have been contained in a reasonable way.
        errorCorrections[boneIndex].w = -dot(boneOrigins[boneIndex], errorCorrections[boneIndex].xyz);
    }
    
    // Calculate what error is left over after this bone absorbed what it could
    leftoverError = errorToHandle - errorCorrections[boneIndex];
    
    // Distribute leftover error to children based on their "capacity"
    totalChildCapacity : vec3 = {0, 0, 0};
    
    // First pass: calculate total capacity of all children for each axis
    for(child in children[boneIndex])
    {
        for(axis = 0; axis < 3; axis++)
        {
            childCoeff = subtreeFunctions[child][axis];
            errorAmount = leftoverError[axis];
            
            // Child can help if its coefficient has same sign as error
            if(childCoeff * errorAmount > 0)
                totalChildCapacity[axis] += abs(childCoeff);
        }
    }
    
    // Second pass: distribute error proportionally
    for(child in children[boneIndex])
    {
        for(axis = 0; axis < 3; axis++)
        {
            childCoeff = subtreeFunctions[child][axis];
            errorAmount = leftoverError[axis];
            
            // Skip if child can't help or no total capacity
            if(childCoeff * errorAmount <= 0 || totalChildCapacity[axis] == 0)
                continue;
                
            // Give child proportional share: (child_capacity / total_capacity) * error
            proportion = abs(childCoeff) / totalChildCapacity[axis];
            errorCorrections[child][axis] = errorAmount * proportion;
        }
    }
}
```

### What this does in plain English:

 - **Accumulate up**: Roll up all volume functions so each bone knows its total subtree contribution
 - **Find the problem**: The root should have spatial coefficients {0,0,0} but doesn't - this is our error
 - **Absorb what you can**: Each bone says "I'll take responsibility for fixing as much error as I can, but only up to my own magnitude and only if the signs match" The sign matching is crucial - you can't fix a +X error by making a bone more -X responsible, that would make things worse.
 - **Pass the buck**: Whatever error a bone can't fix gets distributed to its children based on how much "capacity" each child has to help

The sign matching is crucial - you can't fix a "+X error" by making a bone more "-X responsible", that would make the problem worse. It's like trying to fix a budget deficit by spending more money.

---

## Final Volume Extraction

After error correction, we compute the actual volumes and separate individual bone contributions from subtree totals.

```
// Final step: compute actual volumes after error correction
subtreeVolumes : double[numBones] = {0};
individualBoneVolumes : double[numBones] = {0};
correctedVolumeFunctions = originalVolumeFunctions.copy();

// STEP 1: Apply error corrections and accumulate up the tree
for(boneIndex = numBones-1; boneIndex >= 0; boneIndex--)  // traverse from leaves to root
{
    parentIndex = parents[boneIndex];
    
    // Apply the spatial error correction we computed earlier
    correctedVolumeFunctions[boneIndex] -= backpropCorrection[boneIndex];
    
    // Add this bone's corrected function to its parent (building subtree totals)
    if(parentIndex >= 0)
        correctedVolumeFunctions[parentIndex] += correctedVolumeFunctions[boneIndex];
    
    // Evaluate the volume function at this bone's origin to get actual volume
    // volume = a*x + b*y + c*z + d, where [x,y,z] is the bone origin
    subtreeVolumes[boneIndex] = dot(correctedVolumeFunctions[boneIndex], [origins[boneIndex], 1]);
}

// STEP 2: Extract individual bone volumes from subtree totals
for(boneIndex = numBones-1; boneIndex >= 0; boneIndex--)  // traverse from leaves to root again
{
    parentIndex = parents[boneIndex];
    
    // Individual volume = my subtree volume - what I've already accumulated from children
    individualBoneVolumes[boneIndex] = subtreeVolumes[boneIndex] - individualBoneVolumes[boneIndex];
    
    // Add my subtree volume to my parent's accumulator
    if(parentIndex >= 0)
        individualBoneVolumes[parentIndex] += subtreeVolumes[boneIndex];
}

```

 - **Step 1 is like doing accounting** - we apply our error corrections and then roll up the volume functions so that each bone knows its total subtree volume (itself + all descendants).
 - **Step 2 is the tricky part** - we're doing a "reverse accumulation" to figure out how much volume belongs to each individual bone versus its children. Think of it like this: if bone A has subtree volume 10, and we've already accumulated 7 from its children, then bone A itself contributes 3.

The reason we go backwards (leaves to root) both times is because children need to be processed before their parents - classic dependency ordering for tree algorithms.

Now `individualBoneVolumes[i]` contains the volume that belongs specifically to bone `i`, not including its children.

---

## Extensions

**The beautiful part**: The same backpropagation can be reused for other physical quantities, so long as you can put them in such a way that they vary linearly across the 3-space of the mesh.

**Why this works**: The backprop represents the geometric inconsistency in how your mesh is distributed across the bone hierarchy. This same inconsistency affects any spatially-integrated quantity.

**Other applications**:
- Moments of inertia distribution
- Mass property calculations  
- Center of mass computation
- Any integral over the mesh volume

**Performance notes**: 
- Precompute the backprop once, reuse for multiple physical quantities
- The volume functions are tiny (4 floats per bone), cache-friendly
- No iterations needed - everything is exact linear algebra

**When to use this**:
- You have a skinned mesh with bone weights
- You need volume/mass distribution across bones
- You want mathematically exact results, not approximations
- You're doing physics simulation or procedural animation

**When not to use this**:
- Your mesh isn't skinned (no bone weights available)
- You only care about total volume, not per-bone distribution
- You're working with voxel data (different problem domain)

This algorithm is original work as far as I know - haven't seen this exact linear algebra approach to bone volume distribution in the literature. The mathematical elegance suggests it could be a valuable addition to physics simulation and procedural animation toolkits.
