# Inverse Kinematics Guide

0. [Terms](#0--terms)
1. [Advice](#1--advice)
2. [Fabrik Solver](#2--fabrik-solver)  (forwards and backwards reaching inverse kinematics)
3. [Cyclic Coordinate Descent Solver](#3--cyclic-coordinate-descent-solver)
4. [Jacobian Solver](#4--jacobian-solver)
5. [Mass Spring Solver](#5--mass-spring-solver)
6. [Scale Target Solver](#6--scale-target-solver)

# 0.  Terms

- Effector: this is the part that tries to hit the goal, your hand is the effector that gets the soda bottle.
- Nodes: the roots of the bones in the armature, make sure you have a node for the last tip too, such as by adding a leaf bone before exporting
- normalize: converting a vector to unit length so vec / length(vec)
- path independent: the pose of the armature is not dependent on the previous frame, there is no chain of dependencies
	this ensures that the joints can't get twisted into knots.
- path dependent: the pose is dependent on the previous frame.
	this is generally faster because we need less iterations to converage (because the previous frame is a good hint)
- swizzling: rearranging variables; xyz is a swizzle of yxz
- colinear: when 3 points fall on a line they are colinear; when i say a joint is colinear i mean that the parent of the joint, the joint, and the immediate child of the joint are colinear; the joint is said to be colinear because changing it's angle is what will change this.
- DFS a Depth First Search ordering is the ordering of nodes by exploring as far as possible before backtracking.

# 1.  Advice

- When using IK on an armature; always preprocess the armature so that the nodes are in strict DFS order; this has better caching behavior and ends up simplyifing the code and algorithms a lot more than you would think. I would also advise making the spine the first complete chain in the DFS ordering.

- By extension: sort the effector goals so they apply to the skeleton such that the effectors are in DFS order. For example; say the head is the before the hand in the ordering, and you mark nodes so that effectors further down the line don't affect earlier nodes; the "look" target of the head will bend the spine in a natural way while the arm target reaches for the goal; conversely if you put the arm target first then the spine will twist into the motion and create an overextended robotic look.

- When working with an algorithm that makes angle adjustments, store the angle that would make the joint colinear with it's parent and bump the result if it is equal to this value, don't let it become colinear, this is a degenerate state the solver will have issues with.

- Make sure the rest position of the armature is naturalistic; for example, for a human you want it to be closer to riding a bike than a T-pose; this is because the parts are all bent in a way that hints at the primary way they move; so a simple unconstrained solver is much more likely to get a plausible result, than if the starting position has straight arms and legs.

- use libraries for things like computing the quat between two vectors or matrix inverses; its not worth the trouble.  i recommend GLM or eigen.

- All solver algorithms are sub polynomial if memoized correctly; the speed is not that much of a big deal as compared to the cache optimization.

- If you have scale targets solve them first instead of dealing with them at the same time as position/rotation constraints; as they are less inter-related.



# 2.  Fabrik Solver

### Overview:

The fabrik solver is fast, lightweight and easy for beginners. It excels at situations where there are many targets and many effectors, such as if a few effectors have many targets, or if there are effectors in the middle of the chain etc. It has issues with situations where there are rotation constraints; this is fairly doable in 2D environments, but in 3D the fact that fabrik has no concept of roll means that its not a good choice, and trying to get around that is more trouble than it's worth.


### Let:

	// total joints, 0 is the root, N-1 is the effector, so the chain goes 0-1-2-3 etc.
	  N : int;
	// position of the node in world space
	  positions: vec3[];
	// rotation of the node in world space
	  rotations: quat[];
	// the length of each bone so length[i] = length(position[i+1] - positon[i])
	  lengths : float[];

	// the target we are reaching for; for multiple targets i would suggest using the average
	// position, its easy and works well.
	  goal: vec3;

	// to compute the final rotations we always need the initial transform, regardless of path dependence
	// if we use path dependent rotations then errors accumulate and the solver goes nuts.
	// the initial values MUST be in local space
	  initial_rotations: quat[];
	  initial_positions: vec3[];

### Algorithm:

	// store the position of the root we'll need it later!
	  root := positions[0];
	// set the effector position to the target we wanna reach
	  positions[N-1] = goal;

	// backwards part comes first
	// here we move each joint closer to the effector by adjusting according to the length constraint.
	  for(i : N-2..0)
	  {
	// note that if there is no change then this simplifies to: a + (b - a)
	// you can see the two "a" values cancel and we just get b.  this is what we want
	// when putting things into the normalize function make sure its ordered so it simplifies like this
	// we solve for the position of the parent given the position of the child.
		positons[i] = positions[i+1] + normalize(positions[i] - positions[i+1]) * lengths[i];
	  }

	// forwards part comes net!
	// here we move the root back into place and adjust as we did before
	// the result will be that all the joints bent a bit to get closer to the target
	// multiple iterations are needed for a path independent solver,
	// but for a path dependent solver just do it once, the previous frame works well enough.
	  positions[0] = root;

	  for(i : 1..N)
	  {
		// same as before; we want to simplify it so that if the length changes a cancels: a + (b - a)
		// but we solve for the position of the child given the parent.
		positions[i] = positions[i-1] + normalize(positions[i] - positions[i-1]) * lengths[i-1];
	  }

	// now we need the rotations; we need a child to compute the rotation
	// so the last node can't be computed
	  for(i : 1..N-1)
	  {
	// get the vector of the default transform, we want this normalized so that we can compute rotations.
	// remember the initial values are in local space not world space like the finals!!
		initialVector := initial_positions[i] / lengths[i];
		finalVector  := (positons[i] - positions[i-1]) / lengths[i];

	// the initial vector is in local space so convert the final vector to local space.
	// conjugation is fast, really fast.  conjugate(x, y, z, w)  = {-x, -y, -z, w}
		finalVector =  conjugate(rotations[i-1]) * finalVector;

	// get the rotation needed to convert the initial vector to the final vector
		deltaRotation  := GetRotationBetween(initialVector, finalVector);

		rotations[i] = rotations[i-1] * deltaRotation;
	  }


# 3.  Cyclic Coordinate Descent Solver

### Overview:

The Cyclic Coordinate Descent solver is a solver based on considering each joint in the chain individually and solving it analitically.  Because it operates on arrays of angles and so does the jacobian solver; the output of this solver can be used as input to the jacobian solver and vice versa for refinement.

This is an A-tier solver, its a good all-rounder that produces realistic-ish motion and is fast enough for real time.  This was the solver used for foot placement in many games like Wind Waker!


### Let:

	  struct transform { quat rotation; vec3 position; };

	// total joints, 0 is the root, N-1 is the effector, so the chain goes 0-1-2-3 etc.
	  N : int;
	// the target we are reaching for; for multiple targets run the solver for each target
	// individually and average the output angles
	  goal : vec3;
	// the angles of the joints in the solver, this is for simplicity, translations would
	// also work in this algorithm
	  angles : float[]:

	// local space transforms
	  localSpace : transform[];
	// get the transform of the joint
	// joints are considered to be in between local space transforms, so on the local axises
	  getJoint  : function(node, angles[], first, last) -> rotation;
	// world space transforms from 0 to i, so:
	// _0ToI[i] = localSpace[0] * joint[0] * ... * joint[i-1] * localSpace[i]
	  _0ToI : transform[];
	// world space transforms from i to N, so:
	// _iToN[i] = localSpace[i] * joint[i] * ... * joint[N-1] * localSpace[N]
	  _iToN : transform[];

	// accumulate this and return it, if its too high we need to iterate again
	  lambda: float = 0;

### Algorithm:

	  GetTheta : function = (effector_local_space : vec2, target_local_space : vec2, _default : float) ->
	  {
		length_effector := length(effector);
		length_target 	:= length(target);
		length 			:= length_effector*length_target;

		if(!length)
			return _default;
		else
		{
	// dot product is a cosine projection: a.x * b.x + a.y * b.y
			cosine 	:= dot(effector, target) / length;
	// 2d cross product is a sine projection: a.x * b.y - a.y * b.x
			sine 	:= cross(effector, target) / length;

			return atan2(sine, cosine);
		}
	  }


	  for(node : N-1..0)
	  {
		for(axis : 0..2)
		{
	// we're only considering one axis at a time; so we want to include the joint of I in the computation
	// so the target local space needs to consider all joints preceeding this one
			target := inverse(_0ToI[node] * getJoint(i, angles, 0, axis-1)) * goal;
	// and the effector local space needs to consider all joints after this one!
			effector = getJoint(node, angles, axis+1, 2) * _iToN[node].positition;

	// for each axis get the theta by projecting onto the 2D plane;
	// this can be done just by swizzling because we're already in local space.
			if(axis == 0)
				angles[node] = GetTheta(effector.yz, target.yz, angles[node*3+axis]);
			if(axis == 1)
				angles[node] = GetTheta(effector.xz, target.xz, angles[node*3+axis]);
			if(axis == 2)
				angles[node] = GetTheta(effector.xy, target.xy, angles[node*3+axis]);

	// if you need to apply min/max constraints do it here
	// also make sure that angles[node] isn't a value that can cause gimble lock or become colinear
	// if this happens add a small value to it.
		}

	// update our memo with the new orientation
		_0ToI[node+1] = _0ToI[node] * getJoint(node, angles, 0, 2) * localSpace[node+1];
	  }


# 4.  Jacobian solver

The jacobian algorithm is a C tier solver that is fairly slow but also produces movement that isn't very life-like.  Rather, it will move all bones in the chain equally to try to reach the target; which is not how animals move.

The jacobian solver is extremely hard to get your head around because the language used to describe it is very obtuse, if you didn't major in math its all basically meaningless terminology.  However the algorithm itself isn't hard. lets start by defining a jacobian matrix.

	struct JacobianMatrix
	{
		X : float[];
		Y : float[];
		Z : float[];
	};

In a jacobian matrix each column corresponds to a value for an axis, maybe a delta, maybe a tangent, but a value. And the rows (indexes) correspond to joints (more on this later).  This is basically worthless for computation, you can't do anything useful with this data structure, its junk.  It's arranged like this because of the field of math it came from.

Next we have the jacobian transpose which is this:

	JacobianTranspose : vec3[];

Here the columns are the joints, and the rows are XYZ values; all I've done is turn it from an struct of arrays to an array of structs. Thats what taking a "transpose" means. So while the mathematical papers and algorithm descriptions talk about computing the jacobian and taking the transpose. You never do that, never ever do that: it is far easier to just compute the tranpose directly and work with it.

The jacobian solver has basically 4 functions that need to be explained and defined they are as follows:

a. computing the jacobian transpose matrix
b. multiplying a vector by the jacobian transpose
c. computing the jacobian times the jacobian tranpose (J * J^T)
d. the solver itself.


## 4.a.  Computing the Jacobian Transpose Matrix

As before we consider the jacobian transpose matrix as an array of vec3 objects; each vec3 describes how the end effector will move in response to changes at the current joint. If you've taken calculus before think of it as a secant, or a poor estimate of a tangent.

So the method here is to alter the joint by a small amount, measure what happened and record it. When people say the jacobian is O(N^2) they mean that they need to recompute the armature arm to get how the end effector changed; fact is, you don't. you can just memoize cleverly.


### Let:

	  struct transform { quat rotation; vec3 position; };

	// total joints, 0 is the root, N-1 is the effector, so the chain goes 0-1-2-3 etc.
	  noJoints : int;
	// the angles of the joints in the solver, this is for simplicity,
	// translations would also work in this algorithm
	  angles : float[]:

	// get the transform of the joint
	// joints are considered to be in between local space transforms, so on the local axises
	// if last is overindexed it will be clamped
	  getJoint  : function(node, angles[], first, last = ~0u) -> rotation;

	// total nodes in the armature
	  noNodes  : int;
	// local space transforms
	  localSpace : transform[];

	// get nodeId from joint Id
	  nodeFromJoint : int[];


	// world space transforms from 0 to i, so:
	// _0ToI[i] = localSpace[0] * joint[0] * ... * joint[i-1] * localSpace[i]
	  _0ToI : transform[];
	// world space transforms from i to N, so:
	// _iToN[i] = localSpace[i] * joint[i] * ... * joint[N-1] * localSpace[N]
	  _iToN : transform[];

	// small value to alter joint by to determine what it will do
	  smallValue := 0.01;

	// output, this is the transpose
	  output : vec3[];

### Algorithm:

	  vec3 originalEffectorPos = _0ToI[noNodes].translation;

	  for(i : 0..noJoints-1)
	  {
	// change the value by a small amount
		currentAngle := angles[i];
		angles[i] += smallValue;

	// get new joint transform
		joint : transform = getJoint(bones[i], angles);

	// restore value we don't want to actually alter the angles in this function.
		angles[i] = currentAngle;

	// recompute the effector position
		effector := _0ToI[bones[i]] * joint * _iToN[bones[i]+1].translation;

	// get the secent and store it
		output[i] = (effector - originalEffectorPos) / smallValue;
	  }


## 4.b.  Multiplying a Vector By The Jacobian Transpose

The important thing to understand here is that :

	vec3 = jacobian transpose * vec3

This is unintuitive because the jacobian method adjusts each angle by the same amount as a result of the solver; so if elbows bend on the X axis all joints that can bend on the X axis will; which is kind of odd; it seems like we should get an array of values out of this, one for each joint, but we don't, we get a vec3, a value for each axis that will apply to all joints.

### Algorithm:

	MulByTranspose(transpose : vec3[], vec3 input) -> vec3
		vec3 result = vec3(0);

		for(item : transpose)
			result += input * item;

		return result;

## 4.c.  Computing The Jacobian Times The Jacobian Tranpose (J * J^T)

In most cases matrix multiplications cache really badly and are inherently O(N^2) but because we're multiplying a matrix by it's own transpose we can take some shortcuts.

Now the jacobian matrix has N rows and 3 columns; and the transpose has 3 rows and N columns; so when we multiply them we always get a 3x3 matrix.

Each element is a sum of products, for example an element may represent a sum like:

	shoulder.x * shoulder.x + elbow.x * elbow.x + wrist.x * wrist.x

This represents the cumulative effect of all joints on a particular axis of movement. This captures the interdependencies between the joints, and estimate how collective movements will affect the end effector.

### Let:

	// total joints
	  N : int;
	// jacobian transpose computed above
	  D: vec3[];

	// axis for each joint, 0 = x, 1 = y, 2 = z
	  axis : int[];

### Algorithm version 1:

	result := mat3(0);

	for(i : 0..N-1)
		for(j : 0..N-1)
			result[axis[i]] += D[i] * D[j];

### Algorithm version 2:

	// notice that everything in result [axis[i]] was multiplied by D[i]
	// so we can use the distributive property: a*b + a*c = a * (b + c)
	// (this is something we can do because its a transpose multiplication)
	  result := mat3(0);
	  accumulator := vec3(0);

	  for(i : 0..N-1)
		accumulator += D[i];

	  for(i : 0..N-1)
		result[axis[i]] += D[i] * accumulator;

This is a really good candidate for optimization with SIMD; the compiler will not do this for you because the byte alignment of the input array is not correct. Check the documentation for your language on SIMD and __m128 registers to correctly write this method.  SIMD can be tricky though!

## 4.d.  The jacobian solver

### Let:

	// square of the lambda, this is used by the damped least squares method to better estimate the result
	  lambda2 : float;

	// world space transforms from 0 to i, so:
	// _0ToI[i] = localSpace[0] * joint[0] * ... * joint[i-1] * localSpace[i]
	  _0ToI : transform[];

	  angles : float[]; // current state of the joints
	  axis : int[]; // axis each joint locally roates on
	  min : float[]; // min angle of each joint
	  max : float[]; // max angle of each joint

### Algorithm:

	  error := goal - _0ToI[noNodes].translation;
	  transpose := ComputeTranspose();

	// DLS jacobian solver (damped least squares)
	// this is an optional step that improves solver quality
	// when is it useful?
	// if the matrix is ill conditioned, meaning that there are lots of potential
	// ways to get the end effector to the desired location, or small changes in
	// joint angles can produce large changes in output.

	// When this happens we regularize the system by dampening
	  if(lambda2 > 0)
	  {
	// this step captures the cumulative effects of all joint movements, it
	// essentially forms a basis for understanding how collective joint adjustment
	// will affect the effector.
		JJT : mat3 = ComputeJacobianJacobianTranspose(transpose);

	// next the dampening factor is added to the to the diagonal term before caclulating its inverse.
	// this essentially smooths adjustments leading to more gradual changes in joint positions

	// what?

	// when we increase the diagonal dominance of the matrix it makes it less sensitive
	// to small varations and more numerically stable.
	// this may be called "reducing the condition number"

	// so this creates a tradeoff between accuracy and stability, higher lambda
	// increases stability (less prone to oscilation), but can lead to less
	// accurate tracking by the end effector.  Conversely a lower lambda improves
	// accuracy but can lead to oscilation and overcorrection.
		inverse : mat3 = inverse(JJT + mat3(lambda2));

	// dampen the error term to increase stability.
		error = inverse * error;
	  }

	// here we are effectively distributing the error backwards through the joints
	// mapping the distance from the effector to the goal into joint coordinate space.
	  adjustment := MulByTranspose(error);

	// the adjustment of say adjustment.x applies to all joints that rotate
	// along the X axis this is why the jacobian method leads to unnaturalistic
	// movements; a real animal will say move it's elbow to grab a bean, not
	// move its shoulder, elbow, wrist, and spine equally to grab the bean.

	// this also means that if we have a situation like a dog leg; where we
	// have two hinge joints that bend in opposite directions on the same axis
	// the jacobian solver will produce unnatural movements.
	  for(i : 0..N-1)
		angles[i] = clamp(angles[i] + adjustment[axis[i]], min[i], max[i]);



# 5.  Mass Spring Solver

The spring mass solver will produce a bouncy effect thats good for hair physics and squash and stretch. Don't actually use hook's law to implement this, because the values of K that are stable depend on the value of deltaTime; meaning that inconsistent framerates can lead to blow up; and its not a nice linear function mapping max K to deltaTime that makes for easy clamping either.

### Let:

	// total joints, 0 is the root, N-1 is the effector, so the chain goes 0-1-2-3 etc.
	  N : int;
	// position of the node in world space
	  positions: vec3[];
	// previous position of the node in world space
	  prev_positions: vec3[];
	// mass of each joint
	  mass: float[];

	// rotation of the node in world space
	  rotations: quat[];
	// the length of each bone so length[i] = length(position[i+1] - positon[i])
	  lengths : float[];

	// to compute the final rotations we always need the initial transform, regardless of path dependence
	// if we use path dependent rotations then errors accumulate and the solver goes nuts.
	// the initial values MUST be in local space
	  initial_rotations: quat[];
	  initial_positions: vec3[];

### Algorithm:

	// start at 1 because we assume node 0 to be an anchor.
	  for(i : 1..N-1)
	  {
	// use difference between current and previous positions to approximate velocity;
		prev_positions[i] = position[i] + (position[i] - prev_positions[i]);
		prev_positions[i] += 0.5 * gravity * dt * dt;
	// also add in other forces like collision, bouyancy etc.
	  }
	// swap the two for the next frame
	  prev_positions, positions = positions, prev_positions;

	  for(iteration : noIterations)
	  {
		for(i : 1..N-1)
		{
			invMass2 := 1.0 / (mass[i] * mass[i-1]);
	// get weighted average position of the spring
			center := (positions[i] * mass[i] + positions[i-1] * mass[i-1] ) * invMass2;

	// note that if there is no change then this simplifies to: a + (b - a)
	// you can see the two "a" values cancel and we just get b.  this is what we want
	// when putting things into the normalize function make sure its ordered so it simplifies like this
	// we solve for the position of the ends given the center and lengths
			position[i] = center + normalize(positions[i  ] - center) * lengths[i] * mass[i  ] * invMass2;
			position[i] = center + normalize(positions[i-1] - center) * lengths[i] * mass[i-1] * invMass2;
		}
	  }

	// now we need the rotations; we need a child to compute the rotation
	// so the last node can't be computed
	  for(i : 1..N-1)
	  {

	// get the vector of the default transform, we want this normalized so that we can compute rotations.
	// remember the initial values are in local space not world space like the finals!!
		initialVector := initial_positions[i] / lengths[i];
	// this line is different than in the fabrik rotation finder.
		finalVector  := normalize(positons[i] - positions[i-1]);

	// the initial vector is in local space so convert the final vector to local space.
	// conjugation is fast, really fast.  conjugate(x, y, z, w)  = {-x, -y, -z, w}
		finalVector =  conjugate(rotations[i-1]) * finalVector;

	// get the rotation needed to convert the initial vector to the final vector
		deltaRotation  := GetRotationBetween(initialVector, finalVector);

		rotations[i] = rotations[i-1] * deltaRotation;
	  }

# 6.  Scale Target Solver

The scale target solver is used to distribute scaling information from a goal into the armature.  A very simplified overview would be as follows:


	distribute := (X : float, bucket[], size[])
	{
		for( i : 0..bucket.size-1)
		{
			bucket[i] = min(size[i], X);
			X /= bucket[i];
		}
	}

The unconstrained algorithm will work with negative scales, but the constrained one doesn't.

### Let:

	// number of nodes
	  N : int;
	// scales in armature chain
	  scales : vec3[];
	// goal scale of the effector
	  scaleGoal : vec3;

	// scale constraints
	  minScale : vec3[];
	  maxScale : vec3[];

### Algorithm:

	// solve with constraints
	  if(minScale.size == maxScale.size || minScale.size == scales.size)
	  {
		amountToDistribute := pow(scaleGoal, 1.0 / max(N, 1));
		remaining := scaleGoal;

	// count up from root to leaf so that redistributed scale gets focused towards the effector
		for(i : 0..N-1)
		{
			scales[i] = clamp(amountToDistribute, minScale[i], maxScale[i]);
			remaining /= scales[i];

	// redistribute to put more in next mnode
			if(scales[i] != amountToDistribute)
				amountToDistribute := pow(remaining, 1.0 / max((N-i), 1));
		}
	  }
	  else if(scale.x > 0 && scale.y > 0 && scale.z > 0)
	  {
		amountToDistribute := pow(scaleGoal, 1.0 / max(N, 1));

		for(i : 0..N-1)
			scales[i] = amountToDistribute;
	  }
	  else
	  {
		for(i : 1..N-1)
		{
	// set global space
			scales[i] = lerp(scales[0], scaleGoal, i / float(N));
	// convert to local space
			scales[i] /= lerp(scales[0], scaleGoal, (i+1) / float(N));
	// check for div by 0 etc...
		}
	  }

