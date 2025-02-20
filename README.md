# 5-planar
Repository with all python files for my bachelor's thesis looking for bounds of the edge density in 5-planar graphs.


## ILP-solving for optimal and maximal outer 5-planar graphs
We use an ILP-solver to calculate optimal outer *k*-planar graphs of *n* vertices. The central functionality is in [chord_solver.py](./chord_solver.py).
```python
solver = ChordSolver(k=5, n=12)
solution = solver.solve()
solution.chords         # chords in optimal solution
solution.visualize()    # visualization of solution
```

To get maximal solutions under certain constraints, you can prevent specific edges to be in the solution.

```python
solver.remove_edge((0,2))
solution = solver.solve()
```


## Lower bound for optimal simple 5-planar graphs
Functions to calculate the 12-gon chord fillings we used in our lower bound construction og optimal simple 5-planar graphs can be found in [simple_lower_bound.py](./simple_lower_bound.py).
This construction uses 12-gons wrapped around a cylinder. For most of the 12-gons we start with one of the two optimal outer 5-planar graphs, where one of the main diagonals is not crossed by very short chords. Both optimal outer 5-planar 12-gons can be generated using the function ```optimal_dodecagons()```.

For the top and bottom 12-gon we have a couple of possible multiedges, that we need to prevent. The possible multiedges, that appear in the best possible configuration, are listed in ```top_multiedges``` and ```bottom_multiedges```. Note that this naming does not implie, that the correspoding 12-gon always is at the top or the bottom, both can be both.

With the function ```simple_dodecagon(multiedges)``` chord solutions can be found, which does not contain any edge given in the ```multiedges``` list.
