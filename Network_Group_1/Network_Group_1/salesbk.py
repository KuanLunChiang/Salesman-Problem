
dmc = {}
for i in arange(38):
    for j in arange(38):
        dmc[i,j] = float(DistanceMatrix.loc[i + 1][j])


def subtour(edges):
  visited = [False] * n
  cycles = []
  lengths = []
  selected = [[] for i in range(n)]
  for x,y in edges:
    selected[x].append(y)
  while True:
    current = visited.index(False)
    thiscycle = [current]
    while True:
      visited[current] = True
      neighbors = [x for x in selected[current] if not visited[x]]
      if len(neighbors) == 0:
        break
      current = neighbors[0]
      thiscycle.append(current)
    cycles.append(thiscycle)
    lengths.append(len(thiscycle))
    if sum(lengths) == n:
      break
  print(cycles[lengths.index(min(lengths))])
  return cycles[lengths.index(min(lengths))]

def subtourelim(model, where):
  if where == GRB.callback.MIPSOL:
    selected = []
    # make a list of edges selected in the solution
    for i in range(n):
      sol = model.cbGetSolution([model._vars[i,j] for j in range(n)])
      selected += [(i,j) for j in range(n) if sol[j] > 0.5]
    # find the shortest cycle in the selected edge list
    tour = subtour(selected)
    if len(tour) < n:
      # add a subtour elimination constraint
      expr = 0
      for i in range(len(tour)):
        for j in range(i + 1, len(tour)):
          expr += model._vars[tour[i], tour[j]]
      model.cbLazy(expr <= len(tour) - 1)


n = 38

# Create n random points
m = Model()
# Create variables
vars = {}
for i in range(n):
   for j in range(i + 1):
     vars[i,j] = m.addVar(obj=dmc[i,j], vtype=GRB.BINARY,
                          name='e' + str(i) + '_' + str(j))
     vars[j,i] = vars[i,j]
   m.update()


# Add degree-2 constraint, and forbid loops
for i in range(n):
  m.addConstr(quicksum(vars[i,j] for j in range(n)) == 2)
  vars[i,i].ub = 0

m.update()

m._vars = vars
m.params.LazyConstraints = 1
m.optimize(subtourelim)
solution = m.getAttr('x', vars)
selected = [(i,j) for i in range(n) for j in range(n) if solution[i,j] > 0.5]
assert len(subtour(selected)) == n

selected




x = {(1, 0): (1, 0), (3, 1): (3, 1), (3, 2): (3, 2), (4, 2): (4, 2), (5, 4): (5, 4), (6, 5): (6, 5), (7, 6): (7, 6), (8, 7): (8, 7), (9, 0): (9, 0), (11, 10): (11, 10), (13, 9): (13, 9), (14, 12): (14, 12), (15, 11): (15, 11), (16, 8): (16, 8), (16, 10): (16, 10), (17, 15): (17, 15), (18, 12): (18, 12), (18, 17): (18, 17), (19, 14): (19, 14), (20, 13): (20, 13), (22, 19): (22, 19), (23, 21): (23, 21), (24, 21): (24, 21), (25, 22): (25, 22), (25, 24): (25, 24), (27, 23): (27, 23), (27, 26): (27, 26), (28, 20): (28, 20), (29, 28): (29, 28), (30, 26): (30, 26), (31, 29): (31, 29), (33, 32): (33, 32), (34, 31): (34, 31), (35, 30): (35, 30), (35, 33): (35, 33), (36, 34): (36, 34), (37, 32): (37, 32), (37, 36): (37, 36)}
resnorm =[]
for i in x.keys():
    resnorm.append(i)

n=38

m = Model()
m.setObjective(GRB.NONBASIC_LOWER)
dmc = {}
for i in arange(38):
    for j in arange(38):
        dmc[i,j] = float(DistanceMatrix.loc[i + 1][j])
vars = {}
for i in range(n):
    for j in range(i + 1):
        vars[i,j] = m.addVar(obj=dmc[i,j], vtype=GRB.BINARY,
                            name='e' + str(i) + '_' + str(j))
        vars[j,i] = vars[i,j]
    m.update()
for i in range(n):
    m.addConstr(quicksum(vars[i,j] for j in range(n)) == 2)
    vars[i,i].ub = 0
m.update()
m._vars = vars
m.params.LazyConstraints = 0
res = {}
while True:
    m.optimize()
    selected = {}
    for i in arange(38):
        for j in arange(i+1):
            if vars[i,j].X == 1:
                selected[i,j] = (i,j)
    G = nx.Graph()
    G.add_edges_from(selected)
    cycle = nx.cycle_basis(G)
    print(selected)
    if len(cycle) == 1 and len(cycle[0])==n:
        res = selected
        break
    else:
        m.addConstr(quicksum(vars[l,k] for k in cycle[0] for l in cycle[-1]) >= 2)
resnorm =[]
for i in res.keys():
    resnorm.append(i)

m.optimize()