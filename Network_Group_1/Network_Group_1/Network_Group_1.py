import pandas as pd
import math
import matplotlib.pyplot as plt
from gurobipy import *
import numpy as np
from mpl_toolkits.basemap import Basemap
import networkx as nx



data = pd.read_csv(r'C:\Users\USER\Documents\Imperial College London\Core Module\Network Analytics\Assignment\A2\tempData.txt',  delimiter=r"\s+")

def distance_on_unit_sphere(lat1, long1, lat2, long2):
 
# Convert latitude and longitude to
# spherical coordinates in radians.
    degrees_to_radians = math.pi / 180.0
 
    # phi = 90 - latitude
    phi1 = (90.0 - lat1) * degrees_to_radians
    phi2 = (90.0 - lat2) * degrees_to_radians
 
    # theta = longitude
    theta1 = long1 * degrees_to_radians
    theta2 = long2 * degrees_to_radians
 
    # Compute spherical distance from spherical coordinates.
 
    # For two locations in spherical coordinates
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) =
    # sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
 
    cos = (math.sin(phi1) * math.sin(phi2) * math.cos(theta1 - theta2) + math.cos(phi1) * math.cos(phi2))
    arc = math.acos(cos)
 
    # Remember to multiply arc by the radius of the earth
    # in your favorite set of units to get length.
    return arc



def distance_table_transform(data):
    transData = {'From': '', 
             'To': '',
             'Distance': ''
             }
    fList = []
    toList = []
    dis = []
    for i in arange(len(data)):
        for n in arange(len(data)):
            if (i != n):
                fList.append(i + 1)
                toList.append(n + 1)
                dis.append((6373 *distance_on_unit_sphere((data['latitude'][i] / 1000), (data['longitude'][i] / 1000), (data['latitude'][n] / 1000), (data['longitude'][n] / 1000))))
            else:
                fList.append(i + 1)
                toList.append(n + 1)
                dis.append(0)
    transData['From'] = fList
    transData['To'] = toList
    transData['Distance'] = dis
    return transData




def distance_matrix_transform(transData):
    ### Distance Matrix
    Pdtrans = pd.DataFrame.from_dict(transData)
    DistanceMatrix = pd.pivot_table(Pdtrans, index = ['From'], columns = ['To'])
    pd.DataFrame.to_csv(DistanceMatrix, r'C:\Users\USER\Documents\Imperial College London\Core Module\Network Analytics\Assignment\A2\DistanceMatrix.csv')
    return DistanceMatrix

def draw_map_scatter(data , mapscale=0.5, markersize=5):
    lat = []
    long = []
    for x in data['latitude']: lat.append(x / 1000)
    for x in data['longitude']: long.append(x / 1000)
    
    minlat = min(lat) - mapscale
    maxlat = max(lat) + mapscale
    minlong = min(long) - mapscale
    maxlong = max(long) + mapscale
    meanlat = mean(lat)
    meanlong = mean(long)
    
    plt.figure(figsize= (20,20))
    map = Basemap(projection='lcc', 
                  llcrnrlon = minlong ,
                  llcrnrlat = minlat,
                  urcrnrlon = maxlong,
                  urcrnrlat = maxlat ,
                  lat_0 = meanlat,
                  lon_0 = meanlong,
                  resolution = 'h')
    x,y = map(long,lat)
    map.fillcontinents(color = 'coral', lake_color ='blue')
    map.drawcoastlines(linewidth = 4)
    map.drawcountries(linewidth = 2)
    map.plot(x,y, 'ro', markersize= markersize)
    map.drawlsmask(ocean_color= 'blue',land_color='blue')
    plt.show()


transData = distance_table_transform(data)
DistanceMatrix = distance_matrix_transform(transData)
draw_map_scatter(data,0.2, 7)


### 2c


# Optimize model

def model_setup ():
    m = Model()
    #m.params.LazyConstraints = 1
    m.setObjective(GRB.MINIMIZE)
    dmc = {}
    for i in arange(38):
        for j in arange(38):
            dmc[i,j] = float(DistanceMatrix.loc[i + 1][j])
    vars = {}
    for i in range(38):
        for j in range(i + 1):
            vars[i,j] = m.addVar(obj=dmc[i,j], vtype=GRB.BINARY,
                                name='e' + str(i) + '_' + str(j))
            vars[j,i] = vars[i,j]
        m.update()
    for i in range(38):
        m.addConstr(quicksum(vars[i,j] for j in range(38)) == 2)
        vars[i,i].ub = 0
    m.update()
    m._vars = vars
    res = {}
    return m

def optimize_tour (n = 38):
    m = model_setup()
    while True:
        m.reset()
        m.optimize()
        selected = {}
        for i in arange(38):
            for j in arange(i+1):
                if m._vars[i,j].X >=0.5:
                    selected[i,j] = (i,j)
        G = nx.Graph()
        G.add_edges_from(selected)
        cycle = nx.cycle_basis(G)
        if len(cycle) == 1 and len(cycle[0])==n:
            m.optimize()
            res = selected
            break
        else:
            print(sort(cycle[-1]))
            if len(cycle) ==4:
                m.addConstr(quicksum(m._vars[k,l] for k in sort(cycle[0]) for l in sort(cycle[1]+ cycle [2]+ cycle[3])) >= 2)
            elif len(cycle) == 3:
                m.addConstr(quicksum(m._vars[k,l] for k in sort(cycle[0]) for l in sort(cycle[1]+ cycle [2])) >= 2)
            else:
                m.addConstr(quicksum(m._vars[k,l] for k in sort(cycle[0]) for l in sort(cycle[-1])) >= 2)
    resnorm =[]
    for i in res.keys():
        resnorm.append(i)
    return resnorm

def optimize_tour_map(selected, DistanceMatrix, data,mapscale = 0.5, res = 'f'):
    udv = []
    for i in range(0,len(selected)):
        udv.append((selected[i][0] +1,selected[i][1] +1,DistanceMatrix.loc[selected[i][0] + 1][selected[i][1]] ))
    
    g =nx.Graph()
    node = arange(1,39)
    g.add_nodes_from(node)
    g.add_weighted_edges_from(udv)
    lat = []
    long = []
    for x in data['latitude']: lat.append(x / 1000)
    for x in data['longitude']: long.append(x / 1000)
    minlat = min(lat) - mapscale
    maxlat = max(lat) + mapscale
    minlong = min(long) - mapscale
    maxlong = max(long) + mapscale
    meanlat = mean(lat)
    meanlong = mean(long)
    position = {}

    
    plt.figure(figsize= (20,20))
    map = Basemap(projection='lcc', 
                    llcrnrlon = minlong ,
                    llcrnrlat = minlat,
                    urcrnrlon = maxlong,
                    urcrnrlat = maxlat ,
                    lat_0 = meanlat,
                    lon_0 = meanlong,
                    resolution = res
                    )
    for i in range(0,len(data)):
        position[i] = map(data['longitude'][i]/1000,data['latitude'][i]/1000)
    pos = {}
    for i in range(1,len(data)+1):
        pos[i] = array([position[i-1][0],position[i-1][1]])

    #map.fillcontinents(color = 'coral', alpha = .0)
    map.drawcoastlines(linewidth = 4)
    map.drawcountries(linewidth = 2)
    map.bluemarble()
    #map.drawlsmask(ocean_color= 'blue',land_color='blue')
    nx.draw_networkx_nodes(g, pos = pos, node_size = 50)
    nx.draw_networkx_edges(g, pos= pos, edge_color='lightblue', width = 5, alpha = 1)
    plt.show()




res = optimize_tour()
optimize_tour_map(res, DistanceMatrix, data,mapscale = 0.5, res = '')




