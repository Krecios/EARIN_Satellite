import math
import easygraphics
import copy
import random
import numpy as np

maxiter = 300  # maximum possible number of iterations
# gravity strength
gravConst = 2
# allowed distance from the goal planet
allowedDistance = 25
# maximum allowed velocity of the satellite
vMax = 10
# number of planets in the simulation
planetsCount = 8
# mass of the satellite
satelliteMass = 1
# parameters of the start and goal planet
startPlanetID = 3
destPlanetID = 2
# parameters for the start of the satellite [v, fi, t_0], filled in problemGenerator()
startParameters = []
# copy of the start parameters used for various calculations, copied in the problemGenerator()
baseParameters = []
# mimal value of the distance between the goal and the satellite, calculated during simulation
minDistance = 0
# array of planets, filled in problemGnerator()
planetsBase = []
# copy of the planet array that is used during the simmulation, copied in the problemGenerator()
planets = copy.deepcopy(planetsBase)
# satellite location, filled by the parameters of start planet in porblemGenerator()
satellite = []


def problemGenerator():
    v = random.randrange(0, vMax * 10, 1)
    v = v/10
    print(v)
    fi = random.randrange(0, 360, 1)
    t0 = random.randrange(0, maxiter, 1)
    global baseParameters
    baseParameters = [v, fi, t0]
    global startParameters
    startParameters = copy.deepcopy(baseParameters)
    nPlanets = planetsCount
    while nPlanets != 0:
        r = random.randrange(10, 400, 1)
        fi = random.randrange(0, 360, 1)
        t = random.randrange(360, 1000, 1)
        planet = [r, fi, t]
        planetsBase.append(planet)
        nPlanets -= 1
    global planets
    planets = copy.deepcopy(planetsBase)
    # start coordinates of the satellite equal to the coordinates of the start planet
    global satellite
    satellite = [planets[startPlanetID][0], planets[startPlanetID][1]]


def rotatePlanets():
    # calculating the rotation in one iteration for each planet
    for plrot in planets:
        plrot[1] = round(plrot[1] + (360 / plrot[2]), 2)
        if plrot[1] > 360:
            # adjusting the degree if it passes 360 degrees
            plrot[1] -= 360


def distanceToDest():
    # goal function
    return math.sqrt(planets[destPlanetID][0] ** 2 + satellite[0] ** 2
                     - 2 * planets[destPlanetID][0] * satellite[0]
                     * math.cos(math.radians(satellite[1] - planets[destPlanetID][1])))


def polarToCart(r, fi):
    # conversion from polar to cartesian coordinates
    x = r * math.cos(math.radians(fi))
    y = r * math.sin(math.radians(fi))
    coords = [x, y]
    return coords


def cartToPolar(x, y):
    # conversion from cartesian to polar coordinates
    r = math.sqrt(x*x + y*y)
    fi = math.degrees(math.atan2(y, x))
    coords = [r, fi]
    return coords


def gravPull(r, fi):
    F_sum = [0, 0]
    for pl in planets:
        # gravitational pull of the planets calculation
        F_pl = gravConst * (1/math.sqrt(r**2 + pl[0]**2 + 2 *
                                        r * pl[0] * math.cos(math.radians(fi - pl[1]))))
        currPl = polarToCart(pl[0], pl[1])
        currSa = polarToCart(r, fi)
        temp = [currPl[0] - currSa[0], currPl[1] - currSa[1]]
        fiDifference = math.degrees(math.atan2(temp[1], temp[0]))
        F_vector = polarToCart(F_pl, fiDifference)
        # print(F_vector)
        F_sum[0] += F_vector[0]/satelliteMass
        F_sum[1] += F_vector[1]/satelliteMass
    # print('Vector:')
    # print('Sum: ', F_sum)
    return F_sum


def positionChange(acceleration, satVelocity, satAngle):
    # updating the satellite position after an iteration
    currPosition = polarToCart(satellite[0], satellite[1])
    velVector = polarToCart(satVelocity, satAngle)

    currPosition[0] = currPosition[0] + velVector[0] + \
        (acceleration[0]*acceleration[0])/2
    currPosition[1] = currPosition[1] + velVector[1] + \
        (acceleration[1]*acceleration[1])/2

    newPosition = cartToPolar(currPosition[0], currPosition[1])
    return newPosition


def velocityChange(acceleration, satVelocity, satAngle):
    # updating the satellite velocity vector after an iteration
    velVector = polarToCart(satVelocity, satAngle)
    velVector[0] += acceleration[0]
    velVector[1] += acceleration[1]
    return velVector


def mainloop():
    # visualization
    x = 500
    y = 400
    global minDistance
    minDistance = distanceToDest()
    vectorCollectionBad = []
    vectorCollectionGood = []
    coordsGoal = []
    iterMax = maxiter
    easygraphics.set_render_mode(easygraphics.RenderMode.RENDER_MANUAL)
    while easygraphics.is_run():
        # code from the simmulation function
        while iterMax != 0:
            rotatePlanets()
            if iterMax <= maxiter - startParameters[2]:
                currDistance = distanceToDest()
                if currDistance < minDistance:
                    minDistance = currDistance
                # calculating the gravitational pull
                pull = gravPull(satellite[0], satellite[1])
                # calculating the new position of the satellite
                newSatPosition = positionChange(pull,
                                                startParameters[0], startParameters[1])
                satellite[0] = newSatPosition[0]
                satellite[1] = newSatPosition[1]
                # print(newSatPosition)
                # calculating the new velocity vector of the satellite
                newSatVelocity = velocityChange(pull,
                                                startParameters[0], startParameters[1])
                # print(newSatVelocity)
                newSatVelocity = cartToPolar(
                    newSatVelocity[0], newSatVelocity[1])
                startParameters[0] = newSatVelocity[0]
                startParameters[1] = newSatVelocity[1]
            else:
                satellite[0] = planets[startPlanetID][0]
                satellite[1] = planets[startPlanetID][1]
            if easygraphics.delay_jfps(40):
                # settin the fps limit
                easygraphics.clear_device()
                for pl in planets:
                    # transforming and drawing planets
                    if pl == planets[destPlanetID]:
                        coordsGoal = polarToCart(pl[0], pl[1])
                        easygraphics.set_fill_color(
                            easygraphics.Color.LIGHT_GREEN)
                        easygraphics.draw_circle(
                            coordsGoal[0] + 500, coordsGoal[1] + 400, allowedDistance)
                    easygraphics.circle(x, y, pl[0])
                    coordsPlanet = polarToCart(pl[0], pl[1])
                    easygraphics.set_fill_color(easygraphics.Color.BLUE)
                    easygraphics.draw_circle(
                        coordsPlanet[0] + 500, coordsPlanet[1] + 400, 5)
                # transforming and drawing satellite
                coordsSatellite = polarToCart(satellite[0], satellite[1])
                goal = polarToCart(
                    planets[destPlanetID][0], planets[destPlanetID][1])
                if math.sqrt((coordsSatellite[0] - goal[0])**2 + (coordsSatellite[1] - goal[1])**2) < allowedDistance:
                    vectorCollectionGood.append(coordsSatellite)
                else:
                    vectorCollectionBad.append(coordsSatellite)
                easygraphics.set_fill_color(easygraphics.Color.RED)
                easygraphics.draw_circle(
                    coordsSatellite[0] + 500, coordsSatellite[1] + 400, 3)
                for vec in vectorCollectionGood:
                    easygraphics.set_fill_color(
                        easygraphics.Color.CYAN)
                    easygraphics.draw_circle(vec[0] + 500, vec[1] + 400, 2)
                for vec in vectorCollectionBad:
                    easygraphics.set_fill_color(
                        easygraphics.Color.RED)
                    easygraphics.draw_circle(vec[0] + 500, vec[1] + 400, 2)
                easygraphics.draw_line(
                    coordsSatellite[0] + 500, coordsSatellite[1] + 400, coordsGoal[0] + 500, coordsGoal[1] + 400)
            iterMax -= 1
    easygraphics.close_graph()


def main():
    easygraphics.init_graph(1000, 800)
    mainloop()
    # easygraphics.close_graph()


def simulation(iterMax):
    global minDistance
    minDistance = distanceToDest()
    while iterMax != 0:
        rotatePlanets()
        if iterMax <= maxiter - startParameters[2]:
            currDistance = distanceToDest()
            if currDistance < minDistance:
                minDistance = currDistance
            # calculating the gravitational pull
            pull = gravPull(satellite[0], satellite[1])
            # calculating the new position of the satellite
            newSatPosition = positionChange(pull,
                                            startParameters[0], startParameters[1])
            satellite[0] = newSatPosition[0]
            satellite[1] = newSatPosition[1]
            # print(newSatPosition)
            # calculating the new velocity vector of the satellite
            newSatVelocity = velocityChange(pull,
                                            startParameters[0], startParameters[1])
            # print(newSatVelocity)
            newSatVelocity = cartToPolar(
                newSatVelocity[0], newSatVelocity[1])
            startParameters[0] = newSatVelocity[0]
            startParameters[1] = newSatVelocity[1]
        else:
            satellite[0] = planets[startPlanetID][0]
            satellite[1] = planets[startPlanetID][1]
        iterMax -= 1
    return minDistance


def reset():
    # resets the system to the initial values
    global planets
    planets = copy.deepcopy(planetsBase)
    global satellite
    satellite = [planets[startPlanetID][0], planets[startPlanetID][1]]


def hillClimbing():
    dest = 10000
    global baseParameters
    global startParameters
    sigmaAngle = 1
    sigmaTime = 1
    sigmaVelocity = 0.2
    iterWithNoChange = 0
    steps = 0
    totalIterations = 0
    while dest > allowedDistance:
        if iterWithNoChange % 5 == 0:
            iterWithNoChange = 0
            sigmaAngle *= 1.05
            sigmaTime *= 1.05
            sigmaVelocity *= 1.05
        startParameters = copy.deepcopy(baseParameters)
        iterDest = simulation(maxiter)
        if dest > iterDest:
            dest = iterDest
            sigmaAngle = 1
            sigmaVelocity = 0.2
            sigmaTime = 1
            steps += 1
            totalIterations += 1
        else:
            baseParameters[0] += np.random.normal(0, sigmaVelocity, 1)
            baseParameters[1] += np.random.normal(0, sigmaAngle, 1)
            baseParameters[2] += np.random.normal(0, sigmaTime, 1)
            if baseParameters[0] < 0:
                baseParameters[0] += vMax
            elif baseParameters[0] > vMax:
                baseParameters[0] -= vMax
            if baseParameters[1] < vMax:
                baseParameters[1] += 360
            elif baseParameters[1] > 360:
                baseParameters[1] -= 360
            if baseParameters[2] > maxiter:
                baseParameters[2] -= maxiter
            elif baseParameters[2] < 0:
                baseParameters[2] += maxiter
            iterWithNoChange += 1
            totalIterations += 1
            print(baseParameters)
        reset()
    startParameters = copy.deepcopy(baseParameters)
    print('Start Parameters:', startParameters)
    print('Minimal Distance:', minDistance)
    print('No. Steps: ', steps)
    print('Total Iterations: ', totalIterations)
    easygraphics.easy_run(main)


problemGenerator()
hillClimbing()
