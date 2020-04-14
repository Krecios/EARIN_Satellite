import math
maxiter = 2  # maximum possible number of iterations
# temporary attributes of the satellite for testing
gravConst = 5
# parameters for the start of the satellite [v, fi, t_0]
startParameters = [0.5, 15, 3]
satelliteMass = 0.00001
startPanetID = 2
destPlanetID = 1
planets = []  # array of planets
pl1 = [30, 35, 500]  # planet [r, fi, p]
pl2 = [50, 87, 700]
pl3 = [90, 123, 850]
planets.append(pl1)
planets.append(pl2)
planets.append(pl3)
satellite = [planets[startPanetID][0], planets[startPanetID][1]]


def rotatePlanets():
    # calculating the rotation in one iteration for each planet
    for plrot in planets:
        plrot[1] = round(plrot[1] + (plrot[2] / 360), 2)
        if plrot[1] > 360:
            # adjusting the degree if it passes 360 degrees
            plrot[1] -= 360
        print(plrot[1])


def distanceToDest():
    print("WIP")


def polarToCart(r, fi):
    x = r * math.cos(math.radians(fi))
    y = r * math.sin(math.radians(fi))
    coords = [x, y]
    return coords


def cartToPolar(x, y):
    r = math.sqrt(x*x + y*y)
    fi = math.degrees(math.atan2(y, x))
    coords = [r, fi]
    return coords


def gravPull(r, fi):
    F_sum = [0, 0]
    for pl in planets:
        # grav pull calculation
        F_pl = gravConst * \
            (1/math.sqrt(r*r + pl[0]*pl[0] + 2 *
                         r * pl[0] * math.cos(fi - pl[1])))
        fiDifference = fi - pl[1]
        F_vector = polarToCart(F_pl, fiDifference)
        # print(F_vector)
        F_sum[0] += F_vector[0]
        F_sum[1] += F_vector[1]
    print('Vector:')
    return F_sum


def positionChange(acceleration, satVelocity, satAngle):
    currPosition = polarToCart(satellite[0], satellite[1])
    velVector = polarToCart(satVelocity, satAngle)

    currPosition[0] = currPosition[0] + velVector[0] + \
        (acceleration[0]*acceleration[0])/2
    currPosition[1] = currPosition[1] + velVector[1] + \
        (acceleration[1]*acceleration[1])/2

    newPosition = cartToPolar(currPosition[0], currPosition[1])
    return newPosition


def velocityChange(acceleration, satVelocity, satAngle):
    velVector = polarToCart(satVelocity, satAngle)
    velVector[0] += acceleration[0]
    velVector[1] += acceleration[1]
    return velVector


while maxiter != 0:
    rotatePlanets()
    pull = gravPull(satellite[0], satellite[1])
    print(positionChange(pull, startParameters[0], startParameters[1]))
    print(velocityChange(pull, startParameters[0], startParameters[1]))
    maxiter -= 1
