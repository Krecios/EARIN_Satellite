import math
maxiter = 5  # maximum possible number of iterations
# temporary attributes of the satellite for testing
gravConst = 3
satVelocity = 10
satAngle = 15
startTime = 3
satelliteMass = 0.00001
startPanetID = 2
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


def polarToCart(r, fi):
    x = r * math.cos(fi)
    y = r * math.sin(fi)
    coords = [x, y]
    return coords


def cartToPolar(x, y):
    r = math.sqrt(x*x + y*y)
    fi = math.atan(y/x)
    coords = [r, fi]
    return coords


def gravlPull(v, alpha):
    print("WIP")


while maxiter != 0:
    rotatePlanets()
    maxiter -= 1
