# This is the code for the original model javascript code found at: https:#facultyweb.cs.wwu.edu/~jagodzf/covid-19/original/
import argparse
import datetime
import math
import plotly.express as px
import sys
from functools import reduce
#RT is the same as R0 (R naught)
#RT is the rate of transmission
#If R0 is 2, then that person will give it to two other people
#R0 ranges between 0.5 and 4        
# If R0 < 1.2 there's not enough transmission that it doesn't spread

def strType(xstr):
        try:
            int(xstr)
            return 'int'
        except:
            try:
                float(xstr)
                return 'float'
            except:
                try:
                    complex(xstr)
                    return 'complex'
                except:
                    return 'str'
                    
def toDate(startDate, days):
        # return new date beginning on Jan 1st 2020 + t
        # unsure what 't' is
        return startDate + datetime.timedelta(days=days)

def ReadDecayLine(line):
    lineValues = line.split(" ")
    if (len(lineValues) != 2):
        print(line + ": does not contain exactly 2 values.")
        return None
    
    dayValues = lineValues[0].split("=")
    if (len(dayValues) != 2 or dayValues[0] != "day" or not dayValues[1].isdigit()):
        print(lineValues[0] + ": should be formatted as day={int}")
        return None
    
    metricValues = lineValues[1].split("=")
    if (len(metricValues) != 2 or metricValues[0] != "R0" or strType(metricValues[1]) != "float"):
        print(lineValues[1] + ": should be formatted as R0={float}")
        return None

    #return (day, metric, metricValue)
    return (int(dayValues[1]), "R0", float(metricValues[1]))

def GetR0DecayValues(R0FilePath):
    arrayR0 = []
    XAxis = 0
    YAxis = 2
    
    R0File = open(R0FilePath, "r")
    line1 = R0File.readline().strip()
    if (line1 == ""):
        print("The first line must contain information.")
        return None
    
    line1Values = ReadDecayLine(line1)
    if (line1Values == None):
        return None
    
    if (line1Values[XAxis] != 0):
        print("The first line within the given decay file must have day=0")
        return None

    for line2 in R0File:
        line2 = line2.strip()
        # skip empty lines
        if (line2 == ""):
            continue

        line2Values = ReadDecayLine(line2)
        if (line2Values == None):
            return None
        
        if (line2Values[XAxis] < line1Values[XAxis]):
            print("Cannot go down in days in the given decay file.")
            return None
        
        deltaX = float(line2Values[XAxis] - line1Values[XAxis]) #day values
        deltaY = float(line2Values[YAxis] - line1Values[YAxis]) #R0 Values
        slope = float(deltaY / deltaX)
        b = float(line1Values[YAxis] - slope * line1Values[XAxis])

        arrayR0.append(line1Values[YAxis])
        for j in range(line1Values[XAxis]+1, line2Values[XAxis]):
            arrayR0.append(float(slope * j + b))
        
        line1Values = line2Values
    
    if (line2Values[XAxis] < 365):
        print("The last line within the given decay file must have a day that is greater than or equal to 365.")
        return None

    #currently is not guaranteed to place values in all 365 indexes of the array
    return arrayR0

def BetterCommandLineArgReader():
    parser = argparse.ArgumentParser(description='SEIR Model For Covid 19')
    parser.add_argument("-population", action="store", type=int, dest="N",
                        help="The population size for the model")
    parser.add_argument("-r0", action="store", type=float,
                        help="The rate of transmission to be used for the ENTIRE population")
    parser.add_argument("-i0", action="store", type=int)
    parser.add_argument("-cfr", action="store", type=float)
    parser.add_argument("-psevere", action="store", type=float,
                        help="The probability that the infection is severe")
    parser.add_argument("-hl", "--hospital_lag", action="store", type=int, 
                        help="The hospital lag")
    parser.add_argument("-decay", action="store", 
                        help="Indicate path to .txt file containing age groups and R0 values")
    return parser.parse_args()

def UpdateDefaultValues(defaultValues, arguments):
    argumentsDict = vars(arguments)
    for arg in argumentsDict:
        if argumentsDict[arg] is not None:
            defaultValues[arg.upper()] = argumentsDict[arg]

# f is a func of time t and state y
# y is the initial state, t is the time, h is the timestep
# updated y is returned.
def integrate(m, fn, y, t, h):
        k = []
        for ki in range(0, len(m)):
            _y = y.copy()
            if (ki):
                dt = m[ki-1][0] * h
            else:
                dt = 0
            
            for l in range(0, len(_y)):
                for j in range(0, ki):
                    _y[l] = _y[l] + h * m[ki - 1][j] * k[ki - 1][l]
            #k[ki] = f(t + dt, _y, dt)
            k.append(fn(t + dt, _y))
        
        r = y.copy()
        for l in range(0, len(_y)):
            for j in range(0, len(k)):
                r[l] = r[l] + h * k[j][l] * m[ki][j]
        return r

Integrators = {
    "Euler": [[1]],
    "Midpoint": [
        [0.5, 0.5],
        [0, 1],
    ],
    "Heun": [
        [1, 1],
        [0.5, 0.5],
    ],
    "Ralston": [
        [2 / 3, 2 / 3],
        [0.25, 0.75],
    ],
    "K3": [
        [0.5, 0.5],
        [1, -1, 2],
        [1 / 6, 2 / 3, 1 / 6],
    ],
    "SSP33": [
        [1, 1],
        [0.5, 0.25, 0.25],
        [1 / 6, 1 / 6, 2 / 3],
    ],
    "SSP43": [
        [0.5, 0.5],
        [1, 0.5, 0.5],
        [0.5, 1 / 6, 1 / 6, 1 / 6],
        [1 / 6, 1 / 6, 1 / 6, 1 / 2],
    ],
    "RK4": [
        [0.5, 0.5],
        [0.5, 0, 0.5],
        [1, 0, 0, 1],
        [1 / 6, 1 / 3, 1 / 3, 1 / 6],
    ],
    "RK38": [
        [1 / 3, 1 / 3],
        [2 / 3, -1 / 3, 1],
        [1, 1, -1, 1],
        [1 / 8, 3 / 8, 3 / 8, 1 / 8],
    ],
}

def f(defaultValues):
    #default values for the model
    Time_to_death = 25
    #logN = math.log(7e6) not used
    #N = 7800000
    #I0 = 1
    #R0 = 2.5
    D_incubation = 5.0
    D_infectious = 3.0
    D_recovery_mild = 11.0
    D_recovery_severe = 21.0
    #DHospitalLag = 8
    D_death = Time_to_death - D_infectious
    #CFR = 0.01
    InterventionTime = 10000
    InterventionAmt = 1 / 3
    #Time = 220 not used
    #Xmax = 110000 not used
    #P_SEVERE = 0.04
    duration = 7 * 12 * 1e10
    #seasonal_effect = 0 #not used

    interpolation_steps = 40
    steps = 320 * interpolation_steps
    dt = 1 / interpolation_steps
    sample_step = interpolation_steps

    N = defaultValues["N"]
    I0 = defaultValues["I0"]
    R0 = defaultValues["R0"]
    CFR = defaultValues["CFR"]
    P_SEVERE = defaultValues["PSEVERE"]
    DHospitalLag = defaultValues["HOSPITALLAG"]
    UseDecayingR0 = defaultValues["UseDecayingR0"]
  
    StartDate = datetime.datetime(2020, 1, 15)
    method = Integrators["RK4"]
    def f(t, x):
        nonlocal R0
        nonlocal UseDecayingR0
        if (UseDecayingR0):
            CurrentDate = toDate(StartDate, t)
            DifferenceInDays = CurrentDate - StartDate
            R0 = defaultValues["arrayOfR0s"][DifferenceInDays.days]

        # SEIR ODE
        if (t > InterventionTime and t < InterventionTime + duration):
            beta = (InterventionAmt * R0) / D_infectious
        elif (t > InterventionTime + duration):
            beta = (0.5 * R0) / D_infectious
        else:
            beta = R0 / D_infectious

        alpha = 1 / D_incubation
        gamma = 1 / D_infectious

        S = x[0] # Susectable
        E = x[1] # Exposed
        I = x[2] # Infectious
        Mild = x[3] # Recovering (Mild)
        Severe = x[4] # Recovering (Severe at home)
        Severe_H = x[5] # Recovering (Severe in hospital)
        Fatal = x[6] # Recovering (Fatal)
        R_Mild = x[7] # Recovered
        R_Severe = x[8] # Recovered
        R_Fatal = x[9] # Dead

        p_severe = P_SEVERE
        p_fatal = CFR
        p_mild = 1 - P_SEVERE - CFR

        dS = -beta * I * S
        dE = beta * I * S - alpha * E
        dI = alpha * E - gamma * I
        dMild = p_mild * gamma * I - (1 / D_recovery_mild) * Mild
        dSevere = p_severe * gamma * I - (1 / DHospitalLag) * Severe
        dSevere_H = (1 / DHospitalLag) * Severe - (1 / D_recovery_severe) * Severe_H
        dFatal = p_fatal * gamma * I - (1 / D_death) * Fatal
        dR_Mild = (1 / D_recovery_mild) * Mild
        dR_Severe = (1 / D_recovery_severe) * Severe_H
        dR_Fatal = (1 / D_death) * Fatal

        #      0   1   2   3      4        5          6       7        8          9
        return [
            dS,
            dE,
            dI,         #2
            dMild,
            dSevere,
            dSevere_H,  #5
            dFatal,
            dR_Mild,
            dR_Severe,  #8
            dR_Fatal,
        ]

    v = [1 - I0 / N, 0, I0 / N, 0, 0, 0, 0, 0, 0, 0]
    t = 0

    P = []
    while (steps):
        if ((steps + 1) % sample_step == 0):
            P.append({
            "Time": toDate(StartDate, t),
            "R0": R0,
            "HospitalLag": DHospitalLag,
            "Dead": N * v[9],
            "Susceptible": N * v[0],
            "Hospital": N * (v[5] + v[6]),
            "RecoveredMild": N * v[7],
            "RecoveredSevere": N * v[8],
            "RecoveredTotal": N * (v[7] + v[8]),
            "Infected": N * v[2],
            "Exposed": N * v[1],
            "Sum": N * reduce((lambda a,b: a+b), v)
            })
        v = integrate(method, f, v, t, dt)
        t += dt #this may not be getting updated properly since integrate() is no longer nested within f()
        steps -= 1

    return P

def getTrace(data, name, metric):
    y = []
    x = []
    for i in range (0, len(data)):
        y.append(data[i][metric])
        x.append(data[i]["Time"])

    trace = {
    "x": x,
    'y': y,
    "type": "scatter",
    "name": name,
    }
    return trace

def WeightedAverageR0(values, percentages):
    avg = 0
    for key in values.keys():
        avg += float(values[key] * percentages[key])
    return avg

def main():
    # default values for variables that can be modified with command line arguments go here
    defaultValues = {
        "N": 226387,
        "I0": 1,
        "R0": 2.85,
        "CFR": 0.01,
        "PSEVERE": 0.04,
        "HOSPITALLAG": 8,
        "UseDecayingR0": False,
        "arrayOfR0s": None,
    }

    ageGroupsR0Values = {
        "0-19": 0.8,
        "19-40": 2.0,
        "40-60": 1.3,
        "60-80": 1.1,
        "80+": 1.1
    }

    ageGroupPopulationPercentages = {
        "0-19": 0.24,
        "19-40": 0.3,
        "40-60": 0.22,
        "60-80": 0.2,
        "80+": 0.04
    }

    args = BetterCommandLineArgReader()
    UpdateDefaultValues(defaultValues, args)

    #Will be true if the -decay flag is present with a file path
    defaultValues["UseDecayingR0"] = args.decay is not None
    if (defaultValues["UseDecayingR0"]):
        print("R0 value will decay")
        defaultValues["arrayOfR0s"] = GetR0DecayValues(args.decay)

    print(defaultValues)
    data = f(defaultValues)
    infectedPlotData = getTrace(data, "Infected, seasonal effect = 0", "Infected")
    infectedPlot = px.line(x=infectedPlotData["x"], y=infectedPlotData["y"], title=infectedPlotData["name"])
    infectedPlot.show()


if __name__ == "__main__":
    main()
