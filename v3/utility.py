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

def WeightedAverage(values, percents):
    output = 0
    for key in values.keys():
        output += float(values[key] * percents[key])
    return output