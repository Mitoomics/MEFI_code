def cut(obj, sec):
    return [obj[i:i+sec] for i in range(0,len(obj),sec)]
