
def forward_difference(fun, t, dt):
    from lielab.functions import log
    
    f0 = fun(t)
    f1 = fun(t+dt)
    out = log(f0.inverse()*f1)/dt
    return out

def backward_difference(fun, t, dt):
    from lielab.functions import log
    
    f0 = fun(t-dt)
    f1 = fun(t)
    out = log(f0.inverse()*f1)/dt
    return out

def central_difference(fun, t, dt):
    from lielab.functions import log
    
    f0 = fun(t-dt)
    f1 = fun(t+dt)
    out = log(f0.inverse()*f1)/(2.0*dt)
    return out
