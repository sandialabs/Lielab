from lielab.cppLielab.integrate import ODESolution as _ODESolution

class ODESolution(_ODESolution):
    def __repr__(self):
        # order = ['messages', 'status', 't', 'y', 'ybar', 'p']
        strout = ''
        strout += '{:>14}: {}\n'.format('message', self.message)
        strout += '{:>14}: {}\n'.format('status', self.status)
        strout += '{:>14}: {}\n'.format('t', 'NumPy Array ' + str(self.t.shape))
        strout += '{:>14}: {}\n'.format('y', 'CompositeManifold (' + str(len(self.y)) + ',)')
        strout += '{:>14}: {}\n'.format('ybar', 'NumPy Array ' + str(self.ybar.shape))
        strout += '{:>14}: {}\n'.format('theta', 'CompositeAlgebra (' + str(len(self.theta)) + ',)')
        strout += '{:>14}: {}\n'.format('thetabar', 'NumPy Array ' + str(self.thetabar.shape))
        return strout
