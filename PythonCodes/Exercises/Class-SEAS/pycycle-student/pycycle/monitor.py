import numpy as np
from scipy.integrate import quad

def print_time(t):
    ms    =  round(1000*t)
    y, ms = divmod(ms, 1000*60*60*24*365)
    d, ms = divmod(ms, 1000*60*60*24)
    h, ms = divmod(ms, 1000*60*60)
    m, ms = divmod(ms, 1000*60)
    s, ms = divmod(ms, 1000)
    print(
        "{0:>6}".format(y), 'yr,',
        "{0:>3}".format(d), 'd,',
        "{0:>2}".format(h), 'h,',
        "{0:>2}".format(m), 'm,',
        "{0:>2}".format(s), 's,',
        "{0:>3}".format(ms),'ms.')
    return

class Monitor:
    def __init__(self, thresholds, u_ax, u_fig, v_ax, v_fig):
        self.t          = []
        self.v          = []
        self.sol_stack  = []
        self.plt_stack  = []
        self.thresholds = thresholds
        self.u_ax       = u_ax
        self.u_fig      = u_fig
        self.v_ax       = v_ax
        self.v_fig      = v_fig
        return

    def __call__(self, t, u, v, psi, tau):
        print_time(t)
        print('       v_max = ', v.max(), 'm/s | ', v.max()*100*(60*60*24*365), 'cm/yr')
        
        self.t.append(t)
        self.v.append(v.max())
        
        self.v_ax.cla()
        self.v_ax.plot(self.t, np.log10(self.v), color='#ffcc00',marker='.')
        self.v_fig.canvas.draw()
        
        self.sol_stack.append({'t': t, 'u': u, 'v': v, 'psi': psi, 'tau': tau})
        
        creg = self.thresholds[0]
        for thr in self.thresholds[1:]:
            if v.max() > thr['vthresh']:
                creg = thr
        
        if len(self.plt_stack) == 0:
            self.plt_stack.append({'c': creg['color'], 't': t, 'u': u, 'v': v, 'psi': psi, 'tau': tau})
            self.u_ax.plot(self.plt_stack[-1]['u'],color=creg['color'], linewidth=0.5)
            self.u_fig.canvas.draw()
        
        elif len(self.sol_stack) >= 2:
            t0 = self.plt_stack[-1]['t']
            t1 = self.sol_stack[-2]['t']
            t2 = self.sol_stack[-1]['t']
            tx = t0 + creg['dt']
            while tx < t2:
                if t1 < tx:
                    print('appending to plot stack')
                    
                    self.plt_stack.append({'c': creg['color'], 't': tx, 'u': [], 'v': [], 'psi': [], 'tau': []})
                    
                    weights = [1 - (tx - t1) / (t2 - t1), (tx - t1) / (t2 - t1)]
                    self.plt_stack[-1]['u']   = np.average([self.sol_stack[-2]['u']  , self.sol_stack[-1]['u']  ], axis=0, weights=weights)
                    self.plt_stack[-1]['v']   = np.average([self.sol_stack[-2]['v']  , self.sol_stack[-1]['v']  ], axis=0, weights=weights)
                    self.plt_stack[-1]['psi'] = np.average([self.sol_stack[-2]['psi'], self.sol_stack[-1]['psi']], axis=0, weights=weights)
                    self.plt_stack[-1]['tau'] = np.average([self.sol_stack[-2]['tau'], self.sol_stack[-1]['tau']], axis=0, weights=weights)
                    self.u_ax.plot(self.plt_stack[-1]['u'], color=creg['color'], linewidth=0.5)
                    self.u_fig.canvas.draw()
                
                tx += creg['dt']
                
        if len(self.sol_stack) >= 10:
            self.sol_stack.pop(0)
