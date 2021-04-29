class ODEsolver:
    def __init__(self, fun, t0, y0):
        self.fun = fun
        self.t0 = t0
        self.y0 = y0

    def eulersolver(self, tend, step=None):
        if step is None:
            step = (tend - self.t0) / 1e3
        from ArrayTable import ArrayTable
        result = ArrayTable(2, 0)
        result.setTableHeader(["$t$", "y(t) solved with explicit euler"])
        result.setTableUnit(["/", "/"])
        result.append([self.t0, self.y0])
        t = self.t0
        while t < tend:
            ynext = result.table[1].data[-1] + self.fun(t, result.table[1].data[-1]) * step
            result.append([t + step, ynext])
            t += step
        return result

    def forecastCorrection(self, tend, step=None):
        if step is None:
            step = (tend - self.t0) / 1e3
        from ArrayTable import ArrayTable
        result = ArrayTable(2, 0)
        result.setTableHeader(["$t$", "y(t) solved with forecast correction"])
        result.setTableUnit(["/", "/"])
        result.append([self.t0, self.y0])
        t = self.t0
        while t < tend:
            dydt=self.fun(t, result.table[1].data[-1])
            ynextpridict=result.table[1].data[-1] + dydt * step
            ynext=result.table[1].data[-1]+step/2.*(dydt+self.fun(t+step,ynextpridict))
            result.append([t + step, ynext])
            t += step
        return result

    def Rungekutta4(self, tend, step=None):
        if step is None:
            step = (tend - self.t0) / 1e3
        from ArrayTable import ArrayTable
        result = ArrayTable(2, 0)
        result.setTableHeader(["$t$", "y(t) solved with runge kutta 4"])
        result.setTableUnit(["/", "/"])
        result.append([self.t0, self.y0])
        t = self.t0
        while t < tend:
            s1=self.fun(t,result.table[1].data[-1])
            s2=self.fun(t+step/2.,result.table[1].data[-1]+step/2.*s1)
            s3=self.fun(t+step/2.,result.table[1].data[-1]+step/2.*s2)
            s4=self.fun(t+step,result.table[1].data[-1]+step*s3)
            ynext=result.table[1].data[-1]+step/6.*(s1+2*s2+2*s3+s4)
            result.append([t + step, ynext])
            t += step
        return result

    def implicitEuler(self,tend, step=None):
        if step is None:
            step = (tend - self.t0) / 1e3
        from ArrayTable import ArrayTable
        result = ArrayTable(2, 0)
        result.setTableHeader(["$t$", "y(t) solved with implicit euler"])
        result.setTableUnit(["/", "/"])
        result.append([self.t0, self.y0])
        t = self.t0
        while t<tend:
            temp=self.fun(t,result.table[1].data[-1])
            ynextinit=result.table[1].data[-1]+step/2.*(temp+self.fun(t+step,result.table[1].data[-1] + temp * step))
            def func(w):
                return result.table[1].data[-1]+step*self.fun(t+step,w)-w

            ynext=ynextinit;h=step/10.
            while func(ynext)>1.e-5:
                ynext-=func(ynext)/((func(ynext+h)-func(ynext-h))/2./h)
            result.append([t + step, ynext])
            t += step
        return result


