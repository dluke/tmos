
## 2020-07-30 Pele dependent code moved out of ctfp3d.py
# this file need imports to be executed


def debug_event_function():
    f = open('pele_min_state.out', 'w') 
    def event_function(coords, energy, rms):
        x, y, z, px, py, pz = coords
        fstring = '{} '*8 + '\n'
        f.write(fstring.format(x,y,z,px,py,pz, energy, rms))
    return f, event_function

def spgrad_event_function(bp):
    fp = open('pgrad.out', 'w') 
    fs = open('sgrad.out', 'w') 
    fstring = '{} '*7 + '\n'
    def p_event_function(coords, energy, rms):
        pgrad = bp.cell.pili_grad()
        penergy = bp.cell.energy_pili()
        px, py, pz, ax, ay, az = pgrad
        fp.write(fstring.format(px, py, pz, ax, ay, az, penergy))
    def s_event_function(coords, energy, rms):
        sgrad =  bp.cell.surface_grad()
        senergy = bp.cell.energy_surface()
        px, py, pz, ax, ay, az = sgrad
        fs.write(fstring.format(px, py, pz, ax, ay ,az, senergy))
    return (fp, p_event_function), (fs, s_event_function)


def break_pelefire(maxsteps, bp, FireOptimizer, events):

    for step in range(maxsteps):
        if FireOptimizer.stop_criterion_satisfied():
            return FireOptimizer.get_result()
        FireOptimizer.one_iteration()
        # retrieve axisangle state vector from optimiser an explicitly assign it to the body object
        x = FireOptimizer.get_x()
        origin, p = x[:3], x[3:]
        bp.cell.body.set_pstate( seqtov3d(p) ) # outdated
        for event in events:
            result = FireOptimizer.get_result()
            event(result.coords, result.energy, result.rms)
    return FireOptimizer.get_result()



#import sys; sys.path.insert(0, "/home/dan/.local/lib/python2.7/site-packages")
import pele.optimize
import pele.potentials
pele.optimize import ModifiedFireCPP
def setup_pelefire(args):

    mindir= './minloop/'

    def trypele(cell):
        class pelePotential(pele.potentials.BasePotential):
            def __init__(self, cell):
                self.cell = cell

            def getEnergy(self, state):
                self.cell.set_state(state)
                return self.cell.energy()

            def getGradient(self, state):
                self.cell.set_state(state)
                return self.cell.grad()

            def getEnergyGradient(self, state):
                self.cell.set_state(state)
                #
                return self.cell.energy(), self.cell.grad()

        def make_vtkcellout(cell):
            outform = 'minstep_cell_{:08d}.vtp'
            poutform = 'minstep_pili_{:08d}.vtp'
            interval = tfputils.make_interval(step=10)
            counter = tfputils.make_counter()
            def vtkcellout(coords, energy, rms):
                step = next(counter)
                if next(interval):
                    outfile = outform.format(step)
                    poutfile = poutform.format(step)
                    out = os.path.join(mindir, outfile)
                    piliout = os.path.join(mindir, poutfile)
                    vtkwriter.write_cell3d(cell, out)
                    vtkwriter.write_pili3d(cell, piliout)
            return vtkcellout

        # saved cell just for printing as output
        #savecell = cell.clone()
        originalstate = cell.get_state()

        # max step = 0.001
        #maxstep = args.Pili.d_bound/40.
        maxstep = 1.

        vtkevent = make_vtkcellout(cell)
        bp = pelePotential(cell)

        def pyminimise(bp, dtmax, steps):
            state = bp.cell.get_state()
            Fire = ModifiedFireCPP(state, bp, 
                    dtstart=args.fire.dt,
                    dtmax=dtmax,
                    maxstep=maxstep,
                    tol=args.fire.force_target/np.sqrt(6),
                    stepback=False)
            result = Fire.run(steps)
            return result

        def minimise(bp, dtmax, steps):
            state = bp.cell.get_state()
            Fire = pele.optimize.Fire(state, bp, 
                    dt=args.fire.dt,
                    dtmax=dtmax,
                    maxstep=maxstep)
            tol = args.fire.force_target/np.sqrt(6)
            result = Fire.run(tol, steps)
            return result

        #result = pyminimise(bp, args.fire.dt_max, args.fire.maxsteps)
        #result = minimise(bp, args.fire.dt_max, args.fire.maxsteps)
        def break_cppmin(bp, dtmax, steps, events=[]):
            bp.cell.body.initpstate() # outdated
            Fire = ModifiedFireCPP(bp.cell.get_state(), bp, 
                    dtstart=args.fire.dt,
                    dtmax=dtmax,
                    maxstep=maxstep,
                    tol=args.fire.force_target/np.sqrt(6),
                    stepback=False)
            result = break_pelefire(steps, bp, Fire, events)
            return result
        result = break_cppmin(bp, args.fire.dt_max, args.fire.maxsteps)

        def lower_step_try():
            bp.cell.set_state(originalstate)
            newmaxsteps = 10*args.fire.maxsteps # not set
            new_dt_max = args.fire.dt_max/10.
            print("new maxsteps", newmaxsteps)
            print("new max dt", new_dt_max)
            #result = minimise(bp, dtmax=new_dt_max, steps=newmaxsteps)
            result = break_cppmin(bp, new_dt_max, newmaxsteps)
            return result

        if True: # on
            if result.success==False:
                print('Minimisation failed. Attempt minimisation with lower timestep.')
                print(result)
                result = lower_step_try()
                if result.success:
                    print("lower timestep success. Proceeding")
                else:
                    print("Still failed")

        if result.success==False:
            f, debug_event = debug_event_function()
            (fp, pevent), (fs, sevent) = spgrad_event_function(bp)
            bp.cell.set_state(originalstate)
            thisevents = [debug_event, pevent, sevent]
            if True: # switch for vtk output
                print('Making directory ', mindir)
                if not os.path.exists(mindir):
                    os.mkdir(mindir)
                thisevents.append(vtkevent)
            result = break_cppmin(bp, args.fire.dt_max, args.fire.maxsteps, events=thisevents)
            #result = minimise(bp, args.fire.dt_max, args.fire.maxsteps, events=thisevents)
            f.close()
            fp.close()
            fs.close()

        if args.system.debug:
            print(result)
        if result.success==False:
            print("Minimisation failed")
            print(result)
            print("exiting simulation with seed ", args.system.seed)
            print("pilivar was ", args.ACell.pilivar)
            sys.exit()
        if result.success==True:
            # set the final state
            cell.set_state(result.coords) # check can remove

        return result.nsteps, result.rms
    return trypele


