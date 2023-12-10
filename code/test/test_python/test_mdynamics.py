
import tmos
import parameters
import ctfp3d as tfp
import wrinterval as wr

def TestSummary():
    args = parameters.read().apply()
    cell = tfp.setup_cell(args)
    md = tfp.setup_md3d(args)
    md.step(cell)
    summary = md.get_summary()
    print(summary.rms)
    print(summary.nsteps)
    print(summary.prior_state)
    print(summary.accepted_state)
    # 
    tmp_file = 'mdev.tmp.dat'
    mdev = wr.Track3d.setup_mdevent(tmp_file)
    condition = [True, False, False]
    # TODO update this line
    mdev.track_mdevent(0., 'None', condition, summary)

    with open(tmp_file, 'r') as f:
        header = f.readline()
        print(header)
        line = f.readline()
        print(line)



if __name__=='__main__':

    TestSummary()





