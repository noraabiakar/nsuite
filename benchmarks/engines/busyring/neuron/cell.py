import random
from neuron import h

class cell_parameters:
    def __repr__(self):
        return 'branchy cell parameters: depth {}; branch_probs {}; compartments {}; lengths {}'\
                .format(self.max_depth, self.branch_probs, self.compartments, self.lengths)

    def __init__(self, max_depth, branch_prob, compartment, length, synapses):
        self.max_depth = max_depth          # maximum number of levels
        self.branch_probs = branch_prob      # range of branching probabilities at each level
        self.compartments = compartment      # range of compartment counts at each level
        self.lengths = length                # range of lengths of sections at each level
        self.synapses = synapses             # the nyumber of synapses per cell

def interp(r, i, n):
    p = i * 1.0/(n-1)
    return (1-p)*r[0] + p*r[1]

def printcell(c):
    print('cell with ', len(c.sections), ' levels:')
    s = ''
    for l in range(len(c.sections)):
        s += '  level {} : {} sections\n'.format(l, len(c.sections[l]))
    print(s)

#
#   Branching cell.
#   It branches, and stuff
#
class branchy_cell:
    def __repr__(self):
        s = 'cell_%d\n' % self.gid
        return s

    def __init__(self, gid, params):
        self.pc = h.ParallelContext()
        self.gid = gid

        # generate the soma
        soma = h.Section(name='soma', cell=self)
        soma.L = soma.diam = 12.6157 # Makes a soma of 500 microns squared
        soma.Ra = 100
        soma.cm = 1
        soma.insert('hh')
        for seg in soma:
            seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
            seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
            seg.hh.gl = 0.0003    # Leak conductance in S/cm2
            seg.hh.el = -54.3     # Reversal potential in mV

        self.sections = []
        self.sections.append([soma])

        self.nseg = 0
        self.ncomp = 0

        # build the dendritic tree
        nlev = params.max_depth
        random.seed(gid) # seed the random number generator on gid
        flat_section_list = [soma]
        for i in range(params.max_depth):
            level_secs = []
            count = 0
            # branch prob at this level
            bp = interp(params.branch_probs, i, params.max_depth)
            # length at this level
            l = interp(params.lengths, i, params.max_depth)
            # number of compartments at this level
            nc = round(interp(params.compartments, i, params.max_depth))

            j = 0
            for sec in self.sections[i]:
                # attempt to make some branches
                if random.uniform(0, 1) < bp:
                    for branch  in [0,1]:
                        dend = h.Section(name='dend{}_{}'.format(i, count))
                        dend.L = l      # microns
                        dend.diam = 1   # microns
                        dend.Ra = 100
                        dend.cm = 1
                        dend.nseg = nc
                        dend.insert('pas')
                        for seg in dend:
                            seg.pas.g = 0.001  # Passive conductance in S/cm2
                            seg.pas.e = -65    # Leak reversal potential mV
                        #for seg in dend:
                        #    seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
                        #    seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
                        #    seg.hh.gl = 0.0003    # Leak conductance in S/cm2
                        #    seg.hh.el = -54.3     # Reversal potential in mV

                        dend.connect(sec(1))
                        level_secs.append(dend)
                        flat_section_list.append(dend)
                        count += 1
                        self.ncomp += nc
                j += 1

            self.nseg += count
            if count==0:
                break

            self.sections.append(level_secs)

        self.soma = soma

        # stick a synapse on the soma
        self.synapses = [h.ExpSyn(self.soma(0.5))]
        self.synapses[0].tau = 2

        # add additional synapses that will be connected to the "ghost" network
        for i in range(1, params.synapses):
            seg = random.randint(1, self.nseg-1)
            pos = random.uniform(0, 1)
            self.synapses.append(h.ExpSyn(flat_section_list[seg](pos)))

        self.halfgap_list = []

    def soma_assign_v(self):
        if self.pc.gid_exists(self.gid):
            sec_seg = self.pc.gid2cell(self.gid).soma(0.5)
            self.pc.source_var(sec_seg._ref_v, self.gid, sec=sec_seg.sec)

    # Add a gap_junction between self and other
    # the voltages of 'self' and 'other' need to be visible to the other cell's half gap junction
    # 'source_var' assigns a voltage variable to a unique id
    # 'target_var' attaches a voltage variable (identified using its unique id) to another voltage variable
    # to expose the voltage of 'self' to the half gap_junction at 'other':
    # 1. assign the voltage of a sec on 'self' to a unique id (cell gid) using 'source_var'
    # 2. attach the voltage of the half gap_junction at 'other' to the voltage of a sec of 'self'
    #    using 'target_var' and the unique id (gid) of a sec of 'self'
    def add_point_gap(self, other, ggap, loc1 = 0.5, loc2 = 0.5, name_sec1=None, name_sec2=None):  #ggap in nS
        if self.pc.gid_exists(self.gid):
            self.mk_halfgap(other, ggap, name_sec1, loc1)

        if self.pc.gid_exists(other.gid):
            other.mk_halfgap(self, ggap, name_sec2, loc2)

    # assign the voltage at a sec to the gid of the cell
    # create half gap_junction at a sec, and assign its variables: vgap and g
    # vgap gets the voltage assigned to the gid of the 'other' cell
    # g gets ggap
    def mk_halfgap(self, other, ggap, name_sec, loc):
        # sec seg
        #if name_sec==None:
        sec_seg = self.pc.gid2cell(self.gid).soma(0.5)
        # else:
        #     sec_seg = self.pc.gid2cell(self.gid).sections[name_sec](loc)

        # assign the voltage at the soma to the gid of the cell
        # self.pc.source_var(sec_seg._ref_v, self.gid, sec=sec_seg.sec)

        # create half gap_junction on the soma
        hg = h.HalfGap(sec_seg)

        # attach vgap to the voltage assigned for the 'other' cell's gid
        self.pc.target_var(hg, hg._ref_vgap, other.gid)

        # set the conductance of the half gap_junction
        # must match the second half of the gap_junction
        hg.g = ggap

        # save the state
        self.halfgap_list.append(hg)


    def set_recorder(self):
        """Set soma, dendrite, and time recording vectors on the cell.

        :param cell: Cell to record from.
        :return: the soma, dendrite, and time vectors as a tuple.
        """
        soma_v = h.Vector()   # Membrane potential vector at soma
        dend_v = h.Vector()   # Membrane potential vector at dendrite
        t = h.Vector()        # Time stamp vector
        soma_v.record(self.soma(0.5)._ref_v)
        dend_v.record(self.dend(0.5)._ref_v)
        t.record(h._ref_t)
        return soma_v, dend_v, t
