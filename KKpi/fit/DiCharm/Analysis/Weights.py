#!usr/bin/env python
# ============================================================================
from AnalysisPython.PyRoUts import VE, SE
import DiCharm.Ana as Ana
# ============================================================================


def eff_fix(eff, sigma):
    """
    """
    #
    val = eff. value()
    err = eff. error()
    err = min(err, 1.0)
    #
    nval = val + err * sigma
    #
    nval = max(0.1 * val, nval)
    nval = min(1.9 * val, nval)
    #
    # nval = min (       1.0 , nval )
    #
    return VE(nval, err * err)


# ============================================================================
# the generic cuts
class Cuts  (object):

    def __init__(self,
                 first,
                 second,
                 tos_1,
                 tos_2,
                 ptmin):

        first_name = '_' + first
        second_name = '_' + second
        first_trg = first + '_'
        second_trg = second + '_'

        self.tos1 = lambda s: Ana.isTos(s,  first_trg)
        self.tos2 = lambda s: Ana.isTos(s, second_trg)

        self._good_1 = lambda s: getattr(s, 'good' + first_name)
        self._good_2 = lambda s: getattr(s, 'good' + second_name)
        self._prompt_1 = lambda s: getattr(s, 'prompt' + first_name)
        self._prompt_2 = lambda s: getattr(s, 'prompt' + second_name)

        if tos_1 and tos_2:
            self.tos = lambda s: self.tos1(s) or self.tos2(s)
            ## self.tos = lambda s : self.tos1 ( s ) and self.tos2 ( s )
        elif tos_1:
            self.tos = lambda s: self.tos1(s)
        elif tos_2:
            self.tos = lambda s: self.tos2(s)
        else:
            self.tos = lambda s: True

        self.pt_1 = lambda s: getattr(s, 'pt' + first_name)
        self.pt_2 = lambda s: getattr(s, 'pt' + second_name)
        self.y_1 = lambda s: getattr(s, 'y' + first_name)
        self.y_2 = lambda s: getattr(s, 'y' + second_name)

        pt_min = ptmin
        self.pt_min = pt_min

        print ' CUTS: ptmin = %f GeV ' % pt_min

        if 0 <= first_name.find('psi'):
            self.ok1 = lambda s:      0 <= self.pt_1(
                s) <= 12 and 2 <= self.y_1(s) <= 4.0
        else:
            self.ok1 = lambda s: pt_min <= self.pt_1(
                s) <= 12 and 2 <= self.y_1(s) <= 4.0

        self.ok2 = lambda s: pt_min <= self.pt_2(
            s) <= 12 and 2 <= self.y_2(s) <= 4.0

        self.kine = lambda s: self.ok1(s) and self.ok2(s)

        self.counter_acc_0 = SE()
        self.counter_acc_1 = SE()
        self.counter_acc_2 = SE()
        self.counter_acc_3 = SE()

    #
    def stats(self):
        return (self.counter_acc_0,
                self.counter_acc_1,
                self.counter_acc_2,
                self.counter_acc_3)

    def __call__(self, item):

        ok = True

        ok = ok and self.kine(item)
        #
        ok = ok and self._good_1(item)
        ok = ok and self._good_2(item)
        #
        ok = ok and self._prompt_1(item)
        ok = ok and self._prompt_2(item)
        #

        self.counter_acc_0 += ok

        # prompt
        ok = ok and item.prompt_2c

        self.counter_acc_1 += ok

        # trigger
        ok = ok and self.tos(item)

        self.counter_acc_2 += ok

        ok = ok and item.tmin > 0.001

        self.counter_acc_3 += ok

        return ok

# ============================================================================
# the generic weight for


class Weight (object):

    """
    Weight for 2xCharm
    """

    def __init__(self,
                 first,
                 second,
                 tos_1,
                 tos_2):

        first_name = '_' + first
        second_name = '_' + second
        first_trg = first + '_'
        second_trg = second + '_'

        if 0 <= first_name.find('psi'):
            self._eff_fun_1 = lambda s: Ana.eff_Jpsi(s,  first_name)
            self._eff_trg_1 = lambda s: Ana.trgEff_Jpsi(s,  first_name)
        elif 0 <= first_name.find('D0'):
            self._eff_fun_1 = lambda s: Ana.eff_D0(s,  first_name)
            self._eff_trg_1 = lambda s: Ana.trgEff_D0(s,  first_name)
        elif 0 <= first_name.find('Dp'):
            self._eff_fun_1 = lambda s: Ana.eff_Dp(s,  first_name)
            self._eff_trg_1 = lambda s: Ana.trgEff_Dp(s,  first_name)
        elif 0 <= first_name.find('Ds'):
            self._eff_fun_1 = lambda s: Ana.eff_Ds(s,  first_name)
            self._eff_trg_1 = lambda s: Ana.trgEff_Ds(s,  first_name)
        elif 0 <= first_name.find('Lc'):
            self._eff_fun_1 = lambda s: Ana.eff_Lc(s,  first_name)
            self._eff_trg_1 = lambda s: Ana.trgEff_Lc(s,  first_name)
        else:
            raise AttributeError("Invalid first  name '%s'" % first_name)

        if 0 <= second_name.find('D0'):
            self._eff_fun_2 = lambda s: Ana.eff_D0(s, second_name)
            self._eff_trg_2 = lambda s: Ana.trgEff_D0(s, second_name)
        elif 0 <= second_name.find('Dp'):
            self._eff_fun_2 = lambda s: Ana.eff_Dp(s, second_name)
            self._eff_trg_2 = lambda s: Ana.trgEff_Dp(s, second_name)
        elif 0 <= second_name.find('Ds'):
            self._eff_fun_2 = lambda s: Ana.eff_Ds(s, second_name)
            self._eff_trg_2 = lambda s: Ana.trgEff_Ds(s, second_name)
        elif 0 <= second_name.find('Lc'):
            self._eff_fun_2 = lambda s: Ana.eff_Lc(s, second_name)
            self._eff_trg_2 = lambda s: Ana.trgEff_Lc(s, second_name)
        else:
            raise AttributeError("Invalid second name '%s'" % second_name)

        #
        # note the inverted logic!
        #
        if tos_1 and tos_2:
            self._eff_trg = lambda s : 1 - \
                (1 - self._eff_trg_1(s)) * (1 - self._eff_trg_2(s))
            ## self._eff_trg = lambda s : self._eff_trg_1 ( s ) * self._eff_trg_2 ( s )
        elif tos_1:
            self._eff_trg = lambda s:           self._eff_trg_1(s)
        elif tos_2:
            self._eff_trg = lambda s:           self._eff_trg_2(s)
        else:
            self._eff_trg = lambda s:           VE(1, 0)

        self.counter_eff = SE()
        self.counter_w = SE()

        self.effs = {}
        self.effs['rec_1'] = SE()
        self.effs['rec_2'] = SE()
        self.effs['rec_Tr'] = SE()
        self.effs['pid_Pi'] = SE()
        self.effs['pid_K'] = SE()
        self.effs['pid_P'] = SE()
        self.effs['trigger'] = SE()

        self.pt_1 = lambda s: getattr(s, 'pt' + first_name)
        self.pt_2 = lambda s: getattr(s, 'pt' + second_name)
        self.y_1 = lambda s: getattr(s, 'y' + first_name)
        self.y_2 = lambda s: getattr(s, 'y' + second_name)

    #
    def stats(self):
        return self.counter_eff, self.counter_w, self.effs

    def __call__(self, item):

        ## acceptance & reconstruction & selection
        e1 = self._eff_fun_1(item)
        e2 = self._eff_fun_2(item)

        # PID: pions, kaons & protons
        ePi = Ana.pidEff_pions(item)
        eK = Ana.pidEff_kaons(item)
        eP = Ana.pidEff_protons(item)

        # correct for track recontruction efficiency
        eTr = Ana.recEff_tracks(item)

        # trigger efficiency
        eTrg = self._eff_trg(item)

        self.effs['rec_1'] += e1  . value()
        self.effs['rec_2'] += e2  . value()
        self.effs['pid_Pi'] += ePi . value()
        self.effs['pid_K'] += eK  . value()
        self.effs['pid_P'] += eP  . value()
        self.effs['rec_Tr'] += eTr . value()
        self.effs['trigger'] += eTrg . value()

        eff = VE(1, 0)

        # e1    = eff_fix ( e1 ,  1.0 )
        # e2    = eff_fix ( e2 ,  1.0 )

        eREC = e1 * e2

        ePID = ePi * eK * eP

        eTRG = eTrg

        eTRK = eTr

        # eREC = eff_fix ( eREC ,  1.0 )
        # ePID = eff_fix ( ePID , 1.0 )
        # eTRG = eff_fix ( eTRG ,  1.0 )
        # eTRK = eff_fix ( eTRK , -1.0 )

        eff = VE(1, 0)

        eff *= eREC
        eff *= ePID
        eff *= eTRG
        eff *= eTRK

        ## eff = eff_fix ( eff , +1.0 )
        ## eff = eff_fix ( eff , -1.0 )

        self.counter_eff += eff.value()

        # final result
        weight = 1.0

        if 0 < eff.value():
            weight = 1.0 / eff.value()
        else:
            weight = 0.0

        if 0 == weight or weight > 5.e+4:
            print ' ZERO weight : ', weight, \
                  (  self.pt_1( item ), self.y_1 ( item ) ), \
                  (  self.pt_2( item ), self.y_2 ( item ) ), \
                  (e1, e2, ePi, eK, eP, eTr, eTrg)

        self.counter_w += weight

        return weight

    pass

# ============================================================================
# Weights
# ============================================================================
#
# Psi&D
#


class WpsiD0(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'psi', 'D0', tos_1, tos_2)


class WpsiDp(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'psi', 'Dp', tos_1, tos_2)


class WpsiDs(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'psi', 'Ds', tos_1, tos_2)


class WpsiLc(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'psi', 'Lc', tos_1, tos_2)
#
# D&D
#


class WD0D0(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'D01', 'D02', tos_1, tos_2)


class WD0Dp(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'D0', 'Dp', tos_1, tos_2)


class WD0Ds(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'D0', 'Ds', tos_1, tos_2)


class WD0Lc(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'D0', 'Lc', tos_1, tos_2)


class WDpDp(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'Dp1', 'Dp2', tos_1, tos_2)


class WDpDs(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'Dp', 'Ds', tos_1, tos_2)


class WDpLc(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'Dp', 'Lc', tos_1, tos_2)


class WDsDs(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'Ds1', 'Ds2', tos_1, tos_2)


class WDsLc(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'Ds', 'Lc', tos_1, tos_2)


class WLcLc(Weight):

    def __init__(self, tos_1, tos_2):
        Weight.__init__(self, 'Lc1', 'Lc2', tos_1, tos_2)

# ============================================================================
# Cuts
# ============================================================================
#
# Psi&D
#


class CpsiD0(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'psi', 'D0', tos_1, tos_2, ptmin)


class CpsiDp(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'psi', 'Dp', tos_1, tos_2, ptmin)


class CpsiDs(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'psi', 'Ds', tos_1, tos_2, ptmin)


class CpsiLc(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'psi', 'Lc', tos_1, tos_2, ptmin)
#
# D&D
#


class CD0D0(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'D01', 'D02', tos_1, tos_2, ptmin)


class CD0Dp(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'D0', 'Dp', tos_1, tos_2, ptmin)


class CD0Ds(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'D0', 'Ds', tos_1, tos_2, ptmin)


class CD0Lc(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'D0', 'Lc', tos_1, tos_2, ptmin)


class CDpDp(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'Dp1', 'Dp2', tos_1, tos_2, ptmin)


class CDpDs(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'Dp', 'Ds', tos_1, tos_2, ptmin)


class CDpLc(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'Dp', 'Lc', tos_1, tos_2, ptmin)


class CDsDs(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'Ds1', 'Ds2', tos_1, tos_2, ptmin)


class CDsLc(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'Ds', 'Lc', tos_1, tos_2, ptmin)


class CLcLc(Cuts):

    def __init__(self, tos_1, tos_2, ptmin=2):
        Cuts.__init__(self, 'Lc1', 'Lc2', tos_1, tos_2, ptmin)
