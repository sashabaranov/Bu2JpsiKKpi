#!/usr/bin/env python
# =============================================================================
# @file DiCharm/Samples
#
# Helper module that summarized the samples used for 2xCharm(&Onium)analysis
#
# @author Vanya Belyaev Ivan.Belyaev@cern.ch
# @date   2011-06-19
#
#                   $Revision$
# Last modification $Date$
# by                $Author$
# =============================================================================
__version__ = "$Revision: 124897 $"
__author__ = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__ = "2011-06-07"
# =============================================================================
__all__ = (
    'MC10_samples',  # Monte Carlo data: paths in Bookkeeping DB
    'MC10_from_B',  # Monte Carlo: from B
    'DATA_2011_355pb_samples',  # DATA: LFNs, microDST, tuples & Lumi
    'DATA_2011_1fb_samples',  # DATA: LFNs, microDST, tuples & Lumi
)
# =============================================================================
#
# Patterns for MC10 DB-query paths
mc10_down_pattern = '/MC/MC10/Beam3500GeV-Oct2010-MagDown-Nu2,5/Sim01/Trig0x002e002aFlagged/Reco08/Stripping12Flagged/%d/ALLSTREAMS.DST'
mc10_up_pattern = '/MC/MC10/Beam3500GeV-Oct2010-MagUp-Nu2,5/Sim01/Trig0x002e002aFlagged/Reco08/Stripping12Flagged/%d/ALLSTREAMS.DST'
#
MC10_samples = {
    #
    'MC10-Jpsi/Down': mc10_down_pattern % 24142001,
    'MC10-Jpsi/Up': mc10_up_pattern % 24142001,
    #
    'MC10-Dstar/Down': mc10_down_pattern % 27163003,
    'MC10-Dstar/Up': mc10_up_pattern % 27163003,
    # Ditto:
    'MC10-D*+/Down': mc10_down_pattern % 27163003,
    'MC10-D*+/Up': mc10_up_pattern % 27163003,
    #
    'MC10-Ds/Down': mc10_down_pattern % 23263001,
    'MC10-Ds/Up': mc10_up_pattern % 23263001,
    # Ditto:
    'MC10-Ds+/Down': mc10_down_pattern % 23263001,
    'MC10-Ds+/Up': mc10_up_pattern % 23263001,
    #
    'MC10-D0/Down': mc10_down_pattern % 27163003,
    'MC10-D0/Up': mc10_up_pattern % 27163003,
    #
    'MC10-Dplus/Down': mc10_down_pattern % 21263010,
    'MC10-Dplus/Up': mc10_up_pattern % 21263010,
    # Ditto:
    'MC10-D+/Down': mc10_down_pattern % 21263010,
    'MC10-D+/Up': mc10_up_pattern % 21263010,
    #
    'MC10-LamC/Down': mc10_down_pattern % 25103000,
    'MC10-LamC/Up': mc10_up_pattern % 25103000,
    # Ditto:
    'MC10-Lc+/Down': mc10_down_pattern % 25103000,
    'MC10-Lc+/Up': mc10_up_pattern % 25103000,
}

# ========================================================================
# get the object from read-only db


def _get_data_from_db_(key, dbase='$DICHARMROOT/data/DiCharm.data'):
    """
    Get object from read-only DB 
    """
    import shelve
    import os
    #
    # expand the actual file name
    dbase = os.path.expandvars(dbase)
    dbase = os.path.expanduser(dbase)
    #
    db = shelve.open(dbase, 'r')
    obj = db[key]
    db.close()
    return obj

# ========================================================================
# The samples
DATA_2011_355pb_samples = {
    #
    # original data, collected by  Vladimir...
    #
    # 'CharmMDST/Up'   :  _get_data_from_db_ ( 'DATA:CharmMDST_up_LFNs'   ) ,
    # 'CharmMDST/Down' :  _get_data_from_db_ ( 'DATA:CharmMDST_down_LFNs' ) ,
    #
    # selected microDST for DiCharm modes, Updated...
    #
    'DiCharm/Up':  _get_data_from_db_('DATA:DiCharm_up_LFNs'),
    'DiCharm/Down':  _get_data_from_db_('DATA:DiCharm_down_LFNs'),
    #
    'DiMuon&Charm/Down':  _get_data_from_db_('DATA:DiMuon&Charm_down_LFNs'),
    'DiMuon&Charm/Up':  _get_data_from_db_('DATA:DiMuon&Charm_up_LFNs'),
    #
    # Tuples for trigger studies: Hadrons
    #
    'HadronsTrgTuples/Down':  _get_data_from_db_('DATA:Mesons_Trg_down_tuples'),
    'HadronsTrgTuples/Up':  _get_data_from_db_('DATA:Mesons_Trg_up_tuples'),
    #
    # Tuples for trigger studies: J/psi
    #
    'JpsiTrgTuples/Down':  _get_data_from_db_('DATA:Jpsi_Trg_down_tuples'),
    'JpsiTrgTuples/Up':  _get_data_from_db_('DATA:Jpsi_Trg_up_tuples'),
    #
    # Tuples for mu-PID studies: J/psi
    #
    'JpsiMuPidTuples/Down':  _get_data_from_db_('DATA:Jpsi_muPID_down_LFNs'),
    'JpsiMuPidTuples/Down-Lumi':  _get_data_from_db_('DATA:Jpsi_muPID_down_Lumi'),
    'JpsiMuPidTuples/Up':  _get_data_from_db_('DATA:Jpsi_muPID_up_LFNs'),
    'JpsiMuPidTuples/Up-Lumi':  _get_data_from_db_('DATA:Jpsi_muPID_up_Lumi'),
    #
    # ``signal'' tuples produced with DiCharm.Alg
    #
    'Jpsi&D_Down_Tuples':  _get_data_from_db_('DATA:Jpsi&D_down_tuples'),
    'Jpsi&D_Up_Tuples':  _get_data_from_db_('DATA:Jpsi&D_up_tuples'),
    'D&D_Down_Tuples':  _get_data_from_db_('DATA:D&D_down_tuples'),
    'D&D_Up_Tuples':  _get_data_from_db_('DATA:D&D_up_tuples'),
    #
    # Tuples for J/psi polarization studies: J/psi from B
    'JpsiFromB': _get_data_from_db_('DATA:JpsiB_Polarisation'),
    #
}

# ========================================================================
# The samples
DATA_2011_1fb_samples = {
    #
    '2xDiMuon/Up': _get_data_from_db_('DATA-2011/1fb:2xDiMuon/Up:LFN/1') + _get_data_from_db_('DATA-2011/1fb:2xDiMuon/Up:LFN/2'),
    '2xDiMuon/Down': _get_data_from_db_('DATA-2011/1fb:2xDiMuon/Down:LFN/1') + _get_data_from_db_('DATA-2011/1fb:2xDiMuon/Down:LFN/2'),
    #
    'DiMuon&Charm/Up': _get_data_from_db_('DATA-2011/1fb:DiMuon&Charm/Up:LFN/1') + _get_data_from_db_('DATA-2011/1fb:DiMuon&Charm/Up:LFN/2'),
    'DiMuon&Charm/Down': _get_data_from_db_('DATA-2011/1fb:DiMuon&Charm/Down:LFN/1') + _get_data_from_db_('DATA-2011/1fb:DiMuon&Charm/Down:LFN/2'),
    #
    'DiCharm/Up': _get_data_from_db_('DATA-2011/1fb:DiCharm/Up:LFN/1') + _get_data_from_db_('DATA-2011/1fb:DiCharm/Up:LFN/2'),
    'DiCharm/Down': _get_data_from_db_('DATA-2011/1fb:DiCharm/Down:LFN/1') + _get_data_from_db_('DATA-2011/1fb:DiCharm/Down:LFN/2'),
    #
    #
    # Tuples for trigger studies: J/psi, psi' and Y
    "J/psiTrgTuples":  _get_data_from_db_("DATA-2011/1fb:J/psi_Trg_tuples"),
    "psi'TrgTuples":  _get_data_from_db_("DATA-2011/1fb:psi'_Trg_tuples"),
    "YTrgTuples":  _get_data_from_db_("DATA-2011/1fb:Y_Trg_tuples"),
    #
}
# ========================================================================
# Various MC-tuples
MC10_data = {
    # 'Jpsi_from_B'         :  _get_data_from_db_ ( 'DATA:J/psi_fromB_LFNs'  ) ,
    # 'D0_from_B'           :  _get_data_from_db_ ( 'DATA:D0_fromB_LFNs'     ) ,
    # 'Dp_from_B'           :  _get_data_from_db_ ( 'DATA:Dp_fromB_LFNs'     ) ,
    # 'Ds_from_B'           :  _get_data_from_db_ ( 'DATA:Ds_fromB_LFNs'     ) ,
    # 'Lc_from_B'           :  _get_data_from_db_ ( 'DATA:Lc_fromB_LFNs'     ) ,
    #
    # 'Bd2DX'               :  _get_data_from_db_ ( 'DATA:Bd2DX_Tuples'      ) ,
    # 'Bu2DX'               :  _get_data_from_db_ ( 'DATA:Bu2DX_Tuples'      ) ,
    # 'Bs2DX'               :  _get_data_from_db_ ( 'DATA:Bs2DX_Tuples'      ) ,
    # 'Bi2DX'               :  _get_data_from_db_ ( 'DATA:Bi2DX_Tuples'      ) ,
    # 'Lb2DX'               :  _get_data_from_db_ ( 'DATA:Lb2DX_Tuples'      ) ,
    #
    # 'D0_MCB2'             :  _get_data_from_db_ ( 'DATA:MCB2_D0_Tuples'    ) ,
    # 'Ds_MCB2'             :  _get_data_from_db_ ( 'DATA:MCB2_Ds_Tuples'    ) ,
    # 'Dp_MCB2'             :  _get_data_from_db_ ( 'DATA:MCB2_Dp_Tuples'    ) ,
    # 'Lc_MCB2'             :  _get_data_from_db_ ( 'DATA:MCB2_Lc_Tuples'    ) ,
    #
    'Jpsi_PU':  _get_data_from_db_('MCPU:Jpsi_Tuples'),
    'D0_PU':  _get_data_from_db_('MCPU:D0_Tuples'),
    'Dp_PU':  _get_data_from_db_('MCPU:Dp_Tuples'),
    'Ds_PU':  _get_data_from_db_('MCPU:Ds_Tuples'),
    'Lc_PU':  _get_data_from_db_('MCPU:Lc_Tuples'),
}

# =============================================================================
if '__main__' == __name__:

    print 80 * '*'
    print __doc__
    print ' Author  : ', __author__
    print ' Version : ', __version__
    print ' Date    : ', __date__
    print ' symbols : ', __all__
    print 80 * '*'

# =============================================================================
# The END
# =============================================================================
